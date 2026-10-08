"""Supervised reproduction of the Hopf complex structure on S^{2n−1}×S¹ in the
AMBIENT (chart-free, overlap-free) framework of extra/campaign_unit.py.

The network class is exactly the campaign's CS stage (imported, not copied):

    J = e^{B₂} e^{B₁} J_base e^{−B₁} e^{−B₂},   B_s = P·mlp_s(p)·P   (J² = −P exact)

with J_base a certified ACS (default: the polar-corrected field of a finished
campaign acs unit, --base-from). Supervised Sobolev loss against the closed-form
ambient Hopf structure J_t (src/geometry/hopf_ambient.py):

    L = E_p[ ‖J − J_t‖²_F + λ_D Σ_a ‖D_{Pe_a}(J − J_t)‖²_F ]      (--loss sobolev)
    L = E_p ‖J − J_t‖²_F                                          (--loss value)
    L = E_p ‖N_J‖²   (unsupervised continuation, e.g. from --init-from)  (--loss nij)

D_{Pe_a} = tangential directional derivative along the Parseval frame Pe_a
(Σ_a over it is intrinsic; normal derivatives are extension-dependent and unused).

Reports: relative value / derivative errors (RMS, sup, adversarial sup of the
value error), intrinsic ‖N‖ (mean, sup, adversarial sup), the paper's chart
protocol, and chart-coordinate J / ∂J errors against main's
`hopf_complex_structure_fixed` on the hard-cutoff band — the exact quantities of
main's supervised Phase A, for an apples-to-apples comparison.

Usage:
  python extra/ambient_supervised_hopf.py --verify
  python extra/ambient_supervised_hopf.py --case s3xs1 --base-from DIR --out OUT [--steps 40000] ...
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import os
import sys
import json
import time
import math
import argparse
import torch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from extra.campaign_unit import (Geometry, Net, make_field, load_acs_field, polar_corrected,
                                 acs_residual, nij_norm, adversarial_sup, chart_protocol,
                                 chart_field, dirichlet_energy, DTYPE, BATCH, R_CUT)
from src.geometry.hopf_ambient import hopf_ambient_acs, C_MAIN
from src.geometry.hopf_complex_structure_fixed import hopf_complex_structure
from src.loss.ambient_loss import compute_J_and_DJ
from src.loss.sphere_loss import _nij_tensor

CASES = {'s3xs1': 2, 's5xs1': 3}


# ---------------------------------------------------------------------------
# Target, frame, orientation
# ---------------------------------------------------------------------------

def target_fn(n, c):
    return lambda p: hopf_ambient_acs(p, n, c)


def hopf_frame(p, n):
    """Global ORTHONORMAL frame F_k = Λ₁⁻¹(e_k) = (e_k − q_k q, q_k ∂_θ), k = 1..2n,
    rows = frame vectors in ambient components, shape (B, 2n, 2n+2)."""
    d = 2 * n
    q, w = p[:, :d], p[:, d:]
    eye = torch.eye(d, dtype=p.dtype, device=p.device)
    tang = eye.unsqueeze(0) - q.unsqueeze(2) * q.unsqueeze(1)
    R = torch.stack([-w[:, 1], w[:, 0]], dim=1)
    return torch.cat([tang, q.unsqueeze(2) * R.unsqueeze(1)], dim=2)


def orientation_sign(J, F, n, seed=0):
    """Orientation induced by an ACS, relative to the frame F: sign det[v₁, Jv₁, …, v_n, Jv_n]
    (frame components) for fixed random v_i — well defined for any complex basis."""
    Jf = F @ J @ F.transpose(1, 2)                                # (B, 2n, 2n)
    g = torch.Generator().manual_seed(seed)
    V = torch.randn(2 * n, n, generator=g, dtype=J.dtype)
    cols = []
    for i in range(n):
        v = V[:, i].expand(len(J), -1)
        cols += [v, torch.einsum('bij,bj->bi', Jf, v)]
    return torch.sign(torch.linalg.det(torch.stack(cols, dim=2)))


# ---------------------------------------------------------------------------
# Loss pieces
# ---------------------------------------------------------------------------

def tangential_D(DJ, P):
    """(B,N,N,N) ambient Jacobian [.., k] → (B, a, N, N) directional derivatives along Pe_a."""
    return torch.einsum('bijk,bka->baij', DJ, P)


def sup_loss(geo, J_fn, Jt_fn, x, lam_d, kind):
    P = geo.manifold.projector(x)
    if kind == 'value':
        dJ = J_fn(x) - Jt_fn(x)
        v = (dJ ** 2).sum(dim=(1, 2)).mean()
        return v, {'val': v.item()}
    J, DJ = compute_J_and_DJ(J_fn, x)
    Jt, DJt = compute_J_and_DJ(Jt_fn, x)
    v = ((J - Jt) ** 2).sum(dim=(1, 2)).mean()
    dd = ((tangential_D(DJ - DJt, P)) ** 2).sum(dim=(1, 2, 3)).mean()
    return v + lam_d * dd, {'val': v.item(), 'der': dd.item()}


def unit_loss(geo, J_fn, Jt_fn, x, args):
    if args.loss == 'nij':
        l = (nij_norm(geo, J_fn, x) ** 2).mean()
        return l, {'nij': l.item()}
    return sup_loss(geo, J_fn, Jt_fn, x, args.lam_d, args.loss)


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def ambient_errors(geo, J_fn, Jt_fn, n_pts, chunk=500):
    """Relative value / tangential-derivative errors and intrinsic ‖N‖ on fresh points."""
    acc = {k: [] for k in ('dv', 'tv', 'dd', 'td', 'nij', 'acs')}
    for s in range(0, n_pts, chunk):
        x = geo.sample(min(chunk, n_pts - s))
        P = geo.manifold.projector(x)
        J, DJ = compute_J_and_DJ(J_fn, x)
        Jt, DJt = compute_J_and_DJ(Jt_fn, x)
        J, DJ = J.detach(), DJ.detach()
        acc['dv'].append((J - Jt).norm(dim=(1, 2)))
        acc['tv'].append(Jt.norm(dim=(1, 2)))
        acc['dd'].append(tangential_D(DJ - DJt, P).reshape(len(x), -1).norm(dim=1))
        acc['td'].append(tangential_D(DJt, P).reshape(len(x), -1).norm(dim=1))
        acc['nij'].append(nij_norm(geo, J_fn, x).detach())
        with torch.no_grad():
            acc['acs'].append(acs_residual(geo, J_fn, x))
    a = {k: torch.cat(v) for k, v in acc.items()}
    rv, rd = a['dv'] / a['tv'], a['dd'] / a['td']
    return {
        'rel_val_rms': (a['dv'].pow(2).mean() / a['tv'].pow(2).mean()).sqrt().item(),
        'rel_val_sup': rv.max().item(),
        'rel_der_rms': (a['dd'].pow(2).mean() / a['td'].pow(2).mean()).sqrt().item(),
        'rel_der_sup': rd.max().item(),
        'val_mse': a['dv'].pow(2).mean().item(), 'der_mse': a['dd'].pow(2).mean().item(),
        'nij_mean': a['nij'].mean().item(), 'nij_sup': a['nij'].max().item(),
        'acs_sup': a['acs'].max().item(),
    }


def main_band_points(n, npts=3000):
    """EXACT replica of main's held-out evaluation set (extra/diagnostics/precision_2026_06_16/
    precision_lbfgs_only.py::eval_held_out, used for every main-branch Hopf number): seeds 99999,
    SphereProduct(2n−1, 1).sample_points(npts), keep points with |x| < 0.85 in BOTH charts of each
    factor (the overlap band). Returns patch-(0,0) coordinates (xa, xb)."""
    import numpy as np
    from src.geometry.sphere import EvenSphere
    from src.geometry.sphere_product import SphereProduct
    torch.manual_seed(99999); np.random.seed(99999)
    xa, xb = SphereProduct(2 * n - 1, 1, device='cpu').sample_points(npts)
    xa, xb = xa.to(DTYPE), xb.to(DTYPE)
    xa2, xb2 = EvenSphere.patch_change_coords(xa), EvenSphere.patch_change_coords(xb)
    keep = (xa.norm(dim=1) < 0.85) & (xb.norm(dim=1) < 0.85) & (xa2.norm(dim=1) < 0.85) & (xb2.norm(dim=1) < 0.85)
    return xa[keep], xb[keep]


def chart_comparison(J_fn, n, n_pts=2000, seed=11, band='rcut'):
    """Main-branch Phase A quantities, per patch: J-fit = E‖J_chart − J_main‖²_F,
    ∂J-fit = E‖∂J_chart − ∂J_main‖²_F, flat-chart ‖N‖ mean.
    band='main': main's exact held-out set (overlap band, see main_band_points) — the one to compare
    with main's numbers; band='rcut': uniform in the ball |x_α|,|x_β| ≤ R_CUT of each patch."""
    from src.geometry.sphere import EvenSphere
    m = 2 * n - 1
    g = torch.Generator().manual_seed(seed)
    if band == 'main':
        xa0, xb0 = main_band_points(n)
    out = {}
    for ia in (0, 1):
        for ib in (0, 1):
            if band == 'main':
                xa = xa0 if ia == 0 else EvenSphere.patch_change_coords(xa0)
                xb = xb0 if ib == 0 else EvenSphere.patch_change_coords(xb0)
            else:
                xa = torch.randn(n_pts, m, generator=g, dtype=DTYPE)
                xa = xa / xa.norm(dim=1, keepdim=True) * torch.rand(n_pts, 1, generator=g, dtype=DTYPE) ** (1.0 / m) * R_CUT
                xb = (torch.rand(n_pts, 1, generator=g, dtype=DTYPE) * 2 - 1) * R_CUT
            xc = torch.cat([xa, xb], dim=1)
            n_pts = len(xc)
            J_of = chart_field(J_fn, (m, 1), (ia, ib))
            Jm_of = lambda x, ia=ia, ib=ib: hopf_complex_structure(x[:, :m], x[:, m:], ia, ib)
            Jc, dJc, Jm, dJm = [], [], [], []
            for s in range(0, n_pts, 250):
                xs = xc[s:s + 250]
                Jc.append(J_of(xs).detach())
                dJc.append(torch.func.vmap(torch.func.jacfwd(lambda z: J_of(z.unsqueeze(0))[0]))(xs).detach())
                Jm.append(Jm_of(xs).detach())
                dJm.append(torch.func.vmap(torch.func.jacfwd(lambda z: Jm_of(z.unsqueeze(0))[0]))(xs).detach())
            Jc, dJc, Jm, dJm = map(torch.cat, (Jc, dJc, Jm, dJm))
            Nn = torch.sqrt((_nij_tensor(Jc, dJc) ** 2).sum(dim=(1, 2, 3)) + 1e-20)
            out[f'{ia}{ib}'] = {
                'J_fit': ((Jc - Jm) ** 2).sum(dim=(1, 2)).mean().item(),
                'dJ_fit': ((dJc - dJm) ** 2).sum(dim=(1, 2, 3)).mean().item(),
                'rel_J': ((Jc - Jm).norm(dim=(1, 2)).pow(2).mean() / Jm.norm(dim=(1, 2)).pow(2).mean()).sqrt().item(),
                'rel_dJ': ((dJc - dJm).reshape(n_pts, -1).norm(dim=1).pow(2).mean()
                           / dJm.reshape(n_pts, -1).norm(dim=1).pow(2).mean()).sqrt().item(),
                'nij_flat_mean': Nn.mean().item(), 'nij_flat_max': Nn.max().item(),
                'dJ_main_rms': dJm.reshape(n_pts, -1).norm(dim=1).pow(2).mean().sqrt().item(),
            }
    keys = ['J_fit', 'dJ_fit', 'rel_J', 'rel_dJ', 'nij_flat_mean']
    out['mean'] = {k: sum(out[p][k] for p in ('00', '01', '10', '11')) / 4.0 for k in keys}
    return out


def full_report(geo, J_fn, Jt_fn, n, args):
    rep = {'ambient': ambient_errors(geo, J_fn, Jt_fn, args.n_dense)}
    rel_fn = lambda x: (J_fn(x) - Jt_fn(x)).norm(dim=(1, 2)) / Jt_fn(x).norm(dim=(1, 2))
    rep['adv_rel_val_sup'], _ = adversarial_sup(geo, rel_fn, args.adv_starts, args.adv_steps)
    rep['adv_nij_sup'], _ = adversarial_sup(geo, lambda x: nij_norm(geo, J_fn, x),
                                            args.adv_starts // 2, args.adv_steps // 2)
    rep['chart_cs'], rep['chart_overlap'] = chart_protocol(geo, J_fn, args.n_chart)
    rep['chart_vs_main'] = chart_comparison(J_fn, n, args.n_chart)
    rep['main_protocol'] = chart_comparison(J_fn, n, band='main')
    return rep


# ---------------------------------------------------------------------------
# Verification of the target (Stage 1)
# ---------------------------------------------------------------------------

def verify(c):
    torch.manual_seed(0)
    for case, n in CASES.items():
        geo = Geometry(case)
        Jt = target_fn(n, c)
        x = geo.sample(3000)
        P = geo.manifold.projector(x)
        Jr = hopf_ambient_acs(x, n, c, project=False)
        F = hopf_frame(x, n)
        print(f"== {case} (c = {c:.6f})")
        print(f"   tangency {(P @ Jr @ P - Jr).abs().max():.1e} | ‖J²+P‖ max {acs_residual(geo, Jt, x).max():.1e} "
              f"| intrinsic ‖N‖ max {nij_norm(geo, Jt, x).max():.1e}")
        print(f"   frame orthonormality {(F @ F.transpose(1, 2) - torch.eye(2 * n, dtype=DTYPE)).abs().max():.1e} "
              f"| orientation sign of J_t: {orientation_sign(Jt(x), F, n).unique().tolist()}")
        print(f"   ‖J‖_F mean {Jt(x).norm(dim=(1, 2)).mean():.4f} | Dirichlet energy {dirichlet_energy(geo, Jt, x[:1000]).item():.3f}")
        cmp = chart_comparison(Jt, n, 400)
        for k in ('00', '01', '10', '11'):
            print(f"   chart {k}: rel J {cmp[k]['rel_J']:.1e}  rel ∂J {cmp[k]['rel_dJ']:.1e}  "
                  f"flat ‖N‖ {cmp[k]['nij_flat_mean']:.1e}  (‖∂J_main‖ rms {cmp[k]['dJ_main_rms']:.2f})")


# ---------------------------------------------------------------------------
# Base ACS
# ---------------------------------------------------------------------------

def build_base(geo, n, base_from, base_sign, Jt_fn, tag, mode=None):
    """The acs unit's own certified field J* (Newton–Schulz polar correction), orientation-matched
    to the target. Returns (base_fn, info, mode); `mode` forces the fix (used by --eval-only)."""
    N = geo.N
    K_base = load_acs_field(geo, base_from)
    adv_base, _ = adversarial_sup(geo, lambda x: acs_residual(geo, K_base, x), 256, 150)
    if adv_base >= 1.0:
        raise SystemExit(f"[{tag}] base in {base_from} is not certified (adversarial sup {adv_base:.3f})")
    J_star = polar_corrected(geo, K_base)
    xo = geo.sample(2000)
    F = hopf_frame(xo, n)
    with torch.no_grad():
        o_t = orientation_sign(Jt_fn(xo), F, n)
        o_b = orientation_sign(J_star(xo), F, n)
    if len(o_t.unique()) != 1 or len(o_b.unique()) != 1:
        raise SystemExit(f"[{tag}] orientation sign not constant: target {o_t.unique()}, base {o_b.unique()}")
    same = bool(o_t[0] == o_b[0])
    # Orientation fix. The conjugation family preserves the base's orientation (det e^B > 0), so a
    # base of the wrong orientation can never reach J_t. n odd: −J* flips the orientation.
    # n even (−J* does not): transport J* by the orientation-reversing ISOMETRY S: θ ↦ −θ of the
    # round product, S(q, w₀, w₁) = (q, −w₀, w₁) (w = (sin θ, −cos θ)); S·J*(S p)·S is again a
    # certified ACS (all residuals are S-invariant) with the opposite orientation.
    S = torch.ones(N, dtype=DTYPE)
    S[2 * n] = -1.0
    S = torch.diag(S)
    if mode is None:
        if base_sign == 'auto':
            mode = 'identity' if same else ('negate' if n % 2 == 1 else 'reflect')
        else:
            mode = {'+1': 'identity', '-1': 'negate'}[base_sign]
    if mode == 'identity':
        base_fn = lambda x: J_star(x)
    elif mode == 'negate':
        base_fn = lambda x: -J_star(x)
    else:
        base_fn = lambda x: S @ J_star(x @ S) @ S
    with torch.no_grad():
        o_used = orientation_sign(base_fn(xo), F, n)
    if o_used[0] != o_t[0]:
        raise SystemExit(f"[{tag}] orientation of the base ({mode}) differs from the target's")
    xb = geo.sample(1000)
    info = {'dir': base_from, 'adv_acs_sup_raw': adv_base, 'orientation_fix': mode,
            'acs_sup_used': acs_residual(geo, base_fn, xo).max().item(),
            'orientation_target': o_t[0].item(), 'orientation_base_raw': o_b[0].item(),
            'orientation_base_used': o_used[0].item(),
            'nij_mean': nij_norm(geo, base_fn, xb).mean().item(),
            'rel_dist_to_target': ((base_fn(xb) - Jt_fn(xb)).norm(dim=(1, 2))
                                   / Jt_fn(xb).norm(dim=(1, 2))).mean().item()}
    return base_fn, info, mode


def eval_only(run_dir, args):
    """Re-evaluate a run → eval.json. Uses final.pt (finished run: same base, same orientation fix);
    if absent, the latest state.pt checkpoint of an UNFINISHED run, with --case/--base-from/--hidden/
    --stages/--c taken from the command line (orientation fix re-derived, 'auto')."""
    fpath = os.path.join(run_dir, 'final.pt')
    if os.path.exists(fpath):
        ck = torch.load(fpath, map_location='cpu', weights_only=False)
        partial_step = None
    else:
        st = torch.load(os.path.join(run_dir, 'state.pt'), map_location='cpu', weights_only=False)
        assert args.case and args.base_from, "partial run: pass --case and --base-from"
        d_ = 2 * CASES[args.case]
        ck = {'model': st['model'], 'case': args.case, 'c': args.c, 'base_from': args.base_from,
              'orientation_fix': None, 'stages': args.stages,
              'hidden': [int(h) for h in args.hidden.split(',')] if args.hidden else
              ([256, 256, 256, 256] if d_ >= 6 else [128, 128, 128])}
        partial_step = st['step']
        print(f"[eval] PARTIAL run {run_dir}: state.pt at step {partial_step}", flush=True)
    case, n = ck['case'], CASES[ck['case']]
    geo = Geometry(case, 'ambient')
    Jt_fn = target_fn(n, ck['c'])
    base_fn, base_info, _ = build_base(geo, n, ck['base_from'], 'auto', Jt_fn, run_dir,
                                       mode=ck['orientation_fix'])
    model = Net(geo.N, geo.N, ck['stages'], ck['hidden'], 0.5).to(DTYPE)
    model.load_state_dict(ck['model'])
    model.eval()
    J_fn, _ = make_field(geo, 'cs', model, base_fn=base_fn)
    rep = full_report(geo, J_fn, Jt_fn, n, args)
    name = 'eval.json' if partial_step is None else f'eval_partial_step{partial_step}.json'
    with open(os.path.join(run_dir, name), 'w') as f:
        json.dump({'run': run_dir, 'partial_step': partial_step, 'base': base_info, **rep}, f, indent=2)
    print(json.dumps({'ambient': rep['ambient'], 'adv_rel_val_sup': rep['adv_rel_val_sup'],
                      'adv_nij_sup': rep['adv_nij_sup'], 'chart_cs': rep['chart_cs'],
                      'main_protocol_mean': rep['main_protocol']['mean']}, indent=1))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    ap.add_argument('--case', choices=list(CASES))
    ap.add_argument('--c', type=float, default=C_MAIN, help='Hopf log-period constant (default 1/2π = main)')
    ap.add_argument('--base-from', default=None, help='finished campaign acs unit dir (certified base ACS)')
    ap.add_argument('--base-sign', default='auto', choices=['auto', '+1', '-1'],
                    help="auto: match the target's orientation (flip −J* when n is odd; error when impossible)")
    ap.add_argument('--init-from', default=None, help='final.pt of a previous run of this script (warm start)')
    ap.add_argument('--loss', default='sobolev', choices=['sobolev', 'value', 'nij'])
    ap.add_argument('--lam-d', type=float, default=1.0)
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--out', default=None)
    ap.add_argument('--steps', type=int, default=None)
    ap.add_argument('--hidden', default=None)
    ap.add_argument('--stages', type=int, default=2)
    ap.add_argument('--out-scale', type=float, default=0.5, help='last-layer init scale (campaign CS: 0.5)')
    ap.add_argument('--batch', type=int, default=BATCH)
    ap.add_argument('--lr', type=float, default=1e-3)
    ap.add_argument('--lr-min', type=float, default=1e-5)
    ap.add_argument('--lbfgs', type=int, default=400)
    ap.add_argument('--lbfgs-chunk', type=int, default=100)
    ap.add_argument('--lbfgs-batch', type=int, default=2048)
    ap.add_argument('--eval-every', type=int, default=1000)
    ap.add_argument('--ckpt-every', type=int, default=1000)
    ap.add_argument('--n-dense', type=int, default=4000)
    ap.add_argument('--n-chart', type=int, default=1500)
    ap.add_argument('--adv-starts', type=int, default=128)
    ap.add_argument('--adv-steps', type=int, default=100)
    ap.add_argument('--eval-only', default=None, help='run dir of a finished run: re-evaluate → eval.json')
    args = ap.parse_args()

    if args.verify:
        verify(args.c)
        return
    if args.eval_only:
        torch.manual_seed(args.seed)
        eval_only(args.eval_only, args)
        return
    assert args.case and args.out and args.base_from, "--case, --out and --base-from are required"
    os.makedirs(args.out, exist_ok=True)
    torch.manual_seed(args.seed)
    n = CASES[args.case]
    geo = Geometry(args.case, 'ambient')
    N, d = geo.N, geo.d
    steps = args.steps or (40000 if d == 4 else 80000)
    hidden = [int(h) for h in args.hidden.split(',')] if args.hidden else \
        ([256, 256, 256, 256] if d >= 6 else [128, 128, 128])
    tag = f"{args.case}/{args.loss}/s{args.seed}"
    Jt_fn = target_fn(n, args.c)

    base_fn, base_info, mode = build_base(geo, n, args.base_from, args.base_sign, Jt_fn, tag)
    print(f"[{tag}] base {base_info}", flush=True)

    model = Net(N, N, args.stages, hidden, args.out_scale).to(DTYPE)
    J_fn, meta = make_field(geo, 'cs', model, base_fn=base_fn)
    if args.init_from:
        ck = torch.load(args.init_from, map_location='cpu', weights_only=False)
        model.load_state_dict(ck['model'])
        print(f"[{tag}] warm start from {args.init_from}", flush=True)

    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=steps, eta_min=args.lr_min)
    state_path = os.path.join(args.out, 'state.pt')
    start, elapsed_prev = 0, 0.0
    if os.path.exists(state_path):
        st = torch.load(state_path, map_location='cpu', weights_only=False)
        model.load_state_dict(st['model']); opt.load_state_dict(st['opt']); sched.load_state_dict(st['sched'])
        start, elapsed_prev = st['step'], st['elapsed']
        print(f"[{tag}] resumed at step {start}", flush=True)
    t0 = time.time()
    elapsed = lambda: elapsed_prev + time.time() - t0
    def save_state(s):    # atomic: a kill mid-write (timeout, scheduler time limit) never truncates state.pt
        torch.save({'model': model.state_dict(), 'opt': opt.state_dict(),
                    'sched': sched.state_dict(), 'step': s, 'elapsed': elapsed()}, state_path + '.tmp')
        os.replace(state_path + '.tmp', state_path)
    progress = open(os.path.join(args.out, 'progress.jsonl'), 'a')
    if start == 0:
        ev0 = ambient_errors(geo, J_fn, Jt_fn, 1000)
        progress.write(json.dumps({'step': 0, **ev0}) + '\n'); progress.flush()
        print(f"[{tag}] init {ev0}", flush=True)

    x_val = geo.sample(1500)
    for step in range(start, steps):
        x = geo.sample(args.batch)
        loss, info = unit_loss(geo, J_fn, Jt_fn, x, args)
        if not bool(torch.isfinite(loss)):
            print(f"[{tag}] non-finite loss at step {step}; skipped", flush=True)
            continue
        opt.zero_grad(); loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step(); sched.step()
        if (step + 1) % args.ckpt_every == 0:
            save_state(step + 1)
        if (step + 1) % args.eval_every == 0 or step == steps - 1:
            ev = ambient_errors(geo, J_fn, Jt_fn, 1000)
            rec = {'step': step + 1, 'elapsed_min': elapsed() / 60, 'lr': sched.get_last_lr()[0], 'batch': info, **ev}
            progress.write(json.dumps(rec) + '\n'); progress.flush()
            print(f"[{tag}] step {step+1:6d} | {info} | rel val {ev['rel_val_rms']:.2e} (sup {ev['rel_val_sup']:.1e}) "
                  f"| rel der {ev['rel_der_rms']:.2e} | ‖N‖ {ev['nij_mean']:.3e} (sup {ev['nij_sup']:.2e}) "
                  f"| {elapsed()/60:.1f} min", flush=True)
    save_state(steps)

    # ---- chunked L-BFGS finisher, fresh batch per chunk, revert if held-out loss worsens ----
    lb_log = []
    if args.lbfgs > 0:
        best = unit_loss(geo, J_fn, Jt_fn, x_val, args)[0].item()
        done = 0
        while done < args.lbfgs:
            xf = geo.sample(args.lbfgs_batch)
            snap = {k: v.clone() for k, v in model.state_dict().items()}
            lb = torch.optim.LBFGS(model.parameters(), lr=1.0, max_iter=args.lbfgs_chunk,
                                   history_size=100, line_search_fn='strong_wolfe')

            def closure():
                lb.zero_grad()
                l, _ = unit_loss(geo, J_fn, Jt_fn, xf, args)
                l.backward()
                return l
            try:
                lb.step(closure)
            except Exception as e:
                print(f"[{tag}] L-BFGS chunk failed ({type(e).__name__}); reverting", flush=True)
                model.load_state_dict(snap)
                break
            val = unit_loss(geo, J_fn, Jt_fn, x_val, args)[0].item()
            reverted = not (val <= best)
            if reverted:
                model.load_state_dict(snap)
            else:
                best = val
            done += args.lbfgs_chunk
            rec = {'phase': 'lbfgs', 'iters': done, 'val_loss': val, 'best': best, 'reverted': reverted,
                   'elapsed_min': elapsed() / 60}
            lb_log.append(rec)
            progress.write(json.dumps(rec) + '\n'); progress.flush()
            print(f"[{tag}] lbfgs {done}/{args.lbfgs} | held-out loss {val:.3e} (best {best:.3e})"
                  f"{' REVERTED' if reverted else ''}", flush=True)

    torch.save({'model': model.state_dict(), 'case': args.case, 'hidden': hidden, 'stages': args.stages,
                'base_from': args.base_from, 'orientation_fix': mode, 'c': args.c, 'loss': args.loss},
               os.path.join(args.out, 'final.pt'))
    rep = full_report(geo, J_fn, Jt_fn, n, args)
    result = {'case': args.case, 'n': n, 'c': args.c, 'loss': args.loss, 'lam_d': args.lam_d,
              'seed': args.seed, 'hidden': hidden, 'stages': args.stages, 'out_scale': args.out_scale,
              'batch': args.batch, 'lr': args.lr, 'steps': steps, 'init_from': args.init_from,
              'n_params': sum(p.numel() for p in model.parameters()), 'base': base_info,
              'lbfgs': lb_log, 'elapsed_min': elapsed() / 60, **rep}
    with open(os.path.join(args.out, 'result.json'), 'w') as f:
        json.dump(result, f, indent=2)
    a = rep['ambient']
    print(f"[{tag}] RESULT | rel val rms {a['rel_val_rms']:.2e} sup {a['rel_val_sup']:.2e} (adv {rep['adv_rel_val_sup']:.2e}) "
          f"| rel der rms {a['rel_der_rms']:.2e} | intrinsic ‖N‖ mean {a['nij_mean']:.3e} sup {a['nij_sup']:.3e} "
          f"(adv {rep['adv_nij_sup']:.3e}) | chart ‖N‖ {rep['chart_cs']:.3e} | chart-vs-main {rep['chart_vs_main']['mean']}",
          flush=True)


if __name__ == '__main__':
    main()
