"""Uniform chart-free PINN campaign — one (case, task, seed) unit.

Formulation is chosen by topology (src/geometry/explicit_acs.py):
  parallelizable  (s1xs1, s3xs1, s3xs3, s5xs1)  → global-frame Nijenhuis, no charts
  non-parallelizable (s2, s4, s6, s2xs2, s4xs2) → ambient-projector Nijenhuis, no charts

Tasks (unsupervised PINN, identical across formulations):
  acs  free endomorphism field K (frame: K = net(x); ambient: K = P·net(x)·P),
       loss = mean ‖K² + Id_T‖²_F. From scratch, no base point. Certified a
       posteriori by the polar correction J* = K(−K²_T)^{−1/2} when
       sup ‖K²+Id_T‖ < 1.
  cs   exact-ACS parametrization J = e^{B_s}⋯e^{B_1} J_base e^{−B_1}⋯e^{−B_s}
       (B_k tangential net outputs, initialised at O(1) ⇒ a generic
       NON-integrable start), loss = mean ‖N(J)‖² only. J_base is the
       canonical J₀ (frame cases) or a classical explicit ACS (ambient cases);
       s4 has no ACS, so its cs unit runs the free-field two-term loss
       (w·‖K²+P‖² + ‖N‖²) as the "nothing exists" control.

Every unit is restartable (state.pt), time-guarded (--wall-minutes), logs
progress.jsonl, and ends with a dense evaluation + the ORIGINAL chart-protocol
scores (frame fields are bridged to ambient fields, then to ball charts), in
result.json.

Usage:
  python extra/campaign_unit.py --case s3xs3 --task cs --seed 7 --out DIR \
      [--steps 40000] [--hidden 256,256,256,256] [--stages 2] [--rar] \
      [--wall-minutes 4200] [--eval-every 2000] [--ckpt-every 1000]
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import os
import sys
import json
import time
import argparse
import torch
import torch.nn as nn

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from src.geometry.frames import FRAMES
from src.geometry.explicit_acs import (ALL_CASES, factor_dims, is_frame_case,
                                       ambient_manifold, base_acs)
from src.loss.frame_loss import frame_nij_norm
from src.loss.ambient_loss import (ambient_nij_norm, compute_J_and_DJ,
                                   compute_P_and_DP)
from src.geometry.sphere import EvenSphere, RADIAL_MIDPOINT
from src.geometry.sphere_product import SphereProduct
from src.loss.sphere_loss import _nij_tensor, OverlapLoss
from src.loss.product_loss import ProductOverlapLoss

DTYPE = torch.float64
BATCH = 192
OVERLAP_SIZE = 0.2
R_CUT = RADIAL_MIDPOINT + OVERLAP_SIZE

# Results table of the paper (before this work), for side-by-side reporting.
TABLE_BEFORE = {
    's2':    {'acs': 1.63e-4, 'cs': 5.27e-4},
    's4':    {'acs': 2.83,    'cs': 8.50},
    's6':    {'acs': 3.57e-2, 'cs': 9.87},
    's1xs1': {'acs': 5.81e-5, 'cs': 1.90e-4},
    's2xs2': {'acs': 9.27e-3, 'cs': 4.04e-2},
    's3xs1': {'acs': 1.11e-2, 'cs': 2.00e-1},
    's3xs3': {'acs': 1.61e-1, 'cs': 7.63},
    's4xs2': {'acs': 1.91,    'cs': 6.28},
    's5xs1': {'acs': 7.98e-2, 'cs': 6.42},
}


# ---------------------------------------------------------------------------
# Geometry wrapper
# ---------------------------------------------------------------------------

class Geometry:
    def __init__(self, case):
        self.case = case
        self.dims = factor_dims(case)
        self.d = sum(self.dims)
        if is_frame_case(case):
            self.kind = 'frame'
            self.frame = FRAMES[case]()
            self.N = self.frame.ambient_dim
            probe = self.frame.sample(200)
            c = self.frame.structure_functions(probe)
            self.c_const = c[:1] if (c - c.mean(0, keepdim=True)).abs().max() < 1e-10 else None
        else:
            self.kind = 'ambient'
            self.manifold = ambient_manifold(case)
            self.N = self.manifold.ambient_dim
            self.c_const = None

    def sample(self, n):
        return self.frame.sample(n) if self.kind == 'frame' else self.manifold.sample(n)

    def structure(self, x):
        if self.c_const is not None:
            return self.c_const.expand(len(x), self.d, self.d, self.d)
        return self.frame.structure_functions(x)


def canonical_J0(d):
    J0 = torch.zeros(d, d, dtype=DTYPE)
    for k in range(d // 2):
        J0[2 * k, 2 * k + 1] = -1.0
        J0[2 * k + 1, 2 * k] = 1.0
    return J0


# ---------------------------------------------------------------------------
# Network + field construction
# ---------------------------------------------------------------------------

class Net(nn.Module):
    def __init__(self, n_in, n_mat, stages, hidden, out_scale):
        super().__init__()
        layers, prev = [], n_in
        for h in hidden:
            layers += [nn.Linear(prev, h), nn.GELU()]
            prev = h
        layers += [nn.Linear(prev, stages * n_mat * n_mat)]
        self.net = nn.Sequential(*layers)
        with torch.no_grad():
            self.net[-1].weight.mul_(out_scale)
            self.net[-1].bias.normal_(0.0, out_scale)
        self.n_mat, self.stages = n_mat, stages

    def forward(self, x):
        return self.net(x).reshape(-1, self.stages, self.n_mat, self.n_mat)


def make_field(geo, task, model):
    """Returns (J_fn, meta). J_fn maps ambient points → matrix field:
    frame cases: (B, d, d) in frame components; ambient cases: (B, N, N)."""
    case = geo.case
    if geo.kind == 'frame':
        d = geo.d
        if task == 'acs':
            return (lambda x: model(x)[:, 0]), {'exact_acs': False}
        J0 = canonical_J0(d)

        def J_fn(x):
            raw = model(x)
            J = J0.unsqueeze(0).expand(len(x), d, d)
            for s in range(model.stages):
                Bs = raw[:, s]
                J = torch.matrix_exp(Bs) @ J @ torch.matrix_exp(-Bs)
            return J
        return J_fn, {'exact_acs': True}

    proj = geo.manifold.projector
    if task == 'acs' or case == 's4':
        def K_fn(x):
            P = proj(x)
            return P @ model(x)[:, 0] @ P
        return K_fn, {'exact_acs': False}

    def J_fn(x):
        P = proj(x)
        raw = model(x)
        J = base_acs(case, x, proj)
        for s in range(model.stages):
            Bs = P @ raw[:, s] @ P
            J = torch.matrix_exp(Bs) @ J @ torch.matrix_exp(-Bs)
        return J
    return J_fn, {'exact_acs': True}


# ---------------------------------------------------------------------------
# Residuals
# ---------------------------------------------------------------------------

def identity_T(geo, x):
    if geo.kind == 'frame':
        return torch.eye(geo.d, dtype=DTYPE).expand(len(x), geo.d, geo.d)
    return geo.manifold.projector(x)


def acs_residual(geo, J_fn, x):
    K = J_fn(x)
    return torch.sqrt(((K @ K + identity_T(geo, x)) ** 2).sum(dim=(1, 2)) + 1e-20)


def nij_norm(geo, J_fn, x):
    """Per-sample intrinsic Nijenhuis norm of the field (frame or ambient)."""
    if geo.kind == 'frame':
        def single(xs):
            return J_fn(xs.unsqueeze(0)).squeeze(0)
        J = J_fn(x)
        dJ_amb = torch.func.vmap(torch.func.jacfwd(single))(x)            # (B,d,d,N)
        E = geo.frame.frame_fields(x)
        dJ_frame = torch.einsum('beaN,bcN->beac', dJ_amb, E)
        return frame_nij_norm(J, dJ_frame, geo.structure(x))
    J, DJ = compute_J_and_DJ(J_fn, x)
    P, DP = compute_P_and_DP(geo.manifold, x)
    return ambient_nij_norm(J, DJ, P, DP)


def unit_loss(geo, task, J_fn, x, w_c=10.0):
    if task == 'acs':
        r = acs_residual(geo, J_fn, x)
        return (r ** 2).mean(), {'acs': (r ** 2).mean().item()}
    nij = (nij_norm(geo, J_fn, x) ** 2).mean()
    if geo.case == 's4':                    # no ACS exists: two-term control
        r = acs_residual(geo, J_fn, x)
        acs = (r ** 2).mean()
        return w_c * acs + nij, {'acs': acs.item(), 'nij': nij.item()}
    return nij, {'nij': nij.item()}


# ---------------------------------------------------------------------------
# Polar correction (ACS certification)
# ---------------------------------------------------------------------------

def polar_corrected(geo, K_fn):
    def J_star(x):
        K = K_fn(x)
        IdT = identity_T(geo, x)
        eye = torch.eye(K.shape[-1], dtype=K.dtype).expand_as(K)
        S = -K @ K + (eye - IdT)
        w, V = torch.linalg.eigh(S)
        return K @ (V @ torch.diag_embed(w.clamp(min=1e-12) ** -0.5) @ V.transpose(1, 2))
    return J_star


# ---------------------------------------------------------------------------
# Chart-protocol evaluation (the paper's original metric)
# ---------------------------------------------------------------------------

def ambient_field(geo, J_fn):
    """Frame fields are bridged to ambient tangential fields:
    J_amb = Eᵀ J_frame (E Eᵀ)⁻¹ E."""
    if geo.kind == 'ambient':
        return J_fn

    def J_amb(p):
        E = geo.frame.frame_fields(p)                                       # (B,d,N)
        G = E @ E.transpose(1, 2)
        Jf = J_fn(p)
        return E.transpose(1, 2) @ Jf @ torch.linalg.inv(G) @ E
    return J_amb


def factor_embedding(x, chart):
    if chart == 0:
        return EvenSphere.ball_to_sphere(x, from_north=True)
    p = EvenSphere.ball_to_sphere(x, from_north=False)
    return torch.cat([-p[:, :-1], p[:, -1:]], dim=1)


def chart_field(J_amb, dims, charts):
    def emb_single(xc):
        parts, off = [], 0
        for m, ch in zip(dims, charts):
            parts.append(factor_embedding(xc[off:off + m].unsqueeze(0), ch).squeeze(0))
            off += m
        return torch.cat(parts)

    def J_of(x_chart):
        p = torch.func.vmap(emb_single)(x_chart)
        D = torch.func.vmap(torch.func.jacfwd(emb_single))(x_chart)
        Ja = J_amb(p)
        D_pinv = torch.linalg.inv(D.transpose(1, 2) @ D) @ D.transpose(1, 2)
        return D_pinv @ Ja @ D
    return J_of


def chart_nij(J_of, x_chart, dims, chunk=400):
    off, mask = 0, torch.ones(len(x_chart), dtype=torch.bool)
    for m in dims:
        mask &= x_chart[:, off:off + m].norm(dim=1) <= R_CUT
        off += m
    xa = x_chart[mask]
    d = sum(dims)
    norms = []
    for s in range(0, len(xa), chunk):
        xr = xa[s:s + chunk].detach().clone().requires_grad_(True)
        J = J_of(xr)
        dJ = torch.zeros(len(xr), d, d, d, dtype=DTYPE)
        for i in range(d):
            for j in range(d):
                g, = torch.autograd.grad(J[:, i, j].sum(), xr, retain_graph=True)
                dJ[:, i, j, :] = g
        N = _nij_tensor(J.detach(), dJ)
        norms.append(torch.sqrt((N ** 2).sum(dim=(1, 2, 3)) + 1e-20))
    v = torch.cat(norms) if norms else torch.zeros(1, dtype=DTYPE)
    return v.mean().item(), v.max().item()


def chart_protocol(geo, J_fn, n_test):
    """Returns (chart_cs = mean hard-cutoff Nij over patches, chart_acs = overlap)."""
    J_amb = ambient_field(geo, J_fn)
    dims = geo.dims
    if len(dims) == 1:
        m = dims[0]
        sphere = EvenSphere(dimension=m)
        x1 = sphere.sample_points_patch1(n_test).to(DTYPE)
        x2 = EvenSphere.patch_change_coords(x1)
        J1f, J2f = chart_field(J_amb, dims, (0,)), chart_field(J_amb, dims, (1,))
        n1, _ = chart_nij(J1f, x1, dims)
        n2, _ = chart_nij(J2f, x2, dims)
        with torch.no_grad():
            ov, _ = OverlapLoss(overlap_size=OVERLAP_SIZE).to(DTYPE)(
                J1f(x1), J2f(x2), x1, x2, hard_cutoff=True)
        return 0.5 * (n1 + n2), ov.item()
    ma, mb = dims
    sp = SphereProduct(ma, mb)
    xa0, xb0 = sp.sample_points(n_test)
    xa0, xb0 = xa0.to(DTYPE), xb0.to(DTYPE)
    pa = {0: xa0, 1: EvenSphere.patch_change_coords(xa0)}
    pb = {0: xb0, 1: EvenSphere.patch_change_coords(xb0)}
    nijs, J_dict = [], {}
    for ia, ib in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        J_of = chart_field(J_amb, dims, (ia, ib))
        xc = torch.cat([pa[ia], pb[ib]], dim=1)
        n_mean, _ = chart_nij(J_of, xc, dims)
        nijs.append(n_mean)
        with torch.no_grad():
            J_dict[(ia, ib)] = J_of(xc)
    with torch.no_grad():
        ov, _ = ProductOverlapLoss(ma, mb, OVERLAP_SIZE, OVERLAP_SIZE).to(DTYPE)(
            J_dict, xa0, xb0, pa[1], pb[1], hard_cutoff=True)
    return sum(nijs) / 4.0, ov.item()


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def dense_eval(geo, task, J_fn, n, chunk=500):
    accs, nijs = [], []
    for s in range(0, n, chunk):
        x = geo.sample(min(chunk, n - s))
        with torch.no_grad():
            accs.append(acs_residual(geo, J_fn, x))
        nijs.append(nij_norm(geo, J_fn, x).detach())
    acs, nij = torch.cat(accs), torch.cat(nijs)
    return {'acs_mean': acs.mean().item(), 'acs_sup': acs.max().item(),
            'nij_mean': nij.mean().item(), 'nij_sup': nij.max().item()}


def full_report(geo, task, J_fn, meta, n_dense, n_chart):
    rep = {'dense': dense_eval(geo, task, J_fn, n_dense)}
    field_for_chart = J_fn
    if task == 'acs':
        if rep['dense']['acs_sup'] < 1.0:
            J_star = polar_corrected(geo, J_fn)
            rep['corrected'] = dense_eval(geo, task, J_star, n_dense)
            rep['certified_exact_acs'] = rep['corrected']['acs_sup'] < 1e-8
            field_for_chart = J_star
        else:
            rep['certified_exact_acs'] = False
    cs, acs = chart_protocol(geo, field_for_chart, n_chart)
    rep['chart_cs'], rep['chart_acs'] = cs, acs
    return rep


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', required=True, choices=ALL_CASES)
    ap.add_argument('--task', required=True, choices=['acs', 'cs'])
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--out', required=True)
    ap.add_argument('--steps', type=int, default=None)
    ap.add_argument('--hidden', default=None)
    ap.add_argument('--stages', type=int, default=None)
    ap.add_argument('--rar', action='store_true')
    ap.add_argument('--lbfgs', type=int, default=200)
    ap.add_argument('--wall-minutes', type=float, default=4200.0)
    ap.add_argument('--eval-every', type=int, default=2000)
    ap.add_argument('--ckpt-every', type=int, default=1000)
    ap.add_argument('--n-dense', type=int, default=6000)
    ap.add_argument('--n-chart', type=int, default=3000)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    torch.manual_seed(args.seed)
    geo = Geometry(args.case)
    d, N = geo.d, geo.N
    steps = args.steps or (10000 if args.task == 'acs' else 40000)
    hidden = [int(h) for h in args.hidden.split(',')] if args.hidden else \
        ([256, 256, 256, 256] if d >= 6 else [128, 128, 128])
    stages = args.stages or (2 if args.task == 'cs' else 1)
    n_mat = d if geo.kind == 'frame' else N
    out_scale = 0.3 if args.task == 'acs' else 0.5
    model = Net(N, n_mat, stages, hidden, out_scale).to(DTYPE)
    J_fn, meta = make_field(geo, args.task, model)
    tag = f"{args.case}/{args.task}/s{args.seed} [{geo.kind}]"

    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=steps, eta_min=1e-5)
    state_path = os.path.join(args.out, 'state.pt')
    start_step, elapsed_prev, init_info = 0, 0.0, None
    if os.path.exists(state_path):
        st = torch.load(state_path, map_location='cpu', weights_only=False)
        model.load_state_dict(st['model']); opt.load_state_dict(st['opt'])
        sched.load_state_dict(st['sched'])
        start_step, elapsed_prev, init_info = st['step'], st['elapsed'], st.get('init')
        print(f"[{tag}] resumed at step {start_step} ({elapsed_prev/60:.1f} min used)", flush=True)
    if init_info is None:
        x0 = geo.sample(300)
        _, init_info = unit_loss(geo, args.task, J_fn, x0)
        print(f"[{tag}] init: {init_info}", flush=True)

    progress = open(os.path.join(args.out, 'progress.jsonl'), 'a')
    t0 = time.time()
    elapsed = lambda: elapsed_prev + (time.time() - t0)
    wall_s = args.wall_minutes * 60.0
    nonfinite, step = 0, start_step
    stopped_by_wall = False

    def save_state(step_):
        torch.save({'model': model.state_dict(), 'opt': opt.state_dict(),
                    'sched': sched.state_dict(), 'step': step_,
                    'elapsed': elapsed(), 'init': init_info}, state_path)

    for step in range(start_step, steps):
        if elapsed() > 0.85 * wall_s:
            stopped_by_wall = True
            break
        if args.rar and args.task == 'cs':
            xc = geo.sample(2 * BATCH)
            r = nij_norm(geo, J_fn, xc).detach()
            idx = torch.topk(r, BATCH // 2).indices
            x = torch.cat([xc[idx], geo.sample(BATCH - BATCH // 2)], dim=0)
        else:
            x = geo.sample(BATCH)
        loss, info = unit_loss(geo, args.task, J_fn, x)
        if not bool(torch.isfinite(loss)):
            nonfinite += 1
            if nonfinite >= 200:
                print(f"[{tag}] aborting: 200 consecutive non-finite losses", flush=True)
                break
            continue
        nonfinite = 0
        opt.zero_grad(); loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step(); sched.step()
        if (step + 1) % args.ckpt_every == 0:
            save_state(step + 1)
        if (step + 1) % args.eval_every == 0 or step == steps - 1:
            ev = dense_eval(geo, args.task, J_fn, 1500)
            rec = {'step': step + 1, 'elapsed_min': elapsed() / 60.0, 'batch': info, **ev}
            progress.write(json.dumps(rec) + '\n'); progress.flush()
            print(f"[{tag}] step {step+1:6d} | batch {info} | held-out {ev}", flush=True)
    save_state(step + 1 if not stopped_by_wall else step)

    # L-BFGS polish (time-guarded)
    if elapsed() < 0.92 * wall_s and args.lbfgs > 0:
        xf = geo.sample(768)
        lb = torch.optim.LBFGS(model.parameters(), lr=1.0, max_iter=args.lbfgs,
                               history_size=100, line_search_fn='strong_wolfe')

        def closure():
            lb.zero_grad()
            l, _ = unit_loss(geo, args.task, J_fn, xf)
            l.backward()
            return l
        try:
            lb.step(closure)
        except Exception as e:                     # keep the Adam result on failure
            print(f"[{tag}] LBFGS failed ({type(e).__name__}); keeping Adam weights", flush=True)
        save_state(steps)

    torch.save({'model': model.state_dict(), 'case': args.case, 'task': args.task,
                'seed': args.seed, 'hidden': hidden, 'stages': stages, 'kind': geo.kind},
               os.path.join(args.out, 'final.pt'))

    rep = full_report(geo, args.task, J_fn, meta, args.n_dense, args.n_chart)
    result = {'case': args.case, 'task': args.task, 'seed': args.seed,
              'formulation': geo.kind, 'hidden': hidden, 'stages': stages,
              'steps_done': min(step + 1, steps), 'steps_planned': steps,
              'stopped_by_wall': stopped_by_wall, 'elapsed_min': elapsed() / 60.0,
              'init': init_info, 'exact_acs_by_construction': meta['exact_acs'],
              'table_before': TABLE_BEFORE[args.case][args.task], **rep}
    with open(os.path.join(args.out, 'result.json'), 'w') as f:
        json.dump(result, f, indent=2)
    print(f"[{tag}] RESULT: chart_{args.task} = "
          f"{rep['chart_cs'] if args.task == 'cs' else rep['chart_acs']:.4e} "
          f"(table before {TABLE_BEFORE[args.case][args.task]}) | dense {rep['dense']}", flush=True)


if __name__ == '__main__':
    main()
