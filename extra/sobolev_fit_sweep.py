"""Sobolev-fit architecture sweep (Stage A of the derivative-error program).

For each spectral architecture variant, fit one PatchSubNet to the exact
Calabi-Eckmann structure on S^5 x S^1 with JOINT value + derivative
supervision:

    L = mean ||J_net - J*||_F^2  +  LAMBDA_D * mean ||dJ_net - dJ*||_F^2

then evaluate the per-patch Nijenhuis loss under the paper protocol
(hard cutoff). This measures each architecture's ability to represent the
solution in the first-order Sobolev sense — which is what the Nijenhuis
functional actually requires (an 8%-value-accurate L2 fit was measured to
carry 51% derivative error and Nij ~27).

The best variant is then Sobolev-fitted on all four patches and evaluated
under the full ProductTotalLoss paper protocol — that number is the TRUE
ARCHITECTURE CEILING for from-scratch training to be measured against.

No-bias note: the analytic structure is used offline for architecture
SELECTION only (tune where ground truth exists, deploy the identical recipe
on all manifolds); the from-scratch training loop never sees it.

Run:  python extra/sobolev_fit_sweep.py            (~2-3 h on a laptop CPU)
"""
from __future__ import annotations

import os
os.environ.setdefault('PYTORCH_ENABLE_MPS_FALLBACK', '1')

import sys
import time
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.geometry.sphere import EvenSphere
from src.geometry.sphere_product import SphereProduct
from src.geometry.hopf_complex_structure_fixed import hopf_complex_structure
from src.models import FourPatchProductNet
from src.loss.sphere_loss import _compute_A_and_dA, _compute_dJ_from_dA, _nij_tensor
from src.loss.product_loss import ProductTotalLoss

torch.manual_seed(42)

DIM_ALPHA, DIM_BETA = 5, 1
D = DIM_ALPHA + DIM_BETA
R_FIT   = 0.85
R_EVAL  = 0.6142135  # RADIAL_MIDPOINT + 0.2 (hard-cutoff plateau per factor)
LAMBDA_D = 0.05      # balances value (~38) and derivative (~1156) terms at init

FIT_STEPS  = 600
BATCH      = 192
DATA_N     = 8_000
EVAL_N     = 2_000
LR         = 1e-3

HEAVY = [512, 1024, 1024, 1024, 1024, 512]
LIGHT = [256, 512, 512, 256]

# (name, hidden, activation, use_fourier, n_features, fourier_scale)
VARIANTS = [
    ('gelu_f128_s1.0 (baseline)', HEAVY, 'gelu', True, 128, 1.0),
    ('gelu_f128_s0.5',            HEAVY, 'gelu', True, 128, 0.5),
    ('gelu_f128_s0.25',           HEAVY, 'gelu', True, 128, 0.25),
    ('gelu_f32_s1.0',             HEAVY, 'gelu', True,  32, 1.0),
    ('gelu_f32_s0.5',             HEAVY, 'gelu', True,  32, 0.5),
    ('gelu_nofourier',            HEAVY, 'gelu', False, 0, 1.0),
    ('sin_f128_s1.0',             HEAVY, 'sin',  True, 128, 1.0),
    ('sin_nofourier',             HEAVY, 'sin',  False, 0, 1.0),
    ('light_gelu_nofourier',      LIGHT, 'gelu', False, 0, 1.0),
]


def analytic_targets(ia, ib, n, seed):
    """Precompute (x, J*, dJ*) on patch (ia, ib), filtered to R_FIT, float32."""
    g = torch.Generator().manual_seed(seed)
    sp = SphereProduct(DIM_ALPHA, DIM_BETA)
    xs, Js, dJs = [], [], []
    have = 0
    while have < n:
        torch.manual_seed(int(torch.randint(0, 2**31 - 1, (1,), generator=g)))
        xa, xb = sp.sample_points(4000)
        xa_p = xa if ia == 0 else EvenSphere.patch_change_coords(xa)
        xb_p = xb if ib == 0 else EvenSphere.patch_change_coords(xb)
        keep = (xa_p.norm(dim=1) <= R_FIT) & (xb_p.norm(dim=1) <= R_FIT)
        xa_p, xb_p = xa_p[keep].double(), xb_p[keep].double()
        x = torch.cat([xa_p, xb_p], dim=1).requires_grad_(True)
        J = hopf_complex_structure(x[:, :DIM_ALPHA], x[:, DIM_ALPHA:], ia, ib)
        b = len(x)
        dJ = torch.zeros(b, D, D, D, dtype=torch.float64)
        for i in range(D):
            for j in range(D):
                gr, = torch.autograd.grad(J[:, i, j].sum(), x, retain_graph=True)
                dJ[:, i, j, :] = gr
        xs.append(x.detach().float()); Js.append(J.detach().float()); dJs.append(dJ.float())
        have += b
    return (torch.cat(xs)[:n], torch.cat(Js)[:n], torch.cat(dJs)[:n])


def sobolev_fit(subnet, x_d, J_t, dJ_t, steps=FIT_STEPS, batch=BATCH, log_every=150):
    """Joint value+derivative supervised fit of one subnet. Returns loss parts."""
    opt = torch.optim.Adam(subnet.parameters(), lr=LR)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=steps, eta_min=1e-5)
    J0 = subnet.J0
    for step in range(steps):
        idx = torch.randint(0, len(x_d), (batch,))
        xb_, Jt_, dJt_ = x_d[idx], J_t[idx], dJ_t[idx]
        opt.zero_grad()
        A, dA = _compute_A_and_dA(subnet, xb_)
        J  = torch.matrix_exp(A) @ J0 @ torch.matrix_exp(-A)
        dJ = _compute_dJ_from_dA(A, dA, J0)
        l_val = ((J - Jt_) ** 2).sum(dim=(1, 2)).mean()
        l_der = ((dJ - dJt_) ** 2).sum(dim=(1, 2, 3)).mean()
        loss = l_val + LAMBDA_D * l_der
        loss.backward()
        torch.nn.utils.clip_grad_norm_(subnet.parameters(), 1.0)
        opt.step(); sched.step()
        if step % log_every == 0 or step == steps - 1:
            print(f"    step {step:4d} | val {l_val.item():.3e} | der {l_der.item():.3e}", flush=True)
    return l_val.item(), l_der.item()


def eval_patch(subnet, x_eval, J_t, dJ_t):
    """Paper-protocol per-patch Nij (hard cutoff) + relative errors."""
    mask = (x_eval[:, :DIM_ALPHA].norm(dim=1) <= R_EVAL) & \
           (x_eval[:, DIM_ALPHA:].norm(dim=1) <= R_EVAL)
    xe, Jt_, dJt_ = x_eval[mask], J_t[mask], dJ_t[mask]
    with torch.no_grad():
        pass  # A/dA path needs grad machinery internally; run under enable_grad
    A, dA = _compute_A_and_dA(subnet, xe)
    J0 = subnet.J0
    J  = (torch.matrix_exp(A) @ J0 @ torch.matrix_exp(-A)).detach()
    dJ = _compute_dJ_from_dA(A, dA, J0).detach()
    N = _nij_tensor(J, dJ)
    nij = torch.sqrt((N ** 2).sum(dim=(1, 2, 3)) + 1e-20).mean().item()
    rel_val = ((J - Jt_).norm(dim=(1, 2)) / Jt_.norm(dim=(1, 2))).mean().item()
    rel_der = ((dJ - dJt_).reshape(len(xe), -1).norm(dim=1)
               / dJt_.reshape(len(xe), -1).norm(dim=1).clamp(min=1e-12)).mean().item()
    return nij, rel_val, rel_der, int(mask.sum())


def make_model(hidden, act, use_f, n_f, scale):
    return FourPatchProductNet(
        dim_alpha=DIM_ALPHA, dim_beta=DIM_BETA, hidden_layers=hidden,
        activation=act, use_fourier_features=use_f,
        num_fourier_features=max(n_f, 1), fourier_scale=scale,
        initialization='xavier', matrix_exp_terms=None, use_residual_blocks=True,
    ).float()


def main():
    print("Precomputing analytic targets for patch (0,0) [fit + eval sets] ...", flush=True)
    t0 = time.time()
    x_fit, J_fit, dJ_fit = analytic_targets(0, 0, DATA_N, seed=1)
    x_ev,  J_ev,  dJ_ev  = analytic_targets(0, 0, EVAL_N, seed=2)
    print(f"  done in {time.time()-t0:.0f}s", flush=True)

    results = []
    for (name, hidden, act, use_f, n_f, scale) in VARIANTS:
        print(f"\n### Variant: {name}", flush=True)
        torch.manual_seed(7)
        model = make_model(hidden, act, use_f, n_f, scale)
        subnet = model.patches['patch_00']
        n_par = sum(p.numel() for p in subnet.parameters())
        t0 = time.time()
        l_val, l_der = sobolev_fit(subnet, x_fit, J_fit, dJ_fit)
        nij, rel_val, rel_der, n_act = eval_patch(subnet, x_ev, J_ev, dJ_ev)
        dt = time.time() - t0
        print(f"  -> Nij(0,0) {nij:.3e} | rel val {rel_val:.1%} | rel der {rel_der:.1%} "
              f"| {n_par/1e6:.1f}M par | {dt:.0f}s", flush=True)
        results.append((name, nij, rel_val, rel_der, n_par, dt))
        del model, subnet

    results.sort(key=lambda r: r[1])
    print(f"\n{'='*88}")
    print("SOBOLEV-FIT SWEEP — ranked by paper-protocol patch-(0,0) Nijenhuis")
    print(f"{'='*88}")
    print(f"{'variant':<28} {'Nij(0,0)':>10} {'rel val':>8} {'rel der':>8} {'params':>8} {'time':>6}")
    for name, nij, rv, rd, npar, dt in results:
        print(f"{name:<28} {nij:>10.3e} {rv:>7.1%} {rd:>7.1%} {npar/1e6:>7.1f}M {dt:>5.0f}s")

    # ---- Full 4-patch ceiling for the winner ----
    winner = results[0][0]
    cfg = next(v for v in VARIANTS if v[0] == winner)
    print(f"\n### CEILING RUN: 4-patch Sobolev fit of winner '{winner}'", flush=True)
    torch.manual_seed(7)
    model = make_model(*cfg[1:])
    for (ia, ib) in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        print(f"  patch ({ia},{ib}):", flush=True)
        xd, Jd, dJd = analytic_targets(ia, ib, DATA_N, seed=10 + ia * 2 + ib)
        sobolev_fit(model.patches[f'patch_{ia}{ib}'], xd, Jd, dJd, log_every=300)
        del xd, Jd, dJd

    print("\n  Full paper-protocol evaluation (chunked):", flush=True)
    sp = SphereProduct(DIM_ALPHA, DIM_BETA)
    loss_fn = ProductTotalLoss(
        dim_alpha=DIM_ALPHA, dim_beta=DIM_BETA,
        weights={'nijenhuis': 1.0, 'overlap': 100.0,
                 'overlap_size_alpha': 0.2, 'overlap_size_beta': 0.2,
                 'use_round_metric_contraction': False},
    ).float()
    acc, n_chunks = {}, 4
    for _ in range(n_chunks):
        xa, xb = sp.sample_points(1500)
        xa, xb = xa.float(), xb.float()
        with torch.no_grad():
            J_dict, xa, xb, xa2, xb2 = model(xa, xb)
        _, info = loss_fn(J_dict, xa, xb, xa2, xb2, model,
                          hard_cutoff=True, report_all=True)
        for k, v in info.items():
            if isinstance(v, float):
                acc[k] = acc.get(k, 0.0) + v / n_chunks
    nijs = [acc[f'nij_{ia}{ib}'] for ia, ib in [(0,0),(0,1),(1,0),(1,1)]]
    print(f"\n{'='*88}")
    print(f"ARCHITECTURE CEILING (winner: {winner})")
    print(f"{'='*88}")
    for (ia, ib), v in zip([(0,0),(0,1),(1,0),(1,1)], nijs):
        print(f"  Nij ({ia},{ib}):   {v:.4e}")
    print(f"  MEAN NIJ:      {sum(nijs)/4:.4e}")
    print(f"  OVERLAP:       {acc['overlap_norm']:.4e}")
    print(f"\n  (from-scratch table: CS 6.42, ACS 7.98e-2; certified floor ~1e-10)")


if __name__ == '__main__':
    main()
