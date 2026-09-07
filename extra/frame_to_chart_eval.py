"""Frame → chart conversion + paper-protocol evaluation.

Takes an integrable structure found in the global-frame formulation
(see extra/frame_constant_solve.py) — here the frame-constant solutions on
S³×S³, S³×S¹, S⁵×S¹ — converts it to the codebase's 4-patch ball-chart
representation, and evaluates it under the UNCHANGED paper protocol
(hard-cutoff Nijenhuis per patch + cross-patch overlap, exactly as the
results table). This is the number the table reports going to ~0.

Conversion at a physical point p with patch chart coords x:
    chart basis  = columns of Dσ(x)  (σ = per-factor ball_to_sphere with the
                   direction-reversal-consistent patch conventions),
    V(x)         = frame vectors in chart components  (solve Dσ · V = Eᵀ),
    J_chart(x)   = V · J_frame(p) · V⁻¹.
Everything is differentiable torch, so ∂J_chart for the Nijenhuis evaluation
comes from autograd through the whole chain.

Run:  python extra/frame_to_chart_eval.py
"""
from __future__ import annotations

import os
import sys
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.geometry.sphere import EvenSphere, RADIAL_MIDPOINT
from src.geometry.sphere_product import SphereProduct
from src.geometry.frames import FRAMES
from src.loss.frame_loss import frame_nij_norm
from src.loss.sphere_loss import _nij_tensor
from src.loss.product_loss import ProductNijenhuisLoss, ProductOverlapLoss

torch.manual_seed(42)
DTYPE = torch.float64
OVERLAP_SIZE = 0.2                       # paper-protocol value (product configs)
R_CUT = RADIAL_MIDPOINT + OVERLAP_SIZE
N_TEST = 6000
CHUNK = 1500

CASES = {
    's3xs3': dict(dim_a=3, dim_b=3, table_cs=7.63),
    's3xs1': dict(dim_a=3, dim_b=1, table_cs=2.00e-1),
    's5xs1': dict(dim_a=5, dim_b=1, table_cs=6.42),
}


def factor_embedding(x: torch.Tensor, chart: int) -> torch.Tensor:
    """Per-factor chart embedding consistent with the gluing φ: chart 0 is
    ball_to_sphere(·, from_north=True); chart 1 flips the equatorial ambient
    coordinates of ball_to_sphere(·, from_north=False) so that
    σ₁(φ(x)) = σ₀(x) (φ reverses the tangential direction)."""
    if chart == 0:
        return EvenSphere.ball_to_sphere(x, from_north=True)
    p = EvenSphere.ball_to_sphere(x, from_north=False)
    return torch.cat([-p[:, :-1], p[:, -1:]], dim=1)


def product_embedding(xa, xb, ia, ib):
    return torch.cat([factor_embedding(xa, ia), factor_embedding(xb, ib)], dim=1)


def chart_J_fn(frame, J_frame_const, dim_a, ia, ib):
    """Returns f(x_chart (B, d)) -> J_chart (B, d, d), fully differentiable."""
    d = frame.dim

    def emb_single(xc):
        xa = xc[:dim_a].unsqueeze(0)
        xb = xc[dim_a:].unsqueeze(0)
        return product_embedding(xa, xb, ia, ib).squeeze(0)

    def J_of(x_chart):
        B = x_chart.shape[0]
        p = torch.func.vmap(emb_single)(x_chart)                # (B, N)
        D = torch.func.vmap(torch.func.jacfwd(emb_single))(x_chart)   # (B, N, d)
        E = frame.frame_fields(p)                               # (B, d, N)
        V = torch.linalg.lstsq(D, E.transpose(1, 2)).solution   # (B, d, d)
        Jf = J_frame_const.unsqueeze(0).expand(B, d, d)
        return V @ Jf @ torch.linalg.inv(V)

    return J_of


def dJ_autograd(J_of, x, d):
    xr = x.detach().clone().requires_grad_(True)
    J = J_of(xr)
    dJ = torch.zeros(len(x), d, d, d, dtype=x.dtype)
    for i in range(d):
        for j in range(d):
            g, = torch.autograd.grad(J[:, i, j].sum(), xr, retain_graph=True)
            dJ[:, i, j, :] = g
    return J.detach(), dJ


def solve_constant(name):
    """Reproduce the frame-constant solution (canonical J₀ is integrable)."""
    frame = FRAMES[name]()
    d = frame.dim
    J0 = torch.zeros(d, d, dtype=DTYPE)
    for k in range(d // 2):
        J0[2*k, 2*k+1] = -1.0
        J0[2*k+1, 2*k] = 1.0
    # verify integrability in-frame at random points
    x = frame.sample(500)
    c = frame.structure_functions(x)
    nij = frame_nij_norm(J0.unsqueeze(0).expand(len(x), d, d),
                         torch.zeros(len(x), d, d, d, dtype=DTYPE), c)
    return J0, nij.mean().item()


def evaluate_case(name):
    cfg = CASES[name]
    dim_a, dim_b = cfg['dim_a'], cfg['dim_b']
    d = dim_a + dim_b
    frame = FRAMES[name]()
    J_frame, in_frame_nij = solve_constant(name)

    print("\n" + "=" * 78)
    print(f"{name.upper()}  (S^{dim_a} x S^{dim_b})   —  frame-constant solution → chart protocol")
    print("=" * 78)
    print(f"  in-frame ||N|| (500 pts):      {in_frame_nij:.3e}")

    sp = SphereProduct(dim_a, dim_b)
    xa0, xb0 = sp.sample_points(N_TEST)
    xa0, xb0 = xa0.to(DTYPE), xb0.to(DTYPE)

    # sanity: same-physical-point embedding across all four patches
    patch_a = {0: xa0, 1: EvenSphere.patch_change_coords(xa0)}
    patch_b = {0: xb0, 1: EvenSphere.patch_change_coords(xb0)}
    p_ref = product_embedding(xa0, xb0, 0, 0)
    emb_err = max((product_embedding(patch_a[ia], patch_b[ib], ia, ib) - p_ref)
                  .abs().max().item()
                  for ia, ib in [(0, 1), (1, 0), (1, 1)])
    print(f"  cross-patch embedding agreement: {emb_err:.2e}")

    # ---- Nijenhuis, paper protocol (hard cutoff, per patch) ----
    nij_vals = {}
    for (ia, ib) in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        xa_p, xb_p = patch_a[ia], patch_b[ib]
        mask = (xa_p.norm(dim=1) <= R_CUT) & (xb_p.norm(dim=1) <= R_CUT)
        x_act = torch.cat([xa_p[mask], xb_p[mask]], dim=1)
        J_of = chart_J_fn(frame, J_frame, dim_a, ia, ib)
        norms = []
        for s in range(0, len(x_act), CHUNK):
            J, dJ = dJ_autograd(J_of, x_act[s:s+CHUNK], d)
            N = _nij_tensor(J, dJ)
            norms.append(torch.sqrt((N**2).sum(dim=(1, 2, 3)) + 1e-20))
        v = torch.cat(norms)
        nij_vals[(ia, ib)] = v.mean().item()
        print(f"  Nij patch ({ia},{ib}):  {v.mean():.6e}   (max {v.max():.2e}, {len(x_act)} pts)")
    mean_nij = sum(nij_vals.values()) / 4.0

    # ---- Overlap, paper protocol ----
    ov = ProductOverlapLoss(dim_a, dim_b, OVERLAP_SIZE, OVERLAP_SIZE).to(DTYPE)
    J_dict = {}
    for (ia, ib) in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        with torch.no_grad():
            J_dict[(ia, ib)] = chart_J_fn(frame, J_frame, dim_a, ia, ib)(
                torch.cat([patch_a[ia], patch_b[ib]], dim=1))
    with torch.no_grad():
        ov_loss, _ = ov(J_dict, xa0, xb0, patch_a[1], patch_b[1], hard_cutoff=True)

    print(f"\n  MEAN NIJ (paper protocol):     {mean_nij:.6e}   [table: {cfg['table_cs']}]")
    print(f"  OVERLAP  (paper protocol):     {ov_loss.item():.6e}")
    print(f"  TOTAL (w_nij=1, w_ov=100):     {mean_nij + 100*ov_loss.item():.6e}")
    return mean_nij, ov_loss.item()


if __name__ == '__main__':
    results = {}
    for name in CASES:
        results[name] = evaluate_case(name)
    print("\n" + "=" * 78)
    print("SUMMARY — numerically built complex structures, paper-protocol losses")
    print("=" * 78)
    print(f"{'case':<10} {'CS (mean Nij)':>15} {'ACS (overlap)':>15} {'table CS before':>16}")
    for name, (nij, ovl) in results.items():
        print(f"{name:<10} {nij:>15.6e} {ovl:>15.6e} {CASES[name]['table_cs']:>16}")
