"""Decisive tests of the global-frame (chart-free) Nijenhuis formulation.

1. Frame self-tests: structure functions of the S³×S³ Lie frame are constant
   and antisymmetric; the c = 0 kernel reduces exactly to the coordinate
   kernel `_nij_tensor`.
2. S⁵×S¹: the constant matrix J = i (canonical 2×2 blocks in ℂ³ ≅ ℝ⁶) should
   be EXACTLY integrable in the Hopf-trivialization frame — frame-Nij ≈ 0.
   This is the frame-formulation floor calibration.
3. S³×S³ and S³×S¹ (the failing CS cases): optimize a SINGLE constant matrix
   A ∈ ℝ^{d×d} (36 / 16 parameters) on the purely algebraic frame-Nijenhuis
   loss. By Samelson's theorem frame-constant complex structures exist on
   these Lie groups — if the optimization reaches ~0, we have NUMERICALLY
   BUILT integrable complex structures on the failing cases, in seconds,
   with the PINN loss (restricted to constants) and no reference solution.

Run:  python extra/frame_constant_solve.py
"""
from __future__ import annotations

import os
import sys
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.geometry.frames import FRAMES
from src.loss.frame_loss import nijenhuis_frame_tensor, frame_nij_norm
from src.loss.sphere_loss import _nij_tensor

torch.manual_seed(42)
DTYPE = torch.float64


def build_J0(d: int) -> torch.Tensor:
    J0 = torch.zeros(d, d, dtype=DTYPE)
    for k in range(d // 2):
        J0[2 * k, 2 * k + 1] = -1.0
        J0[2 * k + 1, 2 * k] = 1.0
    return J0


def test_frame_selfchecks():
    print("=" * 78)
    print("1. Frame self-tests")
    print("=" * 78)

    # c = 0 reduction == coordinate kernel on random tensors
    B, d = 32, 6
    J = torch.randn(B, d, d, dtype=DTYPE)
    dJ = torch.randn(B, d, d, d, dtype=DTYPE)
    c0 = torch.zeros(B, d, d, d, dtype=DTYPE)
    diff = (nijenhuis_frame_tensor(J, dJ, c0) - _nij_tensor(J, dJ)).abs().max()
    print(f"  c=0 kernel vs coordinate _nij_tensor: max diff = {diff:.2e}")
    assert diff < 1e-12

    for name in ('s3xs3', 's3xs1', 's5xs1'):
        fr = FRAMES[name]()
        x = fr.sample(500)
        E = fr.frame_fields(x)
        # frame Gram: should be well-conditioned (frame spans the tangent space)
        G = E @ E.transpose(1, 2)
        conds = torch.linalg.cond(G)
        c = fr.structure_functions(x)
        anti = (c + c.transpose(2, 3)).abs().max()
        cvar = (c - c.mean(dim=0, keepdim=True)).abs().max()
        print(f"  {name}: Gram cond max {conds.max():.3f} | c antisym err {anti:.2e} "
              f"| c position-variation {cvar:.2e}"
              + ("  (constant ✓)" if cvar < 1e-10 else "  (position-dependent)"))


def test_s5s1_hopf_constant():
    print("\n" + "=" * 78)
    print("2. S⁵×S¹: constant J = i in the Hopf-trivialization frame")
    print("=" * 78)
    fr = FRAMES['s5xs1']()
    x = fr.sample(2000)
    c = fr.structure_functions(x)
    J = build_J0(6).unsqueeze(0).expand(len(x), 6, 6)
    dJ = torch.zeros(len(x), 6, 6, 6, dtype=DTYPE)
    nij = frame_nij_norm(J, dJ, c)
    print(f"  frame ||N|| at constant i-matrix: mean {nij.mean():.3e}, "
          f"max {nij.max():.3e}")
    print(f"  → ~0 means the frame formulation certifies the Hopf CS on S⁵×S¹.")
    return nij.mean().item()


def constant_solve(name: str, steps: int = 400):
    print("\n" + "=" * 78)
    print(f"3. {name}: constant-structure solve (algebraic frame-Nijenhuis)")
    print("=" * 78)
    fr = FRAMES[name]()
    d = fr.dim
    x = fr.sample(200)
    c = fr.structure_functions(x)
    cvar = (c - c.mean(dim=0, keepdim=True)).abs().max()
    c_use = c if cvar > 1e-10 else c[:1]           # constant c: one sample suffices
    B = c_use.shape[0]
    J0 = build_J0(d)

    A = torch.zeros(d, d, dtype=DTYPE, requires_grad=True)
    dJ = torch.zeros(B, d, d, d, dtype=DTYPE)

    def loss_of(A_):
        J = torch.matrix_exp(A_) @ J0 @ torch.matrix_exp(-A_)
        Jb = J.unsqueeze(0).expand(B, d, d)
        N = nijenhuis_frame_tensor(Jb, dJ, c_use)
        return (N ** 2).sum(dim=(1, 2, 3)).mean()

    with torch.no_grad():
        print(f"  loss at A = 0 (J = J₀):        {loss_of(A).item():.6e}")

    # small random init to break symmetry, then LBFGS
    with torch.no_grad():
        A += 0.1 * torch.randn(d, d, dtype=DTYPE)
    opt = torch.optim.LBFGS([A], lr=1.0, max_iter=steps, history_size=100,
                            tolerance_grad=1e-16, tolerance_change=1e-16,
                            line_search_fn='strong_wolfe')

    def closure():
        opt.zero_grad()
        l = loss_of(A)
        l.backward()
        return l

    opt.step(closure)
    final = loss_of(A).item()
    with torch.no_grad():
        J = torch.matrix_exp(A) @ J0 @ torch.matrix_exp(-A)
        Jb = J.unsqueeze(0).expand(B, d, d)
        nij = frame_nij_norm(Jb, dJ, c_use)
        j2err = (J @ J + torch.eye(d, dtype=DTYPE)).abs().max()
    print(f"  after LBFGS ({steps} iters):     loss {final:.6e}")
    print(f"  frame ||N|| at solution:        mean {nij.mean():.3e}")
    print(f"  J² + I:                         {j2err:.2e}")
    if nij.mean().item() < 1e-6:
        print(f"  *** INTEGRABLE STRUCTURE FOUND on {name} "
              f"(constant in the global frame) ***")
    return final, A.detach(), J.detach()


if __name__ == '__main__':
    test_frame_selfchecks()
    test_s5s1_hopf_constant()
    for name in ('s3xs3', 's3xs1', 's5xs1'):
        constant_solve(name)
