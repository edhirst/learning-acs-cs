"""
The Hopf complex structure on S^{2n−1} × S¹ in the AMBIENT-PROJECTOR
representation (src/geometry/ambient.py): a closed-form tangential matrix
field on ℝ^{2n} × ℝ², identical — as a structure on the manifold — to the
chart formula of `hopf_complex_structure_fixed.py` (verified by pulling it
back through the ball charts; see extra/ambient_supervised_hopf.py --verify).

Construction
------------
Point p = (q, w), q ∈ S^{2n−1} ⊂ ℝ^{2n}, w ∈ S¹ ⊂ ℝ², angle θ with
∂_θ = R(w) := (−w₁, w₀). The Hopf manifold (ℂⁿ∖0)/(z ∼ λz), λ = e^{2πc}, is
identified with S^{2n−1} × S¹ by Φ(q, θ) = e^{cθ} q. With I the complex
structure of ℂⁿ ≅ ℝ^{2n} and the linear isomorphism

    Λ : T_p M → ℝ^{2n},   Λ(v, s·R) = v + c·s·q        (v ⊥ q),

the Hopf structure is J = Λ⁻¹ I Λ (the factor e^{cθ} commutes with I and
cancels). With m := I q (unit, ⊥ q), acting on ambient u = (a, b):

    J(a, b) = ( I a − ⟨a,q⟩ m + c⟨b,R⟩ m + ⟨a,m⟩ q ,  −(1/c)⟨a,m⟩ R ),

i.e. the block matrix

    J_qq = I − m qᵀ + q mᵀ,   J_qw = c · m Rᵀ,   J_wq = −(1/c) · R mᵀ,   J_ww = 0.

It is tangential (J P = P J = J on M), J² = −P, and N_J = 0. Every block is a
polynomial of degree ≤ 2 in (q, w): in the ambient picture the Hopf
structure is a LOW-DEGREE field (contrast the chart picture, where the same
structure carries the chart Jacobians).

Conventions matching the rest of the codebase
---------------------------------------------
* `hopf_complex_structure(_fixed)` uses z ∼ e·z, i.e. t = θ/(2π) ∈ ℝ/ℤ, so
  c = 1/(2π) (DEFAULT here). c = 1 gives the round-ORTHOGONAL Hopf structure
  (Λ is then an isometry) — a different point of the moduli of Hopf
  manifolds, homotopic to the c = 1/(2π) one through J_c.
* Ambient ordering of the codebase embedding (`EvenSphere.ball_to_sphere`):
  q = (x_1, …, x_{2n−1}, z) with the pole coordinate z LAST, while
  `hopf_map_F` orders F ∝ (z, x_1, …, x_{2n−1}) with J₀ = blockdiag(
  [[0,−1],[1,0]], …) on F. Hence I = Lᵀ J₀ L with L the permutation
  q ↦ (q_{2n−1}, q_0, …, q_{2n−2}).
* The S¹ angle of `hopf_map_F` is θ = atan2(w₀, −w₁) with w = (s_xy, s_z)
  = ball_to_sphere(x_β); then ∂_θ w = (−w₁, w₀) = R(w). (Only ∂_θ enters J:
  J is invariant under rotations of the circle.)
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import math
import torch

C_MAIN = 1.0 / (2.0 * math.pi)      # z ∼ e·z : the structure used on main


def hopf_I(n: int, dtype=torch.float64, device=None) -> torch.Tensor:
    """Complex structure I on ℝ^{2n} in the codebase's ambient ordering,
    equal to hopf_map_F's J₀ transported by the coordinate permutation."""
    d = 2 * n
    J0 = torch.zeros(d, d, dtype=dtype, device=device)
    for k in range(n):
        J0[2 * k, 2 * k + 1] = -1.0
        J0[2 * k + 1, 2 * k] = 1.0
    L = torch.zeros(d, d, dtype=dtype, device=device)
    L[0, d - 1] = 1.0                       # F_0 = q_{2n−1}  (pole coordinate)
    for k in range(1, d):
        L[k, k - 1] = 1.0                   # F_k = q_{k−1}
    return L.T @ J0 @ L


def hopf_ambient_acs(p: torch.Tensor, n: int, c: float = C_MAIN,
                     project: bool = True) -> torch.Tensor:
    """Hopf complex structure at ambient points p = (q, w) ∈ S^{2n−1}×S¹.

    Args:
        p: (B, 2n+2) ambient points (q ∈ ℝ^{2n}, w ∈ ℝ²).
        n: complex dimension (S³×S¹: n = 2, S⁵×S¹: n = 3).
        c: log-period constant, λ = e^{2πc}; default 1/(2π) (main's z ∼ e·z).
        project: wrap as P·J·P with the smooth projector (exactly tangential
            also off M; identical on M). Tangential derivatives — the only
            ones that enter N_J — do not depend on this choice.

    Returns:
        (B, 2n+2, 2n+2) matrix field.
    """
    d = 2 * n
    q, w = p[:, :d], p[:, d:]
    I = hopf_I(n, dtype=p.dtype, device=p.device)
    m = q @ I.T                                              # I q
    R = torch.stack([-w[:, 1], w[:, 0]], dim=1)              # ∂_θ
    Jqq = I.unsqueeze(0) - m.unsqueeze(2) * q.unsqueeze(1) + q.unsqueeze(2) * m.unsqueeze(1)
    Jqw = c * m.unsqueeze(2) * R.unsqueeze(1)
    Jwq = -(1.0 / c) * R.unsqueeze(2) * m.unsqueeze(1)
    Jww = torch.zeros(len(p), 2, 2, dtype=p.dtype, device=p.device)
    J = torch.cat([torch.cat([Jqq, Jqw], dim=2), torch.cat([Jwq, Jww], dim=2)], dim=1)
    if not project:
        return J
    P = _projector(p, d)
    return P @ J @ P


def _projector(p: torch.Tensor, d: int) -> torch.Tensor:
    """Block projector for S^{d−1} × S¹ (same formula as AmbientManifold)."""
    from .ambient import AmbientManifold
    return AmbientManifold((d - 1, 1)).projector(p)
