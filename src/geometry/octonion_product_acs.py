"""
Explicit almost complex structure on S⁴×S² ⊂ ℝ⁵×ℝ³ = ℝ⁸ ≅ 𝕆 via octonion
left multiplication.

At (x, y) put a = (x, 0), b = (0, y) — orthogonal unit octonions (the two
radial normals of the product). Then u := b·ā is a unit IMAGINARY octonion
(Re(b ā) = ⟨b, a⟩ = 0), so L_u: v ↦ u·v is an orthogonal complex structure on
ℝ⁸ (u² = −1 and alternativity give L_u² = −1). It maps a ↦ b ↦ −a, hence
preserves the normal plane N = span(a, b) and therefore its orthogonal
complement T = T_{(x,y)}(S⁴×S²). The restriction

    J_ref(x, y) = P · L_u · P

is a smooth global ACS on S⁴×S² — the classical octonionic reason
S²×S⁴ carries an almost complex structure.

Roles in this project (diagnostics + base point, never a CS target):
  * certified exact ACS for the S⁴×S² row (chart-protocol ACS loss → 0);
  * reference Nijenhuis value ("octonionic-product basin");
  * base point for the CS search J = e^B J_ref e^{−B} (open problem),
    exactly as J_oct serves on S⁶.
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import torch

from .octonion_acs import cross7


def oct_mul(a: torch.Tensor, b: torch.Tensor) -> torch.Tensor:
    """Octonion product on ℝ⁸ = ℝ ⊕ Im𝕆, batched (B, 8) × (B, 8) → (B, 8)."""
    a0, av = a[:, :1], a[:, 1:]
    b0, bv = b[:, :1], b[:, 1:]
    real = a0 * b0 - (av * bv).sum(dim=1, keepdim=True)
    imag = a0 * bv + b0 * av + cross7(av, bv)
    return torch.cat([real, imag], dim=1)


def oct_conj(a: torch.Tensor) -> torch.Tensor:
    return torch.cat([a[:, :1], -a[:, 1:]], dim=1)


def left_mult_matrix(u: torch.Tensor) -> torch.Tensor:
    """Matrix of L_u: v ↦ u·v, shape (B, 8, 8)."""
    B = u.shape[0]
    eye = torch.eye(8, dtype=u.dtype, device=u.device)
    cols = [oct_mul(u, eye[k].expand(B, 8)) for k in range(8)]
    return torch.stack(cols, dim=2)                      # column k = u·e_k


def octonion_product_acs(p: torch.Tensor, projector) -> torch.Tensor:
    """J_ref at ambient points p = (x, y) ∈ S⁴×S² ⊂ ℝ⁸; `projector` is the
    manifold's tangential projector function. Returns (B, 8, 8)."""
    x, y = p[:, :5], p[:, 5:]
    a = torch.cat([x, torch.zeros_like(y)], dim=1)
    b = torch.cat([torch.zeros_like(x), y], dim=1)
    u = oct_mul(b, oct_conj(a))                          # unit imaginary
    L = left_mult_matrix(u)
    P = projector(p)
    return P @ L @ P
