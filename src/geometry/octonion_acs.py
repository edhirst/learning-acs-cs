"""
The octonionic (G₂-invariant) almost complex structure on S⁶, expressed in
the codebase's two ball charts.

Theory
------
Identify ℝ⁷ ≅ Im 𝕆 (imaginary octonions). For x ∈ S⁶ and v ∈ T_x S⁶ (v ⊥ x),

    J_x(v) = x × v,

the 7-dimensional cross product induced by octonion multiplication. Standard
identities of the cross product give x × (x × v) = ⟨x,v⟩x − |x|²v = −v on the
tangent space, so J² = −I. This is the classical nearly-Kähler ACS on the
round S⁶: G₂-homogeneous, with NONZERO Nijenhuis tensor of constant
round-metric norm (it is famously non-integrable; LeBrun-type results forbid
any round-orthogonal integrable structure).

Purpose in this project (diagnostics only — never a training target):
  * ACS loss-floor calibration for the sphere-side loss stack in d = 6
    (the sphere analogue of the Calabi–Eckmann floor calibration on products);
  * a reference Nijenhuis value for the "octonionic basin" — comparing trained
    S⁶ fields against it tests whether free ACS-guided training converges to
    the round basin, far from the (Alpöge–Claude/Engel 2026) integrable
    structures, which are torus fibrations metrically unrelated to it;
  * a d = 6 sphere-side target for Sobolev-fit architecture diagnostics.

Chart conventions
-----------------
Pullback through `EvenSphere.ball_to_sphere(x, from_north)`: patch 1 uses
from_north=True, patch 2 uses from_north=False (they embed the same physical
point p at x₁ and x₂ = φ(x₁) respectively). With Dσ the chart differential,

    J_chart(x) = (DσᵀDσ)⁻¹ Dσᵀ · M(σ(x)) · Dσ,     M(p)v = p × v.

Both patches are pullbacks of the same global J, so the cross-chart overlap
condition J₂ = F J₁ F⁻¹ holds by construction.

All operations are differentiable torch (einsum + torch.func.vmap(jacfwd) +
a 6×6 inverse), so autograd through J_chart yields correct ∂J/∂x.
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import torch

from .sphere import EvenSphere


# ---------------------------------------------------------------------------
# Octonion structure constants (7-dim cross product)
# ---------------------------------------------------------------------------

def _build_structure_constants() -> torch.Tensor:
    """f[i,j,k] with (x × v)_k = f_{ijk} x_i v_j, from the Fano-plane triples
    (i, i+1, i+3) mod 7 (1-based). Verified to satisfy the composition-algebra
    identities x×(x×v) = ⟨x,v⟩x − |x|²v and |x×v| = |x||v| for x ⊥ v."""
    f = torch.zeros(7, 7, 7, dtype=torch.float64)
    triples = [(1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 7),
               (5, 6, 1), (6, 7, 2), (7, 1, 3)]
    for (i, j, k) in triples:
        i, j, k = i - 1, j - 1, k - 1
        for (a, b, c) in [(i, j, k), (j, k, i), (k, i, j)]:
            f[a, b, c] = 1.0
        for (a, b, c) in [(j, i, k), (i, k, j), (k, j, i)]:
            f[a, b, c] = -1.0
    return f


_F_CONST = _build_structure_constants()


def cross7(x: torch.Tensor, v: torch.Tensor) -> torch.Tensor:
    """Batched 7-dim cross product, shapes (batch, 7) × (batch, 7) → (batch, 7)."""
    f = _F_CONST.to(dtype=x.dtype, device=x.device)
    return torch.einsum('ijk,bi,bj->bk', f, x, v)


def cross7_matrix(p: torch.Tensor) -> torch.Tensor:
    """M(p) with M(p) @ v = p × v; shape (batch, 7) → (batch, 7, 7)."""
    f = _F_CONST.to(dtype=p.dtype, device=p.device)
    return torch.einsum('ijk,bi->bkj', f, p)


# ---------------------------------------------------------------------------
# Chart pullback
# ---------------------------------------------------------------------------

def chart_embedding(x: torch.Tensor, patch: int) -> torch.Tensor:
    """Ambient embedding of a patch's ball coordinates, consistent with the
    codebase's gluing map φ = EvenSphere.patch_change_coords.

    Patch 1 is `ball_to_sphere(x, from_north=True)`. Patch 2 must satisfy
    σ₂(φ(x₁)) = σ₁(x₁) for the same physical point; because φ reverses the
    tangential direction (x' = −x̂·(1−ρ)/(1+ρ)), the raw
    `ball_to_sphere(·, from_north=False)` differs from σ₁∘φ by the ambient
    isometry R = diag(−1,…,−1, 1) (negate the 6 equatorial coordinates, keep
    the pole axis). We therefore take σ₂ := R ∘ ball_to_sphere(·, False),
    which agrees with σ₁∘φ on the overlap yet stays smooth on the whole
    patch-2 ball (σ₁∘φ itself degenerates at the patch-2 centre).
    """
    if patch == 1:
        return EvenSphere.ball_to_sphere(x, from_north=True)
    p = EvenSphere.ball_to_sphere(x, from_north=False)
    return torch.cat([-p[:, :-1], p[:, -1:]], dim=1)


def _embedding_jacobian(x: torch.Tensor, patch: int) -> torch.Tensor:
    """Dσ of chart_embedding, shape (batch, 7, dim) via forward AD."""
    def sigma_single(x_single: torch.Tensor) -> torch.Tensor:
        return chart_embedding(x_single.unsqueeze(0), patch).squeeze(0)
    return torch.func.vmap(torch.func.jacfwd(sigma_single))(x)


def octonion_acs(x: torch.Tensor, patch: int = 1) -> torch.Tensor:
    """The octonionic ACS in ball coordinates of the given patch.

    Both patches are pullbacks of the one global J through embeddings that
    agree on overlaps (see chart_embedding), so J₂ = F J₁ F⁻¹ holds exactly.

    Args:
        x: Ball coordinates, shape (batch, 6). Must lie in the open ball
           (the pullback degenerates as |x| → 1, like every chart quantity).
        patch: 1 (north-pole-removed embedding) or 2 (south-pole-removed).

    Returns:
        J of shape (batch, 6, 6) with J² = −I to machine precision.
    """
    if patch not in (1, 2):
        raise ValueError(f"patch must be 1 or 2, got {patch}")
    p = chart_embedding(x, patch)                               # (batch, 7)
    D = _embedding_jacobian(x, patch)                           # (batch, 7, 6)
    M = cross7_matrix(p)                                        # (batch, 7, 7)
    DtD = D.transpose(1, 2) @ D                                 # (batch, 6, 6)
    D_pinv = torch.linalg.inv(DtD) @ D.transpose(1, 2)          # (batch, 6, 7)
    return D_pinv @ M @ D                                       # (batch, 6, 6)


def octonion_acs_both_patches(x1: torch.Tensor):
    """Convenience: (J1 at patch-1 coords x1, J2 at x2 = φ(x1), x2)."""
    x2 = EvenSphere.patch_change_coords(x1)
    return octonion_acs(x1, patch=1), octonion_acs(x2, patch=2), x2
