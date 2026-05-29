"""
Analytic Hopf complex structure on S^{m_α} × S¹  —  FIXED VERSION.

Drop-in alternative to `hopf_complex_structure.py` (which is left untouched).
The public API is identical (same function names, signatures, semantics), so
swapping consists of changing one import path, e.g. in train.py:

    # original (folded co-latitude — discontinuous):
    from src.geometry.hopf_complex_structure import hopf_complex_structure

    # fixed (signed angle — smooth):
    from src.geometry.hopf_complex_structure_fixed import hopf_complex_structure

----------------------------------------------------------------------------
The single mathematical change vs. `hopf_complex_structure.py`:
how the scalar Hopf log-radius coordinate t(x_β) is computed.
----------------------------------------------------------------------------

The Hopf manifold M = (ℂⁿ \ {0}) / (z ∼ e·z) is diffeomorphic to S^{2n−1} × S¹.
The standard complex structure J₀ on ℂⁿ \ {0} (multiplication by i) is
invariant under the deck transformation z ↦ e·z, so it descends to a SMOOTH
almost complex structure J on M.  In chart coordinates the codebase encodes
this via the smooth immersion

    F : B^{m_α} × B¹  →  ℂⁿ ≅ ℝ^{2n},      F(x_α, x_β) = e^{t(x_β)} · σ(x_α),

where σ : B^{m_α} → S^{m_α} ⊂ ℝ^{2n} is the codebase's inverse stereographic
embedding (the same one used by `EvenSphere.ball_to_sphere`), and t is the
log-radius coordinate on the Hopf circle ℝ/ℤ.  The complex structure on the
chart is then the pullback

    J(x) = (DF(x))⁻¹ · J₀ · DF(x).

Original (folded) version of t:
    t(x_β) = (1/2π) · arccos( sphere_z(x_β²) )
where sphere_z(s) = (−1 + 6s − s²)/(1 + s)² is the latitude (height in ℝ²)
of the β-point on S¹.  This is the *co-latitude* — the unsigned angular
distance from the north pole.  On a circle it is two-to-one in position
(±θ from a pole share a latitude), so:
    - t depends only on x_β² (it is EVEN in x_β),
    - F(x_α, +x_β) = F(x_α, −x_β) identically,
    - the x_β-column of DF flips sign across x_β = 0,
    - J(−x_β) = S · J(+x_β) · S with S = diag(I_{m_α}, −1),
    - hence J has a JUMP DISCONTINUITY at x_β = 0
      (verified Frobenius norm ≈ 8.8 on S³×S¹, machine precision; see
       /tmp/hopf_verify.py from the diagnostic thread).
Moreover t ∈ (0, 1/2) covers only HALF the Hopf circle.  The function
that the original module returns is therefore not the smooth Hopf complex
structure on M but a folded, half-period, discontinuous patchwork.

Fixed version of t  (this module):
    t(x_β) = (1/2π) · atan2( sphere_xy(x_β),  −sphere_z(x_β²) )
where sphere_xy(x_β) = 4 · x_β · (1 − x_β²) / (1 + x_β²)² is the longitude
of the β-point (odd in x_β; the embedding (sphere_xy, sphere_z) ∈ S¹ ⊂ ℝ²
agrees with `EvenSphere.ball_to_sphere` for d=1).  This is the GENUINE
SIGNED angular coordinate on the Hopf circle: smooth and monotone in x_β
throughout the chart interior, covering the full period (−1/2, 1/2).  Its
only discontinuity is the atan2 branch cut at the north pole (s_z=+1),
i.e. at the chart boundary |x_β| → 1, which the pole filter already
excludes; furthermore the branch-cut jump there is Δt = 1, i.e. exactly
one Hopf period, under which J = (DF)⁻¹·J₀·DF is invariant (J₀ is
scale-invariant — multiplication by i commutes with rescaling — so the
deck-symmetry of M acts trivially on J).

With this t, the formula F = e^t · σ yields a smooth diffeomorphism from
the chart onto an open subset of ℂⁿ \ {0}, and J = (DF)⁻¹·J₀·DF is the
actual smooth Hopf complex structure on M expressed in chart coordinates.

Verified on S³×S¹ (and structurally identical on S^{m_α}×S¹ for any odd
m_α; verifications in /tmp/hopf_fix_audit.py from the diagnostic thread):
    - jump at x_β = 0:               1e-15  (machine precision; was 8.79)
    - ‖J² + I‖_F                  :   1e-16
    - Nijenhuis ‖N_J‖_F           :   1e-15  on the sampled band
    - overlap consistency
       ‖F·J_i·F⁻¹ − J_j‖_F        :   1e-15  for all 6 patch pairs.

Convention (matches `hopf_complex_structure.py`):
    x_alpha    : (batch, dim_alpha)  ball coords for S^{dim_alpha}  (α factor)
    x_beta     : (batch, 1)          ball coord  for S¹             (β factor)
    dim_alpha  = m_α = 2n − 1        must be odd
    d          = dim_alpha + 1 = 2n  total real dimension (even)
    pole_α     ∈ {0, 1}              chart index for S^{m_α}
    pole_β     ∈ {0, 1}              chart index for S¹
"""

import torch
from typing import Dict, Tuple

from .sphere import EvenSphere
# Helpers are unchanged from the original module; reuse them rather than
# duplicating, so the only divergence is the t-formula inside hopf_map_F.
from .hopf_complex_structure import _sphere_z, _J0


# ---------------------------------------------------------------------------
# Core map F  (the fix lives in the computation of t)
# ---------------------------------------------------------------------------

def hopf_map_F(
    x_alpha: torch.Tensor,
    x_beta:  torch.Tensor,
    pole_alpha: int,
    pole_beta:  int,
) -> torch.Tensor:
    """Map (x_α, x_β) ∈ B^{m_α} × B¹ to ℂⁿ ≅ ℝ^{2n} via the Hopf construction,
    using the SIGNED-ANGLE Hopf log-radius coordinate (see module docstring).

    Same signature and semantics as `hopf_complex_structure.hopf_map_F`,
    differing only in the computation of t(x_β):
        original :  t = arccos(sphere_z(x_β²)) / (2π)                (folded)
        fixed    :  t = atan2(sphere_xy(x_β), −sphere_z(x_β²)) / (2π) (signed)

    F(x_α, x_β) = exp(t(x_β)) · σ(x_α),  where σ is the codebase's inverse
    stereographic map B^{m_α} → S^{m_α} ⊂ ℝ^{2n}, namely
        σ(x_α) = ( sphere_z(r_α²),  4(1−r_α²)/(1+r_α²)² · x_α )
        sphere_z(s) = (−1 + 6s − s²) / (1+s)².

    Args:
        x_alpha:    Ball coords for S^{m_α}, shape (batch, m_α).  m_α must be odd.
        x_beta:     Ball coord  for S¹,      shape (batch, 1).
        pole_alpha: Chart index n ∈ {0, 1} for S^{m_α}.
        pole_beta:  Chart index m ∈ {0, 1} for S¹.

    Returns:
        F_out: Shape (batch, m_α + 1) = (batch, d).
    """
    # For pole=1 charts, apply patch_change_coords to map to pole-0 coords first,
    # so that F_{(ia,ib)}(x) = F_{(0,0)}(φ^{ia}(xa), φ^{ib}(xb)).  φ is self-inverse,
    # so the complex structures in all four patches are pullbacks of the same
    # global F, automatically satisfying the overlap conditions.
    xb_in = EvenSphere.patch_change_coords(x_beta) if pole_beta == 1 else x_beta
    xa_in = EvenSphere.patch_change_coords(x_alpha) if pole_alpha == 1 else x_alpha

    # --- β factor: signed-angle Hopf log-radius t ---
    xb        = xb_in[:, 0]                                    # (batch,)
    r_beta_sq = xb ** 2
    sz_beta   = _sphere_z(r_beta_sq)                           # latitude  s_z (even in x_β)
    sxy_beta  = 4.0 * xb * (1.0 - r_beta_sq) / (1.0 + r_beta_sq) ** 2
                                                               # longitude s_xy (odd in x_β);
                                                               # (s_xy, s_z) ∈ S¹ ⊂ ℝ²
    # Genuine signed circle angle of the β-point, measured from the south pole.
    # Branch cut of atan2 is at second arg < 0, first arg = 0  ⇒  s_z = +1  ⇒
    # north pole, i.e. |x_β| → 1, the chart boundary (already excluded by the
    # pole filter).  The branch-cut jump there is Δt = 1 = one Hopf period,
    # under which J = (DF)⁻¹·J₀·DF is invariant (deck symmetry z ∼ e·z), so J
    # is smooth on the entire chart interior, all four patches.
    t      = torch.atan2(sxy_beta, -sz_beta) / (2.0 * torch.pi)   # (batch,) ∈ (−1/2, 1/2)
    exp_t  = torch.exp(t)                                          # (batch,)

    # --- α factor: inverse stereographic σ(x_α) ∈ S^{m_α} ⊂ ℝ^{m_α+1}  (unchanged) ---
    r_alpha_sq = (xa_in ** 2).sum(dim=1)                       # (batch,)
    sz_alpha   = _sphere_z(r_alpha_sq)                         # leading coordinate
    leading    = sz_alpha

    denom   = (1.0 + r_alpha_sq) ** 2                          # (batch,)
    v_alpha = 4.0 * (1.0 - r_alpha_sq).unsqueeze(1) / denom.unsqueeze(1) * xa_in
                                                               # (batch, m_α)

    # F = exp(t) · σ(x_α)  where σ(x_α) = (leading, v_alpha)
    F_lead = (exp_t * leading).unsqueeze(1)                    # (batch, 1)
    F_rest = exp_t.unsqueeze(1) * v_alpha                      # (batch, m_α)
    return torch.cat([F_lead, F_rest], dim=1)                  # (batch, d)


# ---------------------------------------------------------------------------
# J = (DF)⁻¹ · J₀ · DF
# ---------------------------------------------------------------------------

def hopf_complex_structure(
    x_alpha: torch.Tensor,
    x_beta:  torch.Tensor,
    pole_alpha: int,
    pole_beta:  int,
) -> torch.Tensor:
    """Analytic Hopf complex structure J = (DF)⁻¹·J₀·DF in ball coordinates.

    Identical signature/semantics to `hopf_complex_structure.hopf_complex_structure`,
    except that DF is the Jacobian of THIS module's `hopf_map_F` (signed-angle t).
    Consequently J is the genuine smooth Hopf complex structure on M expressed
    in chart coordinates, not the folded, half-period, discontinuous patchwork
    produced by the original module.

    The Jacobian DF is computed by forward-mode AD (`torch.func.jacfwd`),
    consistently with the rest of the codebase.

    Args:
        x_alpha:    Ball coords for S^{m_α}, shape (batch, m_α).
        x_beta:     Ball coord  for S¹,      shape (batch, 1).
        pole_alpha: Chart n ∈ {0, 1}.
        pole_beta:  Chart m ∈ {0, 1}.

    Returns:
        J: Shape (batch, d, d) where d = m_α + 1.
    """
    batch     = x_alpha.shape[0]
    dim_alpha = x_alpha.shape[1]
    d         = dim_alpha + 1                                  # even

    # Concatenate input for Jacobian computation
    x_combined = torch.cat([x_alpha, x_beta], dim=1)           # (batch, d)

    # Single-sample F (so vmap can supply the batch axis)
    def _F(x: torch.Tensor) -> torch.Tensor:
        xa = x[:dim_alpha].unsqueeze(0)
        xb = x[dim_alpha:].unsqueeze(0)
        return hopf_map_F(xa, xb, pole_alpha, pole_beta)[0]

    # Batched forward-mode Jacobian: DF[b] = ∂F/∂x at x_combined[b]
    DF = torch.func.vmap(torch.func.jacfwd(_F))(x_combined)    # (batch, d, d)

    # Canonical J₀ on ℝ^d
    J0       = _J0(d, dtype=x_alpha.dtype, device=x_alpha.device)
    J0_batch = J0.unsqueeze(0).expand(batch, -1, -1)            # (batch, d, d)

    # J = (DF)⁻¹ · J₀ · DF
    DF_inv = torch.linalg.inv(DF)                               # (batch, d, d)
    return DF_inv @ J0_batch @ DF                               # (batch, d, d)


# ---------------------------------------------------------------------------
# All-patches interface
# ---------------------------------------------------------------------------

def hopf_j_all_patches(
    x_alpha: torch.Tensor,
    x_beta:  torch.Tensor,
) -> Dict[Tuple[int, int], torch.Tensor]:
    """Compute the Hopf complex structure in all four ball-patch charts.

    Identical signature/semantics to `hopf_complex_structure.hopf_j_all_patches`.

    Args:
        x_alpha: Patch-0 ball coords for S^{m_α}, shape (batch, m_α).
        x_beta:  Patch-0 ball coord  for S¹,      shape (batch, 1).

    Returns:
        Dict mapping (ia, ib) ↦ J of shape (batch, d, d).
    """
    x_a2 = EvenSphere.patch_change_coords(x_alpha)
    x_b2 = EvenSphere.patch_change_coords(x_beta)
    patch_alpha = {0: x_alpha, 1: x_a2}
    patch_beta  = {0: x_beta,  1: x_b2}

    result: Dict[Tuple[int, int], torch.Tensor] = {}
    for ia in range(2):
        for ib in range(2):
            result[(ia, ib)] = hopf_complex_structure(
                patch_alpha[ia], patch_beta[ib], ia, ib
            )
    return result
