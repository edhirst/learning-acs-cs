"""
Analytic Hopf complex structure on S^{m_α} × S¹.

The Hopf manifold construction gives an integrable (N_J = 0) complex structure
on S^{m_α} × S¹ for any odd m_α = 2n−1, n ≥ 1.  Concretely M = (ℂⁿ)* / z ~ e·z
is diffeomorphic to S^{2n-1} × S¹, and the standard complex structure on ℂⁿ
descends to M.

This module implements J in ball coordinates for the product of patches:
    (pole_alpha, pole_beta) ∈ {0,1}²

where pole_alpha is the chart for S^{m_α} and pole_beta is the chart for S¹.

Convention matching the rest of the codebase:
    x_alpha: (batch, dim_alpha)  — ball coords for S^{dim_alpha}  (α factor)
    x_beta:  (batch, 1)          — ball coord  for S¹              (β factor)
    dim_alpha = m_α = 2n−1       — must be odd
    d = dim_alpha + 1 = 2n       — total real dimension

The map F: B^{m_α} × B¹ → ℂⁿ ≅ ℝ^{2n} is:

    F(x_α, x_β) = exp(t(x_β)) · σ(x_α)

where σ: B^{m_α} → S^{m_α} ⊂ ℝ^{2n} is the codebase's inverse stereographic map
(matching ball_to_sphere in sphere.py), and t(x_β) is the Hopf log-scale from S¹:

    r_α    = |x_α|
    σ(x_α) = ( (−1)^n · sphere_z(r_α),  4(1−r_α²)/(1+r_α²)² · x_α )
    sphere_z(r²) = (−1 + 6r² − r⁴) / (1 + r²)²   ∈ [−1, 1]

    r_β    = |x_β|   (scalar, x_β ∈ B¹)
    t      = (1/2π) arccos( (−1)^m · sphere_z(r_β) )  ∈ (0, 1/2)

    F = exp(t) · σ(x_α)                                # shape (2n,) = (d,)

The complex structure is then:
    J = (DF)^{-1} · J₀ · DF

where J₀ is the canonical block-diagonal structure on ℝ^{2n} and DF is the
Jacobian of F with respect to (x_α, x_β).

The Jacobian DF is computed using torch.func.jacfwd (forward-mode AD), which
is consistent with the rest of the codebase (_compute_A_and_dA in sphere_loss.py).
"""

import torch
import torch.nn.functional as F
from typing import Dict, Tuple

from .sphere import EvenSphere


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _sphere_z(r_sq: torch.Tensor) -> torch.Tensor:
    """Latitude component of inverse stereographic map in codebase ball coordinates.

    For ball coords x with r² = |x|², maps to the last component of the sphere
    point in ℝ^{d+1}, matching ball_to_sphere() in sphere.py:

        sphere_z = (−1 + 6r² − r⁴) / (1 + r²)²

    Range: sphere_z(0) = −1 (south pole), sphere_z(r→1) → +1 (north pole).

    Args:
        r_sq: Squared radii, shape (batch,) or scalar.

    Returns:
        sphere_z: shape matching r_sq, values in [−1, 1].
    """
    num = -1.0 + 6.0 * r_sq - r_sq ** 2
    den = (1.0 + r_sq) ** 2
    return num / den


def _J0(d: int, dtype: torch.dtype, device: torch.device) -> torch.Tensor:
    """Canonical block-diagonal J₀ on ℝ^d, d even.

    J₀ = block_diag([[0, -1], [1, 0]], ...).

    Args:
        d:      Even real dimension.
        dtype:  Tensor dtype.
        device: Tensor device.

    Returns:
        J₀ of shape (d, d).
    """
    J0 = torch.zeros(d, d, dtype=dtype, device=device)
    for k in range(d // 2):
        i, j = 2 * k, 2 * k + 1
        J0[i, j] = -1.0
        J0[j, i] =  1.0
    return J0


# ---------------------------------------------------------------------------
# Core map F
# ---------------------------------------------------------------------------

def hopf_map_F(
    x_alpha: torch.Tensor,
    x_beta:  torch.Tensor,
    pole_alpha: int,
    pole_beta:  int,
) -> torch.Tensor:
    """Map (x_α, x_β) ∈ B^{m_α} × B¹ to ℂⁿ ≅ ℝ^{2n} via the Hopf construction.

    F(x_α, x_β) = exp(t(x_β)) · σ(x_α) where σ: B^{m_α} → S^{m_α} ⊂ ℝ^{2n} is
    the codebase inverse stereographic map and t is the Hopf log-scale:

        sphere_z(r²) = (−1 + 6r² − r⁴) / (1 + r²)²   ∈ [−1, 1]
        σ(x_α)  = ( (−1)^{pole_α} · sphere_z(r_α²),  4(1−r_α²)/(1+r_α²)² · x_α )
        t(x_β)  = (1/2π) arccos( (−1)^{pole_β} · sphere_z(r_β²) )
        F       = exp(t) · σ(x_α)

    Args:
        x_alpha:    Ball coords for S^{m_α}, shape (batch, m_α).  m_α must be odd.
        x_beta:     Ball coord  for S¹,      shape (batch, 1).
        pole_alpha: Chart index n ∈ {0, 1} for S^{m_α} (controls sign of latitude).
        pole_beta:  Chart index m ∈ {0, 1} for S¹ (controls sign in arccos arg).

    Returns:
        F_out: Shape (batch, m_α + 1) = (batch, d).
    """
    # For pole=1 charts, apply patch_change_coords to map to pole-0 coords first.
    # This ensures F_{(ia,ib)}(x) = F_{(0,0)}(φ^{ia}(xa), φ^{ib}(xb)) where φ is
    # self-inverse, so the complex structures in all patches are pullbacks of J₀.
    xb_in = EvenSphere.patch_change_coords(x_beta) if pole_beta == 1 else x_beta
    xa_in = EvenSphere.patch_change_coords(x_alpha) if pole_alpha == 1 else x_alpha

    # --- β factor: compute Hopf log-scale t from S¹ latitude ---
    xb = xb_in[:, 0]                               # (batch,)
    r_beta_sq = xb ** 2
    sz_beta = _sphere_z(r_beta_sq)                 # (batch,) in (−1, 1)
    arg_arccos = sz_beta.clamp(-1.0 + 1e-6, 1.0 - 1e-6)
    t = torch.arccos(arg_arccos) / (2.0 * torch.pi)  # (batch,) in (0, 1/2)
    exp_t = torch.exp(t)                              # (batch,)

    # --- α factor: inverse stereo σ(x_α) ∈ S^{m_α} ⊂ ℝ^{m_α+1} ---
    r_alpha_sq = (xa_in ** 2).sum(dim=1)           # (batch,)
    sz_alpha = _sphere_z(r_alpha_sq)               # (batch,) — latitude of σ(x_α)
    leading = sz_alpha                             # (batch,) — first output component

    # Remaining m_α components: sphere_xy = 4(1 − r²)/(1 + r²)² · x_α
    denom = (1.0 + r_alpha_sq) ** 2               # (batch,)
    v_alpha = 4.0 * (1.0 - r_alpha_sq).unsqueeze(1) / denom.unsqueeze(1) * xa_in  # (batch, m_α)

    # F = exp(t) · σ(x_α)  where σ(x_α) = (leading, v_alpha)
    # Shape (batch, 1 + m_α) = (batch, d)
    F_lead = (exp_t * leading).unsqueeze(1)        # (batch, 1)
    F_rest = exp_t.unsqueeze(1) * v_alpha          # (batch, m_α)
    return torch.cat([F_lead, F_rest], dim=1)      # (batch, d)


# ---------------------------------------------------------------------------
# DF = Jacobian of F
# ---------------------------------------------------------------------------

def hopf_map_DF(
    x_alpha: torch.Tensor,
    x_beta:  torch.Tensor,
    pole_alpha: int,
    pole_beta:  int,
) -> torch.Tensor:
    """Jacobian DF of hopf_map_F w.r.t. (x_α, x_β) via forward-mode AD.

    Args:
        x_alpha:    Ball coords for S^{m_α}, shape (batch, m_α).
        x_beta:     Ball coord  for S¹,      shape (batch, 1).
        pole_alpha: Chart n ∈ {0, 1}.
        pole_beta:  Chart m ∈ {0, 1}.

    Returns:
        DF: Shape (batch, d, d) where d = m_α + 1.
            DF[b, i, j] = ∂F_i / ∂x_j at (x_alpha[b], x_beta[b]).
    """
    dim_alpha = x_alpha.shape[1]
    x_combined = torch.cat([x_alpha, x_beta], dim=1)  # (batch, d)

    def _F(x: torch.Tensor) -> torch.Tensor:
        xa = x[:dim_alpha].unsqueeze(0)
        xb = x[dim_alpha:].unsqueeze(0)
        return hopf_map_F(xa, xb, pole_alpha, pole_beta)[0]

    return torch.func.vmap(torch.func.jacfwd(_F))(x_combined)  # (batch, d, d)


# ---------------------------------------------------------------------------
# J = (DF)^{-1} J₀ DF
# ---------------------------------------------------------------------------

def hopf_complex_structure(
    x_alpha: torch.Tensor,
    x_beta:  torch.Tensor,
    pole_alpha: int,
    pole_beta:  int,
) -> torch.Tensor:
    """Analytic Hopf complex structure J = (DF)⁻¹ J₀ DF in ball coordinates.

    Computes the Jacobian DF of hopf_map_F via forward-mode AD (jacfwd), then
    assembles J via matrix inversion and multiplication with the canonical J₀.

    Args:
        x_alpha:    Ball coords for S^{m_α}, shape (batch, m_α).
        x_beta:     Ball coord  for S¹,      shape (batch, 1).
        pole_alpha: Chart n ∈ {0, 1}.
        pole_beta:  Chart m ∈ {0, 1}.

    Returns:
        J: Shape (batch, d, d) where d = m_α + 1.
    """
    batch = x_alpha.shape[0]
    dim_alpha = x_alpha.shape[1]
    d = dim_alpha + 1   # total real dimension (even)

    # Concatenate input for Jacobian computation
    x_combined = torch.cat([x_alpha, x_beta], dim=1)  # (batch, d)

    # Define F as a function of the combined input
    def _F(x: torch.Tensor) -> torch.Tensor:
        # x: (d,)  — single sample (used by vmap/jacfwd)
        xa = x[:dim_alpha].unsqueeze(0)   # (1, m_α)
        xb = x[dim_alpha:].unsqueeze(0)   # (1, 1)
        return hopf_map_F(xa, xb, pole_alpha, pole_beta)[0]   # (d,)

    # Batched forward-mode Jacobian: DF[b] = ∂F/∂x at x_combined[b]
    # Shape: (batch, d, d)  where DF[b, i, j] = ∂F_i/∂x_j
    DF = torch.func.vmap(torch.func.jacfwd(_F))(x_combined)   # (batch, d, d)

    # J₀ on ℝ^d
    J0 = _J0(d, dtype=x_alpha.dtype, device=x_alpha.device)   # (d, d)
    J0_batch = J0.unsqueeze(0).expand(batch, -1, -1)           # (batch, d, d)

    # J = DF^{-1} · J₀ · DF
    DF_inv = torch.linalg.inv(DF)                              # (batch, d, d)
    J = DF_inv @ J0_batch @ DF                                 # (batch, d, d)
    return J


# ---------------------------------------------------------------------------
# All-patches interface
# ---------------------------------------------------------------------------

def hopf_j_all_patches(
    x_alpha: torch.Tensor,
    x_beta:  torch.Tensor,
) -> Dict[Tuple[int, int], torch.Tensor]:
    """Compute the Hopf complex structure in all four ball-patch charts.

    For patch (ia, ib), J is evaluated at patch-specific coordinates obtained
    by applying EvenSphere.patch_change_coords to the pole=1 factor(s).  Since
    hopf_complex_structure internally reapplies patch_change_coords for pole=1,
    the effective input to F is always in pole-0 ball coordinates, so J_{(ia,ib)}
    is the correct pullback of J₀ to the (ia,ib) local chart.

    This is consistent with what ProductOverlapLoss expects: J_dict[(ia,ib)] at
    the (ia,ib)-patch coordinates of the physical points.

    Args:
        x_alpha: Patch-0 ball coords for S^{m_α}, shape (batch, m_α).
        x_beta:  Patch-0 ball coord  for S¹,      shape (batch, 1).

    Returns:
        Dict mapping (ia, ib) → J tensor of shape (batch, d, d).
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
