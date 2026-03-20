"""
Nijenhuis tensor and overlap consistency losses for S^{m_α} × S^{m_β}.

The product is covered by four patches (i_α, i_β) ∈ {0,1}².  This gives:
  - 4 Nijenhuis loss terms, one per patch.
  - C(4,2) = 6 overlap (patch-consistency) loss terms, one per patch pair.

Radial weighting per patch
--------------------------
Each patch (i_α, i_β) inherits a separate radial weight from each factor:

    w^{(ia,ib)}(ρ_α, ρ_β) = w_α(ρ_α^{patch}) · w_β(ρ_β^{patch})

where ρ_α^{patch} = ‖x_α‖ if i_α=0 else ‖φ_α(x_α)‖, and similarly for β.
The Gaussian decay and plateau radii are per-sphere and configurable
(overlap_size_alpha, overlap_size_beta).

Block Jacobian efficiency for overlap losses
--------------------------------------------
The six patch pairs group by which factors are flipped between the two patches:
  β-only flip  (Δi_α=0, Δi_β=1): (00,01), (10,11)  — use F_β only
  α-only flip  (Δi_α=1, Δi_β=0): (00,10), (01,11)  — use F_α only
  both-flip    (Δi_α=1, Δi_β=1): (00,11), (01,10)  — use block(F_α, F_β)

For each pair the transition Jacobian F is evaluated at the *source* patch
coords. The α and β Jacobians (and their inverses) are computed once each
and reused.  The full-block Jacobian is never assembled; instead the efficient
block_jacobian_transform is called (see SphereProduct.block_jacobian_transform).

Note on the diagonal pair (01)↔(10): this pair involves a both-flip but the
α-Jacobian is evaluated at x_α (from patch-00), and the β-Jacobian at x_b2
(from patch-01's β-coords = φ_β(x_β)), because patch-01 uses x_b2 as its
β-coords. A separate Jacobian F_beta_at_xb2 is computed for this pair.
"""

import torch
import torch.autograd.forward_ad as fwdAD
import torch.nn as nn
from typing import Dict, List, Optional, Tuple

from src.geometry import RADIAL_MIDPOINT
from src.geometry.sphere import EvenSphere
from src.geometry.sphere_product import SphereProduct
from src.loss.sphere_loss import NijenhuisTensorLoss, _metric_norm_sq_12tensor, _metric_norm_sq_11tensor


# All 6 ordered patch pairs (source, target, flip_type)
# flip_type: 'alpha', 'beta', or 'both'
_OVERLAP_PAIRS: List[Tuple[Tuple[int,int], Tuple[int,int], str]] = [
    ((0, 0), (0, 1), 'beta'),
    ((0, 0), (1, 0), 'alpha'),
    ((0, 0), (1, 1), 'both'),
    ((0, 1), (1, 1), 'alpha'),
    ((1, 0), (1, 1), 'beta'),
    ((0, 1), (1, 0), 'both'),
]


class ProductNijenhuisLoss(nn.Module):
    """
    Per-patch Nijenhuis tensor loss for S^{m_α} × S^{m_β}.

    Computes the Nijenhuis tensor on each of the four patches, using the
    product radial weight w_α(ρ_α) · w_β(ρ_β) to down-weight near-pole
    regions of either factor.
    """

    def __init__(
        self,
        total_dim: int,
        overlap_size_alpha: float = 0.1,
        overlap_size_beta:  float = 0.1,
        use_round_metric_contraction: bool = False,
    ):
        """
        Args:
            total_dim: m_α + m_β (even).
            overlap_size_alpha: Extra radial margin for S^{m_α} factor.
            overlap_size_beta:  Extra radial margin for S^{m_β} factor.
            use_round_metric_contraction: If True, compute ||N||² via the
                round product metric instead of the flat Frobenius norm.
        """
        super().__init__()
        self.total_dim = total_dim
        self.r_cut_alpha = RADIAL_MIDPOINT + overlap_size_alpha
        self.r_cut_beta  = RADIAL_MIDPOINT + overlap_size_beta
        self.overlap_size_alpha = overlap_size_alpha
        self.overlap_size_beta  = overlap_size_beta
        self.use_round_metric_contraction = use_round_metric_contraction
        # Reuse the per-patch Nijenhuis computation from the single-sphere loss.
        self._nij = NijenhuisTensorLoss(total_dim, overlap_size=0.0)

    def _radial_weight(
        self,
        rho_alpha: torch.Tensor,
        rho_beta:  torch.Tensor,
        hard_cutoff: bool,
    ) -> torch.Tensor:
        """Product radial weight w_α(ρ_α) · w_β(ρ_β)."""
        if hard_cutoff:
            wa = (rho_alpha <= self.r_cut_alpha).float()
            wb = (rho_beta  <= self.r_cut_beta ).float()
        else:
            exc_a = (rho_alpha - self.r_cut_alpha).clamp(min=0.0)
            exc_b = (rho_beta  - self.r_cut_beta ).clamp(min=0.0)
            wa = torch.exp(-(exc_a / self.overlap_size_alpha) ** 2)
            wb = torch.exp(-(exc_b / self.overlap_size_beta ) ** 2)
        return wa * wb

    def forward(
        self,
        x_alpha:   torch.Tensor,
        x_beta:    torch.Tensor,
        x_a2:      torch.Tensor,
        x_b2:      torch.Tensor,
        model:     nn.Module,
        hard_cutoff: bool = False,
    ) -> Tuple[torch.Tensor, Dict]:
        """
        Compute Nijenhuis loss for all four patches.

        Args:
            x_alpha: Patch-0 α-coords, shape (batch, m_α).
            x_beta:  Patch-0 β-coords, shape (batch, m_β).
            x_a2:    φ_α(x_alpha), shape (batch, m_α).
            x_b2:    φ_β(x_beta),  shape (batch, m_β).
            model:   FourPatchProductNet.
            hard_cutoff: If True use binary masks; intended for evaluation.

        Returns:
            loss: Scalar sum of per-patch weighted Nijenhuis norms.
            info: Dict with 'nij_{ia}{ib}' and 'nij_{ia}{ib}_max' entries.
        """
        patch_alpha = {0: x_alpha, 1: x_a2}
        patch_beta  = {0: x_beta,  1: x_b2}

        total = x_alpha.sum() * 0.0   # differentiable zero on correct device
        info: Dict = {}

        for ia, ib in [(0,0),(0,1),(1,0),(1,1)]:
            xa_p = patch_alpha[ia]
            xb_p = patch_beta[ib]
            rho_a = xa_p.detach().norm(dim=1)
            rho_b = xb_p.detach().norm(dim=1)
            w = self._radial_weight(rho_a, rho_b, hard_cutoff)

            subnet = model.patches[f'patch_{ia}{ib}']
            x_combined = torch.cat([xa_p, xb_p], dim=1)

            if self.use_round_metric_contraction:
                g, g_inv = SphereProduct.round_metric_ball_coords(xa_p, xb_p)
            else:
                g, g_inv = None, None

            # Re-evaluate with gradient tracking (required by NijenhuisTensorLoss)
            nij_l, nij_info = self._nij_with_weight(
                x_combined, subnet, w, hard_cutoff, g, g_inv,
            )

            key = f'nij_{ia}{ib}'
            info[key] = nij_l.item()
            total = total + nij_l

        return total, info

    @staticmethod
    def _nij_with_weight(
        x:           torch.Tensor,
        subnet:      nn.Module,
        w:           torch.Tensor,
        hard_cutoff: bool,
        g:           Optional[torch.Tensor] = None,
        g_inv:       Optional[torch.Tensor] = None,
    ) -> Tuple[torch.Tensor, Dict]:
        """
        Compute Nijenhuis norm for a single patch with an externally supplied weight.

        The weight already encodes the product radial structure; we bypass
        NijenhuisTensorLoss's internal weight and apply w directly.

        Args:
            x:      Combined ball coords, shape (batch, d).
            subnet: PatchSubNet.
            w:      Per-sample weight, shape (batch,).
            hard_cutoff: Passed through (not used here; handled in caller).
            g:      Optional round metric, shape (batch, d, d).  When provided
                    together with g_inv, the metric contraction
                    g_{pi} g^{qj} g^{rk} N^p_{qr} N^i_{jk} is used instead of
                    the flat Frobenius norm.
            g_inv:  Optional inverse metric, shape (batch, d, d).

        Returns:
            loss: Scalar weighted mean Nijenhuis norm.
            info: {'nij_norm': float}
        """
        with torch.enable_grad():
            J = subnet(x)          # (batch, d, d)
        batch, d, _ = J.shape

        x_d = x.detach()
        slices = []
        for l in range(d):
            tangent = x_d.new_zeros(batch, d)
            tangent[:, l] = 1.0
            with fwdAD.dual_level():
                dJ_l = fwdAD.unpack_dual(subnet(fwdAD.make_dual(x_d, tangent))).tangent
            if dJ_l is None:
                dJ_l = x_d.new_zeros(batch, d, d)
            slices.append(dJ_l)
        dJ = torch.stack(slices, dim=-1)    # (batch, d, d, d)

        term1 = torch.einsum('bli,bkjl->bkij', J, dJ)
        term2 = -torch.einsum('blj,bkil->bkij', J, dJ)
        anti  = dJ.permute(0, 1, 3, 2) - dJ
        term3 = torch.einsum('bkl,blij->bkij', J, anti)

        N = term1 + term2 - term3
        if g is not None and g_inv is not None:
            nij_norm = torch.sqrt(_metric_norm_sq_12tensor(N, g, g_inv) + 1e-20)
        else:
            nij_norm = torch.sqrt((N ** 2).sum(dim=(1, 2, 3)) + 1e-20)

        w_sum = w.sum().clamp(min=1e-8)
        loss = (w * nij_norm).sum() / w_sum
        return loss, {'nij_norm': loss.item()}


class ProductOverlapLoss(nn.Module):
    """
    Cross-patch consistency loss for S^{m_α} × S^{m_β}.

    Enforces J^{src}(x) = F · J^{tgt}(x') · F⁻¹ (forward) and the
    symmetric backward direction for all 6 patch pairs.  The block-diagonal
    structure of F is exploited via SphereProduct.block_jacobian_transform,
    which inverts F_α and F_β separately rather than the full (m_α+m_β)-square
    system.
    """

    def __init__(
        self,
        dim_alpha: int,
        dim_beta:  int,
        overlap_size_alpha: float = 0.1,
        overlap_size_beta:  float = 0.1,
        use_round_metric_contraction: bool = False,
    ):
        """
        Args:
            dim_alpha:          m_α.
            dim_beta:           m_β.
            overlap_size_alpha: Equatorial band width for α factor.
            overlap_size_beta:  Equatorial band width for β factor.
            use_round_metric_contraction: If True, compute ||Jdiff||² via the
                round product metric instead of the flat Frobenius norm.  The
                metric is evaluated in the patch where the comparison is
                performed (target for the forward diff, source for the backward).
        """
        super().__init__()
        self.dim_alpha = dim_alpha
        self.dim_beta  = dim_beta
        self.use_round_metric_contraction = use_round_metric_contraction

        r_outer_a = RADIAL_MIDPOINT + overlap_size_alpha
        r_outer_b = RADIAL_MIDPOINT + overlap_size_beta
        r_inner_a = (1.0 - r_outer_a) / (1.0 + r_outer_a)
        r_inner_b = (1.0 - r_outer_b) / (1.0 + r_outer_b)

        self.r_outer_alpha  = r_outer_a
        self.r_outer_beta   = r_outer_b
        self.r_inner_alpha  = r_inner_a
        self.r_inner_beta   = r_inner_b
        self.r_floor_alpha  = max(0.05, r_inner_a - 2.0 * overlap_size_alpha)
        self.r_floor_beta   = max(0.05, r_inner_b - 2.0 * overlap_size_beta)
        self.overlap_size_alpha = overlap_size_alpha
        self.overlap_size_beta  = overlap_size_beta

    def _equatorial_weight(
        self,
        rho_alpha: Optional[torch.Tensor],
        rho_beta:  Optional[torch.Tensor],
        hard_cutoff: bool,
    ) -> torch.Tensor:
        """
        Compute per-sample equatorial weight for the active factor(s).

        Args:
            rho_alpha: ‖x_α‖ per sample, or None if α is not being flipped.
            rho_beta:  ‖x_β‖ per sample, or None if β is not being flipped.
            hard_cutoff: Binary mask (eval) vs soft Gaussian (training).

        Returns:
            w: shape (batch,).
        """
        parts = []
        if rho_alpha is not None:
            if hard_cutoff:
                parts.append(
                    ((rho_alpha >= self.r_inner_alpha) & (rho_alpha <= self.r_outer_alpha)).float()
                )
            else:
                d_out = (rho_alpha - self.r_outer_alpha).clamp(min=0.0)
                d_in  = (self.r_inner_alpha - rho_alpha).clamp(min=0.0)
                parts.append(torch.exp(-((d_out + d_in) / self.overlap_size_alpha) ** 2))
        if rho_beta is not None:
            if hard_cutoff:
                parts.append(
                    ((rho_beta >= self.r_inner_beta) & (rho_beta <= self.r_outer_beta)).float()
                )
            else:
                d_out = (rho_beta - self.r_outer_beta).clamp(min=0.0)
                d_in  = (self.r_inner_beta - rho_beta).clamp(min=0.0)
                parts.append(torch.exp(-((d_out + d_in) / self.overlap_size_beta) ** 2))

        w = parts[0]
        for p in parts[1:]:
            w = w * p
        return w

    def forward(
        self,
        J_dict:  Dict,
        x_alpha: torch.Tensor,
        x_beta:  torch.Tensor,
        x_a2:    torch.Tensor,
        x_b2:    torch.Tensor,
        hard_cutoff: bool = False,
    ) -> Tuple[torch.Tensor, Dict]:
        """
        Compute all 6 pairwise overlap consistency losses.

        Args:
            J_dict:  {(i_α,i_β) → J tensor (batch, d, d)}.
            x_alpha: Patch-0 α-coords (batch, m_α).
            x_beta:  Patch-0 β-coords (batch, m_β).
            x_a2:    φ_α(x_alpha) (batch, m_α).
            x_b2:    φ_β(x_beta) (batch, m_β).
            hard_cutoff: If True use binary masks for evaluation.

        Returns:
            loss: Scalar mean of all 6 pair losses.
            info: Dict with 'overlap_{src}_{tgt}' entries and 'overlap_norm'.
        """
        batch = x_alpha.shape[0]
        zero  = x_alpha.sum() * 0.0
        ma = self.dim_alpha

        # Precompute Jacobians at both patch-0 and patch-1 coords for each factor.
        # Each is used by exactly 3 pairs (2 × α-flip + 1 × both-flip etc.).
        Fa,  Fa_inv  = _jac_and_inv(x_alpha)   # α at patch-0 coords
        Fa2, Fa2_inv = _jac_and_inv(x_a2)      # α at patch-1 coords
        Fb,  Fb_inv  = _jac_and_inv(x_beta)    # β at patch-0 coords
        Fb2, Fb2_inv = _jac_and_inv(x_b2)      # β at patch-1 coords

        # Hard floors: exclude near-singular Jacobians
        floor_a  = x_alpha.norm(dim=1) >= self.r_floor_alpha
        floor_a2 = x_a2.norm(dim=1)   >= self.r_floor_alpha
        floor_b  = x_beta.norm(dim=1)  >= self.r_floor_beta
        floor_b2 = x_b2.norm(dim=1)   >= self.r_floor_beta

        patch_alpha = {0: x_alpha, 1: x_a2}
        patch_beta  = {0: x_beta,  1: x_b2}
        Fa_jac  = {0: (Fa,  Fa_inv),  1: (Fa2, Fa2_inv)}
        Fb_jac  = {0: (Fb,  Fb_inv),  1: (Fb2, Fb2_inv)}
        floor_a_map = {0: floor_a, 1: floor_a2}
        floor_b_map = {0: floor_b, 1: floor_b2}

        total_loss = zero
        info: Dict = {}
        pair_count = 0

        for (ia_s, ib_s), (ia_t, ib_t), flip_type in _OVERLAP_PAIRS:
            # Source patch coords
            xa_s = patch_alpha[ia_s]
            xb_s = patch_beta [ib_s]
            rho_a_s = xa_s.norm(dim=1)
            rho_b_s = xb_s.norm(dim=1)

            # Floor mask (apply on whichever coords are used)
            if flip_type == 'alpha':
                mask = floor_a_map[ia_s]
                rho_w_a, rho_w_b = rho_a_s, None
            elif flip_type == 'beta':
                mask = floor_b_map[ib_s]
                rho_w_a, rho_w_b = None, rho_b_s
            else:  # both
                mask = floor_a_map[ia_s] & floor_b_map[ib_s]
                rho_w_a, rho_w_b = rho_a_s, rho_b_s

            if mask.sum() == 0:
                pair_key = f'overlap_{ia_s}{ib_s}_{ia_t}{ib_t}'
                info[pair_key] = 0.0
                continue

            # Slice to valid samples
            J_src = J_dict[(ia_s, ib_s)][mask]
            J_tgt = J_dict[(ia_t, ib_t)][mask]

            # Build the transition Jacobian for this pair.
            # The Jacobian F maps source patch coords to target patch coords;
            # it is evaluated at the *source* patch coordinates.
            if flip_type == 'beta':
                # Only β is flipped: F = block(I_α, F_β at source β-coords)
                Fb_s, Fb_s_inv = Fb_jac[ib_s]
                Fb_sm, Fb_sm_inv = Fb_s[mask], Fb_s_inv[mask]
                Fa_id = torch.eye(self.dim_alpha, dtype=x_alpha.dtype, device=x_alpha.device)
                Fa_id = Fa_id.unsqueeze(0).expand(mask.sum(), -1, -1)
                J_src_to_tgt = SphereProduct.block_jacobian_transform(
                    Fa_id, Fa_id, Fb_sm, Fb_sm_inv, J_src
                )
                J_tgt_to_src = SphereProduct.block_jacobian_transform(
                    Fa_id, Fa_id,
                    Fb_jac[ib_t][0][mask], Fb_jac[ib_t][1][mask],
                    J_tgt
                )
                rho_w_b_m = rho_b_s[mask]
                w = self._equatorial_weight(None, rho_w_b_m, hard_cutoff)

            elif flip_type == 'alpha':
                # Only α is flipped: F = block(F_α at source α-coords, I_β)
                Fa_s, Fa_s_inv = Fa_jac[ia_s]
                Fa_sm, Fa_sm_inv = Fa_s[mask], Fa_s_inv[mask]
                Fb_id = torch.eye(self.dim_beta, dtype=x_beta.dtype, device=x_beta.device)
                Fb_id = Fb_id.unsqueeze(0).expand(mask.sum(), -1, -1)
                J_src_to_tgt = SphereProduct.block_jacobian_transform(
                    Fa_sm, Fa_sm_inv, Fb_id, Fb_id, J_src
                )
                J_tgt_to_src = SphereProduct.block_jacobian_transform(
                    Fa_jac[ia_t][0][mask], Fa_jac[ia_t][1][mask],
                    Fb_id, Fb_id,
                    J_tgt
                )
                rho_w_a_m = rho_a_s[mask]
                w = self._equatorial_weight(rho_w_a_m, None, hard_cutoff)

            else:  # both flipped
                # F = block(F_α at source α-coords, F_β at source β-coords)
                Fa_s, Fa_s_inv = Fa_jac[ia_s]
                Fb_s, Fb_s_inv = Fb_jac[ib_s]
                Fa_sm,  Fa_sm_inv  = Fa_s[mask],  Fa_s_inv[mask]
                Fb_sm,  Fb_sm_inv  = Fb_s[mask],  Fb_s_inv[mask]
                J_src_to_tgt = SphereProduct.block_jacobian_transform(
                    Fa_sm, Fa_sm_inv, Fb_sm, Fb_sm_inv, J_src
                )
                # Backward: F evaluated at target patch coords
                Fa_t, Fa_t_inv = Fa_jac[ia_t]
                Fb_t, Fb_t_inv = Fb_jac[ib_t]
                J_tgt_to_src = SphereProduct.block_jacobian_transform(
                    Fa_t[mask], Fa_t_inv[mask], Fb_t[mask], Fb_t_inv[mask],
                    J_tgt
                )
                rho_w_a_m = rho_a_s[mask]
                rho_w_b_m = rho_b_s[mask]
                w = self._equatorial_weight(rho_w_a_m, rho_w_b_m, hard_cutoff)

            w_sum = w.sum().clamp(min=1e-8)
            if hard_cutoff and w_sum < 1.0:
                pair_key = f'overlap_{ia_s}{ib_s}_{ia_t}{ib_t}'
                info[pair_key] = 0.0
                continue

            diff_fwd = J_src_to_tgt - J_tgt
            diff_bwd = J_tgt_to_src - J_src

            if self.use_round_metric_contraction:
                # Metric in target patch for diff_fwd, source patch for diff_bwd.
                xa_tgt = patch_alpha[ia_t][mask]
                xb_tgt = patch_beta [ib_t][mask]
                xa_src = patch_alpha[ia_s][mask]
                xb_src = patch_beta [ib_s][mask]
                g_tgt, g_tgt_inv = SphereProduct.round_metric_ball_coords(xa_tgt, xb_tgt)
                g_src, g_src_inv = SphereProduct.round_metric_ball_coords(xa_src, xb_src)
                norm_fwd = torch.sqrt(_metric_norm_sq_11tensor(diff_fwd, g_tgt, g_tgt_inv) + 1e-20)
                norm_bwd = torch.sqrt(_metric_norm_sq_11tensor(diff_bwd, g_src, g_src_inv) + 1e-20)
            else:
                norm_fwd = torch.sqrt((diff_fwd ** 2).sum(dim=(1, 2)) + 1e-20)
                norm_bwd = torch.sqrt((diff_bwd ** 2).sum(dim=(1, 2)) + 1e-20)
            pair_loss = ((w * norm_fwd).sum() + (w * norm_bwd).sum()) / (2.0 * w_sum)

            pair_key = f'overlap_{ia_s}{ib_s}_{ia_t}{ib_t}'
            info[pair_key] = pair_loss.item()
            total_loss = total_loss + pair_loss
            pair_count += 1

        if pair_count > 0:
            total_loss = total_loss / pair_count

        info['overlap_norm'] = total_loss.item()
        return total_loss, info


def _jac_and_inv(x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
    """Compute patch-change Jacobian and its inverse for a single factor."""
    F = EvenSphere.patch_change_jacobian(x)
    return F, torch.linalg.inv(F)


class ProductTotalLoss(nn.Module):
    """
    Combined Nijenhuis + overlap loss for S^{m_α} × S^{m_β}.

    J² = −I is enforced exactly by the exp(A) J₀ exp(−A) parametrisation;
    no constraint loss is required.  Smoothness follows from the network being
    smooth in smooth chart coordinates.
    """

    def __init__(
        self,
        dim_alpha: int,
        dim_beta:  int,
        weights:   Dict,
    ):
        """
        Args:
            dim_alpha: m_α.
            dim_beta:  m_β.
            weights:   Dict with keys 'nijenhuis', 'overlap', 'overlap_size_alpha',
                       'overlap_size_beta', 'use_round_metric_contraction'.
        """
        super().__init__()
        osa  = float(weights.get('overlap_size_alpha', 0.1))
        osb  = float(weights.get('overlap_size_beta',  0.1))
        urmc = bool(weights.get('use_round_metric_contraction', False))

        self.nij_loss     = ProductNijenhuisLoss(
            dim_alpha + dim_beta, osa, osb,
            use_round_metric_contraction=urmc,
        )
        self.overlap_loss = ProductOverlapLoss(
            dim_alpha, dim_beta, osa, osb,
            use_round_metric_contraction=urmc,
        )
        self.weights      = weights

    def forward(
        self,
        J_dict:  Dict,
        x_alpha: torch.Tensor,
        x_beta:  torch.Tensor,
        x_a2:    torch.Tensor,
        x_b2:    torch.Tensor,
        model:   nn.Module,
        hard_cutoff: bool = False,
        report_all:  bool = False,
    ) -> Tuple[torch.Tensor, Dict]:
        """
        Compute total weighted loss.

        Args:
            J_dict:      {(i_α,i_β) → J (batch, d, d)}.
            x_alpha:     Patch-0 α-coords (batch, m_α).
            x_beta:      Patch-0 β-coords (batch, m_β).
            x_a2:        φ_α(x_alpha) (batch, m_α).
            x_b2:        φ_β(x_beta) (batch, m_β).
            model:       FourPatchProductNet.
            hard_cutoff: Binary masks for evaluation.
            report_all:  Compute all components regardless of weight.

        Returns:
            total_loss: Weighted sum.
            info:       Dict with all scalar components.
        """
        w_nij     = float(self.weights.get('nijenhuis', 1.0))
        w_overlap = float(self.weights.get('overlap',   1.0))

        zero = x_alpha.sum() * 0.0

        if w_nij != 0.0 or report_all:
            nij_loss, nij_info = self.nij_loss(
                x_alpha, x_beta, x_a2, x_b2, model, hard_cutoff
            )
        else:
            nij_loss = zero
            nij_info = {f'nij_{ia}{ib}': 0.0 for ia, ib in [(0,0),(0,1),(1,0),(1,1)]}
            for k in list(nij_info.keys()):
                nij_info[k + '_max'] = 0.0

        if w_overlap != 0.0 or report_all:
            overlap_loss, overlap_info = self.overlap_loss(
                J_dict, x_alpha, x_beta, x_a2, x_b2, hard_cutoff
            )
        else:
            overlap_loss = zero
            overlap_info = {'overlap_norm': 0.0}

        total = w_nij * nij_loss / 4.0 + w_overlap * overlap_loss

        info: Dict = {'total_loss': total.item()}
        info.update(nij_info)
        info.update(overlap_info)
        return total, info
