"""
Nijenhuis tensor computation and training losses for the two-patch ball formulation.

In local ball coordinates x ∈ B(0,1) ⊂ ℝ^{2n} an almost complex structure is a
smooth map J: B → End(ℝ^{2n}) satisfying J² = −I.  Because the network is smooth
and the chart is smooth, J is smooth by construction; no smoothness penalty is needed.

Nijenhuis tensor in chart coordinates (coordinate-expression):

    N^k_{ij}(x) = J^l_i ∂_l J^k_j − J^l_j ∂_l J^k_i − J^k_l (∂_i J^l_j − ∂_j J^l_i)

All partial derivatives ∂_l J^k_j are computed by automatic differentiation.

Both losses use smooth radial weighting controlled by `overlap_size` (δ):

  Nijenhuis weight:  w_N(ρ) = exp(−max(0, ρ − (r_eq+δ))² / δ²)
    near-constant in [0, r_eq+δ], Gaussian decay beyond.

  Overlap weight:    w_O(ρ) = exp(−(d_outer + d_inner)² / δ²)
    where d_outer = max(0, ρ−r_outer), d_inner = max(0, r_inner−ρ),
    r_outer = r_eq+δ,  r_inner = (1−r_outer)/(1+r_outer).
    r_inner is the exact patch-change preimage of r_outer, making the plateau
    [r_inner, r_outer] geometrically symmetric about the equator.

where F is the Jacobian of the self-inverse patch-change map and x₂ = φ(x₁).
"""

import torch
import torch.autograd.forward_ad as fwdAD
import torch.nn as nn
from typing import Dict, Tuple

from src.geometry import EvenSphere, RADIAL_MIDPOINT


class NijenhuisTensorLoss(nn.Module):
    """
    Per-patch Nijenhuis tensor loss in ball coordinates.

    Re-evaluates the sub-network with gradient-tracked inputs so that
    ∂J/∂x can be obtained via autograd, then assembles the full
    coordinate-Nijenhuis tensor:

        N^k_{ij} = J^l_i ∂_l J^k_j − J^l_j ∂_l J^k_i − J^k_l(∂_i J^l_j − ∂_j J^l_i)

    The per-sample Nijenhuis norms are combined via a radially weighted mean.
    The weight is near-constant in [0, r_eq + overlap_size] (the patch's natural
    hemisphere plus a margin) and decays as a Gaussian beyond that, so the near-pole
    region—where the other patch has better coverage—contributes less.
    """

    def __init__(self, manifold_dim: int, overlap_size: float = 0.1):
        """
        Args:
            manifold_dim: Real dimension of the manifold.
            overlap_size: Extra radial margin beyond the equator (r_eq = √2−1)
                defining the flat plateau; controls where the Gaussian decay begins.
        """
        super().__init__()
        self.manifold_dim = manifold_dim
        self.overlap_size = overlap_size

    def forward(
        self,
        x: torch.Tensor,
        subnet: nn.Module,
        hard_cutoff: bool = False,
    ) -> Tuple[torch.Tensor, Dict]:
        """
        Compute Nijenhuis loss for a single patch.

        Args:
            x: Ball coordinates, shape (batch, manifold_dim).
            subnet: PatchSubNet to differentiate through.
            hard_cutoff: If True, use a binary mask (ρ ≤ r_cut) instead of
                the soft Gaussian weight. Intended for evaluation only.

        Returns:
            loss: Scalar mean Nijenhuis norm.
            info: Dict with 'nij_norm' (mean).
        """
        with torch.enable_grad():
            J = subnet(x)              # (batch, d, d)
        batch, d, _ = J.shape

        # dJ[b, k, j, l] = ∂_{x^l} J^k_j
        # Computed via d forward-mode AD passes (one per input coordinate l).
        # Forward-mode has cost O(d · forward_pass) vs O(d² · backward_pass)
        # for the reverse-mode loop it replaces, giving ~d/2 speedup.
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
        dJ = torch.stack(slices, dim=-1)            # (batch, d, d, d)

        # term1[b,k,i,j] = Σ_l J^l_i * ∂_l J^k_j = Σ_l J[b,l,i] * dJ[b,k,j,l]
        term1 = torch.einsum('bli,bkjl->bkij', J, dJ)

        # term2[b,k,i,j] = −Σ_l J^l_j * ∂_l J^k_i = −Σ_l J[b,l,j] * dJ[b,k,i,l]
        term2 = -torch.einsum('blj,bkil->bkij', J, dJ)

        # anti[b,l,i,j] = ∂_i J^l_j − ∂_j J^l_i = dJ[b,l,j,i] − dJ[b,l,i,j]
        #   dJ.permute(0,1,3,2)[b,l,i,j] = dJ[b,l,j,i]   ← swaps axes 2↔3
        anti = dJ.permute(0, 1, 3, 2) - dJ

        # term3[b,k,i,j] = Σ_l J^k_l * (∂_i J^l_j − ∂_j J^l_i) = Σ_l J[b,k,l] * anti[b,l,i,j]
        term3 = torch.einsum('bkl,blij->bkij', J, anti)

        N = term1 + term2 - term3           # (batch, k, i, j)
        nij_sq   = (N ** 2).sum(dim=(1, 2, 3))
        nij_norm = torch.sqrt(nij_sq + 1e-20)   # (batch,)

        # Radial weight: plateau in [0, r_eq + overlap_size].  During training a
        # soft Gaussian decay is used beyond the cutoff; during evaluation a hard
        # binary mask is used so the metric is purely geometric.
        rho   = x.detach().norm(dim=1)                          # (batch,)
        r_cut = RADIAL_MIDPOINT + self.overlap_size
        if hard_cutoff:
            w = (rho <= r_cut).float()
        else:
            excess = (rho - r_cut).clamp(min=0.0)
            w = torch.exp(-(excess / self.overlap_size) ** 2)   # (batch,)

        w_sum = w.sum().clamp(min=1e-8)
        loss = (w * nij_norm).sum() / w_sum
        return loss, {'nij_norm': loss.item()}


class OverlapLoss(nn.Module):
    """
    Cross-patch consistency loss with equatorial radial weighting.

    For a point x₁ in patch-1 ball coordinates with x₂ = φ(x₁) the
    corresponding patch-2 coordinates, the J tensors must satisfy:

        J₂(x₂) = F(x₁) · J₁(x₁) · F(x₁)⁻¹

    and symmetrically:

        J₁(x₁) = F(x₂) · J₂(x₂) · F(x₂)⁻¹

    Near either pole the patch-change Jacobian F has condition number
    ∝ 1/ρ², so contributions there are down-weighted by a soft radial
    weight that is near-constant within the equatorial band and decays
    as a Gaussian outside it.  A hard floor at ρ ≥ r_floor prevents
    evaluation at nearly-singular Jacobians.

    The equatorial band is defined by:
        r_inner = (1 − r_outer) / (1 + r_outer),  with r_outer = r_eq + overlap_size.
    r_inner is the exact patch-change preimage of r_outer, so the band
    is geometrically symmetric about the equator (r_eq = √2 − 1).
    """

    def __init__(self, overlap_size: float = 0.1):
        """
        Args:
            overlap_size: Extra radial margin beyond the equator (r_eq = √2−1).
                r_outer = r_eq + overlap_size;
                r_inner = (1 − r_outer) / (1 + r_outer).
        """
        super().__init__()
        self.overlap_size = overlap_size
        r_outer = RADIAL_MIDPOINT + overlap_size
        r_inner = (1.0 - r_outer) / (1.0 + r_outer)
        self.r_outer = r_outer
        self.r_inner = r_inner
        # Hard safety floor: exclude near-singular Jacobians without relying
        # entirely on the soft weight.
        self.r_floor = max(0.05, r_inner - 2.0 * overlap_size)

    def forward(
        self,
        J1: torch.Tensor,
        J2: torch.Tensor,
        x1: torch.Tensor,
        x2: torch.Tensor,
        hard_cutoff: bool = False,
    ) -> Tuple[torch.Tensor, Dict]:
        """
        Compute radially-weighted overlap consistency loss.

        Args:
            J1: J in patch-1 coordinates, shape (batch, dim, dim).
            J2: J in patch-2 coordinates, shape (batch, dim, dim).
            x1: Patch-1 ball coords, shape (batch, dim).
            x2: Patch-2 ball coords, shape (batch, dim).
            hard_cutoff: If True, use a binary mask (ρ ∈ [r_inner, r_outer])
                instead of the soft Gaussian weight. Intended for evaluation only.

        Returns:
            loss: Scalar overlap loss (0 if no valid samples).
            info: Dict with 'overlap_norm'.
        """
        rho1 = x1.norm(dim=1)   # (batch,)
        rho2 = x2.norm(dim=1)   # (batch,)

        # Hard floor to avoid near-singular F inversions.
        mask = (rho1 >= self.r_floor) & (rho2 >= self.r_floor)
        if mask.sum() == 0:
            zero = J1.sum() * 0.0
            return zero, {'overlap_norm': 0.0}

        J1_m  = J1[mask]
        J2_m  = J2[mask]
        x1_m  = x1[mask]
        x2_m  = x2[mask]
        rho1_m = rho1[mask]

        # Equatorial weight: plateau in [r_inner, r_outer].  Training uses a soft
        # Gaussian decay outside; evaluation uses a hard binary mask.
        if hard_cutoff:
            w = ((rho1_m >= self.r_inner) & (rho1_m <= self.r_outer)).float()
        else:
            d_outer = (rho1_m - self.r_outer).clamp(min=0.0)
            d_inner = (self.r_inner - rho1_m).clamp(min=0.0)
            w = torch.exp(-((d_outer + d_inner) / self.overlap_size) ** 2)  # (batch_m,)

        # Transform J1 from patch-1 to patch-2 frame
        J1_to_2 = EvenSphere.patch_change_J(x1_m, J1_m)   # F(x1) J1 F(x1)^{-1}

        # Transform J2 from patch-2 to patch-1 frame
        J2_to_1 = EvenSphere.patch_change_J(x2_m, J2_m)   # F(x2) J2 F(x2)^{-1}

        diff1 = J1_to_2 - J2_m
        diff2 = J2_to_1 - J1_m

        norm1 = torch.sqrt((diff1 ** 2).sum(dim=(1, 2)) + 1e-20)
        norm2 = torch.sqrt((diff2 ** 2).sum(dim=(1, 2)) + 1e-20)

        w_sum = w.sum().clamp(min=1e-8)
        if hard_cutoff and w_sum < 1.0:
            zero = J1.sum() * 0.0
            return zero, {'overlap_norm': 0.0}
        loss  = ((w * norm1).sum() + (w * norm2).sum()) / (2.0 * w_sum)

        info = {
            'overlap_norm': loss.item(),
        }
        return loss, info


class TotalLoss(nn.Module):
    """
    Combined training loss: per-patch Nijenhuis + overlap consistency.

    J² = −I is satisfied exactly by the exp(A) J₀ exp(−A) parametrisation.
    Smoothness of J follows from the smoothness of the network in smooth
    chart coordinates; no smoothness penalty is required.
    """

    def __init__(self, manifold_dim: int, weights: Dict[str, float]):
        """
        Args:
            manifold_dim: Real dimension of the manifold.
            weights: Loss weights; expected keys: 'nijenhuis', 'overlap', 'overlap_size'.
        """
        super().__init__()
        overlap_size = float(weights.get('overlap_size', 0.1))
        self.nij_loss     = NijenhuisTensorLoss(manifold_dim, overlap_size=overlap_size)
        self.overlap_loss = OverlapLoss(overlap_size=overlap_size)
        self.weights      = weights

    def forward(
        self,
        J1:    torch.Tensor,
        J2:    torch.Tensor,
        x1:    torch.Tensor,
        x2:    torch.Tensor,
        model: nn.Module,
        hard_cutoff: bool = False,
        report_all: bool = False,
    ) -> Tuple[torch.Tensor, Dict]:
        """
        Compute total weighted loss.

        Args:
            J1: J in patch 1, shape (batch, dim, dim).
            J2: J in patch 2, shape (batch, dim, dim).
            x1: Patch-1 ball coords, shape (batch, dim).
            x2: Patch-2 ball coords, shape (batch, dim).
            model: TwoPatchAlmostComplexNet instance.
            hard_cutoff: If True, apply binary geometric masks instead of
                soft Gaussian weights. Use for validation and test evaluation.
            report_all: If True, compute all components regardless of their
                weight (for test evaluation). Weights still scale total_loss.

        Returns:
            total_loss: Weighted sum of all loss terms.
            info: Dict with all scalar components.
        """
        w_nij     = self.weights.get('nijenhuis', 1.0)
        w_overlap = self.weights.get('overlap', 1.0)

        zero = J1.sum() * 0.0

        if w_nij != 0.0 or report_all:
            nij1, info1 = self.nij_loss(x1, model.patch1, hard_cutoff=hard_cutoff)
            nij2, info2 = self.nij_loss(x2, model.patch2, hard_cutoff=hard_cutoff)
        else:
            nij1 = nij2 = zero
            info1 = info2 = {'nij_norm': 0.0}

        if w_overlap != 0.0 or report_all:
            overlap, info_ov = self.overlap_loss(J1, J2, x1, x2, hard_cutoff=hard_cutoff)
        else:
            overlap = zero
            info_ov = {'overlap_norm': 0.0}

        total = w_nij * (nij1 + nij2) / 2.0 + w_overlap * overlap

        info = {
            'total_loss':   total.item(),
            'nij1':         nij1.item(),
            'nij2':         nij2.item(),
            'overlap_norm': info_ov['overlap_norm'],
        }

        return total, info


if __name__ == "__main__":
    print("Testing loss functions...")

    from src.models import TwoPatchAlmostComplexNet

    dim       = 4
    batch     = 8
    sphere    = EvenSphere(dimension=dim)

    net = TwoPatchAlmostComplexNet(
        manifold_dim=dim,
        hidden_layers=[32, 64, 32],
        activation='tanh',
    )

    weights = {'nijenhuis': 1.0, 'overlap': 1.0}
    loss_fn = TotalLoss(manifold_dim=dim, weights=weights)

    x1 = sphere.sample_points_patch1(batch)
    J1, J2, x1, x2 = net(x1)

    loss, info = loss_fn(J1, J2, x1, x2, net)
    print(f"Total loss: {info['total_loss']:.4f}")
    print(f"Nij1: {info['nij1']:.4f}  Nij2: {info['nij2']:.4f}")
    print(f"Overlap: {info['overlap_norm']:.4f}")
    assert loss.item() >= 0
    print("Loss function test passed!")
