"""
Geometry utilities for products of spheres S^{m_α} × S^{m_β}.

Each factor sphere is covered by two open ball patches derived from the
modified stereographic projection (see sphere.py).  The product is covered by
2 × 2 = 4 patches indexed (i_α, i_β) ∈ {0,1}².

The patch-change map on the product decomposes as a direct product:

    φ_{α×β}(x_α, x_β) = (φ_α(x_α), φ_β(x_β))

so the Jacobian is exactly block-diagonal:

    F_{α×β} = block(F_α, F_β),    F_{α×β}⁻¹ = block(F_α⁻¹, F_β⁻¹)

This exact block structure is exploited throughout for efficiency and
numerical stability: inverting (and multiplying by) blocks of size m_α and
m_β separately instead of a dense (m_α + m_β) × (m_α + m_β) system.

The total real dimension d = m_α + m_β must be even for an almost complex
structure to exist on the product.  Neither factor is required to be even
individually (e.g. S³ × S³ is allowed).
"""

import torch
import numpy as np
from typing import Dict, Tuple

from .sphere import EvenSphere, RADIAL_MIDPOINT, ball_sample


class SphereProduct:
    """
    Product manifold S^{m_α} × S^{m_β} with four open ball patch charts.

    Each factor sphere contributes two stereographic ball patches; the product
    gives four patches (i_α, i_β) ∈ {0,1}².
    """

    def __init__(self, dim_alpha: int, dim_beta: int, device: str = 'cpu'):
        """
        Args:
            dim_alpha: Real dimension of S^{m_α} (≥ 1).
            dim_beta:  Real dimension of S^{m_β} (≥ 1).
            device:    torch device string.

        Raises:
            ValueError: if dim_alpha < 1, dim_beta < 1, or
                        (dim_alpha + dim_beta) is odd.
        """
        if dim_alpha < 1:
            raise ValueError(f"dim_alpha must be ≥ 1, got {dim_alpha}")
        if dim_beta < 1:
            raise ValueError(f"dim_beta must be ≥ 1, got {dim_beta}")
        if (dim_alpha + dim_beta) % 2 != 0:
            raise ValueError(
                f"dim_alpha + dim_beta must be even for an almost complex structure to exist; "
                f"got {dim_alpha} + {dim_beta} = {dim_alpha + dim_beta} (odd)"
            )

        self.dim_alpha  = dim_alpha
        self.dim_beta   = dim_beta
        self.total_dim  = dim_alpha + dim_beta
        self.complex_dim = self.total_dim // 2
        self.device     = device

    # ------------------------------------------------------------------
    # Sampling
    # ------------------------------------------------------------------

    def sample_points(self, num_points: int) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Sample points on S^{m_α} × S^{m_β} using independent product measure.

        Args:
            num_points: Number of sample points.

        Returns:
            x_alpha: Patch-0 ball coords on S^{m_α}, shape (num_points, dim_alpha).
            x_beta:  Patch-0 ball coords on S^{m_β}, shape (num_points, dim_beta).
        """
        xa = torch.tensor(
            ball_sample(num_points, self.dim_alpha), dtype=torch.float32, device=self.device
        )
        xb = torch.tensor(
            ball_sample(num_points, self.dim_beta), dtype=torch.float32, device=self.device
        )
        return xa, xb

    # ------------------------------------------------------------------
    # Patch coordinates
    # ------------------------------------------------------------------

    def get_all_patch_coords(
        self,
        x_alpha: torch.Tensor,
        x_beta: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Compute all four patch coordinate pairs from patch-(0,0) inputs.

        φ_α and φ_β are computed once and reused across all four patches.

        Args:
            x_alpha: Patch-0 α-coords, shape (batch, dim_alpha).
            x_beta:  Patch-0 β-coords, shape (batch, dim_beta).

        Returns:
            x00: concat of (x_α,    x_β),    shape (batch, total_dim).
            x01: concat of (x_α,    φ_β(x_β)), shape (batch, total_dim).
            x10: concat of (φ_α(x_α), x_β),    shape (batch, total_dim).
            x11: concat of (φ_α(x_α), φ_β(x_β)), shape (batch, total_dim).
        """
        xa2 = EvenSphere.patch_change_coords(x_alpha)
        xb2 = EvenSphere.patch_change_coords(x_beta)
        x00 = torch.cat([x_alpha, x_beta], dim=1)
        x01 = torch.cat([x_alpha, xb2],    dim=1)
        x10 = torch.cat([xa2,     x_beta], dim=1)
        x11 = torch.cat([xa2,     xb2],    dim=1)
        return x00, x01, x10, x11

    # ------------------------------------------------------------------
    # Block-diagonal Jacobian transform
    # ------------------------------------------------------------------

    @staticmethod
    def block_jacobian_transform(
        F_alpha:     torch.Tensor,
        F_alpha_inv: torch.Tensor,
        F_beta:      torch.Tensor,
        F_beta_inv:  torch.Tensor,
        J:           torch.Tensor,
    ) -> torch.Tensor:
        """
        Apply the block-diagonal patch-change Jacobian to a (1,1)-tensor J.

        Because the product patch-change map decomposes as
            φ_{α×β}(x_α, x_β) = (φ_α(x_α), φ_β(x_β)),
        its Jacobian is exactly block-diagonal:
            F = block(F_α, F_β),   F⁻¹ = block(F_α⁻¹, F_β⁻¹)

        The transform J ↦ F J F⁻¹ is computed block-row/block-column wise,
        avoiding assembly of the full (m_α+m_β)² system:

            (F J)_αα = F_α J_αα,    (F J)_αβ = F_α J_αβ
            (F J)_βα = F_β J_βα,    (F J)_ββ = F_β J_ββ

            (F J F⁻¹)_αα = (F J)_αα F_α⁻¹,   (F J F⁻¹)_αβ = (F J)_αβ F_β⁻¹
            (F J F⁻¹)_βα = (F J)_βα F_α⁻¹,   (F J F⁻¹)_ββ = (F J)_ββ F_β⁻¹

        Savings vs. full inversion: O((m+m)³) → O(2m³), i.e. 4× for equal dims.

        Args:
            F_alpha:     α-block Jacobian, shape (batch, m_α, m_α).
            F_alpha_inv: Inverse of F_alpha, shape (batch, m_α, m_α).
            F_beta:      β-block Jacobian, shape (batch, m_β, m_β).
            F_beta_inv:  Inverse of F_beta, shape (batch, m_β, m_β).
            J:           (1,1)-tensor to transform, shape (batch, d, d),
                         d = m_α + m_β.

        Returns:
            F J F⁻¹, shape (batch, d, d).
        """
        ma = F_alpha.shape[1]
        mb = F_beta.shape[1]

        # Split J into 2×2 blocks along both axes
        Jaa = J[:, :ma, :ma]
        Jab = J[:, :ma, ma:]
        Jba = J[:, ma:, :ma]
        Jbb = J[:, ma:, ma:]

        # Left-multiply each row-block by its Jacobian factor
        FJaa = F_alpha @ Jaa
        FJab = F_alpha @ Jab
        FJba = F_beta  @ Jba
        FJbb = F_beta  @ Jbb

        # Right-multiply each column-block by the inverse Jacobian factor
        FJFaa = FJaa @ F_alpha_inv
        FJFab = FJab @ F_beta_inv
        FJFba = FJba @ F_alpha_inv
        FJFbb = FJbb @ F_beta_inv

        top    = torch.cat([FJFaa, FJFab], dim=2)
        bottom = torch.cat([FJFba, FJFbb], dim=2)
        return torch.cat([top, bottom], dim=1)

    @staticmethod
    def get_block_jacobians(
        x_alpha: torch.Tensor,
        x_beta:  torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Compute α and β block Jacobians and their inverses.

        Args:
            x_alpha: α-patch coords, shape (batch, m_α).
            x_beta:  β-patch coords, shape (batch, m_β).

        Returns:
            F_alpha, F_alpha_inv: shapes (batch, m_α, m_α).
            F_beta,  F_beta_inv:  shapes (batch, m_β, m_β).
        """
        F_alpha = EvenSphere.patch_change_jacobian(x_alpha)
        F_beta  = EvenSphere.patch_change_jacobian(x_beta)
        F_alpha_inv = torch.linalg.inv(F_alpha)
        F_beta_inv  = torch.linalg.inv(F_beta)
        return F_alpha, F_alpha_inv, F_beta, F_beta_inv
