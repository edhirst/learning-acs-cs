"""
Geometry utilities for even-dimensional spheres.

Represents S^{2n} via two open ball patches derived from the modified
stereographic projection.  The ball coordinate map is:

    x_i = X_i / (1 + sqrt(1 + |X|²))

where X ∈ ℝ^{2n} are standard stereographic coordinates (north-pole chart
for patch 1, south-pole chart for patch 2).  The image is the open ball
B(0, 1) ⊂ ℝ^{2n}; the overlap region sits around the radial midpoint
r_mid = √2 − 1.

The patch-change map (both directions, self-inverse):
    x' = x · (|x| − 1) / (|x|(|x| + 1))

which satisfies  |x'| = (1 − |x|) / (1 + |x|)  for |x| < 1.
"""

import torch
import numpy as np
from typing import Tuple


# Radial midpoint of the overlap region (image of the equator under each chart)
RADIAL_MIDPOINT: float = float(np.sqrt(2.0) - 1.0)


class EvenSphere:
    """
    Even-dimensional sphere S^{2n} with two open ball patch charts.

    Charts are derived from the modified stereographic projection; each
    sub-network receives 2n-dimensional ball coordinates as input.
    """
    
    def __init__(self, dimension: int, device: str = 'cpu'):
        """
        Args:
            dimension: Real manifold dimension (must be even, ≥ 2).
            device: torch device string.
        """
        if dimension % 2 != 0:
            raise ValueError(f"Dimension must be even, got {dimension}")
        if dimension < 2:
            raise ValueError(f"Dimension must be at least 2, got {dimension}")

        self.dim = dimension               # 2n
        self.ambient_dim = dimension + 1   # 2n + 1  (for reference / visualisation)
        self.complex_dim = dimension // 2  # n
        self.device = device

    # ------------------------------------------------------------------
    # Sampling
    # ------------------------------------------------------------------

    def sample_points_patch1(
        self,
        num_points: int,
        patch_width: float = 1.0,
        density_power: float = 1.0,
        overlap_size: float | None = None,
    ) -> torch.Tensor:
        """
        Sample points in patch-1 ball coordinates uniformly over S^{2n}.

        Half the samples have radii drawn from a beta distribution covering
        the patch-1 hemisphere; the other half are patch-2 beta samples
        mapped via r' = (1−r)/(1+r) to represent the patch-2 hemisphere in
        patch-1 coordinates.

        If overlap_size is given, all points are drawn from the equatorial
        annulus [r_inner, r_outer] (r_outer = r_eq + overlap_size;
        r_inner = (1−r_outer)/(1+r_outer)) via inverse-CDF sampling.  This
        concentrates training signal on the overlap region where the
        consistency loss is active while remaining geometrically valid:
        the annulus is self-symmetric under the patch-change map, so the
        patch-2 images φ(x) lie in the same annulus.

        Args:
            num_points:   Number of points to sample.
            patch_width:  Radial cutoff used only when overlap_size is None.
            density_power: Beta shape parameter used only when overlap_size is None.
            overlap_size: If set, sample from the equatorial annulus defined
                by this overlap margin.

        Returns:
            Ball coordinates in patch 1, shape (num_points, dim).
        """
        if overlap_size is not None:
            r_outer = RADIAL_MIDPOINT + overlap_size
            r_inner = (1.0 - r_outer) / (1.0 + r_outer)
            pts_np = equatorial_annulus_sample(num_points, self.dim, r_inner, r_outer)
        else:
            pts_np = ball_sample(num_points, self.dim, patch_width, density_power)
        return torch.tensor(pts_np, dtype=torch.float32, device=self.device)

    # ------------------------------------------------------------------
    # Patch geometry
    # ------------------------------------------------------------------

    @staticmethod
    def patch_change_coords(x: torch.Tensor) -> torch.Tensor:
        """
        Map ball coordinates from one patch to the other.

        The map is its own inverse:
            x' = x · (|x| − 1) / (|x|(|x| + 1))

        This follows from composing the two modified stereographic projections
        via X' = X / |X|² (inversion in stereographic coordinates).

        Args:
            x: Ball coordinates, shape (batch, dim).

        Returns:
            Ball coordinates in the other patch, shape (batch, dim).
        """
        rho = torch.norm(x, dim=1, keepdim=True).clamp(min=1e-8)
        scale = (rho - 1.0) / (rho * (rho + 1.0))
        return x * scale

    @staticmethod
    def patch_change_jacobian(x: torch.Tensor) -> torch.Tensor:
        """
        Compute the Jacobian  F_{ij} = ∂x'_i/∂x_j  of the patch-change map.

        F has two eigenvalues:
            λ_tan = (ρ−1)/(ρ(ρ+1))        (multiplicity dim−1, tangential)
            λ_rad = 2/(ρ+1)²              (multiplicity 1, radial; always > 0)

        Written in terms of the normalised direction x̂ = x/ρ:

            F = λ_tan · (I − x̂⊗x̂) + λ_rad · (x̂⊗x̂)

        This form avoids the catastrophic float32 cancellation that occurs in
        the equivalent φδ_{ij} + φ'x_ix_j/ρ formula when ρ ≪ 1 (where both
        terms become O(1/ρ) with opposite signs, cancelling to zero in float32
        even though the true value is ≈ 2).

        Args:
            x: Ball coordinates, shape (batch, dim).

        Returns:
            Jacobian matrices, shape (batch, dim, dim).
        """
        batch, dim = x.shape
        rho = torch.norm(x, dim=1, keepdim=True).clamp(min=1e-8, max=0.99)  # (batch, 1)

        lam_tan = (rho - 1.0) / (rho * (rho + 1.0))   # tangential eigenvalue φ(ρ)
        lam_rad = 2.0 / (rho + 1.0).pow(2)             # radial eigenvalue; always positive

        x_hat = x / rho                                        # (batch, dim)
        outer_hat = x_hat.unsqueeze(2) * x_hat.unsqueeze(1)   # x̂⊗x̂, (batch, dim, dim)
        eye = torch.eye(dim, dtype=x.dtype, device=x.device).unsqueeze(0)

        F = lam_tan.unsqueeze(2) * (eye - outer_hat) + lam_rad.unsqueeze(2) * outer_hat
        return F

    @staticmethod
    def patch_change_J(x: torch.Tensor, J: torch.Tensor) -> torch.Tensor:
        """
        Transform a (1,1)-tensor J from one patch to the other.

        Under the patch change x → x' with Jacobian F:
            J' = F J F⁻¹

        Args:
            x: Ball coordinates of the source patch, shape (batch, dim).
            J: Almost complex structure in source patch, shape (batch, dim, dim).

        Returns:
            J expressed in the target patch, shape (batch, dim, dim).
        """
        F = EvenSphere.patch_change_jacobian(x)
        F_inv = torch.linalg.inv(F)
        return F @ J @ F_inv

    # ------------------------------------------------------------------
    # Coordinate conversions (for evaluation / visualisation)
    # ------------------------------------------------------------------

    @staticmethod
    def ball_to_sphere(x: torch.Tensor, from_north: bool = True) -> torch.Tensor:
        """
        Map ball coordinates to ambient sphere coordinates in ℝ^{2n+1}.

        Args:
            x: Ball coordinates, shape (batch, 2n).
            from_north: True uses the patch-1 (north-pole-removed) embedding;
                False uses the patch-2 (south-pole-removed) embedding.

        Returns:
            Points on S^{2n} ⊂ ℝ^{2n+1}, shape (batch, 2n+1).
        """
        rho_sq = (x * x).sum(dim=1, keepdim=True)
        X = 2.0 * x / (1.0 - rho_sq).clamp(min=1e-8)
        R_sq = (X * X).sum(dim=1, keepdim=True)
        denom = 1.0 + R_sq
        sphere_xy = 2.0 * X / denom
        sphere_z  = ((R_sq - 1.0) / denom) if from_north else ((1.0 - R_sq) / denom)
        return torch.cat([sphere_xy, sphere_z], dim=1)


# ------------------------------------------------------------------
# Ball sampling  (pure NumPy, ported from AInstein/sampling/ball.py)
# ------------------------------------------------------------------

def ball_sample(
    num_pts: int,
    dimension: int,
    patch_width: float = 1.0,
    density_power: float = 1.0,
) -> np.ndarray:
    """
    Sample the 2n-ball with radial distribution giving uniform sphere measure.

    Args:
        num_pts: Number of sample points.
        dimension: Ball dimension (= real dimension of sphere = 2n).
        patch_width: Radial cutoff (default 1.0).
        density_power: Beta shape parameter.

    Returns:
        Ball coordinates, shape (num_pts, dimension).
    """
    half = num_pts // 2
    beta_b = density_power / (np.sqrt(2.0) - 1.0) - density_power

    radii_p1       = np.random.beta(density_power, beta_b, size=half)
    radii_p2       = np.random.beta(density_power, beta_b, size=num_pts - half)
    radii_p2_in_p1 = (1.0 - radii_p2) / (1.0 + radii_p2)
    radii = np.concatenate([radii_p1, radii_p2_in_p1]) * patch_width

    # Sample directions uniformly on S^{dimension-1} via isotropic Gaussian.
    # This is correct in all dimensions (unlike spherical-coordinate sampling
    # with uniform angles, which is only uniform on S^1).
    gaussian = np.random.randn(num_pts, dimension)
    norms    = np.linalg.norm(gaussian, axis=1, keepdims=True)
    unit_vec = gaussian / norms

    return radii.reshape(-1, 1) * unit_vec


def equatorial_annulus_sample(
    num_pts: int,
    dimension: int,
    r_inner: float,
    r_outer: float,
) -> np.ndarray:
    """
    Sample uniformly (in ball measure) from the annulus r_inner ≤ ρ ≤ r_outer.

    For the patch-change map φ(x) = x(ρ−1)/(ρ(ρ+1)), the equatorial band
    [r_inner, r_outer] is self-symmetric: φ maps it to itself.  Therefore
    sampling here gives points whose patch-2 images φ(x) also lie in the
    same annulus, making all samples active for the overlap loss without
    needing reweighting.

    Radii are drawn via inverse CDF of the d-ball radial density P(r) ∝ r^{d−1}:
        r = (r_inner^d + U · (r_outer^d − r_inner^d))^{1/d},  U ~ Uniform[0,1].

    Args:
        num_pts:   Number of sample points.
        dimension: Ball dimension d = real manifold dimension.
        r_inner:   Inner radius of the annulus.
        r_outer:   Outer radius of the annulus.

    Returns:
        Ball coordinates, shape (num_pts, dimension).
    """
    d = float(dimension)
    u = np.random.uniform(size=num_pts)
    radii = (r_inner ** d + u * (r_outer ** d - r_inner ** d)) ** (1.0 / d)

    gaussian = np.random.randn(num_pts, dimension)
    norms    = np.linalg.norm(gaussian, axis=1, keepdims=True)
    unit_vec = gaussian / norms

    return radii.reshape(-1, 1) * unit_vec
