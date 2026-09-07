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
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

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
        density_power: float = 1.0,
        overlap_size: float | None = None,
    ) -> torch.Tensor:
        """
        Sample points in patch-0 ball coordinates from S^{2n}.

        When overlap_size is given, samples are drawn from the equatorial
        annulus [r_inner, r_outer] (r_outer = r_eq + overlap_size;
        r_inner = (1−r_outer)/(1+r_outer)) with uniform ball measure.
        The annulus is self-symmetric under the patch-change map so the
        patch-1 images φ(x) lie in the same annulus.

        Otherwise, samples are drawn using sphere_uniform_sample with
        the round measure at density_power=1.

        Args:
            num_points:    Number of points to sample.
            density_power: Latitude density parameter p > 0.  p=1: round
                measure.  p>1: equatorial bias.  p<1: polar bias.
            overlap_size:  If set, sample from equatorial annulus only.

        Returns:
            Ball coordinates in patch 0, shape (num_points, dim).
        """
        if overlap_size is not None:
            r_outer = RADIAL_MIDPOINT + overlap_size
            r_inner = (1.0 - r_outer) / (1.0 + r_outer)
            pts_np = equatorial_annulus_sample(num_points, self.dim, r_inner, r_outer)
        else:
            pts_np = sphere_uniform_sample(num_points, self.dim, density_power)
        return torch.tensor(pts_np, dtype=torch.float64, device=self.device)

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

    @staticmethod
    def round_metric_ball_coords(
        x: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Round metric and its inverse in ball coordinates.

        The unit sphere S^{2n} is embedded via
            p_k = 4 x_k (1−ρ²)/(1+ρ²)²   (k = 1,…,2n)
            p_{2n+1} = (−1 + 6ρ² − ρ⁴)/(1+ρ²)²

        Pulling back the flat metric on ℝ^{2n+1} gives a structured metric
        on the open ball B(0,1) ⊂ ℝ^{2n}:

            g_{ij}(x) = A_tan(ρ) δ_{ij} + B(ρ) x_i x_j

        with radial and tangential eigenvalues:
            A_tan = 16(1−ρ²)² / (1+ρ²)⁴
            A_rad = A_tan + B ρ²  = 16 / (1+ρ²)²
            B     = 64 / (1+ρ²)⁴

        The inverse (by Sherman–Morrison) is:
            g^{ij}(x) = inv_A_tan(ρ) δ_{ij} + D(ρ) x_i x_j

        with:
            inv_A_tan = (1+ρ²)⁴ / (16(1−ρ²)²)
            D         = −(1+ρ²)² / (4(1−ρ²)²)

        Args:
            x: Ball coordinates, shape (batch, dim).

        Returns:
            g:     Metric tensor, shape (batch, dim, dim).
            g_inv: Inverse metric tensor, shape (batch, dim, dim).
        """
        batch, d = x.shape
        rho_sq  = (x * x).sum(dim=1)                    # (batch,)
        op2     = 1.0 + rho_sq                           # (1+ρ²)
        om2     = (1.0 - rho_sq).clamp(min=1e-8)         # (1−ρ²), guarded
        op2_4   = op2 ** 4                               # (1+ρ²)⁴

        A_tan = 16.0 * om2 ** 2 / op2_4                 # (batch,)
        B     = 64.0 / op2_4                             # (batch,)

        eye      = torch.eye(d, dtype=x.dtype, device=x.device).unsqueeze(0)
        outer_xx = x.unsqueeze(2) * x.unsqueeze(1)      # (batch, d, d)

        g = A_tan[:, None, None] * eye + B[:, None, None] * outer_xx

        # Inverse: use (1+ρ²)⁴ / (16(1−ρ²)²) for the diagonal and
        #          −(1+ρ²)² / (4(1−ρ²)²) for the x⊗x coefficient.
        om2_sq    = om2 ** 2                              # (1−ρ²)²
        inv_A_tan = op2_4 / (16.0 * om2_sq)              # (batch,)
        D         = -(op2 ** 2) / (4.0 * om2_sq)         # (batch,)

        g_inv = inv_A_tan[:, None, None] * eye + D[:, None, None] * outer_xx

        return g, g_inv

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
# Sphere sampling  (pure NumPy)
# ------------------------------------------------------------------

def sphere_uniform_sample(
    num_pts: int,
    dimension: int,
    density_power: float = 1.0,
) -> np.ndarray:
    """
    Sample from S^{dimension} and return patch-0 ball coordinates.

    The coarea formula gives the marginal density of the height coordinate
    q = p_{d+1} under the round measure on S^d as

        P(q) ∝ (1 − q²)^{(d−2)/2}

    For density_power = p, the family P_p(q) ∝ (1−q²)^{p(d−2)/2} is used.
    Setting u = (1+q)/2 ∈ [0,1] rewrites this as Beta(α+1, α+1) with
    α = p(d−2)/2.  At p=1, this is the round measure on S^d.  At p>1,
    samples concentrate toward the equator (q≈0, ρ≈r_eq in ball coords);
    at p<1, they concentrate toward the poles (|q|≈1).

    The direction v̂ is sampled uniformly on S^{d−1} (normalised Gaussian).
    The ambient point (√(1−q²) v̂, q) is mapped to patch-0 ball coordinates
    via the closed-form inverse of the modified stereographic projection,
    using a hemisphere split for numerical stability:

        southern (q ≤ 0):  x = √(1+q) / (√(1−q) + √2) · v̂
        northern (q > 0):  ρ₂ = √(1−q) / (√(1+q) + √2)
                           x  = −(1−ρ₂)/(1+ρ₂) · v̂

    Args:
        num_pts:       Number of sample points.
        dimension:     Real dimension d of the sphere S^d (≥ 1).
        density_power: Latitude density parameter p > 0.

    Returns:
        Ball coordinates in patch 0, shape (num_pts, dimension).

    Raises:
        ValueError: if density_power is out of the valid range for dimension.
    """
    d = dimension
    alpha = density_power * (d - 2) / 2.0
    beta_param = alpha + 1.0
    if beta_param <= 0.0:
        raise ValueError(
            f"density_power={density_power} is invalid for dimension d={d}: "
            f"requires density_power < {2.0 / (2 - d):.4g}"
        )

    u     = np.random.beta(beta_param, beta_param, size=num_pts)
    q     = 2.0 * u - 1.0  # height q ∈ (−1, 1)
    gauss = np.random.randn(num_pts, d)
    v_hat = gauss / np.linalg.norm(gauss, axis=1, keepdims=True)

    sqrt_1pq = np.sqrt(np.maximum(1.0 + q, 0.0))  # √(1+q)
    sqrt_1mq = np.sqrt(np.maximum(1.0 - q, 0.0))  # √(1−q)
    sqrt2    = np.sqrt(2.0)

    x     = np.empty((num_pts, d))
    south = q <= 0.0
    north = ~south

    # Southern hemisphere: patch-0 projection is well-conditioned near south pole.
    if south.any():
        scale    = sqrt_1pq[south] / (sqrt_1mq[south] + sqrt2)
        x[south] = scale[:, np.newaxis] * v_hat[south]

    # Northern hemisphere: use patch-1 coords (well-conditioned near north pole),
    # then apply the patch-change map φ.  The direction flips sign under φ.
    if north.any():
        rho2     = sqrt_1mq[north] / (sqrt_1pq[north] + sqrt2)
        rho1     = (1.0 - rho2) / (1.0 + rho2)
        x[north] = -rho1[:, np.newaxis] * v_hat[north]

    return x


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
