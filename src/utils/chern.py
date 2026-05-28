"""
First Chern number and Chern curvature form computation for trained almost
complex structures.

Theory
------
For an almost complex manifold (M, J), the first Chern class is represented by
the Berry curvature of the T^{1,0} eigenspace bundle.  Given J(x) in ball
coordinates, define the projector

    P(x) = (I_{2n} + i J(x)) / 2  :  ℝ^{2n} ⊗ ℂ  →  T^{1,0}_x M

where i = √−1.  P(x) is the spectral projector for eigenvalue +i of J_ℂ.

Berry curvature 2-form (Chern–Weil representative of c₁):

    F_{kl}(x) = i · tr_ℂ[ P (∂_k P (I−P) ∂_l P − ∂_l P (I−P) ∂_k P) ]

which equals  (1/2π) c₁  in de Rham cohomology (independent of the choice of
connection).  For complex dimension n = 1 (S²) this simplifies to the standard
Berry phase formula  F_{kl} = 2 Im⟨∂_k ψ | ∂_l ψ⟩  with |ψ⟩ the unit +i
eigenvector of J_ℂ.

First Chern number (integrated):

    c₁[M] = (1/2π) ∫_M tr_ℂ(F) ∧ ω^{n-1}/(n-1)!

For S² (n=1):  c₁[S²] = (1/2π) ∫_{B²} F_{01} √det(g) d²x   (both patches)

Expected values (topological, independent of which J is learned):
    S²             :  c₁ = 2
    S²n (n ≥ 2)    :  c₁ = 0   (H²(S^{2n}) = 0)
    S¹×S¹ (T²)     :  c₁ = 0   (flat)
    S²×S² (prod J) :  c₁|_{S²×{pt}} = 2  (per factor; c₁² = 8, c₂ = 4)

Relation to N_J
---------------
For an almost Kähler manifold the Chern–Ricci form satisfies (Apostolov–Gauduchon):

    ρ^Ch = ρ^{LC} + (1/4)|N_J|² ω + ∇N_J correction terms

So as N_J → 0 the Chern curvature form converges to the Levi-Civita Ricci form
of the round metric.  The spatial standard deviation of the Chern curvature
density is therefore a complementary training diagnostic: it decreases as J
approaches an integrable structure.

Implementation notes
--------------------
The integrated Chern number is a homotopy-class invariant fixed at
initialization (J = exp(A) J₀ exp(−A) cannot change the homotopy class of J).
It serves as a topological checksum.  The spatially-resolved curvature density
carries per-epoch training signal.

Both `_compute_A_and_dA` and `_compute_dJ_from_dA` from `src.loss.sphere_loss`
are reused to obtain ∂J/∂x via the existing AD path.
"""

import math
import torch
import torch.nn as nn
from typing import Dict, Tuple

from src.geometry.sphere_product import SphereProduct
from src.loss.sphere_loss import _compute_A_and_dA, _compute_dJ_from_dA
from src.models.architecture import _matrix_exp


# ─────────────────────────────────────────────────────────────────────────────
# Core Berry-curvature computation
# ─────────────────────────────────────────────────────────────────────────────


def chern_projector(J: torch.Tensor) -> torch.Tensor:
    """
    Spectral projector P = (I + i J) / 2 onto T^{1,0}M.

    For each point x, P(x) maps ℝ^{2n} ⊗ ℂ onto the +i eigenspace of J_ℂ.

    Args:
        J: Real almost complex structure matrices, shape (batch, d, d).

    Returns:
        P: Complex projector, shape (batch, d, d), dtype complex64 / complex128.
    """
    d = J.shape[-1]
    I = torch.eye(d, dtype=J.dtype, device=J.device).unsqueeze(0)  # (1, d, d)
    J_c = J.to(dtype=torch.promote_types(J.dtype, torch.float32))
    I_c = I.expand_as(J_c)
    # (I + i J) / 2 as a complex tensor
    P = (I_c + 1j * J_c) / 2.0
    return P  # complex dtype


def berry_curvature(J: torch.Tensor, dJ: torch.Tensor) -> torch.Tensor:
    """
    Berry curvature F_{kl}(x) of the T^{1,0} eigenspace bundle.

    F_{kl} = i · tr_ℂ[ P (∂_k P (I−P) ∂_l P − ∂_l P (I−P) ∂_k P) ]

    where P = (I + i J) / 2 and ∂_k P = (i / 2) ∂_k J.

    F is real-valued (the imaginary parts cancel exactly) and antisymmetric
    in (k, l).

    Args:
        J:  J matrices in ball coords, shape (batch, d, d),  real.
        dJ: ∂J/∂x,  shape (batch, d, d, d_ball),  real.
            dJ[b, α, β, k]  =  ∂_k J^α_β  at batch sample b.

    Returns:
        F: Berry curvature 2-form, shape (batch, d_ball, d_ball),  real.
           F[b, k, l] = F_{kl}(x_b).
    """
    batch, d, _, d_ball = dJ.shape

    # Complex projector P = (I + i J)/2 and its complement Q = I - P
    P     = chern_projector(J)                                     # (batch, d, d), complex
    eye_c = torch.eye(d, dtype=P.dtype, device=P.device).unsqueeze(0)
    Q     = eye_c - P                                              # (batch, d, d), complex

    # ∂_k P = (i/2) ∂_k J  — shape (batch, d, d, d_ball), complex
    dP = (0.5j) * dJ.to(dtype=P.dtype)

    # dP_perm[k, b, α, β] = dP[b, α, β, k]
    dP_perm = dP.permute(3, 0, 1, 2)                              # (d_ball, batch, d, d)

    # PdPkQ[k, b, :, :] = P[b] @ dP_k[b] @ Q[b]
    PdPkQ = torch.stack([P @ dP_perm[k] @ Q for k in range(d_ball)])  # (d_ball, batch, d, d)

    # F_half[b,k,l] = i · tr(P dP_k Q dP_l)
    #               = i · Σ_{i,j} PdPkQ[k,b,i,j] · dP_perm[l,b,j,i]
    F_half = 1j * torch.einsum('kbij,lbji->bkl', PdPkQ, dP_perm) # (batch, d_ball, d_ball)

    # Full formula is antisymmetric: F_{kl} = i·tr[P dP_k Q dP_l − P dP_l Q dP_k]
    F_complex = F_half - F_half.permute(0, 2, 1)                  # (batch, d_ball, d_ball)

    # F is real by construction; take real part (imaginary part is numerical noise)
    return F_complex.real                                          # (batch, d_ball, d_ball)


# ─────────────────────────────────────────────────────────────────────────────
# Chern curvature density (scalar, for integration and plotting)
# ─────────────────────────────────────────────────────────────────────────────


def chern_density(
    J: torch.Tensor,
    dJ: torch.Tensor,
    g_sqrt_det: torch.Tensor,
) -> torch.Tensor:
    """
    Chern curvature density  (1/2π) F_{01…}  ∧ vol weighted by √det(g).

    For a 2n-real-dimensional manifold, the contribution to the first Chern
    number from a small ball-coordinate patch element d^{2n}x is:

        dc₁ = (1/2π) · (∑_{k<l} F_{kl} ε^{kl…} ) · √det(g) d^{2n}x

    For S² (d_ball = 2):  dc₁ = (1/2π) F_{01} √det(g) dx₀dx₁
    For d_ball > 2:  dc₁ = (1/2π) Σ_{k<l} F_{kl} · (top-degree Hodge dual) √det(g) d^{2n}x

    For general dimension, the first Chern form integrates against the
    fundamental class as  ∫ c₁ ∧ ω^{n-1}/(n-1)!  where ω is the fundamental
    2-form.  Here we compute the simpler trace density:

        density(x) = (1/2π) tr_ℂ(F(x)) · √det(g(x))

    where tr_ℂ(F(x)) is the trace of F viewed as a 2-form contracted with
    the compatible complex structure (see note below).  For d_ball = 2 this
    is just F_{01}.

    Note: tr_ℂ(F) = Σ_{k<l} ω_{kl}(x) F_{kl}(x) where ω_{kl} = g_{ki} J^i_l
    is the fundamental form in ball coordinates.  Integrating this density over
    all patches (both charts) and dividing by 2π gives c₁[M].

    Args:
        J:           J matrices, shape (batch, d, d),  real.
        dJ:          ∂J/∂x, shape (batch, d, d, d_ball),  real.
        g_sqrt_det:  √det(g(x)),  shape (batch,).

    Returns:
        density: Chern curvature density per sample, shape (batch,),  real.
                 Integrate (sum × dx-volume) to get the Chern number contribution.
    """
    F = berry_curvature(J, dJ)              # (batch, d_ball, d_ball)
    d_ball = F.shape[-1]

    if d_ball == 2:
        # F_{01} is the unique independent component; this is the standard Berry phase.
        tr_F = F[:, 0, 1]                                          # (batch,)
    else:
        # Sum all upper-triangle components F_{kl}, k < l.
        # For S^{2n} with n≥2 this integrates to zero (H²(S^{2n})=0),
        # confirming the topological constraint c₁=0.
        upper = [F[:, k, l] for k in range(d_ball) for l in range(k + 1, d_ball)]
        tr_F  = torch.stack(upper, dim=1).sum(dim=1)               # (batch,)

    return (1.0 / (2.0 * math.pi)) * tr_F * g_sqrt_det             # (batch,)


def _sqrt_det_g(x: torch.Tensor) -> torch.Tensor:
    """
    √det(g(x)) for the round metric in ball coordinates.

    From the metric formula  g = A_tan I + B x⊗x  with
        A_tan = 16(1−ρ²)² / (1+ρ²)⁴
        B     = 64 / (1+ρ²)⁴

    and using det(A I + B v⊗v) = A^{d-1} (A + B ρ²), we get:

        det(g) = A_tan^{d-1} (A_tan + B ρ²)
               = [16(1−ρ²)²/(1+ρ²)⁴]^{d-1} · [16(1−ρ²)²/(1+ρ²)⁴ + 64ρ²/(1+ρ²)⁴]
               = [16(1−ρ²)²/(1+ρ²)⁴]^{d-1} · [16(1−ρ²)² + 64ρ²] / (1+ρ²)⁴

    Note: 16(1−ρ²)² + 64ρ² = 16(1 + ρ²)².  So:

        det(g) = [16(1−ρ²)² / (1+ρ²)⁴]^{d-1} · 16(1+ρ²)² / (1+ρ²)⁴
               = 16^d · (1−ρ²)^{2(d-1)} (1+ρ²)^{2-4d}

    and  √det(g) = 4^d · (1−ρ²)^{d-1} · (1+ρ²)^{1-2d}.

    Args:
        x: Ball coordinates, shape (batch, d).

    Returns:
        √det(g(x)), shape (batch,).
    """
    d = x.shape[1]
    rho_sq  = (x * x).sum(dim=1)                            # (batch,)
    op2     = 1.0 + rho_sq                                  # (1+ρ²)
    om2     = (1.0 - rho_sq).clamp(min=1e-8)                # (1−ρ²)
    return (4.0 ** d) * (om2 ** (d - 1)) * (op2 ** (1 - 2 * d))


# ─────────────────────────────────────────────────────────────────────────────
# Monte Carlo integration over S^{2n}
# ─────────────────────────────────────────────────────────────────────────────


def integrate_chern_number(
    subnet: nn.Module,
    manifold_dim: int,
    n_pts: int = 8000,
    device: torch.device | None = None,
) -> Tuple[float, torch.Tensor, torch.Tensor]:
    """
    Estimate the first Chern number c₁[M] via Monte Carlo integration over S^{2n}.

    The integral is split over both patches:

        c₁[M] = (contribution from patch 1) + (contribution from patch 2)

    Each contribution is estimated as:

        ∫_{B^{2n}} (1/2π) · F_{01…}(x) · √det(g(x)) · d^{2n}x
            ≈  (vol(B^{2n}) / N) · Σ_i density(x_i)     [uniform Monte Carlo]

    where  density = (1/2π) tr_ℂ(F) √det(g).

    However, to be exact: since we sample with the round sphere measure (not the
    flat ball measure), the Monte Carlo estimator is the sample mean of the
    density divided by the sampling PDF and times the total volume.  In practice
    we use a self-normalised estimator:

        ∫_M f dVol_g  ≈  (1/N) Σ_i f(x_i)   when x_i ~ dVol_g  (normalised)

    For the Chern number:

        c₁[M] = ∫_M (1/2π) F ∧ ω^{n-1}/(n-1)!

    In ball coordinates, sampling uniformly with the round metric on each patch:
    both patches cover the same hemisphere (the intersection has measure zero),
    and together they cover all of M.  The round-metric sample on patch 1 covers
    the top hemisphere (patch-1 coordinates x with |x| < 1); patch 2 covers the
    bottom hemisphere.  Since the sampling density is d Vol_g / Vol_g, the
    estimate is:

        c₁[M]  ≈  mean_{x ~ dVol_g, patch1}[density(x)] · Vol(S^{2n}) / 2
                  + mean_{x ~ dVol_g, patch2}[density(x)] · Vol(S^{2n}) / 2

    But the density is defined in ball coordinates, where the integration measure
    is d^{2n}x (flat Lebesgue), not dVol_g.  The relationship is:

        ∫_M f dVol_g = ∫_{B^{2n}} f(x) √det(g(x)) d^{2n}x

    We therefore need to use uniform-in-ball sampling (NOT the round-metric
    sampling from sphere_uniform_sample) to avoid introducing a bias.

    Practically, `sphere_uniform_sample` already returns samples with the round
    sphere measure in ball coordinates with radial PDF proportional to the
    pull-back of the round volume form.  To convert to a flat-ball integral
    estimate, we weight by 1/(sampling density).

    For simplicity and correctness we instead use the following exact relation:

        c₁[M] = (1/2π) ∫_{B^{2n}} F_{01} √det(g) d²x  [for S², d_ball=2]

    Estimated via importance-sampled Monte Carlo, sampling x uniformly on S^{2n}
    in ball coords (density ∝ √det(g(x))), so:

        c₁[M] ≈ (1/N) Σ_i (1/2π) F_{01}(x_i)   [importance-sampled mean]

    This means: sample with the round-sphere measure, then average (1/2π) F_{01}(x_i),
    and the result converges to  (1/Area(S^{2n})) ∫ (1/2π) F_{01} dVol_g,
    which is NOT c₁.  We must multiply by Area(S^{2n}) = 2π^{(2n+1)/2} / Γ((2n+1)/2).

    To avoid this complication, we use a **different strategy**: sample uniformly
    in the ball using the flat Lebesgue measure via rejection sampling, and
    compute the mean of (1/2π) F_{01} √det(g).

    Args:
        subnet:       PatchSubNet (single patch) with J₀ and forward_raw.
        manifold_dim: Real dimension d = 2n.
        n_pts:        Number of Monte Carlo samples per patch.
        device:       Torch device.

    Returns:
        c1:      Estimated first Chern number (scalar float).
        x_pts:   Sample points in ball coords, shape (n_pts, manifold_dim).
        density: Chern density at each point, shape (n_pts,).
    """
    if device is None:
        device = next(subnet.parameters()).device
    dtype = next(subnet.parameters()).dtype

    # Sample uniformly inside the unit ball B^{2n} using Muller's method:
    # direction uniform on S^{d-1}, radius with CDF r^d (beta(1, d/2 + 1) on r²).
    # Here we use simple rejection: sample a hypercube, keep points inside ball.
    x_pts = _sample_uniform_ball(n_pts, manifold_dim, device, dtype)

    with torch.enable_grad():
        A, dA = _compute_A_and_dA(subnet, x_pts)
        # Use _matrix_exp (Padé-13 with MPS fallback) for consistency with training
        J  = _matrix_exp(A) @ subnet.J0 @ _matrix_exp(-A)            # (n_pts, d, d)
        dJ = _compute_dJ_from_dA(A, dA, subnet.J0)                   # (n_pts, d, d, d)

    sqrt_det = _sqrt_det_g(x_pts.detach())                            # (n_pts,)

    with torch.no_grad():
        dens = chern_density(J.detach(), dJ.detach(), sqrt_det)       # (n_pts,)

    # Volume of the unit ball B^{2n}:  V = π^n / n!
    n = manifold_dim // 2
    ball_vol = (math.pi ** n) / math.factorial(n)

    c1 = (ball_vol / n_pts) * dens.sum().item()
    return c1, x_pts.detach(), dens.detach()


def compute_chern_numbers(
    model: nn.Module,
    manifold_type: str,
    manifold_dim: int | None = None,
    dim_alpha: int | None = None,
    dim_beta: int | None = None,
    n_pts: int = 8000,
    device: torch.device | None = None,
) -> Dict[str, float]:
    """
    Compute first Chern number(s) for a trained model.

    For 'sphere':
        Integrates over both patches of S^{2n}.
        Returns {'c1': value}.  Expected: 2.0 for S², 0.0 for S^{2n} (n≥2).

    For 'sphere_product':
        Integrates over each factor independently by restricting to each
        patch pair.  Returns {'c1_alpha': ..., 'c1_beta': ..., 'c1_total': ...}.

    Args:
        model:         TwoPatchAlmostComplexNet or FourPatchProductNet.
        manifold_type: 'sphere' or 'sphere_product'.
        manifold_dim:  Real dimension of S^{2n} (sphere only).
        dim_alpha:     m_α (sphere_product only).
        dim_beta:      m_β (sphere_product only).
        n_pts:         Monte Carlo sample count per patch.
        device:        Torch device.

    Returns:
        Dictionary of Chern number estimates and spatial curvature statistics.
    """
    if device is None:
        device = next(model.parameters()).device

    results: Dict[str, float] = {}

    if manifold_type == 'sphere':
        assert manifold_dim is not None
        # Patch 1
        c1_p1, x1, d1 = integrate_chern_number(model.patch1, manifold_dim, n_pts, device)
        # Patch 2
        c1_p2, x2, d2 = integrate_chern_number(model.patch2, manifold_dim, n_pts, device)
        c1 = c1_p1 + c1_p2
        density_all = torch.cat([d1, d2])
        results['c1'] = c1
        results['c1_patch1'] = c1_p1
        results['c1_patch2'] = c1_p2
        results['chern_density_mean'] = density_all.mean().item()
        results['chern_density_std']  = density_all.std().item()

    elif manifold_type == 'sphere_product':
        assert dim_alpha is not None and dim_beta is not None
        # For each patch combination (i_α, i_β), the subnet sees the full
        # d = dim_alpha + dim_beta dimensional input.  Restricting c₁ to each
        # factor requires the full-dimensional integration; here we compute
        # the total c₁ on the product manifold by integrating the Berry
        # curvature of the product ACS.  For a product J = J_α ⊕ J_β:
        #   c₁(M_α × M_β) = c₁(M_α) [M_α×{pt}] + c₁(M_β) [{pt}×M_β]
        # We estimate c₁ for each factor by restricting to a d_ball = 2 slice.
        results.update(_chern_product(model, dim_alpha, dim_beta, n_pts, device))

    return results


def integrate_chern_number_product(
    model: nn.Module,
    sphere_product: SphereProduct,
    n_pts: int = 8000,
    device: torch.device | None = None,
) -> Dict[str, float]:
    """
    Estimate first Chern numbers for each factor in S^{m_α} × S^{m_β}.

    Args:
        model: FourPatchProductNet model.
        sphere_product: Product manifold carrying factor dimensions.
        n_pts: Monte Carlo samples per factor.
        device: Torch device.

    Returns:
        Dictionary with c1_alpha, c1_beta, c1_total and density stats.
    """
    if device is None:
        device = next(model.parameters()).device
    return _chern_product(
        model=model,
        dim_alpha=sphere_product.dim_alpha,
        dim_beta=sphere_product.dim_beta,
        n_pts=n_pts,
        device=device,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Sphere-product helper
# ─────────────────────────────────────────────────────────────────────────────


def _chern_product(
    model: nn.Module,
    dim_alpha: int,
    dim_beta: int,
    n_pts: int,
    device: torch.device,
) -> Dict[str, float]:
    """
    Chern numbers for S^{m_α} × S^{m_β} from the FourPatchProductNet.

    For the product, the first Chern class of (M_α × M_β, J) has two
    components:

        c₁(α) = ∫_{M_α × {pt_β}} c₁     (restriction to α-factor)
        c₁(β) = ∫_{{pt_α} × M_β} c₁     (restriction to β-factor)

    We estimate these by integrating over each factor subnet separately at a
    fixed reference point for the other factor.  For sub-nets that see the
    full concatenated input, we fix x_β = 0 (origin of the β-ball chart, =
    north pole) and integrate over x_α, and vice versa.

    For patch (0,0) the subnet `model.subnets[(0,0)]` takes
    x = cat(x_α, x_β) ∈ B^{dim_alpha} × B^{dim_beta}.  We fix x_β = 0
    and compute the Berry curvature restricted to the x_α directions.

    For the per-factor integration this is exact only when J is the product
    structure J = J_α ⊕ J_β; for a general J the restriction gives an
    approximation.
    """
    dtype = next(model.parameters()).dtype
    results: Dict[str, float] = {}
    d_alpha = dim_alpha
    d_beta  = dim_beta
    d_total = d_alpha + d_beta

    # ── α-factor: integrate over x_α with x_β = 0 ──
    x_alpha = _sample_uniform_ball(n_pts, d_alpha, device, dtype)
    x_beta0 = torch.zeros(n_pts, d_beta, device=device, dtype=dtype)
    x_full  = torch.cat([x_alpha, x_beta0], dim=1)   # (n_pts, d_total)

    subnet_00 = _get_subnet(model, 0, 0)

    with torch.enable_grad():
        A_full, dA_full = _compute_A_and_dA(subnet_00, x_full)
        J_full  = _matrix_exp(A_full) @ subnet_00.J0 @ _matrix_exp(-A_full)
        dJ_full = _compute_dJ_from_dA(A_full, dA_full, subnet_00.J0)   # (n_pts,d,d,d_total)

    # Restrict dJ to the α-directions (first d_alpha columns)
    dJ_alpha = dJ_full[:, :, :, :d_alpha]          # (n_pts, d_total, d_total, d_alpha)
    sqrt_det_a = _sqrt_det_g(x_alpha.detach())

    with torch.no_grad():
        dens_a = chern_density(J_full.detach(), dJ_alpha.detach(), sqrt_det_a)

    n = d_alpha // 2 if d_alpha % 2 == 0 else None
    if n is not None:
        ball_vol_a = (math.pi ** n) / math.factorial(n)
        c1_a = (ball_vol_a / n_pts) * dens_a.sum().item()
    else:
        # Odd-dimensional factor: c₁ contribution requires the full complex structure
        c1_a = float('nan')

    # ── β-factor: integrate over x_β with x_α = 0 ──
    x_alpha0 = torch.zeros(n_pts, d_alpha, device=device, dtype=dtype)
    x_beta   = _sample_uniform_ball(n_pts, d_beta,  device, dtype)
    x_full2  = torch.cat([x_alpha0, x_beta], dim=1)

    with torch.enable_grad():
        A_full2, dA_full2 = _compute_A_and_dA(subnet_00, x_full2)
        J_full2  = _matrix_exp(A_full2) @ subnet_00.J0 @ _matrix_exp(-A_full2)
        dJ_full2 = _compute_dJ_from_dA(A_full2, dA_full2, subnet_00.J0)

    dJ_beta = dJ_full2[:, :, :, d_alpha:]          # β-directions
    sqrt_det_b = _sqrt_det_g(x_beta.detach())

    with torch.no_grad():
        dens_b = chern_density(J_full2.detach(), dJ_beta.detach(), sqrt_det_b)

    n_b = d_beta // 2 if d_beta % 2 == 0 else None
    if n_b is not None:
        ball_vol_b = (math.pi ** n_b) / math.factorial(n_b)
        c1_b = (ball_vol_b / n_pts) * dens_b.sum().item()
    else:
        c1_b = float('nan')

    results['c1_alpha'] = c1_a
    results['c1_beta']  = c1_b
    results['c1_total'] = (c1_a if not math.isnan(c1_a) else 0.0) + \
                          (c1_b if not math.isnan(c1_b) else 0.0)
    results['chern_density_alpha_mean'] = dens_a.mean().item()
    results['chern_density_alpha_std']  = dens_a.std().item()
    results['chern_density_beta_mean']  = dens_b.mean().item()
    results['chern_density_beta_std']   = dens_b.std().item()
    return results


def _get_subnet(model: nn.Module, i_alpha: int, i_beta: int) -> nn.Module:
    """Extract the (i_α, i_β) PatchSubNet from a FourPatchProductNet."""
    return model.patches[f'patch_{i_alpha}{i_beta}']


# ─────────────────────────────────────────────────────────────────────────────
# Sampling helpers
# ─────────────────────────────────────────────────────────────────────────────


def _sample_uniform_ball(
    n: int,
    d: int,
    device: torch.device,
    dtype: torch.dtype,
) -> torch.Tensor:
    """
    Sample n points uniformly inside the open unit ball B^d.

    Uses the standard method: sample uniformly on the sphere S^{d-1} for
    direction, and r^{1/d} for radius where r ~ Uniform[0,1).

    Args:
        n: Number of points.
        d: Ball dimension.
        device, dtype: Tensor device and dtype.

    Returns:
        Points, shape (n, d), all with ‖x‖ < 1.
    """
    # Uniform direction on S^{d-1}
    z = torch.randn(n, d, device=device, dtype=dtype)
    z = z / z.norm(dim=1, keepdim=True).clamp(min=1e-10)
    # Uniform radius (marginal PDF for ‖x‖ in unit ball is d r^{d-1})
    r = torch.rand(n, 1, device=device, dtype=dtype) ** (1.0 / d)
    return z * r * 0.98    # 0.98 to stay strictly inside (avoid |x|→1 singularity)


# ─────────────────────────────────────────────────────────────────────────────
# Quick self-test
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    """
    Self-test on constant J = J₀.

    For A = 0, J is constant in coordinates so dJ = 0 and the local Berry
    curvature representative computed from dJ vanishes pointwise.
    """
    import sys
    sys.path.insert(0, '.')

    from src.models import TwoPatchAlmostComplexNet

    torch.manual_seed(0)
    model = TwoPatchAlmostComplexNet(
        manifold_dim=2,
        hidden_layers=[32, 64, 32],
        activation='gelu',
    )
    # Zero all weights → A = 0 → J = J₀
    for p in model.parameters():
        p.data.zero_()

    result = compute_chern_numbers(
        model,
        manifold_type='sphere',
        manifold_dim=2,
        n_pts=20000,
    )
    print("Self-test: local Chern density with constant J = J₀")
    for k, v in result.items():
        print(f"  {k}: {v:.4f}")
    print("  Expected local-density integral near 0.0 for constant J in this gauge")
