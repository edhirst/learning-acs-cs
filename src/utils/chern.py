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

    F_{kl}(x) = i · tr_ℂ[ P R^TM_{kl}
                          + P (∇_k P (I−P) ∇_l P
                               − ∇_l P (I−P) ∇_k P) ]

where R^TM is the round Levi-Civita curvature of the tangent bundle.  The
ambient curvature term is essential: a locally constant expression for J can
still represent a nontrivial tangent bundle with nonzero c₁.

First Chern number in complex dimension n:

    ∫_M c₁^n = ∫_M (F / 2π)^n

where F is the real Chern-Weil first-Chern curvature 2-form.  In local
coordinates this top-form coefficient is

    n! Pf(F) / (2π)^n,

with Pf the Pfaffian of the antisymmetric matrix F_ij.  Test samples are drawn
from the round volume measure, so the Monte Carlo integrand divides this
coordinate top-form coefficient by sqrt(det(g)).

For sphere products this is computed on the full product manifold, not on
individual sphere factors.

Reference cases:
    S²             :  c₁ = 2
    S²n (n ≥ 2)    :  ∫ c₁^n = 0   (H²(S^{2n}) = 0)
    S¹×S¹ (T²)     :  c₁ = 0   (flat)
    higher products:  reported as ∫ c₁^n over the full product manifold

Implementation notes
--------------------
The reported scalar is computed only for final test/evaluation reporting; no
Chern-related quantity is tracked during training or validation.  In real
dimension 2 it is the usual first Chern number.  In higher real dimension it is
the top-degree Chern number ∫ c₁^n.

Both `_compute_A_and_dA` and `_compute_dJ_from_dA` from `src.loss.sphere_loss`
are reused to obtain ∂J/∂x via the existing AD path.
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import math
import os
import sys
import torch
import torch.func as func
import torch.nn as nn
from typing import Tuple

if __name__ == '__main__' and __package__ is None:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from src.geometry import EvenSphere, RADIAL_MIDPOINT
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

    F_{kl} = i · tr_ℂ[ P (∇_k P (I−P) ∇_l P − ∇_l P (I−P) ∇_k P) ]

    where P = (I + i J) / 2 and ∇_k P = (i / 2) ∇_k J.

    F is real-valued (the imaginary parts cancel exactly) and antisymmetric
    in (k, l).

    Args:
        J:  J matrices in ball coords, shape (batch, d, d),  real.
        dJ: ∇J, shape (batch, d, d, d_ball), real.
            dJ[b, α, β, k] = ∇_k J^α_β at batch sample b.

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


def _projected_ambient_curvature(
    J: torch.Tensor,
    riemann: torch.Tensor,
) -> torch.Tensor:
    """
    Curvature contribution from the ambient tangent-bundle connection.

    Args:
        J:       J matrices, shape (batch, d, d).
        riemann: R^a_{bkl} of the round/product Levi-Civita connection,
                 shape (batch, d, d, d, d).

    Returns:
        Real 2-form i tr(P R^TM), shape (batch, d, d).
    """
    P = chern_projector(J)
    tr_pr = torch.einsum('bij,bjikl->bkl', P, riemann.to(dtype=P.dtype))
    return (1j * tr_pr).real


def chern_curvature(
    J: torch.Tensor,
    dJ: torch.Tensor,
    riemann: torch.Tensor,
) -> torch.Tensor:
    """
    First-Chern curvature 2-form for the projected tangent-bundle connection.

    This is the old local Berry/projector term plus the ambient round
    Levi-Civita curvature of TM.  Omitting the ambient term incorrectly gives
    zero for locally constant chart representatives of the standard S² ACS.
    """
    return berry_curvature(J, dJ) + _projected_ambient_curvature(J, riemann)


# ─────────────────────────────────────────────────────────────────────────────
# Global first Chern number on held-out test samples
# ─────────────────────────────────────────────────────────────────────────────


def _pfaffian(a: torch.Tensor) -> torch.Tensor:
    """
    Batched Pfaffian for small antisymmetric matrices.

    For a 2-form F = Σ_{i<j} F_ij dx_i∧dx_j, F^n has coordinate coefficient
    n! Pf(F).  The dimensions used in this project are small, so a recursive
    expansion is simpler and less error-prone than a custom factorization.
    """
    d = a.shape[-1]
    if d % 2 != 0:
        raise ValueError(f"Pfaffian requires an even matrix dimension, got {d}")
    if d == 0:
        return a.new_ones(a.shape[:-2])
    if d == 2:
        return a[..., 0, 1]

    total = a.new_zeros(a.shape[:-2])
    for j in range(1, d):
        idx = [k for k in range(d) if k not in (0, j)]
        minor = a[..., idx, :][..., :, idx]
        sign = -1.0 if (j % 2 == 0) else 1.0
        total = total + sign * a[..., 0, j] * _pfaffian(minor)
    return total


def _christoffel_symbols(
    x: torch.Tensor,
    metric_fn,
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    Compute round Levi-Civita Christoffel symbols in local coordinates.

    Returns Γ^a_{kb} with shape (batch, d, d, d), where k is the derivative
    coordinate and b is the lower tensor index being transported.
    """
    def single_metric(x_single: torch.Tensor) -> torch.Tensor:
        g_single, _ = metric_fn(x_single.unsqueeze(0))
        return g_single.squeeze(0)

    g, g_inv = metric_fn(x)
    dg = func.vmap(func.jacfwd(single_metric))(x)  # dg[b, i, j, k] = ∂_k g_ij
    gamma_terms = (
        dg.permute(0, 2, 3, 1)  # ∂_k g_cm
        + dg.permute(0, 2, 1, 3)  # ∂_c g_km
        - dg.permute(0, 3, 1, 2)  # ∂_m g_kc
    )
    gamma = 0.5 * torch.einsum('bam,bmkc->bakc', g_inv, gamma_terms)
    return g, g_inv, gamma


def _sphere_metric_connection(
    x: torch.Tensor,
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Round metric, inverse metric, and Christoffel symbols for one sphere chart."""
    return _christoffel_symbols(x, EvenSphere.round_metric_ball_coords)


def _product_metric_connection(
    x_alpha: torch.Tensor,
    x_beta: torch.Tensor,
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Round product metric, inverse metric, and Christoffel symbols."""
    dim_alpha = x_alpha.shape[1]

    def metric_fn(x_full: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        return SphereProduct.round_metric_ball_coords(
            x_full[:, :dim_alpha],
            x_full[:, dim_alpha:],
        )

    return _christoffel_symbols(torch.cat([x_alpha, x_beta], dim=1), metric_fn)


def _covariant_dJ(
    J: torch.Tensor,
    dJ: torch.Tensor,
    gamma: torch.Tensor,
) -> torch.Tensor:
    """
    Convert coordinate partial derivatives of J into covariant derivatives.

    dJ[b, a, c, k] is ∂_k J^a_c and gamma[b, a, k, m] is Γ^a_{km}.
    """
    upper = torch.einsum('bakm,bmc->back', gamma, J)
    lower = torch.einsum('bmkc,bam->back', gamma, J)
    return dJ + upper - lower


def _first_chern_integrand(
    J: torch.Tensor,
    dJ: torch.Tensor,
    g: torch.Tensor,
    riemann: torch.Tensor,
    orientation_sign: float,
) -> torch.Tensor:
    """
    Pointwise scalar whose round-volume integral gives ∫ c₁^n.

    The Chern-Weil first-Chern form is c₁ = F/(2π).  In real dimension 2n,
    c₁^n has coordinate top-form coefficient

        n! Pf(F) / (2π)^n.

    The sample distribution is the round volume measure, so divide by
    sqrt(det(g)) to convert coordinate top-form density into a volume density.
    """
    F = chern_curvature(J, dJ, riemann)                             # (batch, d, d)
    F = 0.5 * (F - F.transpose(-1, -2))
    real_dim = F.shape[-1]
    complex_dim = real_dim // 2
    top_form_coeff = (
        math.factorial(complex_dim)
        * _pfaffian(F)
        / ((2.0 * math.pi) ** complex_dim)
    )
    sqrt_det_g = torch.linalg.det(g).clamp(min=1e-30).sqrt()
    return orientation_sign * top_form_coeff / sqrt_det_g


def _sphere_volume(dim: int) -> float:
    """Volume of the unit round sphere S^dim."""
    return 2.0 * (math.pi ** ((dim + 1) / 2.0)) / math.gamma((dim + 1) / 2.0)


def _round_sphere_riemann(g: torch.Tensor) -> torch.Tensor:
    """
    Riemann curvature R^a_{bkl} for a unit round sphere in local coordinates.

    The sign convention is calibrated so the standard oriented ACS on S² gives
    c₁[S²] = 2 with the projected-connection formula above.
    """
    batch, d, _ = g.shape
    eye = torch.eye(d, dtype=g.dtype, device=g.device)
    riemann = g.new_zeros(batch, d, d, d, d)
    for a in range(d):
        for b in range(d):
            for k in range(d):
                for l in range(d):
                    riemann[:, a, b, k, l] = (
                        eye[a, l] * g[:, k, b] - eye[a, k] * g[:, l, b]
                    )
    return riemann


def _round_product_riemann(
    g: torch.Tensor,
    dim_alpha: int,
    dim_beta: int,
) -> torch.Tensor:
    """Block-diagonal round Riemann curvature for S^alpha x S^beta."""
    batch, d, _ = g.shape
    riemann = g.new_zeros(batch, d, d, d, d)
    eye = torch.eye(d, dtype=g.dtype, device=g.device)

    for start, block_dim in ((0, dim_alpha), (dim_alpha, dim_beta)):
        stop = start + block_dim
        for a in range(start, stop):
            for b in range(start, stop):
                for k in range(start, stop):
                    for l in range(start, stop):
                        riemann[:, a, b, k, l] = (
                            eye[a, l] * g[:, k, b] - eye[a, k] * g[:, l, b]
                        )
    return riemann


def _local_integrand_for_chart(
    subnet: nn.Module,
    x: torch.Tensor,
    g: torch.Tensor,
    gamma: torch.Tensor,
    riemann: torch.Tensor,
    orientation_sign: float,
) -> torch.Tensor:
    """Compute the first-Chern integrand in one local chart."""
    with torch.enable_grad():
        A, dA = _compute_A_and_dA(subnet, x)
        exp_a = _matrix_exp(A, getattr(subnet, 'matrix_exp_terms', None))
        exp_neg_a = _matrix_exp(-A, getattr(subnet, 'matrix_exp_terms', None))
        J = exp_a @ subnet.J0 @ exp_neg_a
        dJ = _compute_dJ_from_dA(A, dA, subnet.J0)
        dJ_cov = _covariant_dJ(J, dJ, gamma)
    return _first_chern_integrand(
        J.detach(),
        dJ_cov.detach(),
        g.detach(),
        riemann.detach(),
        orientation_sign,
    )


def _accumulate_integrand(
    subnet: nn.Module,
    x: torch.Tensor,
    g: torch.Tensor,
    gamma: torch.Tensor,
    riemann: torch.Tensor,
    orientation_sign: float,
    chunk_size: int,
) -> Tuple[float, int]:
    """Sum local first-Chern integrand values in chunks."""
    if x.shape[0] == 0:
        return 0.0, 0

    total = 0.0
    count = 0
    for start in range(0, x.shape[0], chunk_size):
        end = min(start + chunk_size, x.shape[0])
        vals = _local_integrand_for_chart(
            subnet,
            x[start:end],
            g[start:end],
            gamma[start:end],
            riemann[start:end],
            orientation_sign,
        )
        total += vals.double().sum().item()
        count += vals.numel()
    return total, count


def compute_first_chern_number(
    model: nn.Module,
    manifold_type: str,
    *,
    manifold_dim: int | None = None,
    dim_alpha: int | None = None,
    dim_beta: int | None = None,
    x: torch.Tensor | None = None,
    x_alpha: torch.Tensor | None = None,
    x_beta: torch.Tensor | None = None,
    chunk_size: int = 1024,
) -> float:
    """
    Estimate one global first Chern number from the held-out test sample.

    For `sphere`, `x` is the test sample in patch-0 ball coordinates.
    For `sphere_product`, `x_alpha` and `x_beta` are the product test sample
    in patch-(0,0) coordinates.  The reported scalar is
    ∫ c₁^n over the full manifold, where n is the complex dimension.
    """
    if manifold_type == 'sphere':
        if manifold_dim is None:
            raise ValueError("sphere Chern number requires manifold_dim")
        if x is None:
            raise ValueError("sphere Chern number requires x")

        x = x.detach()
        x2 = EvenSphere.patch_change_coords(x)
        mask0 = x.norm(dim=1) <= RADIAL_MIDPOINT
        total = 0.0
        count = 0

        for subnet, x_chart, mask in (
            (model.patch1, x, mask0),
            (model.patch2, x2, ~mask0),
        ):
            if bool(mask.any()):
                g, _, gamma = _sphere_metric_connection(x_chart[mask])
                riemann = _round_sphere_riemann(g)
                orientation_sign = (
                    1.0
                    if subnet is model.patch1
                    else float((-1) ** (manifold_dim - 1))
                )
                subtotal, n_vals = _accumulate_integrand(
                    subnet,
                    x_chart[mask],
                    g,
                    gamma,
                    riemann,
                    orientation_sign,
                    chunk_size,
                )
                total += subtotal
                count += n_vals

        if count == 0:
            raise ValueError("cannot compute Chern number from an empty test sample")
        return _sphere_volume(manifold_dim) * total / count

    if manifold_type == 'sphere_product':
        if dim_alpha is None or dim_beta is None:
            raise ValueError("sphere_product Chern number requires dim_alpha and dim_beta")
        if x_alpha is None or x_beta is None:
            raise ValueError(
                "sphere_product Chern number requires x_alpha and x_beta"
            )

        x_alpha = x_alpha.detach()
        x_beta = x_beta.detach()
        xa2 = EvenSphere.patch_change_coords(x_alpha)
        xb2 = EvenSphere.patch_change_coords(x_beta)
        mask_a0 = x_alpha.norm(dim=1) <= RADIAL_MIDPOINT
        mask_b0 = x_beta.norm(dim=1) <= RADIAL_MIDPOINT

        total = 0.0
        count = 0
        for ia, ma in ((0, mask_a0), (1, ~mask_a0)):
            xa_chart = x_alpha if ia == 0 else xa2
            for ib, mb in ((0, mask_b0), (1, ~mask_b0)):
                mask = ma & mb
                if not bool(mask.any()):
                    continue

                xb_chart = x_beta if ib == 0 else xb2
                x_full = torch.cat([xa_chart[mask], xb_chart[mask]], dim=1)
                g, _, gamma = _product_metric_connection(
                    xa_chart[mask], xb_chart[mask],
                )
                riemann = _round_product_riemann(g, dim_alpha, dim_beta)
                orientation_sign = float(
                    (-1) ** (ia * (dim_alpha - 1) + ib * (dim_beta - 1))
                )
                subtotal, n_vals = _accumulate_integrand(
                    model.patches[f'patch_{ia}{ib}'],
                    x_full,
                    g,
                    gamma,
                    riemann,
                    orientation_sign,
                    chunk_size,
                )
                total += subtotal
                count += n_vals

        if count == 0:
            raise ValueError("cannot compute Chern number from an empty test sample")
        volume = _sphere_volume(dim_alpha) * _sphere_volume(dim_beta)
        return volume * total / count

    raise ValueError(f"Unknown manifold_type: {manifold_type}")


# ─────────────────────────────────────────────────────────────────────────────
# Quick self-test
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    """
    Self-test on constant J = J₀.

    For A = 0, J is constant in coordinates so dJ = 0 and the local Berry
    curvature representative computed from dJ vanishes pointwise.
    """
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

    sphere = EvenSphere(dimension=2)
    x_test = sphere.sample_points_patch1(2000)
    c1 = compute_first_chern_number(
        model,
        manifold_type='sphere',
        manifold_dim=2,
        x=x_test,
    )
    print("Self-test: first Chern number with constant J = J₀")
    print(f"  c1: {c1:.4f}")
    print("  Expected near 2.0 for the standard S² ACS")
