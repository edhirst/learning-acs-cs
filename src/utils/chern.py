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

First Chern number (integrated):

    c₁[M] = (1/2π) ∫_M F ∧ ω^{n-1}/(n-1)!

For real dimension 2 this is the usual integral of the first Chern form.  In
higher real dimension this module reports the standard Chern-Weil pairing with
the compatible fundamental form,

    ∫_M c₁ ∧ ω^{n-1}/(n-1)!,

estimated globally from the held-out test sample.  For sphere products this is
computed on the full product manifold, not on individual sphere factors.

Reference cases:
    S²             :  c₁ = 2
    S²n (n ≥ 2)    :  c₁ class = 0   (H²(S^{2n}) = 0)
    S¹×S¹ (T²)     :  c₁ = 0   (flat)
    higher products:  reported as the full-manifold Chern-Weil pairing above

Implementation notes
--------------------
The reported scalar is computed only for final test/evaluation reporting; no
Chern-related quantity is tracked during training or validation.  In real
dimension 2 it is the usual first Chern number.  In higher real dimension it is
the full-manifold Chern-Weil pairing with the compatible form.

Both `_compute_A_and_dA` and `_compute_dJ_from_dA` from `src.loss.sphere_loss`
are reused to obtain ∂J/∂x via the existing AD path.
"""

import math
import torch
import torch.nn as nn
from typing import Tuple

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


def _first_chern_integrand(
    J: torch.Tensor,
    dJ: torch.Tensor,
    g: torch.Tensor,
    g_inv: torch.Tensor,
    riemann: torch.Tensor,
) -> torch.Tensor:
    """
    Pointwise scalar whose volume integral gives the first Chern number.

    The Berry curvature F is a 2-form.  On a 2n-dimensional almost Hermitian
    manifold, the scalar paired with the volume form is

        (1 / 2π) ΛF,    where    F ∧ ω^(n-1)/(n-1)! = (ΛF) dVol_g.

    With this codebase's J convention, ΛF = -1/2 F_ij ω^ij.  For real
    dimension 2 this reduces to the usual coordinate integral of F_01.
    """
    F = chern_curvature(J, dJ, riemann)                             # (batch, d, d)
    omega = torch.einsum('bik,bkj->bij', g, J)                      # ω_ij = g_ik J^k_j
    omega_up = torch.einsum('bip,bjq,bpq->bij', g_inv, g_inv, omega)
    lambda_f = -0.5 * torch.einsum('bij,bij->b', F, omega_up)
    return (1.0 / (2.0 * math.pi)) * lambda_f


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
    g_inv: torch.Tensor,
    riemann: torch.Tensor,
) -> torch.Tensor:
    """Compute the first-Chern integrand in one local chart."""
    with torch.enable_grad():
        A, dA = _compute_A_and_dA(subnet, x)
        exp_a = _matrix_exp(A, getattr(subnet, 'matrix_exp_terms', None))
        exp_neg_a = _matrix_exp(-A, getattr(subnet, 'matrix_exp_terms', None))
        J = exp_a @ subnet.J0 @ exp_neg_a
        dJ = _compute_dJ_from_dA(A, dA, subnet.J0)
    return _first_chern_integrand(
        J.detach(),
        dJ.detach(),
        g.detach(),
        g_inv.detach(),
        riemann.detach(),
    )


def _accumulate_integrand(
    subnet: nn.Module,
    x: torch.Tensor,
    g: torch.Tensor,
    g_inv: torch.Tensor,
    riemann: torch.Tensor,
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
            g_inv[start:end],
            riemann[start:end],
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
    ∫ c₁ ∧ ω^(n-1)/(n-1)! over the full manifold.
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
                g, g_inv = EvenSphere.round_metric_ball_coords(x_chart[mask])
                riemann = _round_sphere_riemann(g)
                subtotal, n_vals = _accumulate_integrand(
                    subnet, x_chart[mask], g, g_inv, riemann, chunk_size,
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
                g, g_inv = SphereProduct.round_metric_ball_coords(
                    xa_chart[mask], xb_chart[mask],
                )
                riemann = _round_product_riemann(g, dim_alpha, dim_beta)
                subtotal, n_vals = _accumulate_integrand(
                    model.patches[f'patch_{ia}{ib}'],
                    x_full,
                    g,
                    g_inv,
                    riemann,
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
    print("  Expected near 0.0 for constant J in this gauge")
