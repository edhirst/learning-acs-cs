"""
Decomposed Nijenhuis tensor computation: structurally exploit
``J(x) = J_0 · exp(-2 A(x))``.

Mathematical justification
--------------------------
The architecture parametrises ``J(x) = exp(A(x)) · J_0 · exp(-A(x))`` with
``A(x)`` constrained to the anti-commuting subspace
``𝒜 = {A : A J_0 + J_0 A = 0}``.  Anti-commutation gives
``exp(A) · J_0 = J_0 · exp(-A)``, so
::

    J = exp(A) · J_0 · exp(-A) = J_0 · exp(-2 A).

Two matrix exponentials collapse to one.

Why this matters
----------------
The dominant cost in nijenhuis_fast.py is the vmap+jacfwd trace through
the *whole* subnet — which includes ``matrix_exp`` *twice* per evaluation.
With the simplification:

  • Forward needs only **one** matrix_exp.
  • The MLP that produces ``A`` is matrix_exp-free.  We can therefore
    differentiate ``A`` with vmap+jacfwd over the **MLP backbone only**,
    avoiding the cost of carrying matrix_exp's primal+tangent through
    every direction.
  • The promotion ``∂A → ∂J`` is computed as a single matrix_exp JVP per
    direction (or vmapped over directions), at the same primal point ``-2 A``.
    This decouples the matrix_exp's JVP from the larger autograd trace.

Net cost (compared to FAST):
  • FAST:     ``d × subnet`` evaluations (each = MLP + 2 × matrix_exp), with
              dual-number propagation through every op.
  • DECOMP:   ``d × MLP`` evaluations (with dual numbers, MLP only)
              + ``d × matrix_exp_jvp`` (separate, at the same primal A)
              + ``1 × matrix_exp`` (forward).

For ``d ≥ 4`` the saving is roughly ``(d-1) × matrix_exp_cost`` per step.

Public API
----------
NijenhuisTensorLossDecomposed — drop-in replacement for the reference loss.
"""

import torch
import torch.func as func
import torch.nn as nn
from typing import Dict, Tuple

from src.geometry import EvenSphere, RADIAL_MIDPOINT
from src.loss.sphere_loss import _metric_norm_sq_12tensor


# ─────────────────────────────────────────────────────────────────────────────
# Pure kernel functions
# ─────────────────────────────────────────────────────────────────────────────


def compute_A_and_dA_jacfwd(subnet: nn.Module, x: torch.Tensor):
    """
    Compute (A, dA) where:
        A:  shape (batch, d, d)         — values of A(x).
        dA: shape (batch, d, d, d_in)   — partial derivatives, dA[..., l] = ∂A/∂x_l.

    Differentiation is done over ``subnet.forward_raw`` (the MLP path producing
    A directly, without the matrix_exp closure).  This is significantly
    cheaper than vmap+jacfwd over the full subnet.
    """
    def mlp_single(x_single: torch.Tensor) -> torch.Tensor:
        return subnet.forward_raw(x_single.unsqueeze(0)).squeeze(0)

    A = subnet.forward_raw(x)                         # (b, d, d)
    dA = func.vmap(func.jacfwd(mlp_single))(x)        # (b, d, d, d_in); last index is direction
    return A, dA


def compute_dJ_from_dA(A: torch.Tensor, dA: torch.Tensor, J0: torch.Tensor) -> torch.Tensor:
    """
    Compute dJ from (A, dA) using the chain rule on J = J_0 · exp(-2 A):

        dJ_l = J_0 · D[exp(-2 A)] (-2 dA_l)

    The Frechet derivative D[exp(-2 A)] is evaluated via a per-direction
    ``torch.func.jvp`` on torch.matrix_exp at the same primal ``-2 A``, vmapped
    over the directions.

    Args:
        A:  (batch, d, d)
        dA: (batch, d, d, d_in)
        J0: (d, d) reference complex structure (already signed for this patch).

    Returns:
        dJ: (batch, d, d, d_in)
    """
    def expm_minus_2A(A_):
        return torch.matrix_exp(-2 * A_)

    def jvp_one_direction(A_, dA_l):
        # primal: exp(-2 A);  tangent: D[exp(-2 A)] (-2 dA_l)
        _, tangent = func.jvp(expm_minus_2A, (A_,), (dA_l,))
        return tangent  # (b, d, d)

    # vmap over the direction axis of dA (the last axis).
    dA_perm = dA.movedim(-1, 0)                                            # (d_in, b, d, d)
    dexpB_perm = func.vmap(jvp_one_direction, in_dims=(None, 0))(A, dA_perm)  # (d_in, b, d, d)
    dexpB = dexpB_perm.movedim(0, -1)                                      # (b, d, d, d_in)

    # Multiply by J0 on the left over the row index of the matrix.
    dJ = torch.einsum('kp,bpjl->bkjl', J0, dexpB)
    return dJ


def nij_tensor_antisym(J: torch.Tensor, dJ: torch.Tensor) -> torch.Tensor:
    """
    Antisymmetry-fused Nijenhuis assembly (same as nijenhuis_fast.py).

        T = einsum('bli,bkjl->bkij', J, dJ)
        S = einsum('bkl,blji->bkij', J, dJ)
        N = (T - T^T_{ij}) - (S - S^T_{ij})
    """
    T = torch.einsum('bli,bkjl->bkij', J, dJ)
    S = torch.einsum('bkl,blji->bkij', J, dJ)
    return (T - T.transpose(-2, -1)) - (S - S.transpose(-2, -1))


# ─────────────────────────────────────────────────────────────────────────────
# nn.Module wrapper — drop-in replacement for the reference loss
# ─────────────────────────────────────────────────────────────────────────────


class NijenhuisTensorLossDecomposed(nn.Module):
    """
    Drop-in for sphere_loss.NijenhuisTensorLoss using
    ``J = J_0 · exp(-2 A)`` and a split AD path.

    Same constructor + forward(x, subnet, hard_cutoff) signature; same
    return type ``(loss, info)``.
    """

    def __init__(
        self,
        manifold_dim: int,
        overlap_size: float = 0.1,
        use_round_metric_contraction: bool = False,
    ):
        super().__init__()
        self.manifold_dim = manifold_dim
        self.overlap_size = overlap_size
        self.use_round_metric_contraction = use_round_metric_contraction

    def forward(
        self,
        x: torch.Tensor,
        subnet: nn.Module,
        hard_cutoff: bool = False,
    ) -> Tuple[torch.Tensor, Dict]:
        # Step 1: A and dA via cheap MLP-only AD.
        A, dA = compute_A_and_dA_jacfwd(subnet, x)         # (b, d, d), (b, d, d, d)

        # Step 2: J via the simplified product (one matrix_exp).
        J0 = subnet.J0                                      # (d, d), per-patch signed
        expB = torch.matrix_exp(-2 * A)                     # (b, d, d)
        J = torch.einsum('kp,bpj->bkj', J0, expB)           # (b, d, d)

        # Step 3: dJ via per-direction matrix_exp JVP at the same primal.
        dJ = compute_dJ_from_dA(A, dA, J0)                  # (b, d, d, d)

        # Step 4: Nijenhuis assembly + norm + radial weight.
        N = nij_tensor_antisym(J, dJ)
        if self.use_round_metric_contraction:
            x_d = x.detach()
            g, g_inv = EvenSphere.round_metric_ball_coords(x_d)
            nij_sq = _metric_norm_sq_12tensor(N, g, g_inv)
        else:
            nij_sq = (N ** 2).sum(dim=(1, 2, 3))
        nij_norm = torch.sqrt(nij_sq + 1e-20)

        rho = x.detach().norm(dim=1)
        r_cut = RADIAL_MIDPOINT + self.overlap_size
        if hard_cutoff:
            w = (rho <= r_cut).float()
        else:
            excess = (rho - r_cut).clamp(min=0.0)
            w = torch.exp(-(excess / self.overlap_size) ** 2)

        w_sum = w.sum().clamp(min=1e-8)
        loss = (w * nij_norm).sum() / w_sum
        return loss, {'nij_norm': loss.item()}
