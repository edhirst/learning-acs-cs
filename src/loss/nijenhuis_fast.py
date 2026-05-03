"""
Fast Nijenhuis tensor computation — drop-in replacements for the reference
classes in src.loss.sphere_loss and src.loss.product_loss.

Three optimisations vs the reference (mathematically identical results):

  1. Vectorised forward-mode AD with torch.func.jacfwd + vmap.
     Replaces the Python `for l in range(d):` loop over tangent directions
     with a single vmapped call, batching the d directions into one forward
     pass. Eliminates d-1 kernel-launch + Python-dispatch round trips.

  2. Antisymmetry of N in (i, j) used to fuse 3 einsums into 2.
     N^k_{ij} = J^l_i ∂_l J^k_j − J^l_j ∂_l J^k_i − J^k_l (∂_i J^l_j − ∂_j J^l_i)
     is antisymmetric in (i, j); writing
         T = J^l_i ∂_l J^k_j,   S = J^k_l ∂_i J^l_j,
     gives  N = (T - T^⊤_{ij}) - (S - S^⊤_{ij}),
     so the assembly takes 2 contractions + 2 transposes instead of 3
     contractions + 1 permute.  Saves ~33% of post-AD einsum FLOPs.

  3. Optional torch.compile of the einsum kernel (use_compile=True).
     After (1), the assembly is functional and AOTAutograd-compilable.

Public API:

    NijenhuisTensorLossFast        -- replaces sphere_loss.NijenhuisTensorLoss
    ProductNijenhuisLossFast       -- replaces product_loss.ProductNijenhuisLoss
    compute_dJ_jacfwd(subnet, x)   -- raw kernel, returns dJ
    nij_tensor_antisym(J, dJ)      -- raw kernel, returns N
"""

import torch
import torch.func as func
import torch.nn as nn
from typing import Dict, Optional, Tuple

from src.geometry import EvenSphere, RADIAL_MIDPOINT
from src.geometry.sphere_product import SphereProduct
from src.loss.sphere_loss import _metric_norm_sq_12tensor


# ─────────────────────────────────────────────────────────────────────────────
# Pure kernel functions
# ─────────────────────────────────────────────────────────────────────────────


def compute_dJ_jacfwd(subnet: nn.Module, x: torch.Tensor) -> torch.Tensor:
    """
    Compute dJ[b, k, j, l] = ∂_{x^l} J^k_j via vmapped forward-mode AD.

    Equivalent to the reference loop

        for l in range(d):
            with fwdAD.dual_level():
                dJ_l = fwdAD.unpack_dual(subnet(fwdAD.make_dual(x, e_l))).tangent
        dJ = torch.stack(dJ_l, dim=-1)

    but fuses the d directions into a single vmapped pass.  The autograd
    graph back to the subnet's parameters is preserved, so loss.backward()
    works as in the reference.

    Args:
        subnet: callable producing J of shape (batch, d, d) from (batch, d).
        x:      ball coordinates, shape (batch, d).

    Returns:
        dJ:     shape (batch, d, d, d).  Last index is the partial direction,
                matching the reference convention.
    """
    def _fn_single(x_single: torch.Tensor) -> torch.Tensor:
        # subnet expects (batch, d); call with batch=1 then squeeze.
        return subnet(x_single.unsqueeze(0)).squeeze(0)

    # jacfwd: (d,) -> (d, d, d).  vmap over batch.
    return func.vmap(func.jacfwd(_fn_single))(x)


def nij_tensor_antisym(J: torch.Tensor, dJ: torch.Tensor) -> torch.Tensor:
    """
    Assemble the Nijenhuis tensor with antisymmetry-fused einsums.

    Reference uses 3 einsums:
        term1 = einsum('bli,bkjl->bkij', J, dJ)
        term2 = -einsum('blj,bkil->bkij', J, dJ)
        anti  = dJ.permute(0, 1, 3, 2) - dJ
        term3 = einsum('bkl,blij->bkij', J, anti)
        N = term1 + term2 - term3

    Optimised: 2 einsums + 2 transposes:
        T = einsum('bli,bkjl->bkij', J, dJ)       # = J^l_i ∂_l J^k_j
        S = einsum('bkl,blji->bkij', J, dJ)       # = J^k_l ∂_i J^l_j
        N = (T - T^⊤_{ij}) - (S - S^⊤_{ij})

    Correctness: term2 = -T^⊤_{ij} since term2[k,i,j] = -J^l_j ∂_l J^k_i = -T[k,j,i].
    For term3, expand the (∂_i J^l_j - ∂_j J^l_i) factor to S - S^⊤_{ij}.

    Args:
        J:  (batch, d, d) almost complex structure values.
        dJ: (batch, d, d, d) partials, dJ[b, k, j, l] = ∂_{x^l} J^k_j.

    Returns:
        N:  (batch, d, d, d), N[b, k, i, j] = N^k_{ij}.
    """
    T = torch.einsum('bli,bkjl->bkij', J, dJ)   # T[b,k,i,j] = J^l_i ∂_l J^k_j
    S = torch.einsum('bkl,blji->bkij', J, dJ)   # S[b,k,i,j] = J^k_l ∂_i J^l_j
    return (T - T.transpose(-2, -1)) - (S - S.transpose(-2, -1))


# A torch.compile'd variant of the kernel.  Lazy-initialised on first use so
# importing this module doesn't pay compile cost up-front.
_compiled_kernel = None


def _get_compiled_kernel():
    global _compiled_kernel
    if _compiled_kernel is None:
        _compiled_kernel = torch.compile(nij_tensor_antisym, dynamic=False)
    return _compiled_kernel


# ─────────────────────────────────────────────────────────────────────────────
# nn.Module wrappers — drop-in replacements for the reference classes
# ─────────────────────────────────────────────────────────────────────────────


class NijenhuisTensorLossFast(nn.Module):
    """
    Drop-in replacement for src.loss.sphere_loss.NijenhuisTensorLoss.

    Same constructor (with one extra optional `use_compile` flag) and the
    same forward(x, subnet, hard_cutoff) signature.  Returns the same
    (loss, info) pair.
    """

    def __init__(
        self,
        manifold_dim: int,
        overlap_size: float = 0.1,
        use_round_metric_contraction: bool = False,
        use_compile: bool = False,
    ):
        super().__init__()
        self.manifold_dim = manifold_dim
        self.overlap_size = overlap_size
        self.use_round_metric_contraction = use_round_metric_contraction
        self.use_compile = use_compile

    def _kernel(self, J: torch.Tensor, dJ: torch.Tensor) -> torch.Tensor:
        if self.use_compile:
            return _get_compiled_kernel()(J, dJ)
        return nij_tensor_antisym(J, dJ)

    def forward(
        self,
        x: torch.Tensor,
        subnet: nn.Module,
        hard_cutoff: bool = False,
    ) -> Tuple[torch.Tensor, Dict]:
        with torch.enable_grad():
            J = subnet(x)              # (batch, d, d)
        batch, d, _ = J.shape

        # Vectorised forward-mode AD: one vmapped call instead of d sequential.
        # x.detach() matches the reference (the loss differentiates J w.r.t.
        # parameters, not w.r.t. x).
        x_d = x.detach()
        dJ = compute_dJ_jacfwd(subnet, x_d)        # (batch, d, d, d)

        # Antisymmetry-fused assembly.
        N = self._kernel(J, dJ)

        if self.use_round_metric_contraction:
            g, g_inv = EvenSphere.round_metric_ball_coords(x_d)
            nij_sq = _metric_norm_sq_12tensor(N, g, g_inv)
        else:
            nij_sq = (N ** 2).sum(dim=(1, 2, 3))
        nij_norm = torch.sqrt(nij_sq + 1e-20)

        # Radial weight (identical to reference).
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


class ProductNijenhuisLossFast(nn.Module):
    """
    Drop-in replacement for src.loss.product_loss.ProductNijenhuisLoss.

    Same constructor (with one extra optional `use_compile` flag) and the
    same forward(x_alpha, x_beta, x_a2, x_b2, model, hard_cutoff) signature.
    """

    def __init__(
        self,
        total_dim: int,
        overlap_size_alpha: float = 0.1,
        overlap_size_beta: float = 0.1,
        use_round_metric_contraction: bool = False,
        use_compile: bool = False,
    ):
        super().__init__()
        self.total_dim = total_dim
        self.r_cut_alpha = RADIAL_MIDPOINT + overlap_size_alpha
        self.r_cut_beta  = RADIAL_MIDPOINT + overlap_size_beta
        self.overlap_size_alpha = overlap_size_alpha
        self.overlap_size_beta  = overlap_size_beta
        self.use_round_metric_contraction = use_round_metric_contraction
        self.use_compile = use_compile

    def _kernel(self, J: torch.Tensor, dJ: torch.Tensor) -> torch.Tensor:
        if self.use_compile:
            return _get_compiled_kernel()(J, dJ)
        return nij_tensor_antisym(J, dJ)

    def _radial_weight(
        self,
        rho_alpha: torch.Tensor,
        rho_beta: torch.Tensor,
        hard_cutoff: bool,
    ) -> torch.Tensor:
        if hard_cutoff:
            wa = (rho_alpha <= self.r_cut_alpha).float()
            wb = (rho_beta  <= self.r_cut_beta ).float()
        else:
            exc_a = (rho_alpha - self.r_cut_alpha).clamp(min=0.0)
            exc_b = (rho_beta  - self.r_cut_beta ).clamp(min=0.0)
            wa = torch.exp(-(exc_a / self.overlap_size_alpha) ** 2)
            wb = torch.exp(-(exc_b / self.overlap_size_beta ) ** 2)
        return wa * wb

    def _nij_with_weight(
        self,
        x: torch.Tensor,
        subnet: nn.Module,
        w: torch.Tensor,
        hard_cutoff: bool,
        g: Optional[torch.Tensor] = None,
        g_inv: Optional[torch.Tensor] = None,
    ) -> Tuple[torch.Tensor, Dict]:
        with torch.enable_grad():
            J = subnet(x)
        x_d = x.detach()
        dJ = compute_dJ_jacfwd(subnet, x_d)
        N = self._kernel(J, dJ)

        if g is not None and g_inv is not None:
            nij_norm = torch.sqrt(_metric_norm_sq_12tensor(N, g, g_inv) + 1e-20)
        else:
            nij_norm = torch.sqrt((N ** 2).sum(dim=(1, 2, 3)) + 1e-20)

        w_sum = w.sum().clamp(min=1e-8)
        loss = (w * nij_norm).sum() / w_sum
        return loss, {'nij_norm': loss.item()}

    def forward(
        self,
        x_alpha:   torch.Tensor,
        x_beta:    torch.Tensor,
        x_a2:      torch.Tensor,
        x_b2:      torch.Tensor,
        model:     nn.Module,
        hard_cutoff: bool = False,
    ) -> Tuple[torch.Tensor, Dict]:
        patch_alpha = {0: x_alpha, 1: x_a2}
        patch_beta  = {0: x_beta,  1: x_b2}

        total = x_alpha.sum() * 0.0
        info: Dict = {}

        for ia, ib in [(0, 0), (0, 1), (1, 0), (1, 1)]:
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

            nij_l, _ = self._nij_with_weight(
                x_combined, subnet, w, hard_cutoff, g, g_inv,
            )

            info[f'nij_{ia}{ib}'] = nij_l.item()
            total = total + nij_l

        return total, info
