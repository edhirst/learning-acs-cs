"""
Compiled (graph-cached) Nijenhuis tensor loss.

Goal
----
Trace the forward + backward graph once via AOTAutograd / torch.compile and
reuse the cached compiled artefacts on every training step, eliminating the
per-step Python-level autograd-graph rebuild — i.e. the user's hypothesis
of "compute the derivative path once, cache it."

Reality on PyTorch 2.10
-----------------------
The full Inductor backend (``mode='default' | 'reduce-overhead' |
'max-autotune' | 'max-autotune-no-cudagraphs'``) currently fails when AOTAutograd
tries to compile the *backward* of any forward-mode AD trace
(``torch.func.jacfwd`` / ``torch.func.jvp`` / ``torch.autograd.forward_ad``):

    InductorError: RuntimeError: ZeroTensor is not serializable
    RuntimeError: ZeroTensors are immutable

This appears to be an interaction between Inductor's graph-repro serialization
and PyTorch's internal ZeroTensor optimization for unused gradient inputs.
Reproduces with both ``vmap(jacfwd(...))`` and a manual ``for l in range(d):
torch.autograd.forward_ad.dual_level()`` loop, with or without ``x.detach()``.

Working compile path
--------------------
``backend='aot_eager'`` traces the forward + backward graphs into FX *without*
invoking Inductor codegen.  The graph-rebuild overhead is removed, but the
resulting kernels are still eager dispatches — no kernel fusion gains.

For our workload (matmul-bound: each step is dominated by ~7 evaluations of
the MLP under vmap+jacfwd), the per-op autograd overhead is a small fraction
of total wall-clock, so ``aot_eager`` delivers only modest extra speedup
(~0–5 %) on top of nijenhuis_fast.NijenhuisTensorLossFast.  That's an honest
result, not a bug in this module.

Path forward
------------
1. Wait for the Inductor backward compile bug to be fixed upstream
   (filed in PyTorch issue tracker as a recurring concern).  When it is, the
   ``mode='default'`` path here will activate kernel fusion automatically.
2. Implement a hand-coded JVP for ``PatchSubNet`` layer-by-layer to bypass
   forward-mode AD entirely.  Significant code, but compile-friendly.
3. Implement a ``torch.autograd.Function`` with manual mixed-mode backward.
   The forward computes loss in ``no_grad``, and the backward calls a
   functional reverse-mode + forward-mode chain.  Works around the bug at
   the cost of considerable code.

Public API
----------
NijenhuisTensorLossCompiled — drop-in replacement for the reference loss.
                              Selects compile strategy via ``mode`` argument.
"""

import warnings
from typing import Dict, Optional, Tuple

import torch
import torch.func as func
import torch.nn as nn

from src.geometry import EvenSphere, RADIAL_MIDPOINT
from src.loss.nijenhuis_fast import (
    compute_dJ_jacfwd,
    nij_tensor_antisym,
    NijenhuisTensorLossFast,
)
from src.loss.sphere_loss import _metric_norm_sq_12tensor


# ─────────────────────────────────────────────────────────────────────────────
# Compile strategies
# ─────────────────────────────────────────────────────────────────────────────


# Modes that delegate to torch.compile under different backends.
_COMPILE_MODES = {
    'aot_eager':    dict(backend='aot_eager',  dynamic=False),
    'inductor':     dict(backend='inductor',   dynamic=False),
    'reduce-overhead': dict(backend='inductor', dynamic=False, mode='reduce-overhead'),
}


def _build_loss_fn(subnet: nn.Module, manifold_dim: int, overlap_size: float):
    """Return a closure that takes x and returns the (scalar) Nijenhuis loss.

    Captured by closure so torch.compile can inline subnet's forward and
    cache the trace per-subnet.
    """
    r_cut = float(RADIAL_MIDPOINT + overlap_size)

    def fn_single(x_single: torch.Tensor) -> torch.Tensor:
        return subnet(x_single.unsqueeze(0)).squeeze(0)

    def loss_fn(x: torch.Tensor) -> torch.Tensor:
        J = subnet(x)
        dJ = func.vmap(func.jacfwd(fn_single))(x)
        T = torch.einsum('bli,bkjl->bkij', J, dJ)
        S = torch.einsum('bkl,blji->bkij', J, dJ)
        N = (T - T.transpose(-2, -1)) - (S - S.transpose(-2, -1))
        nij_sq = (N ** 2).sum(dim=(1, 2, 3))
        nij_norm = torch.sqrt(nij_sq + 1e-20)

        rho = x.detach().norm(dim=1)
        excess = (rho - r_cut).clamp(min=0.0)
        w = torch.exp(-(excess / overlap_size) ** 2)
        w_sum = w.sum().clamp(min=1e-8)
        return (w * nij_norm).sum() / w_sum

    return loss_fn


def _wrap_compile(loss_fn, mode: str):
    """Wrap loss_fn with torch.compile under the requested mode.

    `mode` is one of:
        'aot_eager'        — AOTAutograd traces fwd+bwd to FX, no Inductor codegen.
                             Currently the only Inductor-free path that works.
        'inductor'         — Full Inductor codegen (CURRENTLY BROKEN; see module
                             docstring for the ZeroTensor bug).
        'reduce-overhead'  — Inductor with CUDA-Graphs (no-op on CPU).
    """
    if mode not in _COMPILE_MODES:
        raise ValueError(
            f"Unknown compile mode {mode!r}. Available: {list(_COMPILE_MODES)}"
        )
    cfg = _COMPILE_MODES[mode]
    return torch.compile(loss_fn, **cfg)


# ─────────────────────────────────────────────────────────────────────────────
# nn.Module wrapper
# ─────────────────────────────────────────────────────────────────────────────


class NijenhuisTensorLossCompiled(nn.Module):
    """Drop-in for sphere_loss.NijenhuisTensorLoss with torch.compile-cached
    forward + backward graphs.

    Same constructor signature as the reference loss, plus:

        mode: 'aot_eager'  — recommended, currently the only working path.
              'inductor'   — broken in PyTorch 2.10 (see module docstring).
              'reduce-overhead' — same as inductor on CPU; broken on PyTorch 2.10.

    First call to forward(x, subnet) traces (forward + backward) and caches
    the FX graph.  Subsequent calls reuse it.  The cache is keyed by
    ``id(subnet)``, so distinct subnet instances get distinct compiles
    (parameter updates within the same instance do *not* invalidate).

    Limitations:
      • hard_cutoff=True falls back to NijenhuisTensorLossFast (the eval
        path is rare and the radial-mask branching is a graph break).
      • use_round_metric_contraction=True falls back similarly.
      • Variable input shapes will cause recompiles.  Pin batch size to
        avoid this in training.
    """

    def __init__(
        self,
        manifold_dim: int,
        overlap_size: float = 0.1,
        use_round_metric_contraction: bool = False,
        mode: str = 'aot_eager',
    ):
        super().__init__()
        self.manifold_dim = manifold_dim
        self.overlap_size = overlap_size
        self.use_round_metric_contraction = use_round_metric_contraction
        self.mode = mode

        # Cache: id(subnet) -> compiled fn.
        self._cache: Dict[int, callable] = {}
        # Fallback path for hard_cutoff / metric contraction.
        self._fallback = NijenhuisTensorLossFast(
            manifold_dim, overlap_size,
            use_round_metric_contraction=use_round_metric_contraction,
            use_compile=False,
        )

    def _get_compiled(self, subnet: nn.Module):
        key = id(subnet)
        if key not in self._cache:
            loss_fn = _build_loss_fn(subnet, self.manifold_dim, self.overlap_size)
            self._cache[key] = _wrap_compile(loss_fn, self.mode)
        return self._cache[key]

    def forward(
        self,
        x: torch.Tensor,
        subnet: nn.Module,
        hard_cutoff: bool = False,
    ) -> Tuple[torch.Tensor, Dict]:
        # Branches that aren't in the compile graph.
        if hard_cutoff or self.use_round_metric_contraction:
            return self._fallback(x, subnet, hard_cutoff)

        compiled_fn = self._get_compiled(subnet)
        loss = compiled_fn(x)
        return loss, {'nij_norm': loss.item()}
