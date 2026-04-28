"""
Full-matrix Gauss-Newton optimizer for physics-informed neural networks.

Implements the rank-1 EMA estimator of the Gauss-Newton (GN) matrix described
in Wang et al. (2025) "Discovery of Unstable Singularities" (arXiv:2509.14185),
§ "Second-order Optimizer".

Algorithm (one step):
  1. Compute residuals r (flat, shape m) via residual_fn().
  2. Compute gradient g = Jᵀ r via VJP with adjoint r.
  3. Sample z ~ N(0, Iₘ); compute u = Jᵀ z via second VJP.
  4. EMA update:  Ĝ ← β Ĝ + (1−β) outer(u, u).
  5. Solve (Ĝ + γ I) δ = g.
  6. Auto-LR (Martens-Grosse): α* = (g · δ) / (δᵀ (Ĝ + γ I) δ).
  7. Update: θ ← θ − α* δ.
  8. [if adaptive_damping] Compute ρ = actual / predicted decrease; adjust γ.

This is feasible only for small networks (≲ tens of thousands of parameters)
because the GN matrix is stored and inverted as a dense n × n matrix.
"""

from __future__ import annotations

import math
from typing import Callable, Dict, Iterable, List, Optional, Tuple

import torch


class GaussNewtonOptimizer:
    """Full-matrix Gauss-Newton optimizer with rank-1 EMA estimator.

    Args:
        params: Iterable of model parameters to optimise.
        damping: Initial Levenberg-Marquardt damping coefficient γ.
        ema_decay: Exponential moving-average decay β for the GN matrix.
        adaptive_damping: If True, apply Levenberg-Marquardt damping
            adjustment based on the ratio of actual to predicted decrease.
        damping_up_factor: Multiply γ by this factor when ρ < rho_threshold_bad.
        damping_down_factor: Divide γ by this factor when ρ > rho_threshold_good.
        rho_threshold_bad: ρ below this triggers damping increase; ρ < 0 also
            causes step rejection.
        rho_threshold_good: ρ above this triggers damping decrease.
    """

    def __init__(
        self,
        params: Iterable[torch.nn.Parameter],
        damping: float = 1.0,
        ema_decay: float = 0.95,
        adaptive_damping: bool = False,
        damping_up_factor: float = 3.0,
        damping_down_factor: float = 0.5,
        rho_threshold_bad: float = 0.25,
        rho_threshold_good: float = 0.75,
        max_step_norm: float = float('inf'),
    ) -> None:
        self._params: List[torch.nn.Parameter] = list(params)
        self._damping = float(damping)
        self._ema_decay = float(ema_decay)
        self._adaptive_damping = bool(adaptive_damping)
        self._damping_up_factor = float(damping_up_factor)
        self._damping_down_factor = float(damping_down_factor)
        self._rho_bad = float(rho_threshold_bad)
        self._rho_good = float(rho_threshold_good)
        self._max_step_norm = float(max_step_norm)
        self._G_ema: Optional[torch.Tensor] = None  # dense (n, n) GN matrix EMA

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    @property
    def damping(self) -> float:
        """Current damping coefficient γ."""
        return self._damping

    def zero_grad(self) -> None:
        """Clear .grad buffers (for interface compatibility; GN does not use them)."""
        for p in self._params:
            if p.grad is not None:
                p.grad.detach_()
                p.grad.zero_()

    def step(
        self,
        residual_fn: Callable[[], Tuple[torch.Tensor, Dict]],
    ) -> Dict:
        """Perform one Gauss-Newton step.

        Args:
            residual_fn: Callable that returns (r, info) where r is a flat
                residual tensor of shape (m,) with gradient tracking through
                the model parameters, and info is the standard loss info dict
                (scalar values for logging). The GN objective is ||r||².

        Returns:
            info dict from the residual_fn call.
        """
        # ---- Step 1: compute residuals --------------------------------
        r, info = residual_fn()

        params = [p for p in self._params if p.requires_grad]

        if r.numel() == 0:
            return info

        # ---- Step 2: gradient g = Jᵀ r --------------------------------
        grads_g = torch.autograd.grad(
            r, params,
            grad_outputs=r,
            retain_graph=True,
            create_graph=False,
        )
        g = self._flatten(grads_g)      # (n,)

        # Guard: skip update if residuals or gradient are non-finite.
        if not torch.isfinite(r).all() or not torch.isfinite(g).all():
            return info

        # ---- Step 3: rank-1 GN estimate u = Jᵀ z, z ~ N(0, I) --------
        z = torch.randn(r.shape[0], dtype=r.dtype, device=r.device)
        grads_z = torch.autograd.grad(
            r, params,
            grad_outputs=z,
            retain_graph=False,
            create_graph=False,
        )
        u = self._flatten(grads_z)      # (n,)

        # Guard: skip update if rank-1 sample is non-finite (e.g. exploding network).
        if not torch.isfinite(u).all():
            return info

        # ---- Step 4: EMA update of GN matrix --------------------------
        u_outer = torch.outer(u, u)
        if self._G_ema is None:
            self._G_ema = u_outer
        else:
            self._G_ema.mul_(self._ema_decay).add_((1.0 - self._ema_decay) * u_outer)

        # ---- Step 5: solve (G_ema + γI) δ = g -------------------------
        # A = G_ema + γI is positive definite for γ > 0, so Cholesky is both
        # correct and O(n³/3) — far cheaper than the O(n³) LAPACK driver used
        # by lstsq.  Fall back to lstsq only if Cholesky fails numerically.
        n = g.shape[0]
        I = torch.eye(n, dtype=g.dtype, device=g.device)
        A = self._G_ema + self._damping * I
        try:
            L = torch.linalg.cholesky(A)
            delta = torch.cholesky_solve(g.unsqueeze(-1), L).squeeze(-1)
        except torch.linalg.LinAlgError:
            delta = torch.linalg.lstsq(A, g.unsqueeze(-1)).solution.squeeze(-1)

        # Guard: skip update if solve produced non-finite values.
        if not torch.isfinite(delta).all():
            return info

        # ---- Step 6: auto-LR (Martens-Grosse) -------------------------
        Gd = self._G_ema @ delta
        numer = (g * delta).sum()
        denom = (delta * Gd).sum() + self._damping * (delta * delta).sum()
        alpha = (numer / denom.clamp(min=1e-30)).clamp(min=0.0, max=10.0)

        # ---- Adaptive damping: snapshot before update -----------------
        if self._adaptive_damping:
            loss_before = float((r.detach() ** 2).sum())
            old_data = [p.data.clone() for p in params]

        # ---- Step 7: apply update θ ← θ − α δ ------------------------
        update = alpha * delta
        update_norm = update.norm()
        if update_norm > self._max_step_norm:
            update = update * (self._max_step_norm / update_norm)
        self._apply_update(params, update)

        # ---- Step 8: Levenberg-Marquardt damping update ---------------
        if self._adaptive_damping:
            r_new, _ = residual_fn()
            loss_after = float((r_new.detach() ** 2).sum())

            actual_decrease = loss_before - loss_after
            # Predicted decrease in h = ||r||² under the quadratic model:
            # h(θ - α δ) ≈ h(θ) - 2α (Jᵀr · δ) + α² δᵀ (G+γI) δ
            # where ∇h = 2Jᵀr = 2g,  Hessian ≈ 2(G+γI).
            predicted_decrease = float(
                2.0 * alpha * numer - alpha ** 2 * denom
            )
            if abs(predicted_decrease) < 1e-30:
                rho = 1.0
            else:
                rho = actual_decrease / predicted_decrease

            if rho < self._rho_bad:
                self._damping = min(self._damping * self._damping_up_factor, 1e10)
                if rho < 0.0:
                    # Reject step: restore parameters
                    for p, old in zip(params, old_data):
                        p.data.copy_(old)
            elif rho > self._rho_good:
                self._damping = max(self._damping / self._damping_down_factor, 1e-10)

        return info

    # ------------------------------------------------------------------
    # Checkpoint support
    # ------------------------------------------------------------------

    def state_dict(self) -> Dict:
        """Return serialisable state for checkpointing."""
        return {
            'G_ema': self._G_ema,
            'damping': self._damping,
        }

    def load_state_dict(self, state: Dict) -> None:
        """Restore state from a checkpoint."""
        self._G_ema = state.get('G_ema')
        if 'damping' in state:
            self._damping = float(state['damping'])

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _flatten(grads: Tuple[Optional[torch.Tensor], ...]) -> torch.Tensor:
        """Concatenate a tuple of gradient tensors into a flat vector."""
        parts = []
        for g in grads:
            if g is None:
                continue
            parts.append(g.reshape(-1))
        if not parts:
            raise RuntimeError(
                "GaussNewtonOptimizer: all gradient entries are None. "
                "Ensure that the residual tensor depends on the model parameters."
            )
        return torch.cat(parts)

    @staticmethod
    def _apply_update(
        params: List[torch.nn.Parameter],
        delta: torch.Tensor,
    ) -> None:
        """Subtract delta (flat) from each parameter in-place."""
        offset = 0
        for p in params:
            sz = p.numel()
            p.data.sub_(delta[offset:offset + sz].reshape(p.shape))
            offset += sz
