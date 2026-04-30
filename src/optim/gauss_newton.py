"""
Full-matrix Gauss-Newton optimizer for physics-informed neural networks.

Implements the rank-1 EMA estimator of the Gauss-Newton (GN) matrix described
in Wang et al. (2025) "Discovery of Unstable Singularities" (arXiv:2509.14185),
§ "Second-order Optimizer", together with the Martens-Grosse optimal learning
rate + momentum scheme from Martens & Grosse (ICML 2015).

Algorithm (one step):
  1. Compute residuals r (flat, shape m) via residual_fn().
  2. Compute gradient g = Jᵀ r via VJP with adjoint r.
  3. For k = 1..K: sample z_k ~ N(0, Iₘ); compute u_k = Jᵀ z_k via VJP.
  4. EMA update:  Ĝ ← β Ĝ + (1−β) mean_k outer(u_k, u_k).
  5. Solve (Ĝ + γ I) δ = g.
  6. Martens-Grosse optimal (α*, μ*): minimise the quadratic model
       q(v) = gᵀv + ½ vᵀ(Ĝ+γI)v  over  v = −α δ + μ v_{t−1}
     by solving the 2×2 linear system (reduces to scalar α when no velocity).
  7. Update: v_t ← −α* δ + μ* v_{t−1};  θ ← θ + v_t.
  8. [if adaptive_damping] Compute ρ = actual / predicted decrease; adjust γ.

This is feasible only for small networks (≲ tens of thousands of parameters)
because the GN matrix is stored and inverted as a dense n × n matrix.
"""

from __future__ import annotations

import math
from typing import Callable, Dict, Iterable, List, Optional, Tuple

import torch


class GaussNewtonOptimizer:
    """Full-matrix Gauss-Newton optimizer with rank-K EMA and Martens-Grosse momentum.

    Args:
        params: Iterable of model parameters to optimise.
        damping: Initial Levenberg-Marquardt damping coefficient γ.
        ema_decay: Exponential moving-average decay β for the GN matrix.
            Should be ≈ 1 − 1/n for an n-parameter network so the EMA window
            spans enough steps to accumulate a full-rank estimate.
        adaptive_damping: If True, apply Levenberg-Marquardt damping
            adjustment based on the ratio of actual to predicted decrease.
        damping_up_factor: Multiply γ by this factor when ρ < rho_threshold_bad.
        damping_down_factor: Divide γ by this factor when ρ > rho_threshold_good.
        rho_threshold_bad: ρ below this triggers damping increase; ρ < 0 also
            causes step rejection.
        rho_threshold_good: ρ above this triggers damping decrease.
        num_gn_samples: Number of random probe vectors K used per step to
            form a rank-K EMA update of the GN matrix.  K=1 matches the paper;
            higher K warms up the estimate faster at the cost of K extra VJPs.
        use_momentum: If True (default), use the Martens-Grosse optimal (α, μ)
            scheme: jointly minimise the quadratic model over the combined
            direction v = −α δ + μ v_{t−1}.  If False, use the scalar auto-LR
            α* = (g·δ)/(δᵀ(Ĝ+γI)δ) with no momentum (μ = 0).
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
        num_gn_samples: int = 1,
        use_momentum: bool = True,
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
        self._num_gn_samples = max(1, int(num_gn_samples))
        self._use_momentum = bool(use_momentum)
        self._G_ema: Optional[torch.Tensor] = None   # dense (n, n) GN matrix EMA
        self._velocity: Optional[torch.Tensor] = None  # v_{t-1}: last Δθ

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

        # ---- Step 3: rank-K GN estimate  u_k = Jᵀ z_k, z_k ~ N(0, I) --
        # Accumulate K rank-1 outer products before the EMA update so the
        # effective per-step rank is K rather than 1.
        u_outer_sum: Optional[torch.Tensor] = None
        for k in range(self._num_gn_samples):
            retain = (k < self._num_gn_samples - 1)  # keep graph for all but last
            z = torch.randn(r.shape[0], dtype=r.dtype, device=r.device)
            grads_z = torch.autograd.grad(
                r, params,
                grad_outputs=z,
                retain_graph=retain,
                create_graph=False,
            )
            u = self._flatten(grads_z)      # (n,)

            # Guard: skip update if any sample is non-finite.
            if not torch.isfinite(u).all():
                return info

            outer_k = torch.outer(u, u)
            if u_outer_sum is None:
                u_outer_sum = outer_k
            else:
                u_outer_sum = u_outer_sum + outer_k

        # Average over K samples so scale is independent of num_gn_samples.
        u_outer_avg = u_outer_sum / self._num_gn_samples

        # ---- Step 4: EMA update of GN matrix --------------------------
        if self._G_ema is None:
            self._G_ema = u_outer_avg
        else:
            self._G_ema.mul_(self._ema_decay).add_((1.0 - self._ema_decay) * u_outer_avg)

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

        # ---- Step 6: Martens-Grosse optimal (α, μ) --------------------
        # We minimise the quadratic model  q(v) = gᵀv + ½vᵀ(Ĝ+γI)v
        # over  v = −α δ + μ v_{t−1}  (two free scalars).
        #
        # Defining F = Ĝ+γI, p = −δ (GN descent direction), q = v_{t−1}:
        #   ∂q/∂α = 0  and  ∂q/∂μ = 0  give the 2×2 system:
        #     [pᵀFp,  pᵀFq] [α]   [−gᵀp]   [g·δ ]
        #     [pᵀFq,  qᵀFq] [μ] = [−gᵀq] = [−g·v]
        # (with p=−δ: pᵀFp = δᵀFδ, pᵀFq = −δᵀFv, gᵀp = −g·δ)
        #
        # When v_{t−1} is absent (first step) this collapses to the scalar
        # α* = (g·δ)/(δᵀFδ) with μ = 0.
        Gd = self._G_ema @ delta          # (Ĝ)δ, needed for quadratic terms
        dFd = (delta * Gd).sum() + self._damping * (delta * delta).sum()  # δᵀFδ
        gd  = (g * delta).sum()           # g·δ  (≥ 0 when δ is a descent dir)

        v_prev = self._velocity
        if self._use_momentum and v_prev is not None:
            Gv  = self._G_ema @ v_prev
            dFv = (delta * Gv).sum() + self._damping * (delta * v_prev).sum()  # δᵀFv
            vFv = (v_prev * Gv).sum() + self._damping * (v_prev * v_prev).sum()  # vᵀFv
            gv  = (g * v_prev).sum()       # g·v

            # 2×2 system (see docstring derivation)
            det = dFd * vFv - dFv * dFv
            # Fall back to no-momentum if the system is ill-conditioned.
            if det.abs() > 1e-14 * (dFd.abs() * vFv.abs()).clamp(min=1e-30):
                alpha = ((gd * vFv - gv * dFv) / det).clamp(min=0.0, max=10.0)
                mu    = ((gd * dFv - gv * dFd) / (-det)).clamp(min=-0.99, max=0.99)
                # Predicted decrease for adaptive damping (full momentum model)
                _pred_denom = alpha
            else:
                alpha = (gd / dFd.clamp(min=1e-30)).clamp(min=0.0, max=10.0)
                mu    = delta.new_zeros(())
        else:
            alpha = (gd / dFd.clamp(min=1e-30)).clamp(min=0.0, max=10.0)
            mu    = delta.new_zeros(())
            v_prev = None  # simplify update below

        # ---- Adaptive damping: snapshot before update -----------------
        if self._adaptive_damping:
            loss_before = float((r.detach() ** 2).sum())
            old_data = [p.data.clone() for p in params]

        # ---- Step 7: apply update  θ ← θ + v_t  where v_t = −α δ + μ v_{t−1}
        if v_prev is not None:
            new_velocity = -alpha * delta + mu * v_prev
        else:
            new_velocity = -alpha * delta

        update_norm = new_velocity.norm()
        if update_norm > self._max_step_norm:
            new_velocity = new_velocity * (self._max_step_norm / update_norm)

        # θ ← θ + v_t  ↔  subtract (−v_t) from θ
        self._apply_update(params, -new_velocity)
        self._velocity = new_velocity.detach()

        # ---- Step 8: Levenberg-Marquardt damping update ---------------
        if self._adaptive_damping:
            r_new, _ = residual_fn()
            loss_after = float((r_new.detach() ** 2).sum())

            actual_decrease = loss_before - loss_after
            # Predicted decrease under the quadratic model h ≈ hᵀh:
            #   Δh ≈ 2 gᵀ v_t + vᵀ (Ĝ+γI) v_t
            # For the pure-α case this is 2α·gd − α²·dFd.
            # We use the full velocity vector for generality.
            v_t = new_velocity
            Gv_t = self._G_ema @ v_t
            vFv_t = (v_t * Gv_t).sum() + self._damping * (v_t * v_t).sum()
            predicted_decrease = float(2.0 * (g * v_t).sum() + vFv_t)
            if abs(predicted_decrease) < 1e-30:
                rho = 1.0
            else:
                rho = actual_decrease / predicted_decrease

            if rho < self._rho_bad:
                self._damping = min(self._damping * self._damping_up_factor, 1e10)
                if rho < 0.0:
                    # Reject step: restore parameters and clear velocity
                    for p, old in zip(params, old_data):
                        p.data.copy_(old)
                    self._velocity = None
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
            'velocity': self._velocity,
        }

    def load_state_dict(self, state: Dict) -> None:
        """Restore state from a checkpoint."""
        self._G_ema = state.get('G_ema')
        if 'damping' in state:
            self._damping = float(state['damping'])
        self._velocity = state.get('velocity')

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
