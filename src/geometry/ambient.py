"""
Ambient-projector representation of tangential structures on embedded sphere
products — the chart-free formulation for NON-parallelizable manifolds
(S⁴×S² ⊂ ℝ⁸, S⁶ ⊂ ℝ⁷), where no global frame exists (χ ≠ 0).

A tangential endomorphism field is J(x) = P(x) J̃(x) P(x) with P(x) the
orthogonal projector onto T_xM (explicit for sphere factors) and J̃ a single
network on ambient coordinates. One global function ⇒ cross-chart consistency
is automatic: there are no charts and no overlap loss. This removes the same
two structural obstacles (overlap-vs-Nijenhuis stalemate, chart-derivative
inflation) that the global-frame formulation removed on the parallelizable
cases.

All maps are smooth differentiable torch on a neighborhood of M.
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import torch


class AmbientManifold:
    """Product of round spheres  S^{m₁} × … × S^{m_k} ⊂ ℝ^{Σ(mᵢ+1)}."""

    def __init__(self, factor_dims):
        """factor_dims: manifold dims of the factors, e.g. (4, 2) or (6,)."""
        self.factor_dims = tuple(factor_dims)
        if not self.factor_dims:
            raise ValueError("a sphere product must contain at least one factor")
        if any(not isinstance(m, int) or isinstance(m, bool) or m < 1
               for m in self.factor_dims):
            raise ValueError(f"sphere dimensions must be positive integers, got {self.factor_dims!r}")
        self.ambient_dims = tuple(m + 1 for m in self.factor_dims)
        self.ambient_dim = sum(self.ambient_dims)
        self.dim = sum(self.factor_dims)

    def sample(self, n: int, dtype=torch.float64) -> torch.Tensor:
        parts = []
        for na in self.ambient_dims:
            v = torch.randn(n, na, dtype=dtype)
            parts.append(v / v.norm(dim=1, keepdim=True))
        return torch.cat(parts, dim=1)

    def projector(self, x: torch.Tensor) -> torch.Tensor:
        """P(x): block-diagonal of per-factor I − x xᵀ/|x|², shape (B, N, N).

        The |x|² division makes P a smooth ambient extension (defined off the
        sphere), so its ambient derivative DP is meaningful and autograd-able.
        """
        B = x.shape[0]
        N = self.ambient_dim
        # Out-of-place block-diagonal assembly (vmap/jacfwd-safe: no in-place
        # writes into a preallocated tensor).
        row_blocks = []
        off = 0
        for na in self.ambient_dims:
            xf = x[:, off:off + na]
            r2 = (xf * xf).sum(dim=1, keepdim=True).unsqueeze(2)
            eye = torch.eye(na, dtype=x.dtype, device=x.device).unsqueeze(0).expand(B, na, na)
            blk = eye - xf.unsqueeze(2) * xf.unsqueeze(1) / r2
            left = torch.zeros(B, na, off, dtype=x.dtype, device=x.device)
            right = torch.zeros(B, na, N - off - na, dtype=x.dtype, device=x.device)
            row_blocks.append(torch.cat([left, blk, right], dim=2))
            off += na
        return torch.cat(row_blocks, dim=1)


S4xS2 = AmbientManifold((4, 2))
S6 = AmbientManifold((6,))
