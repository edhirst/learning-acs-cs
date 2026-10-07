"""The complex structure J_round = F_* J_Y on the ROUND S^4 x T^2 = S^4 x S^1 x S^1 in R^5 x R^2 x R^2, as a function of
the ambient point q (docs/s4t2_explicit_diffeo.md, section 12.7 / 13).

    J, DJ, res = J_round(q)        q: (B, 9) points of S^4 x T^2 (unit S^4, unit circles)
        J   (B, 9, 9)     tangential complex structure (J = P J P, J^2 = -P), the layout of src/loss/ambient_loss.py
        DJ  (B, 9, 9, 9)  DJ[b, a, r, c] = d_c J_ar along tangent directions (zero on normals), exact (nested AD)
        res (B,)          |F(F^-1(q)) - q|_inf of the inverse (about 1e-13)

Evaluation: p = F^-1(q) in a holomorphic chart of Y chosen automatically (tinverse.Finv), J = DF J0 DF^+ there, DJ
from D^2 F (nijenhuis_round.transport). J_round is invariant under the T^3 = S^1 x T^2 action
(Z -> e^{i theta} Z, circles rotated), because F is T^3-equivariant and T^3 acts holomorphically on Y.
"""
from __future__ import annotations

import torch

from extra.s4t2_explicit import nijenhuis_round as NR, tinverse as TI


def J_round(q, with_DJ=True):
    sol = TI.Finv(q)
    B = q.shape[0]
    J = torch.zeros(B, 9, 9)
    DJ = torch.zeros(B, 9, 9, 9) if with_DJ else None
    res = torch.zeros(B)
    for key, (idx, p, r) in sol.items():
        f = TI.CHART_F[key]
        if with_DJ:
            qq, A, H = NR.DF_and_D2F(f, p)
            Jk, DJk = NR.transport(A, H, qq)
            DJ[idx] = DJk
        else:
            A = torch.stack([torch.func.jvp(f, (p,), (NR.EYE6[k].expand_as(p),))[1] for k in range(6)], -1)
            Jk = A @ NR.J0 @ torch.linalg.pinv(A)
        J[idx] = Jk
        res[idx] = r
    return J, DJ, res


def nijenhuis(q, J, DJ):
    """Per-sample intrinsic round-metric |N| and the cancellation scale (sum of the bracket-term norms)."""
    Nn, scale, _ = NR.nij(q, J, DJ)
    return Nn, scale
