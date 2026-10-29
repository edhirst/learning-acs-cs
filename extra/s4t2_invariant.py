"""Unbiased search for a T^2-INVARIANT complex structure on the round S^4 x T^2.

Symmetry ansatz: J is assumed invariant under the translations of the T^2
factor, an ansatz motivated by the supplied period model
(extra/s4t2_reconstruction.py). Whether the round S^4 x T^2 carries an
integrable J with this symmetry is not established. Only the symmetry is used;
the fibration, the periods and the fillings of that model are not.

Reduced ambient model. A T^2-invariant field on S^4 x T^2 is a field on S^4 with
values in End(TS^4 + R^2), R^2 spanned by the unit Killing fields d/dtheta_1,2
(which commute with every invariant field). In R^7 = R^5 + R^2 with
p = (x, 0, 0), P(x) = diag(I_5 - x x^T/|x|^2, I_2), and fields depending on x
only, the ambient bracket [U, V] = DV U - DU V is the bracket of S^4 x T^2
(the derivatives along the flat directions vanish), so the ambient Nijenhuis
code (src/loss/ambient_loss.py) computes N_J exactly, in the intrinsic norm of
the product of unit spheres. The problem becomes a PDE on S^4 (4-dimensional
domain) for a 6 x 6 field: polynomial classes of degree k in x in R^5 have
49 * binom(5 + k, k) coefficients (degree 4: 6174).
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import math

import torch
import torch.nn as nn

from src.models.poly_field import (PolyField, PolyCache, monomial_keys, monomial_index, fit_field,
                                   poly_residual, poly_residual_and_jacobian, poly_split_losses,
                                   _factorised_jacobian)
from src.loss.ambient_loss import ambient_nij_norm, compute_J_and_DJ, compute_P_and_DP

DTYPE = torch.float64


class S4T2Reduced:
    """S^4 in R^5 plus two flat Killing directions: ambient R^7, tangent rank 6."""
    ambient_dim, dim, n_x = 7, 6, 5

    def sample(self, n, dtype=DTYPE):
        v = torch.randn(n, 5, dtype=dtype)
        return torch.cat([v / v.norm(dim=1, keepdim=True), torch.zeros(n, 2, dtype=dtype)], 1)

    def projector(self, p):
        x = p[:, :5]
        B = p.shape[0]
        r2 = (x * x).sum(1, keepdim=True).unsqueeze(2)
        top = torch.eye(5, dtype=p.dtype).expand(B, 5, 5) - x.unsqueeze(2) * x.unsqueeze(1) / r2
        z52 = torch.zeros(B, 5, 2, dtype=p.dtype)
        bot = torch.cat([torch.zeros(B, 2, 5, dtype=p.dtype), torch.eye(2, dtype=p.dtype).expand(B, 2, 2)], 2)
        return torch.cat([torch.cat([top, z52], 2), bot], 1)

    def retract(self, p):
        x = p[:, :5]
        return torch.cat([x / x.norm(dim=1, keepdim=True), torch.zeros_like(p[:, 5:])], 1)

    def tangential(self, p, g):
        x, gx = p[:, :5], g[:, :5]
        return torch.cat([gx - (gx * x).sum(1, keepdim=True) * x, torch.zeros_like(g[:, 5:])], 1)


class XPolyField(PolyField):
    """K = P M P with M a degree-<=k matrix polynomial in x (the S^4 coordinates only)."""

    def __init__(self, degree, init_scale=0.3):
        nn.Module.__init__(self)
        self.n_in, self.degree, self.trilin = 7, degree, False
        self.keys = monomial_keys(5, degree)
        self.register_buffer('idx', monomial_index(self.keys, 7))
        self.n_feat = self.idx.shape[0]
        self.W = nn.Parameter(torch.randn(49, self.n_feat, dtype=DTYPE) * init_scale / math.sqrt(self.n_feat))
        self.n_mat, self.stages = 7, 1


MAN = S4T2Reduced()


def K_of(field):
    return lambda p: MAN.projector(p) @ field(p)[:, 0] @ MAN.projector(p)


def acs_residual(K_fn, p):
    K = K_fn(p)
    return torch.sqrt(((K @ K + MAN.projector(p)) ** 2).sum(dim=(1, 2)) + 1e-20)


def nij_norm(J_fn, p):
    J, DJ = compute_J_and_DJ(J_fn, p)
    P, DP = compute_P_and_DP(MAN, p)
    return ambient_nij_norm(J, DJ, P, DP)


def inv_sqrt_ns(S, iters=30):
    I = torch.eye(S.shape[-1], dtype=S.dtype).expand_as(S)
    Y, Z = S / 2.0, I
    for _ in range(iters):
        T = 0.5 * (3.0 * I - Z @ Y)
        Y, Z = Y @ T, T @ Z
    return Z / math.sqrt(2.0)


def polar(K_fn):
    def J(p):
        K = K_fn(p)
        S = -K @ K + (torch.eye(7, dtype=K.dtype) - MAN.projector(p))
        return K @ inv_sqrt_ns(S)
    return J


def adversarial_sup(r_fn, n_starts=128, steps=80, lr=0.03):
    xs = MAN.sample(32 * n_starts)
    with torch.no_grad():
        rs = torch.cat([r_fn(xs[i:i + 512]).detach() for i in range(0, len(xs), 512)])
    x = xs[torch.topk(rs, n_starts).indices]
    cur = r_fn(x).detach()
    step = torch.full((n_starts, 1), lr, dtype=DTYPE)
    for _ in range(steps):
        xg = x.clone().requires_grad_(True)
        g, = torch.autograd.grad(r_fn(xg).sum(), xg)
        g = MAN.tangential(xg, g)
        xn = MAN.retract(x + step * g / (g.norm(dim=1, keepdim=True) + 1e-12))
        rn = r_fn(xn).detach()
        up = rn > cur
        x = torch.where(up.unsqueeze(1), xn.detach(), x)
        cur = torch.where(up, rn, cur)
        step = torch.where(up.unsqueeze(1), step * 1.2, step * 0.5).clamp(1e-6, 0.3)
    return cur.max().item()


# ---------------------------------------------------------------------------
# orientation index (Section 4.5 of the report) for the reduced model
# ---------------------------------------------------------------------------

def _pfaffian6(A):
    from extra.campaign_unit import pfaffian
    return pfaffian(A)


def orientation_rows(K, E, s, mu):
    om = torch.einsum('bpa,bpq,bqc->bac', E, K.transpose(1, 2) - K, E)
    o = _pfaffian6(om) / 8.0
    return math.sqrt(mu) * torch.relu(1.0 - s * o).unsqueeze(1), o


def oriented_frame(cache_E, p):
    """Fix the sign of the frame so that (x-normal, E) is a positive basis of R^7."""
    nrm = torch.cat([p[:, :5], torch.zeros_like(p[:, 5:])], 1).unsqueeze(2)
    sgn = torch.sign(torch.linalg.det(torch.cat([nrm, cache_E], 2)))
    return torch.cat([cache_E[:, :, :1] * sgn.view(-1, 1, 1), cache_E[:, :, 1:]], 2)
