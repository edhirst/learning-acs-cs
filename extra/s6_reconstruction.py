"""Numerical reconstruction of the complex structure on S^6 (and S^2 x S^4)
from the learned period map (extra/period_map.py), and its certification.

Model of the smooth part [AC26, Rem. 3.20]: X over B° is (h_Z° x R^4/Z^4)/Delta
with the deck action (Z, x) -> (g_j Z, A_j x), A_j = (T_j^-1)^t, and the complex
structure J whose (1,0)-forms are dZ and Pi(Z) dx, Pi = [[6mu, tau, 1, 0],
[beta, mu, 0, 1]]. J is x-invariant; in the real coordinates
(Re Z, Im Z, x1..x4) it is a 6x6 field of Z only.

Certificate (all on fresh points):
  (a) |N(J)| on the fundamental domain (full coordinate Nijenhuis tensor by AD
      through the q-series) - zero up to round-off, the period map being an
      exact q-series;
  (b) deck consistency |J(gZ) - DPhi_g J(Z) DPhi_g^-1| on the gluing pairs;
  (c) |N| of the globally BLENDED field across the domain edges (partition of
      unity between J and the deck-transported J), the honest single number
      for the compact model;
  (d) D < 0 on the domain (every fibre is a lattice quotient, a complex torus).
The multiple-fibre fillings use the same Pi on discs around Z1, Z2 (twisted
translations leave the x-invariant J unchanged); the cusp filling is the toric
quotient of [AC26, §4] built from the same q-series C(q), holomorphic by
construction.
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import math

import numpy as np
import torch

from extra.period_map import TriangleGroup, mobius, A1, A2
from src.loss.sphere_loss import _nij_tensor

DTYPE = torch.float64
C128 = torch.complex128


class TorchPeriodMap:
    """tau, mu, beta as torch functions of real (a, b), Z = a + i b."""

    def __init__(self, path):
        d = np.load(path)
        self.a = torch.tensor(d['a'], dtype=C128)
        self.b = torch.tensor(d['b'], dtype=C128)
        self.c = torch.tensor(d['cb'], dtype=C128)
        self.q0 = float(d['q0'])
        self.n = torch.arange(len(self.a), dtype=DTYPE)

    def _series(self, Z, coef):
        q = torch.exp(2j * math.pi * Z) / self.q0
        return (q.unsqueeze(-1) ** self.n.to(C128)) @ coef

    def tmb(self, Z):
        return Z + self._series(Z, self.a), self._series(Z, self.b), -Z + self._series(Z, self.c)

    def Pi(self, Z):
        tau, mu, beta = self.tmb(Z)
        o, n = torch.ones_like(tau), torch.zeros_like(tau)
        r1 = torch.stack([6 * mu, tau, o, n], dim=-1)
        r2 = torch.stack([beta, mu, n, o], dim=-1)
        return torch.stack([r1, r2], dim=-2)                     # (B, 2, 4)


def J_of_Pi(Pi, lattice=None):
    """Real 6x6 J (coordinates Re Z, Im Z, x) with (1,0)-forms dZ and Pi dx.
    `lattice` (4x4): fibre coordinates x' with x = lattice @ x' (e.g. a superlattice basis)."""
    if lattice is not None:
        Pi = Pi @ lattice.to(C128)
    B = Pi.shape[0]
    top = torch.zeros(B, 1, 6, dtype=C128)
    top = top + torch.tensor([1.0, 1j, 0, 0, 0, 0], dtype=C128)
    bot = torch.cat([torch.zeros(B, 2, 2, dtype=C128), Pi], dim=2)
    Th = torch.cat([top, bot], dim=1)                            # (B, 3, 6)
    A = torch.cat([Th, Th.conj()], dim=1)
    R = torch.cat([1j * Th, -1j * Th.conj()], dim=1)
    return (torch.linalg.inv(A) @ R).real


def J_field(pm, ab, lattice=None):
    Z = torch.complex(ab[:, 0], ab[:, 1])
    return J_of_Pi(pm.Pi(Z), lattice)


def nijenhuis(J_fn, ab):
    """Full coordinate Nijenhuis tensor of an x-invariant field J(a, b)."""
    def single(p):
        return J_fn(p.unsqueeze(0)).squeeze(0)
    J = J_fn(ab)
    dJ2 = torch.func.vmap(torch.func.jacfwd(single))(ab)          # (B, 6, 6, 2)
    dJ = torch.cat([dJ2, torch.zeros(*dJ2.shape[:3], 4, dtype=DTYPE)], dim=3)
    return _nij_tensor(J, dJ)


def mobius_t(M, Z):
    M = torch.tensor(M, dtype=C128)
    return (M[0, 0] * Z + M[0, 1]) / (M[1, 0] * Z + M[1, 1])


def deck_differential(M, Z, A):
    """DPhi of (Z, x) -> (gZ, A x) in real coordinates (a, b, x): the base block
    is the real 2x2 form of the complex derivative g'(Z)."""
    M = torch.tensor(M, dtype=C128)
    gp = 1.0 / (M[1, 0] * Z + M[1, 1]) ** 2
    B = Z.shape[0]
    D = torch.zeros(B, 6, 6, dtype=DTYPE)
    blk = torch.stack([torch.stack([gp.real, -gp.imag], -1), torch.stack([gp.imag, gp.real], -1)], -2)
    D = torch.cat([torch.cat([blk, torch.zeros(B, 2, 4, dtype=DTYPE)], 2),
                   torch.cat([torch.zeros(B, 4, 2, dtype=DTYPE), torch.tensor(A, dtype=DTYPE).expand(B, 4, 4)], 2)], 1)
    return D


# ---------------------------------------------------------------------------
# Global (blended) field across a domain edge, adversarial sup, full report
# ---------------------------------------------------------------------------

def polar6(K, iters=40):
    """J = K (-K^2)^{-1/2} on R^6 (coupled Newton-Schulz)."""
    I = torch.eye(6, dtype=K.dtype).expand_as(K)
    Y, Z = (-K @ K) / 2.0, I
    for _ in range(iters):
        T = 0.5 * (3.0 * I - Z @ Y)
        Y, Z = Y @ T, T @ Z
    return K @ Z / math.sqrt(2.0)


def blended_field(pm, g, A, centre, kappa=30.0, lattice=None, Alat=None):
    """Smooth global description of J across the edge identified by g:
    J_blend = polar(w J + (1 - w) DPhi^-1 J(g .) DPhi), w = (1 + tanh(kappa (a - a0)))/2.
    On the compact model this is what any partition of unity produces; it is
    integrable iff the two descriptions agree."""
    a0 = float(np.real(centre))
    Am = A if Alat is None else Alat

    def J_fn(ab):
        Z = torch.complex(ab[:, 0], ab[:, 1])
        JA = J_of_Pi(pm.Pi(Z), lattice)
        Zg = mobius_t(g, Z)
        D = deck_differential(g, Z, Am)
        JB = torch.linalg.inv(D) @ J_of_Pi(pm.Pi(Zg), lattice) @ D
        w = 0.5 * (1 + torch.tanh(kappa * (ab[:, 0] - a0))).view(-1, 1, 1)
        return polar6(w * JA + (1 - w) * JB)
    return J_fn


def adversarial_sup_2d(r_fn, ab0, steps=60, lr=0.01, valid=None):
    """Gradient ascent of a per-point residual in the base coordinates,
    constrained to the set where the model is used (`valid(ab) -> bool mask`)."""
    x = ab0.clone()
    rc = r_fn(x).detach()
    step = torch.full((len(x), 1), lr, dtype=DTYPE)
    for _ in range(steps):
        xg = x.clone().requires_grad_(True)
        r = r_fn(xg)
        g, = torch.autograd.grad(r.sum(), xg)
        xn = x + step * g / (g.norm(dim=1, keepdim=True) + 1e-30)
        rn = r_fn(xn).detach()
        up = rn > rc
        if valid is not None:
            up = up & valid(xn)
        x = torch.where(up.unsqueeze(1), xn.detach(), x)
        rc = torch.where(up, rn, rc)
        step = torch.where(up.unsqueeze(1), step * 1.2, step * 0.5)
    return rc.max().item()
