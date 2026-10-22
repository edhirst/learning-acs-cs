"""Fibration-adapted representation of S^6 (and S^2 x S^4), stage 1: the
smooth torus part.

The manifold is presented as the compact complex threefold Y of
Alpöge-Claude / Engel (S^6) or Guo-Fang-Lu (S^2 x S^4): a fibration Y -> P^1
whose fibres over U = P^1 minus {0, 1, inf} are real 4-tori R^4/Lambda, with
monodromies T0 (order 3), T1 (order 4) and T_inf = (T0 T1)^-1 (unipotent),
taken from Engel, Thm 5.1 and eq. (6.1). Only this TOPOLOGICAL data enters;
the complex structure (the periods tau, mu, beta of the paper) is not used.

Over U the true structure is invariant under fibre translations (they act
holomorphically), so the search restricts to x-invariant fields: a 6x6 matrix
field J(t) on the base, in coordinates (t1, t2, x1..x4). The PDE becomes
two-dimensional.

Monodromy is built in exactly by a multivalued gauge
    G(t) = exp(arg(t)/2pi L0) exp(arg(t-1)/2pi L1),   L_i = log diag(I_2, T_i),
(branch cuts on (-inf, 0) and (1, inf)): J = G J~ G^-1 with J~ single-valued
jumps by conjugation with diag(I_2, T0) across the first cut and
diag(I_2, T1) across the second, which is exactly the gluing of the torus
bundle. J~ = e^B J0 e^-B (exact ACS), B a polynomial in t plus Laurent terms
in 1/t, 1/(t-1). Loss: the coordinate Nijenhuis tensor (fibre derivatives
vanish). Solver: dense Levenberg-Marquardt. Certificate: sup |N| on a fine
grid of the region and by adversarial ascent in the base.
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import math

import numpy as np
import scipy.linalg as sla
import torch

from src.loss.sphere_loss import _nij_tensor

DTYPE = torch.float64

T0 = np.array([[0, -1, -1, 0], [1, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], dtype=float)
T1 = np.array([[-1, 1, -1, 2], [-1, 0, -1, 1], [1, 1, 2, -1], [0, 0, 0, 1]], dtype=float)


def monodromy_logs(lattice_change=None, sign=1, dual=False):
    """Real logs of diag(I_2, T_i). `lattice_change` (4x4) re-expresses the
    monodromy in another lattice basis (e.g. the index-2 superlattice used for
    S^2 x S^4): T -> C^-1 T C. `sign` = -1 inverts both monodromies (loop
    orientation convention), `dual` uses the contragredient T^-T (transport of
    fibre coordinates vs of homology classes)."""
    out = []
    for T in (T0, T1):
        if lattice_change is not None:
            T = np.linalg.inv(lattice_change) @ T @ lattice_change
        if dual:
            T = np.linalg.inv(T).T
        if sign < 0:
            T = np.linalg.inv(T)
        L = sla.logm(T).real
        assert np.allclose(sla.expm(L), T, atol=1e-10)
        M = np.zeros((6, 6)); M[2:, 2:] = L
        out.append(torch.tensor(M, dtype=DTYPE))
    return out


def _expm(A, s=6, order=12):
    n = A.shape[-1]
    I = torch.eye(n, dtype=A.dtype).expand_as(A)
    X = A / 2.0 ** s
    E, T = I, I
    for k in range(1, order + 1):
        T = T @ X / k
        E = E + T
    for _ in range(s):
        E = E @ E
    return E


class FibreField(torch.nn.Module):
    """J(t) = G(t) e^{B(t)} J0 e^{-B(t)} G(t)^-1 on the base region."""

    def __init__(self, degree=6, laurent=3, scale=2.0, logs=None):
        super().__init__()
        self.L0, self.L1 = logs if logs is not None else monodromy_logs()
        self.degree, self.laurent, self.scale = degree, laurent, scale
        J0 = torch.zeros(6, 6, dtype=DTYPE)
        for k in range(3):
            J0[2 * k, 2 * k + 1], J0[2 * k + 1, 2 * k] = -1.0, 1.0
        self.register_buffer('J0', J0)
        self.n_feat = self.features(torch.zeros(1, 2, dtype=DTYPE) + 0.3).shape[1]
        self.W = torch.nn.Parameter(torch.zeros(self.n_feat, 6, 6, dtype=DTYPE))

    def features(self, t):
        u, v = (t[:, 0] - 0.5) / self.scale, t[:, 1] / self.scale
        f = [u ** a * v ** b for a in range(self.degree + 1) for b in range(self.degree + 1 - a)
             if a + b >= 1]
        z0 = torch.complex(t[:, 0], t[:, 1])
        for m in range(1, self.laurent + 1):
            for zz in (z0, z0 - 1.0):
                w = zz ** (-m)
                f += [w.real, w.imag]
        # log|t|, log|t-1| (the gauge-stripped periods carry exp(i log|t| L / 2pi)
        # factors), times low-degree polynomials
        lt, lt1 = torch.log(z0.abs()), torch.log((z0 - 1.0).abs())
        low = [torch.ones_like(u)] + [u ** a * v ** b for a in range(3) for b in range(3 - a) if a + b >= 1]
        for L in (lt, lt1, lt ** 2, lt1 ** 2, lt * lt1):
            f += [L * q for q in low]
        return torch.stack(f, dim=1)

    def gauge(self, t):
        a0 = torch.atan2(t[:, 1], t[:, 0]) / (2 * math.pi)                 # cut (-inf, 0)
        a1 = torch.remainder(torch.atan2(t[:, 1], t[:, 0] - 1.0), 2 * math.pi) / (2 * math.pi)  # cut (1, inf)
        return _expm(a0.view(-1, 1, 1) * self.L0) @ _expm(a1.view(-1, 1, 1) * self.L1)

    def J(self, t, W=None):
        W = self.W if W is None else W
        B = torch.einsum('bf,fij->bij', self.features(t), W)
        Jt = _expm(B) @ self.J0 @ _expm(-B)
        G = self.gauge(t)
        return G @ Jt @ torch.linalg.inv(G)


def fibre_nijenhuis(field, t, W=None):
    """Coordinate Nijenhuis tensor (B, 6, 6, 6) of the x-invariant field."""
    def single(ts):
        return field.J(ts.unsqueeze(0), W).squeeze(0)
    J = field.J(t, W)
    dJt = torch.func.vmap(torch.func.jacfwd(single))(t)                     # (B, 6, 6, 2)
    dJ = torch.cat([dJt, torch.zeros(*dJt.shape[:3], 4, dtype=DTYPE)], dim=3)
    return _nij_tensor(J, dJ)


def nij_rows(field, t, W=None):
    N = fibre_nijenhuis(field, t, W)
    iu = torch.triu_indices(6, 6, 1)
    return N[:, :, iu[0], iu[1]].reshape(len(t), -1) * math.sqrt(2.0)


def sample_region(n, r0=0.25, r1=0.25, R=3.0, seed=None):
    """Uniform points of {|t| <= R, |t| >= r0, |t-1| >= r1} (the base minus
    discs around the singular fibres; |t| > R is the cusp neighbourhood)."""
    g = torch.Generator().manual_seed(seed) if seed is not None else None
    out = []
    while sum(len(o) for o in out) < n:
        p = (torch.rand(4 * n, 2, generator=g, dtype=DTYPE) * 2 - 1) * R
        r = p.norm(dim=1)
        ok = (r <= R) & (r >= r0) & ((p - torch.tensor([1.0, 0.0], dtype=DTYPE)).norm(dim=1) >= r1)
        out.append(p[ok])
    return torch.cat(out)[:n]


# ---------------------------------------------------------------------------
# Period-map formulation (fibred structures): the unknown is the fibre
# (1,0)-coframe Pi(t) (2x4 complex); J is defined by (1,0) = span{dt, Pi dx}.
# ---------------------------------------------------------------------------

class PeriodField(torch.nn.Module):
    """Pi(t) = Pi~(t) G_F(t)^{-1}, Pi~ single-valued (real polynomial + Laurent
    features, complex coefficients). Exact monodromy: Pi jumps by right
    multiplication with T^{-1} across the cuts, i.e. x -> T x on the fibres."""

    def __init__(self, degree=6, laurent=3, scale=2.0, logs=None):
        super().__init__()
        self.base = FibreField(degree, laurent, scale, logs)      # reuse features + gauge
        nf = self.base.n_feat + 1                                  # + constant
        self.Wr = torch.nn.Parameter(torch.zeros(nf, 2, 4, dtype=DTYPE))
        self.Wi = torch.nn.Parameter(torch.zeros(nf, 2, 4, dtype=DTYPE))
        with torch.no_grad():                                      # constant generic start Pi = [I | iI]
            self.Wr[0, :, :2] = torch.eye(2, dtype=DTYPE)
            self.Wi[0, :, 2:] = torch.eye(2, dtype=DTYPE)
        self.n_params = 2 * self.Wr.numel()

    def feats(self, t):
        return torch.cat([torch.ones_like(t[:, :1]), self.base.features(t)], dim=1)

    def Pi(self, t, theta=None):
        if theta is None:
            Wr, Wi = self.Wr, self.Wi
        else:
            Wr, Wi = theta[: self.n_params // 2].view_as(self.Wr), theta[self.n_params // 2:].view_as(self.Wi)
        F = self.feats(t)
        Pt = torch.complex(torch.einsum('bf,fij->bij', F, Wr), torch.einsum('bf,fij->bij', F, Wi))
        GF = self.base.gauge(t)[:, 2:, 2:]
        return Pt @ torch.linalg.inv(GF).to(torch.complex128)

    def theta(self):
        return torch.cat([self.Wr.detach().reshape(-1), self.Wi.detach().reshape(-1)])


def period_residual(pf, t, theta=None, w_gauge=1.0, mu_bar=10.0, d0=0.3):
    """Rows: [Re, Im of dbar(Pi) (I - Pi^+ Pi)] (integrability), sqrt(w_gauge)
    (Pi Pi^H - I) (GL2(C) gauge) and sqrt(mu_bar) relu(d0 - |det[Pi; conj Pi]|^(1/4))
    (fibres stay complex tori)."""
    def single(ts):
        P = pf.Pi(ts.unsqueeze(0), theta).squeeze(0)
        return torch.view_as_real(P)
    P = pf.Pi(t, theta)
    dP = torch.func.vmap(torch.func.jacfwd(single))(t)            # (B, 2, 4, 2[re,im], 2[t1,t2])
    dPc = torch.complex(dP[..., 0, :], dP[..., 1, :])             # (B, 2, 4, 2)
    dbar = 0.5 * (dPc[..., 0] + 1j * dPc[..., 1])
    PH = P.conj().transpose(1, 2)
    proj = PH @ torch.linalg.inv(P @ PH) @ P                      # (B, 4, 4) onto row space
    I4 = torch.eye(4, dtype=torch.complex128)
    G2 = P @ PH
    Lc = torch.linalg.cholesky(G2)
    R = torch.linalg.solve_triangular(Lc, dbar @ (I4 - proj), upper=False)   # scale-invariant
    gauge = G2 - torch.eye(2, dtype=torch.complex128)
    D = torch.cat([P, P.conj()], dim=1)
    det = torch.linalg.det(D).abs() ** 0.25
    rows = [torch.view_as_real(R).reshape(len(t), -1),
            math.sqrt(w_gauge) * torch.view_as_real(gauge).reshape(len(t), -1),
            (math.sqrt(mu_bar) * torch.relu(d0 - det)).unsqueeze(1)]
    return torch.cat(rows, dim=1)


def J_from_period(pf, t, theta=None):
    """Real 6x6 J with (1,0)-forms dt and Pi dx: [Theta; conj Theta] J = [i Theta; -i conj Theta]."""
    P = pf.Pi(t, theta)
    B = len(t)
    Th = torch.zeros(B, 3, 6, dtype=torch.complex128)
    Th[:, 0, 0], Th[:, 0, 1] = 1.0, 1j                            # dt = dt1 + i dt2
    Th = torch.cat([Th[:, :1], torch.cat([torch.zeros(B, 2, 2, dtype=torch.complex128), P], dim=2)], dim=1)
    A = torch.cat([Th, Th.conj()], dim=1)
    Rhs = torch.cat([1j * Th, -1j * Th.conj()], dim=1)
    return (torch.linalg.inv(A) @ Rhs).real


class ConjPeriodField(PeriodField):
    """Pi(t) = Pi_ref e^{-B(t)} G_F(t)^{-1}, B real 4x4 (single-valued features):
    every fibre is a genuine complex torus (the fibre structure is a GL4(R)
    conjugate of a fixed one), so the degenerate locus det[Pi; conj Pi] = 0
    cannot be reached by a bounded B."""

    def __init__(self, degree=6, laurent=3, scale=2.0, logs=None):
        torch.nn.Module.__init__(self)
        self.base = FibreField(degree, laurent, scale, logs)
        nf = self.base.n_feat
        self.Wb = torch.nn.Parameter(torch.zeros(nf, 4, 4, dtype=DTYPE))
        P0 = torch.zeros(2, 4, dtype=torch.complex128)
        P0[0, 0], P0[0, 1], P0[1, 2], P0[1, 3] = 1, 1j, 1, 1j       # standard structure on R^4
        self.register_buffer('P0', P0)
        self.n_params = self.Wb.numel()

    def Pi(self, t, theta=None):
        W = self.Wb if theta is None else theta.view_as(self.Wb)
        B = torch.einsum('bf,fij->bij', self.base.features(t), W)
        A_inv = _expm(-B)
        GF = self.base.gauge(t)[:, 2:, 2:]
        return (self.P0 @ (A_inv @ torch.linalg.inv(GF)).to(torch.complex128))

    def theta(self):
        return self.Wb.detach().reshape(-1)
