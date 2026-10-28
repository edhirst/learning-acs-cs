"""Spectral PINN for the period map of the (3,4,inf) family of 2-tori
[Alpöge-Claude 2026, §2-3]: holomorphic tau, mu, beta on the upper half plane
h_Z, determined (Theorem 3.4 there) by the functional equations along the
triangle group Delta = <g1, g2 | g1^3 = g2^4 = 1>,

    (tau2)  tau o g1 = (tau - 1)/tau,       tau o g2 = -1/tau          (tau o g0 = tau - 1)
    (mu1)   mu  o g1 = (1 - mu)/tau,        mu  o g2 = 1 + mu/tau      (mu  o g0 = mu)
    (beta1) beta o g1 = beta + 2 - 6(1-mu)^2/tau,
            beta o g2 = beta - 3 - 6 mu^2/tau                           (beta o g0 = beta + 1)

with mu, beta + tau bounded at the distinguished cusp, and D = Im beta -
6 (Im mu)^2 / Im tau < 0 (the fibres are lattices). RECONSTRUCTION, not
discovery: the equations are taken from the paper; no values of tau, mu,
beta are given to the solver.

Representation (exactly holomorphic, hence N = 0 by construction for the
complex structure built from them): in the coordinate Z with the cusp at inf
and g0(Z) = Z - 1, the functions tau - Z, mu, beta + Z are g0-invariant and
holomorphic on h_Z, i.e. power series in q = exp(2 pi i Z) convergent on
|q| < 1. Training: residuals of the g1, g2 equations on discs around their
fixed points Z1, Z2 (analytic continuation does the rest). tau: Levenberg-
Marquardt (nonlinear); mu: linear least squares given tau; beta: linear least
squares given tau, mu, up to its additive constant c0 (the family parameter).
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import math

import numpy as np

CPX = np.complex128


# ---------------------------------------------------------------------------
# The triangle group in the cusp coordinate
# ---------------------------------------------------------------------------

def _rot(p, theta):
    """Counterclockwise rotation about p in h by theta (det-1 matrix)."""
    C = np.array([[1, -p], [1, -np.conj(p)]], dtype=CPX)
    R = np.diag([np.exp(1j * theta / 2), np.exp(-1j * theta / 2)])
    M = np.linalg.inv(C) @ R @ C
    return M / np.sqrt(np.linalg.det(M))


def mobius(M, Z):
    return (M[0, 0] * Z + M[0, 1]) / (M[1, 0] * Z + M[1, 1])


class TriangleGroup:
    """Delta(3,4,inf) as in [AC26, Prop. 2.11]: triangle A=z1 (angle pi/3),
    B=z2 (angle pi/4), C=inf counterclockwise; x, y the ccw rotations about A,
    B by 2pi/3, 2pi/4; g1 = x^-1, g2 = y^-1, g0 = (g1 g2)^-1 parabolic. Then
    conjugated by z -> Z = -1/(c z) so that the fixed point of g0 is Z = inf and
    g0(Z) = Z - 1 (the sign eps0 = -1 of [AC26, Rem. 2.12])."""

    def __init__(self):
        z1, z2 = np.exp(2j * np.pi / 3), np.exp(1j * np.pi / 4)
        g1 = np.linalg.inv(_rot(z1, 2 * np.pi / 3))
        g2 = np.linalg.inv(_rot(z2, 2 * np.pi / 4))
        g0 = np.linalg.inv(g1 @ g2)
        g0 = g0 / g0[1, 1]
        c = g0[1, 0].real                                  # g0 = [[1,0],[c,1]] fixes z = 0
        S = np.array([[0, -1], [c, 0]], dtype=CPX)          # z -> Z = -1/(c z)
        Si = np.linalg.inv(S)
        conj = lambda M: S @ M @ Si
        self.g1, self.g2, self.g0 = conj(g1), conj(g2), conj(g0)
        self.Z1, self.Z2 = mobius(S, z1), mobius(S, z2)
        self.c = c

    def check(self):
        mp = np.linalg.matrix_power
        out = {
            'g1 fixes Z1': abs(mobius(self.g1, self.Z1) - self.Z1),
            'g2 fixes Z2': abs(mobius(self.g2, self.Z2) - self.Z2),
            'g1^3 = +-I': min(np.abs(mp(self.g1, 3) - s * np.eye(2)).max() for s in (1, -1)),
            'g2^4 = +-I': min(np.abs(mp(self.g2, 4) - s * np.eye(2)).max() for s in (1, -1)),
            'g0(Z) - (Z - 1)': abs(mobius(self.g0, 0.3 + 0.7j) - (0.3 + 0.7j - 1)),
        }
        return out


# Monodromy data of [AC26, Setup] in the basis (gamma, u, w, delta) of V (columns = images)
T1 = np.array([[1, 0, -6, 2], [0, -1, 1, 1], [0, -1, 0, 1], [0, 0, 0, 1]], dtype=float)
T2 = np.array([[1, 6, 0, -3], [0, 0, -1, 1], [0, 1, 0, 0], [0, 0, 0, 1]], dtype=float)
A1, A2 = np.linalg.inv(T1).T, np.linalg.inv(T2).T             # action on Lambda = V*


# ---------------------------------------------------------------------------
# q-series period map
# ---------------------------------------------------------------------------

class PeriodMap:
    """tau = Z + sum a_n (q/q0)^n, mu = sum b_n (q/q0)^n, beta = -Z + sum c_n (q/q0)^n."""

    def __init__(self, K=40, q0=0.2):
        self.K, self.q0 = K, q0
        self.a = np.zeros(K + 1, dtype=CPX)
        self.b = np.zeros(K + 1, dtype=CPX)
        self.cb = np.zeros(K + 1, dtype=CPX)

    def basis(self, Z):
        q = np.exp(2j * np.pi * Z) / self.q0
        return q[:, None] ** np.arange(self.K + 1)[None, :]

    def tau(self, Z, a=None):
        return Z + self.basis(Z) @ (self.a if a is None else a)

    def mu(self, Z, b=None):
        return self.basis(Z) @ (self.b if b is None else b)

    def beta(self, Z, c=None):
        return -Z + self.basis(Z) @ (self.cb if c is None else c)


def collocation_full(G, n_disc=64, n_edge=1500, seed=11):
    """Collocation for the g1- and g2-equations: discs around the fixed points
    plus the Ford-domain gluing pairs (Z with g Z, and W = g^-1 Z with Z, near
    the domain), so every coefficient is constrained where the gluing happens."""
    P1, P2 = collocation(G, n=n_disc, radius=0.25, rings=3)
    pairs = edge_pairs(G, n_edge, np.random.default_rng(seed))
    E1 = np.concatenate([pairs['g1'], mobius(np.linalg.inv(G.g1), pairs['g1inv'])])
    E2 = np.concatenate([pairs['g2'], mobius(np.linalg.inv(G.g2), pairs['g2inv'])])
    return np.concatenate([P1, E1]), np.concatenate([P2, E2])


def collocation(G, n=96, radius=0.06, rings=3):
    """Points on concentric circles around Z1 and Z2 (hyperbolically small)."""
    pts = []
    for Zc in (G.Z1, G.Z2):
        for k in range(1, rings + 1):
            r = radius * k / rings * Zc.imag
            th = 2 * np.pi * (np.arange(n) + 0.5 * k) / n
            pts.append(Zc + r * np.exp(1j * th))
    P = np.concatenate(pts)
    m = len(P) // 2
    return P[:m], P[m:]                                          # around Z1, around Z2


def tau_residual(pm, G, P1, P2, a):
    t1, t1g = pm.tau(P1, a), pm.tau(mobius(G.g1, P1), a)
    t2, t2g = pm.tau(P2, a), pm.tau(mobius(G.g2, P2), a)
    r = np.concatenate([t1g - (t1 - 1) / t1, t2g + 1 / t2])
    return r


def tau_jacobian(pm, G, P1, P2, a):
    """Complex Jacobian d r / d a of `tau_residual` (r is holomorphic in a):
    r1 = tau(g1 Z) - 1 + 1/tau(Z), r2 = tau(g2 Z) + 1/tau(Z), tau = Z + B(Z) a."""
    B1, B1g = pm.basis(P1), pm.basis(mobius(G.g1, P1))
    B2, B2g = pm.basis(P2), pm.basis(mobius(G.g2, P2))
    t1, t2 = pm.tau(P1, a), pm.tau(P2, a)
    return np.concatenate([B1g - B1 / t1[:, None] ** 2, B2g - B2 / t2[:, None] ** 2])


def solve_tau(pm, G, P1, P2, iters=60, lam=1e-3, log=None):
    """Levenberg-Marquardt in the real and imaginary parts of the coefficients,
    with the exact Jacobian. (Until 2026-10-25 a central finite difference with
    step 1e-7 was used: the basis (q/q0)^n reaches ~1e6 on the collocation set
    at n ~ 30, and the difference quotient then failed for K >= 28.)"""
    def F(x):
        a = x[: pm.K + 1] + 1j * x[pm.K + 1:]
        r = tau_residual(pm, G, P1, P2, a)
        return np.concatenate([r.real, r.imag])
    x = np.concatenate([pm.a.real, pm.a.imag])
    r = F(x); L = r @ r
    for it in range(iters):
        Jc = tau_jacobian(pm, G, P1, P2, x[: pm.K + 1] + 1j * x[pm.K + 1:])
        J = np.block([[Jc.real, -Jc.imag], [Jc.imag, Jc.real]])     # d(Re r, Im r)/d(Re a, Im a)
        g = J.T @ r; A = J.T @ J
        while lam < 1e12:
            dx = -np.linalg.solve(A + lam * np.diag(np.diag(A) + 1e-30), g)
            r_new = F(x + dx); L_new = r_new @ r_new
            if L_new < L:
                x, r, L, lam = x + dx, r_new, L_new, max(lam / 3, 1e-15)
                break
            lam *= 4
        if log:
            log(it, L)
        if L < 1e-28 or lam >= 1e12:
            break
    pm.a = x[: pm.K + 1] + 1j * x[pm.K + 1:]
    return L


def solve_mu(pm, G, P1, P2):
    """Linear least squares for b: mu(g1 Z) + mu(Z)/tau(Z) = 1/tau(Z),
    mu(g2 Z) - mu(Z)/tau(Z) = 1."""
    B1, B1g = pm.basis(P1), pm.basis(mobius(G.g1, P1))
    B2, B2g = pm.basis(P2), pm.basis(mobius(G.g2, P2))
    t1, t2 = pm.tau(P1), pm.tau(P2)
    A = np.concatenate([B1g + B1 / t1[:, None], B2g - B2 / t2[:, None]])
    rhs = np.concatenate([1 / t1, np.ones_like(t2)])
    pm.b = np.linalg.lstsq(A, rhs, rcond=None)[0]
    return np.abs(A @ pm.b - rhs).max()


def solve_beta(pm, G, P1, P2, c0):
    """Linear least squares for c_1..c_K (c_0 = c0 is the free constant):
    beta(gZ) - beta(Z) = RHS, beta = -Z + sum c_n basis_n."""
    B1, B1g = pm.basis(P1), pm.basis(mobius(G.g1, P1))
    B2, B2g = pm.basis(P2), pm.basis(mobius(G.g2, P2))
    t1, t2, m1, m2 = pm.tau(P1), pm.tau(P2), pm.mu(P1), pm.mu(P2)
    Zg1, Zg2 = mobius(G.g1, P1), mobius(G.g2, P2)
    # (-Zg + Bg c) - (-Z + B c) = rhs  ->  (Bg - B) c = rhs + Zg - Z
    rhs1 = 2 - 6 * (1 - m1) ** 2 / t1 + Zg1 - P1
    rhs2 = -3 - 6 * m2 ** 2 / t2 + Zg2 - P2
    A = np.concatenate([B1g - B1, B2g - B2])[:, 1:]           # the constant column is zero
    rhs = np.concatenate([rhs1, rhs2])
    c = np.linalg.lstsq(A, rhs, rcond=None)[0]
    pm.cb = np.concatenate([[c0], c])
    return np.abs(A @ c - rhs).max()


def functional_residuals(pm, G, Z):
    """Max residuals of all six functional equations at points Z (anywhere in h)."""
    out = {}
    for name, g in (('g1', G.g1), ('g2', G.g2)):
        Zg = mobius(g, Z)
        t, tg, m, mg, b, bg = pm.tau(Z), pm.tau(Zg), pm.mu(Z), pm.mu(Zg), pm.beta(Z), pm.beta(Zg)
        if name == 'g1':
            out['tau_g1'] = np.abs(tg - (t - 1) / t).max()
            out['mu_g1'] = np.abs(mg - (1 - m) / t).max()
            out['beta_g1'] = np.abs(bg - (b + 2 - 6 * (1 - m) ** 2 / t)).max()
        else:
            out['tau_g2'] = np.abs(tg + 1 / t).max()
            out['mu_g2'] = np.abs(mg - (1 + m / t)).max()
            out['beta_g2'] = np.abs(bg - (b - 3 - 6 * m ** 2 / t)).max()
    return out


def D_function(pm, Z):
    t, m, b = pm.tau(Z), pm.mu(Z), pm.beta(Z)
    return b.imag - 6 * m.imag ** 2 / t.imag


# ---------------------------------------------------------------------------
# Fundamental domain (Ford domain at the distinguished cusp) and validation
# ---------------------------------------------------------------------------

def group_words(G, depth=3):
    """Elements of Delta as products of generators up to the given length (matrices)."""
    gens = [G.g1, np.linalg.inv(G.g1), G.g2, np.linalg.inv(G.g2)]
    words = [np.eye(2, dtype=CPX)]
    frontier = [np.eye(2, dtype=CPX)]
    for _ in range(depth):
        new = []
        for M in frontier:
            for g in gens:
                new.append(g @ M)
        words += new
        frontier = new
    return words


def ford_domain_mask(G, Z, words=None, x0=-0.5):
    """Z in the Ford domain: x0 <= Re Z < x0 + 1 and outside every isometric
    circle |c Z + d| < 1 of the non-parabolic words (c != 0)."""
    words = group_words(G, 3) if words is None else words
    ok = (Z.real >= x0) & (Z.real < x0 + 1)
    for M in words:
        c, d = M[1, 0], M[1, 1]
        if abs(c) > 1e-9:
            ok &= np.abs(c * Z + d) >= 1 - 1e-12
    return ok


def sample_ford(G, n, rng, ymax=2.0):
    words = group_words(G, 3)
    out = []
    while sum(len(o) for o in out) < n:
        Z = rng.uniform(-0.5, 0.5, 4 * n) + 1j * rng.uniform(0.2, ymax, 4 * n)
        out.append(Z[ford_domain_mask(G, Z, words)])
    return np.concatenate(out)[:n]


def edge_pairs(G, n, rng, eps=0.03):
    """Points Z near the Ford-domain boundary together with each generator g for
    which g Z is also (near) the domain: the pairs actually used by the gluing."""
    words = group_words(G, 3)
    Z = sample_ford(G, 20 * n, rng)
    pairs = {}
    for name, g in (('g1', G.g1), ('g2', G.g2), ('g1inv', np.linalg.inv(G.g1)), ('g2inv', np.linalg.inv(G.g2))):
        Zg = mobius(g, Z)
        # near the domain: inside the domain after shifting by the cusp translations
        Zs = Zg - np.round(Zg.real)
        near = Zs.imag > min(G.Z1.imag, G.Z2.imag) - eps
        pairs[name] = Z[near][:n]
    return pairs
