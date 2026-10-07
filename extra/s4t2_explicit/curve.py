"""[SUPERSEDED 2026-10-07: Gamma_E of the old (discontinuous) E; moebius is ported to tcore.moebius.
See docs/s4t2_explicit_diffeo.md section 12.]

The zero arc Gamma_E of E in the orbit space B, its image under h_B, and the Moebius map
mu: S^3 -> R^3 u {inf} with h_B(b0) = A = (0, 1) -> 0 and h_B(b_inf) = B' = (-1, 0) -> inf
(stereographic projection from B', then translation; conformal)."""
from __future__ import annotations

import numpy as np

from extra.s4t2_explicit import hopf as H
from extra.s4t2_explicit.cusp import CUSP0, CUSPINF
from extra.s4t2_explicit.middle import E_W

R_C0, R_CI = 0.6, 1 / 0.6


def _newton2(f, z, iters=40, h=1e-7):
    for _ in range(iters):
        f0 = f(*z)
        fa = (f(z[0] + h, z[1]) - f(z[0] - h, z[1])) / (2 * h)
        fb = (f(z[0], z[1] + h) - f(z[0], z[1] - h)) / (2 * h)
        J = np.array([[fa.real, fb.real], [fa.imag, fb.imag]])
        step = np.linalg.solve(J, np.array([f0.real, f0.imag]))
        z = z - step
        if np.abs(step).max() < 1e-15:
            break
    return z, abs(f(*z))


def gamma_point(r, guess=None):
    """The zero of E on the orbit-space torus |t| = r. Returns (rep, data, u) with u = h_B(point)."""
    if r <= R_C0:
        C = CUSP0
        f = lambda a, l: C.E_torus(np.array([np.exp(l) + 0j]), np.array([r * np.exp(1j * a)]))[0] * np.exp(-l)
        z0 = guess if guess is not None else np.array([np.pi, 0.5 * np.log(r)])
        z, res = _newton2(f, z0)
        t = r * np.exp(1j * z[0]); zeta = np.exp(z[1]) + 0j
        return 'c0', (zeta, t, res), H.h_B_c0(np.array([zeta]), np.array([t]))[0]
    if r >= R_CI:
        s_abs = 1 / r
        C = CUSPINF
        f = lambda a, l: C.E_torus(np.array([np.exp(l) + 0j]), np.array([s_abs * np.exp(1j * a)]))[0] * np.exp(-l)
        z0 = guess if guess is not None else np.array([np.pi, 0.5 * np.log(s_abs) + 0.5 * np.log(abs(C.nu))])
        z, res = _newton2(f, z0)
        s = s_abs * np.exp(1j * z[0]); zp = np.exp(z[1]) + 0j
        return 'ci', (zp, s, res), H.h_B_ci(np.array([zp]), np.array([s]))[0]
    f = lambda a, x4: E_W(np.array([r * np.exp(1j * a)]), np.array([0.0]), np.array([x4]),
                          np.array([np.log(r) + 1j * a]))[0]
    z0 = guess if guess is not None else np.array([np.pi - 1e-9, 0.7])
    z, res = _newton2(f, z0)
    t = r * np.exp(1j * z[0])
    return 'mid', (t, z[1], res), H.h_B_mid(np.array([t]), np.array([z[1]]))[0]


def moebius(u):
    """mu(u) = stereographic projection of u in S^3 from B' = (-1, 0), minus the image (0, 1, 0) of A."""
    P = np.stack([u[..., 0].real, u[..., 0].imag, u[..., 1].real, u[..., 1].imag], -1)
    # 1 + Re u1, computed stably near B' (Re u1 ~ -1): = ((Im u1)^2 + |u2|^2) / (1 - Re u1) on |u| = 1
    den = (P[..., 1] ** 2 + P[..., 2] ** 2 + P[..., 3] ** 2) / (1 - P[..., 0])
    w = P[..., 1:] / den[..., None]
    return w - np.array([0.0, 1.0, 0.0])
