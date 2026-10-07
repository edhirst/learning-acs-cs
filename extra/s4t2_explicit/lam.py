"""[SUPERSEDED 2026-10-07 by tlift.py (fitted arc model, two-phase straightening R, exact linear models,
two-chart lift). Kept as a record; see docs/s4t2_explicit_diffeo.md section 12.]

The orbit-space map Lambda: B -> B_S = S^4/S^1 and the equivariant lift Psi_W: W -> S^4
(docs/s4t2_explicit_diffeo.md, sections 9-10).

Lambda = P_final o moebius o R_gamma o h_B:
  R_gamma(u) = conj(phi(r)) u   rotates the Hopf fibres over |T| = |T(r)| by the fibre phase phi(r) of the zero arc
                                 Gamma_E at radius r = |t| (exact: Newton on E), so that R_gamma h_B(Gamma_E) is the
                                 real great semicircle from A = (0, 1) to B' = (-1, 0); phi = 1 near both ends.
  moebius                        A -> 0, B' -> inf (curve.py); the semicircle -> the ray d = (0, 1, 0).
  P_final(w) = (exp(-nu(log|w|)) rho_O(w/|w|), g(log|w|)),  g = sign * sqrt(1 - exp(-nu)),  nu(tau) = |tau| for
                                 |tau| >= 1 (so |v| = |w| near 0 and 1/|w| near inf: linear, conformal),
                                 rho_O in O(3) with rho_O(d) = (0, 0, -1)  (z1 = 0 meridian).
Lift (equivariant, z1 proportional to E):
  Z1 = sqrt((|v| + v3)/2) E/|E|,  Z2 = (v1 - i v2) (E/|E|) / sqrt(2 (|v| + v3)),  Psi = (Z1, Z2, X).
"""
from __future__ import annotations

import numpy as np

from extra.s4t2_explicit import hopf as H
from extra.s4t2_explicit.cusp import CUSP0, CUSPINF, bump_step
from extra.s4t2_explicit.curve import moebius, R_C0, R_CI
from extra.s4t2_explicit.middle import E_W

RHO_SIGN = -1.0     # det of rho_O: makes Lambda orientation-preserving (check_lam.py, item 2)


def rho_O(w):
    """Orthogonal map with (0, 1, 0) -> (0, 0, -1): (w0, w1, w2) -> (w0, s w2, -w1)  (det = s)."""
    return np.stack([w[..., 0], RHO_SIGN * w[..., 2], -w[..., 1]], -1)


def nu_fn(tau):
    a = np.abs(tau)
    b = bump_step((a - 0.5) / 0.5)
    return (1 - b) * tau ** 2 + b * a


def P_final(w):
    rw = np.linalg.norm(w, axis=-1)
    tau = np.log(rw)
    h = np.exp(-nu_fn(tau))
    g = np.sign(tau) * np.sqrt(np.maximum(1 - h, 0.0))
    v = h[..., None] * rho_O(w / rw[..., None])
    return v, g


# ------------------------------------------------------------------ the zero arc, vectorised Newton
def _G(rep, a, l, r):
    """Invariant version of E on the torus over radius r, unknowns (arg t, log|zeta|) or (arg t, x4)."""
    if rep == 'c0':
        return CUSP0.E_torus(np.exp(l) + 0j, r * np.exp(1j * a)) * np.exp(-l)
    if rep == 'ci':
        return CUSPINF.E_torus(np.exp(l) + 0j, (1 / r) * np.exp(1j * a)) * np.exp(-l)
    return E_W(r * np.exp(1j * a), 0 * a, l, np.log(r) + 1j * a)


_ARC = None


def _arc_table():
    """Gamma_E on a fine log-grid of radii, traced by continuation (each Newton solve starts from the previous
    radius), cached: the initial guesses for gamma_arc."""
    global _ARC
    if _ARC is None:
        rs = np.exp(np.linspace(np.log(1e-9), np.log(1e9), 3601))
        reps = np.where(rs <= R_C0, 'c0', np.where(rs >= R_CI, 'ci', 'mid'))
        A, Lv = np.empty(rs.shape), np.empty(rs.shape)
        prev = None
        for i, (r, rep) in enumerate(zip(rs, reps)):
            if prev is None or prev[0] != rep:
                g = (np.pi, 0.5 * np.log(r) if rep == 'c0' else (0.5 * (np.log(1 / r) - 2.0) if rep == 'ci' else 0.66))
            else:
                g = prev[1:]
            a, l = _newton_damped(rep, np.array([g[0]]), np.array([g[1]]), np.array([r]))
            A[i], Lv[i] = a[0], l[0]
            prev = (rep, a[0], l[0])
        _ARC = (np.log(rs), reps, A, Lv)
    return _ARC


def _newton_damped(key, am, lm, rm, iters=60, h=1e-7, maxstep=0.05):
    for _ in range(iters):
        f0 = _G(key, am, lm, rm)
        fa = (_G(key, am + h, lm, rm) - _G(key, am - h, lm, rm)) / (2 * h)
        fl = (_G(key, am, lm + h, rm) - _G(key, am, lm - h, rm)) / (2 * h)
        det = fa.real * fl.imag - fa.imag * fl.real
        da = (f0.real * fl.imag - f0.imag * fl.real) / det
        dl = (fa.real * f0.imag - fa.imag * f0.real) / det
        sc = np.minimum(1.0, maxstep / np.maximum(np.abs(da), np.abs(dl)).clip(1e-300))
        am, lm = am - sc * da, lm - sc * dl
        if max(np.abs(da).max(), np.abs(dl).max()) < 1e-15:
            break
    return am, lm


def gamma_arc(r):
    """Gamma_E at radii r (array): returns (rep, a, l) with a = arg t (resp. arg s), l = log|zeta| or x4.
    Initial guesses interpolated from the continuation table, damped Newton, convergence asserted."""
    r = np.asarray(r, float)
    lr, reps_t, A_t, L_t = _arc_table()
    rep = np.where(r <= R_C0, 'c0', np.where(r >= R_CI, 'ci', 'mid'))
    a, l = np.full(r.shape, np.pi), np.zeros(r.shape)
    for key in ('c0', 'mid', 'ci'):
        m = rep == key
        if not m.any():
            continue
        mt = reps_t == key                                  # interpolate within the same representation only
        a[m] = np.interp(np.log(r[m]), lr[mt], A_t[mt])
        l[m] = np.interp(np.log(r[m]), lr[mt], L_t[mt])
        am, lm = _newton_damped(key, a[m], l[m], r[m])
        scale = np.abs(_G(key, am, lm + 0.25, r[m])) + np.abs(_G(key, am + 0.5, lm, r[m]))
        res = np.abs(_G(key, am, lm, r[m])) / scale
        assert np.all(res < 1e-8) and np.all(np.abs(am - np.pi) < 1e-6), \
            f'Gamma_E Newton failed ({key}): max residual {res.max():.1e}, max |arg - pi| {np.abs(am - np.pi).max():.1e}'
        a[m], l[m] = am, lm
    return rep, a, l


def gamma_u(r):
    rep, a, l = gamma_arc(r)
    u = np.empty(np.shape(r) + (2,), complex)
    for key in ('c0', 'mid', 'ci'):
        m = rep == key
        if not m.any():
            continue
        if key == 'c0':
            u[m] = H.h_B_c0(np.exp(l[m]) + 0j, r[m] * np.exp(1j * a[m]))
        elif key == 'ci':
            u[m] = H.h_B_ci(np.exp(l[m]) + 0j, (1 / r[m]) * np.exp(1j * a[m]))
        else:
            u[m] = H.h_B_mid(r[m] * np.exp(1j * a[m]), l[m])
    return u


def fibre_phase(r):
    """phi(r): R_gamma = multiplication by conj(phi) maps h_B(Gamma_E(r)) onto the real great semicircle A -> B'."""
    u = gamma_u(np.asarray(r, float))
    return np.where(np.asarray(r) <= 1, u[..., 1] / np.abs(u[..., 1]), -u[..., 0] / np.abs(u[..., 0]))


def Lambda(u, r):
    """Orbit-space map, for a point with Hopf image u = h_B(b) over base radius r = |t|."""
    ph = fibre_phase(r)
    w = moebius(np.conj(ph)[..., None] * u)
    return P_final(w)


def lift(v, X, E):
    s = np.linalg.norm(v, axis=-1) + v[..., 2]
    ph = E / np.abs(E)
    Z1 = np.sqrt(s / 2) * ph
    Z2 = (v[..., 0] - 1j * v[..., 1]) * ph / np.sqrt(2 * s)
    return Z1, Z2, X
