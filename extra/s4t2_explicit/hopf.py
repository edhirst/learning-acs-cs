"""Explicit identification h_B: B = W/S^1_a -> S^3 subset C^2 (Hopf), docs/s4t2_explicit_diffeo.md section 9.

Orbit-space representations of a point of B:
  'c0'  : cusp-0 slice point (zeta, t), 0 < |t| < 0.6         (x4_mid = log|zeta| / P0(t) + xi4(t))
  'mid' : (t, x4) with x4 = x4_mid, the lattice coordinate (principal branch of log t)
  'ci'  : cusp-inf slice point (zeta', s), 0 < |s| < 0.6      (x4_mid = -log|zeta'| / Pinf(s))
Fibre coordinates (mod 1):
  a1 = mu_new^0 (|t| <= RA), blended to x4_mid - xi4(t) - 1/2 (|t| >= RB);
  a2 = mu_new^inf (|s| <= RA), blended to -x4_mid - 1/2 (|s| >= RB);   a2' = -a2 + lam(|t|) delta(t),
  delta = arg t / 2 pi - xi4(t) (single-valued: both have monodromy +1).
Hopf point: |t| <= 1: u2 = e^{2 pi i a1} / sqrt(1 + |T|^2), u1 = T u2;  |t| >= 1: u1 = e^{2 pi i a2'} / sqrt(1 + |T|^-2),
u2 = u1 / T;  T(t) = t k(|t|), k = 4 pi near 0, e^2/(4 pi) near infinity (conformal, scale 2 pi, at b0 and b_inf).
"""
from __future__ import annotations

import numpy as np

from extra.s4t2_explicit.cusp import CUSP0, CUSPINF, bump_step
from extra.s4t2_explicit.model import lattice_coords, sigma_tilde

RA, RB = 0.3, 0.5                  # cusp height -> lattice coordinate blend, in |t| (resp. |s|)
NODE_IN, NODE_OUT = 0.02, 0.05     # node correction radius, in |v| = a0 |x0|^2 + |y0|^2
K0, KINF = 4 * np.pi, np.exp(2.0) / (4 * np.pi)
LK0, LK1 = -4.0, 4.0               # log k interpolated in log|t| on [LK0, LK1]; |T| = r k(r) must increase:
                                   # max |d log k / d log r| = 2 log(16 pi^2 e^-2) / (LK1 - LK0) = 0.77 < 1


def wrap(d):
    """Representative of d mod 1 in (-1/2, 1/2]."""
    return d - np.round(d)


def xi4(t):
    """e4-coordinate of the log-transform translation sigma_tilde(t) (principal branch)."""
    sz, sw = sigma_tilde(t)
    return lattice_coords(sz, sw, t)[..., 3]


def k_of(r):
    lam = bump_step((np.log(r) - LK0) / (LK1 - LK0))
    return np.exp((1 - lam) * np.log(K0) + lam * np.log(KINF))


def T_of(t):
    return t * k_of(np.abs(t))


# ------------------------------------------------------------------ node-corrected cusp heights
def a0_of(C):
    return abs(C.nu) ** 2          # mu = a0 |x0|^2 - |y0|^2 + O(4) at the node P_0


def mu_new(C, zeta, t):
    """mu, replaced near every node by exactly k + a0 |x0|^2 - |y0|^2 (psi_hat-equivariantly)."""
    mu = C.mu_torus(zeta, t)
    n = np.round(mu)
    z0 = zeta * C.multiplier(t) ** (-n)                 # psi_hat^{-n}: near P_0 when mu ~ n
    x0, y0 = t / z0, z0
    a0 = a0_of(C)
    m0 = a0 * np.abs(x0) ** 2 - np.abs(y0) ** 2
    rad = a0 * np.abs(x0) ** 2 + np.abs(y0) ** 2
    c = 1 - bump_step((rad - NODE_IN) / (NODE_OUT - NODE_IN))
    return n + (mu - n) + c * (m0 - (mu - n))


# ------------------------------------------------------------------ fibre coordinates
def a1_c0(zeta, t):
    """a1 from the cusp-0 representation (any 0 < |t| < 0.6)."""
    r = np.abs(t)
    mn = mu_new(CUSP0, zeta, t)
    x4sl = np.log(np.abs(zeta)) / CUSP0.P(t)
    lam = bump_step((r - RA) / (RB - RA))
    return mn + lam * wrap(x4sl - 0.5 - mn)


def a1_mid(t, x4):
    return x4 - xi4(t) - 0.5


def a2_ci(zp, s):
    r = np.abs(s)
    mn = mu_new(CUSPINF, zp, s)
    x4sl = np.log(np.abs(zp)) / CUSPINF.P(s)
    lam = bump_step((r - RA) / (RB - RA))
    return mn + lam * wrap(x4sl - 0.5 - mn)


def a2_mid(t, x4):
    return -x4 - 0.5


def delta(t):
    return np.angle(t) / (2 * np.pi) - xi4(t)


def lam_inf(r):
    """1 for |t| <= 1.2, 0 for |t| >= 1.6."""
    return 1 - bump_step((r - 1.2) / (1.6 - 1.2))


# ------------------------------------------------------------------ the Hopf point
def hopf_from_a1(t, a1):
    T = T_of(t)
    u2 = np.exp(2j * np.pi * a1) / np.sqrt(1 + np.abs(T) ** 2)
    return np.stack([T * u2, u2], -1)


def hopf_from_a2(t, a2):
    """|t| >= 1 side (t may be inf: pass s = 1/t through hopf_from_a2_s)."""
    lam = lam_inf(np.abs(t))
    a2p = -a2 + lam * delta(t) + (1 - lam) * 0.5          # b_inf -> u = (-1, 0) (notes, section 10)
    T = T_of(t)
    u1 = np.exp(2j * np.pi * a2p) / np.sqrt(1 + np.abs(T) ** -2)
    return np.stack([u1, u1 / T], -1)


def hopf_from_a2_s(s, a2):
    """Same, written in s = 1/t (valid at s = 0, where lam_inf = 0)."""
    r = 1 / np.maximum(np.abs(s), 1e-300)
    lam = lam_inf(r)
    tt = np.where(np.abs(s) > 0, 1 / np.where(np.abs(s) > 0, s, 1), 1.0)
    a2p = -a2 + np.where(lam > 0, lam * delta(tt), 0.0) + (1 - lam) * 0.5
    invT = s / k_of(r)                                  # 1/T = s / k(|t|)
    u1 = np.exp(2j * np.pi * a2p) / np.sqrt(1 + np.abs(invT) ** 2)
    return np.stack([u1, u1 * invT], -1)


def h_B_mid(t, x4):
    """h_B for a point in the middle representation (any t with 0.5 <= |t| <= 2)."""
    r = np.abs(t)
    out = np.empty(t.shape + (2,), complex)
    lo = r <= 1
    if lo.any():
        out[lo] = hopf_from_a1(t[lo], a1_mid(t[lo], x4[lo]))
    if (~lo).any():
        out[~lo] = hopf_from_a2(t[~lo], a2_mid(t[~lo], x4[~lo]))
    return out


def h_B_c0(zeta, t):
    return hopf_from_a1(t, a1_c0(zeta, t))


def h_B_ci(zp, s):
    return hopf_from_a2_s(s, a2_ci(zp, s))
