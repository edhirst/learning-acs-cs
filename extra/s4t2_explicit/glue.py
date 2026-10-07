"""Conversions between the descriptions of W = Y/T^2 (docs/s4t2_explicit_diffeo.md):
  * regular part: lattice coordinates of Y_mid, the T^2-slice {x2 = x3 = 0}, i.e. (t, x1, x4) with a branch
    of log t (x1 -> x1 - x4 when log t -> log t + 2 pi i; x1, x4 mod 1);
  * cusp 0: slice {w = 0} of the filling U_0, glued to Y_mid by the log transform Phi (q ~ Phi(iota_0 q));
  * cusp inf: slice {w' = 0} of U_inf, glued by Theta (q ~ Theta^{-1} q).
All T^2-moves use the paper's action (translation by c2 e2 + c3 e3 on Y_mid, multiplier exp(c2 al2 + c3 al3)
on the filling); S^1_a is x1 -> x1 + theta / 2 pi throughout.
"""
from __future__ import annotations

import numpy as np

from extra.s4t2_explicit.cusp import CUSP0, CUSPINF, real_coeffs

TWO_PI = 2 * np.pi


def e_vectors(t, logt):
    """(z, w)-components of e1..e4 at t, with the given value of log t."""
    z = np.zeros_like(t)
    ez = [z + TWO_PI, -1j * t, -t, -1j * logt]
    ew = [z, 1 - t ** 2, 1j * (1 + t ** 2), 2 * t]
    return ez, ew


def mid_point(t, x1, x4, logt):
    """Point (z, w) of Y_mid with lattice coordinates (x1, 0, 0, x4)."""
    ez, ew = e_vectors(t, logt)
    return x1 * ez[0] + x4 * ez[3], x1 * ew[0] + x4 * ew[3]


def to_cusp0(t, x1, x4, logt):
    """The cusp-0 slice point (zeta, t) representing the same orbit of T^2 in W.
    q ~ Phi(q): the Y_mid point is Phi(q), so q = Phi^{-1}(point) = point - sigma_tilde (additive), then move
    by T^2 to w = 0."""
    z, w = mid_point(t, x1, x4, logt)
    s = logt / (2j * np.pi)
    z, w = z - np.pi * s * (s - 1), w - 2 * t * s                 # Phi^{-1}
    c2, c3 = real_coeffs(-w, CUSP0.e2(t), CUSP0.e3(t))            # T^2-move to w = 0 (filling = Y_mid here)
    z = z + c2 * (-1j * t) + c3 * (-t)                            # translation by c2 e2 + c3 e3 (z-part)
    return np.exp(1j * z), t


def to_cuspinf(t, x1, x4, logt):
    """The cusp-inf slice point (zeta', s): Theta in additive form, then T^2-move to w' = 0."""
    z, w = mid_point(t, x1, x4, logt)
    s = 1 / t
    z2, w2 = z - 1j * w / t, w / t ** 2                           # Theta: zeta' = zeta e^{w/t}, w' = w/t^2
    c2, c3 = real_coeffs(-w2, CUSPINF.e2(s), CUSPINF.e3(s))
    z2 = z2 - 1j * (c2 * CUSPINF.al2(s) + c3 * CUSPINF.al3(s))     # multiplier exp(c2 al2 + c3 al3) on zeta'
    return np.exp(1j * z2), s


def E0_W(t, x1, x4, logt):
    return CUSP0.E_torus(*to_cusp0(t, x1, x4, logt))


def Einf_W(t, x1, x4, logt):
    return CUSPINF.E_torus(*to_cuspinf(t, x1, x4, logt))
