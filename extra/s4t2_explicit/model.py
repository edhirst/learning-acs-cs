"""Explicit global model of the threefold Y of [HSSS26] (An explicit complex structure on S^4 x T^2),
with the two commuting actions used for the explicit diffeomorphism Y ~ S^4 x T^2:
  T^2   (c2, c3) in R^2/Z^2 : translation by c2 e2 + c3 e3            (paper, Prop. 'torus action')
  S^1_a theta in R/2piZ     : zeta -> e^{i theta} zeta, i.e. z -> z + theta (translation along e1)
See docs/s4t2_explicit_diffeo.md. All maps are the paper's formulas; nothing is fitted.
Coordinates on Y_mid: (zeta, w, t) in C* x C x C*, additive z with zeta = e^{iz}.
"""
from __future__ import annotations

import numpy as np

TWO_PI = 2.0 * np.pi


# ---------------------------------------------------------------- lattice of the fibre X_t
def lattice(t, k=0):
    """Columns e1..e4 of Lambda_t (paper eq. (lattice)), as a (..., 2, 4) complex array; branch log t + 2 pi i k."""
    t = np.asarray(t, dtype=complex)
    lg = np.log(t) + 2j * np.pi * k
    z = np.zeros_like(t)
    r1 = np.stack([z + TWO_PI, -1j * t, -t, -1j * lg], -1)
    r2 = np.stack([z, 1 - t ** 2, 1j * (1 + t ** 2), 2 * t], -1)
    return np.stack([r1, r2], -2)


def real4(v):
    """(..., 2) complex -> (..., 4) real in the complex orientation (Re z, Im z, Re w, Im w)."""
    return np.stack([v[..., 0].real, v[..., 0].imag, v[..., 1].real, v[..., 1].imag], -1)


def lattice_coords(z, w, t, k=0):
    """Real x in R^4 with (z, w) = sum_j x_j e_j(t)."""
    E = lattice(t, k)
    R = np.stack([E[..., 0, :].real, E[..., 0, :].imag, E[..., 1, :].real, E[..., 1, :].imag], -2)  # (...,4,4)
    rhs = real4(np.stack([np.asarray(z, complex), np.asarray(w, complex)], -1))
    return np.linalg.solve(R, rhs[..., None])[..., 0]


# ---------------------------------------------------------------- maps on Y_mid (zeta, w, t)
def g(j, zeta, w, t):
    if j == 2:
        return np.exp(t) * zeta, w + 1 - t ** 2, t
    if j == 3:
        return np.exp(-1j * t) * zeta, w + 1j * (1 + t ** 2), t
    if j == 4:
        return t * zeta, w + 2 * t, t
    raise ValueError(j)


def rot(theta, zeta, w, t):
    """S^1_a: zeta -> e^{i theta} zeta."""
    return np.exp(1j * theta) * zeta, w, t


def torus(c2, c3, zeta, w, t):
    """T^2 action on Y_mid (paper, Prop. torus action)."""
    return np.exp((c2 - 1j * c3) * t) * zeta, w + c2 * (1 - t ** 2) + 1j * c3 * (1 + t ** 2), t


def Phi(zeta, w, t, k=0):
    """Logarithmic transform (paper, Lemma Phi), s = log t / (2 pi i) on branch k."""
    s = (np.log(t) + 2j * np.pi * k) / (2j * np.pi)
    return np.exp(1j * np.pi * s * (s - 1)) * zeta, w + 2 * t * s, t


def Theta(zeta, w, t):
    """Coordinates at infinity (paper, Lemma 'atinfty')."""
    return zeta * np.exp(w / t), w / t ** 2, 1 / t


def sigma_tilde(t, k=0):
    """Additive translation of Phi: (pi s (s - 1), 2 t s)."""
    s = (np.log(t) + 2j * np.pi * k) / (2j * np.pi)
    return np.pi * s * (s - 1), 2 * t * s
