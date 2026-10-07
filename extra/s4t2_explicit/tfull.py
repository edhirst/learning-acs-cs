"""The diffeomorphism F = (Psi_W o pi, phi): Y -> S^4 x T^2 in the three holomorphic charts of Y
(docs/s4t2_explicit_diffeo.md, section 12.6), as maps R^6 -> R^9 = R^5 x R^2 x R^2 (round S^4 x S^1 x S^1):

  'mid' : (z, w, t), the regular part Y_mid (zeta = e^{iz}); a branch of log t is chosen per point
  'c0'  : (x, y, w), the filling U_0 x C_w at t = 0 (zeta = y, t = x y), |t| <= 0.25, |mu| < 1/2 after the T^2-slice
  'ci'  : (x', y', w'), the filling at infinity (zeta' = y', s = x' y'), Y_mid point Theta^{-1}(zeta', w', s)

The complex structure of Y is the standard one in each chart (the paper's gluings are holomorphic).

T^2-trivialisation phi: Y -> R^2/Z^2 (section 11.7 of the notes):
  regular part         phi = (x2, x3)  (lattice coordinates of the point);
  cusp filling 0       phi0 = C - H c(t),  w = C2 e2bar + C3 e3bar,  c = R(t)^-1 e4bar,  H = mu at the T^2-slice point;
  gluing (|t| in [0.3, 0.5]):  phi = (x2, x3) - (1 - chi(|t|)) D0,  D0 = (H - x4_slice) c(t) + (xi2, xi3)(t);
  infinity: the same through Theta (no xi term).
"""
from __future__ import annotations

import math

import torch

from extra.s4t2_explicit import tcore as T, tlift as L
from extra.s4t2_explicit.tcore import CD, RD, _scatter, bump_step, real_coeffs

CHI_PHI = (0.3, 0.5)


def chi_phi(r):
    return bump_step((r - CHI_PHI[0]) / (CHI_PHI[1] - CHI_PHI[0]))


def D0(t, x1, x4, logt):
    zeta, _ = T.to_cusp0(t, x1, x4, logt)
    H = T.CUSP0.mu_torus(zeta, t)
    x4s = torch.log(zeta.abs()) / T.CUSP0.P(t)
    c2, c3 = T.CUSP0.c_tilde(t)
    xi2, xi3 = T.xi23(t)
    return torch.stack([(H - x4s) * c2 + xi2, (H - x4s) * c3 + xi3], -1)


def Dinf(t, x1, x4, logt):
    zeta, s = T.to_cuspinf(t, x1, x4, logt)
    H = T.CUSPINF.mu_torus(zeta, s)
    x4s = torch.log(zeta.abs()) / T.CUSPINF.P(s)
    c2, c3 = T.CUSPINF.c_tilde(s)
    return torch.stack([(H - x4s) * c2, (H - x4s) * c3], -1)


def phi_mid(t, x, logt):
    """phi on the regular part from the lattice coordinates x = (x1, x2, x3, x4) (branch logt)."""
    r = t.abs()
    phi = x[..., 1:3]
    a, b = r < CHI_PHI[1], r > 1 / CHI_PHI[1]
    if a.any():
        corr = (1 - chi_phi(r[a]))[..., None] * D0(t[a], x[a][..., 0], x[a][..., 3], logt[a])
        phi = _scatter(phi, a, phi[a] - corr)
    if b.any():
        corr = (1 - chi_phi(1 / r[b]))[..., None] * Dinf(t[b], x[b][..., 0], x[b][..., 3], logt[b])
        phi = _scatter(phi, b, phi[b] - corr)
    return phi


def torus_embed(phi):
    a = 2 * math.pi * phi
    return torch.stack([torch.cos(a[..., 0]), torch.sin(a[..., 0]), torch.cos(a[..., 1]), torch.sin(a[..., 1])], -1)


def _out(Z1, Z2, X, phi):
    return torch.cat([torch.stack([Z1.real, Z1.imag, Z2.real, Z2.imag, X], -1), torus_embed(phi)], -1)


def F_mid(z, w, t, logt):
    x = T.lattice_coords(z, w, t, logt)
    Z1, Z2, X = L.Psi_W(t, x[..., 0], x[..., 3], logt)
    return _out(Z1, Z2, X, phi_mid(t, x, logt))


def _slice(C, x, y, w):
    """T^2-slice of the filling point (x, y, w): C with w = C2 e2bar + C3 e3bar, and (x_s, y_s) = (x e^{C.al}, y e^{-C.al})."""
    t = x * y
    C2, C3 = real_coeffs(w, C.e2(t), C.e3(t))
    lam = torch.exp(-(C2 * C.al2(t) + C3 * C.al3(t)))
    return t, C2, C3, x / lam, y * lam


def F_chart(C, psi, x, y, w):
    t, C2, C3, xs, ys = _slice(C, x, y, w)
    Z1, Z2, X = psi(xs, ys)
    H = C.mu_chart(xs, ys)
    c2, c3 = C.c_tilde(t)
    return _out(Z1, Z2, X, torch.stack([C2 - H * c2, C3 - H * c3], -1))


def F_c0(x, y, w):
    return F_chart(T.CUSP0, L.Psi_chart0, x, y, w)


def F_ci(x, y, w):
    return F_chart(T.CUSPINF, L.Psi_chartinf, x, y, w)


# ------------------------------------------------------------------ chart transitions (holomorphic)
def c0_to_mid(x, y, w):
    """(x, y, w) in U_0 x C (t = x y != 0) -> the Y_mid point Phi(zeta = y, w, t) as (z, w, t, logt)."""
    t = x * y
    logt = torch.log(t)
    s = logt / (2j * math.pi)
    return -1j * torch.log(y) + math.pi * s * (s - 1), w + 2 * t * s, t, logt


def ci_to_mid(x, y, w):
    """(x', y', w') at infinity (s = x' y' != 0) -> Theta^{-1}: t = 1/s, w = w' t^2, z = z' + i w / t."""
    s = x * y
    t = 1 / s
    wm = w * t ** 2
    return -1j * torch.log(y) + 1j * wm / t, wm, t, -torch.log(s)


# ------------------------------------------------------------------ real-coordinate wrappers (R^6 -> R^9)
def _c(p, i):
    return torch.complex(p[..., i], p[..., i + 1])


def F_mid_real(p):
    """p = (Re z, Im z, Re w, Im w, Re t, Im t); principal branch of log t."""
    t = _c(p, 4)
    return F_mid(_c(p, 0), _c(p, 2), t, torch.log(t))


def F_c0_real(p):
    return F_c0(_c(p, 0), _c(p, 2), _c(p, 4))


def F_ci_real(p):
    return F_ci(_c(p, 0), _c(p, 2), _c(p, 4))
