"""[SUPERSEDED 2026-10-07: the theta/log blend below is DISCONTINUOUS on the blend tori (the exact
log-ratio is direction-dependent at z_ref); the construction uses tlift.E_W (linear blend). Kept as a record;
see docs/s4t2_explicit_diffeo.md section 12.2.]

The weight-1 function E on all of W (docs/s4t2_explicit_diffeo.md, section 7).

In lattice coordinates (t, x1, x4) of the regular part, every weight-1 function is E = e^{2 pi i x1} F(t, x4) with
    F(sigma + 1, x4) = e^{2 pi i x4} F(sigma, x4),  F(sigma, x4 + 1) = F(sigma, x4),   sigma = arg t / 2 pi
(continuous branch). Theta_{z0} below is the standard such section with one simple zero at z0 = (sigma0, x0).
    |t| <= R0:        E = E0_W                (cusp-0 formula, smooth through the central fibre)
    |t| >= 1/R0:      E = Einf_W              (cusp-inf formula)
    R0 < |t| < 1/R0:  F = Theta_{z(r)} exp(2 pi i (a sigma + b x4) + l(r)),  with the zero z(r) and the periodic
                      part l(r) interpolated between the two cusp functions, read at the plateau radii RP, 1/RP.
"""
from __future__ import annotations

import numpy as np

from extra.s4t2_explicit.cusp import bump_step
from extra.s4t2_explicit.glue import E0_W, Einf_W

R0, RP, R1 = 0.6, 0.7, 0.8            # cusp formula below R0; plateau radius RP; blend on [R1, 1/R1]
KAPPA = 0.5


def theta_std(s, x, kappa=KAPPA, N=12):
    """f(s, x) = sum_n exp(-pi (s + n)^2 / kappa) e^{-2 pi i n x}: f(s+1, x) = e^{2 pi i x} f, zero at (1/2, 1/2)."""
    s, x = np.asarray(s, float), np.asarray(x, float)
    sb = np.floor(s + 0.5)                             # reduce s for accuracy, using the twisted periodicity
    sr = s - sb
    out = np.zeros(np.broadcast(sr, x).shape, complex)
    for n in range(-N, N + 1):
        out = out + np.exp(-np.pi * (sr + n) ** 2 / kappa) * np.exp(-2j * np.pi * n * x)
    return out * np.exp(2j * np.pi * sb * x)


def theta_at(s, x, s0, x0):
    """Section with the twisted periodicity and a simple zero at (s0, x0)."""
    return np.exp(2j * np.pi * s * (x0 - 0.5)) * theta_std(s - s0 + 0.5, x - x0 + 0.5)


def F_cusp(F, r, s, x):
    """F(t, x4) = E(t, x1 = 0, x4) along the continuous branch log t = log r + 2 pi i s."""
    return F(r * np.exp(2j * np.pi * s), 0 * s, x, np.log(r) + 2j * np.pi * s)


class Reference:
    """Zero z, windings (a, b) and periodic part l of F / Theta_z on the torus |t| = r (grid n x n)."""

    def __init__(self, F, r, n=256, near=None):
        self.r = r
        S, X = np.meshgrid(np.arange(n) / n, np.arange(n) / n, indexing='ij')
        G = F_cusp(F, r, S, X)
        i0, j0 = np.unravel_index(np.argmin(np.abs(G)), G.shape)
        z = np.array([S[i0, j0], X[i0, j0]])
        for _ in range(30):                            # Newton on the simple zero
            h = 1e-7
            f0 = F_cusp(F, r, z[0], z[1])
            fs = (F_cusp(F, r, z[0] + h, z[1]) - F_cusp(F, r, z[0] - h, z[1])) / (2 * h)
            fx = (F_cusp(F, r, z[0], z[1] + h) - F_cusp(F, r, z[0], z[1] - h)) / (2 * h)
            Jm = np.array([[fs.real, fx.real], [fs.imag, fx.imag]])
            z = z - np.linalg.solve(Jm, np.array([f0.real, f0.imag]))
        if near is not None:                           # integer representative closest to `near`
            z = z - np.round(z - near)
        self.z = z
        self.zero_residual = abs(F_cusp(F, r, z[0], z[1]))
        R = G / theta_at(S, X, z[0], z[1])
        ph = np.angle(R)
        self.a = int(np.round((np.unwrap(ph[:, n // 3])[-1] - ph[0, n // 3] + np.angle(R[0, n // 3] / R[-1, n // 3])) / (2 * np.pi)))
        self.b = int(np.round((np.unwrap(ph[n // 3, :])[-1] - ph[n // 3, 0] + np.angle(R[n // 3, 0] / R[n // 3, -1])) / (2 * np.pi)))
        L = np.log(np.abs(R)) + 1j * ph - 2j * np.pi * (self.a * S + self.b * X)
        # continuous periodic branch of the imaginary part: unwrap along both directions from the origin
        im = np.unwrap(np.unwrap(L.imag, axis=0), axis=1)   # column 0 first, then every row from it
        self.periodicity_defect = max(abs(im[0, :] - im[-1, :]).max(), abs(im[:, 0] - im[:, -1]).max())
        self.l_hat = np.fft.fft2(np.log(np.abs(R)) + 1j * im)  # periodic smooth -> spectral interpolation
        self.n = n
        self.R = R

    def l_exact(self, F, s, x):
        """Exact log of F / (Theta_z e^{2 pi i (a s + b x)}) at (s, x): the spectral part only selects the
        2 pi branch (its error, ~1e-2, is far below pi), so the blend is exactly F at its ends."""
        ratio = F_cusp(F, self.r, s, x) / (theta_at(s, x, self.z[0], self.z[1])
                                            * np.exp(2j * np.pi * (self.a * s + self.b * x)))
        approx = self.l(s, x)
        ang = np.angle(ratio)
        k = np.round((approx.imag - ang) / (2 * np.pi))
        return np.log(np.abs(ratio)) + 1j * (ang + 2 * np.pi * k)

    def l(self, s, x):
        """Spectral evaluation of the periodic part at arbitrary (s, x)."""
        n = self.n
        k = np.fft.fftfreq(n, 1.0 / n)
        es = np.exp(2j * np.pi * np.multiply.outer(np.asarray(s), k))
        ex = np.exp(2j * np.pi * np.multiply.outer(np.asarray(x), k))
        return np.einsum('...i,ij,...j->...', es, self.l_hat, ex) / n ** 2


REF0 = REFINF = None


def references():
    global REF0, REFINF
    if REF0 is None:
        REF0 = Reference(E0_W, RP)
        REFINF = Reference(Einf_W, 1 / RP, near=REF0.z)
    return REF0, REFINF


def h0(r):
    return r + (RP - r) * bump_step((r - R0) / (R1 - R0))


def hinf(r):
    return r + (1 / RP - r) * bump_step((1 / R0 - r) / (1 / R0 - 1 / R1))


def F_mid(r, s, x):
    """The middle section on R0 <= r <= 1/R0."""
    r0, ri = references()
    lam = bump_step((np.log(r) - np.log(R1)) / (np.log(1 / R1) - np.log(R1)))
    out = np.empty(np.broadcast(r, s, x).shape, complex)
    r, s, x, lam = np.broadcast_arrays(r, s, x, lam)
    pure0, pureinf = lam <= 0, lam >= 1
    mixed = ~(pure0 | pureinf)
    if pure0.any():
        out[pure0] = F_cusp(E0_W, h0(r[pure0]), s[pure0], x[pure0])
    if pureinf.any():
        out[pureinf] = F_cusp(Einf_W, hinf(r[pureinf]), s[pureinf], x[pureinf])
    if mixed.any():
        lm, sm, xm = lam[mixed], s[mixed], x[mixed]
        z = (1 - lm)[:, None] * r0.z + lm[:, None] * ri.z
        ell = (1 - lm) * r0.l_exact(E0_W, sm, xm) + lm * ri.l_exact(Einf_W, sm, xm)
        out[mixed] = theta_at(sm, xm, z[:, 0], z[:, 1]) * np.exp(2j * np.pi * (r0.a * sm + r0.b * xm) + ell)
    return out


def E_W(t, x1, x4, logt):
    """The global weight-1 function on W, in lattice coordinates (regular part); continuous branch logt."""
    r = np.abs(t)
    s = logt.imag / (2 * np.pi)
    out = np.empty(np.broadcast(t, x1, x4).shape, complex)
    t, x1, x4, logt, r, s = np.broadcast_arrays(t, x1, x4, logt, r, s)
    a, b = r <= R0, r >= 1 / R0
    m = ~(a | b)
    if a.any():
        out[a] = E0_W(t[a], x1[a], x4[a], logt[a])
    if b.any():
        out[b] = Einf_W(t[b], x1[b], x4[b], logt[b])
    if m.any():
        out[m] = np.exp(2j * np.pi * x1[m]) * F_mid(r[m], s[m], x4[m])
    return out
