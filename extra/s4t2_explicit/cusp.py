"""The two cusps of Y (paper: Prop. 'Mumford-Tate filling'), on the T^2-slice w = 0 of the filling.

Cusp data (paper, Sections 'The cusp at t = 0' and Lemma 'atinfty'), in the base coordinate t (resp. s = 1/t):
    CUSP0:   alpha2 = t,  alpha3 = -i t,  e2bar = 1 - t^2,  e3bar = i(1 + t^2),  e4bar =  2t,  nu = 1
    CUSPINF: alpha2 = s,  alpha3 =  i s,  e2bar = s^2 - 1,  e3bar = i(1 + s^2),  e4bar = -2s,  nu = e^-2
The T^2 element c multiplies the X-factor by exp(c2 alpha2 + c3 alpha3). On the slice, G4 followed by the
T^2-correction restoring w = 0 is
    psi_hat(zeta, t) = (nu t exp(-(c2 alpha2 + c3 alpha3)) zeta, t),   c = R(t)^{-1} e4bar(t),
and log|zeta| o psi_hat = log|zeta| + P(t), with the paper's height P = log|nu| + log|t| - phi_t(e4bar).

Height (docs/s4t2_explicit_diffeo.md, 5a):  mu = sum_n n w_n,  w_n ~ exp(-2 n log|zeta| + n (n + 1) P),
so mu o psi_hat = mu + 1. In the toric chart U_k (zeta = x^k y^(k+1), t = x y), with j = n - k,
    w_n ~ |x|^(j (j + 1)) |y|^(j (j - 1)) exp((2 k j + j^2 + j) (log|nu| - phi(t))),
smooth across the central fibre; mu = k at the node P_k.

Weight-1 function (notes, section 7):  E = sum_n chi(mu + n) E0(psi_hat^n p),  E0 = y0 + conj(x0) in U_0.
"""
from __future__ import annotations

import numpy as np


def bump_step(u):
    """C-infinity step: 0 for u <= 0, 1 for u >= 1."""
    u = np.clip(u, 0.0, 1.0)
    a = np.where(u > 0, np.exp(-1.0 / np.maximum(u, 1e-300)), 0.0)
    b = np.where(u < 1, np.exp(-1.0 / np.maximum(1 - u, 1e-300)), 0.0)
    return a / (a + b)


def chi(m, inner=0.4, outer=0.9):
    """Cutoff in the height: 1 on |m| <= inner, 0 on |m| >= outer (so sum_n chi(m + n) > 0)."""
    return bump_step((outer - np.abs(m)) / (outer - inner))


def real_coeffs(u, a, b):
    """Real (c2, c3) with c2 a + c3 b = u (a, b R-independent)."""
    det = (a.conj() * b).imag
    return (u.conj() * b).imag / det, (a.conj() * u).imag / det


class Cusp:
    def __init__(self, name, al2, al3, e2, e3, e4, nu, e0coef=1.0):
        self.name, self.al2, self.al3, self.e2, self.e3, self.e4, self.nu = name, al2, al3, e2, e3, e4, nu
        self.e0coef = e0coef        # local model E0 = y0 + e0coef * conj(x0); zero set {mu_node = 0, arg t = pi}

    def _c(self, t):
        t = np.asarray(t, complex)
        return real_coeffs(self.e4(t), self.e2(t), self.e3(t))

    def phi(self, t):
        """phi_t(e4bar), phi_t the R-linear functional with phi_t(e_j bar) = Re alpha_j (j = 2, 3)."""
        t = np.asarray(t, complex)
        c2, c3 = self._c(t)
        return c2 * self.al2(t).real + c3 * self.al3(t).real

    def P(self, t):
        t = np.asarray(t, complex)
        return np.log(abs(self.nu)) + np.log(np.abs(t)) - self.phi(t)

    def multiplier(self, t):
        """psi_hat(zeta, t) = (multiplier(t) zeta, t)."""
        t = np.asarray(t, complex)
        c2, c3 = self._c(t)
        return self.nu * t * np.exp(-(c2 * self.al2(t) + c3 * self.al3(t)))

    def psi_hat(self, zeta, t):
        return self.multiplier(t) * zeta, t

    def mu_torus(self, zeta, t, J=12):
        rho = np.log(np.abs(zeta))
        p = self.P(t)
        nstar = np.round(rho / p - 0.5)
        n = nstar[..., None] + np.arange(-J, J + 1)
        ex = -2 * n * rho[..., None] + n * (n + 1) * p[..., None]
        w = np.exp(ex - ex.max(-1, keepdims=True))
        return (n * w).sum(-1) / w.sum(-1)

    def mu_chart(self, x, y, k, J=8):
        t = x * y
        g = np.log(abs(self.nu)) - self.phi(t)
        X, Y = np.abs(x) ** 2, np.abs(y) ** 2
        num, den = np.zeros_like(X), np.zeros_like(X)
        for j in range(-J, J + 1):
            term = X ** (j * (j + 1) // 2) * Y ** (j * (j - 1) // 2) * np.exp((2 * k * j + j * j + j) * g)
            num, den = num + j * term, den + term
        return k + num / den

    def E0_torus(self, zeta, t):
        """y0 + e0coef conj(x0) in U_0 (y0 = zeta, x0 = t / zeta). With e0coef = |nu| the zero set
        {|y0| = |nu| |x0|, arg t = pi} is exactly where the node height a0|x0|^2 - |y0|^2 vanishes."""
        return zeta * (1 + self.e0coef * np.conj(t) / np.abs(zeta) ** 2)

    def E_torus(self, zeta, t):
        mu = self.mu_torus(zeta, t)
        mt = self.multiplier(t)
        out = np.zeros_like(zeta)
        n0 = -np.round(mu).astype(int)
        for dn in (-1, 0, 1):
            n = n0 + dn
            out = out + chi(mu + n) * self.E0_torus(mt ** n * zeta, t)
        return out


CUSP0 = Cusp('0', lambda t: t, lambda t: -1j * t, lambda t: 1 - t ** 2, lambda t: 1j * (1 + t ** 2),
             lambda t: 2 * t, 1.0)
CUSPINF = Cusp('inf', lambda s: s, lambda s: 1j * s, lambda s: s ** 2 - 1, lambda s: 1j * (1 + s ** 2),
               lambda s: -2 * s, np.exp(-2.0), e0coef=np.exp(-2.0))
