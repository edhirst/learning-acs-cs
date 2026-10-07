"""Torch port (float64 / complex128, differentiable to any order) of the verified building blocks of the explicit
diffeomorphism Y ~ S^4 x T^2 (docs/s4t2_explicit_diffeo.md, section 11): model.py, cusp.py, glue.py, hopf.py,
curve.moebius and lam.P_final. Formulas are transcribed one-to-one; check_tcore.py compares every function with the
numpy original to round-off.

Conventions for automatic differentiation:
  * every function takes tensors of any batch shape; complex inputs are complex128;
  * branches are evaluated only on their own index set (_scatter), so that out-of-range formulas (e.g. cusp
    formulas near |t| = 1, 1/x on the central fibre) never enter the graph;
  * integer choices (rounding of heights to the nearest node) are detached: they select a representative, and the
    result does not depend on the choice where it is used.
"""
from __future__ import annotations

import math

import torch

RD, CD = torch.float64, torch.complex128
TWO_PI = 2.0 * math.pi


def cplx(x):
    return x if torch.is_complex(x) else x.to(CD)


def _scatter(out, mask, val):
    """out[mask] = val, out of place for autograd (index_put on a clone)."""
    if not bool(mask.any()):
        return out
    out = out.clone()
    out[mask] = val
    return out


def abs2(z):
    """|z|^2 as Re^2 + Im^2: polynomial, so its AD derivatives are correct to all orders at z = 0 (z.abs() ** 2 gives
    the right first derivative there but a wrong second one: it differentiates through the non-smooth |z|)."""
    return z.real ** 2 + z.imag ** 2


def where_safe(mask, x, safe):
    """x where mask, else the constant `safe` (keeps unused branches finite)."""
    return torch.where(mask, x, torch.as_tensor(safe, dtype=x.dtype))


# ------------------------------------------------------------------ cutoffs
def bump_step(u):
    """C-infinity step: 0 for u <= 0, 1 for u >= 1 (identical to cusp.bump_step: exp(-1/u) underflows below 1e-3)."""
    u = u.clamp(0.0, 1.0)
    pa, pb = u > 1e-3, (1.0 - u) > 1e-3
    a = torch.where(pa, torch.exp(-1.0 / where_safe(pa, u, 1.0)), torch.zeros_like(u))
    b = torch.where(pb, torch.exp(-1.0 / where_safe(pb, 1.0 - u, 1.0)), torch.zeros_like(u))
    return a / (a + b)


def chi(m, inner=0.4, outer=0.9):
    return bump_step((outer - m.abs()) / (outer - inner))


def real_coeffs(u, a, b):
    """Real (c2, c3) with c2 a + c3 b = u."""
    det = (a.conj() * b).imag
    return (u.conj() * b).imag / det, (a.conj() * u).imag / det


# ------------------------------------------------------------------ lattice of the regular fibres (model.py)
def e_vectors(t, logt):
    z = torch.zeros_like(t)
    ez = [z + TWO_PI, -1j * t, -t, -1j * logt]
    ew = [z, 1 - t ** 2, 1j * (1 + t ** 2), 2 * t]
    return ez, ew


def _det3(a, b, c, d, e, f, g, h, i):
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def lattice_coords(z, w, t, logt):
    """Real x in R^4 with (z, w) = sum_j x_j e_j(t) (branch logt).
    e1 = (2 pi, 0) is real, so (Im z, Re w, Im w) is a 3 x 3 real system for (x2, x3, x4), solved by Cramer's rule
    (elementary operations only: torch.linalg.solve has wrong nested forward-mode derivatives in torch 2.10, which
    would corrupt D^2 F); then x1 = (Re z - sum_{j>1} x_j Re e_j^z) / 2 pi."""
    ez, ew = e_vectors(t, logt)
    a = [ez[j].imag for j in (1, 2, 3)]
    b = [ew[j].real for j in (1, 2, 3)]
    c = [ew[j].imag for j in (1, 2, 3)]
    r = [z.imag, w.real, w.imag]
    det = _det3(a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2])
    x2 = _det3(r[0], a[1], a[2], r[1], b[1], b[2], r[2], c[1], c[2]) / det
    x3 = _det3(a[0], r[0], a[2], b[0], r[1], b[2], c[0], r[2], c[2]) / det
    x4 = _det3(a[0], a[1], r[0], b[0], b[1], r[1], c[0], c[1], r[2]) / det
    x1 = (z.real - x2 * ez[1].real - x3 * ez[2].real - x4 * ez[3].real) / TWO_PI
    return torch.stack([x1, x2, x3, x4], -1)


def point_from_lattice(t, logt, x):
    """(z, w) = sum_j x_j e_j(t)."""
    ez, ew = e_vectors(t, logt)
    z = sum(x[..., j] * ez[j] for j in range(4))
    w = sum(x[..., j] * ew[j] for j in range(4))
    return z, w


def sigma_tilde(t, logt):
    s = logt / (2j * math.pi)
    return math.pi * s * (s - 1), 2 * t * s


# ------------------------------------------------------------------ cusps (cusp.py)
class TCusp:
    def __init__(self, name, al2, al3, e2, e3, e4, nu, e0coef):
        self.name, self.al2, self.al3, self.e2, self.e3, self.e4 = name, al2, al3, e2, e3, e4
        self.nu, self.e0coef = nu, e0coef
        self.lognu = math.log(abs(nu))
        self.a0 = abs(nu) ** 2

    def c_tilde(self, t):
        return real_coeffs(self.e4(t), self.e2(t), self.e3(t))

    def phi(self, t):
        c2, c3 = self.c_tilde(t)
        return c2 * self.al2(t).real + c3 * self.al3(t).real

    def P(self, t):
        return self.lognu + torch.log(t.abs()) - self.phi(t)

    def m(self, t):
        """psi_hat(zeta) = t m(t) zeta; m is smooth and nonvanishing at t = 0."""
        c2, c3 = self.c_tilde(t)
        return self.nu * torch.exp(-(c2 * self.al2(t) + c3 * self.al3(t)))

    def multiplier(self, t):
        return t * self.m(t)

    def mu_torus(self, zeta, t, J=12):
        rho = torch.log(zeta.abs())
        p = self.P(t)
        nstar = torch.round(rho / p - 0.5).detach()
        n = nstar[..., None] + torch.arange(-J, J + 1, dtype=RD)
        ex = -2 * n * rho[..., None] + n * (n + 1) * p[..., None]
        return (n * torch.softmax(ex, -1)).sum(-1)

    def mu_chart(self, x, y, J=8):
        """Height in the toric chart U_0 (zeta = y, t = x y); smooth through the central fibre."""
        g = self.lognu - self.phi(x * y)
        X, Y = abs2(x), abs2(y)
        num, den = 0.0, 0.0
        for j in range(-J, J + 1):
            a, b = j * (j + 1) // 2, j * (j - 1) // 2
            term = torch.exp((j * j + j) * g)
            if a:
                term = term * X ** a
            if b:
                term = term * Y ** b
            num, den = num + j * term, den + term
        return num / den

    def E0(self, x0, y0, sign=1.0):
        return y0 + sign * self.e0coef * x0.conj()

    def E_torus(self, zeta, t, sign=1.0):
        mu = self.mu_torus(zeta, t)
        mt = self.multiplier(t)
        n0 = -torch.round(mu).detach()
        out = torch.zeros_like(zeta)
        for dn in (-1, 0, 1):
            n = n0 + dn
            Z = torch.exp(n * torch.log(mt)) * zeta
            out = out + chi(mu + n) * Z * (1 + sign * self.e0coef * t.conj() / abs2(Z))
        return out

    def E_chart(self, x, y, sign=1.0):
        """The same function in U_0 coordinates (valid for |mu| < 0.6, including the central fibre t = 0)."""
        m = self.m(x * y)
        mu = self.mu_chart(x, y)
        c0, cp, cm = chi(mu), chi(mu + 1), chi(mu - 1)
        out = c0 * self.E0(x, y, sign)
        ys = torch.where(cp > 0, y, torch.ones_like(y))          # n = +1: psi_hat p = (1/(m y), m x y^2) in U_0
        out = out + cp * self.E0(1 / (m * ys), m * x * ys ** 2, sign)
        xs = torch.where(cm > 0, x, torch.ones_like(x))          # n = -1: psi_hat^-1 p = (m x^2 y, 1/(m x)) in U_0
        out = out + cm * self.E0(m * xs ** 2 * y, 1 / (m * xs), sign)
        return out

    def reduce(self, zeta, t):
        """Representative zeta0 = zeta * multiplier^-n with n = round(mu) (detached) and U_0 coordinates (x0, y0)."""
        mu = self.mu_torus(zeta, t)
        n = torch.round(mu).detach()
        z0 = zeta * torch.exp(-n * torch.log(self.multiplier(t)))
        return n, t / z0, z0


CUSP0 = TCusp('0', lambda t: t, lambda t: -1j * t, lambda t: 1 - t ** 2, lambda t: 1j * (1 + t ** 2),
              lambda t: 2 * t, 1.0, 1.0)
CUSPINF = TCusp('inf', lambda s: s, lambda s: 1j * s, lambda s: s ** 2 - 1, lambda s: 1j * (1 + s ** 2),
                lambda s: -2 * s, math.exp(-2.0), math.exp(-2.0))


# ------------------------------------------------------------------ conversions (glue.py)
def mid_point(t, x1, x4, logt):
    ez, ew = e_vectors(t, logt)
    return x1 * ez[0] + x4 * ez[3], x1 * ew[0] + x4 * ew[3]


def to_cusp0_additive(t, x1, x4, logt):
    """Additive coordinate z_s of the cusp-0 slice point (zeta = e^{i z_s}) of the W-point (t, x1, x4)."""
    z, w = mid_point(t, x1, x4, logt)
    s = logt / (2j * math.pi)
    z, w = z - math.pi * s * (s - 1), w - 2 * t * s
    c2, c3 = real_coeffs(-w, CUSP0.e2(t), CUSP0.e3(t))
    return z + c2 * (-1j * t) + c3 * (-t)


def to_cusp0(t, x1, x4, logt):
    return torch.exp(1j * to_cusp0_additive(t, x1, x4, logt)), t


def to_cuspinf_additive(t, x1, x4, logt):
    z, w = mid_point(t, x1, x4, logt)
    s = 1 / t
    z2, w2 = z - 1j * w / t, w / t ** 2
    c2, c3 = real_coeffs(-w2, CUSPINF.e2(s), CUSPINF.e3(s))
    return z2 - 1j * (c2 * CUSPINF.al2(s) + c3 * CUSPINF.al3(s)), s


def to_cuspinf(t, x1, x4, logt):
    z2, s = to_cuspinf_additive(t, x1, x4, logt)
    return torch.exp(1j * z2), s


def E0_W(t, x1, x4, logt, sign=1.0):
    return CUSP0.E_torus(*to_cusp0(t, x1, x4, logt), sign=sign)


def Einf_W(t, x1, x4, logt, sign=1.0):
    return CUSPINF.E_torus(*to_cuspinf(t, x1, x4, logt), sign=sign)


# ------------------------------------------------------------------ Hopf identification (hopf.py)
RA, RB = 0.3, 0.5
NODE_IN, NODE_OUT = 0.02, 0.05
K0, KINF = 4 * math.pi, math.exp(2.0) / (4 * math.pi)
LK0, LK1 = -4.0, 4.0


def wrap(d):
    return d - torch.round(d).detach()


def xi4(t):
    logt = torch.log(t)
    sz, sw = sigma_tilde(t, logt)
    return lattice_coords(sz, sw, t, logt)[..., 3]


def xi23(t):
    logt = torch.log(t)
    sz, sw = sigma_tilde(t, logt)
    x = lattice_coords(sz, sw, t, logt)
    return x[..., 1], x[..., 2]


def k_of(r):
    pos = r > 0
    lr = torch.log(where_safe(pos, r, 1.0))
    lam = torch.where(pos, bump_step((lr - LK0) / (LK1 - LK0)), torch.zeros_like(r))
    return torch.exp((1 - lam) * math.log(K0) + lam * math.log(KINF))


def T_of(t):
    return t * k_of(t.abs())


def node_correct(C, mu, n, x0, y0):
    """mu_new: mu replaced near the node by n + a0|x0|^2 - |y0|^2 ((x0, y0) = U_0 coordinates of psi_hat^-n p)."""
    m0 = C.a0 * abs2(x0) - abs2(y0)
    rad = C.a0 * abs2(x0) + abs2(y0)
    c = 1 - bump_step((rad - NODE_IN) / (NODE_OUT - NODE_IN))
    return mu + c * (n + m0 - mu)


def mu_new_torus(C, zeta, t):
    mu = C.mu_torus(zeta, t)
    n, x0, y0 = C.reduce(zeta, t)
    return node_correct(C, mu, n, x0, y0)


def mu_new_chart(C, x, y):
    """Chart input (U_0, |mu| < 1/2, so n = 0)."""
    mu = C.mu_chart(x, y)
    return node_correct(C, mu, torch.zeros_like(mu), x, y)


def a1_c0(zeta, t):
    r = t.abs()
    mn = mu_new_torus(CUSP0, zeta, t)
    x4sl = torch.log(zeta.abs()) / CUSP0.P(t)
    lam = bump_step((r - RA) / (RB - RA))
    return mn + lam * wrap(x4sl - 0.5 - mn)


def a1_mid(t, x4):
    return x4 - xi4(t) - 0.5


def a2_ci(zp, s):
    r = s.abs()
    mn = mu_new_torus(CUSPINF, zp, s)
    x4sl = torch.log(zp.abs()) / CUSPINF.P(s)
    lam = bump_step((r - RA) / (RB - RA))
    return mn + lam * wrap(x4sl - 0.5 - mn)


def a2_mid(t, x4):
    return -x4 - 0.5


def delta(t):
    return torch.angle(t) / TWO_PI - xi4(t)


def lam_inf(r):
    return 1 - bump_step((r - 1.2) / (1.6 - 1.2))


def hopf_from_a1(t, a1):
    T = T_of(t)
    u2 = torch.exp(2j * math.pi * a1) / torch.sqrt(1 + abs2(T))
    return torch.stack([T * u2, u2], -1)


def hopf_from_a2(t, a2):
    lam = lam_inf(t.abs())
    a2p = -a2 + lam * delta(t) + (1 - lam) * 0.5
    T = T_of(t)
    u1 = torch.exp(2j * math.pi * a2p) / torch.sqrt(1 + (1 / abs2(T)))
    return torch.stack([u1, u1 / T], -1)


def hopf_from_a2_s(s, a2):
    """|t| >= 1.6 side written in s = 1/t (valid at s = 0, where lam_inf = 0)."""
    r_inv = s.abs()
    pos = r_inv > 0
    r = 1 / where_safe(pos, r_inv, 1.0)
    a2p = -a2 + 0.5                                     # lam_inf = 0 for |t| >= 1.6
    kk = torch.where(pos, k_of(r), torch.full_like(r, KINF))
    invT = s / kk
    u1 = torch.exp(2j * math.pi * a2p) / torch.sqrt(1 + abs2(invT))
    return torch.stack([u1, u1 * invT], -1)


def h_B_mid(t, x4):
    """Valid for 0.5 <= |t| <= 2."""
    lo = t.abs() <= 1
    out = torch.zeros(t.shape + (2,), dtype=CD)
    out = _scatter(out, lo, hopf_from_a1(t[lo], a1_mid(t[lo], x4[lo])))
    out = _scatter(out, ~lo, hopf_from_a2(t[~lo], a2_mid(t[~lo], x4[~lo])))
    return out


def h_B_c0(zeta, t):
    return hopf_from_a1(t, a1_c0(zeta, t))


def h_B_ci(zp, s):
    """For |s| <= 1/1.6 (where lam_inf = 0)."""
    return hopf_from_a2_s(s, a2_ci(zp, s))


# ------------------------------------------------------------------ orbit-space target maps (curve.py, lam.py)
def moebius(u):
    P = torch.stack([u[..., 0].real, u[..., 0].imag, u[..., 1].real, u[..., 1].imag], -1)
    den = (P[..., 1] ** 2 + P[..., 2] ** 2 + P[..., 3] ** 2) / (1 - P[..., 0])
    w = P[..., 1:] / den[..., None]
    return w - torch.tensor([0.0, 1.0, 0.0], dtype=RD)


RHO_SIGN = -1.0


def rho_O(w):
    return torch.stack([w[..., 0], RHO_SIGN * w[..., 2], -w[..., 1]], -1)


def nu_fn(tau):
    a = tau.abs()
    b = bump_step((a - 0.5) / 0.5)
    return (1 - b) * tau ** 2 + b * a


def P_final(w):
    rw = torch.linalg.norm(w, dim=-1)
    tau = torch.log(rw)
    h = torch.exp(-nu_fn(tau))
    # X = sign(tau) sqrt(1 - h); for |tau| < 1/2 (nu = tau^2) written as tau sqrt((1 - e^{-x}) / x), x = tau^2,
    # which is smooth through tau = 0 (series below x = 1e-8)
    small = tau.abs() < 0.5
    x = where_safe(small, tau, 1.0) ** 2
    tiny = x < 1e-8
    xs = where_safe(~tiny, x, 1.0)
    f = torch.where(tiny, 1 - x / 4 + 5 * x ** 2 / 96, torch.sqrt(-torch.expm1(-xs) / xs))
    X = torch.where(small, tau * f, torch.sign(tau) * torch.sqrt((1 - h).clamp_min(0.0)))
    v = h[..., None] * rho_O(w / rw[..., None])
    return v, X
