"""F^-1: S^4 x T^2 -> Y (docs/s4t2_explicit_diffeo.md, section 12.7), so that J_round can be evaluated at ANY point
of the round S^4 x T^2 (e.g. uniform samples, supervised targets).

Explicit chain (inverting Lambda = P_final o moebius o R o h_B step by step), then a Newton polish on F itself:
  1. v = Hopf(Z1, Z2); w = P_final^-1(v, X)  (radial: invert nu, monotone; direction rho_O^-1);
  2. u' = inverse stereographic projection from (-1, 0) of w + (0, 1, 0);
  3. r from |u'1| / |u'2| = |T(r)| = r k(r)  (monotone, Newton in log r); u = R^-1 u' (R depends on r only);
  4. t = r e^{i arg(u1/u2)}; fibre coordinate a1 = arg(u2)/2 pi (r <= 1) or a2' = arg(u1)/2 pi (r >= 1);
  5. W-point: middle  x4 = a1 + xi4(t) + 1/2,  resp. a2 = -a2' + lam delta + (1 - lam)/2, x4 = -a2 - 1/2;
              cusps   log|zeta| from a1_c0 (resp. a2_ci) by 1D Newton (monotone along the fibre);
     inside the linear balls the lift is linear: q = U^-1 Z / sqrt(s) directly;
  6. S^1 phase from Z / Psi(point with x1 = 0); T^2 part: (x2, x3) = phi + corrections (resp. C = phi + H c);
  7. Newton polish p <- p + DF(p)^+ (q - F(p)) in the output chart (mid for 0.2 <= |t| <= 5, else the fillings).
"""
from __future__ import annotations

import math

import torch

from extra.s4t2_explicit import tcore as T, tlift as L, tfull as Fm
from extra.s4t2_explicit.tcore import CD, RD, bump_step
from torch.func import jvp

EYE6 = torch.eye(6)


def nu_inv(nu):
    """|tau| with nu_fn(|tau|) = nu (nu_fn is increasing in |tau|): sqrt for nu <= 1/4, identity for nu >= 1,
    Newton on [1/2, 1] in between."""
    a = torch.where(nu <= 0.25, torch.sqrt(nu.clamp_min(0)), torch.where(nu >= 1, nu, 0.75 * torch.ones_like(nu)))
    mid = (nu > 0.25) & (nu < 1)
    for _ in range(60):
        aa = a.clone().requires_grad_(True)
        f = T.nu_fn(aa) - nu
        g, = torch.autograd.grad(f.sum(), aa)
        a = torch.where(mid, (a - f.detach() / g).clamp(0.5, 1.0), a)
    return a


def P_final_inv(v, X):
    nv = torch.linalg.norm(v, dim=-1)
    tau = torch.sign(X) * nu_inv(-torch.log(nv))
    d = torch.stack([v[..., 0], -v[..., 2], -v[..., 1]], -1) / nv[..., None]      # rho_O^-1 (RHO_SIGN = -1)
    return torch.exp(tau)[..., None] * d


def moebius_inv(w):
    y = w + torch.tensor([0.0, 1.0, 0.0], dtype=RD)
    n2 = (y ** 2).sum(-1)
    P1 = (1 - n2) / (1 + n2)
    Pr = 2 * y / (1 + n2)[..., None]
    return torch.stack([torch.complex(P1, Pr[..., 0]), torch.complex(Pr[..., 1], Pr[..., 2])], -1)


def r_from_absT(aT):
    """Solve |T(r)| = r k(r) = aT (increasing) by Newton in log r."""
    lr = torch.log(aT / T.K0)
    for _ in range(80):
        x = lr.clone().requires_grad_(True)
        f = torch.log(T.k_of(torch.exp(x))) + x - torch.log(aT)
        g, = torch.autograd.grad(f.sum(), x)
        step = (f / g).detach()
        lr = lr - step.clamp(-1.0, 1.0)
        if step.abs().max() < 1e-15:
            break
    return torch.exp(lr)


def _cusp_logzeta(afun, C, a_target, t):
    """log|zeta| with afun(e^l, t) = a_target mod 1 (afun decreasing in l: x4_slice = l / P with P < 0)."""
    l = (a_target + 0.5) * C.P(t)
    for _ in range(80):
        x = l.clone().requires_grad_(True)
        f = T.wrap(afun(torch.exp(x).to(CD), t) - a_target)
        g, = torch.autograd.grad(f.sum(), x)
        step = (f / g).detach()
        l = l - step.clamp(-0.5 * C.P(t).abs(), 0.5 * C.P(t).abs())
        if step.abs().max() < 1e-14:
            break
    return l


def W_from_orbit(v, X):
    """A point of W (mid representation, x1 = 0) with Lambda(point) = (v, X), via the explicit chain (exact outside the
    linear-blend shells; there it is an initial guess for the Newton polish)."""
    u1 = moebius_inv(P_final_inv(v, X))
    r = r_from_absT(u1[..., 0].abs() / u1[..., 1].abs())
    p1, p2 = L.arc_phases(r)
    u = torch.stack([p1 * u1[..., 0], p2 * u1[..., 1]], -1)
    t = torch.polar(r, torch.angle(u[..., 0] / u[..., 1]))
    logt = torch.log(t)
    x4 = torch.zeros_like(r)
    a1 = torch.angle(u[..., 1]) / T.TWO_PI
    a2p = torch.angle(u[..., 0]) / T.TWO_PI
    m = (r >= 0.5) & (r <= 1)
    x4 = torch.where(m, a1 + T.xi4(t) + 0.5, x4)
    m = (r > 1) & (r <= 2)
    lam = T.lam_inf(r)
    a2 = -a2p + lam * T.delta(t) + (1 - lam) * 0.5
    x4 = torch.where(m, -a2 - 0.5, x4)
    zeta = torch.ones_like(t)
    c0, ci = r < 0.5, r > 2
    if c0.any():
        l = _cusp_logzeta(T.a1_c0, T.CUSP0, a1[c0], t[c0])
        zeta = zeta.clone(); zeta[c0] = torch.exp(l).to(CD)
    if ci.any():
        s = 1 / t[ci]
        l = _cusp_logzeta(T.a2_ci, T.CUSPINF, -a2p[ci] + 0.5, s)        # lam_inf = 0 there: a2 = -a2' + 1/2
        zeta = zeta.clone(); zeta[ci] = torch.exp(l).to(CD)
    return t, logt, x4, zeta, c0, ci


def _lattice_from_cusp0(zeta, t):
    """Mid lattice coordinates (x1, x4) of the cusp-0 slice point (zeta, t), via Phi."""
    logt = torch.log(t)
    s = logt / (2j * math.pi)
    x = T.lattice_coords(-1j * torch.log(zeta) + math.pi * s * (s - 1), 2 * t * s, t, logt)
    return x[..., 0], x[..., 3], logt


def _lattice_from_cuspinf(zeta, s):
    t = 1 / s
    logt = -torch.log(s)
    x = T.lattice_coords(-1j * torch.log(zeta), torch.zeros_like(zeta), t, logt)
    return x[..., 0], x[..., 3], logt


def initial_guess(q):
    """Point of Y in mid lattice form (t, logt, x1..x4) with F approximately q (exact outside the blend shells)."""
    Z1, Z2 = torch.complex(q[..., 0], q[..., 1]), torch.complex(q[..., 2], q[..., 3])
    X = q[..., 4]
    h = 2 * Z1 * Z2.conj()
    v = torch.stack([h.real, h.imag, Z1.abs() ** 2 - Z2.abs() ** 2], -1)
    t, logt, x4, zeta, c0, ci = W_from_orbit(v, X)
    x1 = torch.zeros_like(x4)
    if c0.any():
        a, b, lg = _lattice_from_cusp0(zeta[c0], t[c0])
        x1[c0], x4[c0], logt = a, b, logt.clone().index_put_((c0.nonzero().flatten(),), lg)
    if ci.any():
        a, b, lg = _lattice_from_cuspinf(zeta[ci], 1 / t[ci])
        x1[ci], x4[ci], logt = a, b, logt.clone().index_put_((ci.nonzero().flatten(),), lg)
    W1, W2, _ = L.Psi_W(t, x1, x4, logt)
    th = torch.angle(Z1 * W1.conj() + Z2 * W2.conj())
    x1 = x1 + th / T.TWO_PI
    phi = torch.stack([torch.atan2(q[..., 6], q[..., 5]), torch.atan2(q[..., 8], q[..., 7])], -1) / T.TWO_PI
    x = torch.stack([x1, torch.zeros_like(x1), torch.zeros_like(x1), x4], -1)
    corr = x[..., 1:3] - Fm.phi_mid(t, x, logt)                         # phi = (x2, x3) - corr
    x23 = phi + corr
    return t, logt, torch.stack([x1, x23[..., 0], x23[..., 1], x4], -1)


def _filling_point(C, zz, ww, tt, logtt):
    """Filling coordinates (zeta = e^{i zz}, w, t) -> U_0 chart point (x, y, w) whose T^2-slice lies in the fundamental
    slab |mu| <= 1/2, using the cusp generator G4 = (nu t zeta, w + e4bar(t)) (its powers act as psi_hat^n)."""
    zeta = torch.exp(1j * zz)
    C2, C3 = T.real_coeffs(ww, C.e2(tt), C.e3(tt))
    zs = zeta * torch.exp(-(C2 * C.al2(tt) + C3 * C.al3(tt)))
    n = torch.round(C.mu_torus(zs, tt)).detach()
    zeta = zeta * torch.exp(-n * (C.lognu + logtt))
    ww = ww - n * C.e4(tt)
    return torch.stack([(tt / zeta).real, (tt / zeta).imag, zeta.real, zeta.imag, ww.real, ww.imag], -1)


def to_chart(t, logt, x):
    """The Y-point with mid lattice coordinates x, in the output charts: {'mid': (mask, p6), 'c0': ..., 'ci': ...};
    mid for 0.2 <= |t| <= 5, the fillings below / above."""
    z, w = T.point_from_lattice(t, logt, x)
    r = t.abs()
    s = logt / (2j * math.pi)
    return {'mid': ((r >= 0.2) & (r <= 5), torch.stack([z.real, z.imag, w.real, w.imag, t.real, t.imag], -1)),
            'c0': (r < 0.2, _filling_point(T.CUSP0, z - math.pi * s * (s - 1), w - 2 * t * s, t, logt)),
            'ci': (r > 5, _filling_point(T.CUSPINF, z - 1j * w / t, w / t ** 2, 1 / t, -logt))}


CHART_F = {'mid': Fm.F_mid_real, 'c0': Fm.F_c0_real, 'ci': Fm.F_ci_real}


def Finv(q, iters=8):
    """F^-1 for points q of S^4 x T^2 in R^9: returns {chart: (indices, p6, residual |F(p) - q|_inf)}."""
    t, logt, x = initial_guess(q)
    out = {}
    for key, val in to_chart(t, logt, x).items():
        mask, P = val
        if mask.any():
            p, res = newton_polish(CHART_F[key], P[mask], q[mask], iters)
            out[key] = (mask.nonzero().flatten(), p, res)
    return out


def newton_polish(f, p, q, iters=8):
    for _ in range(iters):
        Fp = f(p)
        A = torch.stack([jvp(f, (p,), (EYE6[k].expand_as(p),))[1] for k in range(6)], -1)
        d = q - Fp
        d[..., 5:9] = q[..., 5:9] - Fp[..., 5:9]
        p = p + torch.einsum('bkc,bc->bk', torch.linalg.pinv(A), d)
    return p, (f(p) - q).abs().max(-1).values
