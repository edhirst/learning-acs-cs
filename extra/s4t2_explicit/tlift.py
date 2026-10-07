"""The S^1-equivariant diffeomorphism Psi_W: W -> S^4 and the diffeomorphism F = (Psi_W o pi, phi): Y -> S^4 x T^2,
in torch (float64, differentiable to any order). Design and verification: docs/s4t2_explicit_diffeo.md, section 12.

W-points are passed in the 'mid' representation (t, x1, x4, logt): lattice coordinates of the T^2-slice point
(x2 = x3 = 0) over t in C*, with a branch logt of log t (changing it by 2 pi i changes x1 by -x4 and nothing else).
Near the central fibres the chart representations are used (functions with suffix _chart): U_0-coordinates (x, y)
of the slice point of the filling, |mu| < 1/2.
"""
from __future__ import annotations

import json
import math
import os

import torch

from extra.s4t2_explicit import tcore as T
from extra.s4t2_explicit.tcore import CD, RD, _scatter, bump_step, where_safe

# ------------------------------------------------------------------ the weight-1 function E (section 12.1)
R0, RP, R1 = 0.6, 0.7, 0.8      # cusp formulas below R0 (above 1/R0); plateau radius RP; linear blend on [R1, 1/R1]


def h0(r):
    return r + (RP - r) * bump_step((r - R0) / (R1 - R0))


def hinf(r):
    return r + (1 / RP - r) * bump_step((1 / R0 - r) / (1 / R0 - 1 / R1))


def lam_mid(r):
    return bump_step((torch.log(r) - math.log(R1)) / (math.log(1 / R1) - math.log(R1)))


def _at_radius(rr, th):
    """t and its continuous log at radius rr and angle th (th = Im logt of the original point)."""
    return torch.polar(rr, th), torch.complex(torch.log(rr), th)


def E_W(t, x1, x4, logt):
    """E on W (mid representation):
         |t| <= R0:        cusp-0 formula E0_W;      R0 < |t| <= R1: E0_W read at the radius h0(|t|) (-> RP);
         R1 < |t| < 1/R1:  (1 - lam) E0_W(at RP) + lam Einf_W(at 1/RP)   (a linear combination: smooth);
         symmetric on the infinity side."""
    r, th = t.abs(), logt.imag
    out = torch.zeros(t.shape, dtype=CD)
    a, b = r <= R0, r >= 1 / R0
    m1 = (r > R0) & (r <= R1)
    m3 = (r >= 1 / R1) & (r < 1 / R0)
    m2 = (r > R1) & (r < 1 / R1)
    out = _scatter(out, a, T.E0_W(t[a], x1[a], x4[a], logt[a]))
    out = _scatter(out, b, T.Einf_W(t[b], x1[b], x4[b], logt[b]))
    if m1.any():
        tt, lt = _at_radius(h0(r[m1]), th[m1])
        out = _scatter(out, m1, T.E0_W(tt, x1[m1], x4[m1], lt))
    if m3.any():
        tt, lt = _at_radius(hinf(r[m3]), th[m3])
        out = _scatter(out, m3, T.Einf_W(tt, x1[m3], x4[m3], lt))
    if m2.any():
        lam = lam_mid(r[m2])
        t0, l0 = _at_radius(torch.full_like(r[m2], RP), th[m2])
        ti, li = _at_radius(torch.full_like(r[m2], 1 / RP), th[m2])
        out = _scatter(out, m2, (1 - lam) * T.E0_W(t0, x1[m2], x4[m2], l0) + lam * T.Einf_W(ti, x1[m2], x4[m2], li))
    return out


def E_chart0(x, y):
    """E near the central fibre at 0, U_0 coordinates (equals E_W through Phi; check_tlift.py)."""
    return T.CUSP0.E_chart(x, y)


def E_chartinf(x, y):
    return T.CUSPINF.E_chart(x, y)


# ------------------------------------------------------------------ Hopf identification on mid representation
def hB_chart0(x, y):
    """h_B for a chart point near the central fibre at 0 (|t| <= RA, |mu| < 1/2): a1 = mu_new exactly."""
    return T.hopf_from_a1(x * y, T.mu_new_chart(T.CUSP0, x, y))


def hB_chartinf(x, y):
    """h_B for a chart point near the central fibre at infinity (|s| <= RA)."""
    return T.hopf_from_a2_s(x * y, T.mu_new_chart(T.CUSPINF, x, y))


def hB_W(t, x1, x4, logt):
    """h_B: B -> S^3 for a W-point in the mid representation (dispatch c0 / mid / ci; check_hopf.py overlaps)."""
    r = t.abs()
    out = torch.zeros(t.shape + (2,), dtype=CD)
    a, b = r < 0.55, r > 1.8
    m = ~(a | b)
    if a.any():
        out = _scatter(out, a, T.h_B_c0(*T.to_cusp0(t[a], x1[a], x4[a], logt[a])))
    if b.any():
        out = _scatter(out, b, T.h_B_ci(*T.to_cuspinf(t[b], x1[b], x4[b], logt[b])))
    if m.any():
        out = _scatter(out, m, T.h_B_mid(t[m], x4[m]))
    return out


# ------------------------------------------------------------------ the orbit-space map Lambda (section 12.2)
def _x4_zero0(rho):
    """Mid lattice coordinate x4 of the zero of the cusp-0 formula on the torus |t| = rho (rho <= RP): |x| = |y|,
    arg t = pi (exact: there E = y + conj(x))."""
    t = torch.polar(rho, torch.full_like(rho, math.pi))
    return 0.5 * torch.log(rho) / T.CUSP0.P(t) + T.xi4(t)


def _x4_zeroinf(rho):
    """Same at infinity (rho >= 1/RP): |y'| = |nu| |x'|, arg s = -pi; +1 selects the representative continuing the
    cusp-0 branch through the middle (x4 = 0.7246 at 1/RP, versus 0.6655 at RP)."""
    t = torch.polar(rho, torch.full_like(rho, math.pi))
    s = 1 / t
    return -0.5 * (T.CUSPINF.lognu + torch.log(s.abs())) / T.CUSPINF.P(s) + 1.0


def x4_base(r):
    """1/2 (the exact zero |x| = |y| of the cusp-0 formula for r <= 0.366) blended into the exact zero of the
    cusp-infinity formula (valid for r >= 1.65); the blend sits on [1.25, 1.6], away from |t| = 1 where the
    cusp-infinity slice degenerates."""
    B = bump_step((torch.log(r) - math.log(1.25)) / (math.log(1.6) - math.log(1.25)))
    out = torch.full_like(r, 0.5)
    m = r > 1.2
    if m.any():
        out = _scatter(out, m, (1 - B[m]) * 0.5 + B[m] * _x4_zeroinf(r[m]))
    return out


ARC_INT1 = (math.log(0.78), math.log(1.28))     # fit interval of p1 (residual arg t - pi supported in [0.813, 1.232])
ARC_WIN1 = (math.log(0.80), math.log(1.25))     # window W1 = 1 on this interval
ARC_INT2 = (math.log(0.30), math.log(2.00))     # fit interval of p2 (residual x4 - x4_base supported in [0.366, 1.65])
ARC_WIN2 = (math.log(0.35), math.log(1.75))
_ARC_FIT = None


def _window(rho, inner, outer):
    return bump_step((rho - outer[0]) / (inner[0] - outer[0])) * bump_step((outer[1] - rho) / (outer[1] - inner[1]))


def cheb_eval(c, x):
    """Clenshaw: sum_k c_k T_k(x)."""
    b1 = torch.zeros_like(x)
    b2 = torch.zeros_like(x)
    for ck in reversed(c[1:]):
        b1, b2 = 2 * x * b1 - b2 + ck, b1
    return x * b1 - b2 + c[0]


def _arc_fit():
    global _ARC_FIT
    if _ARC_FIT is None:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'arc_fit.json')
        _ARC_FIT = json.load(open(path)) if os.path.exists(path) else {}
    return _ARC_FIT


def _fitted(key, rho, interval, window):
    fit = _arc_fit()
    out = torch.zeros_like(rho)
    if key not in fit:
        return out
    a, b = interval
    m = (rho > a) & (rho < b)
    if m.any():
        xs = 2 * (rho[m] - a) / (b - a) - 1
        out = _scatter(out, m, _window(rho[m], window, interval) * cheb_eval(fit[key]['coef'], xs))
    return out


def arc_model(r):
    """(arg t, x4) of the explicit model of the zero arc Gamma_E at radius r (fit_arc.py)."""
    rho = torch.log(r)
    return math.pi + _fitted('p1', rho, ARC_INT1, ARC_WIN1), x4_base(r) + _fitted('p2', rho, ARC_INT2, ARC_WIN2)


GAM_A, GAM_B = 0.005, 26.0          # the straightening R is the identity for |t| <= GAM_A and |t| >= GAM_B


def arc_phases(r):
    """(phi1, phi2) in U(1)^2: R(u) = (conj(phi1) u1, conj(phi2) u2) sends h_B of the model arc to the real great
    semicircle (0, 1) -> (-1, 0) (u1 <= 0 <= u2). R is a diffeomorphism of S^3 (r is a function of |u1|/|u2|);
    windowed so that R = id exactly near both fixed points (there the arc is exactly real anyway)."""
    out1 = torch.ones(r.shape, dtype=CD)
    out2 = torch.ones(r.shape, dtype=CD)
    m = (r > GAM_A) & (r < GAM_B)
    if m.any():
        rm = r[m]
        th, x4 = arc_model(rm)
        u = hB_W(torch.polar(rm, th), torch.zeros_like(rm), x4, torch.complex(torch.log(rm), th))
        lr = torch.log(rm)
        win = _window(lr, (math.log(0.01), math.log(13.0)), (math.log(GAM_A), math.log(GAM_B)))
        a1 = win * torch.angle(-u[..., 0])
        a2 = win * torch.angle(u[..., 1])
        out1 = _scatter(out1, m, torch.polar(torch.ones_like(rm), a1))
        out2 = _scatter(out2, m, torch.polar(torch.ones_like(rm), a2))
    return out1, out2


def Lam_from_u(u, r):
    """(v, X) = P_final(moebius(R_gamma u)) in B_S = {|v| + X^2 = 1}."""
    p1, p2 = arc_phases(r)
    return T.P_final(T.moebius(torch.stack([p1.conj() * u[..., 0], p2.conj() * u[..., 1]], -1)))


# linear models at the fixed points (section 12.2): Lambda = M v_q exactly for |v_q| <= EPS1, where
#   v_q = Hopf(q), q = (sqrt(a0) conj(x0), y0) in U_0 coordinates of the reduced point (a0 = |nu|^2);
#   M0 = 2 pi Q0, Minf = pi Qinf with Q0, Qinf signed permutation matrices (derived and verified in check_tlift.py)
EPS1, EPS2 = 0.004, 0.008
Q0 = torch.tensor([[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]])
QINF = torch.tensor([[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]])
S0, SINF = 2 * math.pi, math.pi


def v_q(C, x0, y0):
    q1, q2 = math.sqrt(C.a0) * x0.conj(), y0
    h = 2 * q1 * q2.conj()
    return torch.stack([h.real, h.imag, T.abs2(q1) - T.abs2(q2)], -1)


def lin_blend(v, X, vq, rad, M, scale, sgn):
    """Lambda_mod near a fixed point: v -> (1 - c) M vq + c v with c = 0 for rad <= EPS1, 1 for rad >= EPS2
    (rad = |vq|, smooth upstairs); X = sgn sqrt(1 - |v|). Unchanged outside rad < EPS2."""
    c = bump_step((rad - EPS1) / (EPS2 - EPS1))
    vl = vq @ M.T
    vm = (1 - c)[..., None] * vl + c[..., None] * v
    inner = rad <= EPS1
    nv = torch.linalg.norm(torch.where(inner[..., None], torch.ones_like(vm), vm), dim=-1)
    nv = torch.where(inner, scale * rad, nv)
    shell = rad < EPS2
    return (torch.where(shell[..., None], vm, v),
            torch.where(shell, sgn * torch.sqrt((1 - nv).clamp_min(0.0)), X))


M0, MINF = S0 * Q0, SINF * QINF


def node_coords_mid(t, x1, x4, logt):
    """Reduced U_0 coordinates at cusp 0 and at infinity (for the linear models), mid representation."""
    z0 = T.CUSP0.reduce(*T.to_cusp0(t, x1, x4, logt))
    zi = T.CUSPINF.reduce(*T.to_cuspinf(t, x1, x4, logt))
    return z0, zi


def Lam_W(t, x1, x4, logt):
    """Lambda_mod on the mid representation."""
    r = t.abs()
    v, X = Lam_from_u(hB_W(t, x1, x4, logt), r)
    a, b = r < 0.01, r > 25.0          # linear zones: rad >= 2 sqrt(a0)|t| <= EPS2 => |t| <= 0.004, |s| <= 0.0296
    if a.any():
        _, x0, y0 = T.CUSP0.reduce(*T.to_cusp0(t[a], x1[a], x4[a], logt[a]))
        rad = T.CUSP0.a0 * T.abs2(x0) + T.abs2(y0)
        va, Xa = lin_blend(v[a], X[a], v_q(T.CUSP0, x0, y0), rad, M0, S0, -1.0)
        v, X = _scatter(v, a, va), _scatter(X, a, Xa)
    if b.any():
        _, x0, y0 = T.CUSPINF.reduce(*T.to_cuspinf(t[b], x1[b], x4[b], logt[b]))
        rad = T.CUSPINF.a0 * T.abs2(x0) + T.abs2(y0)
        vb, Xb = lin_blend(v[b], X[b], v_q(T.CUSPINF, x0, y0), rad, MINF, SINF, 1.0)
        v, X = _scatter(v, b, vb), _scatter(X, b, Xb)
    return v, X


def _lam_chart(C, hb, x, y, invr, M, scale, sgn):
    t = x * y
    rad = C.a0 * T.abs2(x) + T.abs2(y)
    vq = v_q(C, x, y)
    v = torch.zeros(t.shape + (3,), dtype=RD)
    X = torch.full(t.shape, sgn, dtype=RD)
    gen = rad > EPS1                                   # Lam_from_u is singular exactly at the fixed point
    if gen.any():
        r = t[gen].abs()
        vg, Xg = Lam_from_u(hb(x[gen], y[gen]), 1 / r if invr else r)
        v, X = _scatter(v, gen, vg), _scatter(X, gen, Xg)
    return lin_blend(v, X, vq, rad, M, scale, sgn)


def Lam_chart0(x, y):
    return _lam_chart(T.CUSP0, hB_chart0, x, y, False, M0, S0, -1.0)


def Lam_chartinf(x, y):
    return _lam_chart(T.CUSPINF, hB_chartinf, x, y, True, MINF, SINF, 1.0)


# ------------------------------------------------------------------ the gauge G_A near Lambda^-1(M1) (section 12.3)
GA0, GAI = (0.4, 0.5), (2.0, 2.5)      # glue radii: cusp-0 local model -> e^{2 pi i x1} -> cusp-inf local model


def canonical_branch(x1, x4, logt):
    """Branch of log t with arg t in [0, 2 pi): logt - 2 pi i k, x1 + k x4 (k = floor(arg / 2 pi))."""
    k = torch.floor(logt.imag / T.TWO_PI).detach()
    return x1 + k * x4, logt - 2j * math.pi * k


def _log_local0(t, x1, x4, logt):
    """log(y0 - conj(x0)) at cusp 0 (reduced point), with consistent branches (z_s and log t as given)."""
    zs = T.to_cusp0_additive(t, x1, x4, logt)
    zeta = torch.exp(1j * zs)
    n = torch.round(T.CUSP0.mu_torus(zeta, t)).detach()
    c2, c3 = T.CUSP0.c_tilde(t)
    logmult = T.CUSP0.lognu + logt - (c2 * T.CUSP0.al2(t) + c3 * T.CUSP0.al3(t))
    logy0 = 1j * zs - n * logmult
    y0 = torch.exp(logy0)
    return logy0 + torch.log(1 - t.conj() / T.abs2(y0))


def _log_localinf(t, x1, x4, logt):
    """log(y0' - |nu| conj(x0')) at cusp infinity (reduced point), branches as given (log s = -log t)."""
    zs, s = T.to_cuspinf_additive(t, x1, x4, logt)
    zeta = torch.exp(1j * zs)
    n = torch.round(T.CUSPINF.mu_torus(zeta, s)).detach()
    c2, c3 = T.CUSPINF.c_tilde(s)
    logmult = T.CUSPINF.lognu - logt - (c2 * T.CUSPINF.al2(s) + c3 * T.CUSPINF.al3(s))
    logy0 = 1j * zs - n * logmult
    y0 = torch.exp(logy0)
    return logy0 + torch.log(1 - T.CUSPINF.e0coef * s.conj() / T.abs2(y0))


def GA_hat_W(t, x1, x4, logt):
    """Nonvanishing weight-1 function on the tube {theta_v <= THETA_B} (mid representation):
         = y0 - conj(x0) (cusp 0), e^{2 pi i x1} (canonical branch), y0' - |nu| conj(x0') (cusp infinity),
       glued by exp of a weighted log-ratio on |t| in GA0 and GAI."""
    x1c, lc = canonical_branch(x1, x4, logt)
    r = t.abs()
    expo = 2j * math.pi * x1c.to(CD)                       # log of the weight-1 factor, to be corrected below
    a = r < GA0[1]
    if a.any():
        beta = bump_step((r[a] - GA0[0]) / (GA0[1] - GA0[0]))
        L0 = 2j * math.pi * x1c[a] - _log_local0(t[a], x1c[a], x4[a], lc[a])
        expo = _scatter(expo, a, expo[a] - (1 - beta) * L0)
    b = r > GAI[0]
    if b.any():
        beta = bump_step((r[b] - GAI[0]) / (GAI[1] - GAI[0]))
        Li = 2j * math.pi * x1c[b] - _log_localinf(t[b], x1c[b], x4[b], lc[b])
        expo = _scatter(expo, b, expo[b] - beta * Li)
    return torch.exp(expo)


def GA_hat_chart0(x, y):
    return y - x.conj()


def GA_hat_chartinf(x, y):
    return y - T.CUSPINF.e0coef * x.conj()


# ------------------------------------------------------------------ the lift Psi_W (section 12.4)
THETA_A, THETA_B = math.radians(15.0), math.radians(45.0)
_PSI = None


PSI_R = (0.002, 60.0)       # psi = PSI_LO exactly for |t| <= 0.002 and PSI_HI exactly for |t| >= 60 (there the shell
PSI_LO = -math.pi / 2       # lies in the linear balls, where the two lifts must coincide: g = e^{i(beta - psi)} = 1)


def psi_of_r(r):
    """Branch normaliser: G_A = e^{i psi(|t|)} G_A_hat. psi = -pi/2 (= beta_0) near p0, PSI_HI (= beta_inf + 2 pi k)
    near p_inf; in between a fitted smooth function of log|t| (fit_psi.py) tracking the circular mean of the
    transition phase, so that |Arg g| < pi on the whole shell (verified in check_lift.py)."""
    global _PSI
    if _PSI is None:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'psi_fit.json')
        _PSI = json.load(open(path)) if os.path.exists(path) else {}
    if not _PSI:
        return torch.zeros_like(r)
    a, b = math.log(PSI_R[0]), math.log(PSI_R[1])
    pos = r > 0
    rho = torch.log(where_safe(pos, r, 1.0))
    rho = torch.where(pos, rho, torch.full_like(rho, a - 1.0))
    step = bump_step((rho - a) / (b - a))
    out = PSI_LO + (_PSI['psi_hi'] - PSI_LO) * step
    m = (rho > a) & (rho < b)
    if m.any():
        xs = 2 * (rho[m] - a) / (b - a) - 1
        win = bump_step((rho[m] - a) / _PSI['edge']) * bump_step((b - rho[m]) / _PSI['edge'])
        out = _scatter(out, m, out[m] + win * cheb_eval(_PSI['coef'], xs))
    return out


def shell_weight(v):
    """c = 0 on {theta_v <= THETA_A} (lift A), 1 on {theta_v >= THETA_B} (lift B); theta_v = angle(v, (0,0,-1))."""
    s = -v[..., 2] / torch.linalg.norm(v, dim=-1)
    return bump_step((math.cos(THETA_A) - s) / (math.cos(THETA_A) - math.cos(THETA_B)))


def assemble(v, X, E, GAh, r, use_psi=True):
    """Psi = (Z1, Z2, X) from Lambda = (v, X), E (lift B) and G_A = e^{i psi(r)} GAh (lift A), glued on the shell;
    r = |t| of the point."""
    c = shell_weight(v)
    nv = torch.linalg.norm(v, dim=-1)
    w = torch.complex(v[..., 0], v[..., 1])
    A, B = c < 1, c > 0
    Z1 = torch.zeros(X.shape, dtype=CD)
    Z2 = torch.zeros(X.shape, dtype=CD)
    if A.any():
        GA = GAh[A] * (torch.polar(torch.ones_like(X[A]), psi_of_r(r[A])) if use_psi else 1.0)
        om = GA / GA.abs()
        dm = torch.sqrt(2 * (nv[A] - v[A][..., 2]))
        z1a, z2a = w[A] * om / dm, dm / 2 * om
        sh = B[A]
        if sh.any():                                      # shell: multiply lift A by g^c, g = lift B / lift A
            ob = E[A][sh] / E[A][sh].abs()
            g = ob * om[sh].conj() * (w[A][sh] / w[A][sh].abs()).conj()
            ph = torch.polar(torch.ones_like(c[A][sh]), c[A][sh] * torch.angle(g))
            z1a, z2a = _scatter(z1a, sh, z1a[sh] * ph), _scatter(z2a, sh, z2a[sh] * ph)
        Z1, Z2 = _scatter(Z1, A, z1a), _scatter(Z2, A, z2a)
    onlyB = ~A
    if onlyB.any():
        ob = E[onlyB] / E[onlyB].abs()
        dp = torch.sqrt(2 * (nv[onlyB] + v[onlyB][..., 2]))
        Z1, Z2 = _scatter(Z1, onlyB, dp / 2 * ob), _scatter(Z2, onlyB, w[onlyB].conj() * ob / dp)
    return Z1, Z2, X


def g_shell(v, X, E, GAh, r, use_psi=True):
    """The transition phase g = lift B / lift A (diagnostic)."""
    GA = GAh * (torch.polar(torch.ones_like(X), psi_of_r(r)) if use_psi else 1.0)
    w = torch.complex(v[..., 0], v[..., 1])
    return (E / E.abs()) * (GA / GA.abs()).conj() * (w / w.abs()).conj()


# exact linear models (section 12.2): Psi = sqrt(s) U q, q = (sqrt(a0) conj(x0), y0), X = sgn sqrt(1 - s |q|^2)
SQ2 = math.sqrt(2.0)


def linear_model0(x, y):
    q1, q2 = x.conj(), y
    k = math.sqrt(S0)
    return k * (q1 + q2) / SQ2, k * 1j * (q1 - q2) / SQ2, -torch.sqrt(1 - S0 * (T.abs2(q1) + T.abs2(q2)))


def linear_modelinf(x, y):
    q1, q2 = T.CUSPINF.e0coef * x.conj(), y
    k = math.sqrt(SINF)
    return k * (q1 + q2) / SQ2, k * (q1 - q2) / SQ2, torch.sqrt(1 - SINF * (T.abs2(q1) + T.abs2(q2)))


def _psi_chart(C, lam, Efn, Gfn, lin, x, y, invr, use_psi=True):
    rad = C.a0 * T.abs2(x) + T.abs2(y)
    inner = rad <= EPS1
    Z1 = torch.zeros(x.shape, dtype=CD)
    Z2 = torch.zeros(x.shape, dtype=CD)
    X = torch.zeros(x.shape, dtype=RD)
    if inner.any():
        a, b, c = lin(x[inner], y[inner])
        Z1, Z2, X = _scatter(Z1, inner, a), _scatter(Z2, inner, b), _scatter(X, inner, c)
    o = ~inner
    if o.any():
        v, Xo = lam(x[o], y[o])
        rr = (x[o] * y[o]).abs()
        a, b, c = assemble(v, Xo, Efn(x[o], y[o]), Gfn(x[o], y[o]), 1 / rr if invr else rr, use_psi)
        Z1, Z2, X = _scatter(Z1, o, a), _scatter(Z2, o, b), _scatter(X, o, c)
    return Z1, Z2, X


def Psi_chart0(x, y, use_psi=True):
    """Psi_W near the central fibre at 0 (U_0 chart, |t| <= 0.25, |mu| < 1/2)."""
    return _psi_chart(T.CUSP0, Lam_chart0, E_chart0, GA_hat_chart0, linear_model0, x, y, False, use_psi)


def Psi_chartinf(x, y, use_psi=True):
    return _psi_chart(T.CUSPINF, Lam_chartinf, E_chartinf, GA_hat_chartinf, linear_modelinf, x, y, True, use_psi)


def Psi_W(t, x1, x4, logt, use_psi=True):
    """Psi_W on the mid representation (any t != 0; points inside the linear balls are sent to the chart formulas)."""
    r = t.abs()
    Z1 = torch.zeros(t.shape, dtype=CD)
    Z2 = torch.zeros(t.shape, dtype=CD)
    X = torch.zeros(t.shape, dtype=RD)
    done = torch.zeros(t.shape, dtype=torch.bool)
    a, b = r < 0.01, r > 25.0
    if a.any():
        _, x0, y0 = T.CUSP0.reduce(*T.to_cusp0(t[a], x1[a], x4[a], logt[a]))
        inn = (T.abs2(x0) + T.abs2(y0)) <= EPS1
        if inn.any():
            idx = a.nonzero().flatten()[inn]
            m = torch.zeros_like(done); m[idx] = True
            p, q, s = linear_model0(x0[inn], y0[inn])
            Z1, Z2, X, done = _scatter(Z1, m, p), _scatter(Z2, m, q), _scatter(X, m, s), done | m
    if b.any():
        _, x0, y0 = T.CUSPINF.reduce(*T.to_cuspinf(t[b], x1[b], x4[b], logt[b]))
        inn = (T.CUSPINF.a0 * T.abs2(x0) + T.abs2(y0)) <= EPS1
        if inn.any():
            idx = b.nonzero().flatten()[inn]
            m = torch.zeros_like(done); m[idx] = True
            p, q, s = linear_modelinf(x0[inn], y0[inn])
            Z1, Z2, X, done = _scatter(Z1, m, p), _scatter(Z2, m, q), _scatter(X, m, s), done | m
    o = ~done
    if o.any():
        v, Xo = Lam_W(t[o], x1[o], x4[o], logt[o])
        p, q, s = assemble(v, Xo, E_W(t[o], x1[o], x4[o], logt[o]), GA_hat_W(t[o], x1[o], x4[o], logt[o]),
                           r[o], use_psi)
        Z1, Z2, X = _scatter(Z1, o, p), _scatter(Z2, o, q), _scatter(X, o, s)
    return Z1, Z2, X
