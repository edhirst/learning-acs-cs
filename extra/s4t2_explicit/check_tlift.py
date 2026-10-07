"""Checks of tlift.py: well-definedness on W (branch of log t, lattice translation x4 -> x4 + 1, S^1 weights),
agreement of the chart and mid representations, the exact linear models at the fixed points.
Run from the repo root: python3 extra/s4t2_explicit/check_tlift.py"""
from __future__ import annotations
import math, os, sys
import torch
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
torch.set_default_dtype(torch.float64)
from extra.s4t2_explicit import tcore as T, tlift as L

g = torch.Generator().manual_seed(21)
U = lambda n, a, b: torch.empty(n).uniform_(a, b, generator=g)
fails = []


def rep(name, err, tol):
    ok = err < tol
    if not ok:
        fails.append(name)
    print(f'{name:<62s} {err:.1e}  {"ok" if ok else "FAIL"}')


def rel(a, b):
    return ((a - b).abs() / (1 + b.abs())).max().item()


def chart0_to_mid(x, y):
    """U_0 slice point (zeta = y, t = xy, w = 0) of the cusp-0 filling -> its W-point Phi(q) in mid lattice coords."""
    t = x * y
    logt = torch.log(t)
    s = logt / (2j * math.pi)
    z = -1j * torch.log(y) + math.pi * s * (s - 1)
    w = 2 * t * s
    xl = T.lattice_coords(z, w, t, logt)
    return t, xl[..., 0], xl[..., 3], logt


def chartinf_to_mid(x, y):
    s = x * y
    t = 1 / s
    logt = -torch.log(s)
    z = -1j * torch.log(y)
    xl = T.lattice_coords(z, torch.zeros_like(z), t, logt)
    return t, xl[..., 0], xl[..., 3], logt


n = 3000
# ---------------------------------------------------------------- 1. well-definedness on W (mid representation)
r = torch.exp(U(n, math.log(1e-3), math.log(1e3)))
th = U(n, math.pi / 2 + 0.05, 3 * math.pi / 2 - 0.05)          # the tube side (G_A is only used there)
t = torch.polar(r, th); lg = torch.complex(torch.log(r), th)
x1, x4 = U(n, 0, 1), U(n, 0, 1)
lg2, x12 = lg + 2j * math.pi, x1 - x4                            # another branch of log t
for name, f in (('E_W', L.E_W), ('GA_hat_W', L.GA_hat_W)):
    f0 = f(t, x1, x4, lg)
    rep(f'1. {name}: branch change of log t', rel(f(t, x12, x4, lg2), f0), 1e-11)
    rep(f'   {name}: x4 -> x4 + 1 (lattice translation e4)', rel(f(t, x1, x4 + 1, lg), f0), 1e-11)
    rep(f'   {name}: weight 1 under S^1_a (x1 -> x1 + 0.3)', rel(f(t, x1 + 0.3, x4, lg), f0 * math.e ** 0 * torch.polar(torch.ones(n), torch.full((n,), 2 * math.pi * 0.3))), 1e-11)
v0, X0 = L.Lam_W(t, x1, x4, lg)
for name, args in (('branch change', (t, x12, x4, lg2)), ('x4 -> x4 + 1', (t, x1, x4 + 1, lg)), ('S^1_a (x1 + 0.3)', (t, x1 + 0.3, x4, lg))):
    v1, X1 = L.Lam_W(*args)
    rep(f'   Lambda_W invariant under {name}', max(rel(v1, v0), rel(X1, X0)), 1e-11)
rep('   Lambda_W lands in B_S = {|v| + X^2 = 1}', (v0.norm(dim=-1) + X0 ** 2 - 1).abs().max().item(), 1e-13)
# ---------------------------------------------------------------- 2. chart vs mid representations
for nm, toMid, Echart, Lchart, Gchart, C in (('c0', chart0_to_mid, L.E_chart0, L.Lam_chart0, L.GA_hat_chart0, T.CUSP0),
                                             ('ci', chartinf_to_mid, L.E_chartinf, L.Lam_chartinf, L.GA_hat_chartinf, T.CUSPINF)):
    m = 3000
    rr = torch.exp(U(m, math.log(1e-3), math.log(0.25)))        # |t| (resp. |s|) of the chart point
    al = U(m, math.pi / 2 + 0.05, 3 * math.pi / 2 - 0.05)        # arg t (resp. arg s): tube side
    ratio = torch.exp(U(m, -1.5, 1.5))                            # |y| / |x|
    bx = U(m, -math.pi, math.pi)
    x = torch.polar(torch.sqrt(rr / ratio), bx); y = torch.polar(torch.sqrt(rr * ratio), al - bx)
    keep = C.mu_chart(x, y).abs() < 0.45
    x, y = x[keep], y[keep]
    tm, x1m, x4m, lgm = toMid(x, y)
    rep(f'2. {nm}: E chart = E mid', rel(Echart(x, y), L.E_W(tm, x1m, x4m, lgm)), 1e-10)
    va, Xa = Lchart(x, y); vb, Xb = L.Lam_W(tm, x1m, x4m, lgm)
    rep(f'   {nm}: Lambda chart = Lambda mid', max(rel(va, vb), rel(Xa, Xb)), 1e-10)
    rep(f'   {nm}: G_A_hat chart = G_A_hat mid (|t| <= 0.25: local model)', rel(Gchart(x, y), L.GA_hat_W(tm, x1m, x4m, lgm)), 1e-10)
# ---------------------------------------------------------------- 3. linear models at the fixed points
for nm, C, Lchart, lin in (('p0  ', T.CUSP0, L.Lam_chart0, L.linear_model0), ('pinf', T.CUSPINF, L.Lam_chartinf, L.linear_modelinf)):
    m = 2000
    x = 0.04 * torch.complex(torch.randn(m, generator=g), torch.randn(m, generator=g))
    y = 0.04 * torch.complex(torch.randn(m, generator=g), torch.randn(m, generator=g))
    rad = C.a0 * x.abs() ** 2 + y.abs() ** 2
    keep = (rad > 1e-6) & (rad < L.EPS1)
    x, y = x[keep], y[keep]
    v, X = Lchart(x, y)
    Z1, Z2, Xl = lin(x, y)
    vl = torch.stack([(2 * Z1 * Z2.conj()).real, (2 * Z1 * Z2.conj()).imag, Z1.abs() ** 2 - Z2.abs() ** 2], -1)
    rep(f'3. {nm}: Hopf(linear model) = Lambda_mod (inner ball)', max(rel(vl, v), rel(Xl, X)), 1e-12)
    Ech = (L.E_chart0 if nm == 'p0  ' else L.E_chartinf)(x, y)
    rep(f'   {nm}: Z1 of the linear model has the phase of E', ((Z1 / Z1.abs()) - Ech / Ech.abs()).abs().max().item(), 1e-12)
print('ALL PASS' if not fails else f'FAILURES: {fails}')
