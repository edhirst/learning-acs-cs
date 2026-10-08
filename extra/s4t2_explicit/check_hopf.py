"""Checks of the Hopf identification h_B (hopf.py)."""
from __future__ import annotations
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from extra.s4t2_explicit import hopf as H
from extra.s4t2_explicit.cusp import CUSP0, CUSPINF

rng = np.random.default_rng(8)
n = 3000
d = lambda A, B: np.abs(A - B).max()
# 1. delta single-valued and continuous across the negative real axis
r = np.exp(rng.uniform(np.log(0.5), np.log(2.0), 500))
print(f'1. delta across arg t = +-pi: max jump {np.abs(H.delta(r * np.exp(1j * (np.pi - 1e-9))) - H.delta(r * np.exp(-1j * (np.pi - 1e-9)))).max():.1e}')
# 2. overlaps of representations
r = np.exp(rng.uniform(np.log(0.5), np.log(0.6), n)); t = r * np.exp(1j * rng.uniform(-np.pi, np.pi, n))
zeta = np.exp(rng.uniform(-3, 3, n) + 1j * rng.uniform(-np.pi, np.pi, n))
x4 = np.log(np.abs(zeta)) / CUSP0.P(t) + H.xi4(t)
print(f'2. c0 vs mid on 0.5 <= |t| < 0.6: {d(H.h_B_c0(zeta, t), H.h_B_mid(t, x4)):.1e}')
s = 1 / (np.exp(rng.uniform(np.log(1 / 0.6), np.log(2.0), n)) * np.exp(1j * rng.uniform(-np.pi, np.pi, n)))
zp = np.exp(rng.uniform(-3, 3, n) + 1j * rng.uniform(-np.pi, np.pi, n))
x4b = -np.log(np.abs(zp)) / CUSPINF.P(s)
print(f'   ci vs mid on 0.5 <= |s| < 0.6: {d(H.h_B_ci(zp, s), H.h_B_mid(1 / s, x4b)):.1e}')
# 3. the two Hopf formulas agree at |t| = 1
t1 = np.exp(1j * rng.uniform(-np.pi, np.pi, n)); x41 = rng.uniform(0, 1, n)
A = H.hopf_from_a1(t1, H.a1_mid(t1, x41)); B = H.hopf_from_a2(t1, H.a2_mid(t1, x41))
print(f'3. Hopf charts agree on |t| = 1 (clutching): {d(A, B):.1e}; |u| = 1: {np.abs(np.linalg.norm(A, axis=-1) - 1).max():.1e}')
# 4. fibrewise monotonicity (bijectivity): a1 increasing along every fibre (in x4), incl. cusp zones
def fib_deriv(f, *args, h=1e-6):
    return (f(*args, h) - f(*args, -h)) / (2 * h)
r = np.exp(rng.uniform(np.log(1e-8), np.log(0.6), n)); t = r * np.exp(1j * rng.uniform(-np.pi, np.pi, n))
lz = rng.uniform(-3, 3, n) * np.abs(np.log(r))
g = lambda tt, l, h: H.a1_c0(np.exp(l + h * CUSP0.P(tt)) + 0j, tt)   # log|zeta| += h P  <=>  x4_slice += h
der = fib_deriv(g, t, lz)
print(f'4. d a1 / d x4 along cusp-0 fibres in [{der.min():.3e}, {der.max():.3e}] (> 0)')
gi = lambda ss, l, h: H.a2_ci(np.exp(l + h * CUSPINF.P(ss)) + 0j, ss)
deri = fib_deriv(gi, t, lz)              # s = t: small |s|
print(f'   d a2 / d x4_inf along cusp-inf fibres in [{deri.min():.3e}, {deri.max():.3e}] (> 0)')
# 5. fixed points and conformality (v-coordinates; chart U_0 at each cusp: zeta = y, t = x y)
for nm, C, f in (('b0  ', CUSP0, lambda z, tt: H.h_B_c0(z, tt)), ('binf', CUSPINF, lambda z, ss: H.h_B_ci(z, ss))):
    a0 = abs(C.nu) ** 2
    for eps in (1e-3, 1e-4):
        x = eps * (rng.normal(size=400) + 1j * rng.normal(size=400)); y = eps * (rng.normal(size=400) + 1j * rng.normal(size=400))
        u = f(y, x * y)
        if nm == 'b0  ':
            loc = np.stack([u[:, 0].real, u[:, 0].imag, np.angle(u[:, 1])], 1)
            v = np.stack([2 * (x * y).real, -2 * (x * y).imag, np.abs(x) ** 2 - np.abs(y) ** 2], 1)        # (2 t-bar, m)
        else:
            loc = np.stack([u[:, 1].real, u[:, 1].imag, np.angle(-u[:, 0])], 1)   # b_inf -> (-1, 0) since the 1/2 shift
            sb = np.conj(x * y)
            v = np.stack([2 * np.sqrt(a0) * sb.real, 2 * np.sqrt(a0) * sb.imag, a0 * np.abs(x) ** 2 - np.abs(y) ** 2], 1)
        L, *_ = np.linalg.lstsq(v, loc, rcond=None)
        sv = np.linalg.svd(L, compute_uv=False)
        res = np.abs(v @ L - loc).max() / np.abs(loc).max()
        print(f'5. {nm} eps = {eps:.0e}: |u - u_fixed| max {np.abs(u - (np.array([0, 1]) if nm == "b0  " else np.array([-1, 0]))).max():.1e}; '
              f'linear fit in v: singular values {np.round(sv, 5)}, det sign {np.sign(np.linalg.det(L)):+.0f}, nonlinear residual {res:.1e}')

# 6. T7: compare the simplified formula with the former wrapped formula.
# Include translated half-cell boundaries, both blend endpoints, and their neighbours.
# This tests the finite-series code; the infinite-series proof is C7 in the programme.
rr = np.concatenate([np.linspace(0.3, 0.6, 61), [0.3 - 1e-9, 0.3 + 1e-9, 0.5 - 1e-9, 0.5 + 1e-9]])
cc = np.concatenate([np.linspace(-4, 4, 161), [-0.5 - 1e-10, -0.5 + 1e-10, 0.5 - 1e-10, 0.5 + 1e-10]])
rr, cc = np.meshgrid(rr, cc, indexing='ij')
tt = rr.flatten() * np.exp(0.7j)
centre = cc.flatten()
for name, C, afun in (('c0', CUSP0, H.a1_c0), ('ci', CUSPINF, H.a2_ci)):
    zz = np.exp((centre + 0.5) * C.P(tt) + 0.2j)
    mn = H.mu_new(C, zz, tt)
    diff = np.log(np.abs(zz)) / C.P(tt) - 0.5 - mn
    lam = H.bump_step((np.abs(tt) - H.RA) / (H.RB - H.RA))
    former = mn + lam * H.wrap(diff)
    current = afun(zz, tt)
    err = d(former, current)
    active = lam > 0
    assert np.max(np.abs(diff[active])) < 0.5
    assert err < 1e-13
    shifted = afun(zz * C.multiplier(tt), tt)
    equiv_error = d(shifted, current + 1)
    assert equiv_error < 1e-12
    print(f'6. {name}: {len(tt)} points, max active |sigma-mu| {np.max(np.abs(diff[active])):.6f}; '
          f'old/new error {err:.1e}; shift error {equiv_error:.1e}: PASS')
