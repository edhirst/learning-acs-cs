"""Checks of Lambda (orbit-space map) and of the lift."""
from __future__ import annotations
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from extra.s4t2_explicit import lam as L
from extra.s4t2_explicit import hopf as H
from extra.s4t2_explicit.cusp import CUSP0, CUSPINF

rng = np.random.default_rng(9)
# 1. the zero arc goes exactly onto the z1 = 0 meridian {v1 = v2 = 0, v3 < 0}
rs = np.exp(np.linspace(np.log(1e-5), np.log(1e5), 241))
u = L.gamma_u(rs)
v, X = L.Lambda(u, rs)
print(f'1. Lambda(Gamma_E): max |(v1, v2)| = {np.abs(v[:, :2]).max():.1e}; v3 < 0: {bool((v[:, 2] < 0).all())}; '
      f'X from {X.min():+.4f} to {X.max():+.4f} (min step {np.diff(X).min():+.1e}, max step {np.diff(X).max():+.1e})')
# 2. conformality + orientation at b0 (v0 = (2 t-bar, |x|^2 - |y|^2)) and b_inf (v_inf = (2|nu| s-bar, |nu|^2|x|^2 - |y|^2))
for nm, C in (('b0  ', CUSP0), ('binf', CUSPINF)):
    a0 = abs(C.nu) ** 2
    eps = 1e-4
    x = eps * (rng.normal(size=300) + 1j * rng.normal(size=300)); y = eps * (rng.normal(size=300) + 1j * rng.normal(size=300))
    if nm == 'b0  ':
        uu, r = H.h_B_c0(y, x * y), np.abs(x * y)
        vv = np.stack([2 * (x * y).real, -2 * (x * y).imag, np.abs(x) ** 2 - np.abs(y) ** 2], 1)
    else:
        uu, r = H.h_B_ci(y, x * y), 1 / np.abs(x * y)
        sb = np.conj(x * y)
        vv = np.stack([2 * np.sqrt(a0) * sb.real, 2 * np.sqrt(a0) * sb.imag, a0 * np.abs(x) ** 2 - np.abs(y) ** 2], 1)
    vS, XS = L.Lambda(uu, r)
    M, *_ = np.linalg.lstsq(vv, vS, rcond=None)
    sv = np.linalg.svd(M, compute_uv=False)
    print(f'2. {nm}: X -> {XS.min():+.6f}..{XS.max():+.6f}; v_S = M v: singular values {np.round(sv, 4)}, det M sign {np.sign(np.linalg.det(M)):+.0f} '
          f'(must be +1: both v-coordinates have the same orientation character), residual {np.abs(vv @ M - vS).max() / np.abs(vS).max():.1e}')
# 3. Jacobian of Lambda in the middle representation (t, x4): sign constant?
n = 4000
r = np.exp(rng.uniform(np.log(0.5), np.log(2.0), n)); a = rng.uniform(-np.pi, np.pi, n); x4 = rng.uniform(0, 1, n)
def Lm(r, a, x4):
    t = r * np.exp(1j * a)
    v, X = L.Lambda(H.h_B_mid(t, x4), r)
    return np.concatenate([v, X[:, None]], 1)
h = 1e-6
base = Lm(r, a, x4)
cols = [(Lm(r * np.exp(h), a, x4) - Lm(r * np.exp(-h), a, x4)) / (2 * h), (Lm(r, a + h, x4) - Lm(r, a - h, x4)) / (2 * h),
        (Lm(r, a, x4 + h) - Lm(r, a, x4 - h)) / (2 * h)]
# B_S is 3-dimensional inside R^4: project the Jacobian onto the tangent space via the outward normal of |v| + X^2 = 1
vv, XX = base[:, :3], base[:, 3]
nrm = np.concatenate([vv / np.linalg.norm(vv, axis=1, keepdims=True), 2 * XX[:, None]], 1)
J4 = np.stack(cols + [nrm], -1)
dets = np.linalg.det(J4)
print(f'3. middle: det[dLambda | normal] sign: +{(dets > 0).sum()} / -{(dets < 0).sum()}; min |det| = {np.abs(dets).min():.2e}')
