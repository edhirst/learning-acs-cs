"""Checks of the global weight-1 function E_W (middle.py)."""
from __future__ import annotations
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from extra.s4t2_explicit import middle as M

r0, ri = M.references()
for nm, R in (('cusp 0  ', r0), ('cusp inf', ri)):
    print(f'{nm} reference |t| = {R.r:.4f}: zero z = ({R.z[0]:+.6f}, {R.z[1]:+.6f}), |F(z)| = {R.zero_residual:.1e}, '
          f'windings (a, b) = ({R.a}, {R.b}), periodicity defect of log-ratio {R.periodicity_defect:.1e}')
rng = np.random.default_rng(5)
s, x = rng.uniform(-1, 2, 2000), rng.uniform(-1, 2, 2000)
for nm, R, F in (('cusp 0', r0, M.E0_W), ('cusp inf', ri, M.Einf_W)):
    ex_rec = M.theta_at(s, x, R.z[0], R.z[1]) * np.exp(2j * np.pi * (R.a * s + R.b * x) + R.l_exact(F, s, x))
    print(f'{nm}: exact-branch reconstruction: max rel. error {np.max(np.abs(ex_rec - M.F_cusp(F, R.r, s, x)) / np.abs(M.F_cusp(F, R.r, s, x)).max()):.1e}; '
          f'branch selection margin: max |spectral - exact| imaginary part = '
          f'{np.max(np.abs((R.l(s, x) - R.l_exact(F, s, x)).imag)):.2e} (must be << pi)')
    rec = M.theta_at(s, x, R.z[0], R.z[1]) * np.exp(2j * np.pi * (R.a * s + R.b * x) + R.l(s, x))
    ex = M.F_cusp(F, R.r, s, x)
    print(f'{nm}: spectral reconstruction of F on its reference torus: max rel. error {np.max(np.abs(rec - ex) / np.abs(ex).max()):.1e}')

# continuity of E_W across the radii where the definition changes, and twisted periodicities
def EW(r, s, x4, x1=0.0):
    t = r * np.exp(2j * np.pi * s)
    return M.E_W(np.asarray(t), np.asarray(x1 + 0 * s), np.asarray(x4), np.asarray(np.log(r) + 2j * np.pi * s))
s, x = rng.uniform(0, 1, 500), rng.uniform(0, 1, 500)
for rb in (M.R0, M.R1, 1 / M.R1, 1 / M.R0):
    lo, hi = EW(rb * (1 - 1e-9), s, x), EW(rb * (1 + 1e-9), s, x)
    print(f'continuity at |t| = {rb:.4f}: max |E(r-) - E(r+)| / max|E| = {np.max(np.abs(lo - hi)) / np.abs(lo).max():.1e}')
for r in (0.65, 0.9, 1.0, 1.1, 1.5):
    E0 = EW(r, s, x)
    print(f'|t| = {r}: twisted periodicity sigma -> sigma + 1: '
          f'{np.max(np.abs(EW(r, s + 1, x) - np.exp(2j*np.pi*x) * E0)) / np.abs(E0).max():.1e}; x4 -> x4 + 1: '
          f'{np.max(np.abs(EW(r, s, x + 1) - E0)) / np.abs(E0).max():.1e}')
# one zero per orbit-space torus across the middle
for r in np.exp(np.linspace(np.log(0.55), np.log(1 / 0.55), 13)):
    S, X = np.meshgrid(np.linspace(0, 1, 400, endpoint=False), np.linspace(0, 1, 400, endpoint=False), indexing='ij')
    G = EW(r, S.ravel(), X.ravel()).reshape(S.shape)
    g = np.abs(G) / np.abs(G).max()
    i0, j0 = np.unravel_index(np.argmin(g), g.shape)
    # winding around the boundary of the fundamental domain = number of zeros (twisted periodicity accounted:
    # going around [0,1]^2, the factor e^{2 pi i x} contributes -1, so zeros = winding + 1 ... count locally instead)
    small = g < 4 * np.sort(g.ravel())[40]
    lab = np.zeros(small.shape, int); nl = 0
    for (i, j) in zip(*np.nonzero(small)):
        if lab[i, j]: continue
        nl += 1; st = [(i, j)]
        while st:
            a, b = st.pop()
            a, b = a % 400, b % 400
            if small[a, b] and not lab[a, b]:
                lab[a, b] = nl; st += [(a + 1, b), (a - 1, b), (a, b + 1), (a, b - 1)]
    lam = M.bump_step((np.log(r) - np.log(M.R1)) / (np.log(1 / M.R1) - np.log(M.R1)))
    if 0 < lam < 1:
        th = M.theta_at(S.ravel(), X.ravel(), *((1 - lam) * r0.z + lam * ri.z)).reshape(S.shape)
        q = np.abs(G / th)
        extra = f' | blended: |E / Theta_z(r)| in [{q.min():.2e}, {q.max():.2e}] (> 0: exactly one zero)'
    else:
        extra = ''
    print(f'|t| = {r:.3f}: zero clusters {nl}, deepest at (sigma, x4) = ({S[i0, j0]:.3f}, {X[i0, j0]:.3f}), '
          f'min/max |E| = {g.min():.1e}, max |E| = {np.abs(G).max():.2e}' + extra)
