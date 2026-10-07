"""Step 2 checks for both cusps: the height mu and the weight-1 function E (extra/s4t2_explicit/cusp.py)."""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from extra.s4t2_explicit.cusp import CUSP0, CUSPINF


def check(C, seed):
    print(f'================ cusp {C.name}  (nu = {C.nu:.4g})')
    rng = np.random.default_rng(seed)
    n = 4000
    t = np.exp(rng.uniform(np.log(1e-6), np.log(0.5), n) + 1j * rng.uniform(-np.pi, np.pi, n))
    rho = rng.uniform(-3, 3, n) * np.abs(np.log(np.abs(t)))
    zeta = np.exp(rho + 1j * rng.uniform(-np.pi, np.pi, n))
    m0 = C.mu_torus(zeta, t)
    print(f'mu: equivariance {np.abs(C.mu_torus(*C.psi_hat(zeta, t)) - m0 - 1).max():.1e}', end='')
    dev = 0.0
    k = np.floor(np.log(np.abs(zeta)) / np.log(np.abs(t))).astype(int)
    for kk in np.unique(k):
        s = k == kk
        for kc in (kk - 1, kk, kk + 1):
            x, y = zeta[s] ** -1 * t[s] ** (kc + 1), zeta[s] * t[s] ** (-kc)
            ok = (np.abs(x) < 2) & (np.abs(y) < 2)
            if ok.any():
                dev = max(dev, np.abs(C.mu_chart(x[ok], y[ok], kc) - m0[s][ok]).max())
    print(f' | torus = chart {dev:.1e}', end='')
    h = 1e-6
    dmu = (C.mu_torus(zeta * np.exp(h), t) - C.mu_torus(zeta * np.exp(-h), t)) / (2 * h)
    print(f' | d mu/d log|zeta| in [{dmu.min():.2e}, {dmu.max():.2e}]', end='')
    d = m0 - (np.log(np.abs(zeta)) / C.P(t) - 0.5)
    print(f' | mu - (x4 - 1/2) in [{d.min():+.3f}, {d.max():+.3f}]')
    nu2 = abs(C.nu) ** 2
    r = 10 ** rng.uniform(-4, -1.5, 1000)
    x = r * np.exp(1j * rng.uniform(0, 2 * np.pi, 1000)) * rng.uniform(0, 1, 1000)
    y = r * np.exp(1j * rng.uniform(0, 2 * np.pi, 1000)) * rng.uniform(0, 1, 1000)
    loc = C.mu_chart(x, y, 0) - (nu2 * np.abs(x) ** 2 - np.abs(y) ** 2)
    print(f'mu near the node: mu = |nu|^2 |x|^2 - |y|^2 + O(4), remainder/(|x|^2+|y|^2)^2 <= '
          f'{np.max(np.abs(loc) / (np.abs(x) ** 2 + np.abs(y) ** 2) ** 2):.2f}')

    E = C.E_torus(zeta, t)
    lam = np.exp(1j * rng.uniform(0, 2 * np.pi, n))
    print(f'E: weight 1 {np.max(np.abs(C.E_torus(lam * zeta, t) - lam * E) / np.abs(E)):.1e}'
          f' | psi_hat-invariant {np.max(np.abs(C.E_torus(*C.psi_hat(zeta, t)) - E) / np.abs(E)):.1e}')
    for r in (1e-4, 1e-2, 0.3, 0.6):
        na, nm = 360, 720
        A = np.linspace(-np.pi, np.pi, na, endpoint=False) + np.pi / na
        mus = np.linspace(-0.5, 0.5, nm, endpoint=False)
        rhos = np.linspace(-2.5, 1.5, 6000) * max(abs(np.log(r)) + abs(np.log(abs(C.nu))), 1.0)
        G = np.zeros((na, nm), complex)
        for i, ti in enumerate(r * np.exp(1j * A)):
            mm = C.mu_torus(np.exp(rhos), np.full(rhos.shape, ti))
            o = np.argsort(mm)
            z = np.exp(np.interp(mus, mm[o], rhos[o]))
            G[i] = C.E_torus(z + 0j, np.full(z.shape, ti)) * np.conj(z) / np.abs(z)
        absG = np.abs(G) / np.abs(G).max()
        # connected components of the small set (count zero clusters)
        small = absG < 3 * np.sort(absG.ravel())[min(60, absG.size - 1)]
        lab = np.zeros(small.shape, int); nlab = 0
        for (i, j) in zip(*np.nonzero(small)):
            if lab[i, j]:
                continue
            nlab += 1; stack = [(i, j)]
            while stack:
                a, b = stack.pop()
                if small[a % na, b % nm] and not lab[a % na, b % nm]:
                    lab[a % na, b % nm] = nlab
                    stack += [(a + 1, b), (a - 1, b), (a, b + 1), (a, b - 1)]
        i0, j0 = np.unravel_index(np.argmin(absG), absG.shape)
        loop = [(i0 - 8, j0 - 8), (i0 - 8, j0 + 8), (i0 + 8, j0 + 8), (i0 + 8, j0 - 8)]
        path = np.array([G[int(round(a1 + s * (a2 - a1))) % na, int(round(b1 + s * (b2 - b1))) % nm]
                         for (a1, b1), (a2, b2) in zip(loop, loop[1:] + loop[:1]) for s in np.linspace(0, 1, 64, endpoint=False)])
        wind = (np.sum(np.angle(path[1:] / path[:-1])) + np.angle(path[0] / path[-1])) / (2 * np.pi)
        print(f'E zeros on |t| = {r:g}: {nlab} cluster(s); min at arg t = {A[i0]:+.3f}, mu = {mus[j0]:+.3f}; '
              f'winding {wind:+.2f}; min/max |G| = {absG.min():.1e}')
    print('E smooth across the central fibre (chart U_0: zeta = y, t = x y), values and derivatives as x -> 0, y -> 0:')
    f = lambda x, y: C.E_torus(np.array([y + 0j]), np.array([x * y + 0j]))[0]
    for vary in ('x', 'y'):
        out = []
        for sm in (1e-3, 1e-5, 1e-7, 1e-9):
            x0, y0 = (sm * (1 + 1j), 0.7 + 0.2j) if vary == 'x' else (0.7 - 0.3j, sm * (1 - 1j))
            hh = 1e-6
            dx = (f(x0 + hh, y0) - f(x0 - hh, y0)) / (2 * hh)
            dy = (f(x0, y0 + hh) - f(x0, y0 - hh)) / (2 * hh)
            out.append((f(x0, y0), dx, dy))
        dv = max(abs(a - b) for u, v in zip(out[-2], out[-1]) for a, b in [(u, v)])
        print(f'   {vary} -> 0: E = {out[-1][0]:+.5f}, dE/dx = {out[-1][1]:+.5f}, dE/dy = {out[-1][2]:+.5f}; '
              f'change between 1e-7 and 1e-9: {dv:.1e}')


check(CUSP0, 2)
check(CUSPINF, 3)
