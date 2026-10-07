"""Checks of the conversions in glue.py (lattice coordinates <-> cusp slices)."""
from __future__ import annotations
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from extra.s4t2_explicit.glue import E0_W, Einf_W

rng = np.random.default_rng(4)
n = 3000
for name, F, rr in (('E0_W  ', E0_W, (0.05, 0.95)), ('Einf_W', Einf_W, (1.05, 20.0))):
    r = np.exp(rng.uniform(np.log(rr[0]), np.log(rr[1]), n))
    a = rng.uniform(-np.pi, np.pi, n)
    t, logt = r * np.exp(1j * a), np.log(r) + 1j * a
    x1, x4 = rng.uniform(0, 1, n), rng.uniform(0, 1, n)
    th = rng.uniform(0, 1, n)
    E = F(t, x1, x4, logt)
    rel = lambda A: np.max(np.abs(A - E) / np.abs(E))
    print(f'{name}: weight 1 {np.max(np.abs(F(t, x1 + th, x4, logt) - np.exp(2j*np.pi*th) * E) / np.abs(E)):.1e} | '
          f'e1-period {rel(F(t, x1 + 1, x4, logt)):.1e} | e4-period {rel(F(t, x1, x4 + 1, logt)):.1e} | '
          f'branch change {rel(F(t, x1 - x4, x4, logt + 2j*np.pi)):.1e} | min |E| {np.abs(E).min():.1e}')
# zero of x4 -> F(t, 0, x4) on the circles near |t| = 1 from both sides: one zero per torus (sigma, x4)?
for name, F, rs in (('E0_W', E0_W, (0.7, 0.9, 0.97)), ('Einf_W', Einf_W, (1.03, 1.1, 1.4))):
    for r in rs:
        A = np.linspace(-np.pi, np.pi, 720, endpoint=False) + np.pi / 720
        X4 = np.linspace(0, 1, 720, endpoint=False)
        AA, XX = np.meshgrid(A, X4, indexing='ij')
        G = F(r * np.exp(1j * AA), 0 * AA, XX, np.log(r) + 1j * AA)
        g = np.abs(G) / np.abs(G).max()
        i0, j0 = np.unravel_index(np.argmin(g), g.shape)
        print(f'{name} |t| = {r:4.2f}: min |F|/max = {g.min():.1e} at arg t = {A[i0]:+.3f}, x4 = {X4[j0]:.3f}; '
              f'max |F| = {np.abs(G).max():.2e}, cells < 1e-2: {(g < 1e-2).sum()}')
