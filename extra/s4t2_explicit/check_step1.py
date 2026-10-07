"""Step 1 checks (docs/s4t2_explicit_diffeo.md, sections 1-2), on random points of Y_mid."""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from extra.s4t2_explicit.model import (lattice_coords, g, rot, torus, Phi, Theta, sigma_tilde)

rng = np.random.default_rng(0)
n = 2000
t = np.exp(rng.uniform(-3, 3, n) + 1j * rng.uniform(-np.pi + 0.05, np.pi - 0.05, n))
zeta = np.exp(rng.uniform(-2, 2, n) + 1j * rng.uniform(-np.pi, np.pi, n))
w = rng.normal(size=n) + 1j * rng.normal(size=n)
th = rng.uniform(0, 2 * np.pi, n)
c2, c3 = rng.uniform(-1, 1, n), rng.uniform(-1, 1, n)
z = -1j * np.log(zeta)

def dist(a, b):
    """max RELATIVE deviation (entries reach e^60 through e^{w/t})."""
    return max((np.abs(x - y) / np.maximum(np.abs(x), 1e-300)).max() for x, y in zip(a, b))

print('(a) S^1_a commutes with every identification (max relative deviation):')
for j in (2, 3, 4):
    print(f'    g{j}: {dist(rot(th, *g(j, zeta, w, t)), g(j, *rot(th, zeta, w, t))):.1e}')
print(f'    T^2: {dist(rot(th, *torus(c2, c3, zeta, w, t)), torus(c2, c3, *rot(th, zeta, w, t))):.1e}')
print(f'    Phi: {dist(rot(th, *Phi(zeta, w, t)), Phi(*rot(th, zeta, w, t))):.1e}')
print(f'    Theta (rotation of zeta\' at infinity): {dist(rot(th, *Theta(zeta, w, t)), Theta(*rot(th, zeta, w, t))):.1e}')

print('(b) action on the real lattice coordinates x = (x1..x4):')
x = lattice_coords(z, w, t)
xr = lattice_coords(z + th, w, t)
print(f'    S^1_a : dx = (theta/2pi, 0, 0, 0) to {np.abs(xr - x - np.stack([th / (2 * np.pi), 0 * th, 0 * th, 0 * th], -1)).max():.1e}')
zt, wt, _ = torus(c2, c3, zeta, w, t)
xt = lattice_coords(z + (-1j) * (c2 - 1j * c3) * t, wt, t)        # additive version of the T^2 multiplier
print(f'    T^2   : dx = (0, c2, c3, 0) to {np.abs(xt - x - np.stack([0 * c2, c2, c3, 0 * c2], -1)).max():.1e}')
z4 = z - 1j * np.log(t)                                            # g4 in additive coords: zeta -> t zeta
print(f'    g4    : dx = (0, 0, 0, 1) to {np.abs(lattice_coords(z4, w + 2 * t, t) - x - np.array([0, 0, 0, 1.0])).max():.1e}')
xk = lattice_coords(z, w, t, k=1)
print(f'    branch log t -> log t + 2 pi i: x4 unchanged to {np.abs(xk[:, 3] - x[:, 3]).max():.1e}, '
      f'x1 shift = -x4 to {np.abs(xk[:, 0] - x[:, 0] + x[:, 3]).max():.1e}')

print('(c) clutching of the log transform: xi = lattice coords of sigma_tilde, continued once around |t| = r')
for r in (1e-2, 0.3, 3.0):
    phis = np.linspace(-np.pi, np.pi, 2001)
    xi4 = []
    for p in phis:                                    # continuous branch: log t = ln r + i p
        tt = r * np.exp(1j * p)
        k = 0 if p < np.pi else 0
        s = (np.log(r) + 1j * p) / (2j * np.pi)
        sz, sw = np.pi * s * (s - 1), 2 * tt * s
        # lattice with the same continuous log
        from extra.s4t2_explicit.model import lattice as _lat
        E = _lat(tt)
        E[0, 3] = -1j * (np.log(r) + 1j * p)
        R = np.stack([E[0].real, E[0].imag, E[1].real, E[1].imag])
        xi4.append(np.linalg.solve(R, np.array([sz.real, sz.imag, sw.real, sw.imag]))[3])
    xi4 = np.array(xi4)
    print(f'    r = {r:g}: xi4(arg=pi) - xi4(arg=-pi) = {xi4[-1] - xi4[0]:+.12f}  (expected +1: Euler number of beta = +-1)')

print('(d) coordinates at infinity: x4 at infinity = -x4, and (x2, x3) unchanged')
zp, wp, s = Theta(zeta, w, t)
z2 = z - 1j * w / t                                               # additive z'' with zeta' = e^{i z''}
E = np.zeros((n, 2, 4), complex)
E[:, 0, 0] = 2 * np.pi
E[:, 0, 1], E[:, 1, 1] = -1j * s, s ** 2 - 1                      # (-i alpha2, e2bar), alpha2 = s
E[:, 0, 2], E[:, 1, 2] = -1j * (1j * s), 1j * (1 + s ** 2)        # alpha3 = i s
E[:, 0, 3], E[:, 1, 3] = -1j * (np.log(s) - 2.0), -2 * s          # -i log(nu s), nu = e^-2
R = np.stack([E[:, 0].real, E[:, 0].imag, E[:, 1].real, E[:, 1].imag], 1)
xinf = np.linalg.solve(R, np.stack([z2.real, z2.imag, wp.real, wp.imag], 1)[..., None])[..., 0]
d4 = xinf[:, 3] + x[:, 3]
d23 = xinf[:, 1:3] - x[:, 1:3]
print(f'    x4_inf + x4 is an integer: max dist to Z = {np.abs(d4 - np.round(d4)).max():.1e}; '
      f'(x2, x3)_inf - (x2, x3) integer: {np.abs(d23 - np.round(d23)).max():.1e}')
print(f'    zero section [1, 0, t]: x = {np.abs(lattice_coords(0 * t, 0 * t, t)).max():.1e}')
