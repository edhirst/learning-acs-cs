"""Trace Gamma_E through all scales; image under Moebius o h_B; radial-graph property."""
from __future__ import annotations
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from extra.s4t2_explicit.curve import gamma_point, moebius

rs = np.exp(np.linspace(np.log(1e-6), np.log(1e6), 601))
out, guess, prev_rep = [], None, None
for r in rs:
    rep, data, u = gamma_point(r, guess if True else None)
    if rep != prev_rep:
        guess = None
        rep, data, u = gamma_point(r, None)
    res = data[-1]
    if rep == 'mid':
        guess = np.array([np.angle(data[0]) % (2 * np.pi), data[1]])
    else:
        guess = np.array([np.angle(data[1]) % (2 * np.pi), np.log(abs(data[0]))])
    prev_rep = rep
    out.append((r, rep, res, u, np.angle(data[1] if rep != 'mid' else data[0])))
U = np.array([o[3] for o in out]); W = moebius(U); rho = np.linalg.norm(W, axis=1)
res = np.array([o[2] for o in out]); arg = np.array([o[4] for o in out])
print(f'zero residuals: max {res.max():.1e};  |arg t| - pi: max {np.abs(np.abs(arg) - np.pi).max():.1e}')
print(f'|u| - 1: {np.abs(np.linalg.norm(U, axis=1) - 1).max():.1e};  endpoints u(1e-6) = {np.round(U[0], 6)}, u(1e6) = {np.round(U[-1], 6)}')
d = np.diff(np.log(rho))
print(f'radial graph: d log|w| / d log r in [{(d / np.diff(np.log(rs))).min():.3e}, {(d / np.diff(np.log(rs))).max():.3e}] (>0 => meets each sphere once)')
dirs = W / rho[:, None]
ang = np.degrees(np.arccos(np.clip(dirs @ np.array([0, 1.0, 0]), -1, 1)))
for lo, hi in ((1e-6, 1e-3), (1e-3, 0.3), (0.3, 3), (3, 1e3), (1e3, 1e6)):
    m = (rs >= lo) & (rs <= hi)
    print(f'r in [{lo:g}, {hi:g}]: |w| in [{rho[m].min():.2e}, {rho[m].max():.2e}], angle to the ray d = (0,1,0): max {ang[m].max():.3f} deg')
