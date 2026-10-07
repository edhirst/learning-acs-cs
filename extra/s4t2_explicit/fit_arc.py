"""Explicit smooth model of the zero arc Gamma_E of E (docs/s4t2_explicit_diffeo.md, section 12.2), written to
arc_fit.json and used by tlift.arc_model:

    arg t = pi + W1(rho) p1(rho),     x4 = x4_base(r) + W2(rho) p2(rho),     rho = log r,

  x4_base = (1 - B) * 1/2 + B * x4_zeroinf(r)   (the exact local-model zeros at the two cusps, B a smooth step),
  p1, p2 = Chebyshev fits of the residuals of the Newton-traced arc, W1, W2 = C-infinity windows equal to 1 on the
  supports of the residuals (which vanish identically outside [0.81, 1.24], resp. [0.366, 1.65]).

The model only enters Lambda through the straightening rotation R (any smooth model gives a diffeomorphism); its
accuracy decides how close Lambda(Gamma_E) is to the meridian M1 (check_tlift.py). Run from the repo root.
"""
from __future__ import annotations
import json, math, os, sys
import numpy as np
import torch
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
torch.set_default_dtype(torch.float64)
from torch.func import jvp
from extra.s4t2_explicit import tlift as L

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'arc_fit.json')


def E_on(r, th, x4):
    return L.E_W(torch.polar(r, th), torch.zeros_like(r), x4, torch.complex(torch.log(r), th))


def newton(r, th, x4, iters=60):
    for _ in range(iters):
        f = E_on(r, th, x4)
        _, fa = jvp(lambda a: E_on(r, a, x4), (th,), (torch.ones_like(th),))
        _, fb = jvp(lambda b: E_on(r, th, b), (x4,), (torch.ones_like(x4),))
        det = fa.real * fb.imag - fa.imag * fb.real
        da = (f.real * fb.imag - f.imag * fb.real) / det
        db = (fa.real * f.imag - fa.imag * f.real) / det
        sc = torch.clamp(0.05 / torch.maximum(da.abs(), db.abs()).clamp_min(1e-300), max=1.0)
        th, x4 = th - sc * da, x4 - sc * db
        if max(da.abs().max().item(), db.abs().max().item()) < 1e-15:
            break
    res = E_on(r, th, x4).abs() / (E_on(r, th, x4 + 0.25).abs() + E_on(r, th + 0.5, x4).abs())
    return th, x4, res


def trace(lr):
    """Newton on a log-grid by continuation in blocks (each block starts from the previous block's last point)."""
    r = torch.as_tensor(np.exp(lr))
    th = torch.full_like(r, math.pi)
    x4 = torch.where(r <= 1, torch.full_like(r, 0.5), L._x4_zeroinf(r))
    # initial guesses: exact zeros outside [0.36, 1.7]; inside, continuation in 40 blocks
    inside = ((r > 0.36) & (r < 1.7)).nonzero().flatten()
    blocks = np.array_split(inside.numpy(), 40)
    prev = None
    for blk in blocks:
        idx = torch.as_tensor(blk)
        if prev is not None:
            th[idx], x4[idx] = prev
        t_, x_, _ = newton(r[idx], th[idx], x4[idx])
        th[idx], x4[idx] = t_, x_
        prev = (t_[-1], x_[-1])
    th, x4, res = newton(r, th, x4)
    return r, th, x4, res


def cheb_fit(x, y, a, b, deg):
    xs = 2 * (x - a) / (b - a) - 1
    V = np.polynomial.chebyshev.chebvander(xs, deg)
    c, *_ = np.linalg.lstsq(V, y, rcond=None)
    return c, np.abs(V @ c - y).max()


if __name__ == '__main__':
    lr = np.linspace(math.log(0.2), math.log(3.0), 12001)
    r, th, x4, res = trace(lr)
    print(f'arc traced: {len(lr)} radii, max relative residual {res.max():.1e}')
    base = L.x4_base(r).numpy()
    f1, f2 = th.numpy() - math.pi, x4.numpy() - base
    print(f'max |arg t - pi| = {np.abs(f1).max():.4f}, max |x4 - x4_base| = {np.abs(f2).max():.4f}')
    out = {'provenance': 'extra/s4t2_explicit/fit_arc.py (Newton-traced zero arc of tlift.E_W, Chebyshev least squares)'}
    for key, f, (a, b) in (('p1', f1, L.ARC_INT1), ('p2', f2, L.ARC_INT2)):
        m = (lr >= a) & (lr <= b)
        outside = np.abs(f[~m]).max() if (~m).any() else 0.0
        for deg in (40, 60, 80, 100, 120):
            c, err = cheb_fit(lr[m], f[m], a, b, deg)
            print(f'{key}: degree {deg:3d}: max fit error {err:.2e} (residual outside the interval {outside:.1e})')
            if err < 1e-9:
                break
        out[key] = {'interval': [a, b], 'coef': c.tolist(), 'max_fit_error': err}
    json.dump(out, open(OUT, 'w'), indent=1)
    print('written', OUT)
