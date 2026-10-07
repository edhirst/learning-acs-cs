"""Branch normaliser psi(|t|) of the gauge G_A = e^{i psi} G_A_hat (docs/s4t2_explicit_diffeo.md, section 12.4),
written to psi_fit.json and used by tlift.psi_of_r.

On the shell THETA_A <= theta_v <= THETA_B the two lifts differ by the phase g = lift_B / lift_A, which has a
continuous logarithm (winding 0 around the tube); the glued lift is e^{i c Arg g} lift_A. Arg must be the
continuous branch, which holds iff |Arg(g)| < pi on the shell. psi is chosen to track the circular mean of g along
the tube, with the exact values -pi/2 near p0 and pi + 2 pi k near p_inf (k from continuation), where the shell
lies in the linear balls and g must be exactly 1. Only the branch matters: any psi within the margin gives a
smooth lift; the margin is reported. Run from the repo root.
"""
from __future__ import annotations
import json, math, os, sys
import numpy as np
import torch
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
torch.set_default_dtype(torch.float64)
from extra.s4t2_explicit import tcore as T, tlift as L
from extra.s4t2_explicit.cusp import bump_step

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'psi_fit.json')
gen = torch.Generator().manual_seed(1234)
U = lambda n, a, b: torch.empty(n).uniform_(a, b, generator=gen)


def shell_samples_chart(C, Lam, Ef, Gf, n, lo, hi, inv):
    rr = torch.exp(U(n, math.log(lo), math.log(hi)))
    al = math.pi + U(n, -1.3, 1.3)
    ratio = torch.exp(U(n, -2.5, 2.5)) * math.sqrt(C.a0)
    bx = U(n, -math.pi, math.pi)
    x = torch.polar(torch.sqrt(rr / ratio), bx)
    y = torch.polar(torch.sqrt(rr * ratio), al - bx)
    k = (C.mu_chart(x, y).abs() < 0.45) & ((C.a0 * x.abs() ** 2 + y.abs() ** 2) > L.EPS1)
    x, y = x[k], y[k]
    v, X = Lam(x, y)
    r = (x * y).abs()
    r = 1 / r if inv else r
    return v, L.g_shell(v, X, Ef(x, y), Gf(x, y), r, use_psi=False), torch.log(r)


def shell_samples_mid(n, lo, hi):
    r = torch.exp(U(n, math.log(lo), math.log(hi)))
    th0, x40 = L.arc_model(r)
    th = th0 + U(n, -1.3, 1.3)
    x4 = x40 + U(n, -0.3, 0.3)
    t = torch.polar(r, th)
    lg = torch.complex(torch.log(r), th)
    x1 = torch.zeros_like(r)
    v, X = L.Lam_W(t, x1, x4, lg)
    return v, L.g_shell(v, X, L.E_W(t, x1, x4, lg), L.GA_hat_W(t, x1, x4, lg), r, use_psi=False), torch.log(r)


def collect(batches=8, n=400000):
    out = []
    for _ in range(batches):
        for v, gg, lr in (shell_samples_chart(T.CUSP0, L.Lam_chart0, L.E_chart0, L.GA_hat_chart0, n, 1e-3, 0.25, False),
                          shell_samples_mid(n, 0.2, 5.0),
                          shell_samples_chart(T.CUSPINF, L.Lam_chartinf, L.E_chartinf, L.GA_hat_chartinf, n, 1e-2, 0.25, True)):
            thv = torch.acos((-v[:, 2] / v.norm(dim=1)).clamp(-1, 1))
            sh = (thv >= L.THETA_A) & (thv <= L.THETA_B)
            out.append((gg[sh], lr[sh]))
    return torch.cat([o[0] for o in out]), torch.cat([o[1] for o in out])


def psi_base(rho, a, b, psi_hi):
    return L.PSI_LO + (psi_hi - L.PSI_LO) * bump_step((rho - a) / (b - a))


if __name__ == '__main__':
    gg, lr = collect()
    print(f'shell samples: {len(gg)}')
    a, b = math.log(L.PSI_R[0]), math.log(L.PSI_R[1])
    grid = np.arange(a, b + 1e-9, 0.01)
    h = 0.02
    means, prev, counts = [], None, []
    for c in grid:
        m = (lr - c).abs() < h
        counts.append(int(m.sum()))
        if m.sum() < 5:
            means.append(np.nan)
            continue
        cm = torch.angle(gg[m].mean()).item()
        if prev is not None:
            cm += 2 * math.pi * round((prev - cm) / (2 * math.pi))
        prev = cm
        means.append(cm)
    means = np.array(means)
    ok = np.isfinite(means)
    print(f'grid points without samples: {int((~ok).sum())}; min samples per window {min(counts)}')
    jumps = np.abs(np.diff(means[ok])).max()
    print(f'mean phase: {means[ok][0]:+.4f} -> {means[ok][-1]:+.4f}; largest step between grid points {jumps:.3f} rad')
    k = round((means[ok][-1] - math.pi) / (2 * math.pi))
    psi_hi = math.pi + 2 * math.pi * k
    print(f'psi_hi = pi + 2 pi * {k}')
    edge = 0.5
    res = means[ok] - psi_base(grid[ok], a, b, psi_hi)
    xs = 2 * (grid[ok] - a) / (b - a) - 1
    win = bump_step((grid[ok] - a) / edge) * bump_step((b - grid[ok]) / edge)
    best = None
    for deg in (20, 30, 40, 60):
        V = np.polynomial.chebyshev.chebvander(xs, deg) * win[:, None]
        coef = np.linalg.solve(V.T @ V + 1e-8 * np.eye(deg + 1), V.T @ res)
        L._PSI = {'psi_hi': psi_hi, 'coef': coef.tolist(), 'edge': edge}
        dev = torch.angle(gg * torch.polar(torch.ones(1), -L.psi_of_r(torch.exp(lr)))).abs()
        fit_err = np.abs(V @ coef - res).max()
        print(f'degree {deg}: fit error of the mean {fit_err:.3f} rad; max |Arg(g e^-i psi)| on the shell = {dev.max():.3f} rad '
              f'(99.9% quantile {torch.quantile(dev[:1000000], 0.999):.3f})')
        if best is None or dev.max() < best[0]:
            best = (dev.max().item(), deg, coef)
    margin, deg, coef = best
    json.dump({'provenance': 'extra/s4t2_explicit/fit_psi.py', 'psi_hi': psi_hi, 'edge': edge, 'degree': deg,
               'coef': coef.tolist(), 'max_abs_Arg_on_shell': margin, 'shell_samples': len(gg)}, open(OUT, 'w'), indent=1)
    print(f'written {OUT}: degree {deg}, margin pi - {margin:.3f} = {math.pi - margin:.3f} rad')
