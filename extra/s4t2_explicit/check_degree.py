"""Degree of Psi_W: W -> S^4 as a volume integral: deg = (1/vol S^4) int_W Psi^* vol.
W minus the two central fibres is covered once by (log r, arg t, x1, x4) in R x (-pi, pi] x [0,1)^2; the integrand
det[Psi, d Psi/d(log r, arg t, x1, x4)] is S^1-invariant (independent of x1), so the integral is a 3D quadrature
(trapezoid in the periodic variables, trapezoid in log r on [-L, L]; the integrand decays like r^2 at both cusps).
Run from the repo root: python3 extra/s4t2_explicit/check_degree.py [n_rho n_theta n_x4]"""
from __future__ import annotations
import math, os, sys, time
import torch
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
torch.set_default_dtype(torch.float64)
from extra.s4t2_explicit import check_lift as CL

VOL_S4 = 8 * math.pi ** 2 / 3


def degree(n_rho, n_th, n_x4, lo=-16.0, hi=12.0, chunk=200000):
    rho = torch.linspace(lo, hi, n_rho)
    th = -math.pi + 2 * math.pi * (torch.arange(n_th) + 0.5) / n_th
    x4 = (torch.arange(n_x4) + 0.5) / n_x4
    R, TH, X4 = torch.meshgrid(rho, th, x4, indexing='ij')
    P = torch.stack([R.flatten(), TH.flatten(), torch.zeros(R.numel()), X4.flatten()], -1)
    vals = []
    for i in range(0, len(P), chunk):
        d, _ = CL.det5(CL.psi_mid_real, P[i:i + chunk])
        vals.append(d)
    D = torch.cat(vals).reshape(n_rho, n_th, n_x4)
    prof = D.mean(dim=(1, 2)) * 2 * math.pi                     # integral over theta and x4 at each rho
    w = torch.full((n_rho,), (hi - lo) / (n_rho - 1)); w[0] = w[-1] = w[0] / 2
    return (prof * w).sum().item() / VOL_S4, D, prof, rho


if __name__ == '__main__':
    args = [int(a) for a in sys.argv[1:4]] or [561, 64, 64]
    t0 = time.time()
    deg, D, prof, rho = degree(*args)
    print(f'grid {args}: degree of Psi_W = {deg:.6f}  (det > 0 at {int((D > 0).sum())}/{D.numel()} nodes; '
          f'{time.time() - t0:.0f} s)')
    print(f'   integrand at the ends of the log r range: {prof[0]:.1e}, {prof[-1]:.1e}')
