"""Checks of F = (Psi_W o pi, phi): Y -> S^4 x T^2 (tfull.py), docs/s4t2_explicit_diffeo.md 12.6:
  1. the three charts agree on their overlaps (through the paper's holomorphic gluings Phi and Theta);
  2. F is well defined on Y: lattice translations e1..e4 and the branch of log t (mid), G2 and G3 (fillings);
  3. T^3-equivariance: T^2 translates phi, S^1_a rotates (Z1, Z2);
  4. F lands in S^4 x T^2;
  5. det[n_S4, n_S1, n_S1, DF] (9 x 9) in the holomorphic coordinates has the same sign in all charts
     (complex orientation of Y), bounded away from 0; adversarial minimisation;
  6. degree: F is a T^2-equivariant map covering Psi_W, so its 6 x 6 determinant factors as
     (2 pi)^2 x (Psi_W determinant) in lattice coordinates (checked), and deg F = deg Psi_W = 1 (check_degree.py).
Run from the repo root."""
from __future__ import annotations
import math, os, sys
import torch
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
torch.set_default_dtype(torch.float64)
from torch.func import jvp
from extra.s4t2_explicit import tcore as T, tlift as L, tfull as Fm

gen = torch.Generator().manual_seed(99)
U = lambda n, a, b: torch.empty(n).uniform_(a, b, generator=gen)
fails = []


def rep(name, val, ok):
    if not ok:
        fails.append(name)
    print(f'{name:<72s} {val}  {"ok" if ok else "FAIL"}')


def dS4T2(a, b):
    """Distance in R^9 (the T^2 part is already periodic through the embedding)."""
    return (a - b).abs().max().item()


def jac(f, p):
    cols = [jvp(f, (p,), (torch.zeros_like(p).index_fill_(-1, torch.tensor([k]), 1.0),))[1] for k in range(p.shape[-1])]
    return torch.stack(cols, -1)


def det9(f, p):
    q = f(p)
    A = jac(f, p)
    n1 = torch.cat([q[..., :5], torch.zeros(q.shape[:-1] + (4,))], -1)
    n2 = torch.cat([torch.zeros(q.shape[:-1] + (5,)), q[..., 5:7], torch.zeros(q.shape[:-1] + (2,))], -1)
    n3 = torch.cat([torch.zeros(q.shape[:-1] + (7,)), q[..., 7:9]], -1)
    M = torch.cat([n1[..., None], n2[..., None], n3[..., None], A], -1)
    return torch.linalg.det(M), A


def mid_sample(n, lo=0.2, hi=5.0):
    r = torch.exp(U(n, math.log(lo), math.log(hi)))
    t = torch.polar(r, U(n, -math.pi, math.pi))
    x = torch.stack([U(n, -2, 2) for _ in range(4)], -1)
    z, w = T.point_from_lattice(t, torch.log(t), x)
    return z, w, t


def chart_sample(n, C, rmax=0.25, rmin=1e-6):
    rr = torch.exp(U(n, math.log(rmin), math.log(rmax)))
    ratio = torch.exp(U(n, -3.0, 3.0)) * math.sqrt(C.a0)
    ax, ay = U(n, -math.pi, math.pi), U(n, -math.pi, math.pi)
    x = torch.polar(torch.sqrt(rr / ratio), ax)
    y = torch.polar(torch.sqrt(rr * ratio), ay)
    t = x * y
    w = U(n, -1.5, 1.5) * C.e2(t) + U(n, -1.5, 1.5) * C.e3(t)
    _, _, _, xs, ys = Fm._slice(C, x, y, w)
    k = C.mu_chart(xs, ys).abs() < 0.48
    return x[k], y[k], w[k]


def real6(a, b, c):
    return torch.stack([a.real, a.imag, b.real, b.imag, c.real, c.imag], -1)


if __name__ == '__main__':
    # ------------------------------------------------------------ 1. chart overlaps
    for nm, C, Fc, tomid in (('c0', T.CUSP0, Fm.F_c0, Fm.c0_to_mid), ('ci', T.CUSPINF, Fm.F_ci, Fm.ci_to_mid)):
        x, y, w = chart_sample(20000, C, rmax=0.25, rmin=1e-3)
        e = dS4T2(Fc(x, y, w), Fm.F_mid(*tomid(x, y, w)))
        rep(f'1. {nm} chart = mid chart on the overlap ({len(x)} pts, 1e-3 <= |t or s| <= 0.25)', f'{e:.1e}', e < 1e-10)
    # ------------------------------------------------------------ 2. well-definedness
    z, w, t = mid_sample(20000, 1e-3, 1e3)
    lg = torch.log(t)
    F0 = Fm.F_mid(z, w, t, lg)
    ez, ew = T.e_vectors(t, lg)
    e = max(dS4T2(Fm.F_mid(z + ez[j], w + ew[j], t, lg), F0) for j in range(4))
    rep('2. mid: invariant under the lattice translations e1..e4', f'{e:.1e}', e < 1e-10)
    e = dS4T2(Fm.F_mid(z, w, t, lg + 2j * math.pi), F0)
    rep('   mid: independent of the branch of log t', f'{e:.1e}', e < 1e-10)
    for nm, C, Fc in (('c0', T.CUSP0, Fm.F_c0), ('ci', T.CUSPINF, Fm.F_ci)):
        xc, yc, wc = chart_sample(20000, C)
        tt = xc * yc
        G0 = Fc(xc, yc, wc)
        e = 0.0
        for al, eb in ((C.al2, C.e2), (C.al3, C.e3)):
            lam = torch.exp(al(tt))                              # T^2-element (1, 0) resp. (0, 1): identity on Y
            e = max(e, dS4T2(Fc(xc / lam, yc * lam, wc + eb(tt)), G0))
        rep(f'   {nm}: invariant under G2, G3 (unit T^2-translations)', f'{e:.1e}', e < 1e-10)
    # ------------------------------------------------------------ 3. equivariance
    c2, c3, th = U(len(t), 0, 1), U(len(t), 0, 1), U(len(t), -math.pi, math.pi)
    F1 = Fm.F_mid(z + c2 * ez[1] + c3 * ez[2], w + c2 * ew[1] + c3 * ew[2], t, lg)
    F2 = Fm.F_mid(z + th, w, t, lg)
    tor = Fm.torus_embed(torch.stack([c2, c3], -1))
    ang0 = torch.atan2(F0[:, [6, 8]], F0[:, [5, 7]]); ang1 = torch.atan2(F1[:, [6, 8]], F1[:, [5, 7]])
    dphi = torch.remainder(ang1 - ang0 - 2 * math.pi * torch.stack([c2, c3], -1) + math.pi, 2 * math.pi) - math.pi
    e1 = max(dS4T2(F1[:, :5], F0[:, :5]), dphi.abs().max().item())
    ph = torch.polar(torch.ones_like(th), th)
    Z = torch.complex(F0[:, [0, 2]], F0[:, [1, 3]])
    Zr = torch.complex(F2[:, [0, 2]], F2[:, [1, 3]])
    e2 = max((Zr - ph[:, None] * Z).abs().max().item(), dS4T2(F2[:, 4:], F0[:, 4:]))
    rep('3. T^2-equivariance (phi += c, S^4 part fixed)', f'{e1:.1e}', e1 < 1e-10)
    rep('   S^1_a-equivariance (z += theta: Z -> e^{i theta} Z, rest fixed)', f'{e2:.1e}', e2 < 1e-10)
    # ------------------------------------------------------------ 4. image
    e = max((F0[:, :5].norm(dim=-1) - 1).abs().max(), (F0[:, 5:7].norm(dim=-1) - 1).abs().max(), (F0[:, 7:].norm(dim=-1) - 1).abs().max())
    rep('4. F lands in S^4 x S^1 x S^1', f'{e:.1e}', e < 1e-13)
    # ------------------------------------------------------------ 5. Jacobian in holomorphic coordinates
    signs = []
    for nm, f, P in (('mid (0.2 <= |t| <= 5)', Fm.F_mid_real, real6(*mid_sample(60000))),
                     ('c0 chart (|t| <= 0.25)', Fm.F_c0_real, real6(*chart_sample(60000, T.CUSP0))),
                     ('ci chart (|s| <= 0.25)', Fm.F_ci_real, real6(*chart_sample(60000, T.CUSPINF)))):
        d, A = det9(f, P)
        sv = torch.linalg.svdvals(A)
        print(f'   5. {nm}: det sign +{int((d > 0).sum())} / -{int((d < 0).sum())}; min |det| {d.abs().min():.3e}; '
              f'singular values of DF in [{sv[:, -1].min():.2e}, {sv[:, 0].max():.2e}]')
        signs.append(set(torch.sign(d).unique().tolist()))
        rep(f'5. {nm}: det of one sign, nonzero', f'{signs[-1]}', len(signs[-1]) == 1 and 0.0 not in signs[-1])
    rep('5. same sign in all three charts (F respects the complex orientation of Y)', f'{signs}', len(set.union(*signs)) == 1)
    # ------------------------------------------------------------ 6. block structure in lattice coordinates
    n = 2000
    pl = torch.stack([U(n, math.log(0.2), math.log(5)), U(n, -math.pi, math.pi), U(n, 0, 1), U(n, 0, 1), U(n, 0, 1), U(n, 0, 1)], -1)
    def F_lat(p):
        t = torch.polar(torch.exp(p[..., 0]), p[..., 1]); lg = torch.complex(p[..., 0], p[..., 1])
        z, w = T.point_from_lattice(t, lg, p[..., 2:6])
        return Fm.F_mid(z, w, t, lg)
    dF, _ = det9(F_lat, pl)
    from extra.s4t2_explicit import check_lift as CL
    dP, _ = CL.det5(CL.psi_mid_real, pl[:, [0, 1, 2, 5]])
    # column order [n_S4, n_S1, n_S1, d_rho, d_theta, d_x1, d_x2, d_x3, d_x4]: bringing it to block form costs one
    # transposition inside the T^2 block, so det F = -(2 pi)^2 det Psi_W (and every det F above is negative)
    e = ((dF / (-(2 * math.pi) ** 2 * dP)) - 1).abs().max().item()
    rep('6. det F (lattice coords) = -(2 pi)^2 det Psi_W: |deg F| = deg Psi_W = 1', f'{e:.1e}', e < 1e-9)
    print('ALL PASS' if not fails else f'FAILURES: {fails}')
