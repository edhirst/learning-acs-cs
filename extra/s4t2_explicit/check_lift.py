"""Checks of the lift Psi_W: W -> S^4 (tlift.Psi_W / Psi_chart0 / Psi_chartinf), docs/s4t2_explicit_diffeo.md 12.5.
  1. Psi lands in S^4, is S^1-equivariant (Z -> e^{i theta} Z, X fixed), well defined on W, covers Lambda (Hopf);
  2. the explicit linear models agree with the general formula where both apply (seamless switch at rad = EPS1);
  3. smoothness: D Psi (exact, autograd) is continuous across the gluing loci (shell edges, Lambda^-1(M1) = the
     z1 = 0 sphere, Gamma_E, the glue annuli of G_A, the representation switches) - pairs of points at distance h;
  4. local diffeomorphism: det[Psi, D Psi] (5 x 5) has constant sign, bounded away from 0, on random samples of
     each chart, plus an adversarial minimisation of |det| / scale;
  5. the margin of the shell gluing: max |Arg g| over shell samples (must be < pi).
Run from the repo root."""
from __future__ import annotations
import math, os, sys
import torch
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
torch.set_default_dtype(torch.float64)
from torch.func import jvp
from extra.s4t2_explicit import tcore as T, tlift as L

gen = torch.Generator().manual_seed(77)
U = lambda n, a, b: torch.empty(n).uniform_(a, b, generator=gen)
fails = []


def rep(name, val, ok):
    if not ok:
        fails.append(name)
    print(f'{name:<70s} {val}  {"ok" if ok else "FAIL"}')


def as5(Z1, Z2, X):
    return torch.stack([Z1.real, Z1.imag, Z2.real, Z2.imag, X], -1)


# ---------------------------------------------------------------- representations as real coordinates
def psi_mid_real(p):
    """p = (log r, arg t, x1, x4) -> R^5."""
    t = torch.polar(torch.exp(p[..., 0]), p[..., 1])
    lg = torch.complex(p[..., 0], p[..., 1])
    return as5(*L.Psi_W(t, p[..., 2], p[..., 3], lg))


def psi_c0_real(p):
    return as5(*L.Psi_chart0(torch.complex(p[..., 0], p[..., 1]), torch.complex(p[..., 2], p[..., 3])))


def psi_ci_real(p):
    return as5(*L.Psi_chartinf(torch.complex(p[..., 0], p[..., 1]), torch.complex(p[..., 2], p[..., 3])))


def jac(f, p):
    """(B, 5, 4) Jacobian by forward AD."""
    cols = [jvp(f, (p,), (torch.zeros_like(p).index_fill_(-1, torch.tensor([k]), 1.0),))[1] for k in range(p.shape[-1])]
    return torch.stack(cols, -1)


def det5(f, p):
    J = jac(f, p)
    q = f(p)
    M = torch.cat([q[..., None], J], -1)
    return torch.linalg.det(M), J


def mid_points(n, lo=0.2, hi=5.0):
    return torch.stack([U(n, math.log(lo), math.log(hi)), U(n, -math.pi, math.pi), U(n, 0, 1), U(n, 0, 1)], -1)


def chart_points(n, C, rmax=0.25):
    rr = torch.exp(U(n, math.log(1e-6), math.log(rmax)))
    ratio = torch.exp(U(n, -3.0, 3.0)) * math.sqrt(C.a0)
    ax, ay = U(n, -math.pi, math.pi), U(n, -math.pi, math.pi)
    x = torch.polar(torch.sqrt(rr / ratio), ax)
    y = torch.polar(torch.sqrt(rr * ratio), ay)
    k = C.mu_chart(x, y).abs() < 0.48
    return torch.stack([x.real, x.imag, y.real, y.imag], -1)[k]


if __name__ == '__main__':
    # ------------------------------------------------------------ 1. basic properties
    n = 20000
    lr = U(n, math.log(1e-3), math.log(1e3)); th = U(n, -math.pi, math.pi); x1 = U(n, 0, 1); x4 = U(n, 0, 1)
    t = torch.polar(torch.exp(lr), th); lg = torch.complex(lr, th)
    Z1, Z2, X = L.Psi_W(t, x1, x4, lg)
    rep('1. |Z|^2 + X^2 = 1', f'{(Z1.abs() ** 2 + Z2.abs() ** 2 + X ** 2 - 1).abs().max():.1e}',
        (Z1.abs() ** 2 + Z2.abs() ** 2 + X ** 2 - 1).abs().max() < 1e-13)
    v, Xl = L.Lam_W(t, x1, x4, lg)
    h = 2 * Z1 * Z2.conj()
    hv = torch.stack([h.real, h.imag, Z1.abs() ** 2 - Z2.abs() ** 2], -1)
    rep('   Hopf(Z) = v and X = X of Lambda (orbit map is Lambda)', f'{max((hv - v).abs().max(), (X - Xl).abs().max()):.1e}',
        max((hv - v).abs().max(), (X - Xl).abs().max()) < 1e-13)
    a1, a2, a3 = L.Psi_W(t, x1 + 0.37, x4, lg)
    ph = torch.polar(torch.ones(n), torch.full((n,), 2 * math.pi * 0.37))
    e = max((a1 - ph * Z1).abs().max(), (a2 - ph * Z2).abs().max(), (a3 - X).abs().max())
    rep('   S^1-equivariance (x1 + 0.37)', f'{e:.1e}', e < 1e-12)
    b1, b2, b3 = L.Psi_W(t, x1 - x4, x4, lg + 2j * math.pi)
    c1, c2, c3 = L.Psi_W(t, x1, x4 + 1, lg)
    e = max((b1 - Z1).abs().max(), (b2 - Z2).abs().max(), (c1 - Z1).abs().max(), (c2 - Z2).abs().max())
    rep('   well defined on W (branch change, x4 + 1)', f'{e:.1e}', e < 1e-11)
    # ------------------------------------------------------------ 2. linear models vs general formula
    for nm, C, lam, Ef, Gf, lin, inv in (('p0  ', T.CUSP0, L.Lam_chart0, L.E_chart0, L.GA_hat_chart0, L.linear_model0, False),
                                         ('pinf', T.CUSPINF, L.Lam_chartinf, L.E_chartinf, L.GA_hat_chartinf, L.linear_modelinf, True)):
        p = chart_points(200000, C, rmax=0.01)
        x, y = torch.complex(p[:, 0], p[:, 1]), torch.complex(p[:, 2], p[:, 3])
        rad = C.a0 * x.abs() ** 2 + y.abs() ** 2
        k = (rad > 0.5 * L.EPS1) & (rad <= L.EPS1)
        x, y = x[k], y[k]
        v, X = lam(x, y)
        rr = (x * y).abs()
        G = L.assemble(v, X, Ef(x, y), Gf(x, y), 1 / rr if inv else rr)
        Lm = lin(x, y)
        e = max((G[0] - Lm[0]).abs().max(), (G[1] - Lm[1]).abs().max(), (G[2] - Lm[2]).abs().max())
        rep(f'2. {nm}: general formula = linear model on EPS1/2 < rad <= EPS1 ({len(x)} pts)', f'{e:.1e}', e < 1e-12)
    # ------------------------------------------------------------ 3. continuity of D Psi across the gluing loci
    # jump(h) = |DPsi(p + h d) - DPsi(p - h d)| / |DPsi(p)|: ~ 2 h |D^2 Psi| for a C^2 map, h-independent at a kink.
    # Test: jump(h)/h stays bounded as h -> 0 (ratio of the maxima over h = 1e-4, 1e-6 close to 1).
    for nm, f, pts, loc in (('mid', psi_mid_real, mid_points(60000, 0.03, 30.0), lambda p: f'|t| = {p[0].exp():.4f}'),
                            ('c0', psi_c0_real, chart_points(60000, T.CUSP0), lambda p: f'rad = {(p ** 2).sum():.4f}'),
                            ('ci', psi_ci_real, chart_points(60000, T.CUSPINF), lambda p: f'rad = {T.CUSPINF.a0 * (p[:2] ** 2).sum() + (p[2:] ** 2).sum():.4f}')):
        dirs = torch.randn(pts.shape, generator=gen)
        dirs = dirs / dirs.norm(dim=-1, keepdim=True)
        Jc = jac(f, pts)
        scale = Jc.flatten(1).norm(dim=-1)
        q = {}
        for hh in (1e-4, 1e-5, 1e-6):
            jump = (jac(f, pts + hh * dirs) - jac(f, pts - hh * dirs)).flatten(1).norm(dim=-1) / scale
            q[hh] = jump / hh
        i = int(q[1e-6].argmax())
        print(f'   3. {nm}: max jump/h at h = 1e-4, 1e-5, 1e-6: {q[1e-4].max():.3e}, {q[1e-5].max():.3e}, {q[1e-6].max():.3e} '
              f'(worst point: {loc(pts[i])})')
        ratio = (q[1e-6].max() / q[1e-4].max()).item()
        rep(f'3. {nm}: D Psi is C^1 (jump/h bounded as h -> 0)', f'ratio {ratio:.3f}', 0.5 < ratio < 2.0)
    # ------------------------------------------------------------ 4. Jacobian determinant
    for nm, f, pts in (('mid (0.2 <= |t| <= 5)', psi_mid_real, mid_points(100000)),
                       ('c0 chart', psi_c0_real, chart_points(100000, T.CUSP0)),
                       ('ci chart', psi_ci_real, chart_points(100000, T.CUSPINF))):
        d, J = det5(f, pts)
        sv = torch.linalg.svdvals(J)
        rel = d.abs() / sv.prod(-1).clamp_min(1e-300)
        print(f'   4. {nm}: det sign +{int((d > 0).sum())} / -{int((d < 0).sum())}; min |det| {d.abs().min():.3e}; '
              f'min sigma_4 {sv[:, -1].min():.3e}; max sigma_1 {sv[:, 0].max():.3e}')
        rep(f'4. {nm}: det of constant sign, nonzero', f'{int((d > 0).sum())}/{int((d < 0).sum())}',
            bool((d > 0).all() or (d < 0).all()))
        # adversarial: projected gradient descent on log|det| from the 64 smallest samples
        idx = d.abs().argsort()[:64]
        q = pts[idx].clone().requires_grad_(True)
        opt = torch.optim.Adam([q], lr=2e-3)
        for _ in range(150):
            opt.zero_grad()
            dd, _ = det5(f, q)
            loss = torch.log(dd.abs()).sum()
            loss.backward()
            opt.step()
            with torch.no_grad():
                if nm.startswith('mid'):
                    q[:, 0].clamp_(math.log(0.2), math.log(5.0))
                else:
                    C = T.CUSP0 if nm.startswith('c0') else T.CUSPINF
                    x, y = torch.complex(q[:, 0], q[:, 1]), torch.complex(q[:, 2], q[:, 3])
                    bad = (C.mu_chart(x, y).abs() > 0.48) | ((x * y).abs() > 0.25)
                    q[bad] = pts[idx][bad]
        dd, _ = det5(f, q.detach())
        rep(f'4. {nm}: adversarial min |det| (64 starts x 150 Adam steps)', f'{dd.abs().min():.3e} (sign kept: {bool((torch.sign(dd) == torch.sign(d[idx])).all())})',
            bool((dd.abs().min() > 1e-3) and (torch.sign(dd) == torch.sign(d[idx])).all()))
    # ------------------------------------------------------------ 5. margin of the shell gluing
    from extra.s4t2_explicit import fit_psi as FP
    gg, lr = FP.collect(batches=4)
    dev = torch.angle(gg * torch.polar(torch.ones(1), -L.psi_of_r(torch.exp(lr)))).abs()
    rep(f'5. max |Arg(g e^-i psi)| on {len(gg)} shell samples (< pi)', f'{dev.max():.3f} rad', dev.max() < math.pi - 0.5)
    # adversarial ascent of |Arg| inside the shell (mid representation, the only place with large deviations)
    pts = mid_points(400000, 0.25, 4.0)
    def shell_quant(p):
        t = torch.polar(torch.exp(p[..., 0]), p[..., 1]); lg = torch.complex(p[..., 0], p[..., 1])
        v, X = L.Lam_W(t, p[..., 2], p[..., 3], lg)
        thv = torch.acos((-v[..., 2] / v.norm(dim=-1)).clamp(-1, 1))
        gq = L.g_shell(v, X, L.E_W(t, p[..., 2], p[..., 3], lg), L.GA_hat_W(t, p[..., 2], p[..., 3], lg), t.abs())
        return thv, torch.angle(gq).abs()
    with torch.no_grad():
        thv, ang = shell_quant(pts)
    insh = (thv > L.THETA_A) & (thv < L.THETA_B)
    idx = torch.where(insh, ang, torch.zeros_like(ang)).argsort(descending=True)[:64]
    q = pts[idx].clone().requires_grad_(True)
    opt = torch.optim.Adam([q], lr=2e-3)
    best = ang[idx].max().item()
    for _ in range(200):
        opt.zero_grad()
        thv, ang = shell_quant(q)
        pen = torch.relu(L.THETA_A - thv) + torch.relu(thv - L.THETA_B)
        (-(ang - 50 * pen).sum()).backward()
        opt.step()
        with torch.no_grad():
            thv, ang = shell_quant(q)
            ok = (thv >= L.THETA_A) & (thv <= L.THETA_B)
            if ok.any():
                best = max(best, ang[ok].max().item())
    rep('5. adversarial max |Arg(g e^-i psi)| in the shell (64 starts x 200 steps)', f'{best:.3f} rad', best < math.pi - 0.5)
    print('ALL PASS' if not fails else f'FAILURES: {fails}')
