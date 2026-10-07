"""J_round as a black-box field on the round S^4 x T^2 (jround.py), tested the way an independent user would:
  1. derivatives: DJ from nested AD vs central finite differences of q -> J_round(q) along an orthonormal tangent
     frame (each evaluation goes through F^-1, so neighbouring points may be evaluated in different charts);
     points: uniform on S^4 x T^2, plus points on the chart switches |t| = 0.2 and |t| = 5;
  2. N with the ambient code from the finite-difference DJ (no D^2 F used at all);
  3. control: a non-integrable J_eps = R J_round R^+ (R = P (1 + eps sin(3 q_0) E) P), derivatives also by finite
     differences: N = O(eps);
  4. T^3-invariance: J_round(g q) = dg J_round(q) dg^T for g = (e^{i theta} on Z, rotations of both circles);
  5. the infinity chart: J and DJ via the filling chart at infinity and via the middle chart agree.
Run from the repo root: python3 extra/s4t2_explicit/check_jround.py [n]"""
from __future__ import annotations
import math, os, sys, time
import torch
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
torch.set_default_dtype(torch.float64)
from extra.s4t2_explicit import jround as JR, nijenhuis_round as NR, check_full as CF, tfull as Fm, tcore as T

fails = []


def rep(name, val, ok):
    if not ok:
        fails.append(name)
    print(f'{name:<78s} {val}  {"ok" if ok else "FAIL"}')


def retract(q):
    return torch.cat([q[:, :5] / q[:, :5].norm(dim=1, keepdim=True), q[:, 5:7] / q[:, 5:7].norm(dim=1, keepdim=True),
                      q[:, 7:] / q[:, 7:].norm(dim=1, keepdim=True)], 1)


def tangent_frame(q):
    """(B, 6, 9) orthonormal frame of T_q(S^4 x S^1 x S^1)."""
    B = q.shape[0]
    P1 = torch.eye(5) - q[:, :5, None] * q[:, None, :5]
    _, V = torch.linalg.eigh(P1)
    S4 = V[:, :, 1:].transpose(1, 2)                                   # eigenvalue 1: the 4 tangent directions
    fr = torch.zeros(B, 6, 9)
    fr[:, :4, :5] = S4
    fr[:, 4, 5], fr[:, 4, 6] = -q[:, 6], q[:, 5]
    fr[:, 5, 7], fr[:, 5, 8] = -q[:, 8], q[:, 7]
    return fr


def fd_DJ(Jfun, q, h):
    """DJ[b, a, r, c] = sum_k (D_{f_k} J)[a, r] f_k[c] by central differences along the frame (retracted)."""
    fr = tangent_frame(q)
    DJ = torch.zeros(q.shape[0], 9, 9, 9)
    for k in range(6):
        Jp = Jfun(retract(q + h * fr[:, k]))
        Jm = Jfun(retract(q - h * fr[:, k]))
        DJ = DJ + ((Jp - Jm) / (2 * h))[..., None] * fr[:, k][:, None, None, :]
    return DJ


def switch_points(n, radius):
    gen = CF.gen
    t = torch.polar(torch.full((n,), radius), torch.empty(n).uniform_(-math.pi, math.pi, generator=gen))
    x = torch.empty(n, 4).uniform_(0, 1, generator=gen)
    z, w = T.point_from_lattice(t, torch.log(t), x)
    return Fm.F_mid(z, w, t, torch.log(t))


if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    t0 = time.time()
    torch.manual_seed(17)
    Jonly = lambda qq: JR.J_round(qq, with_DJ=False)[0]
    groups = {'uniform on S^4 x T^2': NR.S4T2.sample(n),
              'on the chart switch |t| = 0.2': switch_points(n // 3, 0.2),
              'on the chart switch |t| = 5': switch_points(n // 3, 5.0)}
    for nm, q in groups.items():
        J, DJ, res = JR.J_round(q)
        P = NR.S4T2.projector(q)
        print(f'{nm}: {len(q)} points, |F(F^-1 q) - q| max {res.max():.1e}, |J^2 + P| max {(J @ J + P).flatten(1).norm(dim=1).max():.1e}')
        for h in (1e-5, 1e-6):
            DJf = fd_DJ(Jonly, q, h)
            err = (DJf - DJ).flatten(1).norm(dim=1) / DJ.flatten(1).norm(dim=1)
            Nn, scale = JR.nijenhuis(q, J, DJf)
            print(f'   h = {h:.0e}: |DJ_fd - DJ_AD| / |DJ| max {err.max():.1e} (median {err.median():.1e});  '
                  f'N from DJ_fd: |N| / scale max {(Nn / scale).max():.1e}')
        rep(f'1-2. {nm}: finite-difference DJ = AD DJ, and N(J_round) from it ~ 0', f'{err.max():.1e}, {(Nn / scale).max():.1e}',
            err.max() < 1e-4 and (Nn / scale).max() < 1e-5)
        Nn, scale = JR.nijenhuis(q, J, DJ)
        rep(f'     {nm}: N from the exact DJ', f'{(Nn / scale).max():.1e}', (Nn / scale).max() < 1e-8)
    # 3. control: non-integrable J_eps through the same black box, finite differences only
    q = NR.S4T2.sample(n)
    E = torch.randn(9, 9, generator=torch.Generator().manual_seed(3))
    for eps in (1e-1, 1e-3):
        def Jeps(qq):
            P = NR.S4T2.projector(qq)
            R = P + eps * torch.sin(3 * qq[:, 0])[:, None, None] * (P @ E @ P)
            return R @ Jonly(qq) @ torch.linalg.pinv(R)
        Je = Jeps(q)
        Nn, scale = JR.nijenhuis(q, Je, fd_DJ(Jeps, q, 1e-5))
        print(f'3. control J_eps, eps = {eps:.0e}: |J^2 + P| {(Je @ Je + NR.S4T2.projector(q)).flatten(1).norm(dim=1).max():.1e}; '
              f'|N| / scale median {(Nn / scale).median():.1e}, max {(Nn / scale).max():.1e}')
    # 4. T^3-invariance
    q = NR.S4T2.sample(n)
    th, a, b = (torch.empty(n).uniform_(-math.pi, math.pi) for _ in range(3))
    def rot(ang):
        c, s = torch.cos(ang), torch.sin(ang)
        return torch.stack([torch.stack([c, -s], -1), torch.stack([s, c], -1)], -2)
    dg = torch.zeros(n, 9, 9)
    dg[:, 0:2, 0:2], dg[:, 2:4, 2:4], dg[:, 4, 4] = rot(th), rot(th), 1.0
    dg[:, 5:7, 5:7], dg[:, 7:9, 7:9] = rot(a), rot(b)
    gq = torch.einsum('bij,bj->bi', dg, q)
    J1, J2 = Jonly(q), Jonly(gq)
    e = ((J2 - dg @ J1 @ dg.transpose(1, 2)).flatten(1).norm(dim=1) / J1.flatten(1).norm(dim=1)).max().item()
    rep('4. T^3-invariance: J(g q) = dg J(q) dg^T', f'{e:.1e}', e < 1e-10)
    # 5. infinity chart vs middle chart. Where they hand over (|t| in [4, 5]: F^-1 switches at |t| = 5) they must agree
    #    to round-off; for larger |t| the middle chart (z, w, t) is badly scaled (w ~ t^2 w', cond(DF) up to 6e4,
    #    versus <= 2e2 in the infinity chart), so the agreement degrades with cond (reported, not tested).
    for lo, hi, tested in ((0.2, 0.25, True), (0.05, 0.2, False)):
        x, y, w = CF.chart_sample(n, T.CUSPINF, rmax=hi, rmin=lo)
        q1, A1, H1 = NR.DF_and_D2F(Fm.F_ci_real, CF.real6(x, y, w))
        zm, wm, tm, _ = Fm.ci_to_mid(x, y, w)
        q2, A2, H2 = NR.DF_and_D2F(Fm.F_mid_real, CF.real6(zm, wm, tm))
        Ja, DJa = NR.transport(A1, H1, q1)
        Jb, DJb = NR.transport(A2, H2, q2)
        eJ = ((Ja - Jb).flatten(1).norm(dim=1) / Ja.flatten(1).norm(dim=1)).max().item()
        eD = ((DJa - DJb).flatten(1).norm(dim=1) / DJa.flatten(1).norm(dim=1)).max().item()
        val = f'{(q1 - q2).abs().max():.1e}, {eJ:.1e}, {eD:.1e}'
        if tested:
            rep('5. infinity vs middle chart at the hand-over |t| in [4, 5]: q, J, DJ', val, eJ < 1e-10 and eD < 1e-8)
        else:
            print(f'   (|t| in [5, 20], middle chart outside its range: q, J, DJ {val})')
    print(f'{time.time() - t0:.0f} s')
    print('ALL PASS' if not fails else f'FAILURES: {fails}')
