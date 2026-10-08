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


def check_fixed_point_balls():
    """Exercise exact poles and nearby axes, which uniform target sampling misses."""
    from extra.s4t2_explicit import tlift as L
    groups = []
    for sign, scale in ((-1, L.S0), (1, L.SINF)):
        radii = torch.tensor([0., 0., 0., 0., L.EPS1 / 4, L.EPS1 / 4, L.EPS1 / 4, L.EPS1 / 4])
        Z = torch.eye(4).repeat(2, 1) * torch.sqrt(scale * radii[:, None])
        height = sign * torch.sqrt(1 - scale * radii)
        phi = torch.tensor([[0., 0.], [0.25, -0.25], [-0.5, 0.125], [0.5, -0.125]]).repeat(2, 1)
        groups.append(torch.cat([Z, height[:, None], Fm.torus_embed(phi)], -1))
    q = torch.cat(groups)
    J, DJ, res = JR.J_round(q)
    assert torch.isfinite(J).all() and torch.isfinite(DJ).all() and torch.isfinite(res).all()
    P = NR.S4T2.projector(q)
    assert (J @ J + P).abs().max() < 1e-11
    assert (P @ J @ P - J).abs().max() < 1e-11
    assert res.max() < 1e-12
    # Independent symbolic pole values (C8f-g), in the ambient tangent embedding.
    pole_indices = torch.tensor([0, 1, 2, 3, 8, 9, 10, 11])
    expected = torch.zeros(len(pole_indices), 9, 9)
    expected[:4, 0, 2] = expected[:4, 1, 3] = -1
    expected[:4, 2, 0] = expected[:4, 3, 1] = 1
    expected[4:, 0, 3] = expected[4:, 2, 1] = 1
    expected[4:, 1, 2] = expected[4:, 3, 0] = -1
    qp = q[pole_indices]
    e2, e3 = torch.zeros_like(qp), torch.zeros_like(qp)
    e2[:, 5], e2[:, 6] = -qp[:, 6], qp[:, 5]
    e3[:, 7], e3[:, 8] = -qp[:, 8], qp[:, 7]
    sign = -qp[:, 4]
    expected += sign[:, None, None] * (e3[:, :, None] * e2[:, None, :] - e2[:, :, None] * e3[:, None, :])
    assert (J[pole_indices] - expected).abs().max() < 1e-12
    N, scale = JR.nijenhuis(q, J, DJ)
    assert (N / scale).max() < 1e-10
    DJf = fd_DJ(lambda qq: JR.J_round(qq, with_DJ=False)[0], q, 1e-6)
    fd_error = ((DJf - DJ).flatten(1).norm(dim=1) / DJ.flatten(1).norm(dim=1).clamp_min(1)).max()
    assert fd_error < 1e-6
    print(f'Fixed-point balls: {len(q)} exact-pole/axis points; residual {res.max():.1e}, '
          f'|J^2+P| {(J @ J + P).abs().max():.1e}, |N|/scale {(N / scale).max():.1e}, '
          f'FD derivative error {fd_error:.1e}: PASS', flush=True)


def check_symbolic_linear_balls():
    """Independent coframe expression, compared with the established forward-map path."""
    from extra.s4t2_explicit import tlift as L
    gen = torch.Generator().manual_seed(91)
    for key, sign, scale, a0 in (('c0',-1,L.S0,1),('ci',1,L.SINF,math.exp(-4))):
        n = 256
        radii = torch.rand(n,generator=gen)*L.EPS1
        radii[:4], radii[4:8], radii[8:12], radii[12:16] = 0, 1e-20, L.EPS1/2, L.EPS1
        normal = torch.randn(n,4,generator=gen)
        normal[:16] = torch.eye(4).repeat(4,1)
        normal = normal/normal.norm(dim=-1,keepdim=True)*radii.sqrt()[:,None]
        q1,q2 = torch.complex(normal[:,0],normal[:,1]), torch.complex(normal[:,2],normal[:,3])
        Z1 = math.sqrt(scale/2)*(q1+q2)
        Z2 = math.sqrt(scale/2)*(q1-q2)*(1j if key == 'c0' else 1)
        Z = torch.stack([Z1.real,Z1.imag,Z2.real,Z2.imag],-1)
        phi = torch.rand(n,2,generator=gen)-0.5
        phi[:4] = torch.tensor([[0.,0.],[0.5,0.125],[-0.5,-0.125],[0.25,-0.25]])
        q = torch.cat([Z,sign*torch.sqrt(1-scale*radii)[:,None],Fm.torus_embed(phi)],-1)
        J = JR.J_linear_ball(q,key)
        reference = JR.J_round(q,with_DJ=False)[0]
        relative = ((J-reference).flatten(1).norm(dim=1)/reference.flatten(1).norm(dim=1)).max()
        P = NR.S4T2.projector(q)
        assert torch.isfinite(J).all() and relative < 1e-11
        assert (J@J+P).abs().max() < 1e-11 and (P@J@P-J).abs().max() < 1e-11
        M = JR.linear_ball_coframe(q,key)
        t = q1.conj()*q2/math.sqrt(a0)
        determinant = sign*(1-T.abs2(t)**2)/(scale**2*a0)
        det_error = ((torch.linalg.det(M)-determinant)/determinant).abs().max()
        assert det_error < 1e-12
        # Interior points, including the exact poles and both source-chart axes.
        # Finite differences of the new formula are independent of D^2 F.
        interior = q[:12]
        Jr,DJr,_ = JR.J_round(interior)
        DJf = fd_DJ(lambda qq: JR.J_linear_ball(qq,key),interior,1e-6)
        derivative_error = ((DJf-DJr).flatten(1).norm(dim=1)/DJr.flatten(1).norm(dim=1).clamp_min(1)).max()
        N,ns = JR.nijenhuis(interior,Jr,DJf)
        assert derivative_error < 1e-6 and (N/ns).max() < 1e-6
        print(f'Symbolic cap {key}: {n} points; J relative error {relative:.1e}, determinant relative error {det_error:.1e}, '
              f'12-point FD DJ error {derivative_error:.1e}, normalized N {(N/ns).max():.1e}: PASS',flush=True)


def check_elementary_caps():
    """New cap tensor vs direct elementary G, including pole/axis/torus seams."""
    from extra.s4t2_explicit import tlift as L, tinverse as TI
    from torch.func import jvp
    gen = torch.Generator().manual_seed(710)
    for key,sign,scale,a0 in (('c0',-1,L.S0,1),('ci',1,L.SINF,math.exp(-4))):
        n = 256
        old_r = .01*torch.rand(n,generator=gen)
        old_r[:16] = torch.tensor([0.,1e-20,.005,.01]).repeat_interleave(4)
        w = torch.randn(n,4,generator=gen)
        w[:16] = torch.eye(4).repeat(4,1)
        w *= torch.sqrt(old_r/scale)[:,None]/w.norm(dim=-1,keepdim=True)
        q1,q2 = torch.complex(w[:,0],w[:,1]),torch.complex(w[:,2],w[:,3])
        Z1 = math.sqrt(scale/2)*(q1+q2)
        Z2 = math.sqrt(scale/2)*(q1-q2)*(1j if key=='c0' else 1)
        Z = torch.stack([Z1,Z2],-1)*torch.sqrt(2/(1+old_r**2))[:,None]
        X = sign*(1-old_r)/torch.sqrt(1+old_r**2)
        xy = TI.algebraic_cap_slice(Z,X,infinity=key=='ci')
        slice_error = max((xy[0]-q1.conj()/math.sqrt(a0)).abs().max().item(),(xy[1]-q2).abs().max().item())
        assert slice_error < 1e-14
        phi = torch.rand(n,2,generator=gen)-.5
        phi[:4] = torch.tensor([[0.,0.],[.5,.125],[-.5,-.125],[.25,-.25]])
        zr = torch.stack([Z[:,0].real,Z[:,0].imag,Z[:,1].real,Z[:,1].imag],-1)
        u = torch.cat([zr,phi],-1)
        q = torch.cat([zr,X[:,None],Fm.torus_embed(phi)],-1)
        M = JR.elementary_cap_coframe(q,key)
        D = torch.stack([jvp(lambda v:JR.elementary_cap_map(v,key),(u,),
                             (torch.eye(6)[j].expand_as(u),))[1] for j in range(6)],-1)
        reference = torch.linalg.inv(D)@NR.J0@D
        J6 = torch.linalg.inv(M)@NR.J0@M
        tensor_error = (J6-reference).abs().max().item()
        assert tensor_error < 1e-11
        # This independently tests all coframe columns at phi=0.
        uzero = torch.cat([zr,torch.zeros_like(phi)],-1)
        Dzero = torch.stack([jvp(lambda v:JR.elementary_cap_map(v,key),(uzero,),
                             (torch.eye(6)[j].expand_as(u),))[1] for j in range(6)],-1)
        assert (M-Dzero).abs().max() < 1e-12
        r = (zr*zr).sum(-1)
        s = torch.sqrt(1-r*r)
        gain = torch.rsqrt(1+s)
        t = xy[0]*xy[1]
        expected_det = sign*(1-T.abs2(t)**2)*gain**4/(scale**2*a0*s)
        det_error = ((torch.linalg.det(M)-expected_det)/expected_det).abs().max().item()
        assert det_error < 1e-12
        J = JR.J_elementary_cap(q,key)
        P = NR.S4T2.projector(q)
        assert (J@J+P).abs().max() < 1e-11 and (P@J@P-J).abs().max() < 1e-11
        # Integer torus changes are holomorphic transitions, hence the same J6.
        shifted = u.clone(); shifted[:,4:6] += torch.tensor([1.,-2.])
        Ds = torch.stack([jvp(lambda v:JR.elementary_cap_map(v,key),(shifted,),
                            (torch.eye(6)[j].expand_as(u),))[1] for j in range(6)],-1)
        assert (torch.linalg.inv(Ds)@NR.J0@Ds-J6).abs().max() < 1e-11
        # Finite-difference N is independent of the coframe's symbolic derivation.
        interior = torch.cat([q[:12],q[16:48]],0)
        Ji = JR.J_elementary_cap(interior,key)
        DJ = fd_DJ(lambda qq:JR.J_elementary_cap(qq,key),interior,1e-6)
        N,ns = JR.nijenhuis(interior,Ji,DJ)
        nij_error = (N/ns.clamp_min(1)).max().item()
        assert nij_error < 1e-7
        perturbation = torch.randn(9,9,generator=gen)
        def nonintegrable(qq):
            projector = NR.S4T2.projector(qq)
            R = projector+.1*torch.sin(3*qq[:,0])[:,None,None]*(projector@perturbation@projector)
            return R@JR.J_elementary_cap(qq,key)@torch.linalg.pinv(R)
        control_points = q[16:24]
        control_J = nonintegrable(control_points)
        control_N,control_scale = JR.nijenhuis(control_points,control_J,fd_DJ(nonintegrable,control_points,1e-6))
        control_error = (control_N/control_scale.clamp_min(1)).median().item()
        assert control_error > 1e-5
        print(f'Elementary cap {key}: {n} points; slice {slice_error:.1e}, direct-G J {tensor_error:.1e}, '
              f'determinant {det_error:.1e}, {len(interior)}-point FD normalized N {nij_error:.1e}, '
              f'nonintegrable control {control_error:.1e}: PASS',flush=True)


if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    t0 = time.time()
    torch.manual_seed(17)
    check_fixed_point_balls()
    check_symbolic_linear_balls()
    check_elementary_caps()
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
