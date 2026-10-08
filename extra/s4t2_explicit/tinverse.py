"""F^-1: S^4 x T^2 -> Y (docs/s4t2_explicit_diffeo.md, section 12.7), so that J_round can be evaluated at ANY point
of the round S^4 x T^2 (e.g. uniform samples, supervised targets).

Explicit chain (inverting Lambda = P_final o moebius o R o h_B step by step), then a Newton polish on F itself:
  1. v = Hopf(Z1, Z2); w = P_final^-1(v, X)  (radial: invert nu, monotone; direction rho_O^-1);
  2. u' = inverse stereographic projection from (-1, 0) of w + (0, 1, 0);
  3. r from |u'1| / |u'2| = |T(r)| = r k(r)  (monotone, Newton in log r); u = R^-1 u' (R depends on r only);
  4. t = r e^{i arg(u1/u2)}; fibre coordinate a1 = arg(u2)/2 pi (r <= 1) or a2' = arg(u1)/2 pi (r >= 1);
  5. W-point: middle  x4 = a1 + xi4(t) + 1/2,  resp. a2 = -a2' + lam delta + (1 - lam)/2, x4 = -a2 - 1/2;
              cusps   log|zeta| from a1_c0 (resp. a2_ci) by 1D Newton (monotone along the fibre);
     inside the linear balls the lift is linear: q = U^-1 Z / sqrt(s) directly;
  6. S^1 phase from Z / Psi(point with x1 = 0); T^2 part: (x2, x3) = phi + corrections (resp. C = phi + H c);
  7. Newton polish p <- p + DF(p)^+ (q - F(p)) in the output chart (mid for 0.2 <= |t| <= 5, else the fillings).
"""
from __future__ import annotations

import math

import torch

from extra.s4t2_explicit import tcore as T, tlift as L, tfull as Fm
from extra.s4t2_explicit.tcore import CD, RD, bump_step
from torch.func import jvp

EYE6 = torch.eye(6)


def transport_shell_inverse(v, infinity=False, steps=96,*,algebraic_profile=False):
    """Contraction inverse for the extended local base map (programme 8.1).

    Return p and |p|<=EPS2. A false flag means use the raw global inverse.
    The exact expression is the limit; its certified contraction ratio is <.68
    for the original profile and <.7 for the algebraic profile (programme 9.2).
    This independent function does not alter Finv or its Newton path.
    """
    scale = L.SINF if infinity else L.S0
    M = scale*(L.QINF if infinity else L.Q0).to(v)
    seed = v@M/scale**2
    p = seed
    for _ in range(steps):
        radius = torch.linalg.vector_norm(p,dim=-1)
        clipped = p*(L.EPS2/radius.clamp_min(L.EPS2))[...,None]
        defect = L.transport_base_local(clipped,infinity,algebraic_profile=algebraic_profile)-clipped@M.T
        p = seed-defect@M/scale**2
    return p, torch.linalg.vector_norm(p,dim=-1) <= L.EPS2


def _transport_axis_radius(a,finite_height=False):
    """Central-fibre radius, including the node correction, for |a|<=1/2.

    Solve a=sign(a)[c(r)r+(1-c(r))r/(1+r)], c=1 on r<=.02,
    c=0 on r>=.05. Strict monotonicity is I1; endpoints r=0,1 bracket.
    """
    lo,hi = torch.zeros_like(a),torch.ones_like(a)
    for _ in range(64):
        r = (lo+hi)/2
        c = 1-bump_step((r-T.NODE_IN)/(T.NODE_OUT-T.NODE_IN))
        height = T.finite_height_seed(r) if finite_height else r/(1+r)
        val = c*r+(1-c)*height
        lower = val < a.abs()
        lo,hi = torch.where(lower,r,lo),torch.where(lower,hi,r)
    return torch.where(a == 0,torch.zeros_like(a),(lo+hi)/2)


def P_algebraic_inv(v, X):
    """Inverse of L.P_algebraic away from the two poles; no root solve."""
    h = torch.linalg.vector_norm(v,dim=-1)
    a = X*torch.sqrt(2-X*X)
    # Rationalize only near the north pole, where 1-a loses precision.
    gain = torch.where(X>0,(1+a)/(h*h),1/(1-a))
    return torch.stack([v[...,0],-v[...,2],-v[...,1]],-1)*gain[...,None]


def algebraic_cap_slice(Z,X,infinity=False):
    """Closed slice inverse of section 9.2's radial pole model.

    Caller supplies a point in the appropriate cap |Z_old|^2<=1/100.
    The expression includes the exact pole and avoids division by |Z|.
    """
    h = (Z.real*Z.real+Z.imag*Z.imag).sum(-1)/(1+X.abs()*torch.sqrt(2-X*X))
    old_Z = Z*torch.sqrt((1+h*h)/2)[...,None]
    scale = L.SINF if infinity else L.S0
    second = old_Z[...,1] if infinity else -1j*old_Z[...,1]
    q1 = (old_Z[...,0]+second)/math.sqrt(2*scale)
    q2 = (old_Z[...,0]-second)/math.sqrt(2*scale)
    C = T.CUSPINF if infinity else T.CUSP0
    return q1.conj()/math.sqrt(C.a0),q2


def algebraic_r_from_absT(absT):
    """Positive quadratic root for |T|=r*algebraic_k(r), no iteration."""
    A,B = (4*math.pi)**2,math.exp(4)
    u = absT*absT
    b = A*(1-u)
    root = torch.sqrt(b*b+4*B*u)
    # Choose a noncancelling root form on each side of u=1.
    denominator = torch.where(u<=1,b+root,torch.ones_like(u))
    y = torch.where(u<=1,2*u/denominator,(root-b)/(2*B))
    return torch.sqrt(y)


def transport_base_inverse(v, X,*,algebraic_radial=False,algebraic_profile=False,finite_height=False):
    """Reference W-point over the new, unfitted base inverse (programme 8.1).

    Return {key: (indices,data)}. c0/ci data=(Re x,Im x,Re y,Im y);
    mid data=(Re t,Im t,x1=0,x4). Charts are local gauges, not a continuous
    global section. Both central fibres and both fixed points are included.
    This is a numerical evaluator of certified roots; exact definitions use
    the infinite height series and the contraction/bisection limits. With
    algebraic_radial=True the nu inverse is eliminated; algebraic_profile=True
    also eliminates the k-profile inverse (programme 9.2). finite_height=True
    eliminates the Gaussian series, but still uses a scalar bisection here.
    """
    out = {}
    assigned = torch.zeros_like(X,dtype=torch.bool)
    C0,Ci = T.height_cusps(finite_height)

    def put(key,idx,data):
        if len(idx):
            out[key] = (torch.cat([out[key][0],idx]),torch.cat([out[key][1],data])) if key in out else (idx,data)

    def xy_data(x,y):
        return torch.stack([x.real,x.imag,y.real,y.imag],-1)

    for key,C,infinity in (('c0',C0,False),('ci',Ci,True)):
        idx = ((X>0) if infinity else (X<0)).nonzero().flatten()
        if not len(idx):
            continue
        local_v = v[idx]
        if algebraic_radial:
            local_v = local_v/(1+X[idx].abs()*torch.sqrt(2-X[idx]*X[idx]))[:,None]
        p,inside = transport_shell_inverse(local_v,infinity,algebraic_profile=algebraic_profile)
        idx,p = idx[inside],p[inside]
        if not len(idx):
            continue
        radius = torch.linalg.vector_norm(p,dim=-1)
        t = torch.complex(p[:,0],-p[:,1])/(2*math.sqrt(C.a0))
        use_y = p[:,2] <= 0
        yr = torch.sqrt(((radius-p[:,2])/2).clamp_min(0))
        xr = torch.sqrt(((radius+p[:,2])/(2*C.a0)).clamp_min(0))
        # Choose the nonzero coordinate as positive, and assign the node separately.
        ys = torch.where(use_y,yr,torch.ones_like(yr)).to(CD)
        xs = torch.where(~use_y,xr,torch.ones_like(xr)).to(CD)
        x = torch.where(use_y,t/torch.where(yr>0,ys,torch.ones_like(ys)),xs)
        y = torch.where(use_y,ys,t/xs)
        put(key,idx,xy_data(x,y))
        assigned[idx] = True

    rest = (~assigned).nonzero().flatten()
    if not len(rest):
        return out
    vr,xr = v[rest],X[rest]
    if algebraic_radial:
        u = moebius_inv(P_algebraic_inv(vr,xr))
    else:
        nv = torch.linalg.vector_norm(vr,dim=-1)
        # Use X near the equator: -log|v| loses the small X^2 by subtraction.
        level = torch.where(xr.abs()<0.5,-torch.log1p(-xr*xr),-torch.log(nv))
        tau = torch.sign(xr)*nu_inv(level)
        direction = torch.stack([vr[:,0],-vr[:,2],-vr[:,1]],-1)/nv[:,None]
        u = moebius_inv(torch.exp(tau)[:,None]*direction)
    zero,infty = u[:,0].abs() == 0,u[:,1].abs() == 0
    for key,C,mask,a in (
            ('c0',C0,zero,T.wrap(torch.angle(u[:,1])/T.TWO_PI)),
            ('ci',Ci,infty,T.wrap(-torch.angle(u[:,0])/T.TWO_PI+0.5))):
        if mask.any():
            aa = a[mask]
            r = _transport_axis_radius(aa,finite_height)
            x = torch.sqrt(torch.where(aa>0,r/C.a0,torch.zeros_like(r))).to(CD)
            y = torch.sqrt(torch.where(aa<0,r,torch.zeros_like(r))).to(CD)
            put(key,rest[mask],xy_data(x,y))

    reg = ~(zero|infty)
    rest,u = rest[reg],u[reg]
    if len(rest):
        r = (algebraic_r_from_absT if algebraic_profile else r_from_absT)(u[:,0].abs()/u[:,1].abs())
        t = torch.polar(r,torch.angle(u[:,0]/u[:,1]))
        a1 = T.wrap(torch.angle(u[:,1])/T.TWO_PI)
        a2p = torch.angle(u[:,0])/T.TWO_PI
        middle = (r>=0.5)&(r<=2)
        if middle.any():
            tt,rr = t[middle],r[middle]
            a2 = -a2p[middle]+T.lam_inf(rr)*T.delta(tt)+(1-T.lam_inf(rr))/2
            x4 = torch.where(rr<=1,a1[middle]+T.xi4(tt)+0.5,-a2-0.5)
            put('mid',rest[middle],torch.stack([tt.real,tt.imag,torch.zeros_like(x4),x4],-1))
        for key,C,mask in (('c0',C0,r<0.5),('ci',Ci,r>2)):
            if not mask.any():
                continue
            tc = t[mask] if key == 'c0' else 1/t[mask]
            aa = a1[mask] if key == 'c0' else T.wrap(-a2p[mask]+0.5)
            pp = C.P(tc)
            lo,hi = (aa+2.5)*pp,(aa-1.5)*pp
            for _ in range(64):
                rho = (lo+hi)/2
                value = T.cusp_angle(C,torch.exp(rho).to(CD),tc)
                lower = value > aa
                lo,hi = torch.where(lower,rho,lo),torch.where(lower,hi,rho)
            y = torch.exp((lo+hi)/2).to(CD)
            put(key,rest[mask],xy_data(tc/y,y))
    if sum(len(idx) for idx,_ in out.values()) != len(X):
        raise RuntimeError('transport base inverse did not cover all input points')
    return out


def nu_inv(nu):
    """|tau| with nu_fn(|tau|) = nu (nu_fn is increasing in |tau|): sqrt for nu <= 1/4, identity for nu >= 1,
    Newton on [1/2, 1] in between."""
    a = torch.where(nu <= 0.25, torch.sqrt(nu.clamp_min(0)), torch.where(nu >= 1, nu, 0.75 * torch.ones_like(nu)))
    mid = (nu > 0.25) & (nu < 1)
    for _ in range(60):
        aa = a.clone().requires_grad_(True)
        f = T.nu_fn(aa) - nu
        g, = torch.autograd.grad(f.sum(), aa)
        a = torch.where(mid, (a - f.detach() / g).clamp(0.5, 1.0), a)
    return a


def P_final_inv(v, X):
    nv = torch.linalg.norm(v, dim=-1)
    tau = torch.sign(X) * nu_inv(-torch.log(nv))
    d = torch.stack([v[..., 0], -v[..., 2], -v[..., 1]], -1) / nv[..., None]      # rho_O^-1 (RHO_SIGN = -1)
    return torch.exp(tau)[..., None] * d


def moebius_inv(w):
    y = w + torch.tensor([0.0, 1.0, 0.0], dtype=RD)
    n2 = (y ** 2).sum(-1)
    P1 = (1 - n2) / (1 + n2)
    Pr = 2 * y / (1 + n2)[..., None]
    return torch.stack([torch.complex(P1, Pr[..., 0]), torch.complex(Pr[..., 1], Pr[..., 2])], -1)


def r_from_absT(aT):
    """Solve |T(r)| = r k(r) = aT (increasing) by Newton in log r."""
    lr = torch.log(aT / T.K0)
    for _ in range(80):
        x = lr.clone().requires_grad_(True)
        f = torch.log(T.k_of(torch.exp(x))) + x - torch.log(aT)
        g, = torch.autograd.grad(f.sum(), x)
        step = (f / g).detach()
        lr = lr - step.clamp(-1.0, 1.0)
        if step.abs().max() < 1e-15:
            break
    return torch.exp(lr)


def _cusp_logzeta(afun, C, a_target, t):
    """log|zeta| with afun(e^l, t) = a_target mod 1 (afun decreasing in l: x4_slice = l / P with P < 0)."""
    l = (a_target + 0.5) * C.P(t)
    for _ in range(80):
        x = l.clone().requires_grad_(True)
        f = T.wrap(afun(torch.exp(x).to(CD), t) - a_target)
        g, = torch.autograd.grad(f.sum(), x)
        step = (f / g).detach()
        l = l - step.clamp(-0.5 * C.P(t).abs(), 0.5 * C.P(t).abs())
        if step.abs().max() < 1e-14:
            break
    return l


def W_from_orbit(v, X):
    """A point of W (mid representation, x1 = 0) with Lambda(point) = (v, X), via the explicit chain (exact outside the
    linear-blend shells; there it is an initial guess for the Newton polish)."""
    u1 = moebius_inv(P_final_inv(v, X))
    r = r_from_absT(u1[..., 0].abs() / u1[..., 1].abs())
    p1, p2 = L.arc_phases(r)
    u = torch.stack([p1 * u1[..., 0], p2 * u1[..., 1]], -1)
    t = torch.polar(r, torch.angle(u[..., 0] / u[..., 1]))
    logt = torch.log(t)
    x4 = torch.zeros_like(r)
    a1 = torch.angle(u[..., 1]) / T.TWO_PI
    a2p = torch.angle(u[..., 0]) / T.TWO_PI
    m = (r >= 0.5) & (r <= 1)
    x4 = torch.where(m, a1 + T.xi4(t) + 0.5, x4)
    m = (r > 1) & (r <= 2)
    lam = T.lam_inf(r)
    a2 = -a2p + lam * T.delta(t) + (1 - lam) * 0.5
    x4 = torch.where(m, -a2 - 0.5, x4)
    zeta = torch.ones_like(t)
    c0, ci = r < 0.5, r > 2
    if c0.any():
        l = _cusp_logzeta(T.a1_c0, T.CUSP0, a1[c0], t[c0])
        zeta = zeta.clone(); zeta[c0] = torch.exp(l).to(CD)
    if ci.any():
        s = 1 / t[ci]
        l = _cusp_logzeta(T.a2_ci, T.CUSPINF, -a2p[ci] + 0.5, s)        # lam_inf = 0 there: a2 = -a2' + 1/2
        zeta = zeta.clone(); zeta[ci] = torch.exp(l).to(CD)
    return t, logt, x4, zeta, c0, ci


def _lattice_from_cusp0(zeta, t):
    """Mid lattice coordinates (x1, x4) of the cusp-0 slice point (zeta, t), via Phi."""
    logt = torch.log(t)
    s = logt / (2j * math.pi)
    x = T.lattice_coords(-1j * torch.log(zeta) + math.pi * s * (s - 1), 2 * t * s, t, logt)
    return x[..., 0], x[..., 3], logt


def _lattice_from_cuspinf(zeta, s):
    t = 1 / s
    logt = -torch.log(s)
    x = T.lattice_coords(-1j * torch.log(zeta), torch.zeros_like(zeta), t, logt)
    return x[..., 0], x[..., 3], logt


def initial_guess(q):
    """Regular-part guess in mid lattice form (t, logt, x1..x4), exact outside the blend shells.
    Finv bypasses this polar chain in the fixed-point caps; central-fibre inputs elsewhere remain T14.
    """
    Z1, Z2 = torch.complex(q[..., 0], q[..., 1]), torch.complex(q[..., 2], q[..., 3])
    X = q[..., 4]
    h = 2 * Z1 * Z2.conj()
    v = torch.stack([h.real, h.imag, Z1.abs() ** 2 - Z2.abs() ** 2], -1)
    t, logt, x4, zeta, c0, ci = W_from_orbit(v, X)
    x1 = torch.zeros_like(x4)
    if c0.any():
        a, b, lg = _lattice_from_cusp0(zeta[c0], t[c0])
        x1[c0], x4[c0], logt = a, b, logt.clone().index_put_((c0.nonzero().flatten(),), lg)
    if ci.any():
        a, b, lg = _lattice_from_cuspinf(zeta[ci], 1 / t[ci])
        x1[ci], x4[ci], logt = a, b, logt.clone().index_put_((ci.nonzero().flatten(),), lg)
    W1, W2, _ = L.Psi_W(t, x1, x4, logt)
    th = torch.angle(Z1 * W1.conj() + Z2 * W2.conj())
    x1 = x1 + th / T.TWO_PI
    phi = torch.stack([torch.atan2(q[..., 6], q[..., 5]), torch.atan2(q[..., 8], q[..., 7])], -1) / T.TWO_PI
    x = torch.stack([x1, torch.zeros_like(x1), torch.zeros_like(x1), x4], -1)
    corr = x[..., 1:3] - Fm.phi_mid(t, x, logt)                         # phi = (x2, x3) - corr
    x23 = phi + corr
    return t, logt, torch.stack([x1, x23[..., 0], x23[..., 1], x4], -1)


def _filling_point(C, zz, ww, tt, logtt):
    """Filling coordinates (zeta = e^{i zz}, w, t) -> U_0 chart point (x, y, w) whose T^2-slice lies in the fundamental
    slab |mu| <= 1/2, using the cusp generator G4 = (nu t zeta, w + e4bar(t)) (its powers act as psi_hat^n)."""
    zeta = torch.exp(1j * zz)
    C2, C3 = T.real_coeffs(ww, C.e2(tt), C.e3(tt))
    zs = zeta * torch.exp(-(C2 * C.al2(tt) + C3 * C.al3(tt)))
    n = torch.round(C.mu_torus(zs, tt)).detach()
    zeta = zeta * torch.exp(-n * (C.lognu + logtt))
    ww = ww - n * C.e4(tt)
    return torch.stack([(tt / zeta).real, (tt / zeta).imag, zeta.real, zeta.imag, ww.real, ww.imag], -1)


def to_chart(t, logt, x):
    """The Y-point with mid lattice coordinates x, in the output charts: {'mid': (mask, p6), 'c0': ..., 'ci': ...};
    mid for 0.2 <= |t| <= 5, the fillings below / above."""
    z, w = T.point_from_lattice(t, logt, x)
    r = t.abs()
    s = logt / (2j * math.pi)
    return {'mid': ((r >= 0.2) & (r <= 5), torch.stack([z.real, z.imag, w.real, w.imag, t.real, t.imag], -1)),
            'c0': (r < 0.2, _filling_point(T.CUSP0, z - math.pi * s * (s - 1), w - 2 * t * s, t, logt)),
            'ci': (r > 5, _filling_point(T.CUSPINF, z - 1j * w / t, w / t ** 2, 1 / t, -logt))}


CHART_F = {'mid': Fm.F_mid_real, 'c0': Fm.F_c0_real, 'ci': Fm.F_ci_real}


def linear_ball_inverse(q, key):
    """Closed inverse in a fixed-point ball, including its pole (programme section 3.4, C8).

    First invert the unitary linear lift on the T^2 slice, then recover C from
    phi=C-H*c_tilde and undo the slice. No logarithm of a vanishing coordinate.
    """
    C, scale = (T.CUSP0, L.S0) if key == 'c0' else (T.CUSPINF, L.SINF)
    Z1, Z2 = torch.complex(q[..., 0], q[..., 1]), torch.complex(q[..., 2], q[..., 3])
    second = -1j * Z2 if key == 'c0' else Z2
    q1 = (Z1 + second) / math.sqrt(2 * scale)
    ys = (Z1 - second) / math.sqrt(2 * scale)
    xs = q1.conj() / math.sqrt(C.a0)
    t = xs * ys
    H = C.mu_chart(xs, ys)
    c2, c3 = C.c_tilde(t)
    phi2 = torch.atan2(q[..., 6], q[..., 5]) / T.TWO_PI
    phi3 = torch.atan2(q[..., 8], q[..., 7]) / T.TWO_PI
    C2, C3 = phi2 + H * c2, phi3 + H * c3
    exponent = C2 * C.al2(t) + C3 * C.al3(t)
    x, y = xs * torch.exp(-exponent), ys * torch.exp(exponent)
    w = C2 * C.e2(t) + C3 * C.e3(t)
    return torch.stack([x.real, x.imag, y.real, y.imag, w.real, w.imag], -1)


def Finv(q, iters=8):
    """F^-1 for points q of S^4 x T^2 in R^9: returns {chart: (indices, p6, residual |F(p) - q|_inf)}."""
    out = {}
    z2 = (q[..., :4] ** 2).sum(-1)
    balls = {'c0': (q[..., 4] < 0) & (z2 <= L.S0 * L.EPS1),
             'ci': (q[..., 4] > 0) & (z2 <= L.SINF * L.EPS1)}
    for key, mask in balls.items():
        if mask.any():
            idx = mask.nonzero().flatten()
            p = linear_ball_inverse(q[mask], key)
            res = (CHART_F[key](p) - q[mask]).abs().max(-1).values
            out[key] = (idx, p, res)
    remaining = (~(balls['c0'] | balls['ci'])).nonzero().flatten()
    if not len(remaining):
        return out
    t, logt, x = initial_guess(q[remaining])
    for key, val in to_chart(t, logt, x).items():
        mask, P = val
        if mask.any():
            idx = remaining[mask]
            p, res = newton_polish(CHART_F[key], P[mask], q[idx], iters)
            entry = (idx, p, res)
            out[key] = tuple(torch.cat([old, new]) for old, new in zip(out[key], entry)) if key in out else entry
    if sum(len(entry[0]) for entry in out.values()) != len(q):
        raise RuntimeError('F inverse failed to assign a chart to every input point')
    return out


def newton_polish(f, p, q, iters=8):
    for _ in range(iters):
        Fp = f(p)
        A = torch.stack([jvp(f, (p,), (EYE6[k].expand_as(p),))[1] for k in range(6)], -1)
        d = q - Fp
        d[..., 5:9] = q[..., 5:9] - Fp[..., 5:9]
        p = p + torch.einsum('bkc,bc->bk', torch.linalg.pinv(A), d)
    return p, (f(p) - q).abs().max(-1).values
