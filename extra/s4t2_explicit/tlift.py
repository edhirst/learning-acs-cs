"""The S^1-equivariant diffeomorphism Psi_W: W -> S^4 and the diffeomorphism F = (Psi_W o pi, phi): Y -> S^4 x T^2,
in torch (float64, differentiable to any order). Design and verification: docs/s4t2_explicit_diffeo.md, section 12.

W-points are passed in the 'mid' representation (t, x1, x4, logt): lattice coordinates of the T^2-slice point
(x2 = x3 = 0) over t in C*, with a branch logt of log t (changing it by 2 pi i changes x1 by -x4 and nothing else).
Near the central fibres the chart representations are used (functions with suffix _chart): U_0-coordinates (x, y)
of the slice point of the filling, |mu| < 1/2.
"""
from __future__ import annotations

import json
import math
import os

import torch

from extra.s4t2_explicit import tcore as T
from extra.s4t2_explicit.tcore import CD, RD, _scatter, bump_step, where_safe

# ------------------------------------------------------------------ the weight-1 function E (section 12.1)
R0, RP, R1 = 0.6, 0.7, 0.8      # cusp formulas below R0 (above 1/R0); plateau radius RP; linear blend on [R1, 1/R1]


def finite_phase_log(g, n, derivative_bound, chart=None):
    """Finite real logarithm of g:S2->S1 with a CERTIFIED ||Dg||<=L.

    n has shape (batch,3); g accepts arbitrary leading batch dimensions.
    The proof and fixed count floor(8L)+1 are in programme 10.6. No sampled
    estimate is a certificate for L. This helper does not construct the
    campaign's missing finite circle lift or a bound for its mismatch.
    Optional north/south selection is for overlap checks (stereographic |a|<=2).
    """
    if not math.isfinite(derivative_bound) or derivative_bound < 0:
        raise ValueError('a finite nonnegative certified derivative bound is required')
    if chart not in (None,'north','south'):
        raise ValueError('chart must be north, south, or None')
    count = math.floor(8*derivative_bound)+1
    times = torch.arange(count+1,dtype=n.dtype,device=n.device)/count
    def primitive(a, north):
        if (a.square().sum(-1)>4).any():
            raise ValueError('the certified stereographic domain is |a|<=2')
        p = a[...,None,:]*times[:,None]
        r2 = p.square().sum(-1,keepdim=True)
        third = (1-r2) if north else (r2-1)
        sphere = torch.cat([2*p,third],-1)/(1+r2)
        values = g(sphere)
        increments = values[...,1:]*values[...,:-1].conj()
        if (increments.real<=0).any():
            raise ValueError('phase increments violate the certified branch domain')
        return torch.atan2(increments.imag,increments.real).sum(-1)
    pole = n.new_tensor([[0.,0.,1.]])
    theta0 = torch.angle(g(pole))[0]
    equator = n.new_tensor([[1.,0.]])
    offset = primitive(equator,True)[0]-primitive(equator,False)[0]
    north = n[:,2]>=0 if chart is None else torch.full_like(n[:,2],chart=='north',dtype=torch.bool)
    out = torch.zeros_like(n[:,0])
    for mask,sign in ((north,1),(~north,-1)):
        if mask.any():
            a = n[mask,:2]/(1+sign*n[mask,2:3])
            value = theta0+primitive(a,sign==1)
            if sign==-1:
                value = value+offset
            out = _scatter(out,mask,value)
    return out


def finite_projection_lift(projectors, z):
    """Exact circle isomorphism from a finite chain of nearby rank-one projectors.

    projectors has shape (...,steps+1,d,d), z shape (...,d), with unit z in
    im P0. A proof that ||Pj-Pj-1||<1 is required at every step; the function
    checks only nonvanishing of the actual evaluations. Reverse the chain to
    invert. This is not an approximation to parallel transport.
    """
    for j in range(1,projectors.shape[-3]):
        z = (projectors[...,j,:,:]@z[...,None])[...,0]
        norm = z.norm(dim=-1,keepdim=True)
        if (norm==0).any():
            raise ValueError('a projection step vanished')
        z = z/norm
    return z


def _bundle_height_weight(h):
    """One on |h|<=1/4, zero on |h|>=2/5, with flat endpoints."""
    return 1-bump_step((h.abs()-.25)/.15)


def _bundle_cusp_torus(C,zeta,t):
    """Four complex components on a punctured cusp; deck invariant."""
    height = C.mu_torus(zeta,t)
    pairs = []
    for centre in (0.,.5):
        shift = torch.floor(centre+.5-height).detach()
        weight = _bundle_height_weight(height+shift-centre)
        active = weight>0
        pair = torch.zeros(zeta.shape+(2,),dtype=CD,device=zeta.device)
        if active.any():
            yy = zeta[active]*torch.exp(shift[active]*torch.log(C.multiplier(t[active])))
            xx = t[active]/yy
            rad = C.a0*T.abs2(xx)+T.abs2(yy)
            frame = torch.stack([math.sqrt(C.a0)*xx.conj(),yy],-1)/torch.sqrt(rad)[...,None]
            pair = _scatter(pair,active,weight[active,None]*frame)
        pairs.append(pair)
    return torch.cat(pairs,-1)


def _bundle_cusp_chart(C,x,y):
    """Central-fibre-safe version on |H|<=.6, excluding the fixed point."""
    height = C.mu_chart(x,y)
    rad = C.a0*T.abs2(x)+T.abs2(y)
    if (rad<=0).any() or (height.abs()>.6).any():
        raise ValueError('bundle chart requires rad>0 and |H|<=.6')
    frame = torch.stack([math.sqrt(C.a0)*x.conj(),y],-1)/torch.sqrt(rad)[...,None]
    first = _bundle_height_weight(height)[...,None]*frame
    shifted = torch.where(height<0,height+.5,height-.5)
    weight = _bundle_height_weight(shifted)
    second = weight[...,None]*frame
    negative = (height<0)&(weight>0)
    if negative.any():
        xx,yy = x[negative],y[negative]
        m = C.m(xx*yy)
        xn,yn = 1/(m*yy),m*xx*yy**2
        rn = C.a0*T.abs2(xn)+T.abs2(yn)
        fn = torch.stack([math.sqrt(C.a0)*xn.conj(),yn],-1)/torch.sqrt(rn)[...,None]
        second = _scatter(second,negative,weight[negative,None]*fn)
    return torch.cat([first,second],-1)


def _bundle_normalize(v,with_norm):
    norm2 = T.abs2(v).sum(-1)
    e = v/torch.sqrt(norm2)[...,None]
    return (e,norm2) if with_norm else e


def bundle_embedding_middle(t,x1,x4,logt,*,finite_height=False,with_norm=False):
    """Explicit equivariant embedding of W minus its nodes into S19.

    Programme 11.1 proves 1/2<=unnormalized norm^2<=2, deck/branch
    independence and central extension. This is a concrete bundle embedding,
    not yet a complete meridian lift. Gaussian height is the legacy default;
    finite_height=True makes this forward block finite elementary.
    """
    r = t.abs()
    b0 = 1-bump_step((r-.3)/.2)
    bi = 1-bump_step((1/r-.3)/.2)
    bm = 1-b0-bi
    v = torch.zeros(t.shape+(10,),dtype=CD,device=t.device)
    C0,Ci = T.height_cusps(finite_height)
    for C,convert,beta,start in ((C0,T.to_cusp0,b0,0),(Ci,T.to_cuspinf,bi,4)):
        active = beta>0
        if active.any():
            z,tc = convert(t[active],x1[active],x4[active],logt[active])
            block = _bundle_cusp_torus(C,z,tc)*beta[active,None]
            pad = torch.zeros(block.shape[:-1]+(10,),dtype=CD,device=t.device)
            pad[...,start:start+4] = block
            v = _scatter(v,active,v[active]+pad)
    active = bm>0
    if active.any():
        tt,one,four,lt = t[active],x1[active],x4[active],logt[active]
        components = []
        for centre in (0.,math.pi):
            angle = torch.angle(tt if centre==0 else -tt)
            weight = 1-bump_step((angle.abs()-math.pi/2)/(math.pi/4))
            branch = angle+centre
            phase = torch.exp(1j*(2*math.pi*one-(branch-lt.imag)*four))
            components.append(weight*phase)
        block = bm[active,None]*torch.stack(components,-1)
        pad = torch.zeros(block.shape[:-1]+(10,),dtype=CD,device=t.device)
        pad[...,8:10] = block
        v = _scatter(v,active,v[active]+pad)
    return _bundle_normalize(v,with_norm)


def bundle_embedding_chart(C,x,y,*,with_norm=False):
    """Same embedding on a pure cusp |xy|<=.3, |H|<=.6, off the node.

    Includes either central axis. On rad<=.004 the only nonzero components
    are exactly (sqrt(a0)*conj(x),y)/sqrt(rad), in the corresponding cusp block.
    """
    if ((x*y).abs()>.3).any():
        raise ValueError('bundle_embedding_chart requires |xy|<=.3')
    block = _bundle_cusp_chart(C,x,y)
    start = 4 if C.name=='inf' else 0
    v = torch.zeros(x.shape+(10,),dtype=CD,device=x.device)
    v[...,start:start+4] = block
    return _bundle_normalize(v,with_norm)


def connection_torus(C, zeta, t, dzeta, dt):
    """Explicit cusp circle connection applied to a tangent vector, t*zeta != 0.

    This independent connection supports the quadrature-lift design; it does not
    change the existing Psi functions. Pass T.FINITE_CUSP0/INF to use the
    finite height. The locally finite chi-average has denominator >= 1/2.
    """
    mu = C.mu_torus(zeta,t)
    n = -torch.round(mu).detach()[...,None]+torch.arange(-1,2,dtype=RD,device=t.device)
    weights = T.chi(mu[...,None]+n)
    rho = torch.log(zeta.abs())[...,None]+n*C.P(t)[...,None]
    f = torch.sigmoid(math.log(C.a0)+torch.log(T.abs2(t))[...,None]-4*rho)
    mean = (weights*(n-f)).sum(-1)/weights.sum(-1)
    return (dzeta/zeta).imag+mean*(dt/t).imag


def connection_chart(C, x, y, dx, dy):
    """The same connection in U_0, including either axis, away from the node.

    This three-term chart expression is used on |mu|<=0.6. Each noncentral
    coordinate divided by below is nonzero on its coefficient's support.
    """
    mu = C.mu_chart(x,y)
    if (mu.abs() > 0.6).any():
        raise ValueError('connection_chart requires |mu| <= 0.6')
    def local(xx,yy,dxx,dyy):
        rad = C.a0*T.abs2(xx)+T.abs2(yy)
        if (rad == 0).any():
            raise ValueError('the circle connection is undefined at a fixed point')
        return (-C.a0*xx.conj()*dxx+yy.conj()*dyy).imag/rad
    w0,wp,wm = T.chi(mu),T.chi(mu+1),T.chi(mu-1)
    out = w0*local(x,y,dx,dy)
    t,dt = x*y,dx*y+x*dy
    tt = T.abs2(t)
    m = C.nu*torch.exp(-2*tt/(1-tt))
    dm = -m*4*(t.conj()*dt).real/(1-tt)**2
    plus = wp > 0
    if plus.any():
        xx,yy,dxx,dyy,mm,dmm = (v[plus] for v in (x,y,dx,dy,m,dm))
        xn,yn = 1/(mm*yy),mm*xx*yy**2
        dxn = -xn*(dmm/mm+dyy/yy)
        dyn = dmm*xx*yy**2+mm*dxx*yy**2+2*mm*xx*yy*dyy
        out = _scatter(out,plus,out[plus]+wp[plus]*local(xn,yn,dxn,dyn))
    minus = wm > 0
    if minus.any():
        xx,yy,dxx,dyy,mm,dmm = (v[minus] for v in (x,y,dx,dy,m,dm))
        xn,yn = mm*xx**2*yy,1/(mm*xx)
        dxn = dmm*xx**2*yy+2*mm*xx*yy*dxx+mm*xx**2*dyy
        dyn = -yn*(dmm/mm+dxx/xx)
        out = _scatter(out,minus,out[minus]+wm[minus]*local(xn,yn,dxn,dyn))
    return out/(w0+wp+wm)


def connection_middle(t, x1, x4, logt, dt, dx1, dx4,*,finite_height=False):
    """Global connection on the regular part, including the two cusp collars.

    Branch changes act on both x1 and dx1. The weighted connection extends
    through the central fibres by connection_chart; no fitted phase is used.
    """
    from torch.func import jvp
    r = t.abs()
    b0 = 1-bump_step((r-0.3)/0.2)
    bi = 1-bump_step((1/r-0.3)/0.2)
    out = (1-b0-bi)*(2*math.pi*dx1+logt.imag*dx4)
    C0,Ci = T.height_cusps(finite_height)
    for C, fun, weight in ((C0,T.to_cusp0,b0),(Ci,T.to_cuspinf,bi)):
        mask = weight > 0
        if mask.any():
            vals = (t[mask],x1[mask],x4[mask],logt[mask])
            tans = (dt[mask],dx1[mask],dx4[mask],dt[mask]/t[mask])
            (z,tc),(dz,dc) = jvp(fun,vals,tans)
            out = _scatter(out,mask,out[mask]+weight[mask]*connection_torus(C,z,tc,dz,dc))
    return out


def algebraic_k(r):
    """Balanced algebraic base profile, programme 9.2 (both endpoint scales)."""
    A = (4*math.pi)**2
    return torch.sqrt((A+math.exp(4)*r*r)/(1+A*r*r))


def suspended_hopf(q):
    """Algebraic Hopf suspension S4->S3; q=(Re z1,Im z1,Re z2,Im z2,X).

    This literature-inspired target model is independent of the maps from W.
    It is not an identification of the repo's W with the standard S4.
    """
    a,b,c,d,h = q.unbind(-1)
    return torch.stack([2*(a*c+b*d),2*(b*c-a*d),a*a+b*b-c*c-d*d,
                        h*torch.sqrt(2-h*h)],-1)


def sphere3_boost(p,inverse=False):
    """Explicit conformal diffeomorphism of S3, denominator >=sqrt(2)-1."""
    sign = -1 if inverse else 1
    den = math.sqrt(2)+sign*p[...,0]
    return torch.cat([((math.sqrt(2)*p[...,0]+sign)/den)[...,None],
                      p[...,1:]/den[...,None]],-1)


def algebraic_torus_fibration(q):
    """h o k o Sigma(h): algebraic S4->S2 with critical values (0,+/-1,0).

    A target for the alternative fibration strategy of programme 10.4, not a
    complex structure and not yet matched to the given complex threefold.
    """
    a,b,c,d = sphere3_boost(suspended_hopf(q)).unbind(-1)
    return torch.stack([2*(a*c+b*d),2*(b*c-a*d),a*a+b*b-c*c-d*d],-1)


def algebraic_hopf_cusp(tc,a,infinity=False):
    """h_B with the algebraic k profile, including the central fibre."""
    r2 = tc.real**2+tc.imag**2
    A = (4*math.pi)**2
    if infinity:
        ratio = tc*torch.sqrt((r2+A)/(A*r2+math.exp(4)))
        u1 = torch.exp(2j*math.pi*(-a+.5))/torch.sqrt(1+T.abs2(ratio))
        return torch.stack([u1,u1*ratio],-1)
    # Use the polynomial |tc|^2: abs(tc)**2 has incorrect higher AD jets at 0.
    ratio = tc*torch.sqrt((A+math.exp(4)*r2)/(1+A*r2))
    u2 = torch.exp(2j*math.pi*a)/torch.sqrt(1+T.abs2(ratio))
    return torch.stack([ratio*u2,u2],-1)


def transport_base_local(p, infinity=False, modified=True,*,algebraic_profile=False):
    """Certified local base map of programme section 8, |p|<=0.008.

    R is identity for this new construction. The raw formula is I2a--d;
    optional blending makes it exactly linear for |p|<=0.004.
    """
    xx = 2*math.pi*p
    if algebraic_profile:
        rho2 = p[...,0]**2+p[...,1]**2
        A = (4*math.pi)**2
        beta = torch.sqrt((1+math.exp(4)*rho2/(4*A))/(1+A*rho2/4))
        xx = xx*torch.stack([beta,-beta if infinity else beta,torch.ones_like(beta)],-1)
    elif infinity:
        rho2 = p[...,0]**2+p[...,1]**2
        outer = rho2 > (2*math.exp(-6))**2
        rho = torch.sqrt(torch.where(outer,rho2,torch.ones_like(rho2)))
        beta = torch.where(outer,T.KINF/T.k_of(2*math.exp(-2)/rho),torch.ones_like(rho))
        xx = xx*torch.stack([beta,-beta,torch.ones_like(beta)],-1)
    a,b,c = xx.unbind(-1)
    den = torch.sqrt(1+a*a+b*b)+a*torch.cos(c)+b*torch.sin(c)
    raw = torch.stack([a*torch.sin(c)-b*torch.cos(c),-torch.sin(c),den-torch.cos(c)],-1)/den[...,None]
    if infinity:
        raw = torch.stack([-raw[...,1],-raw[...,0],raw[...,2]],-1)/2
    if not modified:
        return raw
    linear = p@((SINF if infinity else S0)*(QINF if infinity else Q0).to(p)).T
    # The cutoff is constant on the inner ball. Avoid differentiating |p| at
    # zero even in an inactive expression: nested AD must see a smooth graph.
    radius = torch.sqrt((p*p).sum(-1).clamp_min(EPS1**2))
    cutoff = bump_step((radius-EPS1)/(EPS2-EPS1))
    return linear+cutoff[...,None]*(raw-linear)


def P_algebraic(w):
    """Algebraic replacement for P_final, on 0<|w|<infinity (programme 9.2).

    The fixed points use transport_base_chart's polynomial-radius formulas.
    This block and its inverse contain no profile root or cutoff.
    """
    r2 = (w*w).sum(-1)
    r = torch.sqrt(r2)
    direction = torch.stack([w[...,0],-w[...,2],-w[...,1]],-1)
    return 2*direction/(1+r2)[...,None], (r-1)/torch.sqrt(1+r2)


def transport_base_chart(C,x,y,*,algebraic_radial=False,algebraic_profile=False):
    """Unfitted base in a reduced cusp chart; optional algebraic radial profile.

    False retains section 8 verbatim. True is section 9.2's conjugate base;
    the local shell and all height choices are unchanged.
    """
    # Own the storage before complex real/imag views in nested forward AD.
    x = x.clone(memory_format=torch.contiguous_format)
    y = y.clone(memory_format=torch.contiguous_format)
    xr,xi,yr,yi = x.real,x.imag,y.real,y.imag
    nx,ny = xr*xr+xi*xi,yr*yr+yi*yi
    # Real polynomial form also avoids a nested-AD layout bug after scaled conj.
    p = torch.stack([2*math.sqrt(C.a0)*(xr*yr-xi*yi),
                     -2*math.sqrt(C.a0)*(xr*yi+xi*yr),C.a0*nx-ny],-1)
    rad_up = C.a0*nx+ny
    small = rad_up <= EPS2
    v = torch.zeros(x.shape+(3,),dtype=RD,device=x.device)
    X = torch.zeros(x.shape,dtype=RD,device=x.device)
    infinity = C.name == 'inf'
    if small.any():
        vv = transport_base_local(p[small],infinity,algebraic_profile=algebraic_profile)
        rad = rad_up[small]
        inner = rad <= EPS1
        nv = torch.linalg.vector_norm(torch.where(inner[...,None],torch.ones_like(vv),vv),dim=-1)
        nv = torch.where(inner,(SINF if infinity else S0)*rad,nv)
        if algebraic_radial:
            vv = 2*vv/(1+nv*nv)[...,None]
            hh = (1 if infinity else -1)*(1-nv)/torch.sqrt(1+nv*nv)
        else:
            hh = (1 if infinity else -1)*torch.sqrt(1-nv)
        v,X = _scatter(v,small,vv),_scatter(X,small,hh)
    if (~small).any():
        xx,yy = x[~small],y[~small]
        tc = xx*yy
        # The chart form also covers the central axes. Away from the cusp core,
        # use the original height/log-radius blend, which equals it at r=.3.
        inner = tc.abs() <= 0.3
        aa = T.mu_new_chart(C,xx,yy)
        if (~inner).any():
            aa = _scatter(aa,~inner,T.cusp_angle(C,yy[~inner],tc[~inner]))
        if algebraic_profile:
            u = algebraic_hopf_cusp(tc,aa,infinity)
        else:
            u = T.hopf_from_a2_s(tc,aa) if infinity else T.hopf_from_a1(tc,aa)
        vv,hh = (P_algebraic if algebraic_radial else T.P_final)(transport_stereographic(u))
        v,X = _scatter(v,~small,vv),_scatter(X,~small,hh)
    return v,X


def transport_base_middle(t,x1,x4,logt,*,algebraic_radial=False,algebraic_profile=False,finite_height=False):
    """The same new base map in a middle chart, with R=Id."""
    r = t.abs()
    cusp0,cuspi = r<0.5,r>2
    mid = ~(cusp0|cuspi)
    v = torch.zeros(t.shape+(3,),dtype=RD,device=t.device)
    X = torch.zeros(t.shape,dtype=RD,device=t.device)
    if mid.any():
        tt,aa = t[mid],x4[mid]
        if algebraic_profile:
            lo = tt.abs()<=1
            ratio = tt*algebraic_k(tt.abs())
            a1 = T.a1_mid(tt,aa)
            lam = T.lam_inf(tt.abs())
            a2p = -T.a2_mid(tt,aa)+lam*T.delta(tt)+(1-lam)*.5
            u2 = torch.exp(2j*math.pi*a1)/torch.sqrt(1+T.abs2(ratio))
            u1 = torch.exp(2j*math.pi*a2p)/torch.sqrt(1+1/T.abs2(ratio))
            u = torch.where(lo[:,None],torch.stack([ratio*u2,u2],-1),torch.stack([u1,u1/ratio],-1))
        else:
            u = T.h_B_mid(tt,aa)
        vv,hh = (P_algebraic if algebraic_radial else T.P_final)(transport_stereographic(u))
        v,X = _scatter(v,mid,vv),_scatter(X,mid,hh)
    C0,Ci = T.height_cusps(finite_height)
    for C,mask,convert in ((C0,cusp0,T.to_cusp0),(Ci,cuspi,T.to_cuspinf)):
        if mask.any():
            z,tc = convert(t[mask],x1[mask],x4[mask],logt[mask])
            _,xx,yy = C.reduce(z,tc)
            vv,hh = transport_base_chart(C,xx,yy,algebraic_radial=algebraic_radial,algebraic_profile=algebraic_profile)
            v,X = _scatter(v,mask,vv),_scatter(X,mask,hh)
    return v,X


def transport_stereographic(u):
    """Direct stereographic formula, regular also at u=(1,0).

    The opposite point u=(-1,0) is the infinity node and uses the local base
    formula. This avoids the removable 0/0 of the legacy rationalization at +1.
    """
    vec = torch.stack([u[...,0].imag,u[...,1].real,u[...,1].imag],-1)
    return vec/(1+u[...,0].real)[...,None]-u.real.new_tensor([0.,1.,0.])


def elementary_torus_shift_chart(C,x,y):
    """phi_E-phi in a REDUCED cusp chart (programme 9.3).

    The correction is supported in rad<.008; near its support |H|<.01,
    so the representative is unique. Outside this chart use the middle form.
    """
    radius = C.a0*T.abs2(x)+T.abs2(y)
    active = radius<EPS2
    out = torch.zeros(x.shape+(2,),dtype=RD,device=x.device)
    if active.any():
        xx,yy = x[active],y[active]
        H = C.mu_chart(xx,yy)
        c2,c3 = C.c_tilde(xx*yy)
        cutoff = 1-bump_step((radius[active]-EPS1)/(EPS2-EPS1))
        out = _scatter(out,active,(cutoff*H)[...,None]*torch.stack([c2,c3],-1))
    return out


def elementary_torus_shift_middle(t,x1,x4,logt,*,finite_height=False):
    """Global real torus correction on the regular part, invariant under deck/branch."""
    out = torch.zeros(t.shape+(2,),dtype=RD,device=t.device)
    C0,Ci = T.height_cusps(finite_height)
    for C,mask,convert in ((C0,t.abs()<.01,T.to_cusp0),(Ci,t.abs()>25,T.to_cuspinf)):
        if mask.any():
            z,tc = convert(t[mask],x1[mask],x4[mask],logt[mask])
            _,x,y = C.reduce(z,tc)
            out = _scatter(out,mask,elementary_torus_shift_chart(C,x,y))
    return out


def h0(r):
    return r + (RP - r) * bump_step((r - R0) / (R1 - R0))


def hinf(r):
    return r + (1 / RP - r) * bump_step((1 / R0 - r) / (1 / R0 - 1 / R1))


def lam_mid(r):
    return bump_step((torch.log(r) - math.log(R1)) / (math.log(1 / R1) - math.log(R1)))


def _at_radius(rr, th):
    """t and its continuous log at radius rr and angle th (th = Im logt of the original point)."""
    return torch.polar(rr, th), torch.complex(torch.log(rr), th)


def E_W(t, x1, x4, logt):
    """E on W (mid representation):
         |t| <= R0:        cusp-0 formula E0_W;      R0 < |t| <= R1: E0_W read at the radius h0(|t|) (-> RP);
         R1 < |t| < 1/R1:  (1 - lam) E0_W(at RP) + lam Einf_W(at 1/RP)   (a linear combination: smooth);
         symmetric on the infinity side."""
    r, th = t.abs(), logt.imag
    out = torch.zeros(t.shape, dtype=CD)
    a, b = r <= R0, r >= 1 / R0
    m1 = (r > R0) & (r <= R1)
    m3 = (r >= 1 / R1) & (r < 1 / R0)
    m2 = (r > R1) & (r < 1 / R1)
    out = _scatter(out, a, T.E0_W(t[a], x1[a], x4[a], logt[a]))
    out = _scatter(out, b, T.Einf_W(t[b], x1[b], x4[b], logt[b]))
    if m1.any():
        tt, lt = _at_radius(h0(r[m1]), th[m1])
        out = _scatter(out, m1, T.E0_W(tt, x1[m1], x4[m1], lt))
    if m3.any():
        tt, lt = _at_radius(hinf(r[m3]), th[m3])
        out = _scatter(out, m3, T.Einf_W(tt, x1[m3], x4[m3], lt))
    if m2.any():
        lam = lam_mid(r[m2])
        t0, l0 = _at_radius(torch.full_like(r[m2], RP), th[m2])
        ti, li = _at_radius(torch.full_like(r[m2], 1 / RP), th[m2])
        out = _scatter(out, m2, (1 - lam) * T.E0_W(t0, x1[m2], x4[m2], l0) + lam * T.Einf_W(ti, x1[m2], x4[m2], li))
    return out


def E_chart0(x, y):
    """E near the central fibre at 0, U_0 coordinates (equals E_W through Phi; check_tlift.py)."""
    return T.CUSP0.E_chart(x, y)


def E_chartinf(x, y):
    return T.CUSPINF.E_chart(x, y)


# ------------------------------------------------------------------ Hopf identification on mid representation
def hB_chart0(x, y):
    """h_B for a chart point near the central fibre at 0 (|t| <= RA, |mu| < 1/2): a1 = mu_new exactly."""
    return T.hopf_from_a1(x * y, T.mu_new_chart(T.CUSP0, x, y))


def hB_chartinf(x, y):
    """h_B for a chart point near the central fibre at infinity (|s| <= RA)."""
    return T.hopf_from_a2_s(x * y, T.mu_new_chart(T.CUSPINF, x, y))


def hB_W(t, x1, x4, logt):
    """h_B: B -> S^3 for a W-point in the mid representation (dispatch c0 / mid / ci; check_hopf.py overlaps)."""
    r = t.abs()
    out = torch.zeros(t.shape + (2,), dtype=CD)
    a, b = r < 0.55, r > 1.8
    m = ~(a | b)
    if a.any():
        out = _scatter(out, a, T.h_B_c0(*T.to_cusp0(t[a], x1[a], x4[a], logt[a])))
    if b.any():
        out = _scatter(out, b, T.h_B_ci(*T.to_cuspinf(t[b], x1[b], x4[b], logt[b])))
    if m.any():
        out = _scatter(out, m, T.h_B_mid(t[m], x4[m]))
    return out


# ------------------------------------------------------------------ the orbit-space map Lambda (section 12.2)
def _x4_zero0(rho):
    """Mid lattice coordinate x4 of the zero of the cusp-0 formula on the torus |t| = rho (rho <= RP): |x| = |y|,
    arg t = pi (exact: there E = y + conj(x))."""
    t = torch.polar(rho, torch.full_like(rho, math.pi))
    return 0.5 * torch.log(rho) / T.CUSP0.P(t) + T.xi4(t)


def _x4_zeroinf(rho):
    """Same at infinity (rho >= 1/RP): |y'| = |nu| |x'|, arg s = -pi; +1 selects the representative continuing the
    cusp-0 branch through the middle (x4 = 0.7246 at 1/RP, versus 0.6655 at RP)."""
    t = torch.polar(rho, torch.full_like(rho, math.pi))
    s = 1 / t
    return -0.5 * (T.CUSPINF.lognu + torch.log(s.abs())) / T.CUSPINF.P(s) + 1.0


def x4_base(r):
    """1/2 (the exact zero |x| = |y| of the cusp-0 formula for r <= 0.366) blended into the exact zero of the
    cusp-infinity formula (valid for r >= 1.65); the blend sits on [1.25, 1.6], away from |t| = 1 where the
    cusp-infinity slice degenerates."""
    B = bump_step((torch.log(r) - math.log(1.25)) / (math.log(1.6) - math.log(1.25)))
    out = torch.full_like(r, 0.5)
    m = r > 1.2
    if m.any():
        out = _scatter(out, m, (1 - B[m]) * 0.5 + B[m] * _x4_zeroinf(r[m]))
    return out


ARC_INT1 = (math.log(0.78), math.log(1.28))     # fit interval of p1 (residual arg t - pi supported in [0.813, 1.232])
ARC_WIN1 = (math.log(0.80), math.log(1.25))     # window W1 = 1 on this interval
ARC_INT2 = (math.log(0.30), math.log(2.00))     # fit interval of p2 (residual x4 - x4_base supported in [0.366, 1.65])
ARC_WIN2 = (math.log(0.35), math.log(1.75))
_ARC_FIT = None


def _window(rho, inner, outer):
    return bump_step((rho - outer[0]) / (inner[0] - outer[0])) * bump_step((outer[1] - rho) / (outer[1] - inner[1]))


def cheb_eval(c, x):
    """Clenshaw: sum_k c_k T_k(x)."""
    b1 = torch.zeros_like(x)
    b2 = torch.zeros_like(x)
    for ck in reversed(c[1:]):
        b1, b2 = 2 * x * b1 - b2 + ck, b1
    return x * b1 - b2 + c[0]


def _arc_fit():
    global _ARC_FIT
    if _ARC_FIT is None:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'arc_fit.json')
        _ARC_FIT = json.load(open(path)) if os.path.exists(path) else {}
    return _ARC_FIT


def _fitted(key, rho, interval, window):
    fit = _arc_fit()
    out = torch.zeros_like(rho)
    if key not in fit:
        return out
    a, b = interval
    m = (rho > a) & (rho < b)
    if m.any():
        xs = 2 * (rho[m] - a) / (b - a) - 1
        out = _scatter(out, m, _window(rho[m], window, interval) * cheb_eval(fit[key]['coef'], xs))
    return out


def arc_model(r):
    """(arg t, x4) of the explicit model of the zero arc Gamma_E at radius r (fit_arc.py)."""
    rho = torch.log(r)
    return math.pi + _fitted('p1', rho, ARC_INT1, ARC_WIN1), x4_base(r) + _fitted('p2', rho, ARC_INT2, ARC_WIN2)


GAM_A, GAM_B = 0.005, 26.0          # the straightening R is the identity for |t| <= GAM_A and |t| >= GAM_B


def arc_phases(r):
    """(phi1, phi2) in U(1)^2: R(u) = (conj(phi1) u1, conj(phi2) u2) sends h_B of the model arc to the real great
    semicircle (0, 1) -> (-1, 0) (u1 <= 0 <= u2). R is a diffeomorphism of S^3 (r is a function of |u1|/|u2|);
    windowed so that R = id exactly near both fixed points (there the arc is exactly real anyway)."""
    out1 = torch.ones(r.shape, dtype=CD)
    out2 = torch.ones(r.shape, dtype=CD)
    m = (r > GAM_A) & (r < GAM_B)
    if m.any():
        rm = r[m]
        th, x4 = arc_model(rm)
        u = hB_W(torch.polar(rm, th), torch.zeros_like(rm), x4, torch.complex(torch.log(rm), th))
        lr = torch.log(rm)
        win = _window(lr, (math.log(0.01), math.log(13.0)), (math.log(GAM_A), math.log(GAM_B)))
        a1 = win * torch.angle(-u[..., 0])
        a2 = win * torch.angle(u[..., 1])
        out1 = _scatter(out1, m, torch.polar(torch.ones_like(rm), a1))
        out2 = _scatter(out2, m, torch.polar(torch.ones_like(rm), a2))
    return out1, out2


def Lam_from_u(u, r):
    """(v, X) = P_final(moebius(R_gamma u)) in B_S = {|v| + X^2 = 1}."""
    p1, p2 = arc_phases(r)
    return T.P_final(T.moebius(torch.stack([p1.conj() * u[..., 0], p2.conj() * u[..., 1]], -1)))


# linear models at the fixed points (section 12.2): Lambda = M v_q exactly for |v_q| <= EPS1, where
#   v_q = Hopf(q), q = (sqrt(a0) conj(x0), y0) in U_0 coordinates of the reduced point (a0 = |nu|^2);
#   M0 = 2 pi Q0, Minf = pi Qinf with Q0, Qinf signed permutation matrices (derived and verified in check_tlift.py)
EPS1, EPS2 = 0.004, 0.008
Q0 = torch.tensor([[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]])
QINF = torch.tensor([[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]])
S0, SINF = 2 * math.pi, math.pi


def v_q(C, x0, y0):
    q1, q2 = math.sqrt(C.a0) * x0.conj(), y0
    h = 2 * q1 * q2.conj()
    return torch.stack([h.real, h.imag, T.abs2(q1) - T.abs2(q2)], -1)


def lin_blend(v, X, vq, rad, M, scale, sgn):
    """Lambda_mod near a fixed point: v -> (1 - c) M vq + c v with c = 0 for rad <= EPS1, 1 for rad >= EPS2
    (rad = |vq|, smooth upstairs); X = sgn sqrt(1 - |v|). Unchanged outside rad < EPS2."""
    c = bump_step((rad - EPS1) / (EPS2 - EPS1))
    vl = vq @ M.T
    vm = (1 - c)[..., None] * vl + c[..., None] * v
    inner = rad <= EPS1
    nv = torch.linalg.norm(torch.where(inner[..., None], torch.ones_like(vm), vm), dim=-1)
    nv = torch.where(inner, scale * rad, nv)
    shell = rad < EPS2
    return (torch.where(shell[..., None], vm, v),
            torch.where(shell, sgn * torch.sqrt((1 - nv).clamp_min(0.0)), X))


M0, MINF = S0 * Q0, SINF * QINF


def node_coords_mid(t, x1, x4, logt):
    """Reduced U_0 coordinates at cusp 0 and at infinity (for the linear models), mid representation."""
    z0 = T.CUSP0.reduce(*T.to_cusp0(t, x1, x4, logt))
    zi = T.CUSPINF.reduce(*T.to_cuspinf(t, x1, x4, logt))
    return z0, zi


def Lam_W(t, x1, x4, logt):
    """Lambda_mod on the mid representation."""
    r = t.abs()
    v, X = Lam_from_u(hB_W(t, x1, x4, logt), r)
    a, b = r < 0.01, r > 25.0          # linear zones: rad >= 2 sqrt(a0)|t| <= EPS2 => |t| <= 0.004, |s| <= 0.0296
    if a.any():
        _, x0, y0 = T.CUSP0.reduce(*T.to_cusp0(t[a], x1[a], x4[a], logt[a]))
        rad = T.CUSP0.a0 * T.abs2(x0) + T.abs2(y0)
        va, Xa = lin_blend(v[a], X[a], v_q(T.CUSP0, x0, y0), rad, M0, S0, -1.0)
        v, X = _scatter(v, a, va), _scatter(X, a, Xa)
    if b.any():
        _, x0, y0 = T.CUSPINF.reduce(*T.to_cuspinf(t[b], x1[b], x4[b], logt[b]))
        rad = T.CUSPINF.a0 * T.abs2(x0) + T.abs2(y0)
        vb, Xb = lin_blend(v[b], X[b], v_q(T.CUSPINF, x0, y0), rad, MINF, SINF, 1.0)
        v, X = _scatter(v, b, vb), _scatter(X, b, Xb)
    return v, X


def _lam_chart(C, hb, x, y, invr, M, scale, sgn):
    t = x * y
    rad = C.a0 * T.abs2(x) + T.abs2(y)
    vq = v_q(C, x, y)
    v = torch.zeros(t.shape + (3,), dtype=RD)
    X = torch.full(t.shape, sgn, dtype=RD)
    gen = rad > EPS1                                   # Lam_from_u is singular exactly at the fixed point
    if gen.any():
        r = t[gen].abs()
        vg, Xg = Lam_from_u(hb(x[gen], y[gen]), 1 / r if invr else r)
        v, X = _scatter(v, gen, vg), _scatter(X, gen, Xg)
    return lin_blend(v, X, vq, rad, M, scale, sgn)


def Lam_chart0(x, y):
    return _lam_chart(T.CUSP0, hB_chart0, x, y, False, M0, S0, -1.0)


def Lam_chartinf(x, y):
    return _lam_chart(T.CUSPINF, hB_chartinf, x, y, True, MINF, SINF, 1.0)


# ------------------------------------------------------------------ the gauge G_A near Lambda^-1(M1) (section 12.3)
GA0, GAI = (0.4, 0.5), (2.0, 2.5)      # glue radii: cusp-0 local model -> e^{2 pi i x1} -> cusp-inf local model


def canonical_branch(x1, x4, logt):
    """Branch of log t with arg t in [0, 2 pi): logt - 2 pi i k, x1 + k x4 (k = floor(arg / 2 pi))."""
    k = torch.floor(logt.imag / T.TWO_PI).detach()
    return x1 + k * x4, logt - 2j * math.pi * k


def _log_local0(t, x1, x4, logt):
    """log(y0 - conj(x0)) at cusp 0 (reduced point), with consistent branches (z_s and log t as given)."""
    zs = T.to_cusp0_additive(t, x1, x4, logt)
    zeta = torch.exp(1j * zs)
    n = torch.round(T.CUSP0.mu_torus(zeta, t)).detach()
    c2, c3 = T.CUSP0.c_tilde(t)
    logmult = T.CUSP0.lognu + logt - (c2 * T.CUSP0.al2(t) + c3 * T.CUSP0.al3(t))
    logy0 = 1j * zs - n * logmult
    y0 = torch.exp(logy0)
    return logy0 + torch.log(1 - t.conj() / T.abs2(y0))


def _log_localinf(t, x1, x4, logt):
    """log(y0' - |nu| conj(x0')) at cusp infinity (reduced point), branches as given (log s = -log t)."""
    zs, s = T.to_cuspinf_additive(t, x1, x4, logt)
    zeta = torch.exp(1j * zs)
    n = torch.round(T.CUSPINF.mu_torus(zeta, s)).detach()
    c2, c3 = T.CUSPINF.c_tilde(s)
    logmult = T.CUSPINF.lognu - logt - (c2 * T.CUSPINF.al2(s) + c3 * T.CUSPINF.al3(s))
    logy0 = 1j * zs - n * logmult
    y0 = torch.exp(logy0)
    return logy0 + torch.log(1 - T.CUSPINF.e0coef * s.conj() / T.abs2(y0))


def GA_hat_W(t, x1, x4, logt):
    """Nonvanishing weight-1 function on the tube {theta_v <= THETA_B} (mid representation):
         = y0 - conj(x0) (cusp 0), e^{2 pi i x1} (canonical branch), y0' - |nu| conj(x0') (cusp infinity),
       glued by exp of a weighted log-ratio on |t| in GA0 and GAI."""
    x1c, lc = canonical_branch(x1, x4, logt)
    r = t.abs()
    expo = 2j * math.pi * x1c.to(CD)                       # log of the weight-1 factor, to be corrected below
    a = r < GA0[1]
    if a.any():
        beta = bump_step((r[a] - GA0[0]) / (GA0[1] - GA0[0]))
        L0 = 2j * math.pi * x1c[a] - _log_local0(t[a], x1c[a], x4[a], lc[a])
        expo = _scatter(expo, a, expo[a] - (1 - beta) * L0)
    b = r > GAI[0]
    if b.any():
        beta = bump_step((r[b] - GAI[0]) / (GAI[1] - GAI[0]))
        Li = 2j * math.pi * x1c[b] - _log_localinf(t[b], x1c[b], x4[b], lc[b])
        expo = _scatter(expo, b, expo[b] - beta * Li)
    return torch.exp(expo)


def GA_hat_chart0(x, y):
    return y - x.conj()


def GA_hat_chartinf(x, y):
    return y - T.CUSPINF.e0coef * x.conj()


# ------------------------------------------------------------------ the lift Psi_W (section 12.4)
THETA_A, THETA_B = math.radians(15.0), math.radians(45.0)
_PSI = None


PSI_R = (0.002, 60.0)       # psi = PSI_LO exactly for |t| <= 0.002 and PSI_HI exactly for |t| >= 60 (there the shell
PSI_LO = -math.pi / 2       # lies in the linear balls, where the two lifts must coincide: g = e^{i(beta - psi)} = 1)


def psi_of_r(r):
    """Branch normaliser: G_A = e^{i psi(|t|)} G_A_hat. psi = -pi/2 (= beta_0) near p0, PSI_HI (= beta_inf + 2 pi k)
    near p_inf; in between a fitted smooth function of log|t| (fit_psi.py) tracking the circular mean of the
    transition phase, so that |Arg g| < pi on the whole shell (verified in check_lift.py)."""
    global _PSI
    if _PSI is None:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'psi_fit.json')
        _PSI = json.load(open(path)) if os.path.exists(path) else {}
    if not _PSI:
        return torch.zeros_like(r)
    a, b = math.log(PSI_R[0]), math.log(PSI_R[1])
    pos = r > 0
    rho = torch.log(where_safe(pos, r, 1.0))
    rho = torch.where(pos, rho, torch.full_like(rho, a - 1.0))
    step = bump_step((rho - a) / (b - a))
    out = PSI_LO + (_PSI['psi_hi'] - PSI_LO) * step
    m = (rho > a) & (rho < b)
    if m.any():
        xs = 2 * (rho[m] - a) / (b - a) - 1
        win = bump_step((rho[m] - a) / _PSI['edge']) * bump_step((b - rho[m]) / _PSI['edge'])
        out = _scatter(out, m, out[m] + win * cheb_eval(_PSI['coef'], xs))
    return out


def shell_weight(v):
    """c = 0 on {theta_v <= THETA_A} (lift A), 1 on {theta_v >= THETA_B} (lift B); theta_v = angle(v, (0,0,-1))."""
    s = -v[..., 2] / torch.linalg.norm(v, dim=-1)
    return bump_step((math.cos(THETA_A) - s) / (math.cos(THETA_A) - math.cos(THETA_B)))


def assemble(v, X, E, GAh, r, use_psi=True):
    """Psi = (Z1, Z2, X) from Lambda = (v, X), E (lift B) and G_A = e^{i psi(r)} GAh (lift A), glued on the shell;
    r = |t| of the point."""
    c = shell_weight(v)
    nv = torch.linalg.norm(v, dim=-1)
    w = torch.complex(v[..., 0], v[..., 1])
    A, B = c < 1, c > 0
    Z1 = torch.zeros(X.shape, dtype=CD)
    Z2 = torch.zeros(X.shape, dtype=CD)
    if A.any():
        GA = GAh[A] * (torch.polar(torch.ones_like(X[A]), psi_of_r(r[A])) if use_psi else 1.0)
        om = GA / GA.abs()
        dm = torch.sqrt(2 * (nv[A] - v[A][..., 2]))
        z1a, z2a = w[A] * om / dm, dm / 2 * om
        sh = B[A]
        if sh.any():                                      # shell: multiply lift A by g^c, g = lift B / lift A
            ob = E[A][sh] / E[A][sh].abs()
            g = ob * om[sh].conj() * (w[A][sh] / w[A][sh].abs()).conj()
            ph = torch.polar(torch.ones_like(c[A][sh]), c[A][sh] * torch.angle(g))
            z1a, z2a = _scatter(z1a, sh, z1a[sh] * ph), _scatter(z2a, sh, z2a[sh] * ph)
        Z1, Z2 = _scatter(Z1, A, z1a), _scatter(Z2, A, z2a)
    onlyB = ~A
    if onlyB.any():
        ob = E[onlyB] / E[onlyB].abs()
        dp = torch.sqrt(2 * (nv[onlyB] + v[onlyB][..., 2]))
        Z1, Z2 = _scatter(Z1, onlyB, dp / 2 * ob), _scatter(Z2, onlyB, w[onlyB].conj() * ob / dp)
    return Z1, Z2, X


def g_shell(v, X, E, GAh, r, use_psi=True):
    """The transition phase g = lift B / lift A (diagnostic)."""
    GA = GAh * (torch.polar(torch.ones_like(X), psi_of_r(r)) if use_psi else 1.0)
    w = torch.complex(v[..., 0], v[..., 1])
    return (E / E.abs()) * (GA / GA.abs()).conj() * (w / w.abs()).conj()


# exact linear models (section 12.2): Psi = sqrt(s) U q, q = (sqrt(a0) conj(x0), y0), X = sgn sqrt(1 - s |q|^2)
SQ2 = math.sqrt(2.0)


def linear_model0(x, y):
    q1, q2 = x.conj(), y
    k = math.sqrt(S0)
    return k * (q1 + q2) / SQ2, k * 1j * (q1 - q2) / SQ2, -torch.sqrt(1 - S0 * (T.abs2(q1) + T.abs2(q2)))


def linear_modelinf(x, y):
    q1, q2 = T.CUSPINF.e0coef * x.conj(), y
    k = math.sqrt(SINF)
    return k * (q1 + q2) / SQ2, k * (q1 - q2) / SQ2, torch.sqrt(1 - SINF * (T.abs2(q1) + T.abs2(q2)))


def _psi_chart(C, lam, Efn, Gfn, lin, x, y, invr, use_psi=True):
    rad = C.a0 * T.abs2(x) + T.abs2(y)
    inner = rad <= EPS1
    Z1 = torch.zeros(x.shape, dtype=CD)
    Z2 = torch.zeros(x.shape, dtype=CD)
    X = torch.zeros(x.shape, dtype=RD)
    if inner.any():
        a, b, c = lin(x[inner], y[inner])
        Z1, Z2, X = _scatter(Z1, inner, a), _scatter(Z2, inner, b), _scatter(X, inner, c)
    o = ~inner
    if o.any():
        v, Xo = lam(x[o], y[o])
        rr = (x[o] * y[o]).abs()
        a, b, c = assemble(v, Xo, Efn(x[o], y[o]), Gfn(x[o], y[o]), 1 / rr if invr else rr, use_psi)
        Z1, Z2, X = _scatter(Z1, o, a), _scatter(Z2, o, b), _scatter(X, o, c)
    return Z1, Z2, X


def Psi_chart0(x, y, use_psi=True):
    """Psi_W near the central fibre at 0 (U_0 chart, |t| <= 0.25, |mu| < 1/2)."""
    return _psi_chart(T.CUSP0, Lam_chart0, E_chart0, GA_hat_chart0, linear_model0, x, y, False, use_psi)


def Psi_chartinf(x, y, use_psi=True):
    return _psi_chart(T.CUSPINF, Lam_chartinf, E_chartinf, GA_hat_chartinf, linear_modelinf, x, y, True, use_psi)


def Psi_W(t, x1, x4, logt, use_psi=True):
    """Psi_W on the mid representation (any t != 0; points inside the linear balls are sent to the chart formulas)."""
    r = t.abs()
    Z1 = torch.zeros(t.shape, dtype=CD)
    Z2 = torch.zeros(t.shape, dtype=CD)
    X = torch.zeros(t.shape, dtype=RD)
    done = torch.zeros(t.shape, dtype=torch.bool)
    a, b = r < 0.01, r > 25.0
    if a.any():
        _, x0, y0 = T.CUSP0.reduce(*T.to_cusp0(t[a], x1[a], x4[a], logt[a]))
        inn = (T.abs2(x0) + T.abs2(y0)) <= EPS1
        if inn.any():
            idx = a.nonzero().flatten()[inn]
            m = torch.zeros_like(done); m[idx] = True
            p, q, s = linear_model0(x0[inn], y0[inn])
            Z1, Z2, X, done = _scatter(Z1, m, p), _scatter(Z2, m, q), _scatter(X, m, s), done | m
    if b.any():
        _, x0, y0 = T.CUSPINF.reduce(*T.to_cuspinf(t[b], x1[b], x4[b], logt[b]))
        inn = (T.CUSPINF.a0 * T.abs2(x0) + T.abs2(y0)) <= EPS1
        if inn.any():
            idx = b.nonzero().flatten()[inn]
            m = torch.zeros_like(done); m[idx] = True
            p, q, s = linear_modelinf(x0[inn], y0[inn])
            Z1, Z2, X, done = _scatter(Z1, m, p), _scatter(Z2, m, q), _scatter(X, m, s), done | m
    o = ~done
    if o.any():
        v, Xo = Lam_W(t[o], x1[o], x4[o], logt[o])
        p, q, s = assemble(v, Xo, E_W(t[o], x1[o], x4[o], logt[o]), GA_hat_W(t[o], x1[o], x4[o], logt[o]),
                           r[o], use_psi)
        Z1, Z2, X = _scatter(Z1, o, p), _scatter(Z2, o, q), _scatter(X, o, s)
    return Z1, Z2, X
