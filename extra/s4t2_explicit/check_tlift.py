"""Checks of tlift.py: well-definedness on W (branch of log t, lattice translation x4 -> x4 + 1, S^1 weights),
agreement of the chart and mid representations, the exact linear models at the fixed points.
Run from the repo root: python3 extra/s4t2_explicit/check_tlift.py"""
from __future__ import annotations
import math, os, sys
import torch
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
torch.set_default_dtype(torch.float64)
from extra.s4t2_explicit import tcore as T, tlift as L

g = torch.Generator().manual_seed(21)
U = lambda n, a, b: torch.empty(n).uniform_(a, b, generator=g)
fails = []


def rep(name, err, tol):
    ok = err < tol
    if not ok:
        fails.append(name)
    print(f'{name:<62s} {err:.1e}  {"ok" if ok else "FAIL"}')


def rel(a, b):
    return ((a - b).abs() / (1 + b.abs())).max().item()


def chart0_to_mid(x, y):
    """U_0 slice point (zeta = y, t = xy, w = 0) of the cusp-0 filling -> its W-point Phi(q) in mid lattice coords."""
    t = x * y
    logt = torch.log(t)
    s = logt / (2j * math.pi)
    z = -1j * torch.log(y) + math.pi * s * (s - 1)
    w = 2 * t * s
    xl = T.lattice_coords(z, w, t, logt)
    return t, xl[..., 0], xl[..., 3], logt


def chartinf_to_mid(x, y):
    s = x * y
    t = 1 / s
    logt = -torch.log(s)
    z = -1j * torch.log(y)
    xl = T.lattice_coords(z, torch.zeros_like(z), t, logt)
    return t, xl[..., 0], xl[..., 3], logt


n = 3000
# ---------------------------------------------------------------- 1. well-definedness on W (mid representation)
r = torch.exp(U(n, math.log(1e-3), math.log(1e3)))
th = U(n, math.pi / 2 + 0.05, 3 * math.pi / 2 - 0.05)          # the tube side (G_A is only used there)
t = torch.polar(r, th); lg = torch.complex(torch.log(r), th)
x1, x4 = U(n, 0, 1), U(n, 0, 1)
lg2, x12 = lg + 2j * math.pi, x1 - x4                            # another branch of log t
for name, f in (('E_W', L.E_W), ('GA_hat_W', L.GA_hat_W)):
    f0 = f(t, x1, x4, lg)
    rep(f'1. {name}: branch change of log t', rel(f(t, x12, x4, lg2), f0), 1e-11)
    rep(f'   {name}: x4 -> x4 + 1 (lattice translation e4)', rel(f(t, x1, x4 + 1, lg), f0), 1e-11)
    rep(f'   {name}: weight 1 under S^1_a (x1 -> x1 + 0.3)', rel(f(t, x1 + 0.3, x4, lg), f0 * math.e ** 0 * torch.polar(torch.ones(n), torch.full((n,), 2 * math.pi * 0.3))), 1e-11)
v0, X0 = L.Lam_W(t, x1, x4, lg)
for name, args in (('branch change', (t, x12, x4, lg2)), ('x4 -> x4 + 1', (t, x1, x4 + 1, lg)), ('S^1_a (x1 + 0.3)', (t, x1 + 0.3, x4, lg))):
    v1, X1 = L.Lam_W(*args)
    rep(f'   Lambda_W invariant under {name}', max(rel(v1, v0), rel(X1, X0)), 1e-11)
rep('   Lambda_W lands in B_S = {|v| + X^2 = 1}', (v0.norm(dim=-1) + X0 ** 2 - 1).abs().max().item(), 1e-13)
# ---------------------------------------------------------------- 2. chart vs mid representations
for nm, toMid, Echart, Lchart, Gchart, C in (('c0', chart0_to_mid, L.E_chart0, L.Lam_chart0, L.GA_hat_chart0, T.CUSP0),
                                             ('ci', chartinf_to_mid, L.E_chartinf, L.Lam_chartinf, L.GA_hat_chartinf, T.CUSPINF)):
    m = 3000
    rr = torch.exp(U(m, math.log(1e-3), math.log(0.25)))        # |t| (resp. |s|) of the chart point
    al = U(m, math.pi / 2 + 0.05, 3 * math.pi / 2 - 0.05)        # arg t (resp. arg s): tube side
    ratio = torch.exp(U(m, -1.5, 1.5))                            # |y| / |x|
    bx = U(m, -math.pi, math.pi)
    x = torch.polar(torch.sqrt(rr / ratio), bx); y = torch.polar(torch.sqrt(rr * ratio), al - bx)
    keep = C.mu_chart(x, y).abs() < 0.45
    x, y = x[keep], y[keep]
    tm, x1m, x4m, lgm = toMid(x, y)
    rep(f'2. {nm}: E chart = E mid', rel(Echart(x, y), L.E_W(tm, x1m, x4m, lgm)), 1e-10)
    va, Xa = Lchart(x, y); vb, Xb = L.Lam_W(tm, x1m, x4m, lgm)
    rep(f'   {nm}: Lambda chart = Lambda mid', max(rel(va, vb), rel(Xa, Xb)), 1e-10)
    rep(f'   {nm}: G_A_hat chart = G_A_hat mid (|t| <= 0.25: local model)', rel(Gchart(x, y), L.GA_hat_W(tm, x1m, x4m, lgm)), 1e-10)
# ---------------------------------------------------------------- 3. linear models at the fixed points
for nm, C, Lchart, lin in (('p0  ', T.CUSP0, L.Lam_chart0, L.linear_model0), ('pinf', T.CUSPINF, L.Lam_chartinf, L.linear_modelinf)):
    m = 2000
    x = 0.04 * torch.complex(torch.randn(m, generator=g), torch.randn(m, generator=g))
    y = 0.04 * torch.complex(torch.randn(m, generator=g), torch.randn(m, generator=g))
    rad = C.a0 * x.abs() ** 2 + y.abs() ** 2
    keep = (rad > 1e-6) & (rad < L.EPS1)
    x, y = x[keep], y[keep]
    v, X = Lchart(x, y)
    Z1, Z2, Xl = lin(x, y)
    vl = torch.stack([(2 * Z1 * Z2.conj()).real, (2 * Z1 * Z2.conj()).imag, Z1.abs() ** 2 - Z2.abs() ** 2], -1)
    rep(f'3. {nm}: Hopf(linear model) = Lambda_mod (inner ball)', max(rel(vl, v), rel(Xl, X)), 1e-12)
    Ech = (L.E_chart0 if nm == 'p0  ' else L.E_chartinf)(x, y)
    rep(f'   {nm}: Z1 of the linear model has the phase of E', ((Z1 / Z1.abs()) - Ech / Ech.abs()).abs().max().item(), 1e-12)
# ---------------------------------------------------------------- 4. independent connection for quadrature transport
for C in (T.CUSP0,T.CUSPINF):
    m = 4000
    t = torch.polar(torch.exp(U(m,math.log(1e-5),math.log(0.5))),U(m,-math.pi,math.pi))
    zeta = torch.exp(U(m,-2,2)*C.P(t)).to(torch.complex128)*torch.exp(1j*U(m,-math.pi,math.pi))
    dt = t*torch.complex(U(m,-1,1),U(m,-1,1))
    dz = zeta*torch.complex(U(m,-1,1),U(m,-1,1))
    old = L.connection_torus(C,zeta,t,dz,dt)
    tt = T.abs2(t)
    mm = C.nu*torch.exp(-2*tt/(1-tt))
    dmm = -mm*4*(t.conj()*dt).real/(1-tt)**2
    zn = mm*t*zeta
    dzn = (dmm*t+mm*dt)*zeta+mm*t*dz
    rep(f'4. {C.name}: connection descends through psi_hat (4000)',
        rel(L.connection_torus(C,zn,t,dzn,dt),old),1e-11)
    rep(f'   {C.name}: connection on circle generator',
        (L.connection_torus(C,zeta,t,1j*zeta,torch.zeros_like(t))-1).abs().max().item(),1e-12)
    x = torch.polar(U(m,0.02,0.7),U(m,-math.pi,math.pi))
    y = torch.polar(U(m,0.02,0.7),U(m,-math.pi,math.pi))
    keep = C.mu_chart(x,y).abs() <= 0.55
    x,y = x[keep],y[keep]
    dx,dy = x*complex(0.3,-0.7),y*complex(-0.2,0.4)
    av = L.connection_chart(C,x,y,dx,dy)
    rep(f'   {C.name}: chart and torus connection agree ({len(x)})',
        rel(av,L.connection_torus(C,y,x*y,dy,dx*y+x*dy)),1e-11)
    # Arbitrary transverse vectors at the central axes; no log or 1/t is used.
    rr = torch.linspace(0.01,0.9,64)
    xa = torch.cat([rr/math.sqrt(C.a0),torch.zeros_like(rr)]).to(torch.complex128)
    ya = torch.cat([torch.zeros_like(rr),rr]).to(torch.complex128)
    vals = L.connection_chart(C,xa,ya,torch.full_like(xa,complex(0.2,-0.3)),torch.full_like(ya,complex(-0.4,0.1)))
    assert torch.isfinite(vals).all()
    rep(f'   {C.name}: central-axis circle generator (128)',
        (L.connection_chart(C,xa,ya,-1j*xa,1j*ya)-1).abs().max().item(),1e-12)
    xc,yc = x[:128]*0.01,y[:128]*0.01
    dxc,dyc = torch.full_like(xc,complex(0.2,-0.3)),torch.full_like(yc,complex(-0.4,0.1))
    expected = (-C.a0*xc.conj()*dxc+yc.conj()*dyc).imag/(C.a0*T.abs2(xc)+T.abs2(yc))
    rep(f'   {C.name}: exactly Hopf connection near pole',rel(L.connection_chart(C,xc,yc,dxc,dyc),expected),1e-12)

# ---------------------------------------------------------------- 5. connection gluing and quadrature gauge covariance
from torch.func import jvp
dx1,dx4 = torch.full_like(x1,0.13),torch.full_like(x4,-0.27)
# Use the independent middle sample from section 1, not the cusp loop's t.
tm = torch.polar(r,th)
dm = tm*complex(0.2,-0.3)
am = L.connection_middle(tm,x1,x4,lg,dm,dx1,dx4)
rep('5. global connection: log branch including tangent change',
    rel(L.connection_middle(tm,x1-x4,x4,lg+2j*math.pi,dm,dx1-dx4,dx4),am),1e-11)
rep('   global connection: x4 deck translation',
    rel(L.connection_middle(tm,x1,x4+1,lg,dm,dx1,dx4),am),1e-11)
for C,toMid in ((T.CUSP0,chart0_to_mid),(T.CUSPINF,chartinf_to_mid)):
    xx = torch.complex(U(512,.05,.3),U(512,-.1,.1))
    yy = torch.complex(U(512,.05,.3),U(512,-.1,.1))
    ddx,ddy = xx*complex(.2,-.3),yy*complex(-.1,.4)
    vals,tans = jvp(toMid,(xx,yy),(ddx,ddy))
    tt,one,four,lt = vals
    dtt,done,dfour,_ = tans
    rep(f'   {C.name}: global middle connection equals cusp form',
        rel(L.connection_middle(tt,one,four,lt,dtt,done,dfour),L.connection_chart(C,xx,yy,ddx,ddy)),1e-11)

import numpy as np
nodes,weights = np.polynomial.legendre.leggauss(96)
def quad(fun,a,b):
    u = torch.tensor((nodes+1)*(b-a)/2+a)
    return (fun(u)*torch.tensor(weights)).sum()*(b-a)/2

for C in (T.CUSP0,T.CUSPINF):
    for axis in (0,1):
        def curve(u):
            xx = (0.5+0.1*u+0.2j*torch.sin(2*u))/math.sqrt(C.a0)
            yy = 0.15*u+0.2j*torch.sin(u)
            return (xx,yy) if axis == 0 else (yy/math.sqrt(C.a0),xx*math.sqrt(C.a0))
        def form(u,gauged=False):
            (xx,yy),(ddx,ddy) = jvp(curve,(u,),(torch.ones_like(u),))
            if gauged:
                phase = 0.4*torch.sin(3*u)+2*u
                dphase = 1.2*torch.cos(3*u)+2
                dxg = torch.exp(-1j*phase)*(ddx-1j*dphase*xx)
                dyg = torch.exp(1j*phase)*(ddy+1j*dphase*yy)
                xx,yy = torch.exp(-1j*phase)*xx,torch.exp(1j*phase)*yy
                ddx,ddy = dxg,dyg
            return L.connection_chart(C,xx,yy,ddx,ddy)
        total = quad(form,-.3,.3)
        split = quad(form,-.3,.07)+quad(form,.07,.3)
        reverse = quad(form,.3,-.3)
        gauge = quad(lambda u: form(u,True),-.3,.3)
        phase_jump = 0.8*math.sin(.9)+1.2
        err = max(abs(total-split).item(),abs(total+reverse).item(),abs(gauge-total-phase_jump).item())
        rep(f'   {C.name}: transport composition/reversal/gauge across axis {axis}',err,1e-11)

# A phase map with range larger than 2*pi checks (8.8), including both poles.
def sphere_map(a,north):
    r2 = (a*a).sum(-1)
    return torch.stack([2*a[:,0],2*a[:,1],(1-r2) if north else (r2-1)],-1)/(1+r2)[:,None]
def phase_on_sphere(n):
    return 7*n[:,0]+4*n[:,1]*n[:,2]+5*n[:,2]**2
def primitive(a,north):
    # eta=Im(conj(g)dg), evaluated through g, rather than differentiating a
    # chosen logarithm. The known phase is only the independent comparison.
    def integrand(u):
        b = u[:,None,None]*a[None,:,:]
        d = a[None,:,:].expand_as(b)
        shape = b.shape
        def gg(aa):
            return torch.exp(1j*phase_on_sphere(sphere_map(aa,north)))
        val,der = jvp(gg,(b.reshape(-1,2),),(d.reshape(-1,2),))
        return (val.conj()*der).imag.reshape(shape[:2])
    un = torch.tensor((nodes+1)/2)
    return (integrand(un)*torch.tensor(weights/2)[:,None]).sum(0)
aa = torch.randn(128,2,generator=g)
aa = aa/(1+aa.norm(dim=-1,keepdim=True))
aa[0] = 0  # The two charts evaluate the north and south poles here.
eq = torch.tensor([[1.,0.]])
theta0 = math.atan2(math.sin(5),math.cos(5))
thetaN = theta0+primitive(aa,True)
thetaS = theta0+primitive(eq,True)-primitive(eq,False)+primitive(aa,False)
err = max((torch.exp(1j*thetaN)-torch.exp(1j*phase_on_sphere(sphere_map(aa,True)))).abs().max().item(),
          (torch.exp(1j*thetaS)-torch.exp(1j*phase_on_sphere(sphere_map(aa,False)))).abs().max().item())
rep('5. two-chart real logarithm: 256 samples, phase range >2pi',err,1e-11)
assert torch.cat([thetaN,thetaS]).max()-torch.cat([thetaN,thetaS]).min() > 2*math.pi
angles = torch.linspace(0,2*math.pi,64)
equator = torch.stack([torch.cos(angles),torch.sin(angles)],-1)
real_overlap = primitive(equator,True)-primitive(eq,True)+primitive(eq,False)-primitive(equator,False)
rep('   real logarithms agree (not only their exponentials): 64 overlaps',real_overlap.abs().max().item(),1e-11)
# 6. Global torus-gauge change used by the genuinely elementary pole tensor.
for C,to_mid in ((T.CUSP0,chart0_to_mid),(T.CUSPINF,chartinf_to_mid)):
    z = torch.randn(1024,4,generator=g)
    radius = torch.linspace(1e-6,.01,len(z))
    z *= radius.sqrt()[:,None]/z.norm(dim=-1,keepdim=True)
    x = torch.complex(z[:,0],z[:,1])/math.sqrt(C.a0)
    y = torch.complex(z[:,2],z[:,3])
    t,x1,x4,lg = to_mid(x,y)
    bc = L.elementary_torus_shift_chart(C,x,y)
    bm = L.elementary_torus_shift_middle(t,x1,x4,lg)
    label = 'inf' if C is T.CUSPINF else '0'
    rep(f'6. torus correction cusp {label}: 1024 chart/middle comparisons',(bc-bm).abs().max().item(),1e-12)
    shifted = L.elementary_torus_shift_middle(t,x1-x4,x4,lg+2j*math.pi)
    deck = L.elementary_torus_shift_middle(t,x1,x4+1,lg)
    rep(f'   correction branch/deck invariance, cusp {label}',
        max((shifted-bm).abs().max().item(),(deck-bm).abs().max().item()),1e-12)
    options = dict(algebraic_radial=True,algebraic_profile=True)
    vc,hc = L.transport_base_chart(C,x,y,**options)
    vm,hm = L.transport_base_middle(t,x1,x4,lg,**options)
    vb,hb = L.transport_base_middle(t,x1-x4,x4,lg+2j*math.pi,**options)
    vd,hd = L.transport_base_middle(t,x1,x4+1,lg,**options)
    rep(f'   revised base chart/middle/branch/deck, cusp {label}',
        max((vc-vm).abs().max().item(),(hc-hm).abs().max().item(),
            (vb-vm).abs().max().item(),(hb-hm).abs().max().item(),
            (vd-vm).abs().max().item(),(hd-hm).abs().max().item()),1e-11)
    inner = radius<L.EPS1
    H = C.mu_chart(x[inner],y[inner])
    c2,c3 = C.c_tilde(x[inner]*y[inner])
    rep(f'   exact cancellation b=H*c on the cap, cusp {label}',
        (bc[inner]-H[:,None]*torch.stack([c2,c3],-1)).abs().max().item(),1e-14)
    assert (bc[radius>=L.EPS2] == 0).all()
    axis = torch.linspace(0,.008,128).sqrt().to(torch.complex128)
    zero = torch.zeros_like(axis)
    assert (L.elementary_torus_shift_chart(C,axis/math.sqrt(C.a0),zero) == 0).all()
    assert (L.elementary_torus_shift_chart(C,zero,axis) == 0).all()
    print(f'   cusp {label}: support and 256 exact pole/axis values: PASS')
# 7. Finite height: full cusp collars, base/connection/torus-gauge gluing.
for C,to_mid in ((T.FINITE_CUSP0,chart0_to_mid),(T.FINITE_CUSPINF,chartinf_to_mid)):
    tc = torch.polar(torch.exp(U(2048,math.log(.003),math.log(.5))),U(2048,-math.pi,math.pi))
    zeta = torch.exp(C.P(tc)*U(2048,-2,3)+1j*U(2048,-math.pi,math.pi))
    _,xx,yy = C.reduce(zeta,tc)
    t,one,four,lt = to_mid(xx,yy)
    opts = dict(algebraic_radial=True,algebraic_profile=True)
    vc,hc = L.transport_base_chart(C,xx,yy,**opts)
    vm,hm = L.transport_base_middle(t,one,four,lt,finite_height=True,**opts)
    vd,hd = L.transport_base_middle(t,one,four+1,lt,finite_height=True,**opts)
    vb,hb = L.transport_base_middle(t,one-four,four,lt+2j*math.pi,finite_height=True,**opts)
    rep(f'7. finite {C.name}: 2048 base chart/middle/deck/branch comparisons',
        max(rel(vc,vm),rel(hc,hm),rel(vd,vm),rel(hd,hm),rel(vb,vm),rel(hb,hm)),1e-11)
    dx,dy = xx*complex(.13,.21),yy*complex(-.19,.17)
    vals,tans = jvp(to_mid,(xx,yy),(dx,dy))
    tt,oo,ff,ll = vals
    dt,do,df,_ = tans
    ac = L.connection_chart(C,xx,yy,dx,dy)
    am = L.connection_middle(tt,oo,ff,ll,dt,do,df,finite_height=True)
    # On .3<|tc|<.5 the global connection also has its explicit middle term.
    pure = tc.abs()<=.3
    rep(f'   finite {C.name}: connection chart/middle',rel(ac[pure],am[pure]),1e-11)
    rep(f'   finite {C.name}: global connection branch/deck',
        max(rel(am,L.connection_middle(tt,oo-ff,ff,ll+2j*math.pi,dt,do-df,df,finite_height=True)),
            rel(am,L.connection_middle(tt,oo,ff+1,ll,dt,do,df,finite_height=True))),1e-11)
    rep(f'   finite {C.name}: torus correction chart/middle',
        rel(L.elementary_torus_shift_chart(C,xx,yy),
            L.elementary_torus_shift_middle(tt,oo,ff,ll,finite_height=True)),1e-12)
    # A genuine derivative check at each artificial logarithmic representative seam.
    ts = torch.polar(torch.linspace(.01,.5,81),torch.full((81,),.37))
    rs = torch.cat([k*C.P(ts) for k in (-2,-1,0,1,2)])
    ts = ts.repeat(5)
    def torus_h(rho):
        return C.mu_torus(torch.exp(rho).to(T.CD),ts)
    def chart_h(rho):
        yy = torch.exp(rho).to(T.CD)
        return C.mu_chart(ts/yy,yy)
    fa,fb = torus_h,chart_h
    for degree in range(3):
        rep(f'   finite {C.name}: 405 deck-seam jets, order {degree}',rel(fa(rs),fb(rs)),1e-11)
        pa,pb = fa,fb
        fa = lambda u,pa=pa:jvp(pa,(u,),(torch.ones_like(u),))[1]
        fb = lambda u,pb=pb:jvp(pb,(u,),(torch.ones_like(u),))[1]

# 8. Independent algebraic round-S4 fibration: values, boost inverse and ranks.
qq = torch.randn(4096,5,generator=g)
qq /= qq.norm(dim=-1,keepdim=True)
ss = L.suspended_hopf(qq)
kk = L.sphere3_boost(ss)
ff = L.algebraic_torus_fibration(qq)
rep('8. algebraic fibration: 4096 sphere identities and boost inverses',
    max((ss.square().sum(-1)-1).abs().max().item(),(kk.square().sum(-1)-1).abs().max().item(),
        (ff.square().sum(-1)-1).abs().max().item(),
        (L.sphere3_boost(kk,True)-ss).abs().max().item()),1e-12)
df = torch.stack([jvp(L.algebraic_torus_fibration,(qq,),
                      (torch.eye(5)[i:i+1].expand_as(qq),))[1] for i in range(5)],-1)
proj = torch.eye(5)[None]-qq[:,:,None]*qq[:,None,:]
sv = torch.linalg.svdvals(df@proj)
margin = (sv[:,1]/torch.sqrt(1-qq[:,4]**2)).min().item()
assert margin>1 and sv[:,2].max()<1e-12
print(f'   4096 tangent-rank checks: min sigma2/sqrt(1-X^2)={margin:.6f}>1: PASS')
pp = torch.tensor([[0.,0.,0.,0.,-1.],[0.,0.,0.,0.,1.]])
assert torch.allclose(L.algebraic_torus_fibration(pp),torch.tensor([[0.,1.,0.],[0.,-1.,0.]]),atol=1e-14,rtol=0)
for i in range(4):
    assert (jvp(L.algebraic_torus_fibration,(pp,),
                (torch.eye(5)[i:i+1].expand_as(pp),))[1] == 0).all()
print('   both exact poles: distinct critical values and zero tangential differential: PASS')

# 9. Finite logarithm: a known smooth phase with range exceeding three full turns.
# ||D phase||<=12+9+5=26 follows term by term on the unit sphere.
def phase9(n):
    return 12*n[...,0]+9*n[...,1]*n[...,2]+5*n[...,2]
def circle9(n):
    return torch.exp(1j*phase9(n))
nn = torch.randn(256,3,generator=g)
nn /= nn.norm(dim=-1,keepdim=True)
nn = torch.cat([nn,torch.tensor([[0.,0.,1.],[0.,0.,-1.],[1.,0.,0.]])])
theta9 = L.finite_phase_log(circle9,nn,26)
shift9 = 5-math.atan2(math.sin(5),math.cos(5))
rep('9. finite phase: 259 exact real values, including both poles',
    (theta9-(phase9(nn)-shift9)).abs().max().item(),1e-12)
rep('   finite phase: 259 exponentiation identities',
    (torch.exp(1j*theta9)-circle9(nn)).abs().max().item(),1e-12)
eq9 = torch.stack([torch.cos(torch.linspace(-math.pi,math.pi,65)),
                   torch.sin(torch.linspace(-math.pi,math.pi,65)),torch.zeros(65)],-1)
rep('   finite phase: 65 equator overlaps as REAL logarithms',
    (L.finite_phase_log(circle9,eq9,26,'north')-
     L.finite_phase_log(circle9,eq9,26,'south')).abs().max().item(),1e-12)
# Equatorial tangent jets approach from the two independently defined charts.
v9 = torch.zeros_like(eq9); v9[:,2] = 1
def chart9(n,which):
    return L.finite_phase_log(circle9,n/n.norm(dim=-1,keepdim=True),26,which)
dn9 = jvp(lambda n:chart9(n,'north'),(eq9,),(v9,))[1]
ds9 = jvp(lambda n:chart9(n,'south'),(eq9,),(v9,))[1]
ref9 = jvp(phase9,(eq9,),(v9,))[1]
rep('   finite phase: 65 equator tangent derivatives',
    max((dn9-ds9).abs().max().item(),(dn9-ref9).abs().max().item()),1e-11)

# Rank-one projector chain: all distances are checked, and reversal is an exact
# inverse. This is a test of the finite operation, not a bundle embedding of W.
times9 = torch.linspace(0,1,65)
aa9,bb9 = U(128,-1,1),U(128,-2,2)
vectors9 = torch.stack([torch.cos(aa9[:,None]*times9).to(T.CD),
                       torch.sin(aa9[:,None]*times9)*torch.exp(1j*bb9[:,None]*times9)],-1)
projectors9 = vectors9[:,:,:,None]*vectors9[:,:,None,:].conj()
steps9 = torch.linalg.matrix_norm(projectors9[:,1:]-projectors9[:,:-1],ord=2)
assert steps9.max()<.05
start9 = vectors9[:,0]
end9 = L.finite_projection_lift(projectors9,start9)
back9 = L.finite_projection_lift(projectors9.flip(1),end9)
rot9 = torch.exp(1j*U(128,-math.pi,math.pi))
rotend9 = L.finite_projection_lift(projectors9,rot9[:,None]*start9)
rep('   finite projections: 128 chains, reversal, unit norm, image and equivariance',
    max((back9-start9).abs().max().item(),(end9.norm(dim=-1)-1).abs().max().item(),
        ((projectors9[:,-1]@end9[:,:,None])[:,:,0]-end9).abs().max().item(),
        (rotend9-rot9[:,None]*end9).abs().max().item()),1e-12)

# 10. Concrete ten-component bundle embedding, with independent gauges.
for finite10 in (False,True):
    tt = torch.polar(torch.exp(U(4096,-10,10)),U(4096,-math.pi,math.pi))
    lt = torch.log(tt)
    one,four = U(4096,-2,2),U(4096,-3,3)
    e10,w10 = L.bundle_embedding_middle(tt,one,four,lt,finite_height=finite10,with_norm=True)
    assert w10.min()>=.5-1e-14 and w10.max()<=2+1e-14
    ed10 = L.bundle_embedding_middle(tt,one,four+1,lt,finite_height=finite10)
    eb10 = L.bundle_embedding_middle(tt,one-four,four,lt+2j*math.pi,finite_height=finite10)
    er10 = L.bundle_embedding_middle(tt,one+.137,four,lt,finite_height=finite10)
    rep(f'10. bundle finite={finite10}: 4096 norms/deck/branch/equivariance tests',
        max((T.abs2(e10).sum(-1)-1).abs().max().item(),rel(e10,ed10),rel(e10,eb10),
            rel(er10,e10*complex(math.cos(2*math.pi*.137),math.sin(2*math.pi*.137)))),1e-11)
    for C,to_mid in zip(T.height_cusps(finite10),(chart0_to_mid,chartinf_to_mid)):
        tc = torch.polar(torch.exp(U(1024,-12,math.log(.29))),U(1024,-math.pi,math.pi))
        zz = torch.exp(C.P(tc)*U(1024,-2,3)+1j*U(1024,-math.pi,math.pi))
        _,xx,yy = C.reduce(zz,tc)
        args = to_mid(xx,yy)
        ec = L.bundle_embedding_chart(C,xx,yy)
        em = L.bundle_embedding_middle(*args,finite_height=finite10)
        rep(f'    bundle {C.name}: 1024 chart/middle comparisons',rel(ec,em),1e-11)
        rr = torch.linspace(.0001,.8,128)
        ax = torch.cat([rr.sqrt().to(T.CD)/math.sqrt(C.a0),torch.zeros(128,dtype=T.CD)])
        ay = torch.cat([torch.zeros(128,dtype=T.CD),rr.sqrt().to(T.CD)])
        ac = L.bundle_embedding_chart(C,ax,ay)
        assert torch.isfinite(ac).all()
        assert (T.abs2(ac).sum(-1)-1).abs().max()<1e-14
        # Independent finite differences through either central axis.
        dx,dy = torch.full_like(ax,.13+.07j),torch.full_like(ay,-.11+.19j)
        da = jvp(lambda x,y:L.bundle_embedding_chart(C,x,y),(ax,ay),(dx,dy))[1]
        h10 = 1e-6
        fd = (L.bundle_embedding_chart(C,ax+h10*dx,ay+h10*dy)
              -L.bundle_embedding_chart(C,ax-h10*dx,ay-h10*dy))/(2*h10)
        rep(f'    bundle {C.name}: 256 exact-axis tangent comparisons',rel(da,fd),1e-6)
        # Exact homogeneous pole block, including points on either axis.
        rad10 = torch.linspace(.00001,.004,128)
        frac10 = torch.linspace(0,1,128)
        xp = torch.sqrt(rad10*frac10/C.a0).to(T.CD)
        yp = torch.sqrt(rad10*(1-frac10)).to(T.CD)
        ep = L.bundle_embedding_chart(C,xp,yp)
        expected = torch.zeros_like(ep)
        offset10 = 4 if C.name=='inf' else 0
        expected[:,offset10:offset10+2] = torch.stack([math.sqrt(C.a0)*xp.conj(),yp],-1)/rad10.sqrt()[:,None]
        rep(f'    bundle {C.name}: 128 exact homogeneous pole models',rel(ep,expected),1e-13)

print('ALL PASS' if not fails else f'FAILURES: {fails}')
assert not fails, fails
