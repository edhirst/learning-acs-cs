"""Symbolic certificate and numerical margins for I2 (programme section 4.3).
The certificate reduces both cusps to one rational-trigonometric map and proves
|D Lambda-M| <= Cprime*rad, |Lambda-M*vq| <= Cprime*rad^2/2,
with Cprime=178 at zero and 89 at infinity (including the k transition).
The rational endpoint comparisons below implement the written type-(c) proof.
The separate numerical measurements on the
linear-blend shell EPS1 <= rad <= EPS2 of each cusp, measure
   C_L  = max |Lambda(vq) - M vq| / |vq|^2          (Lambda = the unblended P_final o moebius o R o h_B, in vq),
   C_L' = max |D Lambda(vq) - M| / |vq|             (Jacobian in vq by autograd through the chart point),
and a sample-based estimate for the smallest singular value of D Lambda_mod from the identity N4:
   sigma_min >= scale - c |D Lambda - M| - |grad c| |Lambda - M vq|,   |grad c| <= 2/(EPS2 - EPS1) = 500.
Also the actual sigma_min of D Lambda_mod (autograd) on the shell.  Sanity for the explicit-constant certificate.
Run from the repo root: python3 extra/s4t2_explicit/sym/i2_margins.py
"""
import math, os, sys
import torch
import sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
torch.set_default_dtype(torch.float64)
from extra.s4t2_explicit import tcore as T, tlift as L
from torch.func import jvp

EYE4 = torch.eye(4)


def exact(label, condition):
    assert condition, label
    print(label + ': PASS', flush=True)


# I2a: local orbit map, with normalized coordinates (a,b,c)=2*pi*vq.
a, b, c = sp.symbols('a b c', real=True)
d = sp.sqrt(1 + a*a + b*b)
p = a*sp.cos(c) + b*sp.sin(c)
q = a*sp.sin(c) - b*sp.cos(c)
D = d + p
N = sp.Matrix([q, -sp.sin(c), D-sp.cos(c)])
B = N / D
Q = sp.Matrix([[0,-1,0], [0,0,-1], [1,0,0]])
Oin = sp.diag(1,-1,1)
Oout = sp.Matrix([[0,-1,0], [-1,0,0], [0,0,1]])
u0 = sp.Matrix([p,q,sp.cos(c),sp.sin(c)]) / d
ui = sp.Matrix([-sp.cos(c),sp.sin(c),-a*sp.cos(c)+b*sp.sin(c),a*sp.sin(c)+b*sp.cos(c)]) / d
rho = lambda w: sp.Matrix([w[0],-w[2],-w[1]])
for name, u in (('zero',u0),('infinity',ui)):
    exact('I2a '+name+' sphere parametrization', sp.trigsimp((u.T*u)[0]-1) == 0)
stereo = sp.Matrix([u0[1],u0[2]-1-u0[0],u0[3]])/(1+u0[0])
exact('I2b cusp zero is B', all(sp.simplify(x) == 0 for x in rho(stereo)-B))
P0,P1,P2,P3 = sp.symbols('P0 P1 P2 P3', real=True)
numer = sp.Matrix([P1,P2-1-P0,P3])
exact('I2c inverted shifted stereographic denominator',
      sp.expand((numer.T*numer)[0]-2*(1-P2)*(1+P0)) == P0**2+P1**2+P2**2+P3**2-1)
invstereo = sp.Matrix([ui[1],ui[2]-1-ui[0],ui[3]])/(2*(1-ui[2]))
exact('I2d cusp infinity is Oout B(Oin x)/2',
      all(sp.simplify(x) == 0 for x in rho(invstereo)-Oout*B.subs(b,-b)/2))
variables = sp.Matrix([a,b,c])
DB = B.jacobian(variables)
exact('I2e quotient derivative and linear model',
      all(sp.simplify(x) == 0 for x in DB-N.jacobian(variables)/D+N*sp.Matrix([D]).jacobian(variables)/D**2)
      and DB.subs({a:0,b:0,c:0}) == Q)
Rot = sp.Matrix([[sp.sin(c),-sp.cos(c)],[sp.cos(c),sp.sin(c)]])
Rot0 = Rot.subs(c,0)
exact('I2f rotation and gradient norms used in the bound',
      sp.simplify((Rot-Rot0).T*(Rot-Rot0)) == 2*(1-sp.cos(c))*sp.eye(2)
      and sp.trigsimp(sum(sp.diff(p,v)**2 for v in variables)-1-q*q) == 0
      and sp.simplify(d-1-(a*a+b*b)/(d+1)) == 0)

# Rational elementary constants. For x>0 and n+2>x, the positive exponential
# series tail is at most its first term times 1/(1-x/(n+2)).
def exp_upper(x, n=40):
    return sum(x**k/sp.factorial(k) for k in range(n+1)) + x**(n+1)/sp.factorial(n+1)/(1-x/(n+2))


def atan_bounds(x):
    upper = x-x**3/3+x**5/5
    return upper-x**7/7, upper


lo5,hi5 = atan_bounds(sp.Rational(1,5))
lo239,hi239 = atan_bounds(sp.Rational(1,239))
tan2 = lambda t: 2*t/(1-t*t)
tan4 = tan2(tan2(sp.Rational(1,5)))
exact('I2g Machin identity and 157/50 < pi < 22/7',
      (tan4-sp.Rational(1,239))/(1+tan4/sp.Integer(239)) == 1
      and 16*lo5-4*hi239 > sp.Rational(157,50)
      and 16*hi5-4*lo239 < sp.Rational(22,7))
exact('I2h exponential bounds for the profile domains',
      exp_upper(sp.Integer(1)) < 3 and exp_upper(sp.Integer(2)) < 8
      and exp_upper(sp.Integer(4)) < 55 and exp_upper(sp.Integer(6)) < 500
      and exp_upper(sp.Rational(11,2)) < 250
      and sum(sp.Rational(224,15)**k/sp.factorial(k) for k in range(65)) > 3000000
      and sum(sp.Rational(13,5)**k/sp.factorial(k) for k in range(11)) > sp.Rational(88,7))
S = sp.Rational(51,1000)
factor = (2+S/2)/(1-S)+(1+S/2)/(1-S)+(1+S+S*S/6)*(1+S+S*S/2)/(1-S)**2
exact('I2i |DB-Q| <= (9/2)|x| on |x|<=51/1000', factor < sp.Rational(9,2))
exact('I2i rational decimal enclosure used in the proof', factor < sp.Rational(4443503,1000000))
exact('I2j radial zones and positive denominator',
      2*sp.Rational(22,7)*sp.Rational(1,125) < S
      and (S+S*S+S**3/6)/(1-S) < sp.Rational(1,16)
      and sp.Rational(1,16) < sp.Rational(1,3)
      and sp.Rational(250,8) > 26 and 55 < 250)
u = sp.Symbol('u', positive=True)
exact('I2k monotonic exponential-tail factor',
      sp.simplify(sp.diff(sp.exp(-1/u)/u**2,u)-sp.exp(-1/u)*(1-2*u)/u**4) == 0)
delta = sp.Rational(21,500000)
exact('I2l |1-beta|<=1/500000, |d beta/d log rho|<=1/25000',
      sp.Rational(16,5)/3000000 < sp.Rational(1,500000)
      and sp.Rational(2,5)*258/3000000 < sp.Rational(1,25000)
      and sp.Rational(1,500000)+sp.Rational(1,25000) == delta)
pi_hi = sp.Rational(22,7)
exact('I2m derivative constants: 178 and 89',
      18*pi_hi**2 < 178 and 9*pi_hi**2*(1+delta)+250*pi_hi*delta < 89)
R = sp.Rational(1,125)
exact('I2n blended singular-value margins exceed 2 and 1',
      2*sp.Rational(157,50)-178*R-500*89*R**2 > 2
      and sp.Rational(157,50)-89*R-500*sp.Rational(89,2)*R**2 > 1)
pi_lo = 16*lo5-4*hi239
L0 = 178*R+500*89*R**2
exact('I2o projected-defect contraction ratio < 17/25 at both cusps',
      L0/(2*pi_lo) < sp.Rational(17,25)
      and (89*R+500*sp.Rational(89,2)*R**2)/pi_lo < sp.Rational(17,25))
exact('I2p epsilon=1/100 lies strictly inside both linear caps',
      sp.Rational(1,100) < pi_lo*sp.Rational(1,250))

# Section 9: balanced algebraic profile; the original bounds above are retained.
exact('I2q algebraic beta: 0<a<b<40; beta>40/41',
      55 < (16*sp.Rational(157,50)**2)**2
      and 4*pi_hi**2 < 40
      and 1/(1+40*R**2) > sp.Rational(40,41)**2)
exact('I2r new raw derivative constants 182 and 91',
      18*pi_hi**2+162*pi_hi*R < 182)
new_lip = 182*R+500*91*R**2
exact('I2s algebraic-profile shell: contraction<7/10, singular margins>19/10 and 19/20',
      new_lip/(2*pi_lo) < sp.Rational(7,10)
      and 2*pi_lo-new_lip > sp.Rational(19,10)
      and pi_lo-new_lip/2 > sp.Rational(19,20))


def reduced_orbit(vq, infinity=False):
    """Closed local formula for the unblended Lambda on rad<=0.008.
    No theta series, fitted phase, stereographic singularity or radial inversion.
    """
    xx = 2*math.pi*vq
    if infinity:
        rho2 = vq[...,0]**2+vq[...,1]**2
        outside = rho2 > (2*math.exp(-6))**2
        rho = torch.sqrt(torch.where(outside,rho2,torch.ones_like(rho2)))
        rr = 2*math.exp(-2)/rho
        beta = torch.where(outside,T.KINF/T.k_of(rr),torch.ones_like(rho))
        xx = xx*torch.stack([beta,-beta,torch.ones_like(beta)],-1)
    aa,bb,cc = xx.unbind(-1)
    dd = torch.sqrt(1+aa*aa+bb*bb)+aa*torch.cos(cc)+bb*torch.sin(cc)
    out = torch.stack([aa*torch.sin(cc)-bb*torch.cos(cc),-torch.sin(cc),dd-torch.cos(cc)],-1)/dd[...,None]
    return torch.stack([-out[...,1],-out[...,0],out[...,2]],-1)/2 if infinity else out


def jac(f, p):
    """Batched Jacobian (n, out, 4) of f: R^4 -> R^out by forward-mode jvp along the basis (complex-safe)."""
    return torch.stack([jvp(f, (p,), (EYE4[k].expand_as(p),))[1] for k in range(4)], -1)


torch.manual_seed(5)
n = 4000
for name, C, hb, lam_chart, M, scale, invr in (('cusp 0', T.CUSP0, L.hB_chart0, L.Lam_chart0, L.M0, L.S0, False),
                                               ('cusp inf', T.CUSPINF, L.hB_chartinf, L.Lam_chartinf, L.MINF, L.SINF, True)):
    rad = torch.empty(n).uniform_(L.EPS1, L.EPS2)
    frac = torch.empty(n).uniform_(0, 1)
    a, b = torch.empty(n).uniform_(-math.pi, math.pi), torch.empty(n).uniform_(-math.pi, math.pi)
    x = torch.polar(torch.sqrt(rad * frac / C.a0), a)
    y = torch.polar(torch.sqrt(rad * (1 - frac)), b)

    def lam_unblended(x, y):                       # v of Lambda (no blend) as a function of the chart point
        r = (x * y).abs()
        v, X = L.Lam_from_u(hb(x, y), 1 / r if invr else r)
        return v

    def vq_of(x, y):
        return L.v_q(C, x, y)

    def real4(x, y):
        return torch.stack([x.real, x.imag, y.real, y.imag], -1)

    def from_real4(p):
        return torch.complex(p[..., 0], p[..., 1]), torch.complex(p[..., 2], p[..., 3])

    p = real4(x, y)
    v = lam_unblended(x, y)
    vq = vq_of(x, y)
    err = (v - vq @ M.T).norm(dim=-1)
    nq = vq.norm(dim=-1)
    CL = (err / nq ** 2).max().item()
    # Jacobians in vq: D_p v (3x4) = D_vq v (3x3) . D_p vq (3x4); D_p vq has rank 3 off the origin, so D_vq v = D_p v . pinv(D_p vq)
    Jv = jac(lambda pp: lam_unblended(*from_real4(pp)), p)
    Jq = jac(lambda pp: vq_of(*from_real4(pp)), p)
    DL = Jv @ torch.linalg.pinv(Jq)
    reduced = lambda pp: reduced_orbit(vq_of(*from_real4(pp)), invr)
    reduced_error = (reduced(p)-v).abs().max().item()
    reduced_derivative_error = (jac(reduced,p)-Jv).abs().max().item()
    assert reduced_error < 1e-12 and reduced_derivative_error < 1e-10
    new_raw = L.transport_base_local(vq,invr,modified=False)
    new_blend = lambda pp: L.transport_base_local(vq_of(*from_real4(pp)),invr)
    new_error = max((new_raw-reduced(p)).abs().max().item(),
                    (new_blend(p)-lam_chart(x,y)[0]).abs().max().item())
    new_derivative_error = (jac(new_blend,p)-jac(lambda pp:lam_chart(*from_real4(pp))[0],p)).abs().max().item()
    assert new_error < 1e-12 and new_derivative_error < 1e-10
    dev = torch.linalg.matrix_norm(DL - M, ord=2)
    CLp = (dev / nq).max().item()
    # the blended map and its smallest singular value in vq
    Jm = jac(lambda pp: lam_chart(*from_real4(pp))[0], p)
    DLm = Jm @ torch.linalg.pinv(Jq)
    smin = torch.linalg.svdvals(DLm)[:, -1]
    contraction_sample = (torch.linalg.matrix_norm(DLm-M,ord=2)/scale).max().item()
    assert contraction_sample < 17/25
    bound = scale - dev.max().item() - 500 * err.max().item()
    print(f'{name}: shell rad in [{L.EPS1}, {L.EPS2}], {n} points')
    print(f'   max |Lambda - M vq| / |vq|^2 = {CL:.3f};  max |D Lambda - M|_2 / |vq| = {CLp:.3f};  scale = {scale:.4f}')
    print(f'   max |Lambda - M vq| = {err.max():.2e}, max |D Lambda - M|_2 = {dev.max():.3e}')
    print(f'   sample-based lower estimate (NOT a certificate) = {bound:.4f};  measured min sigma_min = {smin.min():.4f}')
    print(f'   reduced formula vs original: value {reduced_error:.1e}, derivative {reduced_derivative_error:.1e}; '
          f'proven singular-value lower bound: {1 if invr else 2}')
    print(f'   measured defect derivative / scale {contraction_sample:.6f}; proven contraction bound < 0.68')
    print(f'   new base vs independent old shell: value {new_error:.1e}, derivative {new_derivative_error:.1e}')

# The infinity shell is not entirely on the constant-k plateau. At its equator,
# |s|=EPS2*exp(2)/2 and |t|=250*exp(-2)<exp(4), since exp(6)>250.
# A positive Taylor partial sum proves the last strict inequality.
assert sum(sp.Integer(6) ** k / sp.factorial(k) for k in range(8)) > 250
s_equator = torch.tensor(L.EPS2 / (2 * math.sqrt(T.CUSPINF.a0)))
t_equator = 1 / s_equator
print(f'I2 audit (one infinity-shell point): |s|={s_equator:.12f}, |t|={t_equator:.12f}, '
      f'k/KINF-1={T.k_of(t_equator) / T.KINF - 1:.3e}. The certificate includes this transition.')
# Independent composition check for the new profile, before radial compactification.
from extra.s4t2_explicit import tinverse as TI
for C,infinity in ((T.CUSP0,False),(T.CUSPINF,True)):
    scale = L.SINF if infinity else L.S0
    M = L.MINF if infinity else L.M0
    z = torch.randn(4000,4)
    z *= torch.sqrt(torch.linspace(L.EPS1,L.EPS2,len(z)))[:,None]/z.norm(dim=-1,keepdim=True)
    z[:,0:2] /= math.sqrt(C.a0)
    def orbit(z):
        x,y = torch.complex(z[:,0],z[:,1]),torch.complex(z[:,2],z[:,3])
        return L.v_q(C,x,y)
    def direct(z):
        x,y = torch.complex(z[:,0],z[:,1]),torch.complex(z[:,2],z[:,3])
        m = C.a0*T.abs2(x)-T.abs2(y)
        u = L.algebraic_hopf_cusp(x*y,m,infinity)
        w = L.transport_stereographic(u)
        if infinity:
            w = w/(w*w).sum(-1)[:,None]
        return torch.stack([w[:,0],-w[:,2],-w[:,1]],-1)
    p = orbit(z)
    formula = lambda zz:L.transport_base_local(orbit(zz),infinity,modified=False,algebraic_profile=True)
    value_error = (direct(z)-formula(z)).abs().max().item()
    derivative_error = (jac(direct,z)-jac(formula,z)).abs().max().item()
    assert value_error < 1e-12 and derivative_error < 1e-9
    eye3 = torch.eye(3)
    Jnew = torch.stack([jvp(lambda pp:L.transport_base_local(pp,infinity,algebraic_profile=True),
                           (p,),(eye3[j].expand_as(p),))[1] for j in range(3)],-1)
    min_singular = torch.linalg.svdvals(Jnew)[:,-1].min().item()
    contraction = (torch.linalg.matrix_norm(Jnew-M,ord=2)/scale).max().item()
    assert min_singular > (19/20 if infinity else 19/10) and contraction < .7
    print(f'I2 algebraic profile cusp {"inf" if infinity else "0"}: 4000 points; '
          f'composition {value_error:.1e}, derivative {derivative_error:.1e}; '
          f'measured min singular {min_singular:.6f} vs proved {19/20 if infinity else 19/10}; '
          f'defect ratio {contraction:.6f} vs proved <.7')
rr = torch.logspace(-12,12,4097)
back = TI.algebraic_r_from_absT(rr*L.algebraic_k(rr))
profile_error = ((back-rr)/rr).abs().max().item()
assert profile_error < 1e-14
w = torch.randn(4097,3)
w *= torch.logspace(-12,12,4097)[:,None]/w.norm(dim=-1,keepdim=True)
v,h = L.P_algebraic(w)
radial_error = ((TI.P_algebraic_inv(v,h)-w).norm(dim=-1)/w.norm(dim=-1)).max().item()
assert radial_error < 1e-14 and (v.norm(dim=-1)+h*h-1).abs().max() < 1e-14
print(f'Algebraic blocks: 4097 targets each over 24 orders of magnitude; '
      f'profile relative inverse error {profile_error:.1e}, radial {radial_error:.1e}')
# Inverse-first pole certificate, programme 11.2. These are rational endpoint
# checks for the explicit inequalities derived there, not sampled bounds.
Rinv = sp.Rational(1,250)
qinv = Rinv**2/4
sinv = 55*qinv
assert 220/(1-sinv)<221 and 220/(1-sinv)**2<222
assert 1+222*qinv<sp.Rational(1001,1000)
assert 222*Rinv**2<sp.Rational(1,250)
lin_u = (sp.Rational(444,4)+sp.Rational(1001,1000))/(1-qinv)+35*sp.Rational(1001,1000)*Rinv**2/(1-qinv)**2
lin_v = (222+sp.Rational(1001,1000)+sp.Rational(1,4))/(1-qinv)+35*sp.Rational(1001,1000)*Rinv**2/(1-qinv)**2
assert max(lin_u,lin_v)<225
atan_u = sp.Rational(1001,1000)/(1-qinv**3)+sp.Rational(1001,1000)*Rinv**2*qinv**2/(1-qinv**3)**2
atan_v = (1+sp.Rational(1,250))/(1-qinv**3)+sp.Rational(1001,1000)*Rinv**2*qinv**2/(1-qinv**3)**2
assert max(atan_u,atan_v)<2
assert 1+227*Rinv**2<2 and 6/(1-qinv)**2<7 and 225+2+2*7<300
assert 38*Rinv**2<sp.Rational(1,1000)
assert 30*sp.Rational(1001,1000)*(1+700*Rinv**2)+700*Rinv<40
assert 15*sp.Rational(1001,1000)**2+38*Rinv<20
print('I2t inverse-first height: |R|<=38r^3, |R_u|,|R_v|<300r^2; combined defect <=20r^2, derivative <=40r: PASS')
rmax = sp.Rational(1,10000)
assert 2*rmax<Rinv and 40*Rinv<1
assert 80/(1-80*rmax)<100
amin,amax,hmin = 1-80*rmax,1+80*rmax,1-10000*rmax**2
assert 160/amin<200 and (100+400+80)/amin<600
cmax = sp.sqrt(amax/hmin)
cdiff = 80/((1+sp.sqrt(amin))*sp.sqrt(hmin))+10000*rmax/(sp.sqrt(hmin)*(1+sp.sqrt(hmin)))
assert cmax<sp.Rational(201,200) and cdiff<42
assert 1/(2*sp.sqrt(amin*hmin))<sp.Rational(51,100)
assert sp.sqrt(amax)/(2*hmin**sp.Rational(3,2))<sp.Rational(51,100)
assert sp.Rational(51,100)*(400+60000*rmax)<210
assert 42+sp.Rational(201,200)*100<200
assert 42+sp.Rational(201,200)*100+420*(1+100*rmax)+sp.Rational(201,200)*600<1500
epsinv = sp.Rational(1,10**8)
assert 13200*epsinv<sp.Rational(1,1000)
print('I2u inverse-first Hopf lift: value defect <=200r^(3/2), derivative defect <=5000r; blended cap defect <1/1000: PASS')

# Explicit constants for the closed inverse-first projector (programme 11.4).
# These checks verify the endpoint arithmetic of the written derivative proof;
# they do not replace the chain-rule estimates in that proof by sampling.
assert sp.Rational(55)<sp.Rational(9,4)**5
assert sp.Rational(4,9)-sp.Rational(9,32)-sp.Rational(1,10)==sp.Rational(91,1440)
assert 4+10<=14
assert 144/((1+160*16)**sp.Rational(3,2)*sp.sqrt(160+55*16))>sp.Rational(1,10**6)
assert 144/(95**sp.Rational(3,2)*sp.sqrt(sp.Rational(641,4)))>sp.Rational(1,100)
assert sp.sqrt(160)<13 and sp.Rational(641,4)/144<3
assert 4*128*100*10**6<10**12
assert sum(sp.Integer(4)**k/sp.factorial(k) for k in range(6))>sp.Rational(1000,27)
assert 4/(sp.Rational(3,4)**2)*55<392
assert (1600+sp.Rational(9,4)*10000)*sp.Rational(9,8)<10**5
assert (6300+sp.Rational(9,4)*10000)*sp.Rational(9,8)<10**5
assert sp.Rational(10,3)+sp.Rational(32,9)<7
assert sp.Rational(8,9)/(2*sp.sqrt(sp.Rational(3,2)))>sp.Rational(1,3)
assert 14*sp.Rational(55,4)<200
assert (3025+sp.Rational(9,4))*sp.Rational(55,4)<42000
assert 2*42000+7+40<100000
assert 3000+300000*4<2*10**6
assert 2000*(2*10**6+2)+1<10**16
assert 2*sp.Rational(44,7)*(1+20*Rinv)<14
assert sp.Rational(1,1000)*sp.Rational(1,100)<Rinv
assert 2*55*2000+2*55**2<10**6
assert 6*20**3<50000 and 6*3*20**2*1000<10**7
assert 10**7*10+50000*10**7*100<6*10**13
assert (1000+3*(6*10**13*20+500000*1000))/6<10**15
assert 7*(10**15+6*10**13)<10**16
assert 2*(14*10**12*5+3*10**8)<2*10**14
assert sp.sqrt(2)*(2*10**14+10**18)<10**20
delta_closed = sp.Rational(1,10**9)
L_closed = sp.Integer(10)**60*delta_closed**-6
assert L_closed==sp.Integer(10)**114
assert 10**51*delta_closed**-5<L_closed
N_closed = 4*L_closed+1
M_closed = 32*(N_closed+1)*L_closed
K_closed = 8*M_closed+1
assert 2*L_closed/N_closed<sp.Rational(1,2)
assert (9+8*N_closed)*L_closed<M_closed
assert 4*M_closed/K_closed<sp.Rational(1,2)
assert 100*sp.Rational(2,3)*sp.Rational(1,10**5)<sp.Rational(1,2)
assert 14*epsinv<sp.Rational(1,2*10**5)
from extra.s4t2_explicit.cusp import elementary_closed_parameters
parameters_closed = elementary_closed_parameters()
assert parameters_closed['projection_steps']==int(N_closed)
assert parameters_closed['phase_bound']==int(M_closed)
assert parameters_closed['phase_steps']==int(K_closed)
assert parameters_closed['delta']==delta_closed and parameters_closed['epsilon']==epsinv
print('I2v elementary height and global projector constants: L=10^114, N=4L+1, L_g=32(N+1)L, K=8L_g+1; implementation agrees: PASS')

print('ALL PASS')
