"""Fixed-point facts for the symbolic programme (docs/s4t2_symbolic_programme.md, section 4).

N1  phi(t) = Re(c2 al2 + c3 al3) with c = R(t)^-1 e4bar, exactly, for both cusps: a rational function of (Re t, Im t)
    with denominator 1 - |t|^4; its expansion is 2|t|^2 + O(|t|^4) (cusp 0) and 2|s|^2 + O(|s|^4) (cusp infinity),
    so e^{2g} = a0 e^{-2 phi} = a0 (1 - 4|t|^2 + ...).
N2  in the node variables Xh = a0 |x|^2, Y = |y|^2 (rad = Xh + Y, |t|^2 = Xh Y / a0): the chart series is
    mu = sum_j j T_j / sum_j T_j with T_j = e^{-(j^2+j) phi} Xh^{(j^2+j)/2} Y^{(j^2-j)/2}; the node height m0 = Xh - Y;
    mu - m0 = Y^2 - Xh^2 + (terms of total degree >= 3 in (Xh, Y, phi)).  The degree-3 part is printed.
N3  the fibre-derivative identity for the node-corrected height along a torus fibre (rho = log|y|, t fixed):
    d/d rho [mu + c(rad)(m0 - mu)] = (1 - c) mu' + c m0' - 2 c'(rad) m0 (m0 - mu),  with m0' = -2 rad, mu' = -2 Var.
N4  the Jacobian identity for the linear blend Lambda_mod = (1 - c(rad)) M0 vq + c(rad) Lambda(vq):
    D Lambda_mod = (1 - c) M0 + c D Lambda + (Lambda - M0 vq) (grad c)^T.
N5  numerical margins on the node ball rad <= 0.05 (both cusps; sanity, not proof): max |mu - m0| / rad^2,
    min Var / rad, and the minimum of -a1'/(2 rad) with the code's node correction (NODE_IN = 0.02, NODE_OUT = 0.05).
Run from the repo root: python3 extra/s4t2_explicit/sym/node_remainder.py
"""
import math
import sympy as sp

ok = True


def rep(name, cond):
    global ok
    ok &= bool(cond)
    print(f'{name:<84s} {"ok" if cond else "FAIL"}')


# ---------------------------------------------------------------- N1
tr, ti = sp.symbols('t_r t_i', real=True)
t = tr + sp.I * ti


def phi_of(al2, al3, e2, e3, e4):
    R = sp.Matrix([[sp.re(e2), sp.re(e3)], [sp.im(e2), sp.im(e3)]])
    c = R.inv() * sp.Matrix([sp.re(e4), sp.im(e4)])
    return sp.simplify(sp.re(sp.expand(c[0] * al2 + c[1] * al3))), sp.factor(R.det())


phi0, det0 = phi_of(t, -sp.I * t, 1 - t ** 2, sp.I * (1 + t ** 2), 2 * t)
phii, deti = phi_of(t, sp.I * t, t ** 2 - 1, sp.I * (1 + t ** 2), -2 * t)
print('cusp 0:   phi =', sp.factor(phi0), '  det R =', det0)
print('cusp inf: phi =', sp.factor(phii), '  det R =', deti)
r2 = tr ** 2 + ti ** 2
rep("N1a  cusp 0:   phi = 2|t|^2 (1 - Re t^2 ... )/(1 - |t|^4) has expansion 2|t|^2 + O(|t|^4)",
    sp.simplify(sp.series(phi0.subs({tr: sp.Symbol('e') * tr, ti: sp.Symbol('e') * ti}), sp.Symbol('e'), 0, 4).removeO()
                .subs(sp.Symbol('e'), 1) - 2 * r2) == 0)
rep("N1b  cusp inf: phi = 2|s|^2 + O(|s|^4)",
    sp.simplify(sp.series(phii.subs({tr: sp.Symbol('e') * tr, ti: sp.Symbol('e') * ti}), sp.Symbol('e'), 0, 4).removeO()
                .subs(sp.Symbol('e'), 1) - 2 * r2) == 0)
rep("N1c  det R = 1 - |t|^4 (cusp 0), -(1 - |t|^4)... printed above; both nonzero for |t| < 1",
    sp.simplify(det0 - (1 - r2 ** 2)) == 0 or sp.simplify(det0 + (1 - r2 ** 2)) == 0)
rep("N1d  phi = 2|t|^2/(1-|t|^2), exactly at BOTH cusps",
    all(sp.cancel(ph - 2 * r2 / (1 - r2)) == 0 for ph in (phi0, phii)))
rep("N1e  the two signed slice determinants",
    sp.expand(det0 - (1 - r2 ** 2)) == 0 and sp.expand(deti + (1 - r2 ** 2)) == 0)

# ---------------------------------------------------------------- N2
Xh, Y, ph = sp.symbols('Xh Y phi', real=True)
J = 3
T = {j: sp.exp(-(j * j + j) * ph) * Xh ** sp.Rational(j * j + j, 2) * Y ** sp.Rational(j * j - j, 2) for j in range(-J, J + 1)}
mu = sum(j * Tj for j, Tj in T.items()) / sum(T.values())
m0 = Xh - Y
eps = sp.Symbol('eps')
diff = (mu - m0).subs({Xh: eps * Xh, Y: eps * Y, ph: eps * ph})
ser = sp.series(diff, eps, 0, 4).removeO()
deg2 = sp.simplify(ser.coeff(eps, 2))
deg3 = sp.simplify(ser.coeff(eps, 3))
rep("N2a  mu - m0 has no degree-0, 1 terms in (Xh, Y, phi)", sp.simplify(ser.coeff(eps, 0)) == 0 and sp.simplify(ser.coeff(eps, 1)) == 0)
# phi itself is second order (phi = 2|t|^2/(1-|t|^2), |t|^2 = Xh Y/a0), so in (Xh, Y) alone mu - m0 = Y^2 - Xh^2 + O(3)
rep("N2b  degree-2 part of mu - m0 in (Xh, Y, phi) is Y^2 - Xh^2 - 2 Xh phi", sp.simplify(deg2 - (Y ** 2 - Xh ** 2 - 2 * Xh * ph)) == 0)
print('     degree-3 part in (Xh, Y, phi):', sp.factor(deg3))
rep("N2b' with phi = 2 Xh Y/(a0 - Xh Y) (|t|^2 = Xh Y/a0): mu - m0 = Y^2 - Xh^2 + O(3) in (Xh, Y)",
    sp.simplify(sp.series((mu - m0).subs(ph, 2 * Xh * Y / (a0 - Xh * Y) if False else 2 * eps ** 2 * Xh * Y / (sp.Symbol('a0', positive=True) - eps ** 2 * Xh * Y))
                          .subs({Xh: eps * Xh, Y: eps * Y}), eps, 0, 3).removeO().coeff(eps, 2) - (Y ** 2 - Xh ** 2)) == 0)
# the two-term truncation mu1 = (e^{-2phi} Xh - Y)/(1 + e^{-2phi} Xh + Y) differs from mu by terms with |j| >= 2:
mu1 = (sp.exp(-2 * ph) * Xh - Y) / (1 + sp.exp(-2 * ph) * Xh + Y)
d2 = (mu - mu1).subs({Xh: eps * Xh, Y: eps * Y, ph: eps * ph})
rep("N2c  mu - mu_two-term = O(eps^4)", all(sp.simplify(sp.series(d2, eps, 0, 4).removeO().coeff(eps, k)) == 0 for k in range(4)))

# N2d: the exact algebraic identity behind the I1 bound chain (notes section 4): with N = sum j T_j, D = sum T_j,
#   N - m0 D = Xh (e^{-2 phi} - 1) - (Xh - Y)(e^{-2 phi} Xh + Y) + R1 - (Xh - Y) R2,  R1 = sum_{|j|>=2} j T_j, R2 = sum_{|j|>=2} T_j
R1 = sum(j * Tj for j, Tj in T.items() if abs(j) >= 2)
R2 = sum(Tj for j, Tj in T.items() if abs(j) >= 2)
Nn, Dd = sum(j * Tj for j, Tj in T.items()), sum(T.values())
rep("N2d  N - m0 D = Xh(e^{-2phi} - 1) - (Xh - Y)(e^{-2phi} Xh + Y) + R1 - (Xh - Y) R2",
    sp.simplify(Nn - m0 * Dd - (Xh * (sp.exp(-2 * ph) - 1) - (Xh - Y) * (sp.exp(-2 * ph) * Xh + Y) + R1 - (Xh - Y) * R2)) == 0)
rep("N2e  j(j+1) >= 0 for all integers j, so e^{-(j^2+j) phi} <= 1 for phi >= 0; T_j has total degree j^2 in (Xh, Y)",
    all(j * (j + 1) >= 0 and sp.Poly(T[j].subs(ph, 0), Xh, Y).total_degree() == j * j for j in range(-J, J + 1)))

# ---------------------------------------------------------------- N3
rho = sp.Symbol('rho', real=True)
a0, tt = sp.symbols('a0 t2', positive=True)                       # t2 = |t|^2
Xr, Yr = a0 * tt * sp.exp(-2 * rho), sp.exp(2 * rho)              # Xh = a0 |x|^2 = a0 |t|^2 / |y|^2, Y = |y|^2
cf = sp.Function('c')
rad = Xr + Yr
mu_r = mu.subs({Xh: Xr, Y: Yr})
m0_r = Xr - Yr
a1 = mu_r + cf(rad) * (m0_r - mu_r)
lhs = sp.diff(a1, rho)
rhs = (1 - cf(rad)) * sp.diff(mu_r, rho) + cf(rad) * sp.diff(m0_r, rho) - 2 * sp.diff(cf(rad), rho) / sp.diff(rad, rho) * m0_r * (m0_r - mu_r)
rep("N3a  a1' = (1-c) mu' + c m0' - 2 c'(rad) m0 (m0 - mu)   [c' = dc/d rad; d rad/d rho = -2 m0]",
    sp.simplify(sp.diff(rad, rho) + 2 * m0_r) == 0 and sp.simplify(lhs - rhs) == 0)
rep("N3b  m0' = -2 rad", sp.simplify(sp.diff(m0_r, rho) + 2 * rad) == 0)
j_index = sp.Symbol('j', integer=True)
rep("N3c  scaling flow differentiates every chart term by -2j, including on the axes",
    sp.expand(-j_index * (j_index + 1) + j_index * (j_index - 1) + 2 * j_index) == 0
    and all(sp.simplify(-2 * Xh * sp.diff(Tj, Xh) + 2 * Y * sp.diff(Tj, Y) + 2 * j * Tj) == 0
            for j, Tj in T.items()))

# ---------------------------------------------------------------- N4 (matrix identity, 3 x 3, generic smooth maps)
q1, q2, q3 = sp.symbols('q1 q2 q3', real=True)
qv = sp.Matrix([q1, q2, q3])
Lf = sp.Matrix([sp.Function(f'L{i}')(q1, q2, q3) for i in range(3)])
M0 = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'M{i}{j}'))
cq = sp.Function('c')(q1, q2, q3)
Lmod = (1 - cq) * M0 * qv + cq * Lf
Dmod = Lmod.jacobian(qv)
expect = (1 - cq) * M0 + cq * Lf.jacobian(qv) + (Lf - M0 * qv) * sp.Matrix([[sp.diff(cq, s) for s in (q1, q2, q3)]])
rep("N4   D Lambda_mod = (1-c) M0 + c D Lambda + (Lambda - M0 q) (grad c)^T", sp.simplify(Dmod - expect) == sp.zeros(3, 3))

# ---------------------------------------------------------------- I1: exact constants, NOT sampling
# The infinite-series comparisons are proved in programme section 4.2.
# This block checks their algebra and all endpoint constants using rationals.
R = sp.Rational(1, 20)
d0, d1 = sp.Rational(201, 100), sp.Rational(401, 100)
e4_upper = sum(sp.Integer(4) ** k / sp.factorial(k) for k in range(17))
e4_upper += sp.Integer(4) ** 17 / sp.factorial(17) / (1 - sp.Rational(4, 18))
rep("I1a  exp(4) < 55 (positive Taylor series, geometric bound after k=16)", e4_upper < 55)
k, ell, q = sp.symbols('k ell q', integer=True, nonnegative=True)
rep("I1b  (ell+2)^2 >= 4+5ell: remainder ell(ell-1) >= 0 for integer ell >= 0",
    sp.expand((ell + 2) ** 2 - 4 - 5 * ell - ell * (ell - 1)) == 0)
tail0 = 2 / (1 - R ** 5)
tail1 = 2 * (2 - R ** 5) / (1 - R ** 5) ** 2
rep("I1c  sum (ell+2) q^ell = (2-q)/(1-q)^2 for 0 <= q < 1",
    sp.cancel(q * sp.diff(1 / (1 - q), q) + 2 / (1 - q) - (2 - q) / (1 - q) ** 2) == 0)
rep("I1d  R2 <= 2.01 rad^4 and |R1| <= 4.01 rad^4", tail0 < d0 and tail1 < d1)
rep("I1e  4(Xh+Y)^3/27 - Xh^2 Y = (Xh-2Y)^2(4Xh+Y)/27",
    sp.expand(4 * (Xh + Y) ** 3 - 27 * Xh ** 2 * Y - (Xh - 2 * Y) ** 2 * (4 * Xh + Y)) == 0)
B = sp.Rational(16, 27) * 55 / (1 - 55 * R ** 2 / 4)
C = sp.Rational(27, 10)
C_estimate = 1 + B * R + d1 * R ** 2 + d0 * R ** 3
C2_estimate = 1 + d0 * R ** 3 + B * R + (1 + C * R) ** 2
rep("I1f  2 phi Xh <= B rad^3; numerator of variance bound stays positive", B * R ** 2 < 1)
rep("I1g  |mu-m0| <= (27/10) rad^2 at both cusps", C_estimate < C)
rep("I1h  Var >= rad(1-4rad) at both cusps", C2_estimate < 4)
margin = 1 - 4 * R - 2 * C * R ** 2 / sp.Rational(3, 100)
rep("I1i  -a1'/(2rad) >= 7/20 > 0 for 0 < rad <= 1/20", margin == sp.Rational(7, 20))
print(f'PROVEN I1: C={C}, C2=4, margin={margin}; measured margins are printed separately below.')

# ---------------------------------------------------------------- N5 numerical margins (sanity)
import numpy as np
for name, a0v, phifun in (('cusp 0', 1.0, sp.lambdify((tr, ti), phi0, 'numpy')), ('cusp inf', math.exp(-4), sp.lambdify((tr, ti), phii, 'numpy'))):
    Jn = 6
    g = np.random.default_rng(0)
    N = 200000
    radv = g.uniform(0.0005, 0.05, N)
    frac = g.uniform(0, 1, N)
    Xhv, Yv = radv * frac, radv * (1 - frac)
    th = g.uniform(-np.pi, np.pi, N)
    t2 = Xhv * Yv / a0v
    trv, tiv = np.sqrt(t2) * np.cos(th), np.sqrt(t2) * np.sin(th)
    phv = phifun(trv, tiv)
    num = np.zeros(N); den = np.zeros(N); sq = np.zeros(N)
    for j in range(-Jn, Jn + 1):
        Tj = np.exp(-(j * j + j) * phv) * Xhv ** ((j * j + j) / 2) * Yv ** ((j * j - j) / 2)
        num += j * Tj; den += Tj; sq += j * j * Tj
    muv = num / den
    var = sq / den - muv ** 2
    m0v = Xhv - Yv
    C = np.max(np.abs(muv - m0v) / radv ** 2)
    C2 = np.max((1 - var / radv) / radv)
    # code's node correction: c = 1 - s((rad - 0.02)/0.03), c' = -s'/0.03
    uu = np.clip((radv - 0.02) / 0.03, 0, 1)
    with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
        ea, eb = np.exp(-1 / np.maximum(uu, 1e-12)), np.exp(-1 / np.maximum(1 - uu, 1e-12))
        ea[uu <= 0] = 0; eb[uu >= 1] = 0
        sv = ea / (ea + eb)
        spv = sv * (1 - sv) * (1 / np.maximum(uu, 1e-12) ** 2 + 1 / np.maximum(1 - uu, 1e-12) ** 2)
        spv[(uu <= 0) | (uu >= 1)] = 0
    cv, cpv = 1 - sv, -spv / 0.03
    a1p = (1 - cv) * (-2 * var) + cv * (-2 * radv) - 2 * cpv * m0v * (m0v - muv)
    print(f'{name}: max |mu - m0|/rad^2 = {C:.4f};  max (1 - Var/rad)/rad = {C2:.3f};  min of -a1\'/(2 rad) = {np.min(-a1p / (2 * radv)):.4f}')
    rep(f"N5   {name}: a1' < 0 on the node ball (sampled, 2e5 points)", np.all(a1p < 0))

# ---------------------------------------------------------------- N6: finite-height replacement, independent of the Gaussian estimate
rf = sp.Rational(1,20)
sf = 55*rf**2/4
Cf = (sp.Rational(16,27)*55/(1-sf)+sp.Rational(1,4))/(1-rf**2/4)
lowerf = 1-4*sf/(1-sf)
marginf = lowerf-2*35*rf**3/sp.Rational(3,100)
rep("N6a finite height: |H-m0|<=35 rad^3", Cf < 35)
rep("N6b finite height: V/rad>=265/309; corrected slope >14/25",
    lowerf == sp.Rational(265,309) and marginf > sp.Rational(14,25))
rep("N6c active reduced chart |H|<.055, so the node integer is unique",
    rf+35*rf**3 < sp.Rational(11,200))
Af,Yf,qf,uf,vf,ef = sp.symbols('Af Yf qf uf vf ef',positive=True)
rep("N6d exact node remainder numerator",
    sp.cancel((ef*uf-vf)/(1-qf)-(uf-vf)
              -(-uf*(1-ef)+qf*(uf-vf))/(1-qf)) == 0)
print('PROVEN finite height: cubic remainder constant 35; corrected slope >14/25.')
for label,a0f in (('0',1.),('inf',math.exp(-4))):
    rngf = np.random.default_rng(918)
    nr = 200000
    rr = rngf.uniform(.0005,.05,nr)
    uu = rr*rngf.uniform(0,1,nr)
    vv = rr-uu
    ss = uu*vv/a0f
    ee = np.exp(-4*ss/(1-ss))
    qq = ee*uu*vv
    hf = (ee*uu-vv)/(1-qq)
    vfib = (ee*uu+vv)/(1-qq)
    mm = uu-vv
    ccarg = np.clip((rr-.02)/.03,1e-10,1-1e-10)
    with np.errstate(over='ignore',under='ignore'):
        aa,bb = np.exp(-1/ccarg),np.exp(-1/(1-ccarg))
    st = aa/(aa+bb)
    dst = st*(1-st)*(1/ccarg**2+1/(1-ccarg)**2)
    slope = ((st*vfib+(1-st)*rr)-dst*mm*(mm-hf)/.03)/rr
    rem = np.abs(hf-mm)/rr**3
    print(f'finite cusp {label}, {nr} samples: remainder/rad^3 max {rem.max():.6f}; slope min {slope.min():.6f}')
    rep(f"N6e finite cusp {label}: measured remainder and slope within proven bounds",
        rem.max()<35 and slope.min()>14/25)
    # Independent atan-periodization evaluator for the elliptic height. In this
    # node ball q<1/1600, so the eight-term omitted tail is <2e-27.
    aa,bb = ee*uu,vv
    Lell,Vell = np.zeros(nr),np.zeros(nr)
    for _ in range(8):
        Lell += np.arctan(aa)-np.arctan(bb)
        Vell += aa/(1+aa*aa)+bb/(1+bb*bb)
        aa,bb = qq*aa,qq*bb
    Hell = np.arctan((math.pi/2)*np.tan(2*Lell))/math.pi
    dell = 1/(np.cos(2*Lell)**2+(math.pi/2)**2*np.sin(2*Lell)**2)
    vell = dell*Vell
    eslope = ((st*vell+(1-st)*rr)-dst*mm*(mm-Hell)/.03)/rr
    erem = np.abs(Hell-mm)/rr**3
    print(f'elliptic cusp {label}, {nr} samples: remainder/rad^3 max {erem.max():.6f}; slope min {eslope.min():.6f}')
    rep(f"N7 elliptic cusp {label}: remainder<38 and corrected slope>13/25",
        erem.max()<38 and eslope.min()>13/25)

# ---------------------------------------------------------------- N8: the actual two-way construction on the exact inner pole caps.
# On |z_L|^2 <= 10^-8 the pole conjugation H_j is the identity.  The four
# inputs of scalar_finite_height are then all < 1/16, so its four seed branches
# are their identity pieces.  This check includes the resulting H_f in tau_f^-1;
# the earlier height-free torus gauge of programme 9 is a different map.
round_r8 = sp.symbols('round_r8', nonnegative=True)
# On the cap 0<=r<1, both square roots are nonnegative; their equal squares
# therefore identify the positive branches used in (11.17).
rep('N8a0 round-cap gain: |X|sqrt(2-X^2)=sqrt(1-r^2), r=|Z|^2<1',
    sp.expand((1-round_r8)*(1+round_r8)-(1-round_r8**2)) == 0)
eps8 = sp.Rational(1,10**8)
rep('N8a00 inner cap: A,Y <= epsilon and qA,qY <= epsilon^3 < 1/16',
    eps8 < sp.Rational(1,16) and eps8**3 < sp.Rational(1,16))
Ah8, Yh8 = sp.symbols('Ah8 Yh8', nonnegative=True)
qh8 = Ah8 * Yh8
four_seed8 = Ah8 + qh8*Ah8 - Yh8 - qh8*Yh8 + qh8**2*(Ah8-Yh8)/(1-qh8)
rep('N8a inner-cap four-seed torus height H_f=(A-Y)/(1-AY)',
    sp.cancel(four_seed8-(Ah8-Yh8)/(1-qh8)) == 0)

tr8, ti8, h8, ph28, ph38 = sp.symbols('tr8 ti8 h8 ph28 ph38', real=True)
t8 = tr8 + sp.I*ti8
v8 = tr8**2 + ti8**2
for sigma8 in (1, -1):
    W8 = ph28 - sp.I*sigma8*ph38
    C28 = ph28 + h8*2*tr8/(1-v8)
    C38 = ph38 + h8*2*sigma8*ti8/(1-v8)
    A_direct8 = t8*(C28-sp.I*sigma8*C38)
    w_direct8 = C28*sigma8*(1-t8**2)+C38*sp.I*(1+t8**2)
    A_compact8 = t8*W8+h8*2*v8/(1-v8)
    w_compact8 = sigma8*(sp.conjugate(W8)-W8*t8**2+2*h8*t8)
    rep(f'N8b sigma={sigma8}: compact cap exponent and third holomorphic coordinate',
        sp.cancel(A_direct8-A_compact8) == 0
        and sp.cancel(w_direct8-w_compact8) == 0)

# The Jacobian and all six Hessians below are differentiated from the specific
# six-real-coordinate map h=(Re x,Im x,Re y,Im y,Re w,Im w).  Complex real/imag
# parts are expanded by hand before differentiation to keep SymPy's expression
# DAG small.  The normal-frame reduction after this loop never substitutes
# sampled values or assumes that an arbitrary matrix field is a Jacobian.
z8 = sp.symbols('z8_0:4', real=True)
phi28, phi38 = sp.symbols('phi28 phi38', real=True)
u8 = (*z8,phi28,phi38)
J08 = sp.diag(*([sp.Matrix([[0,-1],[1,0]])]*3))
dim8 = 6
T8 = [[[sp.Symbol(f'T8_{out}_{min(i,k)}{max(i,k)}') for k in range(dim8)]
       for i in range(dim8)] for out in range(dim8)]
dJ8 = []
for k in range(dim8):
    Hk8 = sp.Matrix(dim8,dim8,lambda out,i: T8[out][i][k])
    # J=M^-1 J0 M, so at a point in a frame with M=I, dJ=[J0,Hk].
    dJ8.append(J08*Hk8-Hk8*J08)
square8 = J08*J08+sp.eye(dim8)
N8 = [sp.expand(sum(J08[p,i]*dJ8[p][out,j]
                    -J08[p,j]*dJ8[p][out,i]
                    -J08[out,p]*(dJ8[i][p,j]-dJ8[j][p,i])
                    for p in range(dim8)))
      for out in range(dim8) for i in range(dim8) for j in range(dim8)]
rep('N8c normal-frame algebra: all 36 J^2 and 216 N entries vanish',
    len(square8) == 36 and all(x == 0 for x in square8)
    and len(N8) == 216 and all(x == 0 for x in N8))

import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))

for sigma8, scale8, nu8 in ((1,2*sp.pi,sp.Integer(1)),
                            (-1,sp.pi,sp.exp(-2))):
    a8,b8,c8,d8 = z8
    r8 = sum(v*v for v in z8)
    s8 = sp.sqrt(1-r8*r8)
    gain8 = 1/sp.sqrt(1+s8)
    den8 = sp.sqrt(2*scale8)
    if sigma8 == 1:
        p8,q8 = gain8*(a8+d8)/den8,gain8*(b8-c8)/den8
        m8,n8 = gain8*(a8-d8)/den8,gain8*(b8+c8)/den8
    else:
        p8,q8 = gain8*(a8+c8)/den8,gain8*(b8+d8)/den8
        m8,n8 = gain8*(a8-c8)/den8,gain8*(b8-d8)/den8
    # p+iq and m+in are U_j^*Z/sqrt(s_j(1+s)); x_s=(p-iq)/nu.
    tr_actual8 = (p8*m8+q8*n8)/nu8
    ti_actual8 = (p8*n8-q8*m8)/nu8
    v_actual8 = tr_actual8**2+ti_actual8**2
    e_actual8 = sp.exp(-4*v_actual8/(1-v_actual8))
    A_actual8 = e_actual8*(p8*p8+q8*q8)
    Y_actual8 = m8*m8+n8*n8
    q_actual8 = A_actual8*Y_actual8
    H_actual8 = (A_actual8-Y_actual8)/(1-q_actual8)
    f_actual8 = 2*v_actual8/(1-v_actual8)
    ar8 = tr_actual8*phi28+sigma8*ti_actual8*phi38+H_actual8*f_actual8
    ai8 = ti_actual8*phi28-sigma8*tr_actual8*phi38
    co8,si8 = sp.cos(ai8),sp.sin(ai8)
    ea8 = sp.exp(ar8)
    B8 = tr_actual8**2-ti_actual8**2
    C8 = 2*tr_actual8*ti_actual8
    h_actual8 = sp.Matrix([
        (p8*co8-q8*si8)/(nu8*ea8),
        (-p8*si8-q8*co8)/(nu8*ea8),
        ea8*(m8*co8-n8*si8),
        ea8*(m8*si8+n8*co8),
        sigma8*phi28*(1-B8)-phi38*C8+2*sigma8*H_actual8*tr_actual8,
        phi38*(1+B8)-sigma8*phi28*C8+2*sigma8*H_actual8*ti_actual8])
    M_actual8 = h_actual8.jacobian(u8)
    Hess_actual8 = [sp.hessian(h_actual8[out],u8) for out in range(dim8)]
    rep(f'N8d sigma={sigma8}: actual cap 6x6 Jacobian and all 90 Hessian symmetries',
        M_actual8.shape == (6,6) and len(Hess_actual8) == 6
        and not any(x.has(sp.Derivative) for x in M_actual8)
        and not any(x.has(sp.Derivative) for Hx in Hess_actual8 for x in Hx)
        and all(Hx[i,j] == Hx[j,i] for Hx in Hess_actual8
                for i in range(dim8) for j in range(i+1,dim8)))
    # These actual Hessians, after the invertible linear normalization M=I,
    # belong to the symmetric T8 family checked in N8c.  The determinant is
    # certified by the radial, linear-spinor, torus-shear and slice factors;
    # the direct pole value is an independent exact check of their signs.
    M_pole8 = M_actual8.subs(dict.fromkeys(z8,0))
    expected_pole8 = -sigma8/(4*scale8**2*nu8**2)
    rep(f'N8e sigma={sigma8}: actual pole Jacobian determinant',
        sp.simplify(M_pole8.det()-expected_pole8) == 0)
    # Numerical cross-check only: evaluate the independently differentiated
    # SymPy formula against the production cap and torus restoration at the
    # exact pole and two interior points.  It does not certify the continuum.
    import mpmath as mp
    from extra.s4t2_explicit import cusp as cusp8
    h_eval8 = sp.lambdify(u8,h_actual8,'mpmath',cse=True)
    samples8 = ((mp.mpf(0),)*4,
                tuple(map(mp.mpf,('1e-5','2e-5','-1e-5','3e-5'))),
                tuple(map(mp.mpf,('-2e-5','1e-5','3e-5','-4e-5'))))
    worst8 = mp.mpf(0)
    with mp.workdps(70):
        angles8 = (mp.mpf('.123'),mp.mpf('-.279'))
        for vals8 in samples8:
            X8 = -sigma8*mp.sqrt(1-sum(v*v for v in vals8))
            kind8,production8 = cusp8.elementary_global_coordinates(
                (mp.mpc(vals8[0],vals8[1]),mp.mpc(vals8[2],vals8[3])),
                X8,angles8,dps=65)
            direct8 = h_eval8(*vals8,*angles8)
            target8 = mp.matrix([part for z in production8
                                 for part in (z.real,z.imag)])
            worst8 = max(worst8,mp.norm(direct8-target8))
            assert kind8 == ('c0' if sigma8 == 1 else 'ci')
    rep(f'N8e2 sigma={sigma8}: 3 production inner-cap samples (numerical), '
        f'max error {float(worst8):.1e}',worst8 < mp.mpf('1e-55'))

rad8 = sp.symbols('rad8',positive=True)
gr8 = 1/sp.sqrt(1+sp.sqrt(1-rad8**2))
rep('N8f radial 4x4 determinant factor g^4/sqrt(1-r^2)',
    sp.simplify(gr8+2*rad8*sp.diff(gr8,rad8)
                -gr8/sp.sqrt(1-rad8**2)) == 0)
bad8 = [sp.zeros(dim8) for _ in range(dim8)]
bad8[2][0,0] = 1
bad_dJ8 = [J08*B-B*J08 for B in bad8]
bad_N8 = sp.expand(sum(J08[p,0]*bad_dJ8[p][1,3]
                       -J08[p,3]*bad_dJ8[p][1,0]
                       -J08[1,p]*(bad_dJ8[0][p,3]-bad_dJ8[3][p,0])
                       for p in range(dim8)))
rep('N8g non-Jacobian-frame control: N^1_(0,3)=+1',bad_N8 == 1)

# ---------------------------------------------------------------- N9: exact cap/slab overlap, conditional on §11.5's bounds.
# The type-(c) derivative bound ||Dg_j||<L_g, the cap projector
# identification in (11.10), and d<=7r_L are separate written certificates.
# The real agreement of the north/south phase branches uses the connectedness
# argument in §10.6.  N9 checks their exact algebraic consequence on
# alpha/3 < 1-|X| < alpha/2.  In
# particular, the K increments are legitimate smooth atan branches by that
# bound; the following telescoping is not a sampled phase comparison.
import ast
import inspect


def has_statement9(function, statement):
    tree = ast.parse(inspect.getsource(function))
    expected = ast.dump(ast.parse(statement).body[0], include_attributes=False)
    return any(ast.dump(node, include_attributes=False) == expected
               for node in ast.walk(tree))


def has_if_test9(function, test):
    tree = ast.parse(inspect.getsource(function))
    expected = ast.dump(ast.parse(test, mode='eval').body,
                        include_attributes=False)
    return any(isinstance(node, ast.If)
               and ast.dump(node.test, include_attributes=False) == expected
               for node in ast.walk(tree))


rep('N9a source cutoff, direct mismatch, pole selector and cusp action',
    has_statement9(cusp8._elementary_corrected_section,
                   'direct = _elementary_unit(P*mp.matrix(_elementary_pole_inclusion(xx,inf)))')
    and has_statement9(cusp8._elementary_corrected_section,
                       'return _elementary_inner(direct,section)')
    and has_statement9(cusp8._elementary_corrected_section,
                       'weight = 1-scalar_flat_step((distance-alpha/2)/(alpha/2))')
    and has_statement9(cusp8._elementary_corrected_section,
                       'result = tuple(a*mp.exp(-mp.j*weight*theta) for a in result)')
    and has_if_test9(cusp8.elementary_global_slice,
                     'distance<alpha/2 and not force_slab')
    and has_statement9(cusp8.elementary_global_slice,
                       'xi = tuple(a/mp.sqrt(1-X*X) for a in Z)')
    and has_statement9(cusp8.elementary_global_slice,
                       'den = mp.sqrt(scale*(1+abs(X)*mp.sqrt(2-X*X)))')
    and has_statement9(cusp8.elementary_global_slice,
                       'zz = elementary_pole_patch(tuple(a/den for a in pair),inf,dps=dps)')
    and has_statement9(cusp8.elementary_global_slice,
                       'phase = _elementary_inner(ref,result)')
    and has_statement9(cusp8.elementary_global_slice,
                       'return kind,(mp.conj(phase)*x,phase*y)')
    and has_statement9(cusp8._elementary_pole_radial,
                       'if radius>=c:\n    return z')
    and has_statement9(cusp8.scalar_flat_step,
                       'if u<=0:\n    return mp.mpf(0)')
    and has_statement9(cusp8.elementary_pole_patch,
                       'raw = elliptic_pole_raw_inverse(compressed,infinity,dps=dps,elementary=True)')
    and has_statement9(cusp8.elliptic_pole_raw_inverse,
                       'factor = mp.sqrt(source_radius)/length')
    and has_statement9(cusp8.elementary_restore_torus,
                       'H = scalar_finite_height(x,y,mp.exp(-2) if inf else mp.mpf(1))')
    and has_statement9(cusp8.elementary_restore_torus,
                       'return kind,(x*mp.exp(-A),y*mp.exp(A),C[0]*e2+C[1]*e3)')
    and has_statement9(cusp8.elementary_global_coordinates,
                       'return elementary_restore_torus(kind,data,phi)'))

alpha9, epsilon9 = sp.Rational(1,10**5), sp.Rational(1,10**8)
d9 = sp.Symbol('d9', positive=True)
cutoff_argument9 = sp.cancel((d9-alpha9/2)/(alpha9/2))
rep('N9b given d<=7r_L: overlap cutoff=1, r_L>2epsilon',
    cutoff_argument9.subs(d9,alpha9/3) == -sp.Rational(1,3)
    and cutoff_argument9.subs(d9,alpha9/2) == 0
    and alpha9/21 > 2*epsilon9
    and alpha9/21 > sp.Rational(3,2)*epsilon9
    # By (11.11), |g_j(z)-z|<=200r_L^(3/2).  The lower bound below makes
    # the *output* radial inverse the identity as well as the input radial
    # compression.  We use r_L<2alpha/3, weaker than the overlap's bound.
    and 200**2*2*alpha9/3 < sp.Rational(13,25)**2
    and (alpha9/21)*(1-sp.Rational(13,25))**2
        > sp.Rational(3,2)*epsilon9)

# At the certified image radius < 2e-6, the four finite-height seeds are
# all in their identity pieces and |H_f|<.1.  Therefore the only nonzero
# embedding weight is the first cusp pair, e(w)=w/sqrt(|w|^2).
source_bound9 = sp.Rational(2,10**6)
height_bound9 = source_bound9/(1-source_bound9**2)
rep('N9c certified cusp image: first embedding weight 1, second 0',
    source_bound9 < sp.Rational(1,16)
    and height_bound9 < sp.Rational(1,10)
    and has_statement9(cusp8.elliptic_bundle_embedding,
                       'block = [hw(H)*q for q in f]')
    and has_statement9(cusp8.elliptic_bundle_embedding,
                       "w = hw(H+mp.mpf('.5') if H<0 else H-mp.mpf('.5'))"))

# The north/south K-increment products telescope for arbitrary unit phases;
# induction inserts any further increment, including the enormous fixed K.
gA9,gN9,gE9,gS9,gT9 = sp.symbols('gA9 gN9 gE9 gS9 gT9', nonzero=True)
north_phase9 = gA9*(gN9/gA9)*(gT9/gN9)
south_phase9 = gA9*(gN9/gA9)*(gE9/gN9)*(gS9/gE9)*(gT9/gS9)
rep('N9d finite real-log exponent: north and south phase products equal g_j',
    sp.cancel(north_phase9-gT9) == 0
    and sp.cancel(south_phase9-gT9) == 0)

# E and d_j are unit vectors in the same final rank-one projector line.
# Hence E=e^{i theta}d_j and d_j^*E=e^{i theta}.  With the weight exactly 1,
# the corrected finite section is d_j on the whole open overlap.
theta9 = sp.symbols('theta9', real=True)
q19,q29,q39,q49 = sp.symbols('q19 q29 q39 q49', real=True)
qvec9 = sp.Matrix([q19+sp.I*q29,q39+sp.I*q49])
norm29 = sum(sp.expand(sp.conjugate(q)*q) for q in qvec9)
rho9 = sp.Symbol('rho9', positive=True)  # rho9^2=norm29; nonzero by the cap bound.
dvec9 = qvec9/rho9
Evec9 = sp.exp(sp.I*theta9)*dvec9
mismatch9 = (dvec9.conjugate().T*Evec9)[0]
corrected9 = sp.exp(-sp.I*theta9)*Evec9
rep('N9e finite correction E_c=d_j for either pole projector line',
    sp.expand(sp.cancel((mismatch9-sp.exp(sp.I*theta9))*rho9**2
                        /sp.exp(sp.I*theta9)).subs(rho9**2,norm29)) == 0
    and all(sp.simplify(corrected9[k]-dvec9[k]) == 0 for k in range(2)))
rep('N9e2 unit input phase commutes with every normalized projection',
    sp.simplify(sum(sp.conjugate(z)*z for z in Evec9)-norm29/rho9**2) == 0
    and has_statement9(cusp8.elementary_projection_section,
                       'vector = tuple(a*phase/abs(phase) for a in ref)'))

# z_L=c I_j xi with c>0, since its two denominators are positive on the
# overlap.  Thus normalization commutes with this positive rescaling.
scale9, source_radius9 = sp.symbols('scale9 source_radius9', positive=True)
raw9 = sp.sqrt(source_radius9)*(scale9*qvec9)/(scale9*rho9)
rep('N9f raw pole spinor equals sqrt(source radius) times d_j',
    all(sp.simplify(raw9[k]-sp.sqrt(source_radius9)*dvec9[k]) == 0
        for k in range(2)))

# A local source reference is a phase rotation of the same raw spinor.
# Its normalized embedding is e_ref, and lambda=e_ref^*d_j removes that
# phase.  The production cusp action multiplies (nu*conj(x),y) by lambda.
u9 = sp.symbols('u9_0:4', real=True)
phi9 = sp.symbols('phi9_0:2', real=True)
beta9 = sp.Function('beta9', real=True)(*u9)
xr9 = sp.Function('xr9', real=True)(*u9)
xi9 = sp.Function('xi9', real=True)(*u9)
yr9 = sp.Function('yr9', real=True)(*u9)
yi9 = sp.Function('yi9', real=True)(*u9)
nu9 = sp.symbols('nu9', positive=True)
xraw9,yraw9 = xr9+sp.I*xi9,yr9+sp.I*yi9
lambda9 = sp.exp(-sp.I*beta9)
eref9 = sp.exp(sp.I*beta9)*dvec9
lambda_from_inner9 = (eref9.conjugate().T*dvec9)[0]
xref9 = sp.conjugate(sp.exp(sp.I*beta9)*nu9*sp.conjugate(xraw9))/nu9
yref9 = sp.exp(sp.I*beta9)*yraw9
xslab9 = sp.conjugate(lambda9)*xref9
yslab9 = lambda9*yref9
rep('N9g source circle action gives the identical raw cusp W point',
    sp.expand(sp.cancel((lambda_from_inner9-lambda9)*rho9**2
                        /lambda9).subs(rho9**2,norm29)) == 0
    and sp.simplify(xslab9-xraw9) == 0
    and sp.simplify(yslab9-yraw9) == 0)

def restored9(x,y,sigma):
    t = x*y
    rr = t*sp.conjugate(t)
    AA = nu9**2*sp.exp(-4*rr/(1-rr))*x*sp.conjugate(x)
    YY = y*sp.conjugate(y)
    HH = (AA-YY)/(1-AA*YY)  # Exact four-seed H_f here, by N8a/N9c.
    W = phi9[0]-sp.I*sigma*phi9[1]
    AAcomplex = t*W+HH*2*rr/(1-rr)
    return (x*sp.exp(-AAcomplex),y*sp.exp(AAcomplex),
            sigma*(sp.conjugate(W)-W*t*t+2*HH*t))

# Simplify the equal W inputs before composition.  This is the exact open
# equality argument: all first and second jets of the full torus-restored
# maps then agree, without expanding the astronomical projection product.
# Per cusp, 18 complex first entries and 108 complex second entries are
# exactly 36 and 216 real entries, respectively.
hraw9 = [restored9(xraw9,yraw9,sigma) for sigma in (1,-1)]
hslab9 = [restored9(sp.simplify(xslab9),sp.simplify(yslab9),sigma)
          for sigma in (1,-1)]
jets9 = []
for h0,h1 in zip(hraw9,hslab9):
    for a,b in zip(h0,h1):
        delta = b-a
        jets9.extend([delta]+[sp.diff(delta,v) for v in (*u9,*phi9)]
                     +[sp.diff(delta,v,w) for v in (*u9,*phi9)
                       for w in (*u9,*phi9)])
rep('N9h both cusps: exact torus restoration and all 36/216 first/second jets',
    len(jets9) == 2*3*(1+6+36) and all(v == 0 for v in jets9))

# N9i closes the source-reference premise used in N9g.  The separate cap
# bounds in §11.6 give |z_L|^2<2e-6 and raw source norm^2<2e-6; since
# s_j<=2pi<7, hence
# 0<h=s_j|z_L|^2<1/10, |height|<1/10, and Qz_L != 0.
# The following identities are exact on that open region; they do not
# evaluate the astronomical finite projection loop or prove its bounds.
zr9i,zi9i,wr9i,wi9i = sp.symbols('zr9i zi9i wr9i wi9i', real=True)
z9i = sp.Matrix([zr9i+sp.I*zi9i,wr9i+sp.I*wi9i])
def hopf9i(pair):
    pc = sp.expand(2*pair[0]*sp.conjugate(pair[1]))
    return sp.Matrix([sp.re(pc).expand(),sp.im(pc).expand(),
                      sp.expand(pair[0]*sp.conjugate(pair[0])
                                -pair[1]*sp.conjugate(pair[1]))])

p9i = hopf9i(z9i)
r9i = sum(sp.expand(v*sp.conjugate(v)) for v in z9i)
U09i = sp.Matrix([[1,1],[sp.I,-sp.I]])/sp.sqrt(2)
Ui9i = sp.Matrix([[1,1],[1,-1]])/sp.sqrt(2)
O09i = sp.Matrix([-p9i[1],-p9i[2],p9i[0]])
Oi9i = sp.Matrix([p9i[2],-p9i[1],p9i[0]])
hopf_source9i = (
    sp.simplify(hopf9i(U09i*z9i)-O09i)==sp.zeros(3,1)
    and sp.simplify(hopf9i(Ui9i*z9i)-Oi9i)==sp.zeros(3,1)
    and sp.expand(sum(v*v for v in p9i)-r9i**2)==0)
h9i,scale9i = sp.symbols('h9i scale9i',positive=True)
r_input9i = sp.Symbol('r_input9i',positive=True)
rL_bound9i = sp.Rational(2,10**6)
root9i = (1+h9i)/sp.sqrt(1+h9i*h9i)
round_source9i = True
for sign9i in (-1,1):
    X9i = sign9i*(1-h9i)/sp.sqrt(1+h9i*h9i)
    aa9i = sign9i*(1-h9i*h9i)/(1+h9i*h9i)
    # The positive root is root9i.  For south the source stereo coordinate
    # is L=s_j O_j p; for north it is L/h^2.
    expected_gain9i = scale9i if sign9i==-1 else scale9i/h9i**2
    round_source9i &= (
        sp.simplify(2-X9i**2-root9i**2)==0
        and sp.simplify(X9i*root9i-aa9i)==0
        and sp.simplify(1-X9i**2-2*h9i/(1+h9i*h9i))==0
        and sp.simplify((2*scale9i/(1+h9i*h9i))/(1-aa9i)
                        -expected_gain9i)==0)
round_source9i &= (
    sp.simplify(sp.sqrt(2*scale9i/(1+h9i*h9i))
                /sp.sqrt(2*scale9i*r_input9i/(1+h9i*h9i))
                -1/sp.sqrt(r_input9i))==0
    and sp.simplify(sp.sqrt(2*scale9i/(1+h9i*h9i))
                    /sp.sqrt(scale9i*(1+(1-h9i*h9i)/(1+h9i*h9i)))
                    -1)==0)
round_ast9i = (
    has_statement9(cusp8._elementary_pole_inclusion,
                   'pair = ((a+b)/mp.sqrt(2),(a-b)/mp.sqrt(2)) if infinity else ((a-mp.j*b)/mp.sqrt(2),(a+mp.j*b)/mp.sqrt(2))')
    and has_statement9(cusp8.elementary_global_slice,
                       'den = mp.sqrt(scale*(1+abs(X)*mp.sqrt(2-X*X)))')
    and has_statement9(cusp8.elementary_global_slice,
                       'zz = elementary_pole_patch(tuple(a/den for a in pair),inf,dps=dps)')
    and has_statement9(cusp8.elliptic_pole_raw_inverse,
                       'v = tuple(scale*w for w in ((p[2],-p[1],p[0]) if infinity else (-p[1],-p[2],p[0])))'))
rep('N9i1 both Hopf rotations and round-to-pole stereographic gains',
    hopf_source9i and round_source9i and round_ast9i
    and 7*rL_bound9i<sp.Rational(1,10))

# In the north, substituting w=L/|L|^2 into the ordinary inverse
# stereographic formula produces exactly the source's stable rational u.
# South uses w=L directly.  Write L=(a,b,c), |L|^2=a^2+b^2+c^2.
la9i,lb9i,lc9i = sp.symbols('la9i lb9i lc9i',real=True)
lh29i = la9i**2+lb9i**2+lc9i**2
def stereo9i(y0,y1,y2):
    s2 = y0*y0+y1*y1+y2*y2
    return ((1-s2+2*sp.I*y0)/(1+s2),
            2*(y1+sp.I*y2)/(1+s2))

unorth9i = stereo9i(la9i/lh29i,1+lb9i/lh29i,lc9i/lh29i)
north_den9i = 1+2*lb9i+2*lh29i
unorth_source9i = ((-1-2*lb9i+2*sp.I*la9i)/north_den9i,
                   2*(lb9i+lh29i+sp.I*lc9i)/north_den9i)
usouth9i = stereo9i(la9i,1+lb9i,lc9i)
south_s29i = la9i**2+(lb9i+1)**2+lc9i**2
usouth_source9i = ((1-south_s29i+2*sp.I*la9i)/(1+south_s29i),
                    2*(lb9i+1+sp.I*lc9i)/(1+south_s29i))
raw_base_source9i = (
    has_statement9(cusp8.elliptic_raw_base_reference,
                   'aa = xx*mp.sqrt(2-xx*xx)')
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       'w = (vv[0]/(1-aa),-vv[2]/(1-aa),-vv[1]/(1-aa))')
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       'y = (w[0],w[1]+1,w[2])')
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       'u1 = mp.mpc(1-s2,2*y[0])/(1+s2)')
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       'u2 = mp.mpc(2*y[1],2*y[2])/(1+s2)')
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       'inverse_ratio = abs(u1)>abs(u2)')
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       'ratio = u2/u1 if inverse_ratio else u1/u2')
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       'rr = abs(ratio)**2')
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       'b = A*(1-rr)')
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       'root = mp.sqrt(b*b+4*B*rr)')
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       'small = ratio*mp.sqrt((2*B if inverse_ratio else 2)/(b+root))')
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       'inf = inverse_ratio and abs(small)<mp.mpf(\'.5\')')
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       'tc = small if inf else t')
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       'phase_coordinate = -mp.conj(u1) if inf else u2')
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       'a = mp.arg(phase_coordinate)/(2*mp.pi)')
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       'nu = mp.exp(-2) if inf else mp.mpf(1)')
    and has_if_test9(cusp8.elliptic_raw_base_reference,
                     'abs(tc)<=mp.mpf(\'.3\') and elementary')
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       'xs,yy = elementary_cusp_reference(tc,a,nu)')
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       "return ('ci' if inf else 'c0',(xs,yy,a))")
    and has_statement9(cusp8.elliptic_pole_raw_inverse,
                       'aa,bb,cc = v[0],-v[2],-v[1]')
    and has_statement9(cusp8.elliptic_pole_raw_inverse,
                       'u1 = mp.mpc(-1-2*bb,2*aa)/den')
    and has_statement9(cusp8.elliptic_pole_raw_inverse,
                       'u2 = mp.mpc(2*(bb+s2),2*cc)/den')
    and has_statement9(cusp8.elliptic_pole_raw_inverse,
                       'u1 = mp.mpc(1-s2,2*aa)/(1+s2)')
    and has_statement9(cusp8.elliptic_pole_raw_inverse,
                       'u2 = mp.mpc(2*(bb+1),2*cc)/(1+s2)')
    and has_statement9(cusp8.elliptic_pole_raw_inverse,
                       'ratio = u1/u2')
    and has_statement9(cusp8.elliptic_pole_raw_inverse,
                       'ratio = u2/u1')
    and has_statement9(cusp8.elliptic_pole_raw_inverse,
                       'height = mp.arg(u2)/(2*mp.pi)')
    and has_statement9(cusp8.elliptic_pole_raw_inverse,
                       'height = -mp.arg(-u1)/(2*mp.pi)')
    and has_statement9(cusp8.elliptic_pole_raw_inverse,
                       'nu = mp.exp(-2) if infinity else mp.mpf(1)')
    and has_statement9(cusp8.elliptic_pole_raw_inverse,
                       'uu = abs(ratio)**2')
    and has_statement9(cusp8.elliptic_pole_raw_inverse,
                       'den = A*(1-uu)+mp.sqrt(A*A*(1-uu)**2+4*B*uu)')
    and has_statement9(cusp8.elliptic_pole_raw_inverse,
                       'tc = ratio*mp.sqrt((2*B if infinity else 2)/den)'))
# At h=|L|<1/10, south has 1+B>=1-h and
# 1+|y|^2=2+2B+h^2<=2+2h+h^2<221/100.  Hence
# Re(u2)=2(1+B)/(1+|y|^2)>180/221>1/sqrt(2).
# North has den=1+2B+2h^2>1-2h>4/5, and
# |B+h^2+iC|<=h+h^2<11/100.  Hence |u2|<11/40<1/sqrt(2),
# while Re(-u1)=(1+2B)/den>0.  Both source paths select the same
# dominant ratio.  The south phase coordinate u2 and the north coordinate
# -conj(u1) have positive real part, so the half-centred selector is off.
# On the source image |tc|=|zeta_0 zeta_1|/nu<=rho^2/(2nu).
# Here nu is 1 or e^-2>1/8 (I1a gives e^4<55<64), so
# |tc|<4*(2e-6)<.3.
# Hence the raw reference takes its c0/ci return and the north `inf` flag.
dominance9i = ((sp.Rational(180,221))**2>sp.Rational(1,2)
              and (sp.Rational(11,40))**2<sp.Rational(1,2)
              and 2*(1-sp.Rational(1,10))
                  /(2+2*sp.Rational(1,10)+sp.Rational(1,10)**2)
                  ==sp.Rational(180,221)
              and 2*(sp.Rational(1,10)+sp.Rational(1,10)**2)
                  /(1-2*sp.Rational(1,10))==sp.Rational(11,40)
              and 1-2*sp.Rational(1,10)>0
              and 55<64
              and 4*source_bound9<sp.Rational(3,10))
rep('N9i2 both source u formulas, dominant ratio and positive profile path',
    all(sp.factor(left-right)==0 for left,right in
        zip(unorth9i,unorth_source9i))
    and all(sp.factor(left-right)==0 for left,right in
            zip(usouth9i,usouth_source9i))
    and dominance9i and raw_base_source9i)

# On |height|<1/10, inverse_height_moment(height,q,0)=height.  Both
# reference branches and the raw inverse therefore use the same positive
# A=(R+height)/2, Y=(R-height)/2, with A*Y=q.  The raw inverse writes
# A or Y in rationalized form when the sign makes a subtraction small.
height9i = sp.Symbol('height9i',real=True)
q9i = sp.Symbol('q9i',positive=True)
R9i = sp.sqrt(height9i**2+4*q9i)
AA9i,YY9i = (R9i+height9i)/2,(R9i-height9i)/2
seed_source9i = (
    has_statement9(cusp8.elliptic_pole_raw_inverse,
                   'R = mp.sqrt(height*height+4*q)')
    and has_statement9(cusp8.elliptic_pole_raw_inverse,
                       'AA = 2*q/(R-height) if height<0 else (R+height)/2')
    and has_statement9(cusp8.elliptic_pole_raw_inverse,
                       'yy = 2*q/(R+height) if height>=0 else (R-height)/2')
    and has_statement9(cusp8.elliptic_pole_raw_inverse,
                       'moment,source_radius = xx-yy,xx+yy')
    and has_statement9(cusp8.elementary_cusp_reference,
                       'f = inverse_height_moment(b,q,0)')
    and has_statement9(cusp8.elementary_cusp_reference,
                       'cell = int(mp.floor(a+mp.mpf(\'.5\')))')
    and has_statement9(cusp8._terminal_height_data,
                       'kappa,amp,d = mp.mpf(2),mp.mpf(0),mp.mpf(1)')
    and has_statement9(cusp8.inverse_height_moment,
                       'if lam==1 or aa<=mp.mpf(\'.1\'):\n    return _height_G(a,amp,d)')
    and has_statement9(cusp8._height_G,
                       'if amp==0:\n    return a/d if inverse else d*a')
    and has_statement9(cusp8.elementary_cusp_reference,
                       'A = 2/(R-f) if cell==1 else (R+f)/2')
    and has_statement9(cusp8.elementary_cusp_reference,
                       'Y = 2/(R+f) if cell==-1 else (R-f)/2')
    and has_statement9(cusp8.elementary_cusp_reference,
                       'x = mp.sqrt(A/e)/nu')
    and has_statement9(cusp8.elementary_cusp_reference,
                       'return x,t/x')
    and has_statement9(cusp8.elementary_cusp_reference,
                       'y = mp.sqrt(Y)')
    and has_statement9(cusp8.elementary_cusp_reference,
                       'return t/y,y')
    and has_statement9(cusp8.elliptic_bundle_embedding,
                       'return (nu*mp.conj(x)/mp.sqrt(r),y/mp.sqrt(r))'))
rep('N9i3 the raw inverse and cusp reference have identical A,Y invariants',
    sp.simplify(AA9i-YY9i-height9i)==0
    and sp.simplify(AA9i*YY9i-q9i)==0
    and sp.simplify((R9i-height9i)*AA9i-2*q9i)==0
    and sp.simplify((R9i+height9i)*YY9i-2*q9i)==0
    and -sp.Rational(1,10)+sp.Rational(1,2)>0
    and sp.Rational(1,10)+sp.Rational(1,2)<1
    and seed_source9i)

# The source's w0,w1 numerator is precisely the rank-one projector of the
# active cusp frame zeta=(nu*conj(x),y).  The same projector is P_j in that
# two-component block.  Since the already certified Qz_L is nonzero,
# Qz_L/|Qz_L|=e_ref*(e_ref^*z_L)/|e_ref^*z_L| exactly.
xr9i,xi9i,yr9i,yi9i = sp.symbols('xr9i xi9i yr9i yi9i',real=True)
nu9i = sp.symbols('nu9i',positive=True)
x9i,y9i = xr9i+sp.I*xi9i,yr9i+sp.I*yi9i
tc9i = x9i*y9i
zeta9i = sp.Matrix([nu9i*sp.conjugate(x9i),y9i])
rho29i = sp.expand((zeta9i.conjugate().T*zeta9i)[0])
Q9i = zeta9i*zeta9i.conjugate().T/rho29i
xx9i,yy9i = (nu9i**2*x9i*sp.conjugate(x9i),
              y9i*sp.conjugate(y9i))
source_pc9i = 2*nu9i*sp.conjugate(tc9i)
moment9i = xx9i-yy9i
z09i,z19i = sp.symbols('z09i z19i')
input9i = sp.Matrix([z09i,z19i])
raw_numerator9i = sp.Matrix([
    ((rho29i+moment9i)*z09i+source_pc9i*z19i)/(2*rho29i),
    (sp.conjugate(source_pc9i)*z09i+(rho29i-moment9i)*z19i)/(2*rho29i)])
projector_source9i = (
    has_statement9(cusp8.elliptic_pole_raw_inverse,
                   'source_pc = 2*nu*mp.conj(tc)')
    and has_statement9(cusp8.elliptic_pole_raw_inverse,
                       'w0 = ((source_radius+moment)*z0+source_pc*z1)/(2*source_radius)')
    and has_statement9(cusp8.elliptic_pole_raw_inverse,
                       'w1 = (mp.conj(source_pc)*z0+(source_radius-moment)*z1)/(2*source_radius)')
    and has_statement9(cusp8.elliptic_pole_raw_inverse,
                       'factor = mp.sqrt(source_radius)/length')
    and has_statement9(cusp8.elliptic_meridian_projector,
                       'return mp.matrix([[a*mp.conj(b) for b in e] for a in e])'))
# The actual unit phase is c/|c|, c=zeta* z_L.  We verify its real and
# imaginary decomposition from this very frame, then normalize a formal
# nonzero c=(c_r+i c_i); c=0 is excluded by the certified Qz_L bound.
cactual9i = (zeta9i.conjugate().T*z9i)[0]
cactual_re9i = sp.expand(sp.re(cactual9i))
cactual_im9i = sp.expand(sp.im(cactual9i))
Qz_actual9i = Q9i*z9i
Qz_norm29i = sum(sp.conjugate(v)*v for v in Qz_actual9i)
cr9i,cim9i = sp.symbols('cr9i cim9i',real=True)
cformal9i = cr9i+sp.I*cim9i
cnorm9i = sp.sqrt(cr9i**2+cim9i**2)
lam9i = cformal9i/cnorm9i
rho_symbol9i = sp.Symbol('rho_symbol9i',positive=True)
formal_Qz9i = zeta9i*cformal9i/rho_symbol9i
formal_norm9i = sp.sqrt((cr9i**2+cim9i**2)/rho_symbol9i)
normalized_raw9i = sp.sqrt(rho_symbol9i)*formal_Qz9i/formal_norm9i
rawspin9i = lam9i*zeta9i
phase_from_inner9i = (zeta9i.conjugate().T*rawspin9i)[0]/rho29i
slab_pair9i = (sp.conjugate(lam9i)*x9i,lam9i*y9i)
cap_pair9i = (sp.conjugate(rawspin9i[0])/nu9i,rawspin9i[1])
rep('N9i4 actual raw Qz numerator, active frame and cusp phase action',
    all(sp.cancel(v)==0 for v in Q9i*input9i-raw_numerator9i)
    and all(sp.cancel(v)==0 for v in Q9i*input9i
            -zeta9i*(zeta9i.conjugate().T*input9i)[0]/rho29i)
    and sp.expand(cactual9i-cactual_re9i-sp.I*cactual_im9i)==0
    and sp.expand(cactual9i*sp.conjugate(cactual9i)
                  -cactual_re9i**2-cactual_im9i**2)==0
    and sp.cancel(Qz_norm29i
                  -cactual9i*sp.conjugate(cactual9i)/rho29i)==0
    and sp.simplify(lam9i*sp.conjugate(cformal9i)/cnorm9i-1)==0
    and all(sp.simplify(a-b)==0 for a,b in zip(normalized_raw9i,rawspin9i))
    and sp.simplify(phase_from_inner9i-lam9i)==0
    and all(sp.simplify(a-b)==0 for a,b in zip(slab_pair9i,cap_pair9i))
    and projector_source9i)

print('ALL PASS' if ok else 'FAILURES')
assert ok, 'A symbolic identity, constant bound, or numerical check failed'
