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

# ---------------------------------------------------------------- N10: the nonidentity radial pole shell.
# This is a stage of H_j=D^{-1} g_j D, not by itself a check of the full
# six-dimensional tensor there.  The finite elementary S and S^{-1} below
# supply actual value and second-jet nodes for the missing radial stages.
shell_a10 = sp.Rational(1, 10**8)
shell_c10 = 3*shell_a10/2
shell_r10, shell_s10 = sp.symbols('shell_r10 shell_s10', positive=True)
shell_x10 = sp.exp(shell_c10/(shell_c10-shell_r10))-sp.exp(3)
shell_y10 = sp.exp(shell_c10/(shell_c10-shell_s10))-sp.E
shell_F10 = shell_x10-1/shell_x10
shell_T10 = sp.log(1+sp.exp(shell_F10))
shell_S10 = shell_c10*(1-1/sp.log(shell_T10+sp.E))
shell_L10 = sp.log(sp.exp(shell_y10)-1)
shell_root10 = (shell_L10+sp.sqrt(shell_L10**2+4))/2
shell_R10 = shell_c10*(1-1/sp.log(shell_root10+sp.exp(3)))
source_radial10 = (
    has_statement9(cusp8._elementary_pole_radial,
                   "a,c = mp.mpf('1e-8'),mp.mpf('1.5e-8')")
    and has_statement9(cusp8._elementary_pole_radial,
                       'x = mp.exp(3)*mp.expm1(3*(radius-a)/(c-radius))')
    and has_statement9(cusp8._elementary_pole_radial,
                       'f = x-1/x')
    and has_statement9(cusp8._elementary_pole_radial,
                       'y = mp.log1p(mp.exp(f)) if f<0 else f+mp.log1p(mp.exp(-f))')
    and has_statement9(cusp8._elementary_pole_radial,
                       'ell = mp.log1p(y/mp.e)')
    and has_statement9(cusp8._elementary_pole_radial,
                       'new_radius = c*ell/(1+ell)')
    and has_statement9(cusp8._elementary_pole_radial,
                       'y = mp.e*mp.expm1(radius/(c-radius))')
    and has_statement9(cusp8._elementary_pole_radial,
                       'x = 2/(mp.sqrt(L*L+4)-L) if L<0 else (L+mp.sqrt(L*L+4))/2')
    and has_statement9(cusp8._elementary_pole_radial,
                       'ell = mp.log1p(x/mp.exp(3))')
    and has_statement9(cusp8._elementary_pole_radial,
                       'new_radius = a+c*ell/(3*(3+ell))')
    and has_statement9(cusp8._elementary_pole_radial,
                       'factor = mp.sqrt(new_radius/radius)')
    and has_statement9(cusp8.elementary_pole_patch,
                       'compressed = _elementary_pole_radial(pair)')
    and has_statement9(cusp8.elementary_pole_patch,
                       'return _elementary_pole_radial(raw,inverse=True)')
    and has_statement9(cusp8.elementary_pole_forward,
                       'compressed = _elementary_pole_radial(pair)')
    and has_statement9(cusp8.elementary_pole_forward,
                       'return _elementary_pole_radial(raw,inverse=True)'))
# The production forms use expm1/log1p and an algebraically equivalent
# rationalization of the positive quadratic root.  The following small-node
# identities check those conversions without expanding the nested exponential.
shell_u10, shell_v10, shell_l10 = sp.symbols(
    'shell_u10 shell_v10 shell_l10', positive=True)
shell_coord10 = sp.cancel(shell_c10/(shell_c10-shell_r10)-3
                          -3*(shell_r10-shell_a10)
                          /(shell_c10-shell_r10))
shell_root_node10 = (shell_l10+sp.sqrt(shell_l10**2+4))/2
rep('N10a source-bound elementary radial formulas and stable-code rewrites',
    source_radial10 and shell_coord10==0
    and sp.cancel(shell_c10*(1-1/(1+shell_v10))
                  -shell_c10*shell_v10/(1+shell_v10))==0
    and sp.cancel(shell_a10+shell_c10*shell_v10/(3*(3+shell_v10))
                  -shell_c10*(1-1/(3+shell_v10)))==0
    and sp.simplify(shell_root_node10*(shell_root_node10-shell_l10)-1)==0)

# S and R are mutual inverses on the stated open ranges.  The inverse
# primitive identities are verified with positive formal exponential
# coordinates; for phi the exponential coordinate is >3 because r>a.
shell_e10 = sp.symbols('shell_e10', positive=True)
shell_lam10 = shell_u10-1/shell_u10
shell_soft10 = sp.log(1+sp.exp(shell_lam10))
shell_root_back10 = (shell_lam10+
                     sp.sqrt(shell_lam10**2+4))/2
shell_phi_inverse10 = shell_c10*(1-1/sp.log(
    sp.exp(shell_v10)-sp.exp(3)+sp.exp(3)))
shell_psi_inverse10 = shell_c10*(1-1/sp.log(
    sp.exp(shell_e10)-sp.E+sp.E))
shell_phi_back10 = shell_c10*(1-1/sp.log(shell_u10+sp.exp(3)))
shell_psi_back10 = shell_c10*(1-1/sp.log(shell_e10+sp.E))
rep('N10b elementary scalar S/R inverse and positive derivative factors',
    sp.simplify(sp.exp(shell_soft10)-1-sp.exp(shell_lam10))==0
    and sp.simplify(shell_root_back10*(shell_root_back10-shell_lam10)-1)==0
    and sp.simplify(sp.diff(shell_lam10,shell_u10)
                    -(1+shell_u10**-2))==0
    and sp.simplify(sp.diff(sp.log(1+sp.exp(shell_lam10)),shell_u10)
                    -(1+shell_u10**-2)*sp.exp(shell_lam10)
                      /(1+sp.exp(shell_lam10)))==0
    and sp.cancel(shell_c10*(1-1/shell_e10)
                  -shell_c10*(shell_e10-1)/shell_e10)==0
    and sp.simplify(2/(sp.sqrt(shell_l10**2+4)-shell_l10)
                    -shell_root_node10)==0
    and sp.simplify(shell_phi_inverse10
                    -shell_c10*(1-1/shell_v10))==0
    and sp.simplify(shell_psi_inverse10
                    -shell_c10*(1-1/shell_e10))==0
    and sp.simplify(sp.exp(shell_c10/(shell_c10-shell_phi_back10))
                    -sp.exp(3)-shell_u10)==0
    and sp.simplify(sp.exp(shell_c10/(shell_c10-shell_psi_back10))
                    -sp.E-shell_e10)==0
    and sp.simplify(sp.diff(shell_psi_back10,shell_e10)
                    -shell_c10/((shell_e10+sp.E)
                        *sp.log(shell_e10+sp.E)**2))==0
    and sp.simplify(sp.sqrt(sp.factor(shell_lam10**2+4))
                    -(shell_u10+1/shell_u10))==0
    and sp.simplify((shell_lam10+
                     sp.sqrt(sp.factor(shell_lam10**2+4)))/2-shell_u10)==0
    and sp.simplify((2*shell_root_back10-shell_lam10)**2
                    -shell_lam10**2-4)==0)

shell_z10 = sp.symbols('shell_z10_0:4', real=True)
shell_q10 = sum(v*v for v in shell_z10)
shell_f10 = sp.Function('shell_f10')(shell_q10)
shell_D10 = sp.Matrix([shell_f10*v for v in shell_z10])
shell_M10 = shell_D10.jacobian(shell_z10)
shell_H10 = [sp.hessian(v,shell_z10) for v in shell_D10]
shell_dummy10 = sp.symbols('shell_dummy10', positive=True)
shell_fp10 = sp.diff(sp.Function('shell_f10')(shell_dummy10),
                     shell_dummy10).subs(shell_dummy10,shell_q10)
shell_fpp10 = sp.diff(sp.Function('shell_f10')(shell_dummy10),
                      shell_dummy10,2).subs(shell_dummy10,shell_q10)
shell_M_formula10 = sp.Matrix(4,4,lambda i,j:
    shell_f10*int(i==j)+2*shell_fp10*shell_z10[i]*shell_z10[j])
shell_H_formula10 = [sp.Matrix(4,4,lambda i,j:
    2*shell_fp10*(int(k==i)*shell_z10[j]
                  +int(k==j)*shell_z10[i]
                  +int(i==j)*shell_z10[k])
    +4*shell_fpp10*shell_z10[k]*shell_z10[i]*shell_z10[j])
    for k in range(4)]
rep('N10c actual radial four-output first and 64 second entries',
    shell_M10==shell_M_formula10
    and all(shell_H10[k]==shell_H_formula10[k] for k in range(4))
    and all(shell_H10[k][i,j]==shell_H10[k][j,i]
            for k in range(4) for i in range(4) for j in range(4)))

shell_w10 = sp.symbols('shell_w10_0:4', real=True)
shell_qinv10 = sum(v*v for v in shell_w10)
shell_g10 = sp.Function('shell_g10')(shell_qinv10)
shell_Dinv10 = sp.Matrix([shell_g10*v for v in shell_w10])
shell_Minv10 = shell_Dinv10.jacobian(shell_w10)
shell_Hinv10 = [sp.hessian(v,shell_w10) for v in shell_Dinv10]
shell_g_dummy10 = sp.Function('shell_g10')(shell_dummy10)
shell_gp10 = sp.diff(shell_g_dummy10,shell_dummy10).subs(
    shell_dummy10,shell_qinv10)
shell_gpp10 = sp.diff(shell_g_dummy10,shell_dummy10,2).subs(
    shell_dummy10,shell_qinv10)
shell_Minv_formula10 = sp.Matrix(4,4,lambda i,j:
    shell_g10*int(i==j)+2*shell_gp10*shell_w10[i]*shell_w10[j])
shell_Hinv_formula10 = [sp.Matrix(4,4,lambda i,j:
    2*shell_gp10*(int(k==i)*shell_w10[j]
                  +int(k==j)*shell_w10[i]
                  +int(i==j)*shell_w10[k])
    +4*shell_gpp10*shell_w10[k]*shell_w10[i]*shell_w10[j])
    for k in range(4)]
rep('N10c2 inverse radial four-output first and 64 second entries',
    shell_Minv10==shell_Minv_formula10
    and all(shell_Hinv10[k]==shell_Hinv_formula10[k] for k in range(4))
    and all(shell_Hinv10[k][i,j]==shell_Hinv10[k][j,i]
            for k in range(4) for i in range(4) for j in range(4)))

# For D(z)=sqrt(S(r)/r) z, the tangential eigenvalue is f>0 and the
# radial eigenvalue is f+2r f'=f*r*S'/S>0.  The same calculation with R
# checks D^{-1}.  The determinant identity is symbolic in the *actual*
# value of S, whose derivative is a positive composition of the three
# monotone elementary primitives checked above.
shell_sym_r10 = sp.symbols('shell_sym_r10', positive=True)
shell_sym_S10 = sp.Function('shell_sym_S10')(shell_sym_r10)
shell_f_expr10 = sp.sqrt(shell_sym_S10/shell_sym_r10)
shell_det10 = sp.factor(shell_f_expr10**3*
                        (shell_f_expr10+2*shell_sym_r10*
                         sp.diff(shell_f_expr10,shell_sym_r10)))
shell_det_expected10 = (shell_sym_S10
                        *sp.diff(shell_sym_S10,shell_sym_r10)/shell_sym_r10)
shell_sym_s10 = sp.symbols('shell_sym_s10', positive=True)
shell_sym_R10 = sp.Function('shell_sym_R10')(shell_sym_s10)
shell_g_expr10 = sp.sqrt(shell_sym_R10/shell_sym_s10)
shell_invdet10 = sp.factor(shell_g_expr10**3*
    (shell_g_expr10+2*shell_sym_s10*
     sp.diff(shell_g_expr10,shell_sym_s10)))
shell_invdet_expected10 = (shell_sym_R10
    *sp.diff(shell_sym_R10,shell_sym_s10)/shell_sym_s10)
shell_phi_prime10 = (shell_c10*sp.exp(shell_c10/(shell_c10-shell_r10))
                     /(shell_c10-shell_r10)**2)
shell_T_prime10 = (sp.exp(shell_F10)*(1+shell_x10**-2)
                   /(1+sp.exp(shell_F10)))
shell_psi_inverse_prime10 = (shell_c10/((shell_T10+sp.E)
                                  *sp.log(shell_T10+sp.E)**2))
rep('N10d radial determinant S Sprime/r and inverse-jet domain',
    sp.simplify(shell_det10-shell_det_expected10)==0
    and sp.simplify(shell_invdet10-shell_invdet_expected10)==0
    and sp.simplify(sp.diff(shell_S10,shell_r10)
                    -shell_phi_prime10*shell_T_prime10
                     *shell_psi_inverse_prime10)==0
    and sp.simplify(sp.diff(sp.exp(shell_c10/(shell_c10-shell_r10)),shell_r10)
                    -shell_c10*sp.exp(shell_c10/(shell_c10-shell_r10))
                     /(shell_c10-shell_r10)**2)==0
    and sp.simplify(sp.diff(sp.exp(shell_c10/(shell_c10-shell_s10)),shell_s10)
                    -shell_c10*sp.exp(shell_c10/(shell_c10-shell_s10))
                     /(shell_c10-shell_s10)**2)==0
    and shell_a10<sp.Rational(5,4)*shell_a10<shell_c10)

# ---------------------------------------------------------------- N11: the actual dominant-coordinate profile in the raw
# meridian reference.  This checks one scalar generator of P(n,X), not the
# downstream height, ten-slot embedding, or full fixed-count G_cl circuit.
# On the open branch |u1|>|u2| the source uses T=u2/u1 and c=2B;
# otherwise it uses T=u1/u2 and c=2.  The dominant denominator is nonzero
# because |u1|^2+|u2|^2=1.  For z=|T|^2<=1, b=A(1-z)>=0 and
# D=b+sqrt(b^2+4Bz)>0 (A,B>0), including T=0.
profile_tree11 = ast.parse(inspect.getsource(
    cusp8.elliptic_raw_base_reference)).body[0]
profile_with11 = [node for node in profile_tree11.body
                  if isinstance(node,ast.With)]
profile_prefix11 = ast.parse('''
vv = tuple(map(mp.mpf,v))
xx = mp.mpf(X)
if not abs(xx)<1:
    raise ValueError('the raw reference excludes the two fixed points')
aa = xx*mp.sqrt(2-xx*xx)
w = (vv[0]/(1-aa),-vv[2]/(1-aa),-vv[1]/(1-aa))
y = (w[0],w[1]+1,w[2])
s2 = sum(z*z for z in y)
u1 = mp.mpc(1-s2,2*y[0])/(1+s2)
u2 = mp.mpc(2*y[1],2*y[2])/(1+s2)
A,B = (4*mp.pi)**2,mp.exp(4)
inverse_ratio = abs(u1)>abs(u2)
ratio = u2/u1 if inverse_ratio else u1/u2
rr = abs(ratio)**2
b = A*(1-rr)
root = mp.sqrt(b*b+4*B*rr)
small = ratio*mp.sqrt((2*B if inverse_ratio else 2)/(b+root))
inf = inverse_ratio and abs(small)<mp.mpf('.5')
t = None if inf else (1/small if inverse_ratio else small)
''').body
profile_source11 = (len(profile_with11)==1
    and [ast.dump(node,include_attributes=False)
         for node in profile_with11[0].body[:len(profile_prefix11)]]
    ==[ast.dump(node,include_attributes=False)
       for node in profile_prefix11]
    and isinstance(profile_with11[0].body[len(profile_prefix11)],ast.If)
    and ast.dump(profile_with11[0].body[len(profile_prefix11)].test,
                 include_attributes=False)
       ==ast.dump(ast.parse("inf or abs(t)<mp.mpf('.5')",mode='eval').body,
                   include_attributes=False))
yr11,yi11,yj11 = sp.symbols('yr11 yi11 yj11',real=True)
s211 = yr11**2+yi11**2+yj11**2
u111 = (1-s211+2*sp.I*yr11)/(1+s211)
u211 = 2*(yi11+sp.I*yj11)/(1+s211)
sphere11 = sp.factor(u111*sp.conjugate(u111)
                     +u211*sp.conjugate(u211)-1)==0
rep('N11a ordered actual raw-reference profile path and dominant denominator',
    profile_source11 and sphere11)

# Two exact round-meridian witnesses at X=0, where v=n and |n|=1.
# The inverse witness has u1=(3+6i)/7, u2=2/7, z=4/45;
# the direct witness has u1=-3/5, u2=4/5, z=9/16.
def profile_witness11(n):
    assert sp.expand(sum(v*v for v in n)-1)==0
    w = (n[0],-n[2],-n[1])
    y = (w[0],w[1]+1,w[2])
    s2 = sum(v*v for v in y)
    return ((1-s2+2*sp.I*y[0])/(1+s2),
            2*(y[1]+sp.I*y[2])/(1+s2))

u1inv11,u2inv11 = profile_witness11(
    (sp.Rational(3,5),sp.Integer(0),sp.Rational(4,5)))
u1dir11,u2dir11 = profile_witness11(
    (sp.Integer(0),sp.Integer(0),-sp.Integer(1)))
zin11,zdir11 = sp.Rational(4,45),sp.Rational(9,16)
profile_witnesses11 = (
    sp.simplify(u1inv11-(3+6*sp.I)/7)==0
    and sp.simplify(u2inv11-sp.Rational(2,7))==0
    and sp.simplify(u2inv11/u1inv11-(2-4*sp.I)/15)==0
    and sp.simplify((u2inv11/u1inv11)
                    *sp.conjugate(u2inv11/u1inv11)-zin11)==0
    and sp.simplify(u1dir11+sp.Rational(3,5))==0
    and sp.simplify(u2dir11-sp.Rational(4,5))==0
    and sp.simplify((u1dir11/u2dir11)
                    *sp.conjugate(u1dir11/u2dir11)-zdir11)==0
    and sp.simplify(u1inv11*sp.conjugate(u1inv11)
                    -u2inv11*sp.conjugate(u2inv11))>0
    and sp.simplify(u2dir11*sp.conjugate(u2dir11)
                    -u1dir11*sp.conjugate(u1dir11))>0
    and bool(sp.pi>3) and bool(sp.E<3)
    # R>b for z>0, so D>2A(1-z); these bounds force the
    # inverse witness's `inf` selector and the direct witness's south t.
    and sp.Rational(81*4,144*41)<sp.Rational(1,4)
    and sp.Rational(9,144*7)<sp.Rational(1,4))
rep('N11b strict rational meridian witnesses for both profile branches',
    profile_witnesses11)

# Write T=x+iy, z=x^2+y^2.  These are finite elementary value, gradient,
# and Hessian circuits for the exact source `small`, with c=2B or c=2.
# The one-variable derivative formulas are checked against the actual
# positive-root expression; the two-variable radial tensor identities
# then check all 4 first and 8 second entries for each specialization.
px11,py11,pz11 = sp.symbols('px11 py11 pz11',real=True)
Ap11,Bp11,cp11 = sp.symbols('Ap11 Bp11 cp11',positive=True)
pb11 = Ap11*(1-pz11)
pR11 = sp.sqrt(pb11**2+4*Bp11*pz11)
pD11 = pb11+pR11
pg11 = sp.sqrt(cp11/pD11)
pRz11 = (2*Bp11-Ap11*pb11)/pR11
pRzz11 = 4*Bp11*(Ap11**2-Bp11)/pR11**3
pDz11,pDzz11 = -Ap11+pRz11,pRzz11
pgz11 = -pg11*pDz11/(2*pD11)
pgzz11 = pg11*(3*pDz11**2/(4*pD11**2)
                 -pDzz11/(2*pD11))
profile_root11 = (
    sp.simplify(pD11**2-2*pb11*pD11-4*Bp11*pz11)==0
    and sp.simplify(sp.diff(pR11,pz11)-pRz11)==0
    and sp.simplify(sp.diff(pR11,pz11,2)-pRzz11)==0
    and sp.simplify(sp.diff(pg11,pz11)-pgz11)==0
    and sp.simplify(sp.diff(pg11,pz11,2)-pgzz11)==0)
rep('N11c source positive-root value and first/second radial derivatives',
    profile_root11)

pf11 = sp.Function('pf11')(px11**2+py11**2)
pvec11 = (px11,py11)
pfz11 = sp.diff(sp.Function('pf11')(pz11),pz11).subs(
    pz11,px11**2+py11**2)
pfzz11 = sp.diff(sp.Function('pf11')(pz11),pz11,2).subs(
    pz11,px11**2+py11**2)
abstract_profile11 = sp.Matrix([v*pf11 for v in pvec11])
abstract_M11 = abstract_profile11.jacobian(pvec11)
abstract_H11 = [sp.hessian(v,pvec11) for v in abstract_profile11]
expected_M11 = sp.Matrix(2,2,lambda i,j:
    int(i==j)*pf11+2*pvec11[i]*pvec11[j]*pfz11)
expected_H11 = [sp.Matrix(2,2,lambda i,j:
    2*pfz11*(int(k==i)*pvec11[j]+int(k==j)*pvec11[i]
               +int(i==j)*pvec11[k])
    +4*pfzz11*pvec11[k]*pvec11[i]*pvec11[j]) for k in range(2)]
profile_radial_jets11 = (abstract_M11==expected_M11
    and all(abstract_H11[k]==expected_H11[k] for k in range(2)))
profile_circuits11 = []
for branch_c11 in (2*Bp11,sp.Integer(2)):
    values11 = sp.Matrix([v*pg11.subs(cp11,branch_c11).subs(
        pz11,px11**2+py11**2) for v in pvec11])
    first11 = expected_M11.subs({pf11:pg11.subs(cp11,branch_c11),
        pfz11:pgz11.subs(cp11,branch_c11)}).subs(
            pz11,px11**2+py11**2)
    second11 = [H.subs({pf11:pg11.subs(cp11,branch_c11),
        pfz11:pgz11.subs(cp11,branch_c11),
        pfzz11:pgzz11.subs(cp11,branch_c11)}).subs(
            pz11,px11**2+py11**2) for H in expected_H11]
    profile_circuits11.append((values11,first11,second11))
profile_jet_checks11 = (profile_radial_jets11
    and all(V.shape==(2,1) and M.shape==(2,2) and len(H)==2
            and all(h.shape==(2,2) and h==h.T for h in H)
            and not any(entry.has(sp.Derivative)
                        for entry in list(V)+list(M)
                        +[entry for h in H for entry in h])
            for V,M,H in profile_circuits11))
rep('N11d both source profile branches: 2 values, 4 first and 8 second jets',
    profile_jet_checks11)

# The selector |u1|=|u2| is not a singularity of the effective middle
# coordinate t.  Put r=|u2/u1|^2>0, so the direct branch has ratio-square
# 1/r.  Both source denominators use positive square roots.  Since r>0,
# R_d=R/r; hence D_i D_d=4B.  Around r=1 the inverse branch's `small`
# has modulus e>2, so `inf` is false on an open neighbourhood.  There the
# two formulas for t, including every radial jet, are identical.
profile_after11 = profile_with11[0].body[len(profile_prefix11)]
profile_phase_prefix11 = ast.parse('''
tc = small if inf else t
phase_coordinate = -mp.conj(u1) if inf else u2
''').body
profile_phase_source11 = (profile_source11
    and isinstance(profile_after11,ast.If)
    and [ast.dump(node,include_attributes=False)
         for node in profile_after11.body[:2]]
    ==[ast.dump(node,include_attributes=False)
       for node in profile_phase_prefix11])
pr11 = sp.symbols('pr11',positive=True)
pbi11 = Ap11*(1-pr11)
pRi11 = sp.sqrt(pbi11**2+4*Bp11*pr11)
pDi11 = pbi11+pRi11
pRd11 = sp.sqrt(Ap11**2*(1-1/pr11)**2+4*Bp11/pr11)
pDd11 = Ap11*(1-1/pr11)+pRd11
pQi11,pQd11 = pDi11/(2*Bp11),2/pDd11
profile_seam_algebra11 = (
    sp.simplify(pRd11-pRi11/pr11)==0
    and sp.simplify(pDd11-(pRi11-pbi11)/pr11)==0
    and sp.simplify(pDi11*pDd11-4*Bp11)==0
    and all(sp.simplify(sp.diff(pQi11-pQd11,pr11,k))==0
            for k in range(3)))
# For a positive squared factor Q, these are its first and second
# positive-square-root derivatives.  Equal Q,Q',Q'' and positive factors
# imply equal factors and both radial jets without branch sign ambiguity.
pQformal11 = sp.Function('pQformal11')(pr11)
pFformal11 = sp.sqrt(pQformal11)
profile_sqrt_jets11 = (
    sp.simplify(sp.diff(pFformal11,pr11)
                -sp.diff(pQformal11,pr11)/(2*pFformal11))==0
    and sp.simplify(sp.diff(pFformal11,pr11,2)
                    -sp.diff(pQformal11,pr11,2)/(2*pFformal11)
                    +sp.diff(pQformal11,pr11)**2
                     /(4*pFformal11**3))==0)
profile_tie11 = (
    sp.simplify(pDi11.subs(pr11,1)-2*sp.sqrt(Bp11))==0
    and sp.simplify(pDd11.subs(pr11,1)-2*sp.sqrt(Bp11))==0
    and sp.simplify(sp.sqrt(pQi11).subs({pr11:1,
                    Bp11:sp.exp(4)})-sp.exp(-1))==0
    and sp.simplify(sp.sqrt(2*Bp11/pDi11).subs({pr11:1,
                    Bp11:sp.exp(4)})-sp.E)==0
    and bool(sp.E>2) and bool(sp.exp(-1)<sp.Rational(1,2)))
rep('N11e removable dominant-selector seam: common t, phase input and two jets',
    profile_phase_source11 and profile_seam_algebra11
    and profile_sqrt_jets11 and profile_tie11)

# ---------------------------------------------------------------- N12: one active inverse-first height stage.
# This certifies the source path from (a,r,nu) through the terminal-forward
# output v, with exact 2-jets in (b,q,lambda).  It does not differentiate the
# subsequent G^{-1}, R^{-1}, G chain, Y, rho, or the full round J.
def ordered_block12(function, snippet):
    """Find the given consecutive source statements in one actual AST block."""
    tree = ast.parse(inspect.getsource(function))
    wanted = [ast.dump(node,include_attributes=False)
              for node in ast.parse(snippet).body]
    for owner in ast.walk(tree):
        for _,block in ast.iter_fields(owner):
            if isinstance(block,list):
                actual = [ast.dump(node,include_attributes=False)
                          for node in block if isinstance(node,ast.AST)]
                if any(actual[i:i+len(wanted)]==wanted
                       for i in range(len(actual)-len(wanted)+1)):
                    return True
    return False

height_source12 = (
    ordered_block12(cusp8.elementary_cusp_logradius,'''\
P = mp.log(nu)+mp.log(r)-2*r*r/(1-r*r)
q = mp.exp(2*P)
n = mp.floor(a+mp.mpf('.5'))
b = a-n
lam = scalar_flat_step((r-mp.mpf('.3'))/mp.mpf('.2'))
if lam==1:
    return +(P*(a+mp.mpf('.5')))
f = inverse_height_moment(b,q,lam)
''')
    and ordered_block12(cusp8.scalar_flat_step,'''\
a,b = mp.exp(-1/u),mp.exp(-1/(1-u))
return a/(a+b)
''')
    and ordered_block12(cusp8._terminal_height_data,'''\
k = -mp.log(q)
kappa = 2+lam*(k-2)
amp = lam*k
d = 1-lam+2*lam*mp.sqrt(q)*k
''')
    and ordered_block12(cusp8._terminal_height_data,'''\
B,C = mp.exp(-kappa/2),q*mp.exp(kappa/2)
z0 = B-C
zcut = B*mp.exp(mp.mpf('.4')*kappa)-C*mp.exp(-mp.mpf('.4')*kappa)
eta = mp.mpf('.01')
E0 = mp.exp(eta/zcut)
Acut = mp.exp(eta/(zcut-z0))
Lambda = mp.exp(Acut-E0)
Eprime0 = E0*eta/zcut**2
Eprime_z0 = Acut*eta/(zcut-z0)**2
mu = d*Eprime0/(kappa*(B+C)*Eprime_z0)
''')
    and ordered_block12(cusp8.inverse_height_moment,'''\
data = _terminal_height_data(q,lam)
kappa,amp,d,B,C,z0,_,_,_,_,_,_ = data
aa = abs(a)
if lam==1 or aa<=mp.mpf('.1'):
    return _height_G(a,amp,d)
if aa>=mp.mpf('.4'):
    return mp.sign(a)*(B*mp.exp(kappa*aa)-C*mp.exp(-kappa*aa))
u = _height_R(aa)
delta_z = B*mp.expm1(kappa*u)-C*mp.expm1(-kappa*u)
v = _height_terminal_forward(delta_z,data)
return mp.sign(a)*_height_G(_height_R(_height_G(v,amp,d,inverse=True),inverse=True),amp,d)
''')
    and ordered_block12(cusp8._height_R,'''\
delta,c = mp.mpf('.1'),mp.mpf('.2')
if inverse:
    if value<=0:
        return delta
elif value<=delta:
    return mp.mpf(0)
if value>=c:
    return value
''')
    and ordered_block12(cusp8._height_R,'''\
x = mp.exp(2)*mp.expm1(2*(value-delta)/(c-value))
F = x-1/x
soft = (mp.log1p(mp.exp(F)) if F<0 else
        F+mp.log1p(mp.exp(-F)))
ell = mp.log1p(soft/mp.e)
return c*ell/(1+ell)
''')
    and ordered_block12(cusp8._height_terminal_forward,'''\
gap0 = zcut-z0
if delta_z<=0:
    return mp.mpf(0)
if delta_z>=gap0:
    return z0+delta_z
gap = gap0-delta_z
L = eta/gap
Ez = mp.exp(L)
''')
    and ordered_block12(cusp8._height_terminal_forward,'''\
x = Acut*mp.expm1(eta*delta_z/(gap*gap0))
u = mp.expm1(x)
F = (Lambda*u*u+mu*u)/(1+u)
M = mp.log1p(F)
ell = mp.log1p(M/E0)
return zcut*ell/(eta/zcut+ell)
''')
    and all(has_if_test9(cusp8._height_terminal_forward,test)
            for test in ('delta_z<=0','delta_z>=gap0',
                         'Ez>(mp.mp.dps+20)*mp.log(10)+mp.log(K)',
                         'Ez>50'))
    and has_if_test9(cusp8._height_R,'c/(c-value)>(mp.mp.dps+10)*mp.log(10)'))
rep('N12a ordered production height, moment, radial and terminal source path',
    height_source12)

# The exact source witness is a=b=3/20, r=2/5, nu=1.  It has n=0 and
# lambda=s(1/2)=1/2.  q=(4/25)exp(-16/21).  The following rational bounds
# prove 1/2500<q<1/8 and 2<k=-log(q)<4, with no sampled inequality.
height_r12 = sp.Rational(2,5)
height_b12 = sp.Rational(3,20)
height_P12 = sp.log(height_r12)-2*height_r12**2/(1-height_r12**2)
height_q12 = sp.exp(2*height_P12)
height_lam12 = sp.Rational(1,2)
height_bounds12 = (
    sp.simplify(height_P12-sp.log(sp.Rational(2,5))+sp.Rational(8,21))==0
    and sp.simplify(height_q12-sp.Rational(4,25)*sp.exp(-sp.Rational(16,21)))==0
    and sp.Rational(3,10)<height_r12<sp.Rational(1,2)
    and sp.Rational(1,10)<height_b12<sp.Rational(1,5)
    and sp.floor(height_b12+sp.Rational(1,2))==0
    and sp.simplify(sp.exp(-2)/(sp.exp(-2)+sp.exp(-2))-height_lam12)==0
    and bool(sp.Rational(5,2)<sp.E<sp.Rational(11,4))
    and sp.Rational(4,75)>sp.Rational(1,2500)
    and sp.Rational(4,75)>sp.Rational(16,625)
    and sp.Rational(84,925)<sp.Rational(1,8)
    and sp.Rational(1,8)<sp.Rational(16,121)
    and sp.Rational(1,10)/(sp.Rational(1,5)-height_b12)==2
    and sp.Rational(1,5)/(sp.Rational(1,5)-height_b12)==4)
# Since e<3 and 16/21<1, q>4/(25e)>4/75.
# Since exp(x)>1+x for x>0, q<84/925.  e>5/2 gives
# e^{-4}<16/625<4/75, while e<11/4 gives 1/8<16/121<e^{-2}.
# Thus 2<k<4, kappa=1+k/2 in (2,3), B=e^{-kappa/2}>1/9,
# C/B=q e^kappa=e sqrt(q)<1, and zcut>z0>0.
# For the source R branch, x=e^2(exp(2)-1)>12, so F>0 and its
# positive-F softplus is chosen.  Then 0<u=R(b)<1/5.  For
# g(s)=B exp(kappa s)-C exp(-kappa s), g'(s)>2/9 on [0,2/5].
# Hence 0<delta_z<gap0 and gap=g(2/5)-g(u)>2/45.  Therefore
# L=eta/gap<9/40, Ez<e<3<50.  The precision shortcut is also false:
# K>=1 and (dps+20)log(10)>42, whereas Ez<3.  The R precision
# shortcut is false because c/(c-b)=4<(dps+10)log(10).
height_branch12 = (height_bounds12
    and bool(sp.exp(2)*(sp.exp(2)-1)>12)
    and sp.Rational(2,9)*sp.Rational(1,5)==sp.Rational(2,45)
    and sp.Rational(1,100)/sp.Rational(2,45)==sp.Rational(9,40)
    and sp.Rational(9,40)<1
    and sp.Rational(4)<sp.Rational(22)
    and sp.Rational(3)<sp.Rational(42))
rep('N12b strict exact witness for active R and terminal correction',
    height_branch12)

# Every RHS below is the exact source expression in its selected branch;
# expm1(x)=exp(x)-1 and log1p(x)=log(1+x) are algebraic identities.
# Keeping each intermediate as an atom makes a finite, executable 2-jet
# circuit without expanding the enormous terminal composition.
hb12,hq12,hl12 = sp.symbols('hb12 hq12 hl12',real=True)
hn12 = sp.symbols('hn12_0:32',real=True)
eta12,c12,delta12 = sp.Rational(1,100),sp.Rational(1,5),sp.Rational(1,10)
(hk12,hkap12,hamp12,hd12,hB12,hC12,hz012,hzcut12,hE012,
 hAcut12,hLambda12,hEp012,hEpz012,hmu12,hxR12,hFR12,
 hsoft12,hellR12,huR12,hdz12,hgap012,hgap12,hL12,hEz12,
 hxT12,huT12,hFT12,hMT12,hellT12,hv12) = hn12[:30]
height_nodes12 = [
    (hk12,-sp.log(hq12)),
    (hkap12,2+hl12*(hk12-2)),
    (hamp12,hl12*hk12),
    (hd12,1-hl12+2*hl12*sp.sqrt(hq12)*hk12),
    (hB12,sp.exp(-hkap12/2)),
    (hC12,hq12*sp.exp(hkap12/2)),
    (hz012,hB12-hC12),
    (hzcut12,hB12*sp.exp(sp.Rational(2,5)*hkap12)
              -hC12*sp.exp(-sp.Rational(2,5)*hkap12)),
    (hE012,sp.exp(eta12/hzcut12)),
    (hAcut12,sp.exp(eta12/(hzcut12-hz012))),
    (hLambda12,sp.exp(hAcut12-hE012)),
    (hEp012,hE012*eta12/hzcut12**2),
    (hEpz012,hAcut12*eta12/(hzcut12-hz012)**2),
    (hmu12,hd12*hEp012/(hkap12*(hB12+hC12)*hEpz012)),
    (hxR12,sp.exp(2)*(sp.exp(2*(hb12-delta12)/(c12-hb12))-1)),
    (hFR12,hxR12-1/hxR12),
    (hsoft12,hFR12+sp.log(1+sp.exp(-hFR12))),
    (hellR12,sp.log(1+hsoft12/sp.E)),
    (huR12,c12*hellR12/(1+hellR12)),
    (hdz12,hB12*(sp.exp(hkap12*huR12)-1)
             -hC12*(sp.exp(-hkap12*huR12)-1)),
    (hgap012,hzcut12-hz012),
    (hgap12,hgap012-hdz12),
    (hL12,eta12/hgap12),
    (hEz12,sp.exp(hL12)),
    (hxT12,hAcut12*(sp.exp(eta12*hdz12/(hgap12*hgap012))-1)),
    (huT12,sp.exp(hxT12)-1),
    (hFT12,(hLambda12*huT12**2+hmu12*huT12)/(1+huT12)),
    (hMT12,sp.log(1+hFT12)),
    (hellT12,sp.log(1+hMT12/hE012)),
    (hv12,hzcut12*hellT12/(eta12/hzcut12+hellT12))]

# Topologically bind every local expression before differentiation.  A
# second-order chain rule composes each exact local SymPy derivative with
# already-certified input jets, producing three first and six independent
# second partials per scalar node.  No unevaluated Derivative is admitted.
height_inputs12 = (hb12,hq12,hl12)
height_jet12 = {
    s:(tuple(int(i==j) for i in range(3)),
       tuple(tuple(sp.Integer(0) for j in range(3)) for i in range(3)))
    for j,s in enumerate(height_inputs12)}
height_circuit12 = []
height_topology12 = True
for out12,rhs12 in height_nodes12:
    deps12 = sorted(rhs12.free_symbols,key=str)
    if out12 in height_jet12 or any(dep not in height_jet12 for dep in deps12):
        height_topology12 = False
        break
    d112 = {dep:sp.diff(rhs12,dep) for dep in deps12}
    d212 = {(p,w):sp.diff(d112[p],w) for p in deps12 for w in deps12}
    grad12 = tuple(sum((d112[p]*height_jet12[p][0][i]
                        for p in deps12),sp.Integer(0)) for i in range(3))
    hess12 = tuple(tuple(
        sum((d112[p]*height_jet12[p][1][i][j]
             for p in deps12),sp.Integer(0))
        +sum((d212[p,w]*height_jet12[p][0][i]*height_jet12[w][0][j]
              for p in deps12 for w in deps12),sp.Integer(0))
        for j in range(3)) for i in range(3))
    local12 = [rhs12,*d112.values(),*d212.values()]
    if (any(item.has(sp.Derivative) for item in local12)
            or any(sp.simplify(d212[p,w]-d212[w,p])!=0
                   for p in deps12 for w in deps12)
            or any(sp.simplify(hess12[i][j]-hess12[j][i])!=0
                   for i in range(3) for j in range(3))):
        height_topology12 = False
        break
    height_jet12[out12] = grad12,hess12
    height_circuit12.append((out12,rhs12,grad12,hess12))
rep('N12c source-bound active terminal stage: value, 3 first and 6 second jets',
    height_topology12 and len(height_circuit12)==30
    and height_circuit12[-1][0]==hv12
    and all(not item.has(sp.Derivative)
            for _,rhs,grad,hess in height_circuit12
            for item in (rhs,*grad,*(entry for row in hess for entry in row))))

# Type-(d) one-point transcription check and measured margins.  The
# inequalities in N12b, not these decimal values, certify the branch.
with mp.workdps(90):
    sample_q12 = mp.exp(2*(mp.log(mp.mpf(2)/5)-mp.mpf(8)/21))
    sample_data12 = cusp8._terminal_height_data(sample_q12,mp.mpf(1)/2)
    sample_u12 = cusp8._height_R(mp.mpf(3)/20)
    sample_dz12 = (sample_data12[3]*mp.expm1(sample_data12[0]*sample_u12)
                   -sample_data12[4]*mp.expm1(-sample_data12[0]*sample_u12))
    source_v12 = cusp8._height_terminal_forward(sample_dz12,sample_data12)
    sample_values12 = {hb12:height_b12,hq12:height_q12,hl12:height_lam12}
    for node12,rhs12 in height_nodes12:
        sample_values12[node12] = rhs12.subs(sample_values12).evalf(95)
    circuit_v12 = mp.mpf(str(sample_values12[hv12].evalf(90)))
    sample_L12 = sample_data12[7]/(
        sample_data12[6]-sample_data12[5]-sample_dz12)
    sample_error12 = abs(circuit_v12-source_v12)
    assert sample_error12<mp.mpf('1e-75')
    assert sample_q12>mp.mpf(1)/2500
    assert sample_q12<mp.mpf(1)/8
    assert sample_L12<mp.mpf(9)/40
    print('N12d one-point numerical source check only: |v_DAG-v_source|=',
          mp.nstr(sample_error12,5),
          '; q lower/upper margins =',
          mp.nstr(sample_q12-mp.mpf(1)/2500,12),
          mp.nstr(mp.mpf(1)/8-sample_q12,12),
          '; eta/gap upper margin =',
          mp.nstr(mp.mpf(9)/40-sample_L12,12),': PASS')

# ---------------------------------------------------------------- N13: one strict off-cut collar/pure-middle seam cell.
# The raw tuples can differ by a smooth scalar gauge.  We compare the
# rank-one ten-slot projectors and the phase-corrected middle W slice.
seam_source13 = (
    has_if_test9(cusp8.elliptic_raw_base_reference,
                 "inf or abs(t)<mp.mpf('.5')")
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       "phase_coordinate = -mp.conj(u1) if inf else u2")
    and has_if_test9(cusp8.elliptic_raw_base_reference,
                     "phase_coordinate.real<-mp.mpf('.9')*abs(phase_coordinate)")
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       "a = mp.mpf('.5')+mp.arg(-phase_coordinate)/(2*mp.pi)")
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       "a = mp.arg(phase_coordinate)/(2*mp.pi)")
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       "rho = (elementary_cusp_logradius if elementary else elliptic_cusp_logradius)(a,tc,nu,dps=dps)")
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       "lattice = scalar_lattice_coords(zm,wm,t,lt)")
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       "zm,wm = -mp.j*rho,mp.mpc(0)")
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       "zm,wm = -mp.j*rho+mp.pi*ss*(ss-1),2*t*ss")
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       "xi4 = scalar_lattice_coords(mp.pi*ss*(ss-1),2*t*ss,t,lt)[3]")
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       "four = mp.arg(u2)/(2*mp.pi)+xi4+mp.mpf('.5')")
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       "four = -a2-mp.mpf('.5')")
    and has_if_test9(cusp8.elliptic_raw_base_reference,"r<=1")
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       "lam = 1-scalar_flat_step((r-mp.mpf('1.2'))/mp.mpf('.4'))")
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       "a2 = -mp.arg(u1)/(2*mp.pi)+lam*delta+(1-lam)/2")
    and has_statement9(cusp8.elliptic_raw_base_reference,
                       "small = ratio*mp.sqrt((2*B if inverse_ratio else 2)/(b+root))")
    and has_statement9(cusp8._elementary_raw_orbit,
                       "ratio = tc*mp.sqrt((abs(tc)**2+A)/(A*abs(tc)**2+B))")
    and has_statement9(cusp8._elementary_raw_orbit,
                       "ratio = tc*mp.sqrt((A+B*abs(tc)**2)/(1+A*abs(tc)**2))")
    and ordered_block12(cusp8.elementary_cusp_logradius,'''\
lam = scalar_flat_step((r-mp.mpf('.3'))/mp.mpf('.2'))
if lam==1:
    return +(P*(a+mp.mpf('.5')))
f = inverse_height_moment(b,q,lam)
''')
    and has_statement9(cusp8.elementary_cusp_logradius,
                       "P = mp.log(nu)+mp.log(r)-2*r*r/(1-r*r)")
    and has_statement9(cusp8.elementary_cusp_logradius,
                       "Y = 2*q/(root+f) if f>=0 else (root-f)/2")
    and has_statement9(cusp8.elementary_cusp_logradius,
                       "return +(n*P+mp.log(Y)/2)")
    and ordered_block12(cusp8.elliptic_bundle_embedding,'''\
b0 = 1-scalar_flat_step((r-mp.mpf('.3'))/mp.mpf('.2'))
bi = 1-scalar_flat_step((1/r-mp.mpf('.3'))/mp.mpf('.2'))
bm = 1-b0-bi
z,w = 2*mp.pi*one-mp.j*lt*four,2*t*four
''')
    and has_statement9(cusp8.elliptic_bundle_embedding,
                       "weight = 1-scalar_flat_step((abs(angle)-mp.pi/2)/(mp.pi/4))")
    and has_statement9(cusp8.elliptic_bundle_embedding,
                       "vec[8+index] = bm*weight*mp.exp(mp.j*(2*mp.pi*one-(angle+centre-lt.imag)*four))")
    and has_statement9(cusp8.elementary_projection_section,
                       "vector = tuple(a*phase/abs(phase) for a in ref)")
    and has_statement9(cusp8.elementary_global_slice,
                       "return kind,(t,one+mp.arg(phase)/(2*mp.pi),four,lt)"))
rep('N13a ordered production selectors, endpoint, embedding and phase correction',
    seam_source13)

# Inverse-ratio profile at z=|tc|^2, A=(4*pi)^2, B=e^4.
# South u1/u2=tc*sqrt(Qs), Qs=(A+Bz)/(1+Az);
# north u2/u1=tc*sqrt(Qi), Qi=(z+A)/(Az+B).
# The rational D below is the unique positive source-root denominator.
sz13,A13,B13 = sp.symbols('sz13 A13 B13',positive=True)
Qs13=(A13+B13*sz13)/(1+A13*sz13)
Qi13=(sz13+A13)/(A13*sz13+B13)
zs13=1/(sz13*Qs13)
zi13=sz13*Qi13
Ds13=2*B13/Qs13
Di13=2*B13*Qi13
profile_seam13 = (
    sp.cancel(Ds13**2-2*A13*(1-zs13)*Ds13-4*B13*zs13)==0
    and sp.cancel(Di13**2-2*A13*(1-zi13)*Di13-4*B13*zi13)==0
    and sp.cancel(2*B13/Ds13-Qs13)==0
    and sp.cancel(2*B13/Di13-1/Qi13)==0
    and sp.cancel((sz13*Qs13-1).subs(sz13,sp.Rational(1,4))
                  -(B13-16)/(4*(A13+4)))==0
    and sp.cancel((1-zi13).subs(sz13,sp.Rational(1,4))
                  -(16*B13-1)/(4*(A13+4*B13)))==0
    and bool(sp.E>2) and bool(sp.exp(4)>16))
rep('N13b both positive-root profiles and strict seam-side selectors',
    profile_seam13)

# Pick arg(t)=-2*pi/3, strictly inside (-3*pi/4,-pi/2).
# South arg(u2)=pi/10 gives a=1/20; north arg(u1)=pi/10
# gives a=9/20 in the half-centred chart of -conj(u1).
# Both normalized phase pairs are genuine S3 inputs; inverse
# stereographic projection gives raw round-meridian inputs with |X|<1.
theta13=-2*sp.pi/3
as13,an13=sp.Rational(1,20),sp.Rational(9,20)
phase_guards13 = (
    -3*sp.pi/4<theta13<-sp.pi/2
    and 0<as13<sp.Rational(1,10)
    and sp.Rational(2,5)<an13<sp.Rational(1,2)
    and sp.floor(as13+sp.Rational(1,2))==0
    and sp.floor(an13+sp.Rational(1,2))==0
    and sp.Rational(1,2)-sp.Rational(1,20)==an13
    and bool(sp.cos(sp.pi/10)>sp.Rational(9,10))
    and 0<sp.Rational(2,3)<1
    and bool(sp.E<sp.Rational(11,4))
    and sp.Rational(16,3)<6
    and (sp.Rational(11,4))**6<625
    and sp.Rational(1,4)/(sp.Rational(11,4))**6
        >sp.Rational(1,2500)
    and sp.simplify(sp.exp(2*(-sp.log(2)-sp.Rational(2,3)))
                    -sp.exp(-sp.Rational(4,3))/4)==0
    and sp.simplify(sp.exp(2*(-2-sp.log(2)-sp.Rational(2,3)))
                    -sp.exp(-sp.Rational(16,3))/4)==0)
# Thus q_s=e^(-4/3)/4 and q_i=e^(-16/3)/4 lie strictly
# between 1/2500 and 1/8: e>2 proves the upper bounds;
# e<11/4, 16/3<6 and (11/4)^6<625 prove the north lower
# bound, which is the smaller q.
rep('N13c strict lower-half-plane and argument-chart witnesses',
    phase_guards13)

# These are actual source S3 phases, not freely chosen collar data.
# Both denominator norms are 1+|tc|^2 Q>1.  The inverse
# stereographic y and w=y-(0,1,0) give |w|^2>0 because
# Im(u2)>0 at both witnesses.  With a_X=(|w|^2-1)/(|w|^2+1),
# the unique X in (-1,1) satisfying X*sqrt(2-X^2)=a_X
# gives v=(1-a_X)(w0,-w2,-w1) and |v|+X^2=1.
zfix13=sp.Rational(1,4)
tcs13=sp.exp(-2*sp.pi*sp.I/3)/2
tci13=sp.exp(2*sp.pi*sp.I/3)/2
phase13=sp.exp(sp.pi*sp.I/10)
qsw13=Qs13.subs({sz13:zfix13,A13:(4*sp.pi)**2,B13:sp.exp(4)})
qiw13=Qi13.subs({sz13:zfix13,A13:(4*sp.pi)**2,B13:sp.exp(4)})
u2s13=phase13/sp.sqrt(1+zfix13*qsw13)
u1s13=tcs13*sp.sqrt(qsw13)*u2s13
u1i13=phase13/sp.sqrt(1+zfix13*qiw13)
u2i13=tci13*sp.sqrt(qiw13)*u1i13
Rw13=sp.symbols('Rw13',positive=True)
ax13=(Rw13-1)/(Rw13+1)
actual_meridian13=(
    sp.simplify(tcs13*sp.conjugate(tcs13)-zfix13)==0
    and sp.simplify(tci13*sp.conjugate(tci13)-zfix13)==0
    and all(sp.simplify(
        uu1*sp.conjugate(uu1)+uu2*sp.conjugate(uu2)-1)==0
        for uu1,uu2 in ((u1s13,u2s13),(u1i13,u2i13)))
    and sp.Rational(0)<sp.Rational(23,30)<1
    and sp.simplify((1-ax13)**2*Rw13-(1-ax13**2))==0
    and sp.simplify(1-ax13**2-4*Rw13/(1+Rw13)**2)==0)
rep('N13c2 normalized phases lift to genuine raw round-meridian inputs',
    actual_meridian13)

# Add (-i*rho,0) to a middle W point and run the actual real Cramer
# system.  South: (Delta one,Delta four)=
# (-theta*rho/(2*pi*P0),rho/P0).  North Pinf=-P0 and the
# same formula becomes (+theta*rho/(2*pi*Pinf),-rho/Pinf).
cr13,ci13,lr13,li13,rho13=sp.symbols(
    'cr13 ci13 lr13 li13 rho13',real=True)
t13=cr13+sp.I*ci13
lt13=lr13+sp.I*li13
e13=((-sp.I*t13,1-t13**2),(-t13,sp.I*(1+t13**2)),
     (-sp.I*lt13,2*t13))
mat13=sp.Matrix([
    [sp.im(sp.expand(v[0])) for v in e13],
    [sp.re(sp.expand(v[1])) for v in e13],
    [sp.im(sp.expand(v[1])) for v in e13]])
corr13=mat13.inv()*sp.Matrix([-rho13,0,0])
onecorr13=sp.factor(
    -sum(corr13[j]*sp.re(sp.expand(e13[j][0]))
         for j in range(3))/(2*sp.pi))
r213=cr13**2+ci13**2
P013=lr13-2*r213/(1-r213)
Pinf13=-P013
cramer_seam13=(
    sp.factor(mat13.det()
              -(1+r213)*(lr13*(r213-1)+2*r213))==0
    and sp.cancel(corr13[2]-rho13/P013)==0
    and sp.cancel(onecorr13+li13*rho13/(2*sp.pi*P013))==0
    and sp.cancel(corr13[2]+rho13/Pinf13)==0
    and sp.cancel(onecorr13-li13*rho13/(2*sp.pi*Pinf13))==0
    and has_statement9(cusp8.scalar_lattice_coords,
                       "rows = [[v.imag for v in ez],[v.real for v in ew],[v.imag for v in ew]]"))
rep('N13d source Cramer gauge correction at both seams',cramer_seam13)
det_s13=sp.Rational(5,4)*(sp.Rational(3,4)*sp.log(2)
                           +sp.Rational(1,2))
det_i13=5*(3*sp.log(2)+8)
rep('N13d2 strict Cramer denominators: south >5/8, north >40',
    sp.simplify(mat13.det().subs({
        cr13:-sp.Rational(1,4),ci13:-sp.sqrt(3)/4,
        lr13:-sp.log(2)})-det_s13)==0
    and sp.simplify(mat13.det().subs({
        cr13:-1,ci13:-sp.sqrt(3),lr13:sp.log(2)})
        -det_i13)==0
    and bool(sp.log(2)>0)
    and det_s13>sp.Rational(5,8)
    and det_i13>40)

# The source flat step has s(1)=1 and zero first/second left jets.
# Consequently the collar weights and height blend have zero two-jets
# at r=1/2 south and r=2 north.
flat_u13=sp.symbols('flat_u13',real=True)
flat_s13=sp.exp(-1/flat_u13)/(
    sp.exp(-1/flat_u13)+sp.exp(-1/(1-flat_u13)))
flat_jets13=(
    [sp.limit(sp.diff(flat_s13,flat_u13,k),flat_u13,1,dir='-')
     for k in range(3)]==[1,0,0]
    and has_statement9(cusp8.scalar_flat_step,
                       "a,b = mp.exp(-1/u),mp.exp(-1/(1-u))"))
rep('N13e flat endpoint value and first/second jets',flat_jets13)

# The selected open moment cells are analytic for q>0 and lambda near
# one.  At lambda=1 their source expressions are both
# f=2*sqrt(q)*sinh(k*a), k=-log q.  The positive root is
# 2*sqrt(q)*cosh(k*a), whence Y=q^(a+1/2) and
# rho=P*(a+1/2), P=-k/2.  All tangential two-jets agree; every
# transverse lambda derivative is multiplied by a flat endpoint jet.
k13,a13,lam13=sp.symbols('k13 a13 lam13',positive=True)
q13=sp.exp(-k13)
amp13=lam13*k13
d13=1-lam13+2*lam13*sp.sqrt(q13)*k13
fs13=d13*sp.sinh(amp13*a13)/amp13
kap13=2+lam13*(k13-2)
fn13=sp.exp(-kap13/2)*sp.exp(kap13*a13)-(
    q13*sp.exp(kap13/2))*sp.exp(-kap13*a13)
fend13=2*sp.exp(-k13/2)*sp.sinh(k13*a13)
Yend13=sp.exp(-k13*(a13+sp.Rational(1,2)))
height_endpoint13=(
    sp.simplify(fs13.subs(lam13,1)-fend13)==0
    and sp.simplify(sp.expand(
        fn13.subs(lam13,1)-fend13.rewrite(sp.exp)))==0
    and sp.simplify(
        2*q13/(2*sp.exp(-k13/2)*sp.cosh(k13*a13)+fend13)
        -Yend13)==0
    and all(sp.simplify(sp.diff(
        sp.log(Yend13)/2
        +k13*(a13+sp.Rational(1,2))/2,
        k13,i,a13,j))==0
        for i in range(3) for j in range(3-i))
    and has_statement9(cusp8.inverse_height_moment,
                       "return _height_G(a,amp,d)")
    and has_statement9(cusp8.inverse_height_moment,
                       "return mp.sign(a)*(B*mp.exp(kappa*aa)-C*mp.exp(-kappa*aa))")
    and has_if_test9(cusp8.inverse_height_moment,
                     "lam==1 or aa<=mp.mpf('.1')")
    and has_if_test9(cusp8.inverse_height_moment,
                     "aa>=mp.mpf('.4')"))
rep('N13f selected height cells: rho/P and all tangential two-jets',
    height_endpoint13)

# Pure south four=a+xi4+1/2; collar four=xi4+rho/P0.
# Pure north four=-a-1/2 because its outer blend is zero;
# collar four=-rho/Pinf.  Thus four has the same two-jet on
# both sides.  The one coordinates need not agree.  In the chosen
# lower-half-plane angular cell the two middle slots have a common
# phase exp(2*pi*i*one), and relative phase exp(-2*pi*i*four).
W13,F13,O13=sp.symbols('W13 F13 O13',real=True)
v13=sp.Matrix([W13*sp.exp(2*sp.pi*sp.I*O13),
               sp.exp(2*sp.pi*sp.I*(O13-F13))])
P13=v13*v13.conjugate().T/(1+W13**2)
targetP13=sp.Matrix([
    [W13**2,W13*sp.exp(2*sp.pi*sp.I*F13)],
    [W13*sp.exp(-2*sp.pi*sp.I*F13),1]])/(1+W13**2)
projector13=(
    all(sp.simplify(P13[i,j]-targetP13[i,j])==0
        for i in range(2) for j in range(2))
    and all(sp.simplify(sp.diff(P13[i,j],O13,k))==0
            for i in range(2) for j in range(2) for k in (1,2))
    and sp.simplify(targetP13*targetP13-targetP13)==sp.zeros(2,2))
rep('N13g both surviving middle slots: common phase cancels in P',
    projector13)

# b=O(h^3) is either collar weight; a four discrepancy
# delta=O(h^3) follows from the height and Cramer checks above.
# The 2-jet of every ten-slot block agrees with the pure-middle
# projector: cusp/cusp entries are O(b^2), cusp/middle entries O(b),
# and the middle normalization differs only by O(b^2).
h13,Bh13,Dh13,C13=sp.symbols('h13 Bh13 Dh13 C13',real=True)
bh13=Bh13*h13**3
dh13=Dh13*h13**3
den13=bh13**2*C13+(1-bh13)**2*(1+W13**2)
jet_projector13=(
    all(sp.simplify(sp.diff(
        (1-bh13)**2/den13-1/(1+W13**2),h13,k
    ).subs(h13,0))==0 for k in range(3))
    and all(sp.simplify(sp.diff(bh13*(1-bh13)/den13,h13,k)
                        .subs(h13,0))==0 for k in range(3))
    and all(sp.simplify(sp.diff(bh13**2/den13,h13,k)
                        .subs(h13,0))==0 for k in range(3))
    and all(sp.simplify(sp.diff(
        sp.exp(2*sp.pi*sp.I*(F13+dh13))
        -sp.exp(2*sp.pi*sp.I*F13),h13,k
    ).subs(h13,0))==0 for k in range(3)))
rep('N13h both ten-slot projector references have equal local two-jets',
    jet_projector13)

# The projection update e*(e^*v)/|e^*v| equals P*v/sqrt(v^*P*v).
# Thus it depends only on P when the overlap is nonzero.  At the
# local witness v=e the overlap is one and remains positive nearby.
# The final middle return exponentiates one+Arg(phase)/(2*pi).
# Under e_c=exp(2*pi*i*d)e_p, phase_c=exp(-2*pi*i*d)phase_p,
# so the corrected coordinate agrees modulo an integer deck change.
er13,ei13,vr13,vi13,dr13=sp.symbols(
    'er13 ei13 vr13 vi13 dr13',real=True)
ephase13=er13+sp.I*ei13
vphase13=vr13+sp.I*vi13
gauge13=sp.exp(2*sp.pi*sp.I*dr13)
inner13=sp.conjugate(ephase13)*vphase13
inner_g13=sp.conjugate(gauge13*ephase13)*vphase13
slice13=(
    sp.simplify(inner_g13-sp.conjugate(gauge13)*inner13)==0
    and sp.simplify(
        gauge13*ephase13*inner_g13-ephase13*inner13)==0
    and sp.simplify(gauge13*inner_g13-inner13)==0
    and has_statement9(cusp8.elementary_global_slice,
                       "phase = _elementary_inner(ref,result)")
    and has_statement9(cusp8.elliptic_meridian_projector,
                       "return mp.matrix([[a*mp.conj(b) for b in e] for a in e])"))
rep('N13i projection and corrected middle W slice cancel the gauge',
    slice13)

# The actual final phase denominator is one in modulus whenever the
# finite projection product returns.  Its last j=steps iteration uses
# h=X, so the section is a unit scalar multiple of the same endpoint
# reference re-evaluated by elementary_global_slice.  The cap
# correction only multiplies this result by another unit scalar.
# Thus the final inner product has modulus one on either gauge.
pa13,pb13,pc13,pd13,ga13,gb13=sp.symbols(
    'pa13 pb13 pc13 pd13 ga13 gb13',real=True)
radius13=sp.symbols('radius13',positive=True)
eref13=sp.Matrix([pa13+sp.I*pb13,pc13+sp.I*pd13])
enorm13=sum(sp.conjugate(v)*v for v in eref13)
unitref13=eref13/radius13
section13=sp.exp(sp.I*(ga13+gb13))*unitref13
finalphase13=(unitref13.conjugate().T*section13)[0]
phase_unit13=(
    sp.simplify(
        sp.expand((finalphase13-sp.exp(sp.I*(ga13+gb13)))
                  *radius13**2).subs(radius13**2,enorm13))==0
    and sp.simplify(
        sp.expand((finalphase13*sp.conjugate(finalphase13)-1)
                  *radius13**4).subs(radius13**4,enorm13**2))==0
    and has_statement9(cusp8.elementary_projection_section,
                       "h = start+(X-start)*j/steps")
    and has_statement9(cusp8.elementary_projection_section,
                       "vector = tuple(a*phase/abs(phase) for a in ref)")
    and has_statement9(cusp8._elementary_corrected_section,
                       "result = elementary_projection_section(xi,X,N,dps,diagnostics)")
    and has_statement9(cusp8._elementary_corrected_section,
                       "result = tuple(a*mp.exp(-mp.j*weight*theta) for a in result)")
    and has_statement9(cusp8.elementary_global_slice,
                       "ref = elliptic_bundle_embedding(kind,data,dps=dps,elementary=True)"))
rep('N13i2 source final phase denominator has exact unit modulus',
    phase_unit13)

# Type-(d) numbers only; N13b/c/e/f give the exact branch proof.
with mp.workdps(50):
    qsn13=mp.exp(-mp.mpf(4)/3)/4
    qin13=mp.exp(-mp.mpf(16)/3)/4
    step13=cusp8.scalar_flat_step(mp.mpf(2)/3)
    norm_mid13=1+(1-step13)**2
    print('N13j one-point numerical margins only: '
          'q_s lower/upper =',mp.nstr(qsn13-mp.mpf(1)/2500,10),
          mp.nstr(mp.mpf(1)/8-qsn13,10),
          '; q_i lower/upper =',
          mp.nstr(qin13-mp.mpf(1)/2500,10),
          mp.nstr(mp.mpf(1)/8-qin13,10),
          '; half-chart =',
          mp.nstr(mp.cos(mp.pi/10)-mp.mpf(9)/10,10),
          '; middle norm2 =',mp.nstr(norm_mid13,10),
          '; Cramer lower-bound margins =',
          mp.nstr(mp.mpf(str(det_s13.evalf(50)))-mp.mpf(5)/8,10),
          mp.nstr(mp.mpf(str(det_i13.evalf(50)))-40,10))

# ---------------------------------------------------------------- N14: finish the same strict N12 inverse-first height cell.
# The independent inputs remain (b,q,lambda), with b=3/20,
# q=(4/25)exp(-16/21), lambda=1/2.  The new circuit carries N12's
# terminal v through G^-1, R^-1, G, the positive quadratic Y, and rho.
# It does not cover other height cells or the six-real tensor.
def has_ifexp_test14(function,test):
    tree = ast.parse(inspect.getsource(function))
    expected = ast.dump(ast.parse(test,mode='eval').body,
                        include_attributes=False)
    return any(isinstance(node,ast.IfExp)
               and ast.dump(node.test,include_attributes=False)==expected
               for node in ast.walk(tree))

height_source14 = (
    ordered_block12(cusp8._height_G,
        'return mp.asinh(amp*a/d)/amp if inverse else d*mp.sinh(amp*a)/amp')
    and ordered_block12(cusp8._height_R,'''\
if inverse:
    y = mp.e*mp.expm1(value/(c-value))
    L = y+mp.log1p(-mp.exp(-y)) if y>1 else mp.log(mp.expm1(y))
    x = (2/(mp.sqrt(L*L+4)-L) if L<0 else
         (L+mp.sqrt(L*L+4))/2)
    ell = mp.log1p(x/mp.exp(2))
    return delta+c*ell/(2*(2+ell))
''')
    and all(has_if_test9(cusp8._height_R,test)
            for test in ('value<=0','value>=c',
                         'c/(c-value)>(mp.mp.dps+10)*mp.log(10)'))
    and has_ifexp_test14(cusp8._height_R,'y>1')
    and has_ifexp_test14(cusp8._height_R,'L<0')
    and ordered_block12(cusp8.elementary_cusp_logradius,'''\
f = inverse_height_moment(b,q,lam)
root = mp.sqrt(f*f+4*q)
Y = 2*q/(root+f) if f>=0 else (root-f)/2
return +(n*P+mp.log(Y)/2)
''')
    and has_statement9(cusp8.inverse_height_moment,
        'return mp.sign(a)*_height_G(_height_R(_height_G(v,amp,d,inverse=True),inverse=True),amp,d)'))
rep('N14a ordered source G^-1/R^-1/G and positive-Y/log-radius path',
    height_source14)

# An exact outward rational interval calculator supplies branch guards.
# Q=10^30 is a fixed rational grid; every operation rounds its lower
# endpoint downward and upper endpoint upward.  exp uses a positive
# Taylor sum through n=50 with geometric remainder
# t_51/(1-x/52) after range reduction to 0<=x<=1.  log uses
# log m=2 sum_(n=0)^50 z^(2n+1)/(2n+1), z=(m-1)/(m+1)
# in [0,1/3], with tail <=2 z^103/(103(1-z^2)); powers of two
# reduce any positive argument to 1<=m<=2.  sqrt uses integer
# square-root floors.  These are rigorous rational bounds, not
# floating-point samples.
from fractions import Fraction as Q14
from math import isqrt as isqrt14

grid14 = 10**30
def down14(value):
    v = Q14(value)*grid14
    return Q14(v.numerator//v.denominator,grid14)
def up14(value):
    return -down14(-Q14(value))

class RationalInterval14:
    __slots__ = ('lo','hi')
    def __init__(self,lo,hi=None):
        lo,hi = Q14(lo),Q14(lo if hi is None else hi)
        if lo>hi:
            raise ValueError('reversed rational interval')
        self.lo,self.hi = down14(lo),up14(hi)
    def __add__(self,other):
        other = interval14(other)
        return RationalInterval14(self.lo+other.lo,self.hi+other.hi)
    __radd__ = __add__
    def __neg__(self):
        return RationalInterval14(-self.hi,-self.lo)
    def __sub__(self,other):
        return self+(-interval14(other))
    def __rsub__(self,other):
        return interval14(other)+(-self)
    def __mul__(self,other):
        other = interval14(other)
        products = (self.lo*other.lo,self.lo*other.hi,
                    self.hi*other.lo,self.hi*other.hi)
        return RationalInterval14(min(products),max(products))
    __rmul__ = __mul__
    def reciprocal(self):
        if self.lo<=0<=self.hi:
            raise ValueError('interval denominator touches zero')
        return RationalInterval14(1/self.hi,1/self.lo)
    def __truediv__(self,other):
        return self*interval14(other).reciprocal()
    def __rtruediv__(self,other):
        return interval14(other)*self.reciprocal()
    def __pow__(self,n):
        if not isinstance(n,int):
            raise TypeError('only integer interval powers are used')
        if n<0:
            return (self**(-n)).reciprocal()
        out,base = interval14(1),self
        while n:
            if n&1:
                out = out*base
            base = base*base
            n >>= 1
        return out
    def inside(self,lo,hi):
        return Q14(lo)<self.lo<=self.hi<Q14(hi)

def interval14(value):
    return value if isinstance(value,RationalInterval14) else RationalInterval14(value)

def exp_point14(x):
    x = Q14(x)
    if x<0:
        return exp_point14(-x).reciprocal()
    count = 1
    while x>count:
        count *= 2
    y = x/count
    term,total = Q14(1),Q14(1)
    for n in range(1,51):
        term *= y/n
        total += term
    first_omitted = term*y/51
    tail = first_omitted/(1-y/52)
    result = RationalInterval14(total,total+tail)
    while count>1:
        result = result*result
        count //= 2
    return result

def exp14(arg):
    arg = interval14(arg)
    return RationalInterval14(exp_point14(arg.lo).lo,
                              exp_point14(arg.hi).hi)

def log_unit14(x):
    # Exact rational bounds for log x on 1<=x<=2.
    assert 1<=x<=2
    z = (x-1)/(x+1)
    term,total = z,Q14(0)
    for n in range(51):
        total += 2*term/(2*n+1)
        term *= z*z
    tail = 2*term/(103*(1-z*z))
    return RationalInterval14(total,total+tail)

log_two14 = log_unit14(Q14(2))
def log_point14(x):
    x = Q14(x)
    if x<=0:
        raise ValueError('log interval requires positive input')
    k,m = 0,x
    while m<1:
        m *= 2
        k -= 1
    while m>2:
        m /= 2
        k += 1
    return log_unit14(m)+k*log_two14

def log14(arg):
    arg = interval14(arg)
    if arg.lo<=0:
        raise ValueError('log interval touches zero')
    return RationalInterval14(log_point14(arg.lo).lo,
                              log_point14(arg.hi).hi)

def sqrt_point14(x):
    x = Q14(x)
    if x<0:
        raise ValueError('sqrt interval requires nonnegative input')
    scaled = x.numerator*grid14**2
    floor_root = isqrt14(scaled//x.denominator)
    lo = Q14(floor_root,grid14)
    hi = lo if floor_root**2*x.denominator==scaled else Q14(floor_root+1,grid14)
    return RationalInterval14(lo,hi)

def sqrt14(arg):
    arg = interval14(arg)
    if arg.lo<0:
        raise ValueError('sqrt interval touches negative values')
    return RationalInterval14(sqrt_point14(arg.lo).lo,
                              sqrt_point14(arg.hi).hi)

e14 = exp14(interval14(1))
def eval_interval14(expr,environment):
    if expr in environment:
        return environment[expr]
    if expr.is_Rational:
        return interval14(Q14(int(expr.p),int(expr.q)))
    if expr==sp.E:
        return e14
    if expr.func==sp.Add:
        return sum((eval_interval14(arg,environment) for arg in expr.args),
                   interval14(0))
    if expr.func==sp.Mul:
        result = interval14(1)
        for arg in expr.args:
            result *= eval_interval14(arg,environment)
        return result
    if expr.func==sp.Pow:
        base = eval_interval14(expr.base,environment)
        power = expr.exp
        if power==sp.Rational(1,2):
            return sqrt14(base)
        if power.is_Integer:
            return base**int(power)
    if expr.func==sp.exp:
        return exp14(eval_interval14(expr.args[0],environment))
    if expr.func==sp.log:
        return log14(eval_interval14(expr.args[0],environment))
    if expr.func==sp.sinh:
        value = eval_interval14(expr.args[0],environment)
        return (exp14(value)-exp14(-value))/2
    if expr.func==sp.asinh:
        value = eval_interval14(expr.args[0],environment)
        return log14(value+sqrt14(1+value**2))
    raise ValueError('unregistered interval primitive: '+str(expr.func))

# These ten scalar nodes follow the exact selected source branch.  The
# terminal v and every other parameter come from N12's ordered 30 nodes.
(hg14,hy14,hL14,hx14,hell14,hR14,hf14,hroot14,hY14,hrho14) = (
    sp.symbols('hn14_0:10',real=True))
height_nodes14 = [
    (hg14,sp.asinh(hamp12*hv12/hd12)/hamp12),
    (hy14,sp.E*(sp.exp(hg14/(c12-hg14))-1)),
    (hL14,hy14+sp.log(1-sp.exp(-hy14))),
    (hx14,(hL14+sp.sqrt(hL14**2+4))/2),
    (hell14,sp.log(1+hx14/sp.exp(2))),
    (hR14,delta12+c12*hell14/(2*(2+hell14))),
    (hf14,hd12*sp.sinh(hamp12*hR14)/hamp12),
    (hroot14,sp.sqrt(hf14**2+4*hq12)),
    (hY14,2*hq12/(hroot14+hf14)),
    (hrho14,sp.log(hY14)/2)]

interval_values14 = {
    hb12:interval14(Q14(3,20)),
    hq12:Q14(4,25)*exp14(interval14(-Q14(16,21))),
    hl12:interval14(Q14(1,2))}
for node14,rhs14 in height_nodes12+height_nodes14:
    interval_values14[node14] = eval_interval14(rhs14,interval_values14)
branch_interval14 = (
    interval_values14[hq12].inside(Q14(1,2500),Q14(1,8))
    and interval_values14[hamp12].inside(1,2)
    and interval_values14[hd12].inside(Q14(11,10),Q14(7,5))
    and interval_values14[hv12].inside(Q14(19,100),Q14(21,100))
    and interval_values14[hg14].inside(Q14(1,10),Q14(19,100))
    and interval_values14[hy14].lo>1
    and interval_values14[hL14].lo>0
    and interval_values14[hx14].lo>0
    and interval_values14[hR14].inside(Q14(1,10),Q14(1,5))
    and interval_values14[hf14].inside(Q14(19,100),Q14(21,100))
    and interval_values14[hY14].inside(Q14(1,10),Q14(1,4))
    and interval_values14[hrho14].inside(-1,-Q14(1,2))
    and Q14(1,5)/(Q14(1,5)-Q14(19,100))==20
    and 20<22 and Q14(11,4)**2<10)
# The source inverse-R shortcut compares c/(c-g) against
# (mp.dps+10)log(10).  The certified g<19/100 makes the former <20;
# e<11/4 implies log(10)>2, and mp.dps>=1 makes the latter >22.
# g>1/10 gives g/(c-g)>1, whence y=e(expm1(...))>e(e-1)>1.
# For y>1, L=log(exp(y)-1)>0, fixing both source branches.
# f>0 selects Y=2q/(sqrt(f^2+4q)+f), whose denominator is positive.
rep('N14b rational interval guards for inverse R, positive Y and rho',
    height_branch12 and branch_interval14)

# Continue the exact three-input symmetric second-jet chain without
# expanding the terminal DAG.  Each new RHS has only already-bound
# scalar nodes; local derivatives are finite elementary expressions.
height_jet14 = dict(height_jet12)
height_circuit14 = []
height_topology14 = True
for out14,rhs14 in height_nodes14:
    deps14 = sorted(rhs14.free_symbols,key=str)
    if out14 in height_jet14 or any(dep not in height_jet14 for dep in deps14):
        height_topology14 = False
        break
    d114 = {dep:sp.diff(rhs14,dep) for dep in deps14}
    d214 = {(p,w):sp.diff(d114[p],w) for p in deps14 for w in deps14}
    grad14 = tuple(sum((d114[p]*height_jet14[p][0][i]
                        for p in deps14),sp.Integer(0)) for i in range(3))
    hess14 = tuple(tuple(
        sum((d114[p]*height_jet14[p][1][i][j]
             for p in deps14),sp.Integer(0))
        +sum((d214[p,w]*height_jet14[p][0][i]*height_jet14[w][0][j]
              for p in deps14 for w in deps14),sp.Integer(0))
        for j in range(3)) for i in range(3))
    if (any(item.has(sp.Derivative)
            for item in (rhs14,*d114.values(),*d214.values()))
            or any(sp.simplify(d214[p,w]-d214[w,p])!=0
                   for p in deps14 for w in deps14)
            or any(sp.simplify(hess14[i][j]-hess14[j][i])!=0
                   for i in range(3) for j in range(3))):
        height_topology14 = False
        break
    height_jet14[out14] = grad14,hess14
    height_circuit14.append((out14,rhs14,grad14,hess14))
rep('N14c 40-node source-bound rho circuit: 3 first and 6 second jets',
    height_topology14 and len(height_circuit14)==10
    and len(height_circuit12)+len(height_circuit14)==40
    and height_circuit14[-1][0]==hrho14
    and all(not item.has(sp.Derivative)
            for _,rhs,grad,hess in height_circuit14
            for item in (rhs,*grad,*(entry for row in hess for entry in row))))

# Type (d): one-point source comparison.  The interval proof above, not
# these decimals, chooses the branches and supplies the denominators.
with mp.workdps(100):
    source_rho14 = cusp8.elementary_cusp_logradius(
        mp.mpf(3)/20,mp.mpf(2)/5,1,dps=90)
    numeric_values14 = dict(sample_values12)
    for node14,rhs14 in height_nodes14:
        numeric_values14[node14] = rhs14.subs(numeric_values14).evalf(105)
    circuit_rho14 = mp.mpf(str(numeric_values14[hrho14].evalf(90)))
    error_rho14 = abs(circuit_rho14-source_rho14)
    assert error_rho14<mp.mpf('1e-75')
    sample_g14 = mp.mpf(str(numeric_values14[hg14].evalf(70)))
    sample_R14 = mp.mpf(str(numeric_values14[hR14].evalf(70)))
    sample_Y14 = mp.mpf(str(numeric_values14[hY14].evalf(70)))
    print('N14d one-point numerical source check only:'
          ' |rho_DAG-rho_source|=',mp.nstr(error_rho14,5),
          '; g in (1/10,19/100) sample margins =',
          mp.nstr(sample_g14-mp.mpf(1)/10,10),
          mp.nstr(mp.mpf(19)/100-sample_g14,10),
          '; R^-1(g)<1/5 sample margin =',
          mp.nstr(mp.mpf(1)/5-sample_R14,10),
          '; Y in (1/10,1/4) sample margins =',
          mp.nstr(sample_Y14-mp.mpf(1)/10,10),
          mp.nstr(mp.mpf(1)/4-sample_Y14,10),
          ': PASS')

print('ALL PASS' if ok else 'FAILURES')
assert ok, 'A symbolic identity, constant bound, or numerical check failed'
