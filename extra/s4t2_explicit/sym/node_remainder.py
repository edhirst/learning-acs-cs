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

print('ALL PASS' if ok else 'FAILURES')
assert ok, 'A symbolic identity, constant bound, or numerical check failed'
