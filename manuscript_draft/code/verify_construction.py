"""Symbolic checks for the complex structure on S^4 x S^1 x S^1 (see the paper
"An Explicit Complex Structure on S^4 x T^2", Appendix A).  Pure sympy algebra:
no sampling, no optimisation.

Family over t in C^*:  X_t = (C^*_zeta x C_w) / <g2, g3, g4>,
    g2(zeta, w) = (e^t      zeta, w + 1 - t^2)
    g3(zeta, w) = (e^{-i t} zeta, w + i(1 + t^2))
    g4(zeta, w) = (t        zeta, w + 2 t)
In additive coordinates (z, w), zeta = e^{i z}, the period lattice is
    e1 = (2 pi, 0), e2 = (-i t, 1 - t^2), e3 = (-t, i(1 + t^2)), e4 = (-i log t, 2 t).

Checks
 1. det_R(e1, e2, e3, e4) = 2 pi (1 + r^2) (2 r^2 - (1 - r^2) ln r) > 0 for all t != 0.
 2. Monodromy: t -> e^{2 pi i} t sends e4 -> e4 + e1, fixes e1, e2, e3.
 3. Coordinate change near t = oo, (zeta', w') = (zeta e^{w/t}, w/t^2), s = 1/t,
    turns the generators into g2' = (e^s, s^2 - 1), g3' = (e^{i s}, i(1 + s^2)),
    g4^{-1} = (e^{-2} s, -2 s): the same type-II local model as at t = 0.
 4. The twisting section sigma(s) = (pi s (s - 1), 2 t s), t = e^{2 pi i s},
    satisfies sigma(s + 1) - sigma(s) = e4(s): it is a holomorphic section of the
    family over a punctured disc whose Chern class is [e4].
 5. The T^2-action (translations by c2 e2 + c3 e3, c real) extends across t = oo:
    the zeta'-multiplier tends to 1 and the shift to the point -c2 + i c3 of E_oo.
 6. The relative 2-form dz ^ dw equals s^{-2} dz'' ^ dw' near t = oo, and the
    3-form dz ^ dw ^ dt equals -s^{-4} dz'' ^ dw' ^ ds, so
    omega_{Y/P^1} = f^*O(-2) and K_Y = f^*O(-4).
 7. With multipliers e^{lam t}, e^{-i lam t} and g4 = (e^{i beta} t, 2t), the
    determinant is 2 pi (1 + r^2)(2 lam r^2 + (Im beta - ln r)(1 - r^2)).
 8. Toric surface X_Sigma (rays (n,1)): chart coordinates x_n = zeta^{-1} t^{n+1},
    y_n = zeta t^{-n} satisfy t = x_n y_n, zeta = x_n^n y_n^{n+1}; the shift
    psi(zeta,t) = (t zeta, t) satisfies x_{n+1} o psi = x_n, y_{n+1} o psi = y_n;
    on chart overlaps x_{n+1} = t x_n, y_{n+1} = 1/x_n; and
    (dzeta/zeta) ^ dt = -dx_n ^ dy_n (so the log 3-form is a generator near F_0).
 9. Abelian parts: omega(e2bar, e3bar) = 1 - |t|^4 near t = 0 and |s|^4 - 1 near
    s = 0 (opposite orientations), where omega(u, v) = Im(conj(u) v).
10. Local determinant identity used in the filling lemma (checked for arbitrary
    symbolic data): for the R-linear phi with phi(e2bar) = Re alpha2,
    phi(e3bar) = Re alpha3, det_R = -2 pi P omega(e2bar, e3bar) with
    P = Re log(nu t) - phi(e4bar).
11. Exact heights (Remark "the heights P exactly"): from the filling data,
    P_0(r) = ln r - 2 r^2/(1 - r^2) at t = 0 and P_oo(r) = ln r - 2/(1 - r^2) at
    s = 0.  Hence P/ln r = 1 + 2r^2/((1-r^2)|ln r|), resp. 1 + 2/((1-r^2)|ln r|),
    both increasing on (0,1); at r = 1e-2 they are < 1.000044 and < 1.434339, so
    eps = 1e-2 satisfies P/log|t| in [1/2, 3/2] (Prop. 4.1, Step 0); eps = 1/16
    fails at oo.
12. Periods of dz ^ dw (Lemma "periods"): p_jk = z_j w_k - z_k w_j, the six
    closed forms, and no nonzero integral relation (coefficients of 1, t, t^2 and
    of log t, t log t, t^2 log t); X_1 is the product lattice
    (2 pi Z + i Z) x (2Z + 2iZ) after (z, w) -> (z - i w/2, w), e4 -> e4 - e2.
13. Kodaira-Spencer: d/dt (p13/p12) = 4 i t/(1 - t^2)^2, nonzero on C^* \ {+-1}.
14. Admissible region of the deformation (Remark "deformations"): the critical
    point r = sqrt(lam+1) - sqrt(lam) of g = ln r - 2 lam r^2/(1-r^2), the value
    max g = lam - sqrt(lam(lam+1)) - arsinh(sqrt lam), and the counterexample
    lam = 1, beta = -2i, D(1/2) = -1 + (3/4) ln 2 < 0.
15. Toric dual bases: A^T m_i^{(n+1)} = m_i^{(n)} for A = [[1,1],[0,1]].
Positivity in check 1 is the only inequality in the construction and holds
identically in t; everything else is an exact identity.
"""
import sympy as sp

a, b = sp.symbols('a b', real=True)
t = a + sp.I * b
r2 = a**2 + b**2


def realvec(z1, z2):
    return [sp.re(z1), sp.im(z1), sp.re(z2), sp.im(z2)]


# -- 1. lattice determinant ------------------------------------------------
# log t = ln r + i arg t.  e1 has only a Re z component, so expanding along it
# removes the Re z row: only Im(-i log t) = -ln r enters the determinant.
lnr, arg = sp.symbols('lnr arg', real=True)
logt = lnr + sp.I * arg
e1 = (2 * sp.pi, 0)
e2 = (-sp.I * t, 1 - t**2)
e3 = (-t, sp.I * (1 + t**2))
e4 = (-sp.I * logt, 2 * t)
M = sp.Matrix([realvec(*e) for e in (e1, e2, e3, e4)]).T
D = sp.simplify(sp.expand(M.det()))
target = 2 * sp.pi * (1 + r2) * (2 * r2 - (1 - r2) * lnr)
assert sp.simplify(sp.expand(D - target)) == 0, (D, target)
print('1. det_R(e1..e4) =', sp.factor(D))
print('   = 2*pi*(1+r^2)*(2 r^2 - (1-r^2) ln r); arg(t) drops out:',
      sp.diff(D, arg) == 0)

# positivity of g(r) = 2 r^2 - (1 - r^2) ln r on r > 0: (1 - r^2) and ln r
# have opposite signs on (0,1) and on (1,oo) and both vanish only at r = 1,
# so (1 - r^2) ln r <= 0 and g(r) >= 2 r^2 > 0.
r = sp.symbols('r', positive=True)
h = (1 - r**2) * sp.log(r)
assert h.subs(r, 1) == 0
assert sp.sign(h.subs(r, sp.Rational(1, 2))) == -1 and sp.sign(h.subs(r, 2)) == -1
print('   (1-r^2) ln r <= 0 for all r>0 (sign of each factor flips at r=1),'
      ' so g(r) >= 2 r^2 > 0.')

# -- 2. monodromy ----------------------------------------------------------
e4_after = (-sp.I * (logt + 2 * sp.pi * sp.I), 2 * t)
diff = [sp.simplify(e4_after[k] - e4[k] - e1[k]) for k in range(2)]
assert diff == [0, 0]
print('2. monodromy about t=0: e4 -> e4 + e1 (e1, e2, e3 single valued)')

# -- 3. coordinates near infinity -----------------------------------------
tt = sp.symbols('t', nonzero=True)
s = 1 / tt


def transform(mult, shift):
    """Generator (zeta,w) -> (mult*zeta, w+shift) in (zeta', w') coords."""
    # zeta' = zeta e^{w/t}: multiplier becomes mult * e^{shift/t}
    return sp.simplify(mult * sp.exp(shift / tt)), sp.simplify(shift / tt**2)


g2p = transform(sp.exp(tt), 1 - tt**2)
g3p = transform(sp.exp(-sp.I * tt), sp.I * (1 + tt**2))
g4p = transform(tt, 2 * tt)
assert sp.simplify(g2p[0] - sp.exp(s)) == 0 and sp.simplify(g2p[1] - (s**2 - 1)) == 0
assert sp.simplify(g3p[0] - sp.exp(sp.I * s)) == 0 and sp.simplify(g3p[1] - sp.I * (1 + s**2)) == 0
assert sp.simplify(g4p[0] - sp.exp(2) / s) == 0 and sp.simplify(g4p[1] - 2 * s) == 0
print('3. near t=oo (s=1/t): g2 ->', g2p, ' g3 ->', g3p, ' g4 ->', g4p)
print('   so g4^{-1} = (e^{-2} s zeta\', w\' - 2 s): contraction, E_oo = C/<-1, i>')

# -- 4. twisting section ---------------------------------------------------
S = sp.symbols('S')
tS = sp.exp(2 * sp.pi * sp.I * S)
sigma = lambda x: (sp.pi * x * (x - 1), 2 * sp.exp(2 * sp.pi * sp.I * x) * x)
e4S = (2 * sp.pi * S, 2 * tS)                   # -i log t = 2 pi S
d = [sp.simplify(sp.expand(sigma(S + 1)[k] - sigma(S)[k] - e4S[k])) for k in range(2)]
assert d == [0, 0], d
print('4. sigma(S+1) - sigma(S) = e4(S): twist class [e4] in the coinvariants')

# -- 5. T^2 action at infinity ---------------------------------------------
c2, c3 = sp.symbols('c2 c3', real=True)
mult = sp.exp((c2 - sp.I * c3) * tt)            # zeta multiplier of c2 e2 + c3 e3
shift = c2 * (1 - tt**2) + c3 * sp.I * (1 + tt**2)
mp, sh = transform(mult, shift)
assert sp.simplify(mp - sp.exp((c2 + sp.I * c3) * s)) == 0
assert sp.simplify(sh - (c2 * (s**2 - 1) + sp.I * c3 * (1 + s**2))) == 0
print('5. T^2 translation at oo: multiplier exp((c2+i c3) s) -> 1, shift ->',
      '-c2 + i c3 (a point of E_oo): extends over the filling')

# -- 6. relative 2-form ------------------------------------------------------
# z'' = z - i w/t (so zeta' = e^{i z''}), w' = w/t^2: Jacobian on a fibre
J = sp.Matrix([[1, -sp.I / tt], [0, 1 / tt**2]])
assert sp.simplify(J.det() - s**2) == 0
print('6. dz\'\' ^ dw\' = s^2 dz ^ dw: dz^dw has a double pole along F_oo,'
      ' so omega_{Y/P1} = f*O(-2), K_Y = f*O(-4)')
# 3-form: dz ^ dw ^ dt with z = z'' + i w/t, w = w'/s^2, t = 1/s.
zz, ww, S_ = sp.symbols("z2 w2 s", nonzero=True)
z_of = zz + sp.I * (ww / S_**2) * S_            # z = z'' + i w / t, w = w'/s^2
w_of = ww / S_**2
t_of = 1 / S_
Jac = sp.Matrix([[sp.diff(f_, v) for v in (zz, ww, S_)] for f_ in (z_of, w_of, t_of)])
assert sp.simplify(Jac.det() + S_**-4) == 0
print('   and dz^dw^dt = -s^{-4} dz\'\'^dw\'^ds: div = -4 F_oo')

# -- 7. general parameters ---------------------------------------------------
lam, bre, bim = sp.symbols('lam bre bim', real=True)
beta = bre + sp.I * bim
f2 = (-sp.I * lam * t, 1 - t**2)
f3 = (-lam * t, sp.I * (1 + t**2))
f4 = (beta - sp.I * logt, 2 * t)
Mg = sp.Matrix([realvec(*e) for e in (e1, f2, f3, f4)]).T
Dg = sp.expand(Mg.det())
targ = 2 * sp.pi * (1 + r2) * (2 * lam * r2 + (bim - lnr) * (1 - r2))
assert sp.simplify(Dg - sp.expand(targ)) == 0
print('7. general (lam, beta): det = 2 pi (1+r^2)(2 lam r^2 + (Im beta - ln r)(1-r^2))')
# -- 8. toric charts -----------------------------------------------------------
n = sp.symbols('n', integer=True)
Z, T = sp.symbols('zeta t', positive=True)   # monomial identities; positivity only for simplification
xn = lambda k: Z**(-1) * T**(k + 1)
yn = lambda k: Z * T**(-k)
assert sp.simplify(xn(n) * yn(n) - T) == 0
assert sp.simplify(xn(n)**n * yn(n)**(n + 1) - Z) == 0
psiZ = T * Z                                   # zeta-coordinate of psi(zeta, t)
assert sp.simplify(xn(n + 1).subs(Z, psiZ) - xn(n)) == 0
assert sp.simplify(yn(n + 1).subs(Z, psiZ) - yn(n)) == 0
assert sp.simplify(xn(n + 1) - T * xn(n)) == 0 and sp.simplify(yn(n + 1) - 1 / xn(n)) == 0
X_, Y_ = sp.symbols('x y', positive=True)
logz = n * sp.log(X_) + (n + 1) * sp.log(Y_)
Jt = sp.Matrix([[sp.diff(logz, X_), sp.diff(logz, Y_)],
                [sp.diff(X_ * Y_, X_), sp.diff(X_ * Y_, Y_)]])
assert sp.simplify(Jt.det() + 1) == 0
print('8. toric charts: t = x y, zeta = x^n y^(n+1), psi shifts charts,'
      ' overlaps x\' = t x, y\' = 1/x, dlog(zeta) ^ dt = -dx ^ dy')

# -- 9. orientation of the abelian parts --------------------------------------
om = lambda u, v: sp.im(sp.conjugate(u) * v)
w0 = sp.simplify(sp.expand(om(1 - t**2, sp.I * (1 + t**2))))
assert sp.simplify(w0 - (1 - r2**2)) == 0
sa, sb = sp.symbols('sa sb', real=True)
ss = sa + sp.I * sb
winf = sp.simplify(sp.expand(om(ss**2 - 1, sp.I * (1 + ss**2))))
assert sp.simplify(winf - ((sa**2 + sb**2)**2 - 1)) == 0
print('9. omega(e2,e3) = 1 - |t|^4 near 0 and |s|^4 - 1 near oo: opposite signs')

# -- 10. det = -2 pi P omega for arbitrary data ---------------------------------
syms = sp.symbols('A2r A2i A3r A3i u2r u2i u3r u3i u4r u4i lc', real=True)
A2r, A2i, A3r, A3i, u2r, u2i, u3r, u3i, u4r, u4i, lc = syms
al2, al3 = A2r + sp.I * A2i, A3r + sp.I * A3i
u2, u3, u4 = u2r + sp.I * u2i, u3r + sp.I * u3i, u4r + sp.I * u4i
argct = sp.symbols('argct', real=True)
vecs = [(2 * sp.pi, 0), (-sp.I * al2, u2), (-sp.I * al3, u3),
        (-sp.I * (lc + sp.I * argct), u4)]          # log(c t) = lc + i argct
Md = sp.Matrix([realvec(*v) for v in vecs]).T
p_, q_ = sp.symbols('p q', real=True)              # phi(x + i y) = p x + q y
sol = sp.solve([p_ * u2r + q_ * u2i - A2r, p_ * u3r + q_ * u3i - A3r], [p_, q_], dict=True)[0]
phi4 = (p_ * u4r + q_ * u4i).subs(sol)
P = lc - phi4
omg = u2r * u3i - u2i * u3r
assert sp.simplify(Md.det() - (-2 * sp.pi * P * omg)) == 0
print('10. det_R = -2 pi P(t) omega(e2bar, e3bar) for arbitrary data')
# -- 11. exact heights P at both cusps ---------------------------------------------
def height(al2, al3, u2, u3, u4, lognu_abs):
    """P = log|nu| + ln r - phi(u4), phi R-linear with phi(u_j) = Re(alpha_j)."""
    pp, qq = sp.symbols('pp qq', real=True)
    sol_ = sp.solve([pp * sp.re(u2) + qq * sp.im(u2) - sp.re(al2),
                     pp * sp.re(u3) + qq * sp.im(u3) - sp.re(al3)], [pp, qq], dict=True)[0]
    return lognu_abs + lnr - (pp * sp.re(u4) + qq * sp.im(u4)).subs(sol_)
P0 = height(t, -sp.I * t, 1 - t**2, sp.I * (1 + t**2), 2 * t, 0)
Pinf = height(t, sp.I * t, t**2 - 1, sp.I * (1 + t**2), -2 * t, -2)   # t plays s here
assert sp.simplify(P0 - (lnr - 2 * r2 / (1 - r2))) == 0
assert sp.simplify(Pinf - (lnr - 2 / (1 - r2))) == 0
ratio0 = lambda x: 1 + 2 * x**2 / ((1 - x**2) * sp.Abs(sp.log(x)))
ratioi = lambda x: 1 + 2 / ((1 - x**2) * sp.Abs(sp.log(x)))
# Both ratios increase on (0,1): their numerators 2x^2 and 2 are non-decreasing,
# and the common denominator (1 - x^2)|ln x| has derivative 2x ln x - (1 - x^2)/x < 0.
x = sp.symbols('x', positive=True)
den = (1 - x**2) * (-sp.log(x))
assert sp.simplify(sp.diff(den, x) - (2 * x * sp.log(x) - (1 - x**2) / x)) == 0
e2_ = sp.Rational(1, 100)
assert ratio0(e2_) < sp.Rational(1000044, 10**6) and ratioi(e2_) < sp.Rational(1434339, 10**6)
assert ratioi(sp.Rational(1, 16)) > sp.Rational(3, 2)
print('11. P_0 = ln r - 2r^2/(1-r^2), P_oo = ln r - 2/(1-r^2); at eps=1e-2: P/ln r <',
      '%.6f (t=0), %.6f (oo); eps=1/16 gives %.3f > 3/2 at oo'
      % (float(ratio0(e2_)), float(ratioi(e2_)), float(ratioi(sp.Rational(1, 16)))))

# -- 12. periods of dz^dw and the Neron-Severi argument ------------------------------
tc = sp.symbols('tc')
Lg = sp.symbols('Lg')                                  # stands for log t
E = {1: (2 * sp.pi, 0), 2: (-sp.I * tc, 1 - tc**2), 3: (-tc, sp.I * (1 + tc**2)),
     4: (-sp.I * Lg, 2 * tc)}
per = {(j, k): sp.expand(E[j][0] * E[k][1] - E[k][0] * E[j][1])
       for j in range(1, 5) for k in range(j + 1, 5)}
expected = {(1, 2): 2 * sp.pi * (1 - tc**2), (1, 3): 2 * sp.pi * sp.I * (1 + tc**2),
            (1, 4): 4 * sp.pi * tc, (2, 3): 2 * tc,
            (2, 4): sp.I * (1 - tc**2) * Lg - 2 * sp.I * tc**2,
            (3, 4): -(1 + tc**2) * Lg - 2 * tc**2}
assert all(sp.expand(per[k] - expected[k]) == 0 for k in per)
cs = sp.symbols('c12 c13 c14 c23 c24 c34', integer=True)
F = sp.expand(sum(c * per[k] for c, k in zip(cs, sorted(per))))
coeffs = sp.Poly(F, tc, Lg).as_dict()
# log t is not rational in t, so every coefficient of t^a (log t)^b must vanish;
# split each coefficient into real and imaginary parts and solve over the integers,
# using that 1 and pi are linearly independent over Q.
eqs = []
for c in coeffs.values():
    for part in (sp.re(c), sp.im(c)):
        part = sp.expand(part)
        eqs += [part.coeff(sp.pi, 1), part.coeff(sp.pi, 0)]
solns = sp.solve([q for q in eqs if q != 0], cs, dict=True)
assert solns == [dict(zip(cs, [0] * 6))], solns
lat1 = [(2 * sp.pi, 0), (-sp.I, 0), (-1, 2 * sp.I), (0, 2)]
chg = lambda v: (sp.expand(v[0] - sp.I * v[1] / 2), v[1])
new = [chg(v) for v in lat1]
new[3] = (sp.expand(new[3][0] - new[1][0]), new[3][1] - new[1][1])
assert new == [(2 * sp.pi, 0), (-sp.I, 0), (0, 2 * sp.I), (0, 2)], new
print('12. periods p_jk as in the paper; no nonzero integral relation;',
      'X_1 = C/(2pi Z + iZ) x C/(2Z + 2iZ)')

# -- 13. Kodaira-Spencer -----------------------------------------------------------------
dr = sp.simplify(sp.diff(per[(1, 3)] / per[(1, 2)], tc))
assert sp.simplify(dr - 4 * sp.I * tc / (1 - tc**2)**2) == 0
print('13. d/dt (p13/p12) =', dr, '(nonzero for t != 0, +-1)')

# -- 14. admissible region of the deformation -------------------------------------------
lp, rr_ = sp.symbols('lp rr', positive=True)
g = sp.log(rr_) - 2 * lp * rr_**2 / (1 - rr_**2)
rl = sp.sqrt(lp + 1) - sp.sqrt(lp)
assert sp.simplify(sp.diff(g, rr_).subs(rr_, rl)) == 0
assert sp.simplify(((1 - rr_**2)**2 - 4 * lp * rr_**2).subs(rr_, rl)) == 0
# at r_lam: r^2/(1-r^2) = r/(2 sqrt lam), and r_lam = 1/(sqrt lam + sqrt(lam+1)),
# so ln r_lam = -arsinh(sqrt lam) and max g = lam - sqrt(lam(lam+1)) - arsinh(sqrt lam)
assert sp.simplify(rl**2 / (1 - rl**2) - rl / (2 * sp.sqrt(lp))) == 0
assert sp.simplify(rl * (sp.sqrt(lp) + sp.sqrt(lp + 1)) - 1) == 0
assert sp.simplify(sp.asinh(sp.sqrt(lp)).rewrite(sp.log) - sp.log(sp.sqrt(lp) + sp.sqrt(lp + 1))) == 0
Smax = sp.log(rl) - sp.sqrt(lp) * rl
assert sp.simplify(g.subs(rr_, rl) - Smax) == 0
assert sp.simplify(sp.expand(-sp.sqrt(lp) * rl) - (lp - sp.sqrt(lp * (lp + 1)))) == 0
Dfun = lambda lv, bv, rv: 2 * lv * rv**2 + (bv - sp.log(rv)) * (1 - rv**2)
assert sp.simplify(Dfun(1, -2, sp.Rational(1, 2)) - (-1 + sp.Rational(3, 4) * sp.log(2))) == 0
assert Dfun(1, -2, sp.Rational(1, 2)) < 0 and Dfun(1, -2, 1) > 0
print('14. admissible region |Im beta - lam| < sqrt(lam(lam+1)) + arsinh sqrt(lam);',
      'lam=1, beta=-2i: D(1/2) = -1 + 3/4 ln 2 < 0')

# -- 15. toric dual bases ------------------------------------------------------------------
Amat = sp.Matrix([[1, 1], [0, 1]])
m1 = lambda k: sp.Matrix([-1, k + 1])
m2 = lambda k: sp.Matrix([1, -k])
vv = lambda k: sp.Matrix([k, 1])
for k in (n,):
    assert sp.simplify(Amat.T * m1(k + 1) - m1(k)) == sp.zeros(2, 1)
    assert sp.simplify(Amat.T * m2(k + 1) - m2(k)) == sp.zeros(2, 1)
    assert (m1(k).T * vv(k))[0] == 1 and sp.simplify((m1(k).T * vv(k + 1))[0]) == 0
    assert sp.simplify((m2(k).T * vv(k))[0]) == 0 and sp.simplify((m2(k).T * vv(k + 1))[0]) == 1
print('15. dual basis of (v_n, v_{n+1}) is ((-1,n+1),(1,-n)); A^T m_i^(n+1) = m_i^(n)')
print('all checks passed')
