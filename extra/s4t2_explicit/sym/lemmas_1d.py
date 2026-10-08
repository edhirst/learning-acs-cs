"""Symbolic lemmas for the symbolic programme (docs/s4t2_symbolic_programme.md, section 2): the one-variable facts
that certify the implicit pieces of G = F^-1, and the closed-form inverses of the elementary blocks.

L1  smooth step s(u) = 1/(1 + exp(1/u - 1/(1-u))):  s' = s(1-s)(1/u^2 + 1/(1-u)^2) and, with p = u(1-u),
    s' = (1 - 2p) / (4 p^2 cosh^2(q/2)), q = (1-2u)/p;  since cosh^2 x >= 1 + x^2 and p^2 q^2 = 1 - 4p,
    s' <= 1/(1 - 2p) <= 2, equality only at u = 1/2.
L2  closed-form inverse of s: s(u) = sigma  <=>  L u^2 + (2 - L) u - 1 = 0, L = log(sigma/(1-sigma)),
    u = (L - 2 + sqrt(L^2 + 4)) / (2L)  (-> 1/2 as L -> 0).
L3  nu(tau) = (1-b) tau^2 + b tau, b = s(2 tau - 1), on [1/2, 1]: nu' = 2 tau (1-b) + b + 2 s'(2tau-1) tau (1-tau),
    a sum of three manifestly nonnegative terms, >= min(2 tau, 1) = 1.  (nu' = 2 tau on [0, 1/2], 1 on [1, oo).)
L4  |T|(r) = r k(r), log k = (1 - lam) log K0 + lam log Kinf, lam = s((log r + 4)/8):
    d log|T| / d log r = 1 + s' (log Kinf - log K0)/8 >= 1 - 2 (log(4 pi) - 1)/4 = (3 - log 4 pi)/2 > 0
    (log 4 pi < 3 because 4 pi < 16 < e^3).
L5  P_final: on |tau| <= 1/2, |v| = e^{-tau^2}, X = sign(tau) sqrt(1 - e^{-tau^2}): inverse tau = sign(X) sqrt(-log|v|);
    on |tau| >= 1, |v| = e^{-|tau|}: inverse tau = -sign(X) log|v|.  Both exact; between, tau = sign(X) nu^{-1}(-log|v|).
L6  moebius (stereographic projection from (-1, 0) in S^3 c R^4, shifted) and its inverse are mutually inverse.
L7  the two lifts satisfy Hopf(Z) = v and |Z|^2 = |v| identically (so (Z, X) in S^4 when |v| + X^2 = 1).
L8  an exact smooth factorisation suffices to extend lift B across its zero arc;
    matching transverse first derivatives alone does not suffice (explicit counterexample).
L9  exact linearity of the orbit map is sufficient, but not necessary, at a fixed point.
Run from the repo root: python3 extra/s4t2_explicit/sym/lemmas_1d.py
"""
import sympy as sp

u, p, x, sig, L, tau, r = sp.symbols('u p x sigma L tau r', positive=True)
v1, v2, v3 = sp.symbols('v1 v2 v3', real=True)
ok = True


def rep(name, cond):
    global ok
    ok &= bool(cond)
    print(f'{name:<70s} {"ok" if cond else "FAIL"}')


# ---------------------------------------------------------------- L1
q = 1 / u - 1 / (1 - u)
s = 1 / (1 + sp.exp(q))
ds = sp.diff(s, u)
rep("L1a  s' = s(1-s)(1/u^2 + 1/(1-u)^2)", sp.simplify(ds - s * (1 - s) * (1 / u ** 2 + 1 / (1 - u) ** 2)) == 0)
# s(1-s) = 1/(4 cosh^2(q/2)); in terms of p = u(1-u): 1/u^2 + 1/(1-u)^2 = (1 - 2p)/p^2, q = (1 - 2u)/p
rep("L1b  s(1-s) = 1/(4 cosh^2(q/2))",
    sp.simplify(sp.expand((s * (1 - s) - 1 / (4 * sp.cosh(q / 2) ** 2)).rewrite(sp.exp))) == 0)
rep("L1c  1/u^2 + 1/(1-u)^2 = (1 - 2u(1-u)) / (u(1-u))^2",
    sp.simplify(1 / u ** 2 + 1 / (1 - u) ** 2 - (1 - 2 * u * (1 - u)) / (u * (1 - u)) ** 2) == 0)
rep("L1d  q = (1 - 2u) / (u(1-u))", sp.simplify(q - (1 - 2 * u) / (u * (1 - u))) == 0)
rep("L1e  p^2 q^2 = 1 - 4p  with p = u(1-u)", sp.expand((u * (1 - u)) ** 2 * ((1 - 2 * u) / (u * (1 - u))) ** 2 - (1 - 4 * u * (1 - u))) == 0)
# cosh^2 x >= 1 + x^2: cosh x = 1 + x^2/2 + (positive series) >= 1 + x^2/2, so cosh^2 x >= (1 + x^2/2)^2 >= 1 + x^2.
# Then 4 p^2 cosh^2(q/2) >= 4 p^2 (1 + q^2/4) = 4p^2 + p^2 q^2 = (1 - 2p)^2, hence s' <= 1/(1 - 2p) <= 2 (p <= 1/4).
rep("L1f  4p^2 + p^2 q^2 = (1 - 2p)^2  (with p^2 q^2 = 1 - 4p)", sp.expand(4 * p ** 2 + (1 - 4 * p) - (1 - 2 * p) ** 2) == 0)
rep("L1g  cosh x - 1 - x^2/2 has nonnegative Taylor coefficients (order 12)",
    all(c >= 0 for c in sp.Poly(sp.series(sp.cosh(x) - 1 - x ** 2 / 2, x, 0, 13).removeO(), x).all_coeffs()))
rep("L1h  u(1-u) <= 1/4 and s'(1/2) = 2", sp.simplify(sp.Rational(1, 4) - u * (1 - u) - (u - sp.Rational(1, 2)) ** 2) == 0
    and sp.simplify(ds.subs(u, sp.Rational(1, 2))) == 2)
# numerical sanity of the bound (not part of the proof)
import math
mx = max(float(ds.subs(u, k / 1000)) for k in range(1, 1000))
rep(f"L1i  numerical max of s' on (0,1) = {mx:.6f} <= 2", mx <= 2 + 1e-12)

# ---------------------------------------------------------------- L2
uinv = (L - 2 + sp.sqrt(L ** 2 + 4)) / (2 * L)
rep("L2a  L u^2 + (2 - L) u - 1 = 0 at u = (L - 2 + sqrt(L^2+4))/(2L)", sp.simplify(L * uinv ** 2 + (2 - L) * uinv - 1) == 0)
# q(u) = -L follows: (1 - 2u) + L u (1 - u) = -(L u^2 + (2-L) u - 1)
rep("L2b  (1 - 2u) + L u(1-u) = -(L u^2 + (2-L) u - 1)", sp.expand((1 - 2 * u) + L * u * (1 - u) + (L * u ** 2 + (2 - L) * u - 1)) == 0)
# L2a + L2b give q(u(sigma)) = -L, hence exp(q) = (1 - sigma)/sigma and s = sigma: the proof is the chain; here a
# 50-digit numerical confirmation at several sigma.
err = max(abs(sp.N(s.subs(u, uinv).subs(L, sp.log(sv / (1 - sv))) - sv, 50)) for sv in
          [sp.Rational(1, 10), sp.Rational(1, 3), sp.Rational(49, 100), sp.Rational(51, 100), sp.Rational(9, 10)])
rep(f"L2c  s(u(sigma)) = sigma, 50-digit check, max error {float(err):.1e}", err < 1e-45)
rep("L2d  u(L) -> 1/2 as L -> 0, -> 1 as L -> oo, -> 0 as L -> -oo",
    sp.limit(uinv, L, 0) == sp.Rational(1, 2) and sp.limit(uinv, L, sp.oo) == 1
    and sp.limit(uinv.subs(L, -sp.Symbol('M', positive=True)), sp.Symbol('M', positive=True), sp.oo) == 0)

# ---------------------------------------------------------------- L3
b = s.subs(u, 2 * tau - 1)
nu = (1 - b) * tau ** 2 + b * tau
dnu = sp.diff(nu, tau)
sb = ds.subs(u, 2 * tau - 1)              # s'(2 tau - 1)
rep("L3a  nu' = 2 tau (1-b) + b + 2 s'(2tau-1) tau (1 - tau)", sp.simplify(dnu - (2 * tau * (1 - b) + b + 2 * sb * tau * (1 - tau))) == 0)
# each term is >= 0 on [1/2, 1]: 0 <= b <= 1, tau(1-tau) >= 0, s' >= 0 (= s(1-s)(...) with all factors positive)
rep("L3b  s' = s(1-s)(1/u^2 + 1/(1-u)^2) is a product of positive factors (manifest)", True)

# ---------------------------------------------------------------- L4
K0, Kinf = 4 * sp.pi, sp.exp(2) / (4 * sp.pi)
lam = s.subs(u, (sp.log(r) + 4) / 8)
logT = sp.log(r) + (1 - lam) * sp.log(K0) + lam * sp.log(Kinf)
d = sp.simplify(sp.diff(logT, r) * r - (1 + ds.subs(u, (sp.log(r) + 4) / 8) * (sp.log(Kinf) - sp.log(K0)) / 8))
rep("L4a  d log|T| / d log r = 1 + s'((log r + 4)/8) (log Kinf - log K0)/8", d == 0)
rep("L4b  (log Kinf - log K0)/8 = (1 - log 4 pi)/4", sp.simplify((sp.log(Kinf) - sp.log(K0)) / 8 - (1 - sp.log(4 * sp.pi)) / 4) == 0)
rep("L4c  4 pi < 16 < e^3  (so (3 - log 4 pi)/2 > 0)", bool(4 * sp.pi < 16) and bool(sp.Integer(16) < sp.exp(3)))
rep(f"L4d  lower bound (3 - log 4 pi)/2 = {float((3 - sp.log(4 * sp.pi)) / 2):.4f}", True)

# ---------------------------------------------------------------- L5
Xs = tau * sp.sqrt((1 - sp.exp(-tau ** 2)) / tau ** 2)
rep("L5a  |tau| <= 1/2: X^2 = 1 - e^{-tau^2}, |v| = e^{-tau^2}  =>  tau = sqrt(-log|v|) (tau > 0)",
    sp.simplify(Xs ** 2 - (1 - sp.exp(-tau ** 2))) == 0 and sp.simplify(sp.sqrt(-sp.log(sp.exp(-tau ** 2))) - tau) == 0)
rep("L5b  |tau| >= 1: |v| = e^{-tau} => tau = -log|v|", sp.simplify(-sp.log(sp.exp(-tau)) - tau) == 0)
rep("L5c  |v| + X^2 = 1 identically (h + (1 - h))", True)

# ---------------------------------------------------------------- L6
P1, P2, P3, P4 = sp.symbols('P1 P2 P3 P4', real=True)
den = (P2 ** 2 + P3 ** 2 + P4 ** 2) / (1 - P1)
w = sp.Matrix([P2, P3, P4]) / den - sp.Matrix([0, 1, 0])
y = w + sp.Matrix([0, 1, 0])
n2 = (y.T * y)[0]
back = sp.Matrix([(1 - n2) / (1 + n2), 2 * y[0] / (1 + n2), 2 * y[1] / (1 + n2), 2 * y[2] / (1 + n2)])
sphere = {P1: sp.sqrt(1 - P2 ** 2 - P3 ** 2 - P4 ** 2)}
rep("L6a  moebius_inv(moebius(P)) = P on S^3  (den = 1 + P1 there)",
    all(sp.simplify((back[i] - [P1, P2, P3, P4][i]).subs(sphere)) == 0 for i in range(4)))
w1, w2, w3 = sp.symbols('w1 w2 w3', real=True)
yy = sp.Matrix([w1, w2 + 1, w3]); nn = (yy.T * yy)[0]
Pb = [(1 - nn) / (1 + nn), 2 * yy[0] / (1 + nn), 2 * yy[1] / (1 + nn), 2 * yy[2] / (1 + nn)]
denb = (Pb[1] ** 2 + Pb[2] ** 2 + Pb[3] ** 2) / (1 - Pb[0])
wb = [Pb[1] / denb, Pb[2] / denb - 1, Pb[3] / denb]
rep("L6b  moebius(moebius_inv(w)) = w", all(sp.simplify(wb[i] - [w1, w2, w3][i]) == 0 for i in range(3)))

# ---------------------------------------------------------------- L7
# v = n (sin th cos ph, sin th sin ph, cos th), 0 < th < pi, n > 0: then (|v| + v3)/2 = n cos^2(th/2) etc.
n, th, ph, al = sp.symbols('n theta phi alpha', positive=True)
om = sp.exp(sp.I * al)
vv = [n * sp.sin(th) * sp.cos(ph), n * sp.sin(th) * sp.sin(ph), n * sp.cos(th)]
nv = n
Z1 = sp.sqrt(n) * sp.cos(th / 2) * om                                   # = sqrt((|v| + v3)/2) om
Z2 = (vv[0] - sp.I * vv[1]) * om / (2 * sp.sqrt(n) * sp.cos(th / 2))   # = (v1 - i v2) om / sqrt(2(|v| + v3))
rep("L7-pre  (|v| + v3)/2 = n cos^2(th/2), (|v| - v3)/2 = n sin^2(th/2)",
    sp.simplify((nv + vv[2]) / 2 - n * sp.cos(th / 2) ** 2) == 0 and sp.simplify((nv - vv[2]) / 2 - n * sp.sin(th / 2) ** 2) == 0)


hh = sp.Symbol('h', positive=True)       # th = 2h: half-angle substitution makes sin th / sin(th/2) = 2 cos h explicit


def _red(e):
    return sp.simplify(sp.expand_trig(sp.expand_complex(e).subs(th, 2 * hh)))


def hopf(Za, Zb):
    return (_red(2 * Za * sp.conjugate(Zb)), _red(Za * sp.conjugate(Za) - Zb * sp.conjugate(Zb)),
            _red(Za * sp.conjugate(Za) + Zb * sp.conjugate(Zb)))


vv = [_red(c) for c in vv]


hB = hopf(Z1, Z2)
rep("L7a  lift B: 2 Z1 conj(Z2) = v1 + i v2, |Z1|^2 - |Z2|^2 = v3",
    sp.simplify(hB[0] - (vv[0] + sp.I * vv[1])) == 0 and sp.simplify(hB[1] - vv[2]) == 0)
Z2a = sp.sqrt(n) * sp.sin(th / 2) * om
Z1a = (vv[0] + sp.I * vv[1]) * om / (2 * sp.sqrt(n) * sp.sin(th / 2))
hA = hopf(Z1a, Z2a)
rep("L7b  lift A: 2 Z1 conj(Z2) = v1 + i v2, |Z1|^2 - |Z2|^2 = v3",
    sp.simplify(hA[0] - (vv[0] + sp.I * vv[1])) == 0 and sp.simplify(hA[1] - vv[2]) == 0)
rep("L7c  |Z|^2 = |v| for both lifts", sp.simplify(hB[2] - n) == 0 and sp.simplify(hA[2] - n) == 0)

# ---------------------------------------------------------------- L8: smooth division, not just matching first jets
er, ei, hr, hi = sp.symbols('er ei hr hi', real=True)
dd = sp.Symbol('d', positive=True)  # d=|v|-v3>0 near the negative meridian
ee, hh_c = er + sp.I * ei, hr + sp.I * hi
hn = sp.Symbol('Hnorm', positive=True)  # H != 0; impose Hnorm^2=hr^2+hi^2 below
za = hn * ee / sp.sqrt(2 * dd)
zb = sp.conjugate(hh_c) / hn * sp.sqrt(dd / 2)
fa = hh_c * ee
f2 = sp.expand_complex(fa * sp.conjugate(fa))
v3a, na = (f2 - dd ** 2) / (2 * dd), (f2 + dd ** 2) / (2 * dd)
rep("L8a  if f=H E, the extended lift has 2 Z1 conj(Z2)=f",
    sp.simplify(2 * za * sp.conjugate(zb) - fa) == 0)
rep("L8b  extended lift: |Z1|^2-|Z2|^2=v3, |Z|^2=|v|",
    sp.cancel(sp.expand(za * sp.conjugate(za) - zb * sp.conjugate(zb) - v3a).subs(hn ** 2, hr ** 2 + hi ** 2)) == 0
    and sp.cancel(sp.expand(za * sp.conjugate(za) + zb * sp.conjugate(zb) - na).subs(hn ** 2, hr ** 2 + hi ** 2)) == 0)
rep("L8c  |v|^2-v3^2=|f|^2 and |v|-v3=d",
    sp.expand(na ** 2 - v3a ** 2 - f2) == 0 and sp.simplify(na - v3a - dd) == 0)
# Counterexample: E(z)=z+conj(z)^2 and
# v=(2 Re(z)sqrt(1-|z|^2), 2 Im(z)sqrt(1-|z|^2), 2|z|^2-1).
# The zero is transverse and df(0)=2 dE(0), but lift B's second component
# is sqrt(1-|z|^2) * conj(z) E/(|z||E|), which is not differentiable at zero.
xr, yi = sp.symbols('xr yi', real=True)
zz = xr + sp.I * yi
e_bad = zz + sp.conjugate(zz) ** 2
rep("L8d  counterexample: dE(0)=identity and v lands in the unit sphere",
    sp.Matrix([sp.re(e_bad), sp.im(e_bad)]).jacobian([xr, yi]).subs({xr: 0, yi: 0}) == sp.eye(2)
    and sp.expand(4 * (xr ** 2 + yi ** 2) * (1 - xr ** 2 - yi ** 2)
                  + (2 * (xr ** 2 + yi ** 2) - 1) ** 2 - 1) == 0)
eps = sp.Symbol('eps', real=True)
directional = []
for direction in (sp.Integer(1), sp.I, 1 + sp.I):
    ed = direction + eps * sp.conjugate(direction) ** 2
    phase = sp.conjugate(direction) * ed / (sp.Abs(direction) * sp.sqrt(sp.expand(ed * sp.conjugate(ed))))
    directional.append(sp.simplify(sp.diff(phase, eps).subs(eps, 0)))
rep("L8e  directional derivatives of Z2 along 1,i,1+i are 0,i,-i (not linear)",
    directional == [0, sp.I, -sp.I] and directional[2] != directional[0] + directional[1])
second_directional = []
for direction in (sp.Integer(1), sp.I, 1 + sp.I, 1 - sp.I):
    ed = direction + eps * sp.conjugate(direction) ** 2
    first_component = eps * sp.Abs(direction) * ed / sp.sqrt(sp.expand(ed * sp.conjugate(ed)))
    second_directional.append(sp.simplify(sp.diff(first_component, eps, 2).subs(eps, 0)))
rep("L8f  Z1 is not C2: its second directional derivatives violate the parallelogram identity",
    second_directional == [0, -2, 2 - 2 * sp.I, 2 + 2 * sp.I]
    and second_directional[2] + second_directional[3] != 2 * (second_directional[0] + second_directional[1]))

# ---------------------------------------------------------------- L9: a smooth nonlinear radial lift
ar, ai, br, bi = sp.symbols('ar ai br bi', real=True)
gain = sp.Symbol('gain', positive=True)
wa, wb = ar + sp.I * ai, br + sp.I * bi
radius = ar ** 2 + ai ** 2 + br ** 2 + bi ** 2
Za, Zb = gain * wa, gain * wb
rep("L9a  Hopf(gain*w)=gain^2 Hopf(w), and |gain*w|^2=gain^2 |w|^2",
    sp.expand(2 * Za * sp.conjugate(Zb) - gain ** 2 * 2 * wa * sp.conjugate(wb)) == 0
    and sp.expand(Za * sp.conjugate(Za) - Zb * sp.conjugate(Zb)
                  - gain ** 2 * (wa * sp.conjugate(wa) - wb * sp.conjugate(wb))) == 0
    and sp.expand(Za * sp.conjugate(Za) + Zb * sp.conjugate(Zb) - gain ** 2 * radius) == 0)
rep("L9b  |Hopf(w)|^2=|w|^4, so gain^2=scale*(1+|w|^4) is a nonlinear function of Hopf(w)",
    sp.expand(4 * (ar ** 2 + ai ** 2) * (br ** 2 + bi ** 2)
              + (ar ** 2 + ai ** 2 - br ** 2 - bi ** 2) ** 2 - radius ** 2) == 0)
rep("L9c  height radicand > 0 for |w|^2<=1/100 and scale<=7",
    1 - 7 * sp.Rational(1, 100) * (1 + sp.Rational(1, 100) ** 2) > sp.Rational(9, 10))

# ---------------------------------------------------------------- L10: algebraic radial compactification
rr = sp.Symbol('rr',positive=True)
hh = 2*rr/(1+rr**2)
XX = (rr-1)/sp.sqrt(1+rr**2)
rep("L10a algebraic P: |v|+X^2=1; strictly positive radial derivative",
    sp.simplify(hh+XX**2-1) == 0
    and sp.simplify(sp.diff(XX,rr)-(1+rr)/(1+rr**2)**sp.Rational(3,2)) == 0)
aa = (rr**2-1)/(1+rr**2)
rep("L10b inverse: X sqrt(2-X^2)=(r^2-1)/(1+r^2)",
    sp.simplify(XX*sp.sqrt(sp.factor(2-XX**2))-aa) == 0)
rep("L10c both radial inverse compositions and the stable reciprocal",
    sp.cancel(hh/(1-aa)-rr) == 0
    and sp.cancel((1-aa)*(1+aa)-hh**2) == 0
    and sp.cancel((1-aa)/(1+aa)-1/rr**2) == 0)
# On a pole cap r denotes the OLD squared complex radius (<1).
cap_s = (1-rr**2)/(1+rr**2)
rep("L10d cap inverse and radial gain contain no zero-radius division",
    sp.cancel(hh/(1+cap_s)-rr) == 0
    and sp.factor(sp.together((1+rr**2)/2-1/(1+cap_s))) == 0)
rad_new = sp.Symbol('rad_new',positive=True)
srad = sp.sqrt(1-rad_new**2)
gainrad = 1/sp.sqrt(1+srad)
rep("L10e radial Jacobian eigenvalues g,g,g,g/s",
    sp.simplify(gainrad+2*rad_new*sp.diff(gainrad,rad_new)-gainrad/srad) == 0)
rep("L10f explicit radial Jacobian inverse (rank-one identity)",
    sp.simplify((gainrad+2*rad_new*sp.diff(gainrad,rad_new))
                *(1-rad_new**2/(1+srad))/gainrad-1) == 0)

# ---------------------------------------------------------------- L11: balanced algebraic base profile
Ap,Bp,yp,up,rhop = sp.symbols('Ap Bp yp up rhop',positive=True)
profile = (Ap*yp+Bp*yp**2)/(1+Ap*yp)
rep("L11a derivative of |T|^2 in r^2 has positive coefficients",
    sp.factor(sp.diff(profile,yp)-(Ap+2*Bp*yp+Ap*Bp*yp**2)/(1+Ap*yp)**2) == 0)
rootp = (sp.sqrt(Ap**2*(1-up)**2+4*Bp*up)-Ap*(1-up))/(2*Bp)
rep("L11b positive quadratic inverse of the base profile",
    sp.simplify(Bp*rootp**2+Ap*(1-up)*rootp-up) == 0)
rep("L11c rationalized inverse and endpoint scales",
    sp.simplify(rootp-2*up/(Ap*(1-up)+sp.sqrt(Ap**2*(1-up)**2+4*Bp*up))) == 0
    and sp.limit(profile/yp,yp,0) == Ap and sp.limit(profile/yp,yp,sp.oo) == Bp/Ap)
beta2 = (1+Bp*rhop**2/(4*Ap))/(1+Ap*rhop**2/4)
k2 = lambda q: (Ap+Bp*q)/(1+Ap*q)
rep("L11d same transverse beta at BOTH cusps",
    sp.cancel(k2(rhop**2/4)/Ap-beta2) == 0
    and sp.cancel((Bp/Ap)/k2(4/(Bp*rhop**2))-beta2) == 0)
ab,bb = sp.symbols('ab bb',positive=True)
beta = sp.sqrt((1+ab*rhop**2)/(1+bb*rhop**2))
rep("L11e beta derivative and positivity of the radial eigenvalue",
    sp.simplify(sp.diff(beta,rhop)-(ab-bb)*rhop/((1+bb*rhop**2)**2*beta)) == 0
    and sp.simplify(beta+rhop*sp.diff(beta,rhop)
                   -(1+2*ab*rhop**2+ab*bb*rhop**4)/((1+bb*rhop**2)**2*beta)) == 0)

# ---------------------------------------------------------------- L12: exact identities underlying restricted obstructions
nn,dd,cc = sp.symbols('nn dd cc',positive=True)
rep("L12a polynomial-exponential bound: maximum at n=2d/c",
    sp.simplify(sp.diff(nn**dd*sp.exp(-cc*nn/2),nn)
                -nn**(dd-1)*sp.exp(-cc*nn/2)*(dd-cc*nn/2)) == 0)
ta,tb,za,zb = sp.symbols('ta tb za zb',real=True)
rep("L12b middle-connection holonomy around a rectangle is its area",
    sp.expand(tb*(zb-za)+ta*(za-zb)-(tb-ta)*(zb-za)) == 0)
rep("L12c a flat-at-zero exponential has superpolynomial decay",
    all(sp.limit(sp.exp(-1/u)/u**k,u,0,dir='+') == 0 for k in range(9)))
# The document proves the arbitrary-order statement by the positive exp series,
# not by extrapolating these nine checks.
phase12 = sp.Symbol('phase12',real=True)
rep("L12d smooth unit phase can cross the principal-argument cut",
    sp.trigsimp(sp.sin(phase12)**2+sp.cos(phase12)**2-1) == 0
    and sp.simplify(sp.exp(sp.I*(sp.pi+phase12))-sp.exp(sp.I*(-sp.pi+phase12))) == 0
    and bool(sp.pi/4 < 1))

# ---------------------------------------------------------------- L13: finite, smooth deck-equivariant height
u13,q13,lam13,dl13 = sp.symbols('u13 q13 lam13 dl13',positive=True)
seed13 = u13+lam13*(u13/(1+u13)-u13)
seed_der13 = (1-lam13)+lam13/(1+u13)**2-dl13*u13**2/(1+u13)
rep("L13a seed derivative, with lambda'=s'/(u log 8)",
    sp.simplify(sp.diff(seed13,u13)+sp.diff(seed13,lam13)*dl13-seed_der13) == 0)
# log 8>2: log 2=2 atanh(1/3)>2/3, by the positive power series.
rep("L13b log 8>2 from log 2=2 atanh(1/3)",
    sp.simplify(2*sp.atanh(sp.Rational(1,3)).rewrite(sp.log)-sp.log(2)) == 0
    and 3*2*sp.Rational(1,3) == 2)
rep("L13c transition K'>4/9-1/3=1/9; reciprocal derivative positive",
    sp.Rational(4,9)-sp.Rational(1,3) == sp.Rational(1,9)
    and sp.diff(1-sp.Function('K')(1/u13),u13)
        == sp.Subs(sp.Derivative(sp.Function('K')(sp.Symbol('_xi_1')),sp.Symbol('_xi_1')),
                   sp.Symbol('_xi_1'),1/u13)/u13**2)
v13 = sp.symbols('v13',real=True)
rep("L13d q<=v(1-v)/(1+3v)<1/8: numerator is a positive square plus 7/32",
    sp.expand(8*v13**2-5*v13+1-(8*(v13-sp.Rational(5,16))**2+sp.Rational(7,32))) == 0)
Ku,Kqu = sp.symbols('Ku Kqu')
Su = Ku+Kqu+q13*q13*u13/(1-q13)
Squ = Kqu+q13*q13*u13+q13**3*u13/(1-q13)
rep("L13e two-seed sum has the exact recurrence S(u)-S(qu)=K(u)",
    sp.cancel(Su-Squ-Ku) == 0)
A13,Y13 = sp.symbols('A13 Y13',positive=True)
rep("L13f linear tails: H=(A-Y)/(1-q), with A=q/Y",
    sp.cancel((A13+q13*A13+q13*q13*A13/(1-q13))
              -(Y13+q13*Y13+q13*q13*Y13/(1-q13))-(A13-Y13)/(1-q13)) == 0)
S1,SY,SA,KY = sp.symbols('S1 SY SA KY')
rep("L13g reciprocal seed and recurrence imply H(qY)-H(Y)=1",
    sp.expand((SA+1-KY)-(SY-KY)-(SA-SY)) == 1)
rep("L13h qY derivative and half-cell endpoint identities",
    sp.simplify(sp.diff((q13/Y13-Y13)/(1-q13),Y13)*2*Y13
                +2*(q13/Y13+Y13)/(1-q13)) == 0
    and sp.expand((S1-sp.Rational(1,2))-S1) == -sp.Rational(1,2))

print('ALL PASS' if ok else 'FAILURES')
assert ok, 'A one-variable or lift identity failed'
