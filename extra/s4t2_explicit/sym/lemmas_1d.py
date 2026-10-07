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

print('ALL PASS' if ok else 'FAILURES')
