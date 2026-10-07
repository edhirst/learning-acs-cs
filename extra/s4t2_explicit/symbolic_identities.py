"""Symbolic (sympy) proofs of identities of the explicit construction (docs/s4t2_explicit_diffeo.md, section 12.12):
  (a) the cusp height is strictly monotone along the fibres: d mu / d rho = -2 Var_w(n) (rho = log|zeta|), for
      mu = sum n w_n / sum w_n, w_n = exp(-2 n rho + n (n + 1) P) (any finite truncation; the series by dominated
      convergence);
  (b) the differential of Lambda at the fixed point b0, in Hopf coordinates v_q, is exactly 2 pi Q0,
      Q0 = [[0, -1, 0], [0, 0, -1], [1, 0, 0]] (det = 8 pi^3), from the exact local formulas
      (a1 = v3, T = 4 pi t, R = id, moebius from (-1, 0), P_final = rho_O near 0).
Run from the repo root."""
import sympy as sp

rho, P = sp.symbols('rho P', real=True)
ns = range(-3, 4)
w = {n: sp.exp(-2 * n * rho + n * (n + 1) * P) for n in ns}
Z = sum(w.values())
mu = sum(n * w[n] for n in ns) / Z
var = sum(n ** 2 * w[n] for n in ns) / Z - mu ** 2
print('(a) d mu / d rho + 2 Var_w(n) =', sp.simplify(sp.diff(mu, rho) + 2 * var))

v1, v2, v3 = sp.symbols('v1 v2 v3', real=True)
t = (v1 - sp.I * v2) / 2
T = 4 * sp.pi * t
u2 = sp.exp(2 * sp.pi * sp.I * v3) / sp.sqrt(1 + T * sp.conjugate(T))
u1 = T * u2
P1, P2, P3, P4 = sp.re(u1), sp.im(u1), sp.re(u2), sp.im(u2)
wv = sp.Matrix([P2, P3, P4]) / (1 + P1) - sp.Matrix([0, 1, 0])
vS = sp.Matrix([wv[0], -wv[2], -wv[1]])
Jac = vS.jacobian([v1, v2, v3]).subs({v1: 0, v2: 0, v3: 0})
print('(b) D Lambda(b0) / (2 pi) =', sp.simplify(Jac / (2 * sp.pi)).tolist(), '; det =', sp.simplify(Jac.det()))
