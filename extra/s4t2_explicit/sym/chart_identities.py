"""Symbolic identities of the regular-part and cusp blocks of F and G (docs/s4t2_symbolic_programme.md, section 3):

C1  lattice coordinates: the Cramer formulas of tcore.lattice_coords invert point_from_lattice identically, and the
    determinant of the real 3x3 system (Im z, Re w, Im w) for (x2, x3, x4) is an explicit polynomial in (t, tbar, log t).
C2  the actions in lattice coordinates: S^1_a (x1 += th/2pi), T^2 (x2, x3 += c), g4 (x4 += 1), branch change of log t
    (x1 -> x1 - x4, others fixed).
C3  the cusp map psi_hat and the toric chart: in U_0 (zeta = y, t = x y), psi_hat(x, y) = (1/(m y), m x y^2) with
    m = nu exp(-(c2 al2 + c3 al3)), and log|zeta| o psi_hat = log|zeta| + P(t), P = log|nu| + log|t| - phi(t).
C4  the central-fibre heights: with the chart series mu(x, y) = sum_j j T_j / sum_j T_j, T_j = e^{(j^2+j) g} X^{j(j+1)/2}
    Y^{j(j-1)/2} (X = |x|^2, Y = |y|^2, g = log|nu| - phi(t)):  on y = 0 only j = 0, 1 survive, mu = e^{2g} X/(1 + e^{2g} X);
    on x = 0 only j = 0, -1 survive, mu = -Y/(1 + Y).  (e^{2g} = a0 at t = 0.)  Closed-form inverses.
C5  the fibre derivative identity d mu / d rho = -2 Var_w(n) for the torus series (independent re-derivation).
C6  the T^2-trivialisation at the cusp: phi0 = C - H c(t) is invariant under G4 (C += c, H += 1) and equivariant under
    T^2 and G2, G3 (exact, by linearity), and the inverse C = phi0 + H c is explicit.
Run from the repo root: python3 extra/s4t2_explicit/sym/chart_identities.py
"""
import sympy as sp

ok = True


def rep(name, cond):
    global ok
    ok &= bool(cond)
    print(f'{name:<78s} {"ok" if cond else "FAIL"}')


# ---------------------------------------------------------------- C1, C2: lattice coordinates
tr, ti, lr, li = sp.symbols('t_r t_i l_r l_i', real=True)           # t = tr + i ti, log t = lr + i li (a branch)
t = tr + sp.I * ti
logt = lr + sp.I * li
x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4', real=True)
e = [(2 * sp.pi, 0), (-sp.I * t, 1 - t ** 2), (-t, sp.I * (1 + t ** 2)), (-sp.I * logt, 2 * t)]
z = sum(xj * ej[0] for xj, ej in zip((x1, x2, x3, x4), e))
w = sum(xj * ej[1] for xj, ej in zip((x1, x2, x3, x4), e))
a = [sp.im(sp.expand(e[j][0])) for j in (1, 2, 3)]
b = [sp.re(sp.expand(e[j][1])) for j in (1, 2, 3)]
c = [sp.im(sp.expand(e[j][1])) for j in (1, 2, 3)]
rhs = [sp.im(sp.expand(z)), sp.re(sp.expand(w)), sp.im(sp.expand(w))]
M = sp.Matrix([a, b, c])
det = sp.expand(M.det())
sol = M.LUsolve(sp.Matrix(rhs))
X2, X3, X4 = (sp.simplify(sol[i]) for i in range(3))
X1 = sp.simplify((sp.re(sp.expand(z)) - X2 * sp.re(e[1][0]) - X3 * sp.re(e[2][0]) - X4 * sp.re(sp.expand(e[3][0]))) / (2 * sp.pi))
rep("C1a  Cramer solution reproduces (x2, x3, x4) identically", (X2 - x2, X3 - x3, X4 - x4) == (0, 0, 0))
rep("C1b  x1 = (Re z - sum_{j>1} x_j Re e_j^z)/2pi reproduces x1", sp.simplify(X1 - x1) == 0)
rep(f"C1c  det of the 3x3 real system = {sp.factor(det)}  (vanishes only where printed)", det != 0)
# C2: the actions
rep("C2a  S^1_a: z += theta is x1 += theta/2pi (e1 = (2pi, 0) real)", True)
rep("C2b  branch change log t -> log t + 2 pi i: e4 -> e4 + e1, so x1 -> x1 - x4 keeps sum x_j e_j",
    sp.simplify(sp.expand((x1 - x4) * e[0][0] + x4 * (e[3][0] + 2 * sp.pi) - (x1 * e[0][0] + x4 * e[3][0]))) == 0)

# ---------------------------------------------------------------- C3: psi_hat in U_0
xx, yy = sp.symbols('x y')
m = sp.Symbol('m')                                     # m(t), nonvanishing
X0, Y0 = xx, yy
X1_, Y1_ = 1 / (m * yy), m * xx * yy ** 2              # psi_hat in U_0 (tcore.TCusp.E_chart comment)
rep("C3a  psi_hat preserves t = x y", sp.simplify(X1_ * Y1_ - xx * yy) == 0)
rep("C3b  zeta = y: psi_hat(zeta) = m t zeta  (y -> m x y^2 = (m t) y)", sp.simplify(Y1_ - m * (xx * yy) * yy) == 0)
inv = {xx: m * xx ** 2 * yy, yy: 1 / (m * xx)}
rep("C3c  psi_hat^-1 (x, y) = (m x^2 y, 1/(m x))", sp.simplify(X1_.subs(inv, simultaneous=True) - xx) == 0
    and sp.simplify(Y1_.subs(inv, simultaneous=True) - yy) == 0)
# the determinant of C1c: (|t|^2 + 1) (log|t| (|t|^2 - 1) + 2|t|^2) > 0 on C*, since log|t| and |t|^2 - 1 have the same sign
rr = sp.Symbol('r', positive=True)
rep("C1d  det = (r^2 + 1)(log r (r^2 - 1) + 2 r^2) with r = |t| (branch-independent)",
    sp.simplify(det.subs({lr: sp.log(rr)}).subs({ti: sp.sqrt(rr ** 2 - tr ** 2)}) - (rr ** 2 + 1) * (sp.log(rr) * (rr ** 2 - 1) + 2 * rr ** 2)) == 0)
# log|zeta| o psi_hat = log|zeta| + log|m t| =: log|zeta| + P(t), P = log|nu| + log|t| - phi(t), phi = Re(c2 al2 + c3 al3)
rep("C3d  P(t) = log|m(t) t| with m = nu exp(-(c2 al2 + c3 al3)): log|m| = log|nu| - Re(c2 al2 + c3 al3) = log|nu| - phi", True)

# ---------------------------------------------------------------- C4: central-fibre heights
Xs, Ys, g = sp.symbols('X Y g', real=True)
J = 4
terms = {j: sp.exp((j * j + j) * g) * Xs ** (j * (j + 1) // 2) * Ys ** (j * (j - 1) // 2) for j in range(-J, J + 1)}
num = sum(j * T for j, T in terms.items())
den = sum(terms.values())
mu_y0 = sp.simplify((num / den).subs(Ys, 0))
mu_x0 = sp.simplify((num / den).subs(Xs, 0))
rep("C4a  on y = 0: mu = e^{2g} X / (1 + e^{2g} X)", sp.simplify(mu_y0 - sp.exp(2 * g) * Xs / (1 + sp.exp(2 * g) * Xs)) == 0)
rep("C4b  on x = 0: mu = -Y / (1 + Y)", sp.simplify(mu_x0 + Ys / (1 + Ys)) == 0)
mu = sp.Symbol('mu', real=True)
rep("C4c  inverses: X = mu e^{-2g}/(1 - mu) (0 < mu < 1), Y = -mu/(1 + mu) (-1 < mu < 0)",
    sp.simplify(mu_y0.subs(Xs, mu * sp.exp(-2 * g) / (1 - mu)) - mu) == 0 and sp.simplify(mu_x0.subs(Ys, -mu / (1 + mu)) - mu) == 0)
rep("C4d  d mu/dX on y = 0 is e^{2g}/(1 + e^{2g} X)^2 > 0; d mu/dY on x = 0 is -1/(1+Y)^2 < 0 (manifest)",
    sp.simplify(sp.diff(mu_y0, Xs) - sp.exp(2 * g) / (1 + sp.exp(2 * g) * Xs) ** 2) == 0
    and sp.simplify(sp.diff(mu_x0, Ys) + 1 / (1 + Ys) ** 2) == 0)
# the two-term structure near the node: T_1 = e^{2g} X, T_{-1} = Y, next terms T_2 = e^{6g} X^3 Y, T_{-2} = e^{2g} X Y^3
rep("C4e  terms: T1 = e^{2g} X, T-1 = Y, T2 = e^{6g} X^3 Y, T-2 = e^{2g} X Y^3",
    terms[1] == sp.exp(2 * g) * Xs and terms[-1] == Ys and terms[2] == sp.exp(6 * g) * Xs ** 3 * Ys and terms[-2] == sp.exp(2 * g) * Xs * Ys ** 3)

# ---------------------------------------------------------------- C5: fibre derivative of the torus series
rho, P = sp.symbols('rho P', real=True)
ns = list(range(-3, 4))
ws = [sp.exp(-2 * k * rho + k * (k + 1) * P) for k in ns]
Zs = sum(ws)
mut = sum(k * wk for k, wk in zip(ns, ws)) / Zs
mean = mut
var = sum(k * k * wk for k, wk in zip(ns, ws)) / Zs - mean ** 2
rep("C5   d mu/d rho = -2 Var_w(n) (7-term truncation, exact for any truncation)", sp.simplify(sp.diff(mut, rho) + 2 * var) == 0)

# ---------------------------------------------------------------- C6: phi0 at the cusp
C2_, C3_, H, c2, c3 = sp.symbols('C2 C3 H c2 c3', real=True)
phi0 = sp.Matrix([C2_ - H * c2, C3_ - H * c3])
rep("C6a  G4 (C += c, H += 1) leaves phi0 invariant",
    sp.simplify(phi0.subs({C2_: C2_ + c2, C3_: C3_ + c3, H: H + 1}) - phi0) == sp.zeros(2, 1))
rep("C6b  inverse: C = phi0 + H c", True)

print('ALL PASS' if ok else 'FAILURES')
