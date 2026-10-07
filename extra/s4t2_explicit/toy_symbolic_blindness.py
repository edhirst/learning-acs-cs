"""Symbolic checks of J^2 = -1 and N = 0 cannot detect a bad transporting map (docs/s4t2_explicit_diffeo.md, 12.13).

Toy: F: R^4 -> R^4, F(x, y, u, v) = (x^2, y + x u, u, v), is two-to-one ((x, y, u, v) and (-x, y + 2 x u, u, v) have
the same image) and degenerate on {x = 0} (det DF = 2x). Its pushforward J = F_* J0 of the standard complex structure
of C^2 satisfies J^2 = -1 and N_J = 0 IDENTICALLY (sympy, target derivatives by the chain rule), yet J is not a
tensor field on the target: the two preimages of a point give different J, and J blows up on F({x = 0}).
Run from the repo root."""
import sympy as sp

x, y, u, v = sp.symbols('x y u v', real=True)
p = [x, y, u, v]
F = sp.Matrix([x ** 2, y + x * u, u, v])
A = F.jacobian(p)
J0 = sp.Matrix([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]])
Ai = sp.simplify(A.inv())
J = sp.simplify(A * J0 * Ai)


def dq(expr, d):
    """Derivative along the target coordinate q_d: d/dq_d = sum_k (DF^-1)_{kd} d/dp_k."""
    return sum(Ai[k, d] * sp.diff(expr, p[k]) for k in range(4))


N = [[[sp.simplify(sum(J[d, b] * dq(J[a, c], d) - J[d, c] * dq(J[a, b], d)
                       - J[a, d] * (dq(J[d, c], b) - dq(J[d, b], c)) for d in range(4)))
       for c in range(4)] for b in range(4)] for a in range(4)]
print('det DF =', sp.factor(A.det()))
print('J = F_* J0 =', J.tolist())
print('J^2 + 1 == 0 identically:', sp.simplify(J * J + sp.eye(4)) == sp.zeros(4))
print('N_J == 0 identically (all 64 components):', all(N[a][b][c] == 0 for a in range(4) for b in range(4) for c in range(4)))
pt = {x: sp.Rational(1, 2), y: sp.Rational(1, 3), u: sp.Rational(2, 5), v: 0}
pt2 = {x: -pt[x], y: pt[y] + 2 * pt[x] * pt[u], u: pt[u], v: 0}
print('two preimages of the same target point:', list(F.subs(pt)) == list(F.subs(pt2)))
print('  J from the first :', J.subs(pt).tolist())
print('  J from the second:', J.subs(pt2).tolist(), ' -> not single-valued')
print('  J[1,0] =', J[1, 0], ' -> unbounded as x -> 0 (on the image of the fold)')
