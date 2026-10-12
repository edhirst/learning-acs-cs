"""Symbolic (sympy) proofs of identities of the explicit construction (docs/s4t2_explicit_diffeo.md, section 12.12):
  (a) the cusp height is strictly monotone along the fibres: d mu / d rho = -2 Var_w(n) (rho = log|zeta|), for
      mu = sum n w_n / sum w_n, w_n = exp(-2 n rho + n (n + 1) P) (any finite truncation; the series by dominated
      convergence);
  (b) the differential of Lambda at the fixed point b0, in Hopf coordinates v_q, is exactly 2 pi Q0,
      Q0 = [[0, -1, 0], [0, 0, -1], [1, 0, 0]] (det = 8 pi^3), from the exact local formulas
      (a1 = v3, T = 4 pi t, R = id, moebius from (-1, 0), P_final = rho_O near 0).
  (c) on 1/2 <= |t| <= 2, the pure-middle torus-restoration factor of G_cl
      has the six-real-variable formula h=(z,w,t) and an invertible Jacobian.
      The induced tensor J=(Dh)^-1 J0 Dh satisfies J^2=-Id and all 216
      Nijenhuis identities.  Round-to-W and collar factors are outside scope.
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

# (c) The actual middle-chart tensor of the finite closed construction.
# On 1/2 <= |t| <= 2 both cusp correction weights in
# cusp.elementary_restore_torus('mid', ...) vanish.  Its exact holomorphic
# coordinates are
#   z = 2*pi*x1 - i*log(t)*x4 - i*t*phi2 - t*phi3,
#   w = 2*t*x4 + (1-t**2)*phi2 + i*(1+t**2)*phi3,
#   third coordinate = t.
# Use a smooth local log branch; its derivatives below are branch independent.
a, b, x1, x4, phi2, phi3 = sp.symbols('a b x1 x4 phi2 phi3', real=True)
u = (a, b, x1, x4, phi2, phi3)
r2 = a*a+b*b
ell, theta = sp.symbols('ell theta', real=True)
log_radius, angle = sp.log(r2)/2, sp.atan2(b,a)
middle_h = sp.Matrix([
    2*sp.pi*x1+angle*x4+b*phi2-a*phi3,
    -log_radius*x4-a*phi2-b*phi3,
    2*a*x4+(1-a*a+b*b)*phi2-2*a*b*phi3,
    2*b*x4-2*a*b*phi2+(1+a*a-b*b)*phi3,
    a, b,
])
# Keep ell and theta as point values while differentiating them by their
# actual chain rules.  This prevents an unnecessary expansion of log/atan2
# in the 216 tensor coefficients.
M = sp.Matrix([
    [-b*x4/r2-phi3, a*x4/r2+phi2, 2*sp.pi, theta, b, -a],
    [-a*x4/r2-phi2, -b*x4/r2-phi3, 0, -ell, -a, -b],
    [2*x4-2*a*phi2-2*b*phi3, 2*b*phi2-2*a*phi3, 0, 2*a, 1-a*a+b*b, -2*a*b],
    [-2*b*phi2+2*a*phi3, 2*x4-2*a*phi2-2*b*phi3, 0, 2*b, -2*a*b, 1+a*a-b*b],
    [1, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0],
])
assert all(sp.simplify(v) == 0 for v in middle_h.jacobian(u)-M.subs({ell:log_radius,theta:angle}))
print('(c1) exact middle-coordinate Jacobian: PASS')

middle_det = 2*sp.pi*(1+r2)*(ell*(r2-1)+2*r2)
A = M[:4, 2:]
assert sp.factor(A.det()-middle_det) == 0
# The factor (r^2-1) log r is nonnegative for every r>0, so the
# determinant is strictly positive on the entire punctured middle chart.
print('(c2) det Dh =', middle_det, '; positive for |t|>0: PASS')

def middle_derivative(expr, variable):
    """Partial derivative at a point, including log|t| and Arg(t)."""
    ans = sp.diff(expr, variable)
    if variable == a:
        ans += a*sp.diff(expr, ell)/r2-b*sp.diff(expr, theta)/r2
    elif variable == b:
        ans += b*sp.diff(expr, ell)/r2+a*sp.diff(expr, theta)/r2
    return ans

# H[k][out,i] = partial_k partial_i h_out is the actual second jet.
H = [M.applyfunc(lambda e: middle_derivative(e, variable)) for variable in u]
assert all(sp.cancel(H[k][out,i]-H[i][out,k]) == 0
           for k in range(6) for out in range(6) for i in range(6))
print('(c3) actual middle Hessian has all 216 mixed-partial symmetries: PASS')

J0 = sp.diag(*([sp.Matrix([[0,-1],[1,0]])]*3))
# M = [[B,A],[I,0]], so its inverse needs only the four-dimensional
# lattice block A; in this ordering det M = det A.
B = M[:4, :2]
Ainv = A.inv()
Minv = sp.Matrix.vstack(
    sp.Matrix.hstack(sp.zeros(2,4), sp.eye(2)),
    sp.Matrix.hstack(Ainv, -Ainv*B),
)
assert all(sp.cancel(e) == 0 for e in M*Minv-sp.eye(6))
Jmid = Minv*J0*M
assert all(sp.cancel(e) == 0 for e in Jmid*Jmid+sp.eye(6))
print('(c4) actual middle tensor J=Dh^-1 J0 Dh squares to -Id: PASS')

# Differentiate J using d(M^-1)=-M^-1(dM)M^-1.  Write
# S(v,w)=sum_(k,i) v_k w_i H[k][:,i] and T(v,w)=S(v,w)-S(w,v).
# Substitution into the full coordinate formula for N gives, componentwise,
# M N(v,w) = J0*T(Jv,w)-J0*T(Jw,v)-T(Jv,Jw)+T(v,w).
# This identity follows by distributing the four dJ terms and J0^2=-Id;
# The identity is checked against the unreduced coordinate formula at a
# non-coordinate frame below.  The next 216 checks use this specific h's H.
skew_H = [[sp.Matrix([sp.cancel(H[i][out,j]-H[j][out,i])
                      for out in range(6)]) for j in range(6)] for i in range(6)]
assert all(entry == 0 for row in skew_H for vector in row for entry in vector)

def skew_pair(left, right, skew):
    return sp.Matrix([sum(left[i]*right[j]*skew[i][j][out]
                          for i in range(6) for j in range(6)) for out in range(6)])

Nmid = []
for i in range(6):
    for j in range(6):
        ei, ej = sp.eye(6)[:,i], sp.eye(6)[:,j]
        Ji, Jj = Jmid[:,i], Jmid[:,j]
        target_N = (J0*skew_pair(Ji,ej,skew_H)
                    -J0*skew_pair(Jj,ei,skew_H)
                    -skew_pair(Ji,Jj,skew_H)+skew_pair(ei,ej,skew_H))
        Nmid.extend(Minv*target_N)
assert len(Nmid) == 216 and all(v == 0 for v in Nmid)
print('(c5) all 216 exact middle-chart Nijenhuis coefficients vanish: PASS')

# Independent direct-coordinate control at an exact middle point.  It checks
# the unreduced N formula, including its signs, against the factorized form.
point = {a:1,b:0,x1:0,x4:0,phi2:0,phi3:0,ell:0,theta:0}
Mp = M.subs(point)
Hp = [slice_.subs(point) for slice_ in H]
Jp = Jmid.subs(point)
Mip = Minv.subs(point)
assert Mp*Mip == sp.eye(6)

def direct_nijenhuis(jmat, hessian):
    derivatives = [Mip*(J0*slice_-slice_*jmat) for slice_ in hessian]
    return [sp.cancel(sum(
        jmat[p,i]*derivatives[p][out,j]
        -jmat[p,j]*derivatives[p][out,i]
        -jmat[out,p]*(derivatives[i][p,j]-derivatives[j][p,i])
        for p in range(6)))
        for out in range(6) for i in range(6) for j in range(6)]

assert all(v == 0 for v in direct_nijenhuis(Jp,Hp))
print('(c6) unreduced 216-component N formula at t=1 agrees: PASS')

# A frame perturbation with an antisymmetric second jet must be detected.
# Alter only partial_a M[0,x1], leaving the point value M unchanged.
Hbad = [slice_.copy() for slice_ in Hp]
Hbad[0][0,2] += 1
bad_values = direct_nijenhuis(Jp,Hbad)
bad_skew = [[sp.Matrix([Hbad[i][out,j]-Hbad[j][out,i]
                        for out in range(6)]) for j in range(6)] for i in range(6)]
for i in range(6):
    for j in range(6):
        ei, ej = sp.eye(6)[:,i], sp.eye(6)[:,j]
        Ji, Jj = Jp[:,i], Jp[:,j]
        factored = Mip*(J0*skew_pair(Ji,ej,bad_skew)
                        -J0*skew_pair(Jj,ei,bad_skew)
                        -skew_pair(Ji,Jj,bad_skew)+skew_pair(ei,ej,bad_skew))
        assert all(sp.cancel(bad_values[out*36+i*6+j]-factored[out]) == 0
                   for out in range(6))
bad_first = next(((out,i,j,value) for index,value in enumerate(bad_values)
                  if value != 0 for out,i,j in
                  [(index//36,(index//6)%6,index%6)]), None)
assert bad_first is not None
print('(c7) factorization matches all 216 direct negative-control components;'
      ' first nonzero:', bad_first, ': PASS')

# Prove the factorization used in (c5) as a 216-component polynomial identity,
# even when the frame is not a coordinate Jacobian.  A constant source change
# normalizes M to Id at the point, so J=J0 there; the generic frame derivative
# G[k][out,i] remains completely arbitrary and need not be symmetric.
G = [sp.Matrix(6,6,lambda out,i: sp.Symbol(f'G_{k}_{out}_{i}'))
     for k in range(6)]
generic_dJ = [J0*G[k]-G[k]*J0 for k in range(6)]
generic_skew = [[G[i][:,j]-G[j][:,i] for j in range(6)] for i in range(6)]
for i in range(6):
    for j in range(6):
        ei, ej = sp.eye(6)[:,i], sp.eye(6)[:,j]
        Ji, Jj = J0[:,i], J0[:,j]
        factored = (J0*skew_pair(Ji,ej,generic_skew)
                    -J0*skew_pair(Jj,ei,generic_skew)
                    -skew_pair(Ji,Jj,generic_skew)
                    +skew_pair(ei,ej,generic_skew))
        for out in range(6):
            direct = sum(
                J0[p,i]*generic_dJ[p][out,j]
                -J0[p,j]*generic_dJ[p][out,i]
                -J0[out,p]*(generic_dJ[i][p,j]-generic_dJ[j][p,i])
                for p in range(6))
            assert sp.expand(direct-factored[out]) == 0
print('(c8) all 216 factorization identities for arbitrary frame jets: PASS')
print('ALL PASS')
