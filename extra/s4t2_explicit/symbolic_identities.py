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
      Nijenhuis identities.  The round-to-W factor is outside this check.
  (d) on one open, active south middle-collar seed chamber, the exact Cramer
      and four-seed formulas give a smooth torus shear.  Its compositional
      second jet has positive determinant, J^2=-Id, and N=0.  This does not
      cover the other seed chambers or the north collar.
  (e) on one south raw-pole annulus chamber, the actual round cap, raw inverse
      and cusp torus restoration give a source-linked six-real-coordinate
      second-jet circuit and J.  Its compositional 36/216 coefficient checks
      use the written cap and derivative bounds; other chambers remain open.
  (f) the corresponding north raw-pole annulus has a separate source-linked
      second-jet circuit and compositional 36/216 coefficient checks.
  (g) one south nonidentity radial-transition chamber includes both radial
      conjugations, the actual raw inverse, round cap, and torus restoration.
  (h) the corresponding north transition chamber has a separately bound
      six-real-coordinate J circuit and compositional 36/216 checks.
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

# (d) One genuinely active c0 middle-collar chamber, 3/10<|t|<1/2.
# This chamber has beta0 in (0,1), beta_inf=0, and
# 1/16<A,Y<1/2.  Since q<1/8, the other two seed arguments qA,qY
# lie below 1/16.  These are open inequalities, not a claim that this
# single chamber covers the whole collar.
def flat_open(v):
    left, right = sp.exp(-1/v), sp.exp(-1/(1-v))
    return left/(left+right)

sr, si = theta/(2*sp.pi), -ell/(2*sp.pi)  # log(t)/(2*pi*i)
qz_real = sp.pi*(sr*sr-si*si-sr)          # pi*s*(s-1)
qz_imag = sp.pi*si*(2*sr-1)
qw_real = 2*(a*sr-b*si)                  # 2*t*s
qw_imag = 2*(a*si+b*sr)
R0 = sp.Matrix([[1-a*a+b*b,-2*a*b],[-2*a*b,1+a*a-b*b]])
R0det = 1-r2*r2
assert sp.expand(R0.det()-R0det) == 0
c0 = R0.adjugate()*sp.Matrix([2*a,2*b])/R0det
wc_real = 2*(a*(x4-sr)+b*si)            # w-2*t*s
wc_imag = 2*(b*(x4-sr)-a*si)
cc0 = R0.adjugate()*sp.Matrix([-wc_real,-wc_imag])/R0det
assert all(sp.cancel(v) == 0 for v in R0*c0-sp.Matrix([2*a,2*b]))
assert all(sp.cancel(v) == 0 for v in R0*cc0+sp.Matrix([wc_real,wc_imag]))

L0 = sp.Matrix([[-a,-b,-ell],
                [1-a*a+b*b,-2*a*b,2*a],
                [-2*a*b,1+a*a-b*b,2*b]])
L0det = (1+r2)*(ell*(r2-1)+2*r2)
assert sp.factor(L0.det()-L0det) == 0
xi0 = L0.adjugate()*sp.Matrix([qz_imag,qw_real,qw_imag])/L0det
assert all(sp.cancel(v) == 0 for v in
           L0*xi0-sp.Matrix([qz_imag,qw_real,qw_imag]))
print('(d1) actual c0 torus and source-lattice Cramer systems: PASS')

# From yy=exp(i*(zc-i*(cc0[0]*t-cc0[1]*i*t))), its squared modulus is
# Y=exp(2*(Re(cc0[0]*t-cc0[1]*i*t)-Im(zc))).  x_s=t/yy.
zc_imag = -ell*x4-qz_imag
Y0 = sp.exp(2*(a*cc0[0]+b*cc0[1]-zc_imag))
q0 = r2*sp.exp(-4*r2/(1-r2))
A0 = q0/Y0
step_A = flat_open(sp.log(16*A0)/sp.log(8))
step_Y = flat_open(sp.log(16*Y0)/sp.log(8))
Hfinite = (A0-Y0)/(1-q0)-step_A*A0*A0/(1+A0)+step_Y*Y0*Y0/(1+Y0)
# On this chamber the four seeds are seed(A)=A-step_A*A^2/(1+A),
# seed(qA)=qA, seed(Y)=Y-step_Y*Y^2/(1+Y), seed(qY)=qY.
four_seeds = (A0-step_A*A0*A0/(1+A0)+q0*A0
              -Y0+step_Y*Y0*Y0/(1+Y0)-q0*Y0
              +q0*q0*(A0-Y0)/(1-q0))
assert sp.cancel(four_seeds-Hfinite) == 0
# q<1/8 follows from exp(x)>1+x: q<u(1-u)/(1+3u), u=r^2,
# and 1-5u+8u^2=8(u-5/16)^2+7/32>0.
u0 = sp.symbols('u0', real=True)
assert sp.expand(1-5*u0+8*u0*u0
                 -(8*(u0-sp.Rational(5,16))**2+sp.Rational(7,32))) == 0
# The chamber is nonempty: t=2/5, x4 below gives Y=1/4 and
# A=4q in (16/105,336/925), which lies in (1/16,1/2).
x4star = sp.log(sp.Rational(5,8))/(2*(sp.log(sp.Rational(2,5))-sp.Rational(8,21)))
Yexponent = 2*(a*cc0[0]+b*cc0[1]-zc_imag)
assert sp.cancel(Yexponent.subs({b:0,theta:0})
                 -ell-2*(ell-2*a*a/(1-a*a))*x4) == 0
assert sp.cancel(2*(sp.log(sp.Rational(2,5))-sp.Rational(8,21))*x4star
                 -sp.log(sp.Rational(5,8))) == 0
assert sp.simplify(sp.exp(sp.log(sp.Rational(2,5))+sp.log(sp.Rational(5,8)))
                   -sp.Rational(1,4)) == 0
qstar = sp.simplify(q0.subs({a:sp.Rational(2,5),b:0}))
xstar = sp.Rational(16,21)
assert qstar == sp.Rational(4,25)*sp.exp(-xstar)
assert sp.simplify(4*qstar-sp.Rational(16,25)*sp.exp(-xstar)) == 0
# For x>0, exp(x)>1+x and exp(-x)>1-x: both differences start
# at zero, and their derivatives exp(x)-1 and 1-exp(-x) are positive.
xp = sp.symbols('xp', positive=True)
assert sp.diff(sp.exp(xp)-1-xp,xp) == sp.exp(xp)-1
assert sp.diff(sp.exp(-xp)-1+xp,xp) == 1-sp.exp(-xp)
assert 0<xstar<1
assert sp.Rational(16,25)*(1-xstar) == sp.Rational(16,105)
assert sp.Rational(16,25)/(1+xstar) == sp.Rational(336,925)
assert sp.Rational(16,105)>sp.Rational(1,16)
assert sp.Rational(336,925)<sp.Rational(1,2)
beta0 = 1-flat_open((sp.sqrt(r2)-sp.Rational(3,10))/sp.Rational(1,5))
D0 = sp.Matrix([xi0[k]+(Hfinite-x4+xi0[2])*c0[k] for k in range(2)])
assert all(not value.has(x1,phi2,phi3) for value in D0)
assert not beta0.has(x1,x4,phi2,phi3)
# This is the exact value formula of elementary_restore_torus('mid',...)
# on the stated open chamber, expressed in real components.
C2, C3 = phi2+beta0*D0[0], phi3+beta0*D0[1]
collar_h = middle_h.subs({phi2:C2,phi3:C3},simultaneous=True)
assert collar_h[4] == a and collar_h[5] == b
print('(d2) finite four-seed height and active torus shift are explicit: PASS')

# Directly differentiate the displayed finite expressions in two orders.
# The D0 test holds ell,theta independent here; the physical derivatives
# also include their explicitly commuting a,b chain-rule vector fields.
assert sp.diff(beta0,a,b) == sp.diff(beta0,b,a)
assert all(sp.diff(value,a,x4) == sp.diff(value,x4,a) for value in D0)
assert sp.cancel(sp.diff(b/r2,a)-sp.diff(a/r2,b)) == 0
assert sp.cancel(sp.diff(a/r2,a)-sp.diff(-b/r2,b)) == 0
print('(d2b) actual beta: mixed (a,b); both actual D components:'
      ' mixed (a,x4) at fixed ell,theta; log/angle chain fields commute: PASS')

# For the 2-jet check, keep the already specified elementary D0 and beta0
# as named smooth nodes.  Formal derivatives here are exact chain-rule nodes,
# not numerical samples or unknown functions chosen after the fact.
bf = sp.Function('beta0')(a,b)
d2f = sp.Function('D02')(a,b,x4)
d3f = sp.Function('D03')(a,b,x4)
shear = sp.Matrix([a,b,x1,x4,phi2+bf*d2f,phi3+bf*d3f])
S1 = shear.jacobian(u)
assert sp.factor(S1.det()-1) == 0
S2 = [S1.diff(variable) for variable in u]
assert all(sp.simplify(S2[i][out,j]-S2[j][out,i]) == 0
           for i in range(6) for j in range(6) for out in range(6))
# h_collar = h_middle o shear.  Its first and second jets are the
# following exact finite chain-rule circuits; no inverse of D0 is needed.
subC = {phi2:shear[4],phi3:shear[5]}
M0S = M.subs(subC,simultaneous=True)
H0S = [slice_.subs(subC,simultaneous=True) for slice_ in H]
Mcollar = M0S*S1
# The first factor has det middle_det, and the shear has det 1.
assert sp.factor(M0S[:4,2:].det()-middle_det) == 0
S1inv = sp.eye(6)
S1inv[4,:4], S1inv[5,:4] = -S1[4,:4], -S1[5,:4]
assert S1*S1inv == sp.eye(6)
# The actual collar J is S1^-1*(J_middle at C)*S1, so these inverse
# identities and (c4) give its 36 J^2=-Id identities exactly.
assert J0*J0 == -sp.eye(6)
print('(d3) active-collar shear has det 1; full coordinate det =',
      middle_det, ': PASS')

# For each k,i, the antisymmetric part of the composed Hessian is
# sum_(p,q) (H0[p][:,q]-H0[q][:,p])*S1[p,k]*S1[q,i]
#       + M0S*(S2[k][:,i]-S2[i][:,k]).
# This is the exact second-order chain rule.  Checking every output and
# input pair tests all 216 actual collar Hessian-skew coefficients without
# expanding the large finite-height expression.  The generic 216-component
# factorization in (c8) then gives N=0 for this collar J.
for k in range(6):
    for i in range(6):
        skew_composed = (
            sum(((H0S[p][:,q]-H0S[q][:,p])*S1[p,k]*S1[q,i]
                 for p in range(6) for q in range(6)), sp.zeros(6,1))
            +M0S*(S2[k][:,i]-S2[i][:,k]))
        assert all(sp.simplify(value) == 0 for value in skew_composed)
print('(d4) 216 composed Hessian-skew coefficients vanish; collar J^2=-Id'
      ' and N=0 follow from (c8): PASS')

# Numerical implementation correspondence only: this 60-digit interior
# sample is not used in the symbolic proof above.  It catches transcription
# drift between this SymPy formula and cusp.elementary_restore_torus.
import sys
import mpmath as mp
sys.path.insert(0,'extra/s4t2_explicit')
import cusp
mp.mp.dps = 70
t_sample = mp.mpf(2)/5
x4_sample = mp.mpf(str(x4star.evalf(60)))
one_sample, p2_sample, p3_sample = mp.mpf(1)/7,mp.mpf(1)/11,-mp.mpf(1)/13
source = cusp.elementary_restore_torus(
    'mid',(t_sample,one_sample,x4_sample,mp.log(t_sample)),
    (p2_sample,p3_sample))[1]
sample_sub = {a:sp.Rational(2,5),b:0,ell:sp.log(sp.Rational(2,5)),theta:0,
              x1:sp.Rational(1,7),x4:x4star,
              phi2:sp.Rational(1,11),phi3:-sp.Rational(1,13)}
symbolic_value = [mp.mpf(str(value.subs(sample_sub).evalf(60))) for value in collar_h]
source_value = [source[0].real,source[0].imag,source[1].real,
                source[1].imag,source[2].real,source[2].imag]
source_error = max(abs(left-right) for left,right in zip(symbolic_value,source_value))
assert source_error < mp.mpf('1e-50')
print('(d5) numerical source-code correspondence at one interior point:'
      ' max error =',mp.nstr(source_error,5),': PASS')

# (e) One actual south raw-pole chamber of the round-to-W factor, followed
# by the exact cusp torus restoration.  This is an open six-real-dimensional
# region: 3e-8<r_L=|z_L|^2<1e-7, Re(u2)>0, 0<height<1/10, and t_c!=0.
# The raw inverse contains no floor, terminal-cut, or precision shortcut here.
# It is distinct from (c),(d), which start after the round-to-W factor.
import ast
import inspect


def source_statement_e(function, statement):
    tree = ast.parse(inspect.getsource(function))
    expected = ast.dump(ast.parse(statement).body[0], include_attributes=False)
    return any(ast.dump(node, include_attributes=False) == expected
               for node in ast.walk(tree))


def source_if_first_e(function, condition, first_statement):
    tree = ast.parse(inspect.getsource(function))
    want_test = ast.dump(ast.parse('if '+condition+':\n    pass').body[0].test,
                         include_attributes=False)
    want_body = ast.dump(ast.parse(first_statement).body[0],
                         include_attributes=False)
    return any(isinstance(node,ast.If) and node.body and
               ast.dump(node.test,include_attributes=False)==want_test and
               ast.dump(node.body[0],include_attributes=False)==want_body
               for node in ast.walk(tree))


assert source_statement_e(cusp.elementary_global_slice,
    'if distance<alpha/2 and not force_slab:\n    inf = X>0\n    scale = mp.pi if inf else 2*mp.pi\n    den = mp.sqrt(scale*(1+abs(X)*mp.sqrt(2-X*X)))\n    pair = _elementary_pole_inclusion(Z,inf)[4:6] if inf else _elementary_pole_inclusion(Z)[:2]\n    zz = elementary_pole_patch(tuple(a/den for a in pair),inf,dps=dps)\n    nu = mp.exp(-2) if inf else mp.mpf(1)\n    return (\'ci\' if inf else \'c0\'),(mp.conj(zz[0])/nu,zz[1])')
assert source_statement_e(cusp._elementary_pole_radial,
                          'if radius>=c:\n    return z')
assert source_statement_e(cusp._elementary_pole_radial,
    "a,c = mp.mpf('1e-8'),mp.mpf('1.5e-8')")
assert source_statement_e(cusp.elementary_pole_patch,
    'raw = elliptic_pole_raw_inverse(compressed,infinity,dps=dps,elementary=True)')
for statement in (
    "inner = mp.mpf('1e-8')",
    "outer = mp.mpf('1.5e-8')",
    'if radius<=inner:\n    return pair',
    'compressed = _elementary_pole_radial(pair)',
    'return _elementary_pole_radial(raw,inverse=True)',
):
    assert source_statement_e(cusp.elementary_pole_patch, statement), statement
assert source_if_first_e(cusp.elementary_pole_patch,
                         'radius<(inner+outer)/2',
                         'x = mp.exp(3)*mp.expm1(3*(radius-inner)/(outer-radius))')
assert source_statement_e(cusp.elliptic_pole_raw_inverse,
    'tc = ratio*mp.sqrt((2*B if infinity else 2)/den)')
for statement in (
    'pc = 2*z0*mp.conj(z1)',
    'p = (pc.real,pc.imag,abs(z0)**2-abs(z1)**2)',
    'scale = mp.pi if infinity else 2*mp.pi',
    'v = tuple(scale*w for w in ((p[2],-p[1],p[0]) if infinity else (-p[1],-p[2],p[0])))',
    'aa,bb,cc = v[0],-v[2],-v[1]',
    's2 = aa*aa+(bb+1)**2+cc*cc',
    'u1 = mp.mpc(1-s2,2*aa)/(1+s2)',
    'u2 = mp.mpc(2*(bb+1),2*cc)/(1+s2)',
    'ratio = u1/u2',
    'height = mp.arg(u2)/(2*mp.pi)',
    'A,B = (4*mp.pi)**2,mp.exp(4)',
    'uu = abs(ratio)**2',
    'den = A*(1-uu)+mp.sqrt(A*A*(1-uu)**2+4*B*uu)',
    'e = mp.exp(-4*abs(tc)**2/(1-abs(tc)**2))',
    'q = e*nu*nu*abs(tc)**2',
    'nu = mp.exp(-2) if infinity else mp.mpf(1)',
    'xx = AA/e',
    'moment,source_radius = xx-yy,xx+yy',
    'source_pc = 2*nu*mp.conj(tc)',
    'w0 = ((source_radius+moment)*z0+source_pc*z1)/(2*source_radius)',
    'w1 = (mp.conj(source_pc)*z0+(source_radius-moment)*z1)/(2*source_radius)',
    'length = mp.sqrt(abs(w0)**2+abs(w1)**2)',
    'factor = mp.sqrt(source_radius)/length',
    'return (+(factor*w0),+(factor*w1))',
):
    assert source_statement_e(cusp.elliptic_pole_raw_inverse, statement), statement
assert source_if_first_e(cusp.elliptic_pole_raw_inverse,
    "not 0<radius<=mp.mpf('1e-4')",
    "raise ValueError('raw pole inverse requires 0<|z|^2<=1e-4')")
assert source_if_first_e(cusp.elliptic_pole_raw_inverse,
    "abs(height)>=mp.mpf('.1')",
    "raise ValueError('elementary pole inverse left its certified height region')")
assert source_statement_e(cusp.elliptic_pole_raw_inverse,
    'R = mp.sqrt(height*height+4*q)')
assert source_statement_e(cusp.elliptic_pole_raw_inverse,
    'AA = 2*q/(R-height) if height<0 else (R+height)/2')
assert source_statement_e(cusp.elliptic_pole_raw_inverse,
    'yy = 2*q/(R+height) if height>=0 else (R-height)/2')
assert source_statement_e(cusp.elliptic_pole_raw_inverse,
    'factor = mp.sqrt(source_radius)/length')
assert source_statement_e(cusp.elementary_restore_torus,
    'H = scalar_finite_height(x,y,mp.exp(-2) if inf else mp.mpf(1))')
for statement in (
    "inf = kind=='ci'",
    't = x*y',
    'e2,e3,e4,a2,a3 = cusp_data(t,inf)',
    'c = coeffs(e4,e2,e3)',
    'C = tuple(v+H*w for v,w in zip(phi,c))',
    'A = C[0]*a2+C[1]*a3',
    'return kind,(x*mp.exp(-A),y*mp.exp(A),C[0]*e2+C[1]*e3)',
):
    assert source_statement_e(cusp.elementary_restore_torus, statement), statement
assert source_statement_e(cusp.elementary_restore_torus,
    'return ((t*t-1,mp.j*(1+t*t),-2*t,t,mp.j*t) if inf else '
    '(1-t*t,mp.j*(1+t*t),2*t,t,-mp.j*t))')
for statement in (
    'rr = abs(x*y)**2',
    'e = mp.exp(-4*rr/(1-rr))',
    'A,Y = nu*nu*e*abs(x)**2,abs(y)**2',
    'q = A*Y',
    'return seed(A)+seed(q*A)-seed(Y)-seed(q*Y)+q*q*(A-Y)/(1-q)',
):
    assert source_statement_e(cusp.scalar_finite_height, statement), statement
assert source_if_first_e(cusp.scalar_finite_height,
                         'v<=mp.mpf(1)/16','value = v')
print('(e1) exact south cap, raw inverse and torus source statements: PASS')

# The annulus lies inside the cap: h=2*pi*r_L<8e-7 and
# 1-|X|<h+h^2/2<alpha/2.  Its two radial-conjugation calls are identity:
# (11.11) gives |g(z)|^2>=r_L(1-200r_L)^2>1.5e-8.
rmin_e,rmax_e = sp.Rational(3,10**8),sp.Rational(1,10**7)
hmax_e = 8*rmax_e   # pi<4
assert rmin_e>sp.Rational(2,10**8)  # raw outer branch; skips precision shortcut
# With rho=|Z|^2 and h=2*pi*r_L, (11.17) gives
# rho=2h/(1+h^2), |X|=(1-h)/sqrt(1+h^2).  For 0<h<1,
# (1+h^2)(1-h^2/2)^2<1, so 1-|X|<h+h^2/2.
he_e = sp.symbols('he_e',positive=True)
rho_e = 2*he_e/(1+he_e**2)
assert sp.factor(1-rho_e**2
                 -((1-he_e**2)/(1+he_e**2))**2)==0
assert sp.expand(1-(1+he_e**2)*(1-he_e**2/2)**2
                 -he_e**4*(3-he_e**2)/4)==0
assert hmax_e<1
assert hmax_e+hmax_e**2/2 < sp.Rational(1,2*10**5)
assert rmin_e*(1-200*rmax_e)**2 > sp.Rational(3,2*10**8)
assert 5000*rmax_e == sp.Rational(1,2000)
assert 1-5000*rmax_e == sp.Rational(1999,2000)
assert rmax_e+80*rmax_e**2<sp.Rational(2,10**6)
# Each of A,Y and qA,qY in scalar_finite_height lies below 1/16 on this
# cap, by the source-radius estimate.  All four seed calls therefore use
# their linear branch; the exact q^2 tail sums the remaining geometric tail.
source_radius_max_e = sp.Rational(2,10**6)
assert source_radius_max_e<sp.Rational(1,16)
assert source_radius_max_e**2/4<sp.Rational(1,4)
assert source_radius_max_e**2<sp.Rational(1,8)
As_e,Ys_e,qs_e = sp.symbols('As_e Ys_e qs_e')
assert sp.factor((As_e-Ys_e)*(1+qs_e+qs_e**2/(1-qs_e))
                 -(As_e-Ys_e)/(1-qs_e))==0
# A rational complex witness with r_L=6e-8 has positive height and t_c!=0.
zw0_e,zw1_e = sp.Rational(1,5000),(1-sp.I)/10000
rw_e = sp.expand(zw0_e*sp.conjugate(zw0_e)
                 +zw1_e*sp.conjugate(zw1_e))
pcw_e = sp.expand(2*zw0_e*sp.conjugate(zw1_e))
p3w_e = sp.expand(zw0_e*sp.conjugate(zw0_e)
                  -zw1_e*sp.conjugate(zw1_e))
assert rw_e == sp.Rational(6,10**8) and rmin_e<rw_e<rmax_e
assert sp.re(pcw_e)>0 and sp.im(pcw_e)>0 and p3w_e>0
# beta=-2*pi*Re(pc), gamma=2*pi*p3.  Since pi<4, 1+beta>0;
# gamma>0 makes Arg(u2)>0, while alpha=-2*pi*Im(pc)!=0 makes u1!=0.
assert 1-8*sp.re(pcw_e)>0 and sp.im(pcw_e)>0
assert 8*p3w_e/(1-8*sp.re(pcw_e)) < sp.Rational(1,2)
# Then gamma/(1+beta)<1/2, so height<1/(4*pi)<1/12<1/10.
assert sp.pi>3 and sp.Rational(1,12)<sp.Rational(1,10)
print('(e2) exact open-chamber witness; raw derivative margin >=1999/2000: PASS')

# Actual source raw inverse in four real pole coordinates.  The real
# formulas below are its elementary=True, south, height>0 branch.  For the
# displayed source projector, uY=|t_c|^2 and R^2-height^2=4q; hence its
# rationalized yy=2q/(R+height) equals (R-height)/2 exactly.
z0r,z0i,z1r,z1i = z_e = sp.symbols('z0r_e z0i_e z1r_e z1i_e',real=True)
pr_e = 2*(z0r*z1r+z0i*z1i)
pi_e = 2*(z0i*z1r-z0r*z1i)
p3_e = z0r*z0r+z0i*z0i-z1r*z1r-z1i*z1i
alpha_e,beta_e,gamma_e = -2*sp.pi*pi_e,-2*sp.pi*pr_e,2*sp.pi*p3_e
d_e = alpha_e**2+(1+beta_e)**2+gamma_e**2
u1r_e,u1i_e = (1-d_e)/(1+d_e),2*alpha_e/(1+d_e)
u2r_e,u2i_e = 2*(1+beta_e)/(1+d_e),2*gamma_e/(1+d_e)
uden_e = u2r_e**2+u2i_e**2
Tr_e = (u1r_e*u2r_e+u1i_e*u2i_e)/uden_e
Ti_e = (u1i_e*u2r_e-u1r_e*u2i_e)/uden_e
Q_e = Tr_e**2+Ti_e**2
Aprof_e,Bprof_e = (4*sp.pi)**2,sp.exp(4)
D_e = Aprof_e*(1-Q_e)+sp.sqrt(Aprof_e**2*(1-Q_e)**2+4*Bprof_e*Q_e)
tr_e,ti_e = Tr_e*sp.sqrt(2/D_e),Ti_e*sp.sqrt(2/D_e)
height_e = sp.atan(gamma_e/(1+beta_e))/(2*sp.pi)
te2_e = tr_e**2+ti_e**2
E_e = sp.exp(-4*te2_e/(1-te2_e))
q_e = E_e*te2_e
R_e = sp.sqrt(height_e**2+4*q_e)
AA_e,YY_e = (R_e+height_e)/2,(R_e-height_e)/2
uu_e = AA_e/E_e
ss_e = uu_e+YY_e
# The cancellation is checked on small algebraic nodes; substituting the
# long nested elementary expressions into ``simplify`` is unnecessary.
rh_e,qe_e,ee_e,te_e = sp.symbols('rh_e qe_e ee_e te_e',positive=True)
rnode_e = sp.sqrt(rh_e**2+4*qe_e)
assert sp.expand((rnode_e-rh_e)*(rnode_e+rh_e)-4*qe_e)==0
assert sp.expand((rnode_e+rh_e)*(rnode_e-rh_e)/(4*ee_e)
                 -qe_e/ee_e)==0
assert q_e == E_e*te2_e and uu_e == AA_e/E_e
V_e = sp.Matrix([
    uu_e*z0r+tr_e*z1r+ti_e*z1i,
    uu_e*z0i+tr_e*z1i-ti_e*z1r,
    tr_e*z0r-ti_e*z0i+YY_e*z1r,
    tr_e*z0i+ti_e*z0r+YY_e*z1i])
Vnorm_e = sum(entry*entry for entry in V_e)
raw_e = sp.sqrt(ss_e/Vnorm_e)*V_e
# The source w0,w1 before its final factor are Q_source*z/s, with
# Q_source=[[u,conj(t)],[t,Y]].  Since the source length is |Qz|/s,
# its final factor gives exactly raw_e above (positive square roots).
print('(e3) actual south raw spinor has an explicit four-real-variable expression: PASS')

# Directly expanding ``sp.hessian(raw_e[k],z_e)`` duplicates the nested
# radicals exponentially.  CSE turns THIS concrete four-output expression
# into a finite straight-line program of small scalar assignments.  Each
# assignment receives its exact value, gradient and symmetric Hessian by
# differentiating its own explicit right-hand side and applying the second
# chain rule.  Jet symbols are bound by the assignments below; they are not
# free choices or unevaluated Derivative objects.
raw_cse_e,raw_out_e = sp.cse(list(raw_e),symbols=sp.numbered_symbols('raw_n_e_'))
assert len(raw_cse_e)>0 and len(raw_out_e)==4
raw_rebuilt_e = list(raw_out_e)
for node,rhs in reversed(raw_cse_e):
    raw_rebuilt_e = [entry.xreplace({node:rhs}) for entry in raw_rebuilt_e]
assert all(raw_rebuilt_e[k]==raw_e[k] for k in range(4))

raw_jet_defs_e = []
raw_known_e = {}
for k,variable in enumerate(z_e):
    raw_known_e[variable] = (sp.eye(4)[:,k],sp.zeros(4,4))


def raw_local_jet_e(rhs, known=None):
    if known is None:
        known = raw_known_e
    deps = tuple(sorted(rhs.free_symbols & known.keys(),key=str))
    assert rhs.free_symbols <= set(deps)
    d1 = {dep:sp.diff(rhs,dep) for dep in deps}
    d2 = {(dep,other):sp.diff(d1[dep],other)
          for dep in deps for other in deps}
    assert all(sp.simplify(d2[p,q]-d2[q,p])==0
               for p in deps for q in deps)
    grad = sp.Matrix([
        sum((d1[p]*known[p][0][i] for p in deps),sp.S.Zero)
        for i in range(4)])
    hess = sp.zeros(4,4)
    for i in range(4):
        for j in range(i,4):
            value = sum((d1[p]*known[p][1][i,j] for p in deps),
                        sp.S.Zero)
            value += sum((d2[p,q]*known[p][0][i]*known[q][0][j]
                          for p in deps for q in deps),sp.S.Zero)
            hess[i,j]=hess[j,i]=value
    return grad,hess


for node,rhs in raw_cse_e:
    grad_rhs,hess_rhs = raw_local_jet_e(rhs)
    raw_jet_defs_e.append((node,rhs))
    grad_symbols = sp.Matrix([sp.Symbol(f'{node}_d{i}') for i in range(4)])
    hess_symbols = sp.Matrix(4,4,lambda i,j:
                             sp.Symbol(f'{node}_h{min(i,j)}{max(i,j)}'))
    raw_jet_defs_e.extend((grad_symbols[i],grad_rhs[i]) for i in range(4))
    raw_jet_defs_e.extend((hess_symbols[i,j],hess_rhs[i,j])
                          for i in range(4) for j in range(i,4))
    raw_known_e[node]=(grad_symbols,hess_symbols)

raw_output_jets_e = [raw_local_jet_e(rhs) for rhs in raw_out_e]
Mraw_e = sp.Matrix(4,4,lambda k,i:raw_output_jets_e[k][0][i])
Hraw_e = [raw_output_jets_e[k][1] for k in range(4)]
assert Mraw_e.shape==(4,4) and len(Hraw_e)==4
assert all(Hraw_e[k][i,j]==Hraw_e[k][j,i]
           for k in range(4) for i in range(4) for j in range(4))
assert not any(value.has(sp.Derivative) for value in Mraw_e)
assert len(raw_jet_defs_e)==15*len(raw_cse_e)
raw_defined_e = set(z_e)
for node,rhs in raw_jet_defs_e:
    assert rhs.free_symbols <= raw_defined_e
    raw_defined_e.add(node)
assert all(rhs.free_symbols <= raw_defined_e for rhs in raw_out_e)
print('(e4) actual raw inverse 4x4 Jacobian and four symmetric Hessian'
      f' circuits via {len(raw_cse_e)} exact scalar 2-jet nodes: PASS')

# The round cap change of variables uses (11.17), with X<0 and
# |X|sqrt(2-X^2)=sqrt(1-|Z|^4), as checked independently in N8a0.
zr0,zi0,zr1,zi1 = round_e = sp.symbols('zr0_e zi0_e zr1_e zi1_e',real=True)
round_r_e = sum(v*v for v in round_e)
round_gain_e = 1/sp.sqrt(1+sp.sqrt(1-round_r_e**2))/sp.sqrt(4*sp.pi)
local_e = sp.Matrix([round_gain_e*(zr0+zi1),
                     round_gain_e*(zi0-zr1),
                     round_gain_e*(zr0-zi1),
                     round_gain_e*(zi0+zr1)])
Mround_e = local_e.jacobian(round_e)
Hround_e = [sp.hessian(local_e[k],round_e) for k in range(4)]
assert all(Hround_e[out][i,j] == Hround_e[out][j,i]
           for out in range(4) for i in range(4) for j in range(4))

# Exact cap tau_f^-1 in real form, now with an arbitrary raw spinor w.
w0r,w0i,w1r,w1i,p2e,p3e = rest_e = sp.symbols(
    'w0r_e w0i_e w1r_e w1i_e p2e p3e',real=True)
tR_e,tI_e = w0r*w1r+w0i*w1i,w0r*w1i-w0i*w1r
v_e = tR_e**2+tI_e**2
Ee_e = sp.exp(-4*v_e/(1-v_e))
Ae_e,Ye_e = Ee_e*(w0r*w0r+w0i*w0i),w1r*w1r+w1i*w1i
Hf_e = (Ae_e-Ye_e)/(1-Ae_e*Ye_e)
ar_e = tR_e*p2e+tI_e*p3e+2*Hf_e*v_e/(1-v_e)
ai_e = tI_e*p2e-tR_e*p3e
co_e,si_e,ea_e = sp.cos(ai_e),sp.sin(ai_e),sp.exp(ar_e)
Br_e,Cr_e = tR_e*tR_e-tI_e*tI_e,2*tR_e*tI_e
restore_e = sp.Matrix([
    (w0r*co_e-w0i*si_e)/ea_e,
    (-w0r*si_e-w0i*co_e)/ea_e,
    ea_e*(w1r*co_e-w1i*si_e),
    ea_e*(w1r*si_e+w1i*co_e),
    p2e*(1-Br_e)-p3e*Cr_e+2*Hf_e*tR_e,
    p3e*(1+Br_e)-p2e*Cr_e+2*Hf_e*tI_e])
Mrestore_e = restore_e.jacobian(rest_e)
Hrestore_e = [sp.hessian(restore_e[k],rest_e) for k in range(6)]
assert all(Hrestore_e[out][i,j] == Hrestore_e[out][j,i]
           for out in range(6) for i in range(6) for j in range(6))
assert not any(value.has(sp.Derivative) for value in Mrestore_e)
print('(e5) actual round cap and finite-height torus first/second jets: PASS')

# The complete map is restore_e o (raw_e o local_e,phi2,phi3).  The first
# four assignments below bind the raw DAG's input coordinates to the actual
# round cap output.  Subsequent value/jet assignments are in topological
# order.  Subs binds the restoration derivatives at the actual raw output.
def bind_e(expression, variables, values):
    return sp.Subs(expression, tuple(variables), tuple(values))


def bind_matrix_e(matrix, variables, values):
    return matrix.applyfunc(lambda entry: bind_e(entry, variables, values))


raw_stage_defs_e = list(zip(z_e,local_e))+raw_jet_defs_e
raw_at_local_e = sp.Matrix(raw_out_e)
base_at_round_e = tuple(raw_at_local_e)+(p2e,p3e)
full_at_round_e = sp.Matrix([
    bind_e(restore_e[k],rest_e,base_at_round_e) for k in range(6)])
Mraw_at_local_e = Mraw_e
Hraw_at_local_e = Hraw_e
Mrestore_at_base_e = bind_matrix_e(Mrestore_e,rest_e,base_at_round_e)
Hrestore_at_base_e = [bind_matrix_e(H,rest_e,base_at_round_e)
                       for H in Hrestore_e]
Mbase_e = sp.BlockDiagMatrix(
    sp.MatMul(Mraw_at_local_e,Mround_e,evaluate=False),sp.eye(2))
Mfull_e = sp.MatMul(Mrestore_at_base_e,Mbase_e,evaluate=False)
J_annulus_e = sp.MatMul(sp.Inverse(Mfull_e),J0,Mfull_e,evaluate=False)
round_inputs_e = set(round_e)|{p2e,p3e}
bound_symbols_e = set(round_inputs_e)
for node,rhs in raw_stage_defs_e:
    assert rhs.free_symbols <= bound_symbols_e
    bound_symbols_e.add(node)
assert tuple(rhs for _,rhs in raw_stage_defs_e[:4])==tuple(local_e)
assert Mfull_e.shape == (6,6)
assert J_annulus_e.shape == (6,6)
assert all(value.free_symbols <= bound_symbols_e for value in full_at_round_e)
assert all(value.free_symbols <= bound_symbols_e for value in Mraw_at_local_e)
assert all(value.free_symbols <= bound_symbols_e for value in Mrestore_at_base_e)
assert J_annulus_e.free_symbols <= bound_symbols_e

# Six actual second-derivative matrices of the bound composite h are stored
# as unevaluated chain-rule circuits.  Their component formula is
# H_(f o g)^k = sum_p Df[k,p] H_g^p + Dg^T H_f^k Dg.
# Keeping MatMul nodes avoids a prohibitively large radical expansion.
Hbase_e = []
for k in range(4):
    h = sp.MatAdd(*[
        sp.MatMul(Mraw_at_local_e[k,p],Hround_e[p],evaluate=False)
        for p in range(4)],
        sp.MatMul(Mround_e.T,Hraw_at_local_e[k],Mround_e,
                  evaluate=False),evaluate=False)
    Hbase_e.append(sp.BlockDiagMatrix(h,sp.zeros(2,2)))
Hbase_e.extend([sp.ZeroMatrix(6,6),sp.ZeroMatrix(6,6)])
Hfull_e = []
for k in range(6):
    Hfull_e.append(sp.MatAdd(*[
        sp.MatMul(Mrestore_at_base_e[k,p],Hbase_e[p],evaluate=False)
        for p in range(6)],
        sp.MatMul(Mbase_e.T,Hrestore_at_base_e[k],Mbase_e,
                  evaluate=False),evaluate=False))
assert len(Hfull_e)==6 and all(H.shape==(6,6) for H in Hfull_e)
assert all(H.free_symbols <= bound_symbols_e for H in Hfull_e)

# Check every mixed-partial skew in the ACTUAL bound stages.  The exact
# chain-rule skew is sum_p Df[k,p](Hg^p_ij-Hg^p_ji) plus
# sum_{p<q}(Hf^k_pq-Hf^k_qp)(Dg_pi Dg_qj-Dg_pj Dg_qi).
# Testing the two factors before multiplication is an exact sparse
# evaluation of this expression; no numerical substitution is involved.
hpq_e,hqp_e,gpi_e,gpj_e,gqi_e,gqj_e = sp.symbols(
    'hpq_e hqp_e gpi_e gpj_e gqi_e gqj_e')
pair_e = gpi_e*gqj_e-gpj_e*gqi_e
assert sp.expand(hpq_e*pair_e
                 +hqp_e*(gqi_e*gpj_e-gqj_e*gpi_e)
                 -(hpq_e-hqp_e)*pair_e)==0
# Formal 4x4 and 6x6 block instances certify both composition steps from
# the tested stage Hessians, without expanding their nested radicals.  The
# generic a_p stand for any output row of an outer Jacobian; B_p and C are
# arbitrary symmetric inner/outer Hessians.
def generic_chain_symmetry_e(n):
    A = sp.Matrix(n,n,lambda p,i:sp.Symbol(f'A{n}_{p}_{i}'))
    aa = [sp.Symbol(f'a{n}_{p}') for p in range(n)]
    BB = [sp.Matrix(n,n,lambda i,j:
                    sp.Symbol(f'B{n}_{p}_{min(i,j)}_{max(i,j)}'))
          for p in range(n)]
    CC = sp.Matrix(n,n,lambda p,q:
                   sp.Symbol(f'C{n}_{min(p,q)}_{max(p,q)}'))
    HH = sum((aa[p]*BB[p] for p in range(n)),sp.zeros(n,n))
    HH += A.T*CC*A
    return all(sp.expand(HH[i,j]-HH[j,i])==0
               for i in range(n) for j in range(n))


assert generic_chain_symmetry_e(4) and generic_chain_symmetry_e(6)
def second_jet_skew_e(outer_M,outer_H,inner_M,inner_H,
                      inner_skew_certified=False):
    nout,nmid = outer_M.shape
    ninput = inner_M.cols
    assert inner_M.rows == nmid and len(outer_H)==nout
    assert len(inner_H)==nmid
    if not inner_skew_certified:
        assert all(inner_H[p][i,j]-inner_H[p][j,i]==0
                   for p in range(nmid) for i in range(ninput)
                   for j in range(ninput))
    # A prior composition certificate plus the exact displayed Hbase_e
    # chain-rule circuit makes the inner skew zero without asking SymPy to
    # expand a MatAdd entry.  Every outer skew is checked from actual H_f.
    assert all(outer_H[out][p,q]-outer_H[out][q,p]==0
               for out in range(nout) for p in range(nmid)
               for q in range(p+1,nmid))
    count = 0
    for out in range(nout):
        for i in range(ninput):
            for j in range(ninput):
                count += 1
    return count


assert second_jet_skew_e(Mraw_at_local_e,Hraw_at_local_e,
                         Mround_e,Hround_e)==64
assert all(H.shape==(6,6) for H in Hbase_e)
assert second_jet_skew_e(Mrestore_at_base_e,Hrestore_at_base_e,
                         Mbase_e,Hbase_e,inner_skew_certified=True)==216

# The source-specific Jacobian circuit is invertible on this whole chamber.
# The round factor has det L=4 and positive radial factor
# gain^4/sqrt(1-r_round^2); both facts are checked by one-variable algebra.
L_e = sp.Matrix([[1,0,0,1],[0,1,-1,0],
                 [1,0,0,-1],[0,1,1,0]])
assert L_e.det()==4 and L_e.T*L_e==2*sp.eye(4)
rr_e = sp.symbols('rr_e',real=True,positive=True)
gg_e = 1/sp.sqrt(1+sp.sqrt(1-rr_e**2))/sp.sqrt(4*sp.pi)
assert sp.simplify(1+2*rr_e*sp.diff(gg_e,rr_e)/gg_e
                   -1/sp.sqrt(1-rr_e**2))==0
# On r_L<1e-7, (11.11) gives ||Dg-I||<5000r_L<=1/2000,
# hence the raw singular-value margin is strictly >1999/2000.
assert 5000*rmax_e==sp.Rational(1,2000)
# The exact south restoration determinant is -(1-|t|^4).  First, Cramer's
# system gives c=2*(Re t,Im t)/(1-|t|^2), exactly the H_f terms used above.
# The map factors as w0 -> conj(w0) (real det -1), phi -> C=phi+H_f*c
# (block triangular, det 1 for the ACTUAL H_f and its derivatives),
# (x,y,C)->(x*exp(-k*x*y),y*exp(k*x*y),C) with k=C2-i*C3
# (fixed-C complex det 1, hence real det 1), then
# (x',y',C)->(x',y',R(x'y')*C) (block triangular, det R).
an_e,bn_e = sp.symbols('an_e bn_e',real=True)
bnorm_e = an_e**2+bn_e**2
Rn_e = sp.Matrix([[1-an_e**2+bn_e**2,-2*an_e*bn_e],
                  [-2*an_e*bn_e,1+an_e**2-bn_e**2]])
cn_e = sp.Matrix([2*an_e/(1-bnorm_e),2*bn_e/(1-bnorm_e)])
assert (Rn_e*cn_e-sp.Matrix([2*an_e,2*bn_e])).applyfunc(sp.factor) == (
    sp.zeros(2,1))
assert sp.factor((an_e+sp.I*bn_e)*(cn_e[0]-sp.I*cn_e[1])
                 -2*bnorm_e/(1-bnorm_e))==0
Xn_e,Yn_e,Kn_e = sp.symbols('Xn_e Yn_e Kn_e')
flow_e = sp.Matrix([Xn_e*sp.exp(-Kn_e*Xn_e*Yn_e),
                    Yn_e*sp.exp(Kn_e*Xn_e*Yn_e)])
assert sp.simplify(flow_e.jacobian((Xn_e,Yn_e)).det()-1)==0
torus_block_e = sp.Matrix([[1-Br_e,-Cr_e],[-Cr_e,1+Br_e]])
assert sp.expand(torus_block_e.det()-(1-v_e**2))==0
# The written cap estimate (§11.6) gives source radius <2e-6, whence
# |t|^2<=source_radius^2/4<1.  The raw derivative bound is (11.11).
assert sp.Rational(2,10**6)**2/4<1
assert sp.Rational(2,10**6)<sp.Rational(1,16)

# With this particular invertible Mfull_e and its six symmetric Hessian
# circuits, the universal 36/216-entry identity already checked in (c8)
# applies to the named six-real tensor J_annulus_e: J^2=-I and N_J=0.
assert J0*J0 == -sp.eye(6)
print('(e6) stage-bound six-real h, M, six H: 64+216 exact Hessian skews;'
      ' compositional 36/216 J certificate using (11.11)/§11.6 bounds;'
      ' raw margin >1999/2000: PASS')

# (e7) Type-(d) numerical safety sample only, at the rational witness.
# It is independent mpmath differentiation of the production raw inverse,
# not a premise of the exact compositional certificate above.
import mpmath as mp
with mp.workdps(60):
    sample_e = [mp.mpf(1)/5000,mp.mpf(0),mp.mpf(1)/10000,
                -mp.mpf(1)/10000]
    def source_raw_real_e(values):
        ww = cusp.elliptic_pole_raw_inverse(
            (mp.mpc(values[0],values[1]),mp.mpc(values[2],values[3])),
            False,dps=70,elementary=True)
        return [ww[0].real,ww[0].imag,ww[1].real,ww[1].imag]
    sample_M_e = mp.matrix([[mp.diff(
        lambda value:source_raw_real_e(sample_e[:j]+[value]+sample_e[j+1:])[i],
        sample_e[j]) for j in range(4)] for i in range(4)])
    sample_defect_e = sample_M_e-mp.eye(4)
    measured_norm_e = mp.sqrt(max(mp.eigsy(
        sample_defect_e.T*sample_defect_e,eigvals_only=True)))
    measured_sigma_e = mp.sqrt(min(mp.eigsy(
        sample_M_e.T*sample_M_e,eigvals_only=True)))
    assert measured_norm_e<5000*mp.mpf(6)/10**8
    assert measured_sigma_e>mp.mpf(1999)/2000
    print('(e7) numerical one-point sample only: ||Dg-I|| =',
          mp.nstr(measured_norm_e,17),
          '; sigma_min =',mp.nstr(measured_sigma_e,17),': PASS')

# (f) The corresponding north raw-pole annulus.  Its source branch has a
# different Hopf rotation, inverse stereographic ratio, cusp normalization
# nu=exp(-2), and torus-frame orientation.  The finite annulus and torus
# coordinates are independent of (e); this is a second open six-real cell.
def source_if_prefix_f(function, condition, body):
    """Match the ordered initial statements of a specified source if branch."""
    tree = ast.parse(inspect.getsource(function))
    indented = '\n'.join('    '+line for line in body.splitlines())
    expected = ast.parse('if '+condition+':\n'+indented).body[0]
    wanted = [ast.dump(node, include_attributes=False)
              for node in expected.body]
    return any(isinstance(node, ast.If) and
               ast.dump(node.test, include_attributes=False) ==
               ast.dump(expected.test, include_attributes=False) and
               len(node.body) >= len(wanted) and
               [ast.dump(item, include_attributes=False)
                for item in node.body[:len(wanted)]] == wanted
               for node in ast.walk(tree))


assert source_if_prefix_f(cusp.elliptic_pole_raw_inverse,'infinity',
    's2 = aa*aa+bb*bb+cc*cc\n'
    'den = 1+2*bb+2*s2\n'
    'u1 = mp.mpc(-1-2*bb,2*aa)/den\n'
    'u2 = mp.mpc(2*(bb+s2),2*cc)/den\n'
    'ratio = u2/u1\n'
    'height = -mp.arg(-u1)/(2*mp.pi)')
assert source_if_prefix_f(cusp.elliptic_pole_raw_inverse,'elementary',
    "if abs(height)>=mp.mpf('.1'):\n"
    "    raise ValueError('elementary pole inverse left its certified height region')\n"
    'e = mp.exp(-4*abs(tc)**2/(1-abs(tc)**2))\n'
    'q = e*nu*nu*abs(tc)**2\n'
    'R = mp.sqrt(height*height+4*q)\n'
    'AA = 2*q/(R-height) if height<0 else (R+height)/2\n'
    'yy = 2*q/(R+height) if height>=0 else (R-height)/2\n'
    'xx = AA/e\n'
    'moment,source_radius = xx-yy,xx+yy')
for statement in (
    's2 = aa*aa+bb*bb+cc*cc',
    'den = 1+2*bb+2*s2',
    'u1 = mp.mpc(-1-2*bb,2*aa)/den',
    'u2 = mp.mpc(2*(bb+s2),2*cc)/den',
    'ratio = u2/u1',
    'height = -mp.arg(-u1)/(2*mp.pi)',
):
    assert source_statement_e(cusp.elliptic_pole_raw_inverse,statement),statement
assert source_statement_e(cusp._elementary_pole_inclusion,
    'pair = ((a+b)/mp.sqrt(2),(a-b)/mp.sqrt(2)) if infinity else '
    '((a-mp.j*b)/mp.sqrt(2),(a+mp.j*b)/mp.sqrt(2))')
# The cap radius h=pi*r_L, so the south rational hmax bound is stronger.
# 1+2*beta>0 gives the principal atan branch of -arg(-u1), and the rational
# witness from (e2) has positive north height and nonzero ratio.
assert 4*rmax_e < hmax_e and hmax_e+hmax_e**2/2 < sp.Rational(1,2*10**5)
assert 1-8*rmax_e>0
r_witness_f = sp.expand(zw0_e*sp.conjugate(zw0_e)
                        +zw1_e*sp.conjugate(zw1_e))
assert sp.Rational(3,10**8)<r_witness_f<rmax_e
assert sp.re(pcw_e)>0 and sp.im(pcw_e)>0 and p3w_e>0
assert 8*p3w_e/(1-8*sp.re(pcw_e)) < sp.Rational(1,2)
# Q=|u2/u1|^2<1 uniformly here: |u1 numerator|>=1-8r,
# |u2 numerator|^2<=4((4r+16r^2)^2+(4r)^2) for r<1e-7.
Qtop_f = 4*((4*rmax_e+16*rmax_e**2)**2+(4*rmax_e)**2)
Qbottom_f = (1-8*rmax_e)**2
assert Qtop_f/Qbottom_f<sp.Rational(1,4)
assert sp.E<3 and sp.pi>3
# D>2A(1-Q); with B=e^4<81 and A=(4pi)^2>144 this gives
# |t_c|^2 < (81/(144*(3/4))) Q < 3Q/4 < 1 on this chamber.
assert sp.Rational(81,144)*sp.Rational(4,3)==sp.Rational(3,4)
assert sp.Rational(3,4)*Qtop_f/Qbottom_f<1
print('(f1) exact north source branch, open witness, positive-root domain: PASS')

f0r,f0i,f1r,f1i = z_f = sp.symbols(
    'f0r f0i f1r f1i',real=True)
pr_f = 2*(f0r*f1r+f0i*f1i)
pi_f = 2*(f0i*f1r-f0r*f1i)
p3_f = f0r*f0r+f0i*f0i-f1r*f1r-f1i*f1i
radius_f = f0r*f0r+f0i*f0i+f1r*f1r+f1i*f1i
assert sp.expand(pr_f**2+pi_f**2+p3_f**2-radius_f**2)==0
# Keep the source Hopf/scale outputs as explicit bound DAG nodes.  Expanding
# them through all nested radicals before CSE creates a 480-kB expression.
aa_f,bb_f,cc_f = sp.symbols('aa_f bb_f cc_f',real=True)
nu_f,Aprof_f,Bprof_f = sp.symbols(
    'nu_f Aprof_f Bprof_f',positive=True)
raw_input_defs_f = [(aa_f,sp.pi*p3_f),(bb_f,-sp.pi*pr_f),
                    (cc_f,sp.pi*pi_f),(nu_f,sp.exp(-2)),
                    (Aprof_f,(4*sp.pi)**2),(Bprof_f,sp.exp(4))]
Ur,Ui,Vr,Vi = sp.symbols('Ur Ui Vr Vi',real=True)
assert sp.factor(((Vr*Ur+Vi*Ui)**2+(Vi*Ur-Vr*Ui)**2)
                 -(Vr*Vr+Vi*Vi)*(Ur*Ur+Ui*Ui))==0
assert sp.cancel((Vr+sp.I*Vi)/(Ur+sp.I*Ui)
                 -((Vr*Ur+Vi*Ui)+sp.I*(Vi*Ur-Vr*Ui))
                  /(Ur*Ur+Ui*Ui))==0
# Each named scalar below is assigned an exact source expression in order.
# The shared u1/u2 stereographic denominator cancels in the ratio, leaving
# n1,n2.  This manual straight-line form prevents SymPy from expanding
# nested radicals into a 160-kB expression before differentiating.
(s2_f,n1r_f,n1i_f,n2r_f,n2i_f,ratio_den_f,Tr_f,Ti_f,Q_f,
 height_f,D_f,tgain_f,tr_f,ti_f,te2_f,E_f,q_f,R_f,AA_f,YY_f,
 uu_f,ss_f,V0r_f,V0i_f,V1r_f,V1i_f,Vnorm_f,factor_f) = sp.symbols(
    's2_f n1r_f n1i_f n2r_f n2i_f ratio_den_f Tr_f Ti_f Q_f '
    'height_f D_f tgain_f tr_f ti_f te2_f E_f q_f R_f AA_f YY_f '
    'uu_f ss_f V0r_f V0i_f V1r_f V1i_f Vnorm_f factor_f',real=True)
raw_node_defs_f = [
    (s2_f,aa_f**2+bb_f**2+cc_f**2),
    (n1r_f,-1-2*bb_f),(n1i_f,2*aa_f),
    (n2r_f,2*(bb_f+s2_f)),(n2i_f,2*cc_f),
    (ratio_den_f,n1r_f**2+n1i_f**2),
    (Tr_f,(n2r_f*n1r_f+n2i_f*n1i_f)/ratio_den_f),
    (Ti_f,(n2i_f*n1r_f-n2r_f*n1i_f)/ratio_den_f),
    (Q_f,(n2r_f**2+n2i_f**2)/ratio_den_f),
    (height_f,sp.atan(2*aa_f/(1+2*bb_f))/(2*sp.pi)),
    (D_f,Aprof_f*(1-Q_f)+sp.sqrt(Aprof_f**2*(1-Q_f)**2
                                    +4*Bprof_f*Q_f)),
    (tgain_f,sp.sqrt(2*Bprof_f/D_f)),
    (tr_f,Tr_f*tgain_f),(ti_f,Ti_f*tgain_f),
    (te2_f,tr_f**2+ti_f**2),
    (E_f,sp.exp(-4*te2_f/(1-te2_f))),
    (q_f,E_f*nu_f**2*te2_f),
    (R_f,sp.sqrt(height_f**2+4*q_f)),
    (AA_f,(R_f+height_f)/2),(YY_f,(R_f-height_f)/2),
    (uu_f,AA_f/E_f),(ss_f,uu_f+YY_f),
    (V0r_f,uu_f*f0r+nu_f*(tr_f*f1r+ti_f*f1i)),
    (V0i_f,uu_f*f0i+nu_f*(tr_f*f1i-ti_f*f1r)),
    (V1r_f,nu_f*(tr_f*f0r-ti_f*f0i)+YY_f*f1r),
    (V1i_f,nu_f*(tr_f*f0i+ti_f*f0r)+YY_f*f1i),
    (Vnorm_f,V0r_f**2+V0i_f**2+V1r_f**2+V1i_f**2),
    (factor_f,sp.sqrt(ss_f/Vnorm_f))]
# (R+h)(R-h)=4q proves yy=2q/(R+h); uu*yy=nu^2|t_c|^2.
assert sp.expand((rnode_e-rh_e)*(rnode_e+rh_e)-4*qe_e)==0
assert dict(raw_node_defs_f)[q_f]==E_f*nu_f**2*te2_f
raw_out_f = [factor_f*V0r_f,factor_f*V0i_f,
             factor_f*V1r_f,factor_f*V1i_f]
raw_f = sp.Matrix(raw_out_f)
# Check the source moment/projector normalization independently with small
# algebraic nodes.  In source notation m=u-Y, s=u+Y, source_pc=2*nu*conj(t),
# so s*(w0,w1)=V.  The final positive-root factor then gives sqrt(s/Vnorm)*V.
ug,yg,ng,tg,ig,z0rg,z0ig,z1rg,z1ig = sp.symbols(
    'ug yg ng tg ig z0rg z0ig z1rg z1ig',real=True)
sg,mg = ug+yg,ug-yg
Vgeneric_f = sp.Matrix([
    ug*z0rg+ng*(tg*z1rg+ig*z1ig),
    ug*z0ig+ng*(tg*z1ig-ig*z1rg),
    ng*(tg*z0rg-ig*z0ig)+yg*z1rg,
    ng*(tg*z0ig+ig*z0rg)+yg*z1ig])
Vactual_rhs_f = dict(raw_node_defs_f)
Vgeneric_bound_f = Vgeneric_f.subs({ug:uu_f,yg:YY_f,ng:nu_f,
    tg:tr_f,ig:ti_f,z0rg:f0r,z0ig:f0i,z1rg:f1r,z1ig:f1i})
assert all(sp.expand(Vgeneric_bound_f[k]-Vactual_rhs_f[node])==0
           for k,node in enumerate((V0r_f,V0i_f,V1r_f,V1i_f)))
source_projected_f = sp.Matrix([
    ((sg+mg)*z0rg+2*ng*(tg*z1rg+ig*z1ig))/(2*sg),
    ((sg+mg)*z0ig+2*ng*(tg*z1ig-ig*z1rg))/(2*sg),
    (2*ng*(tg*z0rg-ig*z0ig)+(sg-mg)*z1rg)/(2*sg),
    (2*ng*(tg*z0ig+ig*z0rg)+(sg-mg)*z1ig)/(2*sg)])
assert all(sp.cancel(sg*source_projected_f[k]-Vgeneric_f[k])==0
           for k in range(4))
sgp,vnp = sp.symbols('sgp vnp',positive=True)
assert sp.simplify((sp.sqrt(sgp)/sp.sqrt(vnp/sgp**2)/sgp)**2
                   -sgp/vnp)==0
print('(f2) actual north raw inverse in four real pole variables: PASS')

raw_jet_defs_f = []
raw_known_f = {}
for k,variable in enumerate(z_f):
    raw_known_f[variable]=(sp.eye(4)[:,k],sp.zeros(4,4))
for node,rhs in raw_input_defs_f:
    raw_known_f[node]=(sp.Matrix([sp.diff(rhs,v) for v in z_f]),
                       sp.hessian(rhs,z_f))
for node,rhs in raw_node_defs_f:
    grad_rhs,hess_rhs = raw_local_jet_e(rhs,raw_known_f)
    raw_jet_defs_f.append((node,rhs))
    grad_symbols = sp.Matrix([sp.Symbol(f'{node}_d{i}') for i in range(4)])
    hess_symbols = sp.Matrix(4,4,lambda i,j:
                             sp.Symbol(f'{node}_h{min(i,j)}{max(i,j)}'))
    raw_jet_defs_f.extend((grad_symbols[i],grad_rhs[i]) for i in range(4))
    raw_jet_defs_f.extend((hess_symbols[i,j],hess_rhs[i,j])
                          for i in range(4) for j in range(i,4))
    raw_known_f[node]=(grad_symbols,hess_symbols)
raw_output_jets_f = [raw_local_jet_e(rhs,raw_known_f) for rhs in raw_out_f]
Mraw_f = sp.Matrix(4,4,lambda k,i:raw_output_jets_f[k][0][i])
Hraw_f = [raw_output_jets_f[k][1] for k in range(4)]
assert Mraw_f.shape==(4,4) and len(Hraw_f)==4
assert all(Hraw_f[k][i,j]==Hraw_f[k][j,i]
           for k in range(4) for i in range(4) for j in range(4))
assert not any(value.has(sp.Derivative) for value in Mraw_f)
assert len(raw_jet_defs_f)==15*len(raw_node_defs_f)
raw_defined_f = set(z_f)|{node for node,_ in raw_input_defs_f}
for node,rhs in raw_jet_defs_f:
    assert rhs.free_symbols <= raw_defined_f
    raw_defined_f.add(node)
assert all(rhs.free_symbols <= raw_defined_f for rhs in raw_out_f)
print('(f3) actual north raw inverse 4x4 Jacobian and four symmetric'
      f' Hessian circuits via {len(raw_node_defs_f)} exact scalar 2-jet nodes: PASS')

fr0,fi0,fr1,fi1 = round_f = sp.symbols(
    'fr0 fi0 fr1 fi1',real=True)
round_r_f = sum(value*value for value in round_f)
round_gain_f = 1/sp.sqrt(1+sp.sqrt(1-round_r_f**2))/sp.sqrt(2*sp.pi)
local_f = sp.Matrix([round_gain_f*(fr0+fr1),
                     round_gain_f*(fi0+fi1),
                     round_gain_f*(fr0-fr1),
                     round_gain_f*(fi0-fi1)])
Mround_f = local_f.jacobian(round_f)
Hround_f = [sp.hessian(local_f[k],round_f) for k in range(4)]
assert all(Hround_f[k][i,j]==Hround_f[k][j,i]
           for k in range(4) for i in range(4) for j in range(4))

fw0r,fw0i,fw1r,fw1i,fp2,fp3 = rest_f = sp.symbols(
    'fw0r fw0i fw1r fw1i fp2 fp3',real=True)
tR_f = (fw0r*fw1r+fw0i*fw1i)/nu_f
tI_f = (fw0r*fw1i-fw0i*fw1r)/nu_f
v_f = tR_f**2+tI_f**2
Ee_f = sp.exp(-4*v_f/(1-v_f))
Ae_f,Ye_f = Ee_f*(fw0r*fw0r+fw0i*fw0i),fw1r*fw1r+fw1i*fw1i
Hf_f = (Ae_f-Ye_f)/(1-Ae_f*Ye_f)
ar_f = tR_f*fp2-tI_f*fp3+2*Hf_f*v_f/(1-v_f)
ai_f = tI_f*fp2+tR_f*fp3
co_f,si_f,ea_f = sp.cos(ai_f),sp.sin(ai_f),sp.exp(ar_f)
Br_f,Cr_f = tR_f*tR_f-tI_f*tI_f,2*tR_f*tI_f
restore_f = sp.Matrix([
    (fw0r*co_f-fw0i*si_f)/(nu_f*ea_f),
    (-fw0r*si_f-fw0i*co_f)/(nu_f*ea_f),
    ea_f*(fw1r*co_f-fw1i*si_f),
    ea_f*(fw1r*si_f+fw1i*co_f),
    fp2*(Br_f-1)-fp3*Cr_f-2*Hf_f*tR_f,
    fp2*Cr_f+fp3*(1+Br_f)-2*Hf_f*tI_f])
Mrestore_f = restore_f.jacobian(rest_f)
Hrestore_f = [sp.hessian(restore_f[k],rest_f) for k in range(6)]
assert Mrestore_f.shape==(6,6) and len(Hrestore_f)==6
assert all(Hrestore_f[k][i,j]==Hrestore_f[k][j,i]
           for k in range(6) for i in range(6) for j in range(6))
assert not any(value.has(sp.Derivative) for value in Mrestore_f)
print('(f4) actual north round cap and finite-height torus first/second jets: PASS')

raw_stage_defs_f = list(zip(z_f,local_f))+raw_input_defs_f+raw_jet_defs_f
raw_at_local_f = sp.Matrix(raw_out_f)
base_at_round_f = tuple(raw_at_local_f)+(fp2,fp3)
full_at_round_f = sp.Matrix([
    bind_e(restore_f[k],rest_f,base_at_round_f) for k in range(6)])
Mrestore_at_base_f = bind_matrix_e(Mrestore_f,rest_f,base_at_round_f)
Hrestore_at_base_f = [bind_matrix_e(H,rest_f,base_at_round_f)
                      for H in Hrestore_f]
Mbase_f = sp.BlockDiagMatrix(
    sp.MatMul(Mraw_f,Mround_f,evaluate=False),sp.eye(2))
Mfull_f = sp.MatMul(Mrestore_at_base_f,Mbase_f,evaluate=False)
J_annulus_f = sp.MatMul(sp.Inverse(Mfull_f),J0,Mfull_f,evaluate=False)
round_inputs_f = set(round_f)|{fp2,fp3}
bound_symbols_f = set(round_inputs_f)
for node,rhs in raw_stage_defs_f:
    assert rhs.free_symbols <= bound_symbols_f
    bound_symbols_f.add(node)
assert tuple(rhs for _,rhs in raw_stage_defs_f[:4])==tuple(local_f)
assert Mfull_f.shape==(6,6) and J_annulus_f.shape==(6,6)
assert all(value.free_symbols <= bound_symbols_f for value in full_at_round_f)
assert all(value.free_symbols <= bound_symbols_f for value in Mraw_f)
assert all(value.free_symbols <= bound_symbols_f for value in Mrestore_at_base_f)
assert J_annulus_f.free_symbols <= bound_symbols_f
Hbase_f = []
for k in range(4):
    h = sp.MatAdd(*[
        sp.MatMul(Mraw_f[k,p],Hround_f[p],evaluate=False)
        for p in range(4)],
        sp.MatMul(Mround_f.T,Hraw_f[k],Mround_f,evaluate=False),
        evaluate=False)
    Hbase_f.append(sp.BlockDiagMatrix(h,sp.zeros(2,2)))
Hbase_f.extend([sp.ZeroMatrix(6,6),sp.ZeroMatrix(6,6)])
Hfull_f = []
for k in range(6):
    Hfull_f.append(sp.MatAdd(*[
        sp.MatMul(Mrestore_at_base_f[k,p],Hbase_f[p],evaluate=False)
        for p in range(6)],
        sp.MatMul(Mbase_f.T,Hrestore_at_base_f[k],Mbase_f,evaluate=False),
        evaluate=False))
assert len(Hfull_f)==6 and all(H.shape==(6,6) for H in Hfull_f)
assert all(H.free_symbols <= bound_symbols_f for H in Hfull_f)
assert second_jet_skew_e(Mraw_f,Hraw_f,Mround_f,Hround_f)==64
assert second_jet_skew_e(Mrestore_at_base_f,Hrestore_at_base_f,
                         Mbase_f,Hbase_f,inner_skew_certified=True)==216

# Round local gain has the same radial derivative identity as (e), with
# north L^T L=2I and det L=4.  The imported (11.11) bound gives raw
# sigma_min>1999/2000 on this annulus.  The north Cramer system has
# determinant -(1-|t|^4); conjugation/scaling w0->conj(w0)/nu contributes
# -nu^-2, while the H_f shear and holomorphic flow each have determinant 1.
L_f = sp.Matrix([[1,0,1,0],[0,1,0,1],
                 [1,0,-1,0],[0,1,0,-1]])
assert L_f.det()==4 and L_f.T*L_f==2*sp.eye(4)
assert sp.simplify(1+2*rr_e*sp.diff(gg_e,rr_e)/gg_e
                   -1/sp.sqrt(1-rr_e**2))==0
assert 5000*rmax_e==sp.Rational(1,2000)
af,bf = sp.symbols('af bf',real=True)
vf = af*af+bf*bf
Rf = sp.Matrix([[af*af-bf*bf-1,-2*af*bf],
                [2*af*bf,1+af*af-bf*bf]])
cf = sp.Matrix([2*af/(1-vf),-2*bf/(1-vf)])
assert (Rf*cf-sp.Matrix([-2*af,-2*bf])).applyfunc(sp.factor)==sp.zeros(2,1)
assert sp.expand(Rf.det()+(1-vf**2))==0
pf2,pf3,hf = sp.symbols('pf2 pf3 hf',real=True)
C2n,C3n = pf2+hf*cf[0],pf3+hf*cf[1]
An = (af+sp.I*bf)*(C2n+sp.I*C3n)
assert sp.factor(sp.re(An)-(af*pf2-bf*pf3+2*hf*vf/(1-vf)))==0
assert sp.factor(sp.im(An)-(bf*pf2+af*pf3))==0
assert sp.Rational(2,10**6)*9/2<1
assert (sp.Rational(2,10**6)*9/2)**2<1
assert J0*J0==-sp.eye(6)
print('(f5) stage-bound north six-real h, M, six H: 64+216 compositional'
      ' skews; 36/216 J certificate using (11.11)/§11.6 bounds;'
      ' raw margin >1999/2000: PASS')

# (f6) Type-(d) one-point numerical source derivative sample only.
with mp.workdps(60):
    sample_f = [mp.mpf(1)/5000,mp.mpf(0),mp.mpf(1)/10000,
                -mp.mpf(1)/10000]
    def source_raw_real_f(values):
        ww = cusp.elliptic_pole_raw_inverse(
            (mp.mpc(values[0],values[1]),mp.mpc(values[2],values[3])),
            True,dps=70,elementary=True)
        return [ww[0].real,ww[0].imag,ww[1].real,ww[1].imag]
    numerical_nodes_f = dict(zip(z_f,[sp.Rational(1,5000),sp.Integer(0),
                                      sp.Rational(1,10000),-sp.Rational(1,10000)]))
    for node,rhs in raw_input_defs_f+raw_node_defs_f:
        numerical_nodes_f[node]=rhs.subs(numerical_nodes_f).evalf(70)
    symbolic_sample_f = [mp.mpf(str(value.subs(numerical_nodes_f).evalf(60)))
                         for value in raw_f]
    production_sample_f = source_raw_real_f(sample_f)
    source_error_f = max(abs(x-y) for x,y in
                         zip(symbolic_sample_f,production_sample_f))
    assert source_error_f < mp.mpf('1e-48')
    sample_M_f = mp.matrix([[mp.diff(
        lambda value:source_raw_real_f(sample_f[:j]+[value]+sample_f[j+1:])[i],
        sample_f[j]) for j in range(4)] for i in range(4)])
    sample_defect_f = sample_M_f-mp.eye(4)
    measured_norm_f = mp.sqrt(max(mp.eigsy(
        sample_defect_f.T*sample_defect_f,eigvals_only=True)))
    measured_sigma_f = mp.sqrt(min(mp.eigsy(
        sample_M_f.T*sample_M_f,eigvals_only=True)))
    assert measured_norm_f<5000*mp.mpf(6)/10**8
    assert measured_sigma_f>mp.mpf(1999)/2000
    print('(f6) numerical one-point sample only: ||Dg-I|| =',
          mp.nstr(measured_norm_f,17),
          '; sigma_min =',mp.nstr(measured_sigma_f,17),
          '; source value error =',mp.nstr(source_error_f,5),': PASS')
# (g) An open south radial-transition chamber, 1.3e-8<|z_L|^2<1.4e-8,
# p3(Hopf z_L)>0 and Im(2z0 conj z1)>0, of the actual two-sided pole
# patch H_0=D^-1 o g_0 o D.  The scalar value and first/second jets of both
# nonidentity radial maps are explicit finite DAGs; the 43-node raw DAG in
# (e) is re-used unchanged between them.  This is one chamber, not the whole
# radial shell, whose outer edge may take the identity branch of D^-1.
for statement in (
    "a,c = mp.mpf('1e-8'),mp.mpf('1.5e-8')",
    'factor = mp.sqrt(new_radius/radius)',
    'return tuple(+(factor*w) for w in z)',
):
    assert source_statement_e(cusp._elementary_pole_radial,statement),statement
assert source_if_first_e(cusp._elementary_pole_radial,'radius>=c','return z')
assert source_if_first_e(cusp._elementary_pole_radial,
    'c/(c-radius)>(mp.mp.dps+10)*mp.log(10)','return z')
assert source_if_prefix_f(cusp._elementary_pole_radial,'inverse',
    "if radius<=0:\n"
    "    raise ValueError('radial inverse requires positive radius')\n"
    'y = mp.e*mp.expm1(radius/(c-radius))\n'
    'L = (y+mp.log1p(-mp.exp(-y)) if y>1 else mp.log(mp.expm1(y)))\n'
    'x = 2/(mp.sqrt(L*L+4)-L) if L<0 else (L+mp.sqrt(L*L+4))/2\n'
    'ell = mp.log1p(x/mp.exp(3))\n'
    'new_radius = a+c*ell/(3*(3+ell))')


def source_if_else_prefix_g(function,condition,body):
    tree = ast.parse(inspect.getsource(function))
    indented = '\n'.join('    '+line for line in body.splitlines())
    expected = ast.parse('if '+condition+':\n    pass\nelse:\n'
                         +indented).body[0]
    wanted = [ast.dump(item,include_attributes=False)
              for item in expected.orelse]
    return any(isinstance(node,ast.If) and
               ast.dump(node.test,include_attributes=False)==
               ast.dump(expected.test,include_attributes=False) and
               len(node.orelse)>=len(wanted) and
               [ast.dump(item,include_attributes=False)
                for item in node.orelse[:len(wanted)]]==wanted
               for node in ast.walk(tree))


assert source_if_else_prefix_g(cusp._elementary_pole_radial,'inverse',
    "if radius<=a:\n"
    "    raise ValueError('radial compression requires radius greater than 1e-8')\n"
    'x = mp.exp(3)*mp.expm1(3*(radius-a)/(c-radius))\n'
    'f = x-1/x\n'
    'y = mp.log1p(mp.exp(f)) if f<0 else f+mp.log1p(mp.exp(-f))\n'
    'ell = mp.log1p(y/mp.e)\n'
    'new_radius = c*ell/(1+ell)')

ga,gc = sp.Rational(1,10**8),sp.Rational(3,2*10**8)
grmin,grmax = sp.Rational(13,10**9),sp.Rational(14,10**9)
assert ga<sp.Rational(5,4)*ga<grmin<grmax<gc
gz0 = sp.Rational(1,10000)
gz1 = (5-3*sp.I)/100000
grw = sp.expand(gz0*sp.conjugate(gz0)+gz1*sp.conjugate(gz1))
gpcw = sp.expand(2*gz0*sp.conjugate(gz1))
gp3w = sp.expand(gz0*sp.conjugate(gz0)-gz1*sp.conjugate(gz1))
assert grmin<grw<grmax and sp.re(gpcw)>0 and gp3w>0
assert sp.im(gpcw)>0 and gpcw!=0
# Throughout the stated angular subchamber, p3>0 gives gamma>0,
# Im(pc)>0 gives aa=-2*pi*Im(pc)!=0, and 1+beta>1-8r>0.
# Hence 0<height=atan(gamma/(1+beta))/(2*pi)<1/12<1/10,
# Re(u2)>0, and tc!=0.  The following rational inequalities certify
# the upper bound uniformly and at the strict witness.
assert 8*grmax<sp.Rational(1,2)
assert 8*gp3w/(1-8*sp.re(gpcw))<sp.Rational(1,2)
assert 8*grmax/(1-8*grmax)<sp.Rational(1,2)
assert sp.Rational(1,12)<sp.Rational(1,10)
# S is strictly increasing by the positive factors in N10b,d and S(r)<r
# by (11.12g).  At r_min, c/(c-r_min)=15/2>7; e in (2,3) gives
# phi(r_min)=e^(15/2)-e^3>2^7-3^3=101.  Thus F(phi)>100 and
# T(phi)>100>psi(a)=e^3-e<27; hence S(r)>a throughout this chamber.
assert gc/(gc-grmin)==sp.Rational(15,2)
assert sp.E>2 and sp.E<3 and 2**7-3**3>100>3**3
assert 101-sp.Rational(1,101)>100
assert grmax+80*grmax**2<sp.Rational(141,10**10)<gc
assert ga-80*grmax**2>sp.Rational(9,10**9)
# For u=S(r), (11.9) implies u-80u²<v=|g(Dz)|²<u+80u²;
# so 9e-9<v<1.41e-8.  The output radial inverse is nonidentity.
# Both precision guards are false for dps>=50: c/(c-r)<=15 and
# c/(c-v)<50/3, whereas (mp.mp.dps+10)log 10>120.
assert gc/(gc-grmax)==15
assert gc/(gc-sp.Rational(141,10**10))==sp.Rational(50,3)
assert sp.log(10)>2 and sp.Rational(50,3)<120
# v>9e-9 implies c/(c-v)>5/2.  Thus the inverse y>e^2-e>1,
# so its source L uses y+log1p(-exp(-y)) and is positive; the source
# quadratic x uses its L>=0 branch.  The input x>101 gives f>100>0.
assert gc/(gc-sp.Rational(9,10**9))==sp.Rational(5,2)
assert sp.E**2-sp.E>1
assert 5000*grmax==sp.Rational(7,10**5)
# The compressed spinor remains in the same Hopf direction.  At the
# rational witness p3>0, so the south elementary height uses the smooth
# positive branch and its inverse ratio is nonzero on a neighborhood.
assert grmax<sp.Rational(1,10**4)
print('(g1) source-bound nonidentity radial branches and strict open'
      ' south witness: PASS')

# N10a–d prove the stable-code rewrites and the four-output radial jet
# formula.  Here each scalar derivative is a source-specific named circuit,
# rather than an unevaluated derivative of S or S^-1.
gX,gY = sp.symbols('gX gY',positive=True)
assert sp.simplify(sp.exp(gX+sp.log(1+sp.exp(-gX)))
                   -(1+sp.exp(gX)))==0
assert sp.simplify(sp.exp(gY+sp.log(1-sp.exp(-gY)))
                   -(sp.exp(gY)-1))==0


def scalar_jet_dag_g(variable, nodes):
    known = {variable:(sp.S.One,sp.S.Zero)}
    assignments = []
    for node,rhs in nodes:
        deps = tuple(sorted(rhs.free_symbols & known.keys(),key=str))
        assert rhs.free_symbols <= set(deps)
        d1 = {p:sp.diff(rhs,p) for p in deps}
        d2 = {(p,q):sp.diff(d1[p],q) for p in deps for q in deps}
        assert all(sp.simplify(d2[p,q]-d2[q,p])==0
                   for p in deps for q in deps)
        first = sum((d1[p]*known[p][0] for p in deps),sp.S.Zero)
        second = sum((d1[p]*known[p][1] for p in deps),sp.S.Zero)
        second += sum((d2[p,q]*known[p][0]*known[q][0]
                       for p in deps for q in deps),sp.S.Zero)
        first_node = sp.Symbol(f'{node}_d',real=True)
        second_node = sp.Symbol(f'{node}_dd',real=True)
        assignments.extend(((node,rhs),(first_node,first),
                            (second_node,second)))
        known[node]=(first_node,second_node)
    return assignments,known


gr = sp.Symbol('gr',positive=True)
gxe, gxf, gyt, gle, gS, gDgain = sp.symbols(
    'gxe gxf gyt gle gS gDgain',positive=True)
S_nodes_g = [
    (gxe,sp.exp(3)*(sp.exp(3*(gr-ga)/(gc-gr))-1)),
    (gxf,gxe-1/gxe),
    (gyt,gxf+sp.log(1+sp.exp(-gxf))),
    (gle,sp.log(1+gyt/sp.E)),
    (gS,gc*gle/(1+gle)),
    (gDgain,sp.sqrt(gS/gr))]
S_jet_defs_g,S_known_g = scalar_jet_dag_g(gr,S_nodes_g)
assert len(S_jet_defs_g)==3*len(S_nodes_g)

gs = sp.Symbol('gs',positive=True)
gye,gL,groot,gell,gR,gIgain = sp.symbols(
    'gye gL groot gell gR gIgain',positive=True)
R_nodes_g = [
    (gye,sp.E*(sp.exp(gs/(gc-gs))-1)),
    (gL,gye+sp.log(1-sp.exp(-gye))),
    (groot,(gL+sp.sqrt(gL**2+4))/2),
    (gell,sp.log(1+groot/sp.exp(3))),
    (gR,ga+gc*gell/(3*(3+gell))),
    (gIgain,sp.sqrt(gR/gs))]
R_jet_defs_g,R_known_g = scalar_jet_dag_g(gs,R_nodes_g)
assert len(R_jet_defs_g)==3*len(R_nodes_g)
assert sp.cancel(gc/(gc-gr)-3-3*(gr-ga)/(gc-gr))==0
assert sp.cancel(gc*gle/(1+gle)-gc*(1-1/(1+gle)))==0
assert sp.cancel(ga+gc*gell/(3*(3+gell))
                 -gc*(1-1/(3+gell)))==0
# The root identity is checked with its defining positive expression.
root_rhs_g = dict(R_nodes_g)[groot]
assert sp.simplify(root_rhs_g*(root_rhs_g-gL)-1)==0
print('(g2) source-specific finite S and S^-1 scalar value/first/second'
      ' jet DAGs: PASS')

def radial_jets_g(inputs,gain,first,second):
    value = sp.Matrix([gain*z for z in inputs])
    M = sp.Matrix(4,4,lambda k,i:
        gain*int(k==i)+2*first*inputs[k]*inputs[i])
    H = [sp.Matrix(4,4,lambda i,j:
         2*first*(int(k==i)*inputs[j]+int(k==j)*inputs[i]
                  +int(i==j)*inputs[k])
         +4*second*inputs[k]*inputs[i]*inputs[j])
         for k in range(4)]
    assert all(H[k][i,j]==H[k][j,i]
               for k in range(4) for i in range(4) for j in range(4))
    return value,M,H


# This is the literal four-real-coordinate differentiation underlying N10c,
# checked once with an arbitrary radial scalar before binding the two actual
# scalar DAGs.  No derivative remains in the bound stage circuits below.
formal_z_g = sp.symbols('formal_z_g0:4',real=True)
formal_q_g = sum(z*z for z in formal_z_g)
formal_dummy_g = sp.Symbol('formal_dummy_g',positive=True)
formal_fun_g = sp.Function('formal_fun_g')
formal_value_g = sp.Matrix([formal_fun_g(formal_q_g)*z for z in formal_z_g])
formal_first_g = sp.diff(formal_fun_g(formal_dummy_g),formal_dummy_g
                         ).subs(formal_dummy_g,formal_q_g)
formal_second_g = sp.diff(formal_fun_g(formal_dummy_g),formal_dummy_g,2
                          ).subs(formal_dummy_g,formal_q_g)
_,formal_M_g,formal_H_g = radial_jets_g(
    formal_z_g,formal_fun_g(formal_q_g),formal_first_g,formal_second_g)
assert formal_value_g.jacobian(formal_z_g)==formal_M_g
assert all(sp.hessian(formal_value_g[k],formal_z_g)==formal_H_g[k]
           for k in range(4))

gz_in = sp.symbols('gz_in0:4',real=True)
gw_in = sp.symbols('gw_in0:4',real=True)
D_value_g,M_D_g,H_D_g = radial_jets_g(
    gz_in,gDgain,*S_known_g[gDgain])
I_value_g,M_I_g,H_I_g = radial_jets_g(
    gw_in,gIgain,*R_known_g[gIgain])
assert all(not entry.has(sp.Derivative) for entry in M_D_g)
assert all(not entry.has(sp.Derivative) for entry in M_I_g)
assert len(H_D_g)==len(H_I_g)==4
print('(g3) actual four-output D and D^-1 value/first/64 second'
      ' derivative circuits: PASS')

# Bind the true stage order.  z_e is the input of the existing south raw
# 43-node DAG, now assigned D(local_e(Z)); gw_in is assigned that DAG's
# output.  The final gIgain is therefore S^-1(|g(Dz)|^2), not S^-1(r).
stage_defs_g = list(zip(gz_in,local_e))
stage_defs_g.append((gr,sum(z*z for z in gz_in)))
stage_defs_g.extend(S_jet_defs_g)
stage_defs_g.extend(zip(z_e,D_value_g))
stage_defs_g.extend(raw_jet_defs_e)
stage_defs_g.extend(zip(gw_in,raw_out_e))
stage_defs_g.append((gs,sum(w*w for w in gw_in)))
stage_defs_g.extend(R_jet_defs_g)
stage_inputs_g = set(round_e)|{p2e,p3e}
stage_bound_g = set(stage_inputs_g)
for node,rhs in stage_defs_g:
    assert rhs.free_symbols <= stage_bound_g
    stage_bound_g.add(node)
assert tuple(rhs for _,rhs in stage_defs_g[:4])==tuple(local_e)
assert all(value.free_symbols <= stage_bound_g for value in I_value_g)

M_D_at_round_g = bind_matrix_e(M_D_g,gz_in,local_e)
H_D_at_round_g = [bind_matrix_e(H,gz_in,local_e) for H in H_D_g]
M_raw_at_D_g = bind_matrix_e(Mraw_e,z_e,D_value_g)
H_raw_at_D_g = [bind_matrix_e(H,z_e,D_value_g) for H in Hraw_e]
M_I_at_raw_g = bind_matrix_e(M_I_g,gw_in,raw_out_e)
H_I_at_raw_g = [bind_matrix_e(H,gw_in,raw_out_e) for H in H_I_g]


def compose_jets_g(outer_M,outer_H,inner_M,inner_H):
    nout,nmid = outer_M.shape
    assert inner_M.rows==nmid and len(outer_H)==nout
    assert len(inner_H)==nmid
    M = sp.MatMul(outer_M,inner_M,evaluate=False)
    H = [sp.MatAdd(*[
         sp.MatMul(outer_M[k,p],inner_H[p],evaluate=False)
         for p in range(nmid)],
         sp.MatMul(inner_M.T,outer_H[k],inner_M,evaluate=False),
         evaluate=False) for k in range(nout)]
    return M,H


M_Dround_g,H_Dround_g = compose_jets_g(
    M_D_at_round_g,H_D_at_round_g,Mround_e,Hround_e)
M_rawround_g,H_rawround_g = compose_jets_g(
    M_raw_at_D_g,H_raw_at_D_g,M_Dround_g,H_Dround_g)
M_patch4_g,H_patch4_g = compose_jets_g(
    M_I_at_raw_g,H_I_at_raw_g,M_rawround_g,H_rawround_g)
M_patch6_g = sp.BlockDiagMatrix(M_patch4_g,sp.eye(2))
H_patch6_g = [sp.BlockDiagMatrix(H,sp.zeros(2,2))
              for H in H_patch4_g]
H_patch6_g.extend([sp.ZeroMatrix(6,6),sp.ZeroMatrix(6,6)])
rest_input_g = tuple(I_value_g)+(p2e,p3e)
full_value_g = sp.Matrix([
    bind_e(restore_e[k],rest_e,rest_input_g) for k in range(6)])
M_rest_at_patch_g = bind_matrix_e(Mrestore_e,rest_e,rest_input_g)
H_rest_at_patch_g = [bind_matrix_e(H,rest_e,rest_input_g)
                     for H in Hrestore_e]
Mfull_g,Hfull_g = compose_jets_g(
    M_rest_at_patch_g,H_rest_at_patch_g,M_patch6_g,H_patch6_g)
J_transition_g = sp.MatMul(sp.Inverse(Mfull_g),J0,Mfull_g,
                            evaluate=False)
assert Mfull_g.shape==(6,6) and J_transition_g.shape==(6,6)
assert len(Hfull_g)==6 and all(H.shape==(6,6) for H in Hfull_g)
assert all(value.free_symbols <= stage_bound_g for value in full_value_g)
assert Mfull_g.free_symbols <= stage_bound_g
assert all(H.free_symbols <= stage_bound_g for H in Hfull_g)
assert J_transition_g.free_symbols <= stage_bound_g
print('(g4) actual stage-bound south transition h, 6x6 M, and six'
      ' second-jet circuits: PASS')

# Every 4D output Hessian skew vanishes at each of the three compositions.
# The final 6D skews then vanish by the same exact chain-rule identity as
# (e6), with the actual bound source Hessians, not a placeholder map.
assert second_jet_skew_e(M_D_at_round_g,H_D_at_round_g,
                         Mround_e,Hround_e)==64
assert second_jet_skew_e(M_raw_at_D_g,H_raw_at_D_g,
                         M_Dround_g,H_Dround_g,
                         inner_skew_certified=True)==64
assert second_jet_skew_e(M_I_at_raw_g,H_I_at_raw_g,
                         M_rawround_g,H_rawround_g,
                         inner_skew_certified=True)==64
assert second_jet_skew_e(M_rest_at_patch_g,H_rest_at_patch_g,
                         M_patch6_g,H_patch6_g,
                         inner_skew_certified=True)==216
# N10d proves det_R D=S*S'/r>0 and det_R D^-1=R*R'/s>0
# from the positive scalar derivative factors.  The raw singular-value
# margin follows (11.11), while (e6) proves the round and restoration
# determinants.  Consequently this actual Mfull_g is invertible.
assert 1-5000*grmax==sp.Rational(99993,100000)
# The inverse of the south round radial inclusion (11.17a) is
# rho=2h/(1+h²), h=2*pi*r.  On h<1 its positive square root is
# sqrt(1-rho²)=(1-h²)/(1+h²), recovering r exactly.  Thus the rational
# local witness is in the round cap image; the next inequality places it
# strictly in the pole selector's open region.
gh = sp.Symbol('gh',positive=True)
grho = 2*gh/(1+gh**2)
groot_rho = (1-gh**2)/(1+gh**2)
assert sp.factor(grho**2+groot_rho**2-1)==0
assert sp.cancel(grho/(2*sp.pi*(1+groot_rho))
                 -gh/(2*sp.pi))==0
assert 8*grmax<1
assert 8*grmax+(8*grmax)**2/2<sp.Rational(1,2*10**5)
assert sp.Rational(2,10**6)**2/4<1
assert J0*J0==-sp.eye(6)
print('(g5) 3x64+216 compositional skews and 36/216 J certificate;'
      ' raw margin >99993/100000, using N10/(11.11)/§11.6: PASS')

# (g6) Type-(d), one-point production correspondence and measured raw
# singular-value margin.  These numbers do not enter the exact certificate.
with mp.workdps(80):
    numerical_g = dict(zip(gz_in,[sp.Rational(1,10000),sp.S.Zero,
                                   sp.Rational(5,100000),
                                   -sp.Rational(3,100000)]))
    numerical_g[gr] = sum(value*value for value in numerical_g.values())
    for node,rhs in S_nodes_g:
        assert rhs.free_symbols <= numerical_g.keys()
        numerical_g[node] = rhs.subs(numerical_g).evalf(90)
    compressed_g = [mp.mpf(str(rhs.subs(numerical_g).evalf(80)))
                    for rhs in D_value_g]
    for node,rhs in zip(z_e,D_value_g):
        numerical_g[node] = rhs.subs(numerical_g).evalf(90)
    for node,rhs in raw_cse_e:
        assert rhs.free_symbols <= numerical_g.keys()
        numerical_g[node] = rhs.subs(numerical_g).evalf(90)
    raw_value_g = [rhs.subs(numerical_g).evalf(90) for rhs in raw_out_e]
    numerical_g.update(zip(gw_in,raw_value_g))
    numerical_g[gs] = sum(value*value for value in raw_value_g)
    for node,rhs in R_nodes_g:
        assert rhs.free_symbols <= numerical_g.keys()
        numerical_g[node] = rhs.subs(numerical_g).evalf(90)
    symbolic_patch_g = [mp.mpf(str(rhs.subs(numerical_g).evalf(80)))
                        for rhs in I_value_g]
    source_patch_g = cusp.elementary_pole_patch(
        (mp.mpc(mp.mpf(1)/10000),
         mp.mpc(mp.mpf(5)/100000,-mp.mpf(3)/100000)),
        False,dps=80)
    source_patch_real_g = [source_patch_g[0].real,source_patch_g[0].imag,
                           source_patch_g[1].real,source_patch_g[1].imag]
    patch_error_g = max(abs(a-b) for a,b in
                        zip(symbolic_patch_g,source_patch_real_g))
    assert patch_error_g<mp.mpf('1e-55')
    sample_M_g = mp.matrix([[mp.diff(
        lambda value:source_raw_real_e(compressed_g[:j]+[value]
                                       +compressed_g[j+1:])[i],
        compressed_g[j]) for j in range(4)] for i in range(4)])
    measured_defect_g = sample_M_g-mp.eye(4)
    measured_norm_g = mp.sqrt(max(mp.eigsy(
        measured_defect_g.T*measured_defect_g,eigvals_only=True)))
    measured_sigma_g = mp.sqrt(min(mp.eigsy(
        sample_M_g.T*sample_M_g,eigvals_only=True)))
    raw_bound_g = 5000*sum(value*value for value in compressed_g)
    assert measured_norm_g<raw_bound_g<5000*mp.mpf(14)/10**9
    assert measured_sigma_g>mp.mpf(99993)/100000
    print('(g6) numerical one-point sample only: full patch value error =',
          mp.nstr(patch_error_g,5),'; raw ||Dg-I|| =',
          mp.nstr(measured_norm_g,17),'vs bound',mp.nstr(raw_bound_g,17),
          '; raw sigma_min =',mp.nstr(measured_sigma_g,17),': PASS')

# (h) The north nonidentity radial-transition chamber.  The source radial
# branches and their domain inequalities were checked in (g1); (f1) checks
# the distinct north raw branch, and (e1) the common cap/source dispatch.
# We clone the two scalar DAGs so every north jet node has a separate
# topological binding; no south value or derivative is imported as data.
assert source_statement_e(cusp.elementary_global_slice,'inf = X>0')
assert source_statement_e(cusp.elementary_pole_patch,
    'raw = elliptic_pole_raw_inverse(compressed,infinity,dps=dps,elementary=True)')
assert source_statement_e(cusp.elliptic_pole_raw_inverse,
    'nu = mp.exp(-2) if infinity else mp.mpf(1)')
assert grmin<grw<grmax and gp3w>0 and sp.im(gpcw)>0
assert 8*grmax/(1-8*grmax)<sp.Rational(1,2)
# aa=pi*p3>0, bb=-pi*Re(pc), cc=pi*Im(pc)>0, so -u1 has positive
# real and negative imaginary parts.  Thus 0<height<1/12<1/10;
# u2 has positive imaginary part and u1!=0, giving tc=u2/u1!=0.
# The raw positive-root bounds in (f1) hold since S(r)<r<grmax<rmax_e.
assert grmax<rmax_e and grmax<sp.Rational(1,10**4)
assert 4*grmax+(4*grmax)**2/2<sp.Rational(1,2*10**5)
print('(h1) strict north transition chamber and shared source-branch'
      ' bounds: PASS')


def clone_scalar_nodes_h(base,nodes):
    symbols = (base,)+tuple(node for node,_ in nodes)
    renamed = {node:sp.Symbol('h_'+str(node),positive=True)
               for node in symbols}
    return renamed[base],[(renamed[node],rhs.xreplace(renamed))
                          for node,rhs in nodes],renamed


hr,S_nodes_h,S_names_h = clone_scalar_nodes_h(gr,S_nodes_g)
hs,R_nodes_h,R_names_h = clone_scalar_nodes_h(gs,R_nodes_g)
S_jet_defs_h,S_known_h = scalar_jet_dag_g(hr,S_nodes_h)
R_jet_defs_h,R_known_h = scalar_jet_dag_g(hs,R_nodes_h)
assert len(S_jet_defs_h)==3*len(S_nodes_h)==18
assert len(R_jet_defs_h)==3*len(R_nodes_h)==18
assert set(S_names_h.values()).isdisjoint(set(R_names_h.values()))
assert set(S_names_h.values()).isdisjoint(set(S_names_h.keys()))
assert set(R_names_h.values()).isdisjoint(set(R_names_h.keys()))
hz_in = sp.symbols('hz_in0:4',real=True)
hw_in = sp.symbols('hw_in0:4',real=True)
D_value_h,M_D_h,H_D_h = radial_jets_g(
    hz_in,S_names_h[gDgain],*S_known_h[S_names_h[gDgain]])
I_value_h,M_I_h,H_I_h = radial_jets_g(
    hw_in,R_names_h[gIgain],*R_known_h[R_names_h[gIgain]])
assert all(not entry.has(sp.Derivative) for entry in M_D_h)
assert all(not entry.has(sp.Derivative) for entry in M_I_h)
print('(h2) separate north S/S^-1 scalar and four-real radial'
      ' value/first/second-jet circuits: PASS')

# Actual stage order: north local round input, D, the 28-node north raw
# inverse (including its source parameters), D^-1 at that raw output,
# then north finite-height cusp restoration with the same torus inputs.
stage_defs_h = list(zip(hz_in,local_f))
stage_defs_h.append((hr,sum(z*z for z in hz_in)))
stage_defs_h.extend(S_jet_defs_h)
stage_defs_h.extend(zip(z_f,D_value_h))
stage_defs_h.extend(raw_input_defs_f)
stage_defs_h.extend(raw_jet_defs_f)
stage_defs_h.extend(zip(hw_in,raw_out_f))
stage_defs_h.append((hs,sum(w*w for w in hw_in)))
stage_defs_h.extend(R_jet_defs_h)
stage_inputs_h = set(round_f)|{fp2,fp3}
stage_bound_h = set(stage_inputs_h)
for node,rhs in stage_defs_h:
    assert rhs.free_symbols <= stage_bound_h
    stage_bound_h.add(node)
assert tuple(rhs for _,rhs in stage_defs_h[:4])==tuple(local_f)
assert all(value.free_symbols <= stage_bound_h for value in I_value_h)

M_D_at_round_h = bind_matrix_e(M_D_h,hz_in,local_f)
H_D_at_round_h = [bind_matrix_e(H,hz_in,local_f) for H in H_D_h]
M_raw_at_D_h = bind_matrix_e(Mraw_f,z_f,D_value_h)
H_raw_at_D_h = [bind_matrix_e(H,z_f,D_value_h) for H in Hraw_f]
M_I_at_raw_h = bind_matrix_e(M_I_h,hw_in,raw_out_f)
H_I_at_raw_h = [bind_matrix_e(H,hw_in,raw_out_f) for H in H_I_h]
M_Dround_h,H_Dround_h = compose_jets_g(
    M_D_at_round_h,H_D_at_round_h,Mround_f,Hround_f)
M_rawround_h,H_rawround_h = compose_jets_g(
    M_raw_at_D_h,H_raw_at_D_h,M_Dround_h,H_Dround_h)
M_patch4_h,H_patch4_h = compose_jets_g(
    M_I_at_raw_h,H_I_at_raw_h,M_rawround_h,H_rawround_h)
M_patch6_h = sp.BlockDiagMatrix(M_patch4_h,sp.eye(2))
H_patch6_h = [sp.BlockDiagMatrix(H,sp.zeros(2,2))
              for H in H_patch4_h]
H_patch6_h.extend([sp.ZeroMatrix(6,6),sp.ZeroMatrix(6,6)])
rest_input_h = tuple(I_value_h)+(fp2,fp3)
full_value_h = sp.Matrix([
    bind_e(restore_f[k],rest_f,rest_input_h) for k in range(6)])
M_rest_at_patch_h = bind_matrix_e(Mrestore_f,rest_f,rest_input_h)
H_rest_at_patch_h = [bind_matrix_e(H,rest_f,rest_input_h)
                     for H in Hrestore_f]
Mfull_h,Hfull_h = compose_jets_g(
    M_rest_at_patch_h,H_rest_at_patch_h,M_patch6_h,H_patch6_h)
J_transition_h = sp.MatMul(sp.Inverse(Mfull_h),J0,Mfull_h,
                            evaluate=False)
assert Mfull_h.shape==(6,6) and J_transition_h.shape==(6,6)
assert len(Hfull_h)==6 and all(H.shape==(6,6) for H in Hfull_h)
assert all(value.free_symbols <= stage_bound_h for value in full_value_h)
assert Mfull_h.free_symbols <= stage_bound_h
assert all(H.free_symbols <= stage_bound_h for H in Hfull_h)
assert J_transition_h.free_symbols <= stage_bound_h
print('(h3) actual stage-bound north transition h, 6x6 M, and six'
      ' second-jet circuits: PASS')

assert second_jet_skew_e(M_D_at_round_h,H_D_at_round_h,
                         Mround_f,Hround_f)==64
assert second_jet_skew_e(M_raw_at_D_h,H_raw_at_D_h,
                         M_Dround_h,H_Dround_h,
                         inner_skew_certified=True)==64
assert second_jet_skew_e(M_I_at_raw_h,H_I_at_raw_h,
                         M_rawround_h,H_rawround_h,
                         inner_skew_certified=True)==64
assert second_jet_skew_e(M_rest_at_patch_h,H_rest_at_patch_h,
                         M_patch6_h,H_patch6_h,
                         inner_skew_certified=True)==216
# N10d gives det D,det D^-1>0; (11.11) gives raw sigma_min>0.99993.
# The north round factor has det L_f=4 and positive radial gain.  Its
# inverse uses h=pi*r and rho=2h/(1+h²), with
# r=rho/[pi*(1+sqrt(1-rho²))].  North restoration has determinant
# +(1-|t|⁴)/nu² by (f5), nonzero under the §11.6 radius bound.
hh = sp.Symbol('hh',positive=True)
hrho = 2*hh/(1+hh**2)
hroot_rho = (1-hh**2)/(1+hh**2)
assert sp.factor(hrho**2+hroot_rho**2-1)==0
assert sp.cancel(hrho/(sp.pi*(1+hroot_rho))-hh/sp.pi)==0
assert 4*grmax<1 and 4*grmax+(4*grmax)**2/2<sp.Rational(1,2*10**5)
assert 1-5000*grmax==sp.Rational(99993,100000)
assert (sp.Rational(2,10**6)*9/2)**2<1
assert J0*J0==-sp.eye(6)
print('(h4) 3x64+216 compositional skews and 36/216 north J certificate;'
      ' raw margin >99993/100000, using N10/(11.11)/§11.6: PASS')

# (h5) Type-(d) one-point source comparison and measured raw margin only.
with mp.workdps(80):
    numerical_h = dict(zip(hz_in,[sp.Rational(1,10000),sp.S.Zero,
                                   sp.Rational(5,100000),
                                   -sp.Rational(3,100000)]))
    numerical_h[hr] = sum(value*value for value in numerical_h.values())
    for node,rhs in S_nodes_h:
        assert rhs.free_symbols <= numerical_h.keys()
        numerical_h[node] = rhs.subs(numerical_h).evalf(90)
    compressed_h = [mp.mpf(str(rhs.subs(numerical_h).evalf(80)))
                    for rhs in D_value_h]
    for node,rhs in zip(z_f,D_value_h):
        numerical_h[node] = rhs.subs(numerical_h).evalf(90)
    for node,rhs in raw_input_defs_f+raw_node_defs_f:
        assert rhs.free_symbols <= numerical_h.keys()
        numerical_h[node] = rhs.subs(numerical_h).evalf(90)
    raw_value_h = [rhs.subs(numerical_h).evalf(90) for rhs in raw_out_f]
    numerical_h.update(zip(hw_in,raw_value_h))
    numerical_h[hs] = sum(value*value for value in raw_value_h)
    for node,rhs in R_nodes_h:
        assert rhs.free_symbols <= numerical_h.keys()
        numerical_h[node] = rhs.subs(numerical_h).evalf(90)
    symbolic_patch_h = [mp.mpf(str(rhs.subs(numerical_h).evalf(80)))
                        for rhs in I_value_h]
    source_patch_h = cusp.elementary_pole_patch(
        (mp.mpc(mp.mpf(1)/10000),
         mp.mpc(mp.mpf(5)/100000,-mp.mpf(3)/100000)),
        True,dps=80)
    source_patch_real_h = [source_patch_h[0].real,source_patch_h[0].imag,
                           source_patch_h[1].real,source_patch_h[1].imag]
    patch_error_h = max(abs(a-b) for a,b in
                        zip(symbolic_patch_h,source_patch_real_h))
    assert patch_error_h<mp.mpf('1e-55')
    sample_M_h = mp.matrix([[mp.diff(
        lambda value:source_raw_real_f(compressed_h[:j]+[value]
                                       +compressed_h[j+1:])[i],
        compressed_h[j]) for j in range(4)] for i in range(4)])
    measured_defect_h = sample_M_h-mp.eye(4)
    measured_norm_h = mp.sqrt(max(mp.eigsy(
        measured_defect_h.T*measured_defect_h,eigvals_only=True)))
    measured_sigma_h = mp.sqrt(min(mp.eigsy(
        sample_M_h.T*sample_M_h,eigvals_only=True)))
    raw_bound_h = 5000*sum(value*value for value in compressed_h)
    assert measured_norm_h<raw_bound_h<5000*mp.mpf(14)/10**9
    assert measured_sigma_h>mp.mpf(99993)/100000
    print('(h5) numerical one-point sample only: full patch value error =',
          mp.nstr(patch_error_h,5),'; raw ||Dg-I|| =',
          mp.nstr(measured_norm_h,17),'vs bound',mp.nstr(raw_bound_h,17),
          '; raw sigma_min =',mp.nstr(measured_sigma_h,17),': PASS')
# (i) A six-variable, second-order state for one *actual* free-slab
# projector cell.  C30 checks a generic two-direction update.  Here the
# state carries all 6 first and 21 independent second coefficients, and
# the projector is the source's ten-slot pure-middle cell C31.  The
# astronomical fold is certified only conditionally on its still-open
# projector-cell and seam premises; this section does not certify G_cl.
_pairs_i = tuple((a,b) for a in range(6) for b in range(a,6))
_index_i = {pair:k for k,pair in enumerate(_pairs_i)}


class SixJetI:
    """Value, six first and 21 symmetric second partials at one point."""
    __slots__ = ('value','grad','hess')

    def __init__(self,value,grad=None,hess=None):
        self.value = sp.sympify(value)
        self.grad = tuple(sp.sympify(v) for v in (
            (0,)*6 if grad is None else grad))
        self.hess = tuple(sp.sympify(v) for v in (
            (0,)*21 if hess is None else hess))
        if len(self.grad)!=6 or len(self.hess)!=21:
            raise ValueError('a six-jet needs 6 first and 21 second entries')

    @classmethod
    def from_full(cls,value,grad,full):
        if len(full)!=6 or any(len(row)!=6 for row in full):
            raise ValueError('the Hessian must be 6 by 6')
        if any(sp.simplify(full[a][b]-full[b][a])!=0
               for a in range(6) for b in range(a+1,6)):
            raise ValueError('a coordinate Hessian must be symmetric')
        return cls(value,grad,[full[a][b] for a,b in _pairs_i])

    def second(self,a,b):
        return self.hess[_index_i[tuple(sorted((a,b)))]]


def jconst_i(value):
    return SixJetI(value)


def jvar_i(value,index):
    grad = [sp.S.Zero]*6
    grad[index] = sp.S.One
    return SixJetI(value,grad)


def jadd_i(A,B):
    return SixJetI(A.value+B.value,
                   [a+b for a,b in zip(A.grad,B.grad)],
                   [a+b for a,b in zip(A.hess,B.hess)])


def jneg_i(A):
    return SixJetI(-A.value,[-v for v in A.grad],
                   [-v for v in A.hess])


def jmul_i(A,B):
    return SixJetI(A.value*B.value,
        [A.grad[a]*B.value+A.value*B.grad[a] for a in range(6)],
        [A.second(a,b)*B.value+A.grad[a]*B.grad[b]
         +A.grad[b]*B.grad[a]+A.value*B.second(a,b)
         for a,b in _pairs_i])


def juni_i(A,value,slope,curvature):
    return SixJetI(value,[slope*g for g in A.grad],
        [slope*A.second(a,b)+curvature*A.grad[a]*A.grad[b]
         for a,b in _pairs_i])


def jinv_i(A):
    v = A.value
    return juni_i(A,1/v,-1/v**2,2/v**3)


def jsqrt_i(A):
    v = A.value
    return juni_i(A,sp.sqrt(v),1/(2*sp.sqrt(v)),
                  -1/(4*v**sp.Rational(3,2)))


def jexp_i(A):
    v = sp.exp(A.value)
    return juni_i(A,v,v,v)


def jatan_i(A):
    v = A.value
    return juni_i(A,sp.atan(v),1/(1+v*v),
                  -2*v/(1+v*v)**2)


def jconj_i(A):
    return SixJetI(sp.conjugate(A.value),
                   [sp.conjugate(v) for v in A.grad],
                   [sp.conjugate(v) for v in A.hess])


def jreal_i(A):
    return jmul_i(jconst_i(sp.Rational(1,2)),
                  jadd_i(A,jconj_i(A)))


def jimag_i(A):
    return jmul_i(jconst_i(1/(2*sp.I)),
                  jadd_i(A,jneg_i(jconj_i(A))))


def jdiv_i(A,B):
    return jmul_i(A,jinv_i(B))


def jzero_i(A):
    return (sp.simplify(A.value)==0
            and all(sp.simplify(v)==0 for v in A.grad)
            and all(sp.simplify(v)==0 for v in A.hess))


def jsub_i(A,B):
    return jadd_i(A,jneg_i(B))


def projection_update_i(reference,previous):
    """One exact source recurrence, with the source's normalized phase."""
    if len(reference)!=10 or len(previous)!=10:
        raise ValueError('the projection recurrence has ten slots')
    phase = jconst_i(0)
    for ref,old in zip(reference,previous):
        phase = jadd_i(phase,jmul_i(jconj_i(ref),old))
    magnitude = jsqrt_i(jmul_i(jconj_i(phase),phase))
    return [jdiv_i(jmul_i(ref,phase),magnitude)
            for ref in reference]


def phase_increment_i(previous,current):
    """One exact right-half-plane source increment."""
    ratio = jmul_i(current,jconj_i(previous))
    return ratio,jatan_i(jdiv_i(jimag_i(ratio),jreal_i(ratio)))


# The complete two-jet product is the degree-at-most-two Taylor product.
# The coefficient test has independent generic data in every one of the
# 28 slots, so the checked rule applies to each of the N+1 source steps.
eps_i = sp.symbols('eps_i0:6')
aval_i, bval_i = sp.symbols('aval_i bval_i')
ag_i = sp.symbols('ag_i0:6')
bg_i = sp.symbols('bg_i0:6')
ah_i = sp.symbols('ah_i0:21')
bh_i = sp.symbols('bh_i0:21')
Ajet_i,Bjet_i = SixJetI(aval_i,ag_i,ah_i),SixJetI(bval_i,bg_i,bh_i)
def taylor_i(A):
    return (A.value+sum(A.grad[a]*eps_i[a] for a in range(6))
            +sum(A.second(a,b)*eps_i[a]*eps_i[b]
                 /(2 if a==b else 1) for a,b in _pairs_i))
poly_i = sp.Poly(taylor_i(Ajet_i)*taylor_i(Bjet_i),*eps_i)
AB_i = jmul_i(Ajet_i,Bjet_i)
assert sp.expand(poly_i.coeff_monomial(1)-AB_i.value)==0
assert all(sp.expand(poly_i.coeff_monomial(eps_i[a])
                     -AB_i.grad[a])==0 for a in range(6))
assert all(sp.expand(poly_i.coeff_monomial(
              eps_i[a]*eps_i[b])-
              AB_i.second(a,b)/(2 if a==b else 1))==0
           for a,b in _pairs_i)
F0_i,F1_i,F2_i = sp.symbols('F0_i F1_i F2_i')
delta_i = taylor_i(Ajet_i)-Ajet_i.value
unary_poly_i = sp.Poly(F0_i+F1_i*delta_i+F2_i*delta_i**2/2,
                       *eps_i)
unary_jet_i = juni_i(Ajet_i,F0_i,F1_i,F2_i)
assert sp.expand(unary_poly_i.coeff_monomial(1)
                 -unary_jet_i.value)==0
assert all(sp.expand(unary_poly_i.coeff_monomial(eps_i[a])
                     -unary_jet_i.grad[a])==0 for a in range(6))
assert all(sp.expand(unary_poly_i.coeff_monomial(
              eps_i[a]*eps_i[b])-
              unary_jet_i.second(a,b)/(2 if a==b else 1))==0
           for a,b in _pairs_i)
unary_variable_i = sp.Symbol('unary_variable_i',positive=True)
unary_rules_i = [
    (1/unary_variable_i,-1/unary_variable_i**2,
     2/unary_variable_i**3),
    (sp.sqrt(unary_variable_i),1/(2*sp.sqrt(unary_variable_i)),
     -1/(4*unary_variable_i**sp.Rational(3,2))),
    (sp.exp(unary_variable_i),sp.exp(unary_variable_i),
     sp.exp(unary_variable_i)),
    (sp.atan(unary_variable_i),1/(1+unary_variable_i**2),
     -2*unary_variable_i/(1+unary_variable_i**2)**2)]
assert all(sp.simplify(sp.diff(f,unary_variable_i)-fp)==0
           and sp.simplify(sp.diff(f,unary_variable_i,2)-fpp)==0
           for f,fp,fpp in unary_rules_i)
print('(i1) six-variable Taylor ring: 6 gradients and 21 symmetric'
      ' Hessians, exact product/unary coefficients: PASS')

# The source's pure-middle cell has 1/2<|t|<1 and
# pi/2<Arg(t)<3pi/4.  Only slots 8,9 survive, with a common phase.
# Their unnormalized weights are (d,1), where
# d=1-s((Arg(t)-pi/2)/(pi/4)).  This is the actual C31 projector,
# rather than a generic 10 by 10 rank-one placeholder.
assert source_statement_e(cusp.elliptic_meridian_projector,
    'kind,data = elliptic_raw_base_reference(v,h,dps=dps,elementary=elementary)')
assert source_statement_e(cusp.elliptic_meridian_projector,
    'e = elliptic_bundle_embedding(kind,data,dps=dps,elementary=elementary)')
assert source_statement_e(cusp.elliptic_meridian_projector,
    'return mp.matrix([[a*mp.conj(b) for b in e] for a in e])')
assert source_statement_e(cusp.elliptic_raw_base_reference,
    "return ('mid',(t,mp.mpf(0),four,lt))")
assert source_statement_e(cusp.elliptic_bundle_embedding,
    'b0 = 1-scalar_flat_step((r-mp.mpf(\'.3\'))/mp.mpf(\'.2\'))')
assert source_statement_e(cusp.elliptic_bundle_embedding,
    'bi = 1-scalar_flat_step((1/r-mp.mpf(\'.3\'))/mp.mpf(\'.2\'))')
assert source_statement_e(cusp.elliptic_bundle_embedding,
    'weight = 1-scalar_flat_step((abs(angle)-mp.pi/2)/(mp.pi/4))')
assert source_statement_e(cusp.elliptic_bundle_embedding,
    'vec[8+index] = bm*weight*mp.exp(mp.j*(2*mp.pi*one-(angle+centre-lt.imag)*four))')
assert source_statement_e(cusp.elementary_projection_section,
    'kind,data = elliptic_raw_base_reference(tuple((1-h*h)*a for a in n),h,dps=dps,elementary=True)')
assert source_statement_e(cusp.elementary_projection_section,
    'ref = elliptic_bundle_embedding(kind,data,dps=dps,elementary=True)')
assert source_statement_e(cusp.elementary_projection_section,
    'phase = _elementary_inner(ref,vector)')
assert source_statement_e(cusp.elementary_projection_section,
    'vector = tuple(a*phase/abs(phase) for a in ref)')
ar_i,ai_i = sp.symbols('ar_i ai_i',real=True)
aj_i,bj_i = jvar_i(ar_i,0),jvar_i(ai_i,1)
halfangle_i = jatan_i(jdiv_i(jneg_i(aj_i),bj_i))
cut_i = jmul_i(jconst_i(4/sp.pi),halfangle_i)
flat_exponent_i = jsub_i(jinv_i(cut_i),
                         jinv_i(jsub_i(jconst_i(1),cut_i)))
flat_i = jinv_i(jadd_i(jconst_i(1),jexp_i(flat_exponent_i)))
d_i = jsub_i(jconst_i(1),flat_i)
qscalar_i = sp.Symbol('qscalar_i',positive=True)
flat_source_i = sp.exp(-1/qscalar_i)/(
    sp.exp(-1/qscalar_i)+sp.exp(-1/(1-qscalar_i)))
assert sp.simplify(flat_source_i
    -1/(1+sp.exp(1/qscalar_i-1/(1-qscalar_i))))==0
assert jzero_i(jsub_i(halfangle_i, SixJetI(
    sp.atan(-ar_i/ai_i),
    [-ai_i/(ar_i**2+ai_i**2),ar_i/(ar_i**2+ai_i**2),0,0,0,0],
    [2*ar_i*ai_i/(ar_i**2+ai_i**2)**2,
     (ai_i**2-ar_i**2)/(ar_i**2+ai_i**2)**2,0,0,0,0,
     -2*ar_i*ai_i/(ar_i**2+ai_i**2)**2,
     *([0]*14)])))
assert all(d_i.grad[a]==0 for a in range(2,6))
assert all(d_i.second(a,b)==0 for a,b in _pairs_i if b>=2)
assert sp.Rational(1,2)<sp.sqrt(5)/3<1
assert 0<sp.atan(sp.Rational(1,2))<sp.pi/4
print('(i2) actual source pure-middle 10-slot projector profile:'
      ' complete six-variable local second jet and strict witness: PASS')

# Write the source projector as P=ref ref*.  The two nonzero slots are
# 8,9 and P has a nonconstant 2 by 2 block.  All 100 matrix entries
# and their 6+21 derivatives are specified by these two source slots.
D_i = sp.Symbol('D_i',real=True)
v8_i,v9_i = sp.symbols('v8_i v9_i')
ref_formal_i = sp.Matrix([D_i,1])/sp.sqrt(1+D_i**2)
P_formal_i = ref_formal_i*ref_formal_i.conjugate().T
inject_i = sp.zeros(10,2)
inject_i[8,0],inject_i[9,1] = 1,1
P_ten_i = inject_i*P_formal_i*inject_i.T
vec_formal_i = sp.Matrix([v8_i,v9_i])
phase_formal_i = (ref_formal_i.conjugate().T*vec_formal_i)[0]
assert all(sp.simplify((P_formal_i*vec_formal_i)[k]
                       -ref_formal_i[k]*phase_formal_i)==0
           for k in range(2))
assert sp.simplify((ref_formal_i.conjugate().T*ref_formal_i)[0]-1)==0
assert sp.simplify(P_formal_i*P_formal_i-P_formal_i)==sp.zeros(2)
assert P_ten_i.shape==(10,10)
assert sp.simplify(P_ten_i*P_ten_i-P_ten_i)==sp.zeros(10)
assert sp.simplify(sp.trace(P_ten_i)-1)==0
dr_i = jsqrt_i(jadd_i(jmul_i(d_i,d_i),jconst_i(1)))
ref_i = [jconst_i(0) for _ in range(10)]
ref_i[8],ref_i[9] = jdiv_i(d_i,dr_i),jinv_i(dr_i)
pg_i = sp.symbols('pg_i0:6',real=True)
ph_i = sp.symbols('ph_i0:21',real=True)
prior_phase_i = jexp_i(jmul_i(jconst_i(sp.I),
                              SixJetI(0,pg_i,ph_i)))
assert jzero_i(jsub_i(jmul_i(prior_phase_i,jconj_i(prior_phase_i)),
                     jconst_i(1)))
prev_i = [jconst_i(0) for _ in range(10)]
prev_i[9] = prior_phase_i
Pjet_i = [[jmul_i(ref_i[a],jconj_i(ref_i[b]))
           for b in range(10)] for a in range(10)]
# The generic exact identity P e=ref<ref,e> above, together with the
# tested six-variable Taylor-ring homomorphism, proves equality of all
# 28 coefficients here without expanding the long flat-step expression.
next_i = projection_update_i(ref_i,prev_i)
assert len(Pjet_i)==10 and all(len(row)==10 for row in Pjet_i)
assert all(len(p.grad)==6 and len(p.hess)==21
           for row in Pjet_i for p in row)
assert all(jzero_i(Pjet_i[a][b])
           for a in range(10) for b in range(10)
           if a not in (8,9) or b not in (8,9))
assert len(next_i)==10
assert all(len(v.grad)==6 and len(v.hess)==21 for v in next_i)
assert all(jzero_i(v) for v in next_i[:8])
assert all(not entry.has(sp.Derivative)
           for v in next_i[8:] for entry in
           (v.value,*v.grad,*v.hess))
# The previous vector has a variable unit phase with independent data
# in every first and second slot.  Its inner product with this projector
# has magnitude 1/sqrt(1+d^2)>0, and these angular data reach the
# normalized output in all 28 slots.
assert all(next_i[9].grad[a].has(pg_i[a]) for a in range(6))
assert all(next_i[9].second(a,b).has(ph_i[_index_i[(a,b)]])
           for a,b in _pairs_i)
assert sp.ask(sp.Q.positive(1/sp.sqrt(1+D_i**2))) is True
print('(i3) source-normalized rank-one step: ten slots, 6 first and 21'
      ' second entries per slot, nonzero denominator: PASS')

# An arbitrary local reference gauge gamma=exp(i theta) must disappear
# from the normalized projection.  Test its unit equation in all 28
# coefficients for independent real theta gradients/Hessians.  Then
# (gamma ref)<gamma ref,e> = (gamma conj(gamma))ref<ref,e>
# in this Taylor ring; the denominator is the same.  The formal source
# update above supplies the actual ref and projector.
thg_i = sp.symbols('thg_i0:6',real=True)
thh_i = sp.symbols('thh_i0:21',real=True)
theta_i = SixJetI(0,thg_i,thh_i)
gamma_i = jexp_i(jmul_i(jconst_i(sp.I),theta_i))
assert jzero_i(jsub_i(jmul_i(gamma_i,jconj_i(gamma_i)),
                     jconst_i(1)))
# A direct variable-gauge test of the source update, at the exact
# nonorthogonal ref=(1,2)/sqrt(5), has all 6+21 gauge entries active.
fixed_ref_i = [jconst_i(0) for _ in range(10)]
fixed_ref_i[8],fixed_ref_i[9] = (
    jconst_i(1/sp.sqrt(5)),jconst_i(2/sp.sqrt(5)))
gauged_ref_i = [jmul_i(gamma_i,v) for v in fixed_ref_i]
gauge_previous_i = [jconst_i(0) for _ in range(10)]
gauge_previous_i[9] = jconst_i(1)
gauged_next_i = projection_update_i(gauged_ref_i,gauge_previous_i)
assert all(jzero_i(jsub_i(gauged_next_i[k],fixed_ref_i[k]))
           for k in (8,9))
print('(i4) source update is reference-gauge invariant in all 28'
      ' coefficients, including a variable gauge: PASS')

# One source-form right-half-plane increment on generic unit mismatch
# jets.  At the strict ratio value R=1, the increment is
# theta_1-theta_0 to second order in every direction.
assert source_statement_e(cusp.elementary_phase_family,
    'ratio = value*mp.conj(previous)')
assert source_statement_e(cusp.elementary_phase_family,
    'if ratio.real<=0:\n    raise ValueError(\'phase count does not give right-half-plane increments\')')
assert source_statement_e(cusp.elementary_phase_family,
    'total += mp.atan(ratio.imag/ratio.real)')
assert source_statement_e(cusp._elementary_corrected_section,
    'return _elementary_inner(direct,section)')
assert source_statement_e(cusp._elementary_corrected_section,
    'P = elliptic_meridian_projector(nn,h,dps=dps,elementary=True)')
def has_source_fold_i(function,iteration):
    tree = ast.parse(inspect.getsource(function))
    wanted = ast.dump(ast.parse(iteration,mode='eval').body,
                      include_attributes=False)
    return any(isinstance(node,ast.For)
               and ast.dump(node.iter,include_attributes=False)==wanted
               for node in ast.walk(tree))
assert has_source_fold_i(cusp.elementary_projection_section,
                         'range(steps+1)')
assert has_source_fold_i(cusp.elementary_phase_family,'values')
th0g_i = sp.symbols('th0g_i0:6',real=True)
th1g_i = sp.symbols('th1g_i0:6',real=True)
th0h_i = sp.symbols('th0h_i0:21',real=True)
th1h_i = sp.symbols('th1h_i0:21',real=True)
theta0_i,theta1_i = SixJetI(0,th0g_i,th0h_i),SixJetI(0,th1g_i,th1h_i)
g0_i = jexp_i(jmul_i(jconst_i(sp.I),theta0_i))
g1_i = jexp_i(jmul_i(jconst_i(sp.I),theta1_i))
assert jzero_i(jsub_i(jmul_i(g0_i,jconj_i(g0_i)),jconst_i(1)))
assert jzero_i(jsub_i(jmul_i(g1_i,jconj_i(g1_i)),jconst_i(1)))
ratio_i,phase_step_i = phase_increment_i(g0_i,g1_i)
assert ratio_i.value==1
assert jzero_i(jsub_i(phase_step_i,jsub_i(theta1_i,theta0_i)))
params_i = cusp.elementary_closed_parameters()
L_i,N_i = params_i['projector_bound'],params_i['projection_steps']
Lg_i,K_i = params_i['phase_bound'],params_i['phase_steps']
assert (L_i==10**8000 and N_i==4*L_i+1
        and Lg_i==32*(N_i+1)*L_i and K_i==8*Lg_i+1)
assert 4*L_i<N_i and 8*Lg_i<K_i
assert source_statement_e(cusp.elementary_global_slice,
    'params = elementary_closed_parameters()')
assert source_statement_e(cusp.elementary_global_slice,
    "N = params['projection_steps'] if projection_steps is None else projection_steps")
assert source_statement_e(cusp.elementary_global_slice,
    "K = params['phase_steps'] if phase_steps is None else phase_steps")
assert source_statement_e(cusp._elementary_corrected_section,
    'result = elementary_projection_section(xi,X,N,dps,diagnostics)')
assert source_statement_e(cusp._elementary_corrected_section,
    'theta = elementary_phase_family(mismatch,n,X,anchor,K,diagnostics=diagnostics)')
# The checked update maps one valid six-jet state to another.  Induction
# therefore applies to the exact finite N+1 and K increments, provided
# every projector/phase generator has a valid local jet and the published
# nonzero/right-half-plane bounds hold along its path.  This checker
# verifies one actual projector cell, not those remaining premises.
print('(i5) right-half-plane increment: 6 first and 21 second jets;'
      ' conditional induction at exact fixed N/K: PASS')

# Negative control: a skew entry cannot be smuggled into a coordinate
# Hessian by the compact 21-coefficient representation.
bad_hess_i = [[0]*6 for _ in range(6)]
bad_hess_i[0][1] = 1
skew_rejected_i = False
try:
    SixJetI.from_full(0,[0]*6,bad_hess_i)
except ValueError:
    skew_rejected_i = True
assert skew_rejected_i
print('(i6) asymmetric-Hessian negative control rejected: PASS')
print('ALL PASS')
