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
C7  discrete-Gaussian centring and reflection identities used to remove wrap from the fibre blend;
    the infinite-series argument and the strict half-cell bound are in programme section 3.3.
C8  closed inverses in the two linear fixed-point balls, including the poles.
C9  full complex coframe, determinant and tensor formula on those target caps.
C10 naturality in dimension six, using a completely arbitrary symmetric second jet.
C11 explicit circle-connection identities for the proposed quadrature lift.
C12 branch and deck invariance of the existing T^2 trivialization.
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

# ---------------------------------------------------------------- C7: no wrap is needed
j, centre, a, b = sp.symbols('j centre a b', real=True)
rep("C7a  torus exponent = P(j-centre)^2 - P centre^2, centre=rho/P-1/2",
    sp.expand((P * (j - centre) ** 2 - P * centre ** 2).subs(centre, rho / P - sp.Rational(1, 2))
              - (-2 * j * rho + j * (j + 1) * P)) == 0)
weight = sp.exp(-a * (j - centre) ** 2)
rep("C7b  at a half-integer centre b, reflection j -> 2b-j preserves weights",
    sp.simplify(weight.subs({j: 2 * b - j, centre: b}, simultaneous=True) - weight.subs(centre, b)) == 0)
rep("C7c  reflected contributions to the centred first moment cancel",
    sp.expand((j - b) + (2 * b - j - b)) == 0)
D, N, S = sp.symbols('D N S', positive=True)
rep("C7d  d mean/d centre = 2a Var (a=-P>0), by the quotient rule",
    sp.cancel((2 * a * (S - centre * N) * D - N * 2 * a * (N - centre * D)) / D ** 2
              - 2 * a * (S / D - (N / D) ** 2)) == 0)
rep("C7e  integer translation of centre is the same translation of the index",
    sp.expand(((j + 1) - (centre + 1)) ** 2 - (j - centre) ** 2) == 0)
# For a>0 the Gaussian and all derivatives converge uniformly on compact sets.
# Thus mean'(centre)>0 and mean(b)=b at every b in (1/2)Z.
# Centre and mean lie in the same half-cell, so |mean-centre|<1/2.
# The node correction is absent whenever the blend weight is nonzero:
# |t| <= rad/(2 sqrt(a0)) <= exp(2)/40 < 1/5 < RA=3/10 on its support.
rep("C7f  exp(2)<8 suffices to separate node support from the blend annulus",
    sp.Rational(8, 40) < sp.Rational(3, 10) and bool(sp.exp(2) < 8))
theta, xi = sp.symbols('theta xi', real=True)
aa1 = x4 - xi - sp.Rational(1, 2)
aa2 = -x4 - sp.Rational(1, 2)
aa2p = -aa2 + theta / (2 * sp.pi) - xi
rep("C7g  clutching: a2prime-a1 = arg(t)/(2pi) modulo integers",
    sp.expand(aa2p - aa1 - theta / (2 * sp.pi)) == 1)
delta = sp.Symbol('delta', real=True)
rep("C7h  finite centred truncation: square differences for pairing j,-j and j,1-j",
    sp.expand((-j - delta) ** 2 - (j - delta) ** 2 - 4 * j * delta) == 0
    and sp.expand((j - delta) ** 2 - (1 - j - delta) ** 2 - (2 * j - 1) * (1 - 2 * delta)) == 0)

# ---------------------------------------------------------------- C8: the inverse avoids polar coordinates at Z=0
U0 = sp.Matrix([[1, 1], [sp.I, -sp.I]]) / sp.sqrt(2)
Ui = sp.Matrix([[1, 1], [1, -1]]) / sp.sqrt(2)
rep("C8a  both linear lift matrices have the explicit inverse U*",
    all(sp.simplify(U.conjugate().T * U) == sp.eye(2) for U in (U0, Ui)))
Z1, Z2 = sp.symbols('Z1 Z2')
scale = sp.Symbol('scale', positive=True)
for name, U, second in (('c0', U0, -sp.I * Z2), ('ci', Ui, Z2)):
    recovered = sp.Matrix([Z1 + second, Z1 - second]) / sp.sqrt(2 * scale)
    rep(f"C8b  {name}: forward linear lift composed with inverse is (Z1,Z2)",
        sp.simplify(sp.sqrt(scale) * U * recovered - sp.Matrix([Z1, Z2])) == sp.zeros(2, 1))
xs, ys, exponent = sp.symbols('xs ys exponent')
xp, yp = xs * sp.exp(-exponent), ys * sp.exp(exponent)
rep("C8c  undoing the T^2 slice preserves t and recovers (xs,ys)",
    sp.simplify(xp * yp - xs * ys) == 0
    and sp.simplify(xp * sp.exp(exponent) - xs) == 0
    and sp.simplify(yp * sp.exp(-exponent) - ys) == 0)
phi2, phi3 = sp.symbols('phi2 phi3', real=True)
rep("C8d  C=phi+H c recovers the prescribed torus phases after the slice inverse",
    phi0.subs({C2_: phi2 + H * c2, C3_: phi3 + H * c3}) == sp.Matrix([phi2, phi3]))
rep("C8e  in the linear balls, |t| <= exp(2)*EPS1/2 < 2/125 < 1",
    sp.Rational(8, 1) * sp.Rational(4, 1000) / 2 == sp.Rational(2, 125)
    and sp.Rational(2, 125) < 1)
J_slice = sp.diag(-sp.I, sp.I)  # coordinates (sqrt(a0)*conj(x), y)
rep("C8f  pole values on complex Z: south (-Z2,Z1), north (-i Z2,-i Z1)",
    sp.simplify(U0 * J_slice * U0.conjugate().T) == sp.Matrix([[0, -1], [1, 0]])
    and sp.simplify(Ui * J_slice * Ui.conjugate().T) == sp.Matrix([[0, -sp.I], [-sp.I, 0]]))
Jw = sp.Matrix([[0, -1], [1, 0]])
Ci = sp.diag(-1, 1)  # w=-phi2+i phi3 at the north pole
rep("C8g  pole values on torus angles: south Jw, north -Jw", Ci.inv() * Jw * Ci == -Jw)

# ---------------------------------------------------------------- C9: an explicit tensor on the whole caps
tr, ti = sp.symbols('tr ti', real=True)
tc = tr+sp.I*ti
r2 = tr*tr+ti*ti
Fh = 2*r2/(1-r2)
for sigma in (1,-1):
    e2c, e3c = sigma*(1-tc**2), sp.I*(1+tc**2)
    c2c, c3c = 2*tr/(1-r2), 2*sigma*ti/(1-r2)
    W = phi2-sp.I*sigma*phi3
    Ac = (phi2+H*c2c)*tc-sp.I*sigma*(phi3+H*c3c)*tc
    wc = (phi2+H*c2c)*e2c+(phi3+H*c3c)*e3c
    R = sp.Matrix([[sp.re(e2c),sp.re(e3c)],[sp.im(e2c),sp.im(e3c)]])
    rep(f"C9a sigma={sigma}: closed exponent, w, and slice determinant",
        sp.simplify(Ac-tc*W-H*Fh) == 0
        and sp.simplify(wc-sigma*(sp.conjugate(W)-W*tc**2+2*H*tc)) == 0
        and sp.expand(R.det()-sigma*(1-r2**2)) == 0)
t_, W_, Wbar_, A_, theta_, thetabar_ = sp.symbols('t W Wbar A theta thetabar')
f_, dH_, df_, dt_, dW_, dWbar_, dxs_, dys_ = sp.symbols('f dH df dt dW dWbar dxs dys')
dA_ = W_*dt_+t_*dW_+H*df_+f_*dH_
eta1, eta2 = dxs_-xs*dA_, dys_+ys*dA_
rep("C9b coframe: exp(A) dx and exp(-A) dy; dt=y_s eta1+x_s eta2",
    sp.simplify(sp.exp(A_)*(sp.exp(-A_)*dxs_-xs*sp.exp(-A_)*dA_)-eta1) == 0
    and sp.simplify(sp.exp(-A_)*(sp.exp(A_)*dys_+ys*sp.exp(A_)*dA_)-eta2) == 0
    and sp.expand(ys*eta1+xs*eta2-ys*dxs_-xs*dys_) == 0)
w_ = Wbar_-W_*t_**2+2*H*t_
dw_ = sum(sp.diff(w_,v)*dv for v,dv in ((t_,dt_),(W_,dW_),(Wbar_,dWbar_),(H,dH_)))
rep("C9c differential of w", sp.expand(dw_-(dWbar_-t_**2*dW_+(2*H-2*W_*t_)*dt_+2*t_*dH_)) == 0)
rep("C9d torus translation is holomorphic in (x,y,w), t=xy",
    sp.simplify(xs*sp.exp(-t_*(W_+theta_)-H*f_)-sp.exp(-theta_*t_)*xs*sp.exp(-t_*W_-H*f_)) == 0
    and sp.expand(w_.subs({W_:W_+theta_,Wbar_:Wbar_+thetabar_}, simultaneous=True)-w_-thetabar_+theta_*t_**2) == 0)
lam = sp.Symbol('lambda')
Axy = sp.Matrix([[1-lam*xs*ys,-lam*xs**2],[lam*ys**2,1+lam*xs*ys]])
shear = sp.eye(6)
for row in range(2):
    for col in range(4):
        shear[4+row,col] = sp.Symbol(f'cder{row}{col}')
rep("C9e determinant factors: complex normal block and real shear both have determinant 1",
    sp.expand(Axy.det()) == 1 and shear.det() == 1)
Xh,Yh,fh = sp.symbols('Xh Yh fh', positive=True)
jj = sp.Symbol('j', integer=True)
aa,bb = jj*(jj+1)/2,jj*(jj-1)/2
term = sp.exp(-jj*(jj+1)*fh)*Xh**aa*Yh**bb
rep("C9f termwise derivatives for the height coframe",
    sp.simplify(sp.diff(term,Xh)-aa*sp.exp(-jj*(jj+1)*fh)*Xh**(aa-1)*Yh**bb) == 0
    and sp.simplify(sp.diff(term,Yh)-bb*sp.exp(-jj*(jj+1)*fh)*Xh**aa*Yh**(bb-1)) == 0
    and sp.simplify(sp.diff(term,fh)+jj*(jj+1)*term) == 0)
dNh,dDh,Nh,Dh = sp.symbols('dNh dDh Nh Dh')
rep("C9g quotient derivative for H", sp.cancel((dNh*Dh-Nh*dDh)/Dh**2-(dNh-Nh/Dh*dDh)/Dh) == 0)
J6 = sp.diag(Jw,Jw,Jw)
rep("C9h real complex-structure block squares to -Id", J6*J6 == -sp.eye(6))
# The cap coordinate embedding has differential E4=(I4,-Z^T/X).
# Its left inverse on tangent vectors is R4=(I4-Z Z^T,-X Z).
zr = sp.Matrix(sp.symbols('z0:4', real=True))
height = sp.Symbol('height', nonzero=True, real=True)
E4 = sp.eye(4).col_join(-zr.T/height)
R4 = (sp.eye(4)-zr*zr.T).row_join(-height*zr)
rep("C9i tangent coordinate extraction composed with embedding is Id", sp.simplify(R4*E4) == sp.eye(4))
angle = sp.Symbol('angle', real=True)
Ec = 2*sp.pi*sp.Matrix([-sp.sin(angle),sp.cos(angle)])
Rc = Ec.T/(2*sp.pi)**2
rep("C9j circle-angle coordinate extraction composed with embedding is Id", sp.trigsimp((Rc*Ec)[0]) == 1)
normals = sp.Matrix.vstack(zr,sp.Matrix([height]))
rep("C9k coordinate extraction kills the sphere normal on |Z|^2+X^2=1",
    sp.simplify(R4*normals-zr*(1-(zr.T*zr)[0]-height**2)) == sp.zeros(4,1))
a0_ = sp.Symbol('a0',positive=True)
z1c,z2c = zr[0]+sp.I*zr[1],zr[2]+sp.I*zr[3]
for U in (U0,Ui):
    qlin = U.conjugate().T*sp.Matrix([z1c,z2c])/sp.sqrt(scale)
    xc,yc = sp.conjugate(qlin[0])/sp.sqrt(a0_),qlin[1]
    realmap = sp.Matrix([sp.re(xc),sp.im(xc),sp.re(yc),sp.im(yc)])
    rep("C9l real determinant of Z -> (xs,ys) is -1/(scale^2*a0)",
        sp.simplify(realmap.jacobian(zr).det()+1/(scale**2*a0_)) == 0)

# ---------------------------------------------------------------- C10: generic six-dimensional naturality
# For a local diffeomorphism F, normalize A=DF to Id at the point by a constant
# real linear change of target coordinates. Every symmetric Hessian remains
# possible. Then J=A J6 A^-1 and d_k J=H_k J6-J6 H_k.
# This is a jet identity, not a numerical test and not a check that F is a diffeomorphism.
dim = 6
Hjet = [[[sp.Symbol(f'H{out}_{min(i,k)}{max(i,k)}') for k in range(dim)]
         for i in range(dim)] for out in range(dim)]
dJjet = [sp.Matrix(dim,dim,lambda out,i: Hjet[out][i][k])*J6
          -J6*sp.Matrix(dim,dim,lambda out,i: Hjet[out][i][k]) for k in range(dim)]
rep("C10a differentiated J^2=-Id for an arbitrary second jet",
    all(sp.expand(x) == 0 for DJ in dJjet for x in DJ*J6+J6*DJ))
all_N = []
for out in range(dim):
    for i in range(dim):
        for j in range(dim):
            nij = sum(J6[p,i]*dJjet[p][out,j]-J6[p,j]*dJjet[p][out,i]
                      -J6[out,p]*(dJjet[i][p,j]-dJjet[j][p,i]) for p in range(dim))
            all_N.append(sp.expand(nij))
rep("C10b all 216 Nijenhuis components vanish for arbitrary symmetric Hessians",
    all(value == 0 for value in all_N))
# Control: A=Id+q_2 E_00 is a frame matrix, but is not a coordinate Jacobian.
# Its invertible conjugation still gives J^2=-Id, while N^1_{0,3}=-1 at zero.
bad_H = [sp.zeros(dim) for _ in range(dim)]
bad_H[2][0,0] = 1
bad_dJ = [HH*J6-J6*HH for HH in bad_H]
control_N = sum(J6[p,0]*bad_dJ[p][1,3]-J6[p,3]*bad_dJ[p][1,0]
                -J6[1,p]*(bad_dJ[0][p,3]-bad_dJ[3][p,0]) for p in range(dim))
rep("C10c non-coordinate-frame control has N^1_(0,3)=-1", control_N == -1)

# ---------------------------------------------------------------- C11: connection descent and transport identities
# The source circle acts by (x,y)->(exp(-i theta)x,exp(i theta)y).
xr,xi,yr,yi = sp.symbols('xr xi yr yi',real=True)
weight0 = sp.Symbol('weight0',positive=True)
radc = weight0*(xr*xr+xi*xi)+yr*yr+yi*yi
alpha_local = sp.Matrix([weight0*xi,-weight0*xr,-yi,yr])/radc
circle_generator = sp.Matrix([xi,-xr,-yi,yr])
rep("C11a local Hopf connection evaluates to 1 on the circle generator",
    sp.cancel((alpha_local.T*circle_generator)[0]-1) == 0)
Xc,Yc,dchi,dtheta,nidx = sp.symbols('Xc Yc dchi dtheta nidx',real=True)
rep("C11b torus expression alpha_n=dchi+(n-f_n)dtheta, f_n=X_n/(X_n+Y_n)",
    sp.cancel((-Xc*(dtheta-(dchi+nidx*dtheta))+Yc*(dchi+nidx*dtheta))/(Xc+Yc)
              -dchi-(nidx-Xc/(Xc+Yc))*dtheta) == 0)
av = sp.Symbol('average',real=True)
rep("C11c deck shift dchi->dchi+dtheta, average->average-1 preserves alpha",
    sp.expand(dchi+dtheta+(av-1)*dtheta-dchi-av*dtheta) == 0)
dx1_,dx4_,theta_base = sp.symbols('dx1 dx4 theta_base',real=True)
rep("C11d middle connection 2pi dx1+arg(t) dx4 is branch invariant",
    sp.expand(2*sp.pi*(dx1_-dx4_)+(theta_base+2*sp.pi)*dx4_-2*sp.pi*dx1_-theta_base*dx4_) == 0)
rep("C11e connection partition and exact linear-cap support",
    (sp.Rational(9,10)-sp.Rational(1,2))/(sp.Rational(9,10)-sp.Rational(2,5)) > sp.Rational(1,2)
    and sp.Rational(1,250)+sp.Rational(27,10)*sp.Rational(1,250)**2 < sp.Rational(1,10)
    and sp.Rational(1,20)+sp.Rational(27,10)*sp.Rational(1,20)**2 == sp.Rational(227,4000)
    and sp.Rational(227,4000) < sp.Rational(1,2))
# Gauge coordinates theta_b=theta_a+h(b) require A_b=A_a-dh.
# On a path, integral(A_b)=integral(A_a)-h(end)+h(start).
theta_start,integral_A,ha,hb = sp.symbols('theta_start integral_A h_start h_end',real=True)
rep("C11f quadrature transport is independent of the phase chart",
    sp.expand((theta_start+ha)-(integral_A-hb+ha)-(theta_start-integral_A+hb)) == 0)
Iab,Ibc = sp.symbols('Iab Ibc',real=True)
rep("C11g transport composes and reverses exactly",
    sp.expand((theta_start-Iab)-Ibc-(theta_start-(Iab+Ibc))) == 0
    and sp.expand(theta_start-Iab+Iab-theta_start) == 0)
Theta = sp.Symbol('Theta',real=True)
rep("C11h north-end mismatch is cancelled by the phase correction",
    sp.simplify(sp.exp(-sp.I*Theta)*sp.exp(sp.I*Theta)-1) == 0)
Bmap = sp.Matrix(3,3,sp.symbols('B0:9'))
triangular = sp.eye(6)
triangular[:3,:3] = Bmap
for row in range(3,6):
    for col in range(3):
        triangular[row,col] = sp.Symbol(f'phase_der{row}_{col}')
rep("C11i base map plus circle/torus translations has determinant det(base)",
    sp.expand(triangular.det()-Bmap.det()) == 0)

# C12: the torus trivialization has a global real correction on each gluing annulus.
ss = sp.Symbol('s')
rep("C12a Phi correction under log(t)->log(t)+2pi i is e4",
    sp.expand(sp.pi*(ss+1)*ss-sp.pi*ss*(ss-1)-2*sp.pi*ss) == 0
    and sp.expand(2*t_*(ss+1)-2*t_*ss-2*t_) == 0)
x2g,x3g,x4g,xi2g,xi3g,xi4g,chg2,chg3,Hg = sp.symbols('x2g x3g x4g xi2g xi3g xi4g chg2 chg3 Hg',real=True)
x4slice = x4g-xi4g
Cg = sp.Matrix([x2g-xi2g+x4slice*chg2,x3g-xi3g+x4slice*chg3])
Dg = sp.Matrix([xi2g+(Hg-x4slice)*chg2,xi3g+(Hg-x4slice)*chg3])
rep("C12b phi0=C-H*c=x23-D0 exactly",
    sp.simplify(Cg-Hg*sp.Matrix([chg2,chg3])-sp.Matrix([x2g,x3g])+Dg) == sp.zeros(2,1))
rep("C12c D0 invariant under deck x4,H -> x4+1,H+1",
    sp.simplify(Dg.subs({x4g:x4g+1,Hg:Hg+1},simultaneous=True)-Dg) == sp.zeros(2,1))
rep("C12d D0 invariant under branch xi4,H -> xi4+1,H-1",
    sp.simplify(Dg.subs({xi4g:xi4g+1,Hg:Hg-1},simultaneous=True)-Dg) == sp.zeros(2,1))
Di = sp.Matrix([(Hg+x4g)*chg2,(Hg+x4g)*chg3])
rep("C12e Dinf invariant under x4,H -> x4+1,H-1",
    sp.simplify(Di.subs({x4g:x4g+1,Hg:Hg-1},simultaneous=True)-Di) == sp.zeros(2,1))
u0,v0,ui0,uj0,vi0,vj0 = sp.symbols('u0 v0 ui0 uj0 vi0 vj0',real=True)
wedge = ui0*vj0-uj0*vi0
constraint_i,constraint_j = u0*ui0+v0*vi0,u0*uj0+v0*vj0
rep("C12f unit-circle phase form is closed: its exterior derivative is forced to vanish",
    sp.expand((u0*u0+v0*v0)*wedge-u0*(constraint_i*vj0-constraint_j*vi0)
              -v0*(ui0*constraint_j-uj0*constraint_i)) == 0)

# C13: explicit real logarithm of a smooth S^2 -> S^1 map, by two radial integrals.
u1s,u2s,tau_s = sp.symbols('u1s u2s tau_s',real=True)
usq = u1s*u1s+u2s*u2s
SN = sp.Matrix([2*u1s,2*u2s,1-usq])/(1+usq)
SS = sp.Matrix([2*u1s,2*u2s,usq-1])/(1+usq)
rep("C13a two stereographic parametrizations and their inverses",
    sp.cancel((SN.T*SN)[0]-1) == 0 and sp.cancel((SS.T*SS)[0]-1) == 0
    and sp.cancel(SN[0]/(1+SN[2])-u1s) == 0
    and sp.cancel(SN[1]/(1+SN[2])-u2s) == 0
    and sp.cancel(SS[0]/(1-SS[2])-u1s) == 0
    and sp.cancel(SS[1]/(1-SS[2])-u2s) == 0
    and SN.subs({u1s:1,u2s:0}) == SS.subs({u1s:1,u2s:0}))
# d_j int_0^1 A_i(tu) u_i dt = int_0^1 d_t(t A_j(tu)) dt
# when dA=0. The symmetric first jet is arbitrary at each tu.
a1s,a2s,h11,h12,h22 = sp.symbols('a1s a2s h11 h12 h22',real=True)
Avec = sp.Matrix([a1s,a2s])
Hclosed = sp.Matrix([[h11,h12],[h12,h22]])
uv = sp.Matrix([u1s,u2s])
rep("C13b radial primitive derivative identity for an arbitrary closed one-form",
    Avec+tau_s*Hclosed.T*uv == Avec+tau_s*Hclosed*uv)
gr,gi,dr,di = sp.symbols('gr gi dr di',real=True)
eta = gr*di-gi*dr
constraint = gr*dr+gi*di
rep("C13c d(exp(-i Theta)g)=0 when dTheta=Im(conj(g)dg) and |g|=1",
    sp.expand((dr+sp.I*di)-sp.I*(gr+sp.I*gi)*eta
              -(gr+sp.I*gi)*constraint-(1-gr*gr-gi*gi)*(dr+sp.I*di)) == 0)

# C14: inverse compositions and the two Hopf local reference sections.
nr1,nr2,nr3 = sp.symbols('nr1 nr2 nr3',real=True)
section_y = sp.Matrix([(nr1+sp.I*nr2)/sp.sqrt(2*(1-nr3)),sp.sqrt((1-nr3)/2)])
section_x = sp.Matrix([sp.sqrt((1+nr3)/2),(nr1-sp.I*nr2)/sp.sqrt(2*(1+nr3))])
# Algebra is performed with real squared amplitudes, avoiding unnecessary sqrt branches.
ny = (1-nr3)/2
nx = (1+nr3)/2
rep("C14a Hopf sections have unit norm on n1^2+n2^2+n3^2=1",
    sp.cancel((nr1**2+nr2**2)/(4*ny)+ny-1-(nr1**2+nr2**2+nr3**2-1)/(2*(1-nr3))) == 0
    and sp.cancel(nx+(nr1**2+nr2**2)/(4*nx)-1-(nr1**2+nr2**2+nr3**2-1)/(2*(1+nr3))) == 0)
chi_s,theta_s,bphase,flowphase = sp.symbols('chi_s theta_s bphase flowphase',real=True)
forwardphase = bphase-flowphase-chi_s*theta_s
rep("C14b forward and inverse transport phases compose to the identity",
    sp.expand(forwardphase+chi_s*theta_s+flowphase-bphase) == 0)
radius_s,epsilon_s = sp.symbols('radius_s epsilon_s',positive=True)
rep("C14c endpoint radial rescalings compose to the identity",
    sp.simplify(sp.sqrt(epsilon_s/radius_s)*sp.sqrt(radius_s/epsilon_s)-1) == 0)

# The projection to a closed Euclidean ball is nonexpanding. For p=ru,q=sv,
# u.v=z<=1, the squared-distance difference is a sum of nonnegative terms.
rproj,sproj,Rproj,zproj = sp.symbols('rproj sproj Rproj zproj',real=True)
rep("C14d radial projection: one point outside the ball",
    sp.expand(rproj*rproj+sproj*sproj-2*rproj*sproj*zproj
              -(Rproj*Rproj+sproj*sproj-2*Rproj*sproj*zproj)
              -(rproj-Rproj)*(rproj+Rproj-2*sproj*zproj)) == 0)
rep("C14e radial projection: both points outside the ball",
    sp.expand(rproj*rproj+sproj*sproj-2*rproj*sproj*zproj
              -2*Rproj*Rproj*(1-zproj)
              -(rproj-sproj)**2-2*(rproj*sproj-Rproj*Rproj)*(1-zproj)) == 0)

# C15: orbit-chart and same-Hopf-fibre identities used by the global construction.
rk,xx4,xxi4,tht = sp.symbols('rk xx4 xxi4 tht',real=True)
rep("C15a h_B normalization and the middle clutching phase",
    sp.cancel((1+rk*rk)/(1+rk*rk)-1) == 0
    and sp.expand(-(-xx4-sp.Rational(1,2))+tht/(2*sp.pi)-xxi4
                  -(xx4-xxi4-sp.Rational(1,2))-tht/(2*sp.pi)) == 1)
za,zb,zc,zd = sp.symbols('za zb zc zd',real=True)
Zc = sp.Matrix([za+sp.I*zb,zc+sp.I*zd])
radz = za*za+zb*zb+zc*zc+zd*zd
hz = 2*Zc[0]*sp.conjugate(Zc[1])
vz = za*za+zb*zb-zc*zc-zd*zd
projector_z = sp.Matrix([[radz+vz,hz],[sp.conjugate(hz),radz-vz]])/2
rep("C15b Hopf data determine the complex rank-one projector",
    sp.simplify(Zc*Zc.conjugate().T-projector_z) == sp.zeros(2))
rep("C15c Hopf norm equals the squared complex norm",
    sp.expand(sp.re(hz)**2+sp.im(hz)**2+vz**2-radz**2) == 0)

# The slice-radius equalities used by both the base clutching and the T^2 correction.
logabs,phasearg,ff,lnu,xone,xfour,xi_one,xi_four = sp.symbols(
    'logabs phasearg ff lnu xone xfour xi_one xi_four',real=True)
logcomplex = logabs+sp.I*phasearg
zslice0 = 2*sp.pi*(xone-xi_one)-sp.I*(xfour-xi_four)*(logcomplex-ff)
zslicei = 2*sp.pi*xone+sp.I*xfour*(logcomplex+lnu-ff)
rep("C15d cusp slice radii: rho0=(x4-xi4)P0 and rhoi=-x4 Pi",
    sp.expand(-sp.im(zslice0)-(xfour-xi_four)*(logabs-ff)) == 0
    and sp.expand(-sp.im(zslicei)+xfour*(logabs+lnu-ff)) == 0)
tc_sym,lc_sym = sp.symbols('tc_sym lc_sym',nonzero=True)
ez_mid = [-sp.I*tc_sym,-tc_sym,-sp.I*lc_sym]
ew_mid = [1-tc_sym**2,sp.I*(1+tc_sym**2),2*tc_sym]
ezi = [-sp.I/tc_sym,1/tc_sym,sp.I*(lc_sym+2)]
ewi = [1/tc_sym**2-1,sp.I*(1+1/tc_sym**2),-2/tc_sym]
rep("C15e Theta sends the middle e2,e3 to cusp e2,e3 and e4 to minus cusp e4",
    all(sp.simplify(ez_mid[j]-sp.I*ew_mid[j]/tc_sym-(1 if j<2 else -1)*ezi[j]) == 0
        and sp.simplify(ew_mid[j]/tc_sym**2-(1 if j<2 else -1)*ewi[j]) == 0 for j in range(3)))
rep("C15f explicit height root bracket and central endpoint",
    sp.Rational(1,2)+sp.Rational(27,10)*sp.Rational(1,20)**2 == sp.Rational(2027,4000)
    and sp.Rational(2027,4000) < 1 and sp.Rational(1,1)/(1+sp.Rational(1,1)) == sp.Rational(1,2))

# C16: the height is a standard theta logarithmic derivative, not a new
# finite elementary function. DLMF 20.2.3 uses exp(2 i n z), q^(n^2).
rho16,P16 = sp.symbols('rho16 P16',real=True)
n16 = sp.Symbol('n16',integer=True)
z16 = sp.I*(rho16-P16/2)
exponent16 = P16*n16**2+2*sp.I*n16*z16
rep("C16a theta summand exponent equals the original Gaussian exponent",
    sp.expand(exponent16-(-2*n16*rho16+n16*(n16+1)*P16)) == 0)
rep("C16b theta z-derivative divided by 2i gives the height numerator",
    sp.simplify(sp.diff(sp.exp(P16*n16**2+2*sp.I*n16*ss),ss)/(2*sp.I)
                -n16*sp.exp(P16*n16**2+2*sp.I*n16*ss)) == 0)
rep("C16c deck reindexing: w_n(rho+P)=exp(-2rho) w_(n-1)(rho)",
    sp.expand(exponent16.subs(rho16,rho16+P16)
              -exponent16.subs(n16,n16-1)+2*rho16) == 0)
import mpmath as mp
with mp.workdps(60):
    theta_error = mp.mpf(0)
    for pp in ('-.1','-.5','-1','-3','-12'):
        pp = mp.mpf(pp)
        for centre in ('-2.3','-.49','0','.37','2.4'):
            rr = (mp.mpf(centre)+mp.mpf('.5'))*pp
            z = 1j*(rr-pp/2)
            height_theta = mp.jtheta(3,z,mp.exp(pp),1)/(2j*mp.jtheta(3,z,mp.exp(pp)))
            weights = [(n,mp.exp(-2*n*rr+n*(n+1)*pp)) for n in range(-160,161)]
            height_sum = mp.fsum(n*w for n,w in weights)/mp.fsum(w for n,w in weights)
            theta_error = max(theta_error,abs(height_theta-height_sum))
    rep(f"C16d numerical 60-digit theta/series comparison: 25 points, error {float(theta_error):.1e}",theta_error < mp.mpf('1e-50'))

# C17: remove the height from the cap torus trivialization exactly.
phie,H17,n17,eta17,c17 = sp.symbols('phie H17 n17 eta17 c17',real=True)
correction17 = eta17*(H17-n17)*c17
rep("C17a torus correction is deck invariant on a uniquely reduced node slab",
    sp.expand(correction17.subs({H17:H17+1,n17:n17+1},simultaneous=True)-correction17) == 0)
Cold = phie-correction17+H17*c17
rep("C17b in the cap eta=1,n=0 the new coordinate is C=phi_E",
    sp.expand(Cold.subs({eta17:1,n17:0})-phie) == 0)
sx,sy,st,sW,sWb = sp.symbols('sx sy st sW sWb')
dx,dy,dt,dW,dWb = sp.symbols('dx dy dt dW dWb')
sA = st*sW
rep("C17c elementary cap coframe is the scaled differential of its explicit G",
    sp.simplify(sp.exp(sA)*(sp.diff(sx*sp.exp(-sA),sx)*dx
        +sp.diff(sx*sp.exp(-sA),st)*dt+sp.diff(sx*sp.exp(-sA),sW)*dW)
        -(dx-sx*(sW*dt+st*dW))) == 0
    and sp.expand(sum(sp.diff(sWb-sW*st**2,v)*dv for v,dv in ((sWb,dWb),(sW,dW),(st,dt)))
                  -(dWb-st**2*dW-2*sW*st*dt)) == 0)
for sigma17 in (1,-1):
    t17 = tr+sp.I*ti
    e217,e317 = sigma17*(1-t17**2),sp.I*(1+t17**2)
    M17 = sp.eye(6)
    for ii17 in range(4):
        M17[ii17,4],M17[ii17,5] = sp.symbols(f'b{ii17}_0 b{ii17}_1')
    M17[4:,4:] = sp.Matrix([[sp.re(e217),sp.re(e317)],[sp.im(e217),sp.im(e317)]])
    rep(f"C17d sigma={sigma17}: slice coframe determinant = sigma(1-|t|^4)",
        sp.expand(M17.det()-sigma17*(1-(tr*tr+ti*ti)**2)) == 0)
rep("C17e elementary-cap denominator and torus correction support margins",
    sp.Rational(1,100)**2*55/(4*(sp.Rational(157,50))**2) < sp.Rational(1,1000)
    and sp.Rational(1,125)+sp.Rational(27,10)*sp.Rational(1,125)**2 < sp.Rational(1,100)
    and sp.Rational(1,100) < sp.Rational(157,50)/250)

# ---------------------------------------------------------------- C18: fully specified algebraic Matsumoto--Fukaya target
a18,b18,c18,d18,h18 = sp.symbols('a18 b18 c18 d18 h18',real=True)
root18 = sp.sqrt(2)
v18 = sp.Matrix([a18,b18,c18,d18])
k18 = sp.Matrix([(root18*a18+1)/(root18+a18),b18/(root18+a18),
                 c18/(root18+a18),d18/(root18+a18)])
ki18 = sp.Matrix([(root18*k18[0]-1)/(root18-k18[0]),
                  k18[1]/(root18-k18[0]),k18[2]/(root18-k18[0]),k18[3]/(root18-k18[0])])
rep("C18a algebraic S3 boost and its closed inverse",
    all(sp.cancel(v)==0 for v in ki18-v18)
    and sp.cancel((k18.dot(k18)-1)*(root18+a18)**2-(v18.dot(v18)-1))==0)
zz18 = sp.Matrix([2*(a18*c18+b18*d18),2*(b18*c18-a18*d18),
                  a18*a18+b18*b18-c18*c18-d18*d18])
rad18 = a18*a18+b18*b18+c18*c18+d18*d18
rep("C18b Hopf norm identity and suspension sphere identity",
    sp.expand(zz18.dot(zz18)-rad18**2)==0
    and sp.expand((1-h18*h18)**2+h18*h18*(2-h18*h18)-1)==0)
rep("C18c suspension meridian singular value squared = 4(1-X^2)/(2-X^2)",
    sp.cancel(4*h18*h18*(1-h18*h18)+4*(1-h18*h18)**3/(2-h18*h18)
              -4*(1-h18*h18)/(2-h18*h18))==0)
# For a tangent vector w to S3, |Dk w|^2=|w|^2/(sqrt(2)+a)^2.
w18 = sp.Matrix(sp.symbols('w18_0:4',real=True))
dk18 = k18.jacobian(v18)*w18
res18 = sp.factor((dk18.dot(dk18)-w18.dot(w18)/(root18+a18)**2)*(root18+a18)**4)
rep("C18d conformal factor follows from the sphere and tangency equations",
    sp.expand(res18-(w18[0]**2*(v18.dot(v18)-1)
                   -2*w18[0]*(root18+a18)*(v18.dot(w18))))==0)
def hopf18(v):
    a,b,c,d = v
    return sp.Matrix([2*(a*c+b*d),2*(b*c-a*d),a*a+b*b-c*c-d*d])
rep("C18e the two poles have distinct values (0,-epsilon,0)",
    all(hopf18(k18.subs({a18:0,b18:0,c18:0,d18:sign}))==sp.Matrix([0,-sign,0])
        for sign in (-1,1)))
# At a pole the two tangent target coordinates have quadratic forms
# c+epsilon*b and sqrt(2)*a, where (a,b,c)=Hopf(z).
for sign in (-1,1):
    Q18 = sp.Matrix([zz18[2]+sign*zz18[1],root18*zz18[0]])
    l18,m18 = sp.symbols('l18 m18',real=True)
    Hess18 = sp.hessian(l18*Q18[0]+m18*Q18[1],(a18,b18,c18,d18))
    rep(f"C18f pole {sign}: every nonzero real quadratic pencil is nondegenerate",
        sp.factor(Hess18.det()-64*(l18*l18+m18*m18)**2)==0)
    eps18 = sp.symbols('eps18',real=True)
    pole18 = hopf18(k18.subs({a18:eps18**2*zz18[0],b18:eps18**2*zz18[1],
                             c18:eps18**2*zz18[2],d18:sign},simultaneous=True))
    rep(f"C18g pole {sign}: actual tangent-coordinate quadratic part",
        all(sp.simplify(sp.diff(pole18[j],eps18,2).subs(eps18,0)/2-Q18[i])==0
            for i,j in enumerate((0,2))))

# ---------------------------------------------------------------- C19: invertible elliptic height, with explicit branch bookkeeping
rho19,L19,omega19 = sp.symbols('rho19 L19 omega19',real=True)
Q19 = sp.symbols('Q19',positive=True)
rep("C19a arctangent derivative is the sech kernel",
    sp.simplify(sp.diff(sp.atan(sp.exp(2*rho19)),rho19)
                -sp.sech(2*rho19)).rewrite(sp.exp).simplify()==0)
rep("C19b beta/reflection Fourier transform denominator and Fourier coefficient",
    sp.simplify(sp.sin(sp.pi*(sp.Rational(1,2)-sp.I*omega19/4))
                -sp.cosh(sp.pi*omega19/4))==0
    and sp.simplify(2*Q19/(1+Q19**2)-sp.sech(-sp.log(Q19))).rewrite(sp.exp).simplify()==0)
u19 = sp.symbols('u19',real=True)
T19 = sp.atan(sp.pi*sp.tan(2*u19)/2)/sp.pi
rep("C19c normalized amplitude derivative and cubic Taylor term",
    sp.trigsimp(sp.diff(T19,u19)
                 -1/(sp.cos(2*u19)**2+sp.pi**2*sp.sin(2*u19)**2/4))==0
    and sp.simplify(sp.series(T19,u19,0,4).removeO()
                   -u19-(4-sp.pi**2)*u19**3/3)==0)
R19 = sp.Rational(1,20)
qn19 = R19**2/4
extra19 = 1/(3*(1-qn19**3))+2/(1-qn19)**3
vd19 = sp.Rational(265,309)/(1+R19**2)*(1-6*(R19/(1-qn19))**2)
rep("C19d normalized elliptic node: remainder constant 38, corrected slope >13/25",
    35+extra19<38 and vd19-2*38*R19**3/sp.Rational(3,100)>sp.Rational(13,25))

# Special-function identities are proved in the document by comparing Fourier
# coefficients with DLMF 22.11.3, using the beta integral. These high-precision
# evaluations are independent numerical evidence, not a substitute for that proof.
import mpmath as mp
import os,sys
sys.path.insert(0,os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..','..')))
from extra.s4t2_explicit.cusp import (elliptic_height,elliptic_height_inverse,elliptic_height_chart,
    elliptic_height_reparam,elliptic_axis_radius,elliptic_cusp_angle,elliptic_cusp_logradius)
with mp.workdps(80):
    worst19 = mp.mpf(0)
    cases19 = 0
    for ps in ('-1.1','-2','-5','-20','-60','-100'):
        P19 = mp.mpf(ps)
        q19 = mp.exp(2*P19)
        for ss in ('-2','-1.50001','-1','-.00000001','0','.1','.5','.9','1','2.7'):
            rho = P19*mp.mpf(ss)
            # Independent arctangent periodization, reduced to bounded arguments.
            n19 = mp.floor(rho/P19)
            Y19 = mp.exp(2*(rho-n19*P19))
            A19 = q19/Y19
            raw19 = n19+2/mp.pi*mp.fsum(mp.atan(A19*q19**j)-mp.atan(Y19*q19**j) for j in range(160))
            for norm19 in (False,True):
                h19 = elliptic_height(rho,P19,norm19,dps=80)
                expected19 = elliptic_height_reparam(raw19) if norm19 else raw19
                back19 = elliptic_height_inverse(h19,P19,norm19,dps=80)
                err19 = max(abs(h19-expected19),abs(back19-rho)/(1+abs(rho)))
                worst19 = max(worst19,err19)
                cases19 += 1
    rep(f"C19e {cases19} amplitude/independent-sum/inverse tests, 80 digits, error {float(worst19):.1e}",
        worst19<mp.mpf('1e-65'))
    axis19 = mp.mpf(0)
    for nu19 in (mp.mpf(1),mp.exp(-2)):
        for rr19 in (mp.mpf('.001'),mp.mpf('.05'),mp.mpf('.5'),mp.mpf(1),mp.mpf(3)):
            for axis in (0,1):
                x,y = (mp.sqrt(rr19)/nu19,mp.mpf(0)) if axis==0 else (mp.mpf(0),mp.sqrt(rr19))
                target19 = elliptic_height_reparam((1 if axis==0 else -1)*2/mp.pi*mp.atan(rr19))
                axis19 = max(axis19,abs(elliptic_height_chart(x,y,nu19,dps=80)-target19))
                axis19 = max(axis19,abs(elliptic_axis_radius(target19)-rr19))
    rep("C19f both exact axes (20 cases), plus the exact node",
        axis19<mp.mpf('1e-70') and elliptic_height_chart(0,0,dps=80)==0)
    collar19 = mp.mpf(0)
    for nu19 in (mp.mpf(1),mp.exp(-2)):
        for rr19 in map(mp.mpf,('.29','.3','.32','.4','.48','.5')):
            tc19 = rr19*mp.exp(mp.j*mp.mpf('.7'))
            for a19 in map(mp.mpf,('-2.5','-.4','0','.4','.5','1.3')):
                rh19 = elliptic_cusp_logradius(a19,tc19,nu19,dps=80)
                back19 = elliptic_cusp_angle(rh19,tc19,nu19,dps=80)
                collar19 = max(collar19,abs(back19-a19))
                if rr19 == mp.mpf('.5'):
                    PP = mp.log(nu19)+mp.log(rr19)-2*rr19**2/(1-rr19**2)
                    collar19 = max(collar19,abs(rh19-PP*(a19+mp.mpf('.5'))))
    rep(f"C19g 72 complete collar inverses and exact middle match, error {float(collar19):.1e}",
        collar19<mp.mpf('1e-65'))
    jets19 = mp.mpf(0)
    count_jets19 = 0
    for PP in (mp.mpf(-2),mp.mpf(-20)):
        qq = mp.exp(2*PP)
        def independent_height19(rho,normalized):
            AA,YY = qq*mp.exp(-2*rho),mp.exp(2*rho)
            h = 2/mp.pi*mp.fsum(mp.atan(AA*qq**j)-mp.atan(YY*qq**j)
                                for j in range(160))
            return elliptic_height_reparam(h) if normalized else h
        for cell19 in (-2,0,2):
            for normalized19 in (False,True):
                for order19 in (1,2,3):
                    lhs19 = mp.diff(lambda rr: elliptic_height(rr,PP,normalized19,dps=80),
                                    cell19*PP,order19)
                    rhs19 = mp.diff(lambda rr: independent_height19(rr,normalized19),
                                    cell19*PP,order19)
                    jets19 = max(jets19,abs(lhs19-rhs19))
                    count_jets19 += 1
    rep(f"C19h {count_jets19} seam jets through order three versus independent sum, error {float(jets19):.1e}",
        jets19<mp.mpf('1e-65'))

# ---------------------------------------------------------------- C20: finite phase logarithm and finite projection transport
a20,b20 = sp.symbols('a20 b20',real=True)
stereo20 = sp.Matrix([2*a20,2*b20,1-a20*a20-b20*b20])/(1+a20*a20+b20*b20)
ds20 = stereo20.jacobian([a20,b20])
rep("C20a stereographic metric: ||DS||=2/(1+|a|^2)<=2",
    (ds20.T*ds20-4*sp.eye(2)/(1+a20*a20+b20*b20)**2).applyfunc(sp.factor)==sp.zeros(2))
u20,v20,x20,y20 = sp.symbols('u20 v20 x20 y20',real=True)
rep("C20b unit phases: Re(z*conj(w))=1-|z-w|^2/2",
    sp.expand(2*(u20*x20+v20*y20)+(u20-x20)**2+(v20-y20)**2
              -(u20*u20+v20*v20+x20*x20+y20*y20))==0)
z20 = sp.symbols('z20_0:5',nonzero=True)
rep("C20c finite phase ratios telescope exactly",
    sp.cancel(sp.prod(z20[j]/z20[j-1] for j in range(1,5))-z20[4]/z20[0])==0)
rep("C20d reverse rank-one normalized projection cancels its phase",
    sp.cancel((x20+sp.I*y20)*(x20-sp.I*y20)/(x20*x20+y20*y20)-1)==0)

# ---------------------------------------------------------------- C21: the concrete ten-component embedding
h21 = sp.symbols('h21',real=True)
reduced21 = sp.symbols('reduced21',real=True)
j21 = sp.symbols('j21',integer=True)
m21,z21 = sp.symbols('m21 z21',nonzero=True)
rep('C21a height representative reindexing under the deck action',
    sp.simplify((sp.floor(sp.Rational(1,2)-(h21+1))-
                (sp.floor(sp.Rational(1,2)-h21)-1)).subs(h21,sp.Rational(1,2)-reduced21))==0
    and sp.powsimp(m21**(j21-1)*m21*z21-m21**j21*z21)==0)
b21 = sp.symbols('b21',real=True)
rep('C21b norm lower bound is a positive square',
    sp.expand(b21*b21+(1-b21)**2-sp.Rational(1,2)-2*(b21-sp.Rational(1,2))**2)==0)
th21,x21,y21 = sp.symbols('th21 x21 y21',real=True)
phase21 = 2*sp.pi*x21-2*sp.pi*j21*y21
rep('C21c middle lattice/branch phases are integer multiples of 2pi',
    sp.simplify(sp.exp(sp.I*(phase21.subs(y21,y21+1)-phase21)))==1
    and sp.expand((2*sp.pi*(x21-y21)-(2*sp.pi*j21-2*sp.pi)*y21)-phase21)==0)

# C22: rational inverse stereographic chart at the north pole, and the
# normalized Hopf projection identities used by the inverse-first cap.
a22,b22,c22 = sp.symbols('a22 b22 c22',real=True)
r22 = a22*a22+b22*b22+c22*c22
D22 = 1+2*b22+2*r22
ui22 = sp.Matrix([(-1-2*b22)/D22,2*a22/D22,2*(b22+r22)/D22,2*c22/D22])
rep('C22a stable north inverse lies in S3 and recovers the inverted stereographic vector',
    sp.factor(ui22.dot(ui22)-1)==0
    and all(sp.factor(v)==0 for v in (ui22[1]/(1+ui22[0])-a22/r22,
                    ui22[2]/(1+ui22[0])-1-b22/r22,ui22[3]/(1+ui22[0])-c22/r22)))
s22 = sp.symbols('s22',positive=True)
P22 = sp.Matrix([[s22+c22,a22+sp.I*b22],[a22-sp.I*b22,s22-c22]])/(2*s22)
rep('C22b rank-one projector identity, modulo |p|^2=s^2',
    (P22*P22-P22-sp.eye(2)*(r22-s22*s22)/(4*s22*s22)).applyfunc(sp.factor)==sp.zeros(2))
za22,zb22,zc22,zd22 = sp.symbols('za22 zb22 zc22 zd22',real=True)
zz22 = sp.Matrix([za22+sp.I*zb22,zc22+sp.I*zd22])
norm22 = za22*za22+zb22*zb22+zc22*zc22+zd22*zd22
hp22 = sp.Matrix([2*(za22*zc22+zb22*zd22),2*(zb22*zc22-za22*zd22),
                  za22*za22+zb22*zb22-zc22*zc22-zd22*zd22])
rep('C22c projected squared length equals the Hopf dot-product expression',
    sp.expand((zz22.conjugate().T*P22*zz22)[0]-(norm22+hp22.dot(sp.Matrix([a22,b22,c22]))/s22)/2)==0)
Dh22 = hp22.jacobian([za22,zb22,zc22,zd22])
rep('C22d Hopf horizontal right inverse has norm 1/(2 sqrt(r))',
    (Dh22*Dh22.T-4*norm22*sp.eye(3)).applyfunc(sp.expand)==sp.zeros(3))
stereo22 = sp.Matrix([1-r22,2*a22,2*b22,2*c22])/(1+r22)
Ds22 = stereo22.jacobian([a22,b22,c22])
rep('C22e inverse stereographic S3 derivative norm is 2/(1+|w|^2)',
    (Ds22.T*Ds22-4*sp.eye(3)/(1+r22)**2).applyfunc(sp.factor)==sp.zeros(3))

# C23: elementary inverse-first height, exact seams, and independent pole chain.
from extra.s4t2_explicit.cusp import (inverse_height_moment,elementary_cusp_logradius,
    elementary_cusp_reference,elliptic_pole_raw_inverse,elliptic_pole_patch,scalar_flat_step)
E23,Y23,q23,f23 = sp.symbols('E23 Y23 q23 f23',positive=True)
rep('C23a quadratic inverse: E-q/E=f gives Y=q/E and the endpoint jets',
    sp.factor((Y23**2+f23*Y23-q23).subs({Y23:q23/E23,f23:E23-q23/E23}))==0)
rep('C23a2 open positive-axis inverse, modulo R^2=f^2+4q',
    sp.factor((E23+f23)/2-2*q23/(E23+f23)-f23-(E23**2-f23**2-4*q23)/(2*(E23+f23)))==0)
rep('C23a3 half-centred argument chart has .4<a<.6',sp.cos(sp.pi/5)<sp.Rational(85,100))
rep('C23b transition positivity derivative is a sum of squares',
    sp.factor(2*E23+1/(4*E23)-1-(2*(E23-sp.Rational(1,4))**2+sp.Rational(1,8))/E23)==0)
rep('C23c normalized-projection connection phase propagates with coefficient one',
    sp.cancel((q23+f23*Y23)/Y23-f23-q23/Y23)==0)
ct23,st23 = sp.symbols('ct23 st23',real=True)
aa23,bb23,er23,ei23,fr23,fi23 = sp.symbols('aa23 bb23 er23 ei23 fr23 fi23',real=True)
uu23,vv23 = sp.Matrix([1,0]),sp.Matrix([ct23,st23])
du23 = sp.Matrix([sp.I*aa23,er23+sp.I*ei23])
dv23 = sp.I*bb23*vv23+(fr23+sp.I*fi23)*sp.Matrix([-st23,ct23])
PP23,QQ23 = uu23*uu23.conjugate().T,vv23*vv23.conjugate().T
dPP23 = du23*uu23.conjugate().T+uu23*du23.conjugate().T
dQQ23 = dv23*vv23.conjugate().T+vv23*dv23.conjugate().T
dc23 = (dv23.conjugate().T*uu23+vv23.conjugate().T*du23)[0]
actual23 = bb23+sp.im(dc23)/ct23
claimed23 = aa23+sp.im((uu23.conjugate().T*QQ23*(dQQ23+dPP23)*uu23)[0])/ct23**2
rep('C23c2 full two-line normalized-projection phase derivative, all tangent jets',
    sp.factor(actual23-claimed23).subs(st23**2,1-ct23**2).simplify()==0)
with mp.workdps(85):
    height_error23 = mp.mpf(0)
    slope_min23,slope_max23 = mp.inf,mp.mpf(0)
    height_cases23 = 0
    for nu23 in (mp.mpf(1),mp.exp(-2)):
        for r23 in map(mp.mpf,('1e-12','.03','.3','.4','.499','.5')):
            tc23 = r23 if r23==mp.mpf('.5') else r23*mp.exp(mp.j*mp.mpf('.7'))
            P23 = mp.log(nu23*r23)-2*r23*r23/(1-r23*r23)
            qv23 = mp.exp(2*P23)
            lam23 = scalar_flat_step((r23-mp.mpf('.3'))/mp.mpf('.2'))
            for a23 in map(mp.mpf,('-.5','-.499','-.3','-.2','-.1','0','.1','.2','.3','.499','.5')):
                rho23 = elementary_cusp_logradius(a23,tc23,nu23,dps=65)
                height_error23 = max(height_error23,abs(elementary_cusp_logradius(a23+1,tc23,nu23,dps=65)-rho23-P23))
                derivative23 = mp.diff(lambda aa:elementary_cusp_logradius(aa,tc23,nu23,dps=65),a23)
                assert derivative23<0
                if abs(a23)<mp.mpf('.5'):
                    slope23 = mp.diff(lambda aa:inverse_height_moment(aa,qv23,lam23),a23)
                    slope_min23,slope_max23 = min(slope_min23,slope23),max(slope_max23,slope23)
                if r23==mp.mpf('.5'):
                    height_error23 = max(height_error23,abs(rho23-P23*(a23+mp.mpf('.5'))))
                height_cases23 += 1
            for seam23 in (mp.mpf('-.5'),mp.mpf('.5'),mp.mpf('1.5')):
                for order23 in (1,2,3):
                    derivative23 = mp.diff(lambda aa:elementary_cusp_logradius(aa,tc23,nu23,dps=65),seam23,order23)
                    expected23 = -(1-lam23)+lam23*P23 if order23==1 else 0
                    height_error23 = max(height_error23,abs(derivative23-expected23))
    rep(f'C23d {height_cases23} elementary height/deck/collar cases and 108 seam jets: {float(height_error23):.1e}',
        height_error23<mp.mpf('1e-50') and slope_min23>0)
    print(f'     sampled f_a range [{float(slope_min23):.6f},{float(slope_max23):.6f}] (numerical observation)')
    axis_reference_error23 = mp.mpf(0)
    for nu23 in (mp.mpf(1),mp.exp(-2)):
        for a23 in map(mp.mpf,('-.55','-.5','-.49','.49','.5','.55')):
            tc23 = mp.mpc('.1','.03')
            xx23,yy23 = elementary_cusp_reference(tc23,a23,nu23)
            axis_reference_error23 = max(axis_reference_error23,
                abs(xx23*yy23-tc23),abs(mp.log(abs(yy23))-elementary_cusp_logradius(a23,tc23,nu23,dps=85)))
        for a23 in (mp.mpf('-.5'),mp.mpf('.5')):
            for order23 in (1,2,3):
                coord23 = 0 if a23>0 else 1
                derivative23 = mp.diff(lambda aa:elementary_cusp_reference(0,aa,nu23)[coord23],a23,order23)
                expected23 = 1/nu23 if a23>0 else (-1)**order23
                axis_reference_error23 = max(axis_reference_error23,abs(derivative23-expected23))
    rep(f'C23d2 open central-axis references: 12 inverse comparisons and 12 seam jets, {float(axis_reference_error23):.1e}',
        axis_reference_error23<mp.mpf('1e-70'))
    pole_error23,value_ratio23,derivative_ratio23,patch_defect23 = [mp.mpf(0)]*4
    pole_cases23 = 0
    def forward_local23(pair,inf):
        nu = mp.exp(-2) if inf else mp.mpf(1)
        x,y = mp.conj(pair[0])/nu,pair[1]
        tc = x*y
        u,v = abs(pair[0])**2,abs(pair[1])**2
        e = mp.exp(-4*abs(tc)**2/(1-abs(tc)**2))
        height = e*u-v
        A,B = (4*mp.pi)**2,mp.exp(4)
        if inf:
            ratio = tc*mp.sqrt((abs(tc)**2+A)/(A*abs(tc)**2+B))
            u1 = mp.exp(2*mp.j*mp.pi*(-height+mp.mpf('.5')))/mp.sqrt(1+abs(ratio)**2)
            u2 = u1*ratio
        else:
            ratio = tc*mp.sqrt((A+B*abs(tc)**2)/(1+A*abs(tc)**2))
            u2 = mp.exp(2*mp.j*mp.pi*height)/mp.sqrt(1+abs(ratio)**2)
            u1 = ratio*u2
        w = (u1.imag/(1+u1.real),u2.real/(1+u1.real)-1,u2.imag/(1+u1.real))
        if inf:
            norm2 = sum(a*a for a in w)
            w = tuple(a/norm2 for a in w)
        return mp.matrix([w[0],-w[2],-w[1]])
    for inf23 in (False,True):
        # Stay strictly inside the closed certified domain: rounding a unit
        # vector exactly on its outer boundary can put its mpf norm outside it.
        for r23 in map(mp.mpf,('1e-12','1e-8','1.25e-8','1.5e-8','1.75e-8','2e-8','1e-6','9.99e-5')):
            for direction23 in ((1,0,0,0),(0,0,1,0),(1,2,3,-1),(-2,1,0,3)):
                zz23 = mp.matrix(direction23)*mp.sqrt(r23/sum(a*a for a in direction23))
                pair23 = (mp.mpc(zz23[0],zz23[1]),mp.mpc(zz23[2],zz23[3]))
                raw23 = elliptic_pole_raw_inverse(pair23,inf23,dps=65,elementary=True)
                pc23 = 2*pair23[0]*mp.conj(pair23[1])
                p23 = (pc23.real,pc23.imag,abs(pair23[0])**2-abs(pair23[1])**2)
                target23 = mp.matrix((p23[2],-p23[1],p23[0]) if inf23 else (-p23[1],-p23[2],p23[0]))*(mp.pi if inf23 else 2*mp.pi)
                pole_error23 = max(pole_error23,mp.norm(forward_local23(raw23,inf23)-target23)/r23)
                value_ratio23 = max(value_ratio23,mp.norm(mp.matrix(raw23)-mp.matrix(pair23))/r23**mp.mpf('1.5'))
                if r23<mp.mpf('1e-4'):
                    def real_lift23(v,patch=False):
                        fn = elliptic_pole_patch if patch else elliptic_pole_raw_inverse
                        a,b = fn((mp.mpc(v[0],v[1]),mp.mpc(v[2],v[3])),inf23,dps=65,elementary=True)
                        return mp.matrix([a.real,a.imag,b.real,b.imag])
                    jac23 = mp.matrix(4)
                    patchjac23 = mp.matrix(4)
                    for k23 in range(4):
                        basis23 = mp.eye(4)[:,k23]
                        jac23[:,k23] = mp.diff(lambda h:real_lift23(zz23+h*basis23),0)
                        patchjac23[:,k23] = mp.diff(lambda h:real_lift23(zz23+h*basis23,True),0)
                    derivative_ratio23 = max(derivative_ratio23,mp.norm(jac23-mp.eye(4))/r23)
                    if r23<=mp.mpf('2e-8'):
                        patch_defect23 = max(patch_defect23,mp.norm(patchjac23-mp.eye(4)))
                pole_cases23 += 1
    rep(f'C23e {pole_cases23} independent inverse-first pole compositions: {float(pole_error23):.1e}',pole_error23<mp.mpf('1e-45'))
    rep(f'C23f pole value/derivative constants {float(value_ratio23):.3f}/{float(derivative_ratio23):.3f} vs 200/5000; blend {float(patch_defect23):.1e} vs .001',
        value_ratio23<200 and derivative_ratio23<5000 and patch_defect23<mp.mpf('.001'))

# C24: finite first-jet arithmetic eliminates the formal derivatives in (11.19).
v24,a24 = sp.symbols('v24 a24',positive=True)
w24,b24,t24 = sp.symbols('w24 b24 t24',real=True)
rules24 = [(1/v24,-1/v24**2),(sp.sqrt(v24),1/(2*sp.sqrt(v24))),
           (sp.exp(v24),sp.exp(v24)),(sp.log(v24),1/v24),
           (sp.sin(v24),sp.cos(v24)),(sp.cos(v24),-sp.sin(v24)),
           (sp.atan(v24),1/(1+v24**2))]
rep('C24a elementary unary first-jet rules',all(sp.simplify(sp.diff(f,v24)-df)==0 for f,df in rules24))
rep('C24b addition and multiplication first-jet rules',
    sp.diff((v24+t24*w24)+(a24+t24*b24),t24).subs(t24,0)==w24+b24
    and sp.expand(sp.diff((v24+t24*w24)*(a24+t24*b24),t24).subs(t24,0)-a24*w24-v24*b24)==0)
rep('C24c local argument first jet',
    sp.factor(sp.diff(sp.atan((a24+t24*b24)/(v24+t24*w24)),t24).subs(t24,0)
              -(v24*b24-a24*w24)/(v24*v24+a24*a24))==0)
arg24 = 2*sp.atan(a24/(sp.sqrt(v24*v24+a24*a24)+v24))
rep('C24d half-angle argument chart has the same first jet',
    sp.simplify(sp.diff(arg24,v24)*w24+sp.diff(arg24,a24)*b24
                -(v24*b24-a24*w24)/(v24*v24+a24*a24))==0)

# C25: exact forward reversal on raw caps and the south raw-base obstruction.
za25,zb25,zc25,zd25,wa25,wb25,wc25,wd25 = sp.symbols(
    'za25 zb25 zc25 zd25 wa25 wb25 wc25 wd25',real=True)
z25 = sp.Matrix([za25+sp.I*zb25,zc25+sp.I*zd25])
w25 = sp.Matrix([wa25+sp.I*wb25,wc25+sp.I*wd25])
rz25 = (z25.conjugate().T*z25)[0]
sw25 = (w25.conjugate().T*w25)[0]
Pz25 = z25*z25.conjugate().T/rz25
Pw25 = w25*w25.conjugate().T/sw25
ov25 = (z25.conjugate().T*w25)[0]
rep('C25a reverse raw-cap projection has positive real phase',
    all(sp.simplify(v)==0 for v in
        Pw25*Pz25*w25-sp.conjugate(ov25)*ov25*w25/(rz25*sw25)))
rep('C25a2 explicit nonzero projection margin',
    sp.Rational(40,250)/(1-sp.Rational(20,250))==sp.Rational(4,23)
    and 1-sp.Rational(4,23)**2/4==sp.Rational(525,529))

a25,b25,c25 = sp.symbols('a25 b25 c25',real=True)
d25 = sp.sqrt(1+a25*a25+b25*b25)
D25 = d25+a25*sp.cos(c25)+b25*sp.sin(c25)
Bnum25 = [a25*sp.sin(c25)-b25*sp.cos(c25),-sp.sin(c25),D25-sp.cos(c25)]
rep('C25b exact norm identity for the south raw-base block B',
    sp.trigsimp(sp.expand(sum(v*v for v in Bnum25)-2*D25*(d25-sp.cos(c25))))==0)

tt25,p125,p225,p325,rr25,kk25,bb25 = sp.symbols(
    'tt25 p125 p225 p325 rr25 kk25 bb25',real=True)
rho225 = p125**2+p225**2
u25,v25 = (rr25+p325)/2,(rr25-p325)/2
R325 = -4*u25*u25*v25
ac25 = kk25*tt25*p125*(1+bb25*tt25**2*rho225)
bc25 = kk25*tt25*p225*(1+bb25*tt25**2*rho225)
cc25 = kk25*(tt25*p325+tt25**3*R325)
dc25 = sp.sqrt(1+ac25*ac25+bc25*bc25)
Dc25 = dc25+ac25*sp.cos(cc25)+bc25*sp.sin(cc25)
F425 = (2*bb25*rho225**2+2*p325*R325
        +kk25**2*(rr25**2*(p125**2-rho225/2-p225*p325)-rho225**2/4-p325**4/12))
rawnorm25 = 2*(dc25-sp.cos(cc25))/(Dc25*kk25**2)
target25 = tt25**2*rr25**2-kk25*tt25**3*p125*rr25**2+tt25**4*F425
difference25 = sp.expand(sp.series(rawnorm25,tt25,0,5).removeO()-target25)
unit25 = rr25**2-p125**2-p225**2-p325**2
rep('C25b2 south raw-base norm expansion through order four',
    all(sp.factor(sp.rem(sp.expand(difference25.coeff(tt25,j)),unit25,rr25))==0
        for j in range(2,5)))

ee25,uu25,vv25 = sp.symbols('ee25 uu25 vv25',real=True)
beta25 = sp.sqrt((1+sp.exp(4)*ee25**2/(64*sp.pi**2))/(1+4*sp.pi**2*ee25**2))
beta225 = sp.exp(4)/(128*sp.pi**2)-2*sp.pi**2
height25 = ee25*uu25*(sp.exp(-4*ee25**2*uu25*vv25/(1-ee25**2*uu25*vv25))-1)
rep('C25c algebraic profile and elementary height jets',
    sp.simplify(sp.diff(beta25,ee25,2).subs(ee25,0)/2-beta225)==0
    and sp.simplify(sp.series(height25,ee25,0,4).removeO()
                    +4*ee25**3*uu25**2*vv25)==0)

x25,y25 = sp.symbols('x25 y25',real=True)
N25 = F425-kk25**2*p125**2*rr25**2/4
Nxy25 = sp.expand(N25.subs({p125:2*x25*y25,p225:0,p325:x25**2-y25**2,
                            rr25:x25**2+y25**2,kk25:2*sp.pi,bb25:beta225}))
eta25 = 1+sp.exp(4)/(64*sp.pi**2)-16*sp.pi**2/3
rep('C25d nonpolynomial sixth-order norm term on a real spinor plane',
    sp.simplify(Nxy25.subs(x25,sp.I*y25)-16*eta25*y25**8)==0
    and sp.rem(Nxy25,x25**2+y25**2,x25)!=0)
rep('C25d2 explicit obstruction coefficient bound using e^4<55 and pi>3',
    1+sp.Rational(55,576)-48 < -46)

with mp.workdps(80):
    forward_error25 = mp.mpf(0)
    forward_cases25 = 0
    for inf25 in (False,True):
        scale25 = mp.pi if inf25 else 2*mp.pi
        for radius25 in map(mp.mpf,('3e-8','1e-6','7e-5')):
            for direction25 in ((1,0,0,0),(1,2,3,-1),(-2,1,0,3)):
                unit25 = mp.matrix(direction25)/mp.sqrt(sum(v*v for v in direction25))
                source25 = (mp.mpc(unit25[0],unit25[1])*mp.sqrt(radius25),
                            mp.mpc(unit25[2],unit25[3])*mp.sqrt(radius25))
                base25 = forward_local23(source25,inf25)
                pL25 = mp.matrix([base25[2],-base25[1],base25[0]])/scale25 if inf25 else \
                       mp.matrix([base25[2],-base25[0],-base25[1]])/scale25
                rL25 = mp.norm(pL25)
                assert 2*mp.mpf('1e-8')<rL25<mp.mpf('1e-4')
                nL25 = pL25/rL25
                projector25 = mp.matrix([[1+nL25[2],nL25[0]+mp.j*nL25[1]],
                                         [nL25[0]-mp.j*nL25[1],1-nL25[2]]])/2
                projected25 = projector25*mp.matrix(source25)
                zL25 = projected25*mp.sqrt(rL25)/mp.norm(projected25)
                recovered25 = elliptic_pole_raw_inverse((zL25[0],zL25[1]),inf25,
                                                        dps=70,elementary=True)
                forward_error25 = max(forward_error25,
                    mp.norm(mp.matrix(recovered25)-mp.matrix(source25))/mp.sqrt(radius25))
                forward_cases25 += 1
    rep(f'C25e {forward_cases25} independent raw-cap forward compositions: '
        f'{float(forward_error25):.1e}',forward_error25<mp.mpf('1e-60'))

# C26: the patched shell descends to a three-real-variable Hopf-base map.
alpha26,gamma26 = sp.symbols('alpha26 gamma26',real=True)
h26 = sp.expand(ov25*sp.conjugate(ov25))/(rz25*sw25)
sigma_z26,sigma_w26 = 2*Pz25-sp.eye(2),2*Pw25-sp.eye(2)
rep('C26a rank-one projector anticommutator',
    all(sp.cancel(v)==0 for v in
        Pz25*Pw25+Pw25*Pz25-h26*sp.eye(2)-(sigma_z26+sigma_w26)/2))

mix26 = alpha26*z25+gamma26*Pw25*z25
def hopf26(pair):
    cross = 2*pair[0]*sp.conjugate(pair[1])
    return sp.Matrix([sp.re(cross).expand(complex=True),
                      sp.im(cross).expand(complex=True),
                      sp.expand(pair[0]*sp.conjugate(pair[0])
                                -pair[1]*sp.conjugate(pair[1]))])

rep('C26b Hopf vector of a projected spinor blend',
    all(sp.cancel(v)==0 for v in
        hopf26(mix26)-(alpha26**2+alpha26*gamma26)*hopf26(z25)
        -rz25/sw25*(gamma26**2*h26+alpha26*gamma26)*hopf26(w25)))

r26,s26,hpos26,theta26 = sp.symbols('r26 s26 hpos26 theta26',positive=True)
a26 = 1-theta26
g26 = theta26*sp.sqrt(s26/(r26*hpos26))
rep('C26c phase projection and blended squared norm',
    all(sp.cancel(v)==0 for v in
        Pz25*mix26-(alpha26+gamma26*h26)*z25)
    and sp.cancel((mix26.conjugate().T*mix26)[0]
                  -rz25*(alpha26**2+2*alpha26*gamma26*h26
                          +gamma26**2*h26))==0
    and sp.simplify(a26**2*r26+2*a26*g26*r26*hpos26
                    +g26**2*r26*hpos26
                    -(a26**2*r26+theta26**2*s26
                      +2*a26*theta26*sp.sqrt(r26*s26*hpos26)))==0
    and sp.simplify((a26**2+a26*g26)*r26-a26**2*r26
                    -a26*theta26*sp.sqrt(r26*s26/hpos26))==0
    and sp.simplify(r26/s26*(g26**2*hpos26+a26*g26)
                    -theta26**2-a26*theta26*sp.sqrt(r26/(s26*hpos26)))==0)

with mp.workdps(80):
    def hopf_mp26(pair):
        cross = 2*pair[0]*mp.conj(pair[1])
        return mp.matrix([cross.real,cross.imag,abs(pair[0])**2-abs(pair[1])**2])

    shell_error26 = mp.mpf(0)
    shell_cases26 = 0
    epsilon26 = mp.mpf('1e-8')
    for inf26 in (False,True):
        scale26 = mp.pi if inf26 else 2*mp.pi
        for radius26 in map(mp.mpf,('1e-8','1.1e-8','1.5e-8','1.9e-8','2e-8')):
            for direction26 in ((1,0,0,0),(0,0,1,0),(-2,1,0,3)):
                unit26 = mp.matrix(direction26)/mp.sqrt(sum(v*v for v in direction26))
                z26 = mp.matrix([mp.mpc(unit26[0],unit26[1])*mp.sqrt(radius26),
                                 mp.mpc(unit26[2],unit26[3])*mp.sqrt(radius26)])
                raw26 = mp.matrix(elliptic_pole_raw_inverse(tuple(z26),inf26,
                                                              dps=70,elementary=True))
                output26 = mp.matrix(elliptic_pole_patch(tuple(z26),inf26,
                                                           dps=70,elementary=True))
                p26,q26 = hopf_mp26(z26),hopf_mp26(raw26)
                rL26,rS26 = mp.norm(p26),mp.norm(q26)
                nL26,nS26 = p26/rL26,q26/rS26
                hval26 = (1+sum(nL26[k]*nS26[k] for k in range(3)))/2
                theta_val26 = scalar_flat_step((rL26-epsilon26)/epsilon26)
                base26 = forward_local23(tuple(raw26),inf26)
                forward_p26 = (mp.matrix([base26[2],-base26[1],base26[0]])/scale26
                               if inf26 else
                               mp.matrix([base26[2],-base26[0],-base26[1]])/scale26)
                predicted26 = ((1-theta_val26)**2*p26+theta_val26**2*q26
                               +theta_val26*(1-theta_val26)
                               *mp.sqrt(rL26*rS26/hval26)*(nL26+nS26))
                predicted_norm26 = ((1-theta_val26)**2*rL26+theta_val26**2*rS26
                                    +2*theta_val26*(1-theta_val26)
                                    *mp.sqrt(rL26*rS26*hval26))
                PL26 = mp.matrix([[1+nL26[2],nL26[0]+mp.j*nL26[1]],
                                  [nL26[0]-mp.j*nL26[1],1-nL26[2]]])/2
                projected26 = PL26*output26
                recovered26 = mp.sqrt(rL26)*projected26/mp.norm(projected26)
                shell_error26 = max(shell_error26,
                    mp.norm(forward_p26-p26)/rL26,
                    mp.norm(hopf_mp26(output26)-predicted26)/rL26,
                    abs(mp.norm(hopf_mp26(output26))-predicted_norm26)/rL26,
                    mp.norm(projected26-((1-theta_val26)
                        +theta_val26*mp.sqrt(rS26*hval26/rL26))*z26)/mp.sqrt(rL26),
                    mp.norm(recovered26-z26)/mp.sqrt(rL26))
                shell_cases26 += 1
    rep(f'C26d {shell_cases26} two-cap shell base/norm/phase compositions: '
        f'{float(shell_error26):.1e}',shell_error26<mp.mpf('1e-60'))

# C27: the south raw lift has a non-polynomial homogeneous cubic term.
p127,p227,p327,tau27 = sp.symbols('p127 p227 p327 tau27',real=True)
R27 = sp.symbols('R27',positive=True)
kappa27 = 2*sp.pi
p27 = sp.Matrix([p127,p227,p327])
delta27 = kappa27*sp.Matrix([
    (-p127**2+p227**2+p327**2)/2+p227*p327,
    -p127*(p227+p327),-p127*p327])
a27,b27,c27 = (tau27*kappa27*v for v in p27)
d27 = sp.sqrt(1+a27*a27+b27*b27)
D27 = d27+a27*sp.cos(c27)+b27*sp.sin(c27)
B27 = sp.Matrix([a27*sp.sin(c27)-b27*sp.cos(c27),
                 -sp.sin(c27),D27-sp.cos(c27)])/D27
Q027 = sp.Matrix([[0,-1,0],[0,0,-1],[1,0,0]])
C027 = Q027.T*B27/kappa27
rho227 = p127**2+p227**2
beta27 = sp.sqrt((1+sp.exp(4)*tau27**2*rho227/(64*sp.pi**2))
                 /(1+4*sp.pi**2*tau27**2*rho227))
RH27 = (sp.exp(-4*tau27**2*rho227/(4-tau27**2*rho227))-1)\
       *tau27*(R27+p327)/2
rep('C27a south B/C0 quadratic jet and cubic onset of profile/height terms',
    all(sp.expand(sp.series(C027[k],tau27,0,3).removeO()
                  -tau27*p27[k]-tau27**2*delta27[k])==0 for k in range(3))
    and sp.simplify(sp.series(tau27*(beta27-1),tau27,0,3).removeO())==0
    and sp.simplify(sp.series(RH27,tau27,0,3).removeO())==0)

def sigma27(v):
    return sp.Matrix([[v[2],v[0]+sp.I*v[1]],
                      [v[0]-sp.I*v[1],-v[2]]])

dot27 = p27.dot(delta27)/R27
s_tau27 = sp.sqrt(R27**2-2*tau27*p27.dot(delta27)
                  +tau27**2*delta27.dot(delta27))
Q_tau27 = (sp.eye(2)+sigma27(p27-tau27*delta27)/s_tau27)/2
Qdot27 = -sigma27(delta27)/(2*R27)+dot27*sigma27(p27)/(2*R27**2)
spinor_p27 = hopf26(z25)
spinor_delta27 = delta27.subs(dict(zip(p27,spinor_p27)),simultaneous=True)
spinor_r27 = rz25
spinor_dot27 = spinor_p27.dot(spinor_delta27)/spinor_r27
spinor_Qdot27 = (-sigma27(spinor_delta27)
                 +spinor_dot27*sigma27(spinor_p27)/spinor_r27)/(2*spinor_r27)
norm_dot27 = (z25.conjugate().T*spinor_Qdot27*z25)[0]
formal_dot27 = (-spinor_dot27*z25/(2*spinor_r27)
                +spinor_Qdot27*z25-norm_dot27*z25/(2*spinor_r27))
rep('C27b normalized-projector first variation cancels to the stated cubic',
    all(sp.simplify(v)==0 for v in
        sp.diff(Q_tau27,tau27).subs(tau27,0)-Qdot27)
    and all(sp.cancel(v)==0 for v in
            sigma27(spinor_p27)*z25-spinor_r27*z25)
    and sp.cancel(norm_dot27)==0
    and all(sp.cancel(v)==0 for v in
            formal_dot27+sigma27(spinor_delta27)*z25/(2*spinor_r27)))

x27,y27 = sp.symbols('x27 y27',real=True)
plane_delta27 = delta27.subs({p127:2*x27*y27,p227:0,
                              p327:x27**2-y27**2})
imag_num27 = sp.factor(-sp.im((sigma27(plane_delta27)
                                  *sp.Matrix([x27,y27]))[0])/2)
rep('C27c real-plane imaginary cubic numerator is not divisible by r',
    sp.simplify(imag_num27-kappa27*x27*y27**2*(x27**2-y27**2))==0
    and sp.simplify(imag_num27.subs(x27,sp.I*y27)
                    +2*sp.I*kappa27*y27**5)==0
    and sp.rem(imag_num27,x27**2+y27**2,x27)!=0)

# C28: the new pole patch is a flat radial conjugation with an elementary inverse.
from extra.s4t2_explicit.cusp import elementary_pole_patch, elementary_pole_forward

x28,l28 = sp.symbols('x28 l28',positive=True)
f28 = x28-1/x28
soft28 = sp.log(1+sp.exp(f28))
root28 = (l28+sp.sqrt(l28**2+4))/2
rep('C28a softplus quadratic inverse and strict derivative factor',
    sp.simplify(sp.exp(soft28)-1-sp.exp(f28))==0
    and sp.simplify(root28*(root28-l28)-1)==0
    and sp.simplify(sp.diff(f28,x28)-(1+x28**-2))==0)

coord28 = sp.symbols('coord28',positive=True)
outer_scale28 = sp.symbols('outer_scale28',positive=True)
back28 = outer_scale28*(1-1/coord28)
rep('C28a2 radial endpoint coordinate inverse identities',
    sp.simplify(outer_scale28/(outer_scale28-back28)-coord28)==0
    and sp.simplify(outer_scale28*(1-1/sp.log(sp.exp(coord28)))-back28)==0)

c28,t28,q28,d28 = sp.symbols('c28 t28 q28 d28',positive=True)
ell28 = sp.log(1+q28*d28)
shifted28 = c28*(1-1/(c28/t28+ell28))
A28,B28 = sp.exp(3),sp.E
raw_x28 = 1/q28-A28
outer_d28 = B28-A28-q28/(1-A28*q28)\
            +sp.log(1+sp.exp(-1/q28+A28+q28/(1-A28*q28)))
soft_large28 = raw_x28-1/raw_x28\
               +sp.log(1+sp.exp(-raw_x28+1/raw_x28))
rep('C28b exact outer-seam expression and flat-factor identity',
    sp.simplify(soft_large28+B28-1/q28-outer_d28)==0
    and sp.simplify(shifted28-(c28-t28)-t28**2*ell28/(c28+t28*ell28))==0)

xp28,u28,D28 = sp.symbols('xp28 u28 D28',positive=True)
Fdiff28 = (xp28-1/xp28)-(x28-1/x28)
W28 = sp.symbols('W28',real=True)
xw28 = x28+x28**2*u28*W28
rep('C28c inner-seam conjugation identity',
    sp.simplify(Fdiff28-(xp28-x28)*(1+1/(x28*xp28)))==0
    and sp.simplify(((xw28-1/xw28)-(x28-1/x28))/u28
                    -W28*(1+x28**2*(1+x28*u28*W28))
                    /(1+x28*u28*W28))==0)

eps28 = sp.Rational(1,10**8)
a28,c28r,b28,R28 = eps28,3*eps28/2,2*eps28,sp.Rational(1,10**4)
alpha28 = sp.Rational(1,10**5)
cap28 = sp.Rational(17,10**7)
rho28 = sp.Rational(1,10**7)
rep('C28d raw-collar and source-domain rational margins',
    b28-80*b28**2>c28r
    and sp.simplify(R28**sp.Rational(1,2)-200*R28**sp.Rational(3,2)
                    -sp.Rational(98,10**4))==0
    and 5000*rho28<1
    and rho28**sp.Rational(1,2)-200*rho28**sp.Rational(3,2)>sp.sqrt(c28r)
    and alpha28/(3*(2-alpha28/2))<cap28
    and alpha28/21>b28
    and cap28+80*cap28**2<sp.Rational(2,10**6)
    and sp.Rational(2,10**6)<sp.Rational(98,10**4)**2)

ur28,ui28,vr28,vi28 = sp.symbols('ur28 ui28 vr28 vi28',real=True)
den28 = 1+ur28
stereo28 = (ui28/den28,vr28/den28-1,vi28/den28)
norm28 = sum(v*v for v in stereo28)
rep('C28h rational north stereographic inversion on the unit S3',
    sp.simplify(den28**2*(norm28-2*(1-vr28)/den28)
                -(ur28**2+ui28**2+vr28**2+vi28**2-1))==0)

h28 = sp.symbols('h28',nonnegative=True)
X28 = (1-h28)/sp.sqrt(1+h28**2)
rep('C28i forward cap sphere and inverse denominator identities',
    sp.simplify(X28**2+2*h28/(1+h28**2)-1)==0
    and sp.simplify((2-X28**2)*(1+h28**2)-(1+h28)**2)==0
    and sp.simplify(1+X28*(1+h28)/sp.sqrt(1+h28**2)
                    -2/(1+h28**2))==0)

with mp.workdps(85):
    epsilon28 = mp.mpf('1e-8')
    inner28,outer28 = epsilon28,3*epsilon28/2
    def soft_mp28(x):
        f = x-1/x
        return mp.log1p(mp.exp(f)) if f<0 else f+mp.log1p(mp.exp(-f))
    def soft_inv_mp28(y):
        L = mp.log(mp.expm1(y))
        return 2/(mp.sqrt(L*L+4)-L) if L<0 else (L+mp.sqrt(L*L+4))/2
    def scale_mp28(r):
        if r>=outer28:
            return r
        x = mp.exp(outer28/(outer28-r))-mp.exp(3)
        y = soft_mp28(x)
        ell = mp.log1p(y/mp.e)
        return outer28*ell/(1+ell)
    def scale_inv_mp28(s):
        if s>=outer28:
            return s
        y = mp.exp(outer28/(outer28-s))-mp.e
        x = soft_inv_mp28(y)
        ell = mp.log1p(x/mp.exp(3))
        return inner28+outer28*ell/(3*(3+ell))
    def radial_mp28(pair,inverse=False):
        r = sum(abs(v)**2 for v in pair)
        s = scale_inv_mp28(r) if inverse else scale_mp28(r)
        factor = mp.sqrt(s/r)
        return (factor*pair[0],factor*pair[1])
    def raw_forward_mp28(pair,inf):
        base = forward_local23(pair,inf)
        scale = mp.pi if inf else 2*mp.pi
        pL = (mp.matrix([base[2],-base[1],base[0]])/scale if inf else
              mp.matrix([base[2],-base[0],-base[1]])/scale)
        rL = mp.norm(pL)
        nL = pL/rL
        projector = mp.matrix([[1+nL[2],nL[0]+mp.j*nL[1]],
                               [nL[0]-mp.j*nL[1],1-nL[2]]])/2
        projected = projector*mp.matrix(pair)
        z = mp.sqrt(rL)*projected/mp.norm(projected)
        return (z[0],z[1])
    def patched_mp28(pair,inf,inverse=False):
        r = sum(abs(v)**2 for v in pair)
        if r<=inner28:
            return pair
        compressed = radial_mp28(pair)
        middle = (raw_forward_mp28(compressed,inf) if inverse else
                  elliptic_pole_raw_inverse(compressed,inf,dps=80,elementary=True))
        return radial_mp28(middle,inverse=True)
    inverse_error28 = mp.mpf(0)
    raw_error28 = mp.mpf(0)
    phase_error28 = mp.mpf(0)
    production_error28 = mp.mpf(0)
    max_source28 = mp.mpf(0)
    cases28 = 0
    for inf28 in (False,True):
        for multiplier28 in map(mp.mpf,('0.5','1','1.001','1.01','1.05','1.1',
                                        '1.25','1.4','1.49','1.8','2','3','10','100','170')):
            radius28 = multiplier28*epsilon28
            for direction28 in ((1,0,0,0),(0,0,1,0),(1,2,3,-1),(-2,1,0,3)):
                unit28 = mp.matrix(direction28)/mp.sqrt(sum(v*v for v in direction28))
                z28 = (mp.mpc(unit28[0],unit28[1])*mp.sqrt(radius28),
                       mp.mpc(unit28[2],unit28[3])*mp.sqrt(radius28))
                output28 = patched_mp28(z28,inf28)
                max_source28 = max(max_source28,sum(abs(v)**2 for v in output28))
                recovered28 = patched_mp28(output28,inf28,inverse=True)
                inverse_error28 = max(inverse_error28,
                    mp.norm(mp.matrix(recovered28)-mp.matrix(z28))/mp.sqrt(radius28))
                production28 = elementary_pole_patch(z28,inf28,dps=80)
                try:
                    production_inverse28 = elementary_pole_forward(output28,inf28,dps=80)
                except Exception as exc:
                    raise AssertionError(
                        f'C28g production forward failed at infinity={inf28}, '
                        f'multiplier={multiplier28}, direction={direction28}') from exc
                production_error28 = max(
                    production_error28,
                    mp.norm(mp.matrix(production28)-mp.matrix(output28))/mp.sqrt(radius28),
                    mp.norm(mp.matrix(production_inverse28)-mp.matrix(recovered28))/mp.sqrt(radius28))
                if multiplier28>=2:
                    raw28 = elliptic_pole_raw_inverse(z28,inf28,dps=80,elementary=True)
                    raw_error28 = max(raw_error28,
                        mp.norm(mp.matrix(output28)-mp.matrix(raw28))/mp.sqrt(radius28))
                if multiplier28 in (mp.mpf('1.25'),mp.mpf('1.8'),mp.mpf('100')):
                    phase28 = mp.exp(mp.j*mp.mpf('.37'))
                    phase_output28 = patched_mp28(tuple(phase28*v for v in z28),inf28)
                    phase_error28 = max(phase_error28,
                        mp.norm(mp.matrix(phase_output28)
                                -phase28*mp.matrix(output28))/mp.sqrt(radius28))
                cases28 += 1
    rep(f'C28e {cases28} two-cap elementary patch/inverse compositions: '
        f'{float(inverse_error28):.1e}; raw collar {float(raw_error28):.1e}; '
        f'phase {float(phase_error28):.1e}; max source radius² '
        f'{float(max_source28):.7e}',
        inverse_error28<mp.mpf('1e-55') and raw_error28<mp.mpf('1e-55')
        and phase_error28<mp.mpf('1e-55') and max_source28<mp.mpf('2e-6'))
    rep(f'C28g {cases28} production versus independent H_j evaluations: '
        f'{float(production_error28):.1e}',production_error28<mp.mpf('1e-55'))

    def real_vector28(pair):
        return mp.matrix((pair[0].real,pair[0].imag,pair[1].real,pair[1].imag))
    def complex_pair28(vector):
        return (mp.mpc(vector[0],vector[1]),mp.mpc(vector[2],vector[3]))
    def central_jacobian28(func,pair,step):
        base = real_vector28(pair)
        jac = mp.matrix(4,4)
        for k in range(4):
            direction = mp.matrix([int(k==j) for j in range(4)])
            plus = real_vector28(func(complex_pair28(base+step*direction)))
            minus = real_vector28(func(complex_pair28(base-step*direction)))
            for j in range(4):
                jac[j,k] = (plus[j]-minus[j])/(2*step)
        return jac
    derivative_error28 = mp.mpf(0)
    minimum_det28 = mp.inf
    derivative_cases28 = 0
    for inf28 in (False,True):
        for multiplier28 in map(mp.mpf,('1.1','1.25','1.8','2','100')):
            radius28 = multiplier28*epsilon28
            direction28 = mp.matrix((1,2,3,-1))/mp.sqrt(15)
            z28 = (mp.mpc(direction28[0],direction28[1])*mp.sqrt(radius28),
                   mp.mpc(direction28[2],direction28[3])*mp.sqrt(radius28))
            w28 = patched_mp28(z28,inf28)
            step28 = mp.sqrt(radius28)*mp.mpf('1e-7')
            jac28 = central_jacobian28(lambda v: patched_mp28(v,inf28),z28,step28)
            inv_jac28 = central_jacobian28(
                lambda v: patched_mp28(v,inf28,inverse=True),w28,step28)
            defect28 = inv_jac28*jac28-mp.eye(4)
            derivative_error28 = max(derivative_error28,
                                     max(abs(defect28[j,k]) for j in range(4)
                                         for k in range(4)))
            minimum_det28 = min(minimum_det28,mp.det(jac28))
            derivative_cases28 += 1
    rep(f'C28f {derivative_cases28} numerical two-cap inverse-Jacobian samples: '
        f'{float(derivative_error28):.1e}; min det {float(minimum_det28):.6e}',
        derivative_error28<mp.mpf('1e-7') and minimum_det28>0)

# C29: exact identities for the two-way terminal-cut height in (11.3)--(11.7).
u29,v29,mu29,La29 = sp.symbols('u29 v29 mu29 La29', positive=True)
F29 = (La29*u29**2+mu29*u29)/(1+u29)
Finv29 = 2*v29/(sp.sqrt((mu29-v29)**2+4*La29*v29)+mu29-v29)
rep('C29a quadratic height correction has the stated positive inverse',
    sp.simplify((La29*u29**2+(mu29-v29)*u29-v29).subs(u29,Finv29))==0
    and sp.factor(sp.diff(F29,u29)-(La29+(mu29-La29)/(1+u29)**2))==0)
vflat29 = sp.symbols('vflat29', positive=True)
tail29 = sp.factor(vflat29*(1+F29.subs(u29,1/vflat29-1))/La29)
convex29 = (1-vflat29)**2+(mu29/La29)*vflat29*(1-vflat29)+vflat29/La29
rep('C29b terminal correction is the convex flat-tail expression',
    sp.factor(tail29-convex29)==0)
Ep029,Epz029,gp029,d29 = sp.symbols('Ep029 Epz029 gp029 d29',positive=True)
mu_tuned29 = d29*Ep029/(gp029*Epz029)
rep('C29c tuned endpoint derivative of T composed with g is d',
    sp.cancel((Epz029*mu_tuned29/Ep029)*gp029-d29)==0)
z29,q29,Y29,k29,a29 = sp.symbols('z29 q29 Y29 k29 a29',positive=True)
rad29 = (z29+sp.sqrt(z29*z29+4*q29))/2
rep('C29d positive quadratic root, inverse and open height seam',
    sp.simplify(rad29**2-z29*rad29-q29)==0
    and sp.simplify((z29-sp.sqrt(z29*z29+4*q29))/2+q29/rad29)==0
    and sp.simplify(sp.exp(k29*(a29-sp.Rational(1,2)))
        -q29*sp.exp(-k29*(a29-sp.Rational(1,2)))
        -(q29/(q29*sp.exp(-k29*(a29-sp.Rational(1,2))))
          -q29*sp.exp(-k29*(a29-sp.Rational(1,2)))))==0)
# Derive b_R'=d(R/R')/da using only the elementary differential chain:
# E'=E L^2/c, L'=L^2/c, F'=D E L^2/c,
# M'=D E L^2/(c F), D'=D_x E L^2/c.
E29,L29,Fp29,M29,D29,Dx29,c29 = sp.symbols('E29 L29 Fp29 M29 D29 Dx29 c29',positive=True)
br29 = c29*Fp29*M29*(M29-1)/(D29*E29*L29**2)
chain29 = {E29:E29*L29**2/c29,L29:L29**2/c29,
           Fp29:D29*E29*L29**2/c29,
           M29:D29*E29*L29**2/(c29*Fp29),
           D29:Dx29*E29*L29**2/c29}
deriv29 = sum(sp.diff(br29,key)*val for key,val in chain29.items())
claimed29 = M29**2+M29-1-Fp29*M29*(M29-1)*(1+2/L29)/(E29*D29) \
            -Fp29*M29*(M29-1)*Dx29/D29**2
rep('C29e one-soft radial distortion derivative identity',
    sp.factor(deriv29-claimed29)==0)
ur29,ui29,vr29,vi29,lr29,li29 = sp.symbols(
    'ur29 ui29 vr29 vi29 lr29 li29',real=True)
xi29 = (ur29+sp.I*ui29,vr29+sp.I*vi29)
phase29 = lr29+sp.I*li29
hopf29 = lambda pair: sp.Matrix([sp.re(sp.expand(2*pair[0]*sp.conjugate(pair[1]))),
    sp.im(sp.expand(2*pair[0]*sp.conjugate(pair[1]))),
    sp.expand(pair[0]*sp.conjugate(pair[0])-pair[1]*sp.conjugate(pair[1]))])
scaled29 = hopf29(tuple(phase29*w for w in xi29))
source29 = hopf29(xi29)
rep('C29f forward phase recovery preserves the Hopf base and circle value',
    all(sp.factor(dif)==0 for dif in scaled29-(lr29**2+li29**2)*source29)
    and sp.expand((sp.conjugate(phase29*xi29[0])*phase29**2*xi29[0]
                   +sp.conjugate(phase29*xi29[1])*phase29**2*xi29[1])
                  -(lr29**2+li29**2)*phase29*(abs(xi29[0])**2+abs(xi29[1])**2))==0)

# The positive inverse of S and both logarithmic radial-coordinate inverses
# give the two-way R composition in (11.3b). The intervening h and G_* maps
# then cancel in the stated order in (11.4).
xg29,ellg29,cg29 = sp.symbols('xg29 ellg29 cg29',positive=True)
wg29 = xg29-1/xg29
sg29 = sp.log(1+sp.exp(wg29))
xback29 = (wg29+sp.sqrt(wg29**2+4))/2
rg29 = cg29-cg29/ellg29
rep('C29g positive one-soft inverse and both R radial-coordinate inverses',
    sp.simplify(sp.exp(sg29)-1-sp.exp(wg29))==0
    and sp.factor(wg29**2+4-(xg29+1/xg29)**2)==0
    and sp.simplify(xback29.subs(sp.sqrt(wg29**2+4),xg29+1/xg29)-xg29)==0
    and sp.simplify(cg29-cg29/sp.log(sp.exp(2)+sp.exp(ellg29)-sp.exp(2))-rg29)==0
    and sp.simplify(cg29-cg29/sp.log(sp.E+sp.exp(ellg29)-sp.E)-rg29)==0)

# Cusp torus restoration reverses the real Cramer coefficients and the
# deck-equivariant height shift without choosing a circle argument.
xh29,yh29,c2h29,c3h29,h29,b2h29,b3h29,t2h29,t3h29 = sp.symbols(
    'xh29 yh29 c2h29 c3h29 h29 b2h29 b3h29 t2h29 t3h29')
Ah29 = c2h29*t2h29+c3h29*t3h29
ph2h29,ph3h29 = c2h29-h29*b2h29,c3h29-h29*b3h29
Ch2h29,Ch3h29 = ph2h29+h29*b2h29,ph3h29+h29*b3h29
rep('C29h cusp forward torus coefficients reverse tau_f restoration',
    sp.simplify(Ch2h29-c2h29)==0 and sp.simplify(Ch3h29-c3h29)==0
    and sp.simplify(xh29*sp.exp(Ah29)*sp.exp(-Ah29)-xh29)==0
    and sp.simplify(yh29*sp.exp(-Ah29)*sp.exp(Ah29)-yh29)==0)

# A reproducible numerical margin for the independently proved <200 bound.
# This is a sample, not a proof of the inequality on the continuum.
with mp.workdps(55):
    radial_max29 = mp.mpf(0)
    radial_arg29 = mp.mpf(0)
    for radial_i29 in range(1,2000):
        radial_a29 = mp.mpf('0.1')+mp.mpf('0.1')*radial_i29/2000
        radial_L29 = mp.mpf('0.2')/(mp.mpf('0.2')-radial_a29)
        radial_E29 = mp.exp(radial_L29)
        radial_x29 = radial_E29-mp.exp(2)
        radial_w29 = radial_x29-1/radial_x29
        if radial_w29<0:
            radial_S29 = mp.log1p(mp.exp(radial_w29))
            radial_sigma29 = mp.exp(radial_w29)/(1+mp.exp(radial_w29))
        else:
            radial_S29 = radial_w29+mp.log1p(mp.exp(-radial_w29))
            radial_sigma29 = 1/(1+mp.exp(-radial_w29))
        radial_D29 = radial_sigma29*(1+1/radial_x29**2)
        radial_Dp29 = (radial_sigma29*(1-radial_sigma29)
                       *(1+1/radial_x29**2)**2
                       -2*radial_sigma29/radial_x29**3)
        radial_F29 = mp.e+radial_S29
        radial_M29 = mp.log(radial_F29)
        radial_derivative29 = (radial_M29**2+radial_M29-1
            -radial_F29*radial_M29*(radial_M29-1)
             *(1+2/radial_L29)/(radial_E29*radial_D29)
            -radial_F29*radial_M29*(radial_M29-1)
             *radial_Dp29/radial_D29**2)
        if abs(radial_derivative29)>radial_max29:
            radial_max29 = abs(radial_derivative29)
            radial_arg29 = radial_a29
rep(f'C29i 1999 radial-derivative samples (numerical): max '
    f'{float(radial_max29):.6f} at a={float(radial_arg29):.7f}; '
    f'margin to 200 is {float(200-radial_max29):.6f}', radial_max29<200)

# C30: second-jet calculus for the fixed-count projection and phase recurrences.
# These are local derivative identities and an induction step for any fixed
# positive counts. They do not identify the resulting jets with every branch of
# the specific global G_cl evaluator; that source correspondence and the chart
# seams require their own checks.
from extra.s4t2_explicit.cusp import elementary_closed_parameters
counts30 = elementary_closed_parameters()
L30 = counts30['projector_bound']
N30 = counts30['projection_steps']
Lg30 = counts30['phase_bound']
K30 = counts30['phase_steps']
count_ok30 = (L30 == 10**8000 and N30 == 4*L30+1
              and Lg30 == 32*(N30+1)*L30 and K30 == 8*Lg30+1
              and N30 > 0 and K30 > 0)
rep('C30a exact certified counts for the two fixed-count folds', count_ok30)

s30,t30 = sp.symbols('s30 t30', real=True)
u30 = sp.Symbol('u30', positive=True)
ui30,uj30,uij30 = sp.symbols('ui30 uj30 uij30', real=True)
lift30 = u30+ui30*s30+uj30*t30+uij30*s30*t30
unary30 = [
    (1/u30, -1/u30**2, 2/u30**3),
    (sp.sqrt(u30), 1/(2*sp.sqrt(u30)), -1/(4*u30**sp.Rational(3,2))),
    (sp.exp(u30), sp.exp(u30), sp.exp(u30)),
    (sp.log(u30), 1/u30, -1/u30**2),
    (sp.sin(u30), sp.cos(u30), -sp.sin(u30)),
    (sp.cos(u30), -sp.sin(u30), -sp.cos(u30)),
    (sp.atan(u30), 1/(1+u30**2), -2*u30/(1+u30**2)**2),
]
unary_ok30 = all(
    sp.simplify(sp.diff(f.subs(u30,lift30),s30).subs({s30:0,t30:0})-fp*ui30)==0
    and sp.simplify(sp.diff(f.subs(u30,lift30),s30,t30).subs({s30:0,t30:0})
                    -fp*uij30-fpp*ui30*uj30)==0
    for f,fp,fpp in unary30)
flat30 = sp.exp(-1/u30)/(sp.exp(-1/u30)+sp.exp(-1/(1-u30)))
flat_slope30 = flat30*(1-flat30)*(1/u30**2+1/(1-u30)**2)
flat_curvature30 = flat30*(1-flat30)*(
    (1-2*flat30)*(1/u30**2+1/(1-u30)**2)**2
    -2/u30**3+2/(1-u30)**3)
flat_ok30 = (sp.simplify(sp.diff(flat30,u30)-flat_slope30)==0
             and sp.simplify(sp.diff(flat30,u30,2)-flat_curvature30)==0)
rep('C30b elementary and open flat-step second-jet rules',
    unary_ok30 and flat_ok30)

# Complex conjugation is componentwise real differentiation. The one-component
# modulus identity below sums immediately over the ten embedding components.
a30,b30,ai30,bi30,aj30,bj30,aij30,bij30 = sp.symbols(
    'a30 b30 ai30 bi30 aj30 bj30 aij30 bij30', real=True)
z30 = a30+ai30*s30+aj30*t30+aij30*s30*t30 + sp.I*(
    b30+bi30*s30+bj30*t30+bij30*s30*t30)
z030,z_i30,z_j30,z_ij30 = (
    a30+sp.I*b30, ai30+sp.I*bi30,
    aj30+sp.I*bj30, aij30+sp.I*bij30)
r30 = sp.expand(z30*sp.conjugate(z30))
modulus_ok30 = (
    sp.expand(sp.diff(r30,s30).subs({s30:0,t30:0})
              -2*sp.re(sp.conjugate(z030)*z_i30))==0
    and sp.expand(sp.diff(r30,s30,t30).subs({s30:0,t30:0})
                  -2*sp.re(sp.conjugate(z_j30)*z_i30
                           +sp.conjugate(z030)*z_ij30))==0)
p30,pi30,pj30,pij30,ri30,rj30,rij30 = sp.symbols(
    'p30 pi30 pj30 pij30 ri30 rj30 rij30', real=True)
r030 = sp.Symbol('r030', positive=True)
p_lift30 = p30+pi30*s30+pj30*t30+pij30*s30*t30
r_lift30 = r030+ri30*s30+rj30*t30+rij30*s30*t30
normal_mixed30 = (pij30/sp.sqrt(r030)
    -(pi30*rj30+pj30*ri30+p30*rij30)/(2*r030**sp.Rational(3,2))
    +3*p30*ri30*rj30/(4*r030**sp.Rational(5,2)))
normal_ok30 = sp.simplify(
    sp.diff(p_lift30/sp.sqrt(r_lift30),s30,t30).subs({s30:0,t30:0})
    -normal_mixed30)==0
rep('C30c complex modulus and nonzero normalized-projection second jets',
    modulus_ok30 and normal_ok30)

# One matrix-vector entry of p_j=P_j e_(j-1); matrix multiplication is a sum
# of these entries. The normalized p_j rule above therefore preserves symmetric
# second jets at every projection step, by induction for the exact N30 count.
P30,Pi30,Pj30,Pij30,e30,ei30,ej30,eij30 = sp.symbols(
    'P30 Pi30 Pj30 Pij30 e30 ei30 ej30 eij30', real=True)
P_lift30 = P30+Pi30*s30+Pj30*t30+Pij30*s30*t30
e_lift30 = e30+ei30*s30+ej30*t30+eij30*s30*t30
project_mixed30 = Pij30*e30+Pi30*ej30+Pj30*ei30+P30*eij30
project_ok30 = sp.expand(
    sp.diff(P_lift30*e_lift30,s30,t30).subs({s30:0,t30:0})
    -project_mixed30)==0

# The phase increment is atan(Im R/Re R), R=g_l conjugate(g_(l-1)).
# Its right-half-plane denominator is certified separately in (11.14)--(11.16).
ar30,br30,ari30,bri30,arj30,brj30,arij30,brij30 = sp.symbols(
    'ar30 br30 ari30 bri30 arj30 brj30 arij30 brij30', real=True)
Rreal30 = ar30+ari30*s30+arj30*t30+arij30*s30*t30
Rimag30 = br30+bri30*s30+brj30*t30+brij30*s30*t30
den30 = ar30**2+br30**2
num_i30 = ar30*bri30-br30*ari30
arg_mixed30 = ((arj30*bri30+ar30*brij30-brj30*ari30-br30*arij30)/den30
    -2*num_i30*(ar30*arj30+br30*brj30)/den30**2)
arg_ok30 = sp.factor(
    sp.diff(sp.atan(Rimag30/Rreal30),s30,t30).subs({s30:0,t30:0})
    -arg_mixed30)==0
arg_symmetric30 = sp.factor(
    arg_mixed30-arg_mixed30.xreplace({
        ari30:arj30, arj30:ari30, bri30:brj30, brj30:bri30}))==0
rep('C30d each projection and phase step preserves symmetric second jets',
    project_ok30 and normal_ok30 and arg_ok30 and arg_symmetric30)

# Negative control: an arbitrary invertible frame still squares to -Id but
# need not have zero Nijenhuis tensor. This prevents interpreting C30's
# symmetric-coordinate-jet induction as a claim about arbitrary frame fields.
q30 = sp.Symbol('q30', real=True)
E0030 = sp.zeros(dim)
E0030[0,0] = 1
frame30 = sp.eye(dim)+q30*E0030
frame_J30 = frame30*J6*frame30.inv()
frame_square30 = sp.simplify(frame_J30*frame_J30+sp.eye(dim))==sp.zeros(dim)
frame_dJ30 = [sp.zeros(dim) for _ in range(dim)]
frame_dJ30[2] = frame_J30.diff(q30).subs(q30,0)
frame_N30 = sp.expand(sum(
    J6[p,0]*frame_dJ30[p][1,3]-J6[p,3]*frame_dJ30[p][1,0]
    -J6[1,p]*(frame_dJ30[0][p,3]-frame_dJ30[3][p,0])
    for p in range(dim)))
rep('C30e non-coordinate-frame control: J^2=-Id yet N^1_(0,3)=-1',
    frame_square30 and frame_N30 == -1)

# Bind the recurrence checks above to the actual two large loop sites in the
# scalar evaluator. This is a source-level correspondence check for those sites,
# not a symbolic execution of the entire G_cl call graph or its chart seams.
import ast
import inspect
from extra.s4t2_explicit import cusp as cusp30


def source_tree30(function):
    return ast.parse(inspect.getsource(function)).body[0]


def matching_assignment30(tree,name,expression):
    expected = ast.dump(ast.parse(expression,mode='eval').body,
                        include_attributes=False)
    return any(isinstance(node,ast.Assign)
               and any(isinstance(target,ast.Name) and target.id==name
                       for target in node.targets)
               and ast.dump(node.value,include_attributes=False)==expected
               for node in ast.walk(tree))


def matching_call30(tree,expression):
    expected = ast.dump(ast.parse(expression,mode='eval').body,
                        include_attributes=False)
    return any(isinstance(node,ast.Call)
               and ast.dump(node,include_attributes=False)==expected
               for node in ast.walk(tree))


global_tree30 = source_tree30(cusp30.elementary_global_coordinates)
slice_tree30 = source_tree30(cusp30.elementary_global_slice)
correct_tree30 = source_tree30(cusp30._elementary_corrected_section)
project_tree30 = source_tree30(cusp30.elementary_projection_section)
phase_tree30 = source_tree30(cusp30.elementary_phase_family)
# The global function unpacks the slice return, so inspect that assignment
# separately rather than treating the tuple target as a scalar name.
source_counts30 = (
    matching_call30(global_tree30,'elementary_global_slice(Z,X,**kwargs)')
    and matching_call30(global_tree30,'elementary_restore_torus(kind,data,phi)')
    and matching_assignment30(slice_tree30,'params','elementary_closed_parameters()')
    and matching_assignment30(slice_tree30,'N',
        "params['projection_steps'] if projection_steps is None else projection_steps")
    and matching_assignment30(slice_tree30,'K',
        "params['phase_steps'] if phase_steps is None else phase_steps")
    and matching_call30(slice_tree30,
        '_elementary_corrected_section(xi,X,N,K,dps,diagnostics)')
    and matching_call30(correct_tree30,
        'elementary_projection_section(xi,X,N,dps,diagnostics)')
    and matching_call30(correct_tree30,
        'elementary_phase_family(mismatch,n,X,anchor,K,diagnostics=diagnostics)')
    and matching_call30(correct_tree30,
        'elementary_projection_section(xx,h,N,dps,diagnostics)')
    and count_ok30)
rep('C30f exact default counts reach both source recurrence sites',source_counts30)

project_loop30 = any(
    isinstance(node,ast.For)
    and ast.dump(node.iter,include_attributes=False)
        ==ast.dump(ast.parse('range(steps+1)',mode='eval').body,
                    include_attributes=False)
    for node in ast.walk(project_tree30))
project_source30 = (
    project_loop30
    and matching_assignment30(project_tree30,'phase','_elementary_inner(ref,vector)')
    and matching_assignment30(project_tree30,'vector',
                              'tuple(a*phase/abs(phase) for a in ref)')
    and matching_call30(project_tree30,
        'elliptic_raw_base_reference(tuple((1-h*h)*a for a in n),h,dps=dps,elementary=True)')
    and matching_call30(project_tree30,
        'elliptic_bundle_embedding(kind,data,dps=dps,elementary=True)'))
phase_loop30 = any(
    isinstance(node,(ast.For,ast.comprehension))
    and ast.dump(node.iter,include_attributes=False)
        ==ast.dump(ast.parse('range(steps+1)',mode='eval').body,
                    include_attributes=False)
    for node in ast.walk(phase_tree30))
phase_source30 = (
    phase_loop30
    and matching_assignment30(phase_tree30,'ratio','value*mp.conj(previous)')
    and matching_call30(phase_tree30,'mp.atan(ratio.imag/ratio.real)'))
rep('C30g source projection/phase updates match the checked jet primitives',
    project_source30 and phase_source30)

# If ref is unit and phase=ref*vector, P=ref ref* gives
# P vector=ref phase and |P vector|=|phase|. The source update above is
# therefore exactly normalized projection, not an unrelated phase recurrence.
ar1,ai1,ar2,ai2,vr1,vi1,vr2,vi2 = sp.symbols(
    'ar1 ai1 ar2 ai2 vr1 vi1 vr2 vi2', real=True)
ref30 = sp.Matrix([ar1+sp.I*ai1,ar2+sp.I*ai2])
vec30 = sp.Matrix([vr1+sp.I*vi1,vr2+sp.I*vi2])
phase_value30 = (ref30.conjugate().T*vec30)[0]
Pvec30 = ref30*ref30.conjugate().T*vec30
projection_value30 = all(sp.expand(Pvec30[i]-ref30[i]*phase_value30)==0
                         for i in range(2))
projection_norm30 = sp.expand(
    (Pvec30.conjugate().T*Pvec30)[0]
    -sp.conjugate(phase_value30)*phase_value30
     *(ref30.conjugate().T*ref30)[0])==0
right30 = sp.Symbol('right30', positive=True)
imag30 = sp.Symbol('imag30', real=True)
angle30 = sp.atan(imag30/right30)
argument_value30 = (
    sp.simplify(sp.cos(angle30)-right30/sp.sqrt(right30**2+imag30**2))==0
    and sp.simplify(sp.sin(angle30)-imag30/sp.sqrt(right30**2+imag30**2))==0
    and sp.simplify(sp.diff(angle30,right30)
                    +imag30/(right30**2+imag30**2))==0
    and sp.simplify(sp.diff(angle30,imag30)
                    -right30/(right30**2+imag30**2))==0)
rep('C30h actual projection and right-half-plane phase values/derivatives',
    projection_value30 and projection_norm30 and argument_value30)

# Here h is the inverse-direction coordinate map G_cl and M=Dh. Thus
# J=M^{-1} J0 M, in contrast with C10's forward-map convention. Using
# D(M^{-1})=-M^{-1}(DM)M^{-1} gives, for H_k=partial_k M,
# D_k J=M^{-1}(J0 H_k-H_k J). A constant linear change of the *domain*
# coordinates normalizes M to Id at a point without changing the Nijenhuis
# equations; the holomorphic target coordinates retain their standard J0.
# For a coordinate map H^a_{ij}=H^a_{ji}; test all 36 square entries and
# all 216 N components in this inverse convention.
dJ_inverse30 = [J6*sp.Matrix(dim,dim,lambda out,i:Hjet[out][i][k])
                -sp.Matrix(dim,dim,lambda out,i:Hjet[out][i][k])*J6
                for k in range(dim)]
square_inverse30 = all(value==0 for value in J6*J6+sp.eye(dim))
dsquare_inverse30 = all(
    sp.expand(value)==0
    for DJ in dJ_inverse30 for value in DJ*J6+J6*DJ)
N_inverse30 = [sp.expand(sum(
    J6[p,i]*dJ_inverse30[p][out,j]
    -J6[p,j]*dJ_inverse30[p][out,i]
    -J6[out,p]*(dJ_inverse30[i][p,j]-dJ_inverse30[j][p,i])
    for p in range(dim)))
    for out in range(dim) for i in range(dim) for j in range(dim)]
rep('C30i inverse-map convention: 36 J^2 and 216 N components vanish',
    square_inverse30 and dsquare_inverse30
    and len(N_inverse30)==216 and all(value==0 for value in N_inverse30))

# The cusp part of tau_f^{-1} is a fixed-size block. Check that the actual
# restore and source-recovery functions use the same two cusp bases, Cramer
# formula, finite H_f and exponential coordinates. The source reduction is a
# separate deck-chart transition; these identities are in one open deck cell.
def returned_expr30(tree,expression):
    expected = ast.dump(ast.parse(expression,mode='eval').body,
                        include_attributes=False)
    return any(isinstance(node,ast.Return)
               and ast.dump(node.value,include_attributes=False)==expected
               for node in ast.walk(tree))


restore_tree30 = source_tree30(cusp30.elementary_restore_torus)
source_tree_cusp30 = source_tree30(cusp30._elementary_source_coordinates)
basis_tree30 = source_tree30(cusp30._elementary_cusp_basis)
coeff_tree30 = source_tree30(cusp30._elementary_real_coeffs)
height_tree30 = source_tree30(cusp30.scalar_finite_height)
restore_basis30 = next(node for node in ast.walk(restore_tree30)
                       if isinstance(node,ast.FunctionDef) and node.name=='cusp_data')
restore_coeff30 = next(node for node in ast.walk(restore_tree30)
                       if isinstance(node,ast.FunctionDef) and node.name=='coeffs')
source_cusp_ok30 = (
    ast.dump(restore_basis30.body[-1].value,include_attributes=False)
        ==ast.dump(basis_tree30.body[-1].value,include_attributes=False)
    and ast.dump(restore_coeff30.body[-1].value,include_attributes=False)
        ==ast.dump(coeff_tree30.body[-1].value,include_attributes=False)
    and matching_assignment30(restore_tree30,'H',
        'scalar_finite_height(x,y,mp.exp(-2) if inf else mp.mpf(1))')
    and matching_assignment30(restore_tree30,'C',
        'tuple(v+H*w for v,w in zip(phi,c))')
    and matching_assignment30(restore_tree30,'A','C[0]*a2+C[1]*a3')
    and returned_expr30(restore_tree30,
        'kind,(x*mp.exp(-A),y*mp.exp(A),C[0]*e2+C[1]*e3)')
    and matching_assignment30(source_tree_cusp30,'C',
        '_elementary_real_coeffs(w,e2,e3)')
    and matching_assignment30(source_tree_cusp30,'H',
        'scalar_finite_height(x0,y0,nu)+deck')
    and matching_assignment30(source_tree_cusp30,'phi',
        'tuple(C[j]-H*c[j] for j in range(2))'))
rep('C30j both cusp signs share the exact source/restore finite-height block',
    source_cusp_ok30)

tc_r30,tc_i30,C2_30,C3_30,Hc30 = sp.symbols(
    'tc_r30 tc_i30 C2_30 C3_30 Hc30',real=True)
tc30 = tc_r30+sp.I*tc_i30
qrad30 = tc_r30**2+tc_i30**2
basis_ok30 = True
for inf30 in (False,True):
    E2_30,E3_30,E4_30,A2_30,A3_30 = (
        (tc30**2-1,sp.I*(1+tc30**2),-2*tc30,tc30,sp.I*tc30)
        if inf30 else
        (1-tc30**2,sp.I*(1+tc30**2),2*tc30,tc30,-sp.I*tc30))
    det_c30 = sp.expand(sp.re(E2_30)*sp.im(E3_30)
                         -sp.re(E3_30)*sp.im(E2_30))
    sign_c30 = -1 if inf30 else 1
    W_c30 = C2_30*E2_30+C3_30*E3_30
    def cramer30(value):
        return ((sp.re(value)*sp.im(E3_30)-sp.re(E3_30)*sp.im(value))/det_c30,
                (sp.re(E2_30)*sp.im(value)-sp.re(value)*sp.im(E2_30))/det_c30)
    recovered30 = cramer30(W_c30)
    coeff_c30 = cramer30(E4_30)
    basis_ok30 &= (
        sp.factor(det_c30-sign_c30*(1-qrad30**2))==0
        and all(sp.cancel(recovered30[j]-C)==0
                for j,C in enumerate((C2_30,C3_30)))
        and sp.simplify(E4_30-coeff_c30[0]*E2_30-coeff_c30[1]*E3_30)==0)
    ph_c30 = (C2_30-Hc30*coeff_c30[0],
              C3_30-Hc30*coeff_c30[1])
    A_c30 = C2_30*A2_30+C3_30*A3_30
    xc30,yc30 = sp.symbols('xc30 yc30')
    basis_ok30 &= (
        all(sp.cancel(ph_c30[j]+Hc30*coeff_c30[j]-C)==0
                for j,C in enumerate((C2_30,C3_30)))
        and sp.simplify((xc30*sp.exp(-A_c30))
                        *(yc30*sp.exp(A_c30))-xc30*yc30)==0
        and sp.simplify(xc30*sp.exp(-A_c30)*sp.exp(A_c30)-xc30)==0
        and sp.simplify(yc30*sp.exp(A_c30)*sp.exp(-A_c30)-yc30)==0)
rep('C30k both cusp Cramer determinants and torus maps reverse exactly',
    basis_ok30)

# One genuinely noncap open branch of the actual four-seed H_f. If
# 1/16<A<1/2, 0<Y<1/16, q=AY, then qA<1/64 and qY<1/512; no seed is
# reciprocal, the A seed uses the flat interpolation, and the other three
# seeds are in their specified linear branch. It is outside the pole cap:
# nu^2|x|^2+|y|^2 >= A > 1/16. Explicit interior witnesses with Y=1/100
# are |x|^2=1/8 at zero and |x|^2=8 at infinity. Their respective |xy|^2
# are 1/800 and 2/25, both below the cusp-chart cutoff 9/100. For the first,
# A=(1/8)exp(-4/799) lies between (1/8)(1-4/799) and 1/8. For the second,
# A=8exp(-4-8/23) lies between 32/385 and 128/625, using e^4<55,
# 5/2<e<3, and exp(8/23)<sqrt(e)<7/4. Both intervals lie in (1/16,1/2).
seed_source30 = (
    matching_assignment30(height_tree30,'q','A*Y')
    and matching_assignment30(height_tree30,'weight',
        'scalar_flat_step(mp.log(16*v)/mp.log(8))')
    and matching_assignment30(height_tree30,'value',
        'v+weight*(v/(1+v)-v)')
    and returned_expr30(height_tree30,
        'seed(A)+seed(q*A)-seed(Y)-seed(q*Y)+q*q*(A-Y)/(1-q)'))
Ah30,Yh30 = sp.symbols('Ah30 Yh30',positive=True)
qh30 = Ah30*Yh30
uh30 = sp.log(16*Ah30)/sp.log(8)
step_h30 = sp.exp(-1/uh30)/(sp.exp(-1/uh30)+sp.exp(-1/(1-uh30)))
seed_A30 = Ah30+step_h30*(Ah30/(1+Ah30)-Ah30)
tail_h30 = qh30*Ah30-Yh30-qh30*Yh30+qh30**2*(Ah30-Yh30)/(1-qh30)
height_branch30 = seed_A30+tail_h30
height_jet30 = (sp.diff(height_branch30,Ah30,Yh30)
                -sp.diff(height_branch30,Yh30,Ah30))
height_ok30 = (seed_source30 and sp.simplify(height_jet30)==0
               and sp.simplify(sp.diff(height_branch30,Ah30,Yh30)
                               -sp.diff(tail_h30,Ah30,Yh30))==0
               and sp.Rational(1,2)**2*sp.Rational(1,16)<sp.Rational(1,16)
               and sp.Rational(1,2)*sp.Rational(1,16)**2<sp.Rational(1,16)
               and sp.Rational(1,8)*(1-sp.Rational(4,799))>sp.Rational(1,16)
               and sp.Rational(32,385)>sp.Rational(1,16)
               and sp.Rational(128,625)<sp.Rational(1,2))
rep('C30l source four-seed H_f: explicit noncap branch has symmetric Hessian',
    height_ok30)

# Mixed jets of the exact cusp restoration, for arbitrary symmetric jets of
# its inputs. Each complex identity covers its real and imaginary components.
# The Cramer denominator above is nonzero on |t|<1; hence its coefficients
# have such jets, as does the fixed open branch of H_f just checked.
def lift2_30(name):
    v,vs,vt,vst = sp.symbols(f'{name} {name}s {name}t {name}st')
    return v+vs*s30+vt*t30+vst*s30*t30


def at0_30(expr):
    return expr.subs({s30:0,t30:0})


xj30,yj30,Hj30,ph2j30,ph3j30 = (
    lift2_30(name) for name in ('xj30','yj30','Hj30','ph2j30','ph3j30'))
c2j30,c3j30,a2j30,a3j30,e2j30,e3j30 = (
    lift2_30(name) for name in
    ('c2j30','c3j30','a2j30','a3j30','e2j30','e3j30'))
C2j30,C3j30 = ph2j30+Hj30*c2j30,ph3j30+Hj30*c3j30
Aj30 = C2j30*a2j30+C3j30*a3j30
Z1j30,Z2j30,Wj30 = (xj30*sp.exp(-Aj30),yj30*sp.exp(Aj30),
                     C2j30*e2j30+C3j30*e3j30)
def jet_mixed30(expr):
    return at0_30(sp.diff(expr,s30,t30))


def jet_first30(expr,direction):
    return at0_30(sp.diff(expr,direction))


def product_mixed30(left,right):
    return (jet_mixed30(left)*at0_30(right)
            +jet_first30(left,s30)*jet_first30(right,t30)
            +jet_first30(left,t30)*jet_first30(right,s30)
            +at0_30(left)*jet_mixed30(right))


C2_mixed30 = (jet_mixed30(ph2j30)+product_mixed30(Hj30,c2j30))
C3_mixed30 = (jet_mixed30(ph3j30)+product_mixed30(Hj30,c3j30))
structural_c30 = (
    sp.expand(jet_mixed30(C2j30)-C2_mixed30)==0
    and sp.expand(jet_mixed30(C3j30)-C3_mixed30)==0
    and sp.expand(jet_mixed30(Aj30)
                  -product_mixed30(C2j30,a2j30)
                  -product_mixed30(C3j30,a3j30))==0)
A0j30,Asj30,Atj30,Astj30 = (
    at0_30(Aj30),jet_first30(Aj30,s30),
    jet_first30(Aj30,t30),jet_mixed30(Aj30))
Z1_mixed30 = sp.exp(-A0j30)*(
    jet_mixed30(xj30)-jet_first30(xj30,s30)*Atj30
    -jet_first30(xj30,t30)*Asj30
    +at0_30(xj30)*(Asj30*Atj30-Astj30))
Z2_mixed30 = sp.exp(A0j30)*(
    jet_mixed30(yj30)+jet_first30(yj30,s30)*Atj30
    +jet_first30(yj30,t30)*Asj30
    +at0_30(yj30)*(Asj30*Atj30+Astj30))
structural_c30 &= (
    sp.simplify(jet_mixed30(Z1j30)-Z1_mixed30)==0
    and sp.simplify(jet_mixed30(Z2j30)-Z2_mixed30)==0
    and sp.expand(jet_mixed30(Wj30)
                  -product_mixed30(C2j30,e2j30)
                  -product_mixed30(C3j30,e3j30))==0
    and all(sp.simplify(sp.diff(expr,s30,t30)
                        -sp.diff(expr,t30,s30))==0
            for expr in (Z1j30,Z2j30,Wj30)))
rep('C30m cusp restoration has explicit symmetric second jets in all six outputs',
    structural_c30 and height_ok30 and basis_ok30)
rep('C30n cusp-restoration J^2/N via its coordinate Hessian (round-to-W separate)',
    source_cusp_ok30 and basis_ok30 and height_ok30 and structural_c30
    and square_inverse30 and all(value==0 for value in N_inverse30))

# C31: the actual meridian projector on one pure-middle open chart.
# Assume 1/2<|t|<1 and pi/2<Arg(t)<3pi/4, away from the argument cuts.
# Both cusp radial cutoffs vanish; the raw reference is ('mid',t,0,x4,log t).
# Arg(-t)+pi=Arg(t), so the two surviving embedding entries have the *same*
# phase, independently of x4. Their weights are (d,1), where
# d=1-s((Arg(t)-pi/2)/(pi/4)). The projector has only the last 2x2 block
# [[d^2,d],[d,1]]/(1+d^2). This specializes P_j=P(n,h_j), rather than an
# arbitrary smooth rank-one projector. The path from the south through the
# cusp regions and its subsequent phase correction remain separate.
def matching_statement31(tree,statement):
    expected = ast.dump(ast.parse(statement).body[0],include_attributes=False)
    return any(ast.dump(node,include_attributes=False)==expected
               for node in ast.walk(tree))


meridian_tree31 = source_tree30(cusp30.elliptic_meridian_projector)
raw_tree31 = source_tree30(cusp30.elliptic_raw_base_reference)
bundle_tree31 = source_tree30(cusp30.elliptic_bundle_embedding)
middle_source31 = (
    matching_assignment30(meridian_tree31,'v',
        'tuple((1-h*h)*mp.mpf(a) for a in n)')
    and matching_call30(meridian_tree31,
        'elliptic_raw_base_reference(v,h,dps=dps,elementary=elementary)')
    and matching_call30(meridian_tree31,
        'elliptic_bundle_embedding(kind,data,dps=dps,elementary=elementary)')
    and matching_call30(meridian_tree31,
        'mp.matrix([[a*mp.conj(b) for b in e] for a in e])')
    and matching_assignment30(raw_tree31,'aa','xx*mp.sqrt(2-xx*xx)')
    and matching_assignment30(raw_tree31,'u1',
        'mp.mpc(1-s2,2*y[0])/(1+s2)')
    and matching_assignment30(raw_tree31,'u2',
        'mp.mpc(2*y[1],2*y[2])/(1+s2)')
    and matching_assignment30(raw_tree31,'inverse_ratio',
        'abs(u1)>abs(u2)')
    and matching_assignment30(raw_tree31,'ratio',
        'u2/u1 if inverse_ratio else u1/u2')
    and matching_assignment30(raw_tree31,'small',
        'ratio*mp.sqrt((2*B if inverse_ratio else 2)/(b+root))')
    and matching_assignment30(raw_tree31,'inf',
        "inverse_ratio and abs(small)<mp.mpf('.5')")
    and matching_assignment30(raw_tree31,'t',
        'None if inf else (1/small if inverse_ratio else small)')
    and returned_expr30(raw_tree31,"'mid',(t,mp.mpf(0),four,lt)")
    and matching_assignment30(bundle_tree31,'b0',
        "1-scalar_flat_step((r-mp.mpf('.3'))/mp.mpf('.2'))")
    and matching_assignment30(bundle_tree31,'bi',
        "1-scalar_flat_step((1/r-mp.mpf('.3'))/mp.mpf('.2'))")
    and matching_assignment30(bundle_tree31,'bm','1-b0-bi')
    and matching_assignment30(bundle_tree31,'angle',
        'mp.arg(t if index==0 else -t)')
    and matching_assignment30(bundle_tree31,'weight',
        '1-scalar_flat_step((abs(angle)-mp.pi/2)/(mp.pi/4))')
    and matching_statement31(bundle_tree31,
        'vec[8+index]=bm*weight*mp.exp(mp.j*(2*mp.pi*one'
        '-(angle+centre-lt.imag)*four))')
    and matching_assignment30(bundle_tree31,'norm2',
        'sum(abs(q)**2 for q in vec)')
    and matching_assignment30(bundle_tree31,'result',
        'tuple(+(q/mp.sqrt(norm2)) for q in vec)'))
rep('C31a actual raw-reference, embedding and meridian-projector source path',
    middle_source31)

# On 1/2<|t|<1 the dominant coordinate is u1: indeed
# |u1/u2|^2>1 iff B|t|^4>1, and B=e^4>16. The actual source therefore
# sets ratio=u2/u1, solves for small=1/t by its positive root, then takes
# t=1/small. This positive factor gives Arg(t)=Arg(u1/u2). If
# sigma=|1/t|^2, z=|u2/u1|^2=sigma(sigma+A)/(A sigma+B), and b=A(1-z),
# both sigma and the source's 2Bz/(b+sqrt(b^2+4Bz)) solve the quadratic.
radius31,Aprof31,Bprof31 = sp.symbols('radius31 Aprof31 Bprof31',positive=True)
zprof31 = radius31*(radius31+Aprof31)/(Aprof31*radius31+Bprof31)
bprof31 = Aprof31*(1-zprof31)
rootprof31 = sp.Symbol('rootprof31',positive=True)
inverse_radius31 = 2*Bprof31*zprof31/(bprof31+rootprof31)
profile31 = (
    sp.factor(zprof31-1
              -(radius31**2-Bprof31)
               /(Aprof31*radius31+Bprof31))==0
    and
    sp.factor(radius31**2+bprof31*radius31-Bprof31*zprof31)==0
    and sp.factor((inverse_radius31**2
                   +bprof31*inverse_radius31-Bprof31*zprof31)
                  *(bprof31+rootprof31)**2/(Bprof31*zprof31)
                  -(bprof31**2+4*Bprof31*zprof31-rootprof31**2))==0)

# The raw u1 conj(u2) argument is an explicit rational function of (n,X).
# In the stated upper-left angular sector its half-angle formula is smooth.
ym131,ym231,ym331 = sp.symbols('ym131 ym231 ym331',real=True)
ymnorm31 = ym131**2+ym231**2+ym331**2
u131 = (1-ymnorm31+2*sp.I*ym131)/(1+ymnorm31)
u231 = 2*(ym231+sp.I*ym331)/(1+ymnorm31)
realarg31 = (1-ymnorm31)*ym231+2*ym131*ym331
imagarg31 = 2*ym131*ym231-(1-ymnorm31)*ym331
arg_numerator31 = sp.expand(u131*sp.conjugate(u231)
    -2*(realarg31+sp.I*imagarg31)/(1+ymnorm31)**2)
arg_source31 = sp.simplify(arg_numerator31)==0
n131,n231,n331,X31 = sp.symbols('n131 n231 n331 X31',real=True)
target_den31 = 1-X31*sp.sqrt(2-X31**2)
target_factor31 = (1-X31**2)/target_den31
target_y31 = (target_factor31*n131,
              1-target_factor31*n331,-target_factor31*n231)
target_v31 = tuple((1-X31**2)*v for v in (n131,n231,n331))
target_aa31 = X31*sp.sqrt(2-X31**2)
target_w31 = (target_v31[0]/(1-target_aa31),
              -target_v31[2]/(1-target_aa31),
              -target_v31[1]/(1-target_aa31))
target_y_source31 = (target_w31[0],target_w31[1]+1,target_w31[2])
target_source31 = (
    matching_assignment30(raw_tree31,'w',
        '(vv[0]/(1-aa),-vv[2]/(1-aa),-vv[1]/(1-aa))')
    and matching_assignment30(raw_tree31,'y','(w[0],w[1]+1,w[2])')
    and matching_assignment30(raw_tree31,'aa','xx*mp.sqrt(2-xx*xx)')
    and all(sp.simplify(target_y_source31[j]-target_y31[j])==0
            for j in range(3)))

# An exact point in this open source branch is t=(-1+2i)/3. Write
# s=|t|^2=5/9 and R=|u1/u2|^2=s(A+Bs)/(1+As), with A=(4pi)^2>144
# and 16<B=e^4<64. Here pi>3 follows from an inscribed regular hexagon;
# 2<e<8/3+5/96=87/32<11/4 follows by bounding the series tail after
# 1/3! by (1/4!)/(1-1/5), and (11/4)^4<64. The rational inequalities
# below prove 1<R<2,
# so the source chooses inverse_ratio, obtains small=1/t (hence no `inf`),
# and enters its pure-middle branch. Choose u2>0 and u1/|u1|=t/|t|.
# Then |u1|^2,|u2|^2 lie in (1/3,2/3). For the inverse S3 chart,
# D=1+Re(u1)>1-sqrt(2/15)>3/5, while D<1. The resulting real stereographic
# w=(Im(u1)/D,u2/D-1,0) has 1/2<w[0]<4/3 and -1/2<w[1]<2/3.
# Thus 1/2<|w|<2 and X=(|w|-1)/sqrt(1+|w|^2) lies in (-1/2,1/2),
# since X'(r)=(r+1)/(1+r^2)^(3/2)>0 and its endpoint values are +/-1/sqrt(5).
# This supplies a nonempty open (n,X) cell, with all strict inequalities.
tw31 = (-1+2*sp.I)/3
sw31 = sp.Rational(5,9)
rw31 = sp.Symbol('rw31',positive=True)
Xw31 = (rw31-1)/sp.sqrt(1+rw31**2)
aw31 = Xw31*sp.sqrt(2-Xw31**2)
# Since 2-X^2=(r+1)^2/(1+r^2) and r>0, its positive square root is
# (r+1)/sqrt(1+r^2), giving a=(r^2-1)/(r^2+1).
aw_root31 = (rw31+1)/sp.sqrt(1+rw31**2)
aw_formula31 = (rw31**2-1)/(rw31**2+1)
target_scale_w31 = (1-Xw31**2)/(1-aw_formula31)
witness31 = (
    sp.simplify(tw31*sp.conjugate(tw31)-sw31)==0
    and sp.re(tw31)<0 and sp.im(tw31)>-sp.re(tw31)
    and sp.Rational(8,3)+sp.Rational(5,96)==sp.Rational(87,32)
    and sp.Rational(87,32)<sp.Rational(11,4)
    and sp.Rational(11,4)**4<64
    and sp.Rational(1,4)<sw31<1
    and 16*sw31**2>1
    and 64*sw31**2 < 2+144*sw31
    and sw31*144/(1+144*sw31)>sp.Rational(1,2)
    and sp.Rational(2,15)<sp.Rational(4,25)
    and sp.Rational(4,15)>sp.Rational(1,4)
    and sp.Rational(8,15)<sp.Rational(16,25)
    and sp.Rational(1,3)>sp.Rational(1,4)
    and sp.Rational(4,3)**2+sp.Rational(2,3)**2<4
    and sp.simplify(2-Xw31**2-aw_root31**2)==0
    and sp.simplify(Xw31*aw_root31-aw_formula31)==0
    and sp.simplify(sp.diff(Xw31,rw31)
                    -(rw31+1)/(1+rw31**2)**sp.Rational(3,2))==0
    and sp.simplify(Xw31.subs(rw31,sp.Rational(1,2))
                    +1/sp.sqrt(5))==0
    and sp.simplify(Xw31.subs(rw31,2)-1/sp.sqrt(5))==0
    and sp.Rational(1,5)<sp.Rational(1,4)
    and sp.simplify(target_scale_w31-rw31)==0)
rep('C31b positive profile, target argument and exact open-cell witness',
    profile31 and arg_source31 and target_source31 and witness31)

# On pi/2<theta<3pi/4, d1=1 and d0=d. The local projector is independent
# of the raw reference's x1=0 and x4. Its 2-jet is checked directly below.
d31 = sp.Symbol('d31',real=True)
Pmid31 = sp.Matrix([[d31*d31,d31],[d31,1]])/(1+d31*d31)
Pd31 = Pmid31.diff(d31)
Pdd31 = Pd31.diff(d31)
projector_mid31 = (
    sp.simplify(Pmid31*Pmid31-Pmid31)==sp.zeros(2)
    and Pmid31.T==Pmid31 and sp.simplify(sp.trace(Pmid31)-1)==0
    and sp.simplify(Pmid31*Pd31+Pd31*Pmid31-Pd31)==sp.zeros(2)
    and sp.simplify(Pmid31*Pdd31+Pdd31*Pmid31
                    +2*Pd31*Pd31-Pdd31)==sp.zeros(2))
theta31 = sp.Symbol('theta31',real=True)
stheta31 = (theta31-sp.pi/2)/(sp.pi/4)
flat_theta31 = sp.exp(-1/stheta31)/(
    sp.exp(-1/stheta31)+sp.exp(-1/(1-stheta31)))
dtheta31 = 1-flat_theta31
first_dtheta31 = -4*flat_slope30.subs(u30,stheta31)/sp.pi
second_dtheta31 = -16*flat_curvature30.subs(u30,stheta31)/sp.pi**2
flat_chain31 = (
    sp.simplify(sp.diff(dtheta31,theta31)-first_dtheta31)==0
    and sp.simplify(sp.diff(dtheta31,theta31,2)-second_dtheta31)==0)
theta_i31,theta_j31,theta_ij31 = sp.symbols(
    'theta_i31 theta_j31 theta_ij31',real=True)
d0_31,ds_31,dss_31 = sp.symbols('d0_31 ds_31 dss_31',real=True)
d_lift31 = (d0_31+ds_31*theta_i31*s30+ds_31*theta_j31*t30
            +(dss_31*theta_i31*theta_j31
              +ds_31*theta_ij31)*s30*t30)
Pjet31 = Pmid31.subs(d31,d_lift31)
expected_Pij31 = (Pd31.subs(d31,d0_31)
    *(dss_31*theta_i31*theta_j31+ds_31*theta_ij31)
    +Pdd31.subs(d31,d0_31)*ds_31**2*theta_i31*theta_j31)
middle_jet31 = all(sp.factor(jet_mixed30(Pjet31[i,j])
                                  -expected_Pij31[i,j])==0
                   for i in range(2) for j in range(2))
rep('C31c pure-middle projector and its complete local second-jet chain',
    projector_mid31 and middle_jet31 and flat_chain31)

# Generic exact two-jet product for P=VV*/(V*V). It extends the preceding
# 2x2 calculation entrywise to the actual ten-component embedding; positivity
# of V*V is (11.2), and C30c gives the normalized component second jet.
e031,es31,et31,est31,f031,fs31,ft31,fst31 = sp.symbols(
    'e031 es31 et31 est31 f031 fs31 ft31 fst31')
e_lift31 = e031+es31*s30+et31*t30+est31*s30*t30
f_lift31 = f031+fs31*s30+ft31*t30+fst31*s30*t30
entry31 = e_lift31*sp.conjugate(f_lift31)
entry_mixed31 = (est31*sp.conjugate(f031)
    +es31*sp.conjugate(ft31)+et31*sp.conjugate(fs31)
    +e031*sp.conjugate(fst31))
generic_projector31 = sp.expand(jet_mixed30(entry31)-entry_mixed31)==0
rep('C31d ten-component normalized projector has exact symmetric 2-jets',
    generic_projector31 and modulus_ok30 and normal_ok30)

# C32: a second actual meridian-projector generator, with *both* cusp pairs
# active.  This is the positive-height cusp cell
#       1/10 < H < 1/4,  |xy| < 3/10,  (x,y) != (0,0).
# The first pair has weight 1.  The second, unshifted pair has the nonconstant
# weight w=1-s((1/4-H)/(3/20)) in (0,1).  The source's H<0 deck-shift branch is
# inactive.  This checks a genuine noncap projector cell, for either cusp;
# it does not assert that all meridian steps lie in this one cell.
bundle_tree32 = source_tree30(cusp30.elliptic_bundle_embedding)
reference_tree32 = source_tree30(cusp30.elementary_cusp_reference)
def matching_if_test32(tree,test):
    expected = ast.dump(ast.parse(test,mode='eval').body,
                        include_attributes=False)
    return any(isinstance(node,ast.If)
               and ast.dump(node.test,include_attributes=False)==expected
               for node in ast.walk(tree))

source_cusp32 = (
    matching_assignment30(meridian_tree31,'e',
        'elliptic_bundle_embedding(kind,data,dps=dps,elementary=elementary)')
    and matching_statement31(meridian_tree31,
        'return mp.matrix([[a*mp.conj(b) for b in e] for a in e])')
    and matching_if_test32(bundle_tree32,"kind in ('c0','ci')")
    and matching_if_test32(bundle_tree32,
        "abs(x*y)>mp.mpf('.3') or abs(H)>mp.mpf('.6')")
    and matching_if_test32(bundle_tree32,'H<0 and w>0')
    and matching_assignment30(bundle_tree32,'nu',
        "mp.exp(-2) if kind=='ci' else mp.mpf(1)")
    and matching_assignment30(bundle_tree32,'H',
        'mp.mpf(data[2]) if elementary else elliptic_height_chart(x,y,nu,dps=dps)')
    and matching_statement31(bundle_tree32,
        "return 1-scalar_flat_step((abs(h)-mp.mpf('.25'))/mp.mpf('.15'))")
    and matching_statement31(bundle_tree32,
        'r = nu*nu*abs(x)**2+abs(y)**2')
    and matching_statement31(bundle_tree32,
        'return (nu*mp.conj(x)/mp.sqrt(r),y/mp.sqrt(r))')
    and matching_statement31(bundle_tree32,'block = [hw(H)*q for q in f]')
    and matching_statement31(bundle_tree32,
        "w = hw(H+mp.mpf('.5') if H<0 else H-mp.mpf('.5'))")
    and matching_statement31(bundle_tree32,'block.extend(w*q for q in f)')
    and matching_statement31(bundle_tree32,
        "start = 4 if kind=='ci' else 0")
    and matching_statement31(bundle_tree32,'vec[start:start+4] = block')
    and matching_statement31(bundle_tree32,
        'norm2 = sum(abs(q)**2 for q in vec)')
    and matching_statement31(bundle_tree32,
        'result = tuple(+(q/mp.sqrt(norm2)) for q in vec)')
    and matching_statement31(reference_tree32,
        "cell = int(mp.floor(a+mp.mpf('.5')))")
    and matching_statement31(reference_tree32,
        'f = inverse_height_moment(b,q,0)')
    and matching_statement31(reference_tree32,
        'x = mp.sqrt(A/e)/nu')
    and matching_statement31(reference_tree32,'return x,t/x'))
H32 = sp.Symbol('H32',real=True)
cut32 = (sp.Rational(1,4)-H32)/sp.Rational(3,20)
first_cut32 = (H32-sp.Rational(1,4))/sp.Rational(3,20)
# The strict cell inequalities imply 0<cut32<1, so this is exactly the
# source's open flat-step branch.  A source reference at t=1/10, a=1/5
# realizes the cell: q=nu^2 exp(-4/99)/100<1/8 for nu in {1,e^-2}, and
# the certified increasing moment in (11.4) gives a finite positive frame.
witness32 = (sp.simplify(cut32.subs(H32,sp.Rational(1,10))-1)==0
             and cut32.subs(H32,sp.Rational(1,4))==0
             and sp.simplify(first_cut32+cut32)==0
             and 0<cut32.subs(H32,sp.Rational(1,5))<1
             and sp.Rational(1,10)<sp.Rational(1,5)<sp.Rational(1,4)
             and sp.Rational(1,10)<sp.Rational(3,10)
             and sp.Rational(1,100)<sp.Rational(1,8)
             and bool(sp.exp(-sp.Rational(4,99))<1)
             and bool(sp.exp(-4)<1)
             and sp.floor(sp.Rational(1,5)+sp.Rational(1,2))==0)
# At that witness the certified f(a,q,0) is positive (strictly increasing
# with f(0,q,0)=0).  Its positive value gives the source's a>0 branch below;
# these two identities check the resulting nonzero frame.
f_w32,nu_w32 = sp.symbols('f_w32 nu_w32',positive=True)
t_w32 = sp.Rational(1,10)
e_w32 = sp.exp(-sp.Rational(4,99))
q_w32 = nu_w32**2*e_w32*t_w32**2
root_w32 = sp.sqrt(f_w32**2+4*q_w32)
A_w32,Y_w32 = (root_w32+f_w32)/2,(root_w32-f_w32)/2
x_w32 = sp.sqrt(A_w32/e_w32)/nu_w32
y_w32 = t_w32/x_w32
witness32 &= (sp.simplify(x_w32*y_w32-t_w32)==0
              and sp.simplify(A_w32*Y_w32-q_w32)==0)
rep('C32a both-cusp two-active-pair source path and exact open witness',
    source_cusp32 and witness32)

xr32,xi32,yr32,yi32 = sp.symbols('xr32 xi32 yr32 yi32',real=True)
nu32 = sp.Symbol('nu32',positive=True)  # 1 or exp(-2), constant in jets.
weight32 = sp.Symbol('weight32',real=True)
radius32 = (nu32**2*(xr32**2+xi32**2)+yr32**2+yi32**2)
spinor32 = sp.Matrix([nu32*(xr32-sp.I*xi32),yr32+sp.I*yi32])
norm32 = spinor32*spinor32.conjugate().T
Q32 = norm32/radius32
B32 = sp.Matrix([[1,weight32],[weight32,weight32**2]])/(1+weight32**2)
P32 = sp.kronecker_product(B32,Q32)
# A separate positive radius symbol makes conjugation of the source's
# normalizing square root unambiguous; then impose the manifest sum of
# squares radius.  The four entries occupy 0:4 (south) or 4:8 (north) in C^10.
positive_radius32 = sp.Symbol('positive_radius32',positive=True)
source_vector32 = sp.Matrix([*spinor32,weight32*spinor32[0],
                             weight32*spinor32[1]]) / sp.sqrt(
                                 positive_radius32*(1+weight32**2))
source_projector32 = sp.kronecker_product(
    B32,norm32/positive_radius32)
source_norm32 = (source_vector32.conjugate().T*source_vector32)[0]
source_value32 = (
    all(sp.simplify(v)==0 for v in
        source_vector32*source_vector32.conjugate().T-source_projector32)
    and source_projector32.subs(positive_radius32,radius32)==P32
    and sp.simplify(source_norm32-radius32/positive_radius32)==0
    and all(sp.cancel(v)==0 for v in P32*P32-P32)
    and all(sp.cancel(v)==0 for v in P32.conjugate().T-P32)
    and sp.cancel(sp.trace(P32)-1)==0)
rep('C32b actual four-slot projector B(weight) tensor Q(x,y), rank one',
    source_value32)

flat32 = sp.exp(-1/cut32)/(sp.exp(-1/cut32)
                             +sp.exp(-1/(1-cut32)))
actual_weight32 = 1-flat32
slope32 = flat32*(1-flat32)*(1/cut32**2+1/(1-cut32)**2)
curvature32 = flat32*(1-flat32)*(
    (1-2*flat32)*(1/cut32**2+1/(1-cut32)**2)**2
    -2/cut32**3+2/(1-cut32)**3)
weight_H32 = sp.Rational(20,3)*slope32
weight_HH32 = -sp.Rational(400,9)*curvature32
weight_jet32 = (
    sp.simplify(sp.diff(actual_weight32,H32)-weight_H32)==0
    and sp.simplify(sp.diff(actual_weight32,H32,2)-weight_HH32)==0)

# These are the complete source-coordinate first and second jets of Q.
# The denominator radius32 is positive on the stated nonnode cell.
source_vars32 = (xr32,xi32,yr32,yi32)
Qjet_ok32 = True
for i32,vi32 in enumerate(source_vars32):
    zi32 = spinor32.diff(vi32)
    Ni32 = zi32*spinor32.conjugate().T+spinor32*zi32.conjugate().T
    ri32 = sp.diff(radius32,vi32)
    Qi32 = Ni32/radius32-norm32*ri32/radius32**2
    Qjet_ok32 &= all(sp.cancel(v)==0 for v in Qi32-Q32.diff(vi32))
    for j32,vj32 in enumerate(source_vars32[:i32+1]):
        zj32 = spinor32.diff(vj32)
        Nj32 = zj32*spinor32.conjugate().T+spinor32*zj32.conjugate().T
        rj32 = sp.diff(radius32,vj32)
        Nij32 = (zi32*zj32.conjugate().T+zj32*zi32.conjugate().T)
        rij32 = sp.diff(radius32,vi32,vj32)
        Qij32 = (Nij32/radius32
                 -(Ni32*rj32+Nj32*ri32+norm32*rij32)/radius32**2
                 +2*norm32*ri32*rj32/radius32**3)
        Qjet_ok32 &= all(sp.cancel(v)==0
                         for v in Qij32-Q32.diff(vi32,vj32))
rep('C32c actual normalized cusp frame: all 4 first and 10 second jets',
    Qjet_ok32 and weight_jet32)

# Directly differentiate every entry of the 4x4 projector in its five local
# source variables, with the actual open flat weight substituted.  Check all
# 400 second entries for mixed-partial symmetry; no opaque smooth projector
# or free jet is introduced.
actual_P32 = P32.subs(weight32,actual_weight32)
full_vars32 = (*source_vars32,H32)
full_first32 = [actual_P32.diff(v) for v in full_vars32]
full_second32 = [[actual_P32.diff(vi,vj) for vj in full_vars32]
                 for vi in full_vars32]
full_jets32 = (len(full_first32)==5 and len(full_second32)==5
    and all(len(row)==5 for row in full_second32)
    and all(full_second32[i][j]==full_second32[j][i]
            for i in range(5) for j in range(5))
    and any(full_first32[4][0,j]!=0 for j in range(4))
    and any(full_second32[0][4][0,j]!=0 for j in range(4)))
rep('C32d nonconstant two-pair projector: 80 first and 400 second entries, symmetric',
    full_jets32)

# C33: the negative-height cusp cell has a genuinely different second frame.
# On -1/4<H<-1/10, |xy|<3/10 and y!=0, both pair weights are nonzero,
# while the source's H<0 deck branch replaces (x,y) in the second frame by
# (1/(m*y),m*x*y*y), m=nu*exp(-2|xy|^2/(1-|xy|^2)).  This is a local
# meridian-projector generator for either cusp, not a claim about the full
# meridian path.  The ordered source branch check below binds all assignments
# and guards together, rather than finding their AST nodes independently.
expected_cusp_body33 = ast.parse('''
x,y = map(mp.mpc,data[:2])
nu = mp.exp(-2) if kind=='ci' else mp.mpf(1)
H = mp.mpf(data[2]) if elementary else elliptic_height_chart(x,y,nu,dps=dps)
if abs(x*y)>mp.mpf('.3') or abs(H)>mp.mpf('.6'):
    raise ValueError('pure cusp embedding requires |xy|<=.3 and |H|<=.6')
f = frame(x,y,nu)
block = [hw(H)*q for q in f]
w = hw(H+mp.mpf('.5') if H<0 else H-mp.mpf('.5'))
if H<0 and w>0:
    rr = abs(x*y)**2
    m = nu*mp.exp(-2*rr/(1-rr))
    f = frame(1/(m*y),m*x*y*y,nu)
block.extend(w*q for q in f)
start = 4 if kind=='ci' else 0
vec[start:start+4] = block
''').body
cusp_test33 = ast.dump(ast.parse("kind in ('c0','ci')",mode='eval').body,
                       include_attributes=False)
cusp_nodes33 = [node for node in ast.walk(bundle_tree32)
                if isinstance(node,ast.If)
                and ast.dump(node.test,include_attributes=False)==cusp_test33]
source_branch33 = (len(cusp_nodes33)==1 and
    [ast.dump(node,include_attributes=False) for node in cusp_nodes33[0].body]
    ==[ast.dump(node,include_attributes=False) for node in expected_cusp_body33])
flat_tree33 = source_tree30(cusp30.scalar_flat_step)
expected_flat_body33 = ast.parse('''
import mpmath as mp
u = mp.mpf(u)
if u<=0:
    return mp.mpf(0)
if u>=1:
    return mp.mpf(1)
a,b = mp.exp(-1/u),mp.exp(-1/(1-u))
return a/(a+b)
''').body
source_branch33 &= (
    source_cusp32
    and [ast.dump(node,include_attributes=False) for node in flat_tree33.body[1:]]
        ==[ast.dump(node,include_attributes=False) for node in expected_flat_body33]
    and matching_if_test32(flat_tree33,'u<=0')
    and matching_if_test32(flat_tree33,'u>=1')
    and matching_statement31(flat_tree33,
                             'a,b = mp.exp(-1/u),mp.exp(-1/(1-u))')
    and matching_statement31(flat_tree33,'return a/(a+b)')
    and matching_if_test32(reference_tree32,'a>0')
    and matching_statement31(reference_tree32,
        "Y = 2/(R+f) if cell==-1 else (R-f)/2")
    and matching_statement31(reference_tree32,'y = mp.sqrt(Y)')
    and matching_statement31(reference_tree32,'return t/y,y')
    and matching_assignment30(raw_tree31,'phase_coordinate',
                              '-mp.conj(u1) if inf else u2')
    and matching_if_test32(raw_tree31,
        "phase_coordinate.real<-mp.mpf('.9')*abs(phase_coordinate)")
    and matching_assignment30(raw_tree31,'a',
                              'mp.arg(phase_coordinate)/(2*mp.pi)')
    and matching_statement31(raw_tree31,
        'xs,yy = elementary_cusp_reference(tc,a,nu)')
    and returned_expr30(raw_tree31,
        "('ci' if inf else 'c0'),(xs,yy,a)"))

H33 = sp.Symbol('H33',real=True)
cut33 = (H33+sp.Rational(1,4))/sp.Rational(3,20)
flat33 = sp.exp(-1/cut33)/(sp.exp(-1/cut33)
                             +sp.exp(-1/(1-cut33)))
weight33 = 1-flat33
# At t=1/10, a=-1/6 the reference lies in the strict cell.  Its cell
# integer is zero and the certified odd, increasing moment has f<0.  Write
# F=-f>0; then R=sqrt(F^2+4q), Y=(R+F)/2>0 and x=t/sqrt(Y), y=sqrt(Y).
# This also proves that the source's 1/y deck shift is defined.
t_w33,a_w33 = sp.Rational(1,10),-sp.Rational(1,6)
F_w33,nu_w33 = sp.symbols('F_w33 nu_w33',positive=True)
e_w33 = sp.exp(-sp.Rational(4,99))
q_w33 = e_w33*nu_w33**2*t_w33**2
R_w33 = sp.sqrt(F_w33**2+4*q_w33)
A_w33,Y_w33 = (R_w33-F_w33)/2,(R_w33+F_w33)/2
y_w33 = sp.sqrt(Y_w33)
x_w33 = t_w33/y_w33
witness33 = (
    sp.Rational(-1,4)<a_w33<sp.Rational(-1,10)
    and 0<cut33.subs(H33,a_w33)<1
    and sp.floor(a_w33+sp.Rational(1,2))==0
    and t_w33<sp.Rational(3,10)
    and bool(e_w33<1)
    and all(bool(q_w33.subs(nu_w33,nu)<sp.Rational(1,8))
            for nu in (sp.Integer(1),sp.exp(-2)))
    and sp.simplify(x_w33*y_w33-t_w33)==0
    and sp.simplify(A_w33*Y_w33-q_w33)==0
    and sp.simplify(e_w33*nu_w33**2*x_w33**2-A_w33)==0)
# An exact round-meridian witness reaches this cell, too: X=0 and
# n=(sqrt(3)/4,3sqrt(3)/8,5/8).  The raw stereographic algebra gives
# |u1|^2=13/49<36/49=|u2|^2, Arg(u2)=-pi/3, and
# |ratio|^2=13/36.  In the source's inverse_ratio=False branch,
# |small|^2=2z/(b+sqrt(b^2+4Bz))<z/b=13/(23A)<1/100,
# since A=(4pi)^2>144.  Thus it selects the south cusp, its direct
# argument chart, and a nonzero |tc|<1/10, with height a=-1/6.
n_w33 = (sp.sqrt(3)/4,3*sp.sqrt(3)/8,sp.Rational(5,8))
target_y_w33 = (n_w33[0],1-n_w33[2],-n_w33[1])
s2_w33 = sum(v*v for v in target_y_w33)
u1_w33 = (1-s2_w33+2*sp.I*target_y_w33[0])/(1+s2_w33)
u2_w33 = 2*(target_y_w33[1]+sp.I*target_y_w33[2])/(1+s2_w33)
ratio_w33 = sp.simplify(u1_w33/u2_w33)
z_w33 = sp.Rational(13,36)
Aprof_w33,Bprof_w33 = (4*sp.pi)**2,sp.exp(4)
bprof_w33 = Aprof_w33*(1-z_w33)
rootprof_w33 = sp.sqrt(bprof_w33**2+4*Bprof_w33*z_w33)
round_witness33 = (
    sp.simplify(sum(v*v for v in n_w33)-1)==0
    and sp.simplify(s2_w33-sp.Rational(3,4))==0
    and sp.simplify(u1_w33*sp.conjugate(u1_w33)-sp.Rational(13,49))==0
    and sp.simplify(u2_w33*sp.conjugate(u2_w33)-sp.Rational(36,49))==0
    and sp.arg(u2_w33)/(2*sp.pi)==a_w33
    and sp.re(u2_w33)>0
    and sp.simplify(ratio_w33*sp.conjugate(ratio_w33)-z_w33)==0
    and sp.simplify(bprof_w33-Aprof_w33*sp.Rational(23,36))==0
    and bool(sp.pi>3)
    and sp.Rational(13,23*144)<sp.Rational(1,100))
witness33 &= round_witness33
rep('C33a negative cusp deck branch: ordered source path and strict witness',
    source_branch33 and witness33)

xr33,xi33,yr33,yi33 = sp.symbols('xr33 xi33 yr33 yi33',real=True)
nu33,m33 = sp.symbols('nu33 m33',positive=True)
w33 = sp.Symbol('w33',real=True)
X33,Y33 = xr33+sp.I*xi33,yr33+sp.I*yi33
sx33,sy33 = xr33**2+xi33**2,yr33**2+yi33**2
rho33 = sx33*sy33
z033 = sp.Matrix([nu33*sp.conjugate(X33),Y33])
z133 = sp.Matrix([nu33*Y33/(m33*sy33),m33*X33*Y33**2])
r033 = nu33**2*sx33+sy33
r133 = nu33**2/(m33**2*sy33)+m33**2*sx33*sy33**2
R033,R133 = sp.symbols('R033 R133',positive=True)
V33 = sp.Matrix([z033[0]/sp.sqrt(R033),z033[1]/sp.sqrt(R033),
                 w33*z133[0]/sp.sqrt(R133),
                 w33*z133[1]/sp.sqrt(R133)])/sp.sqrt(1+w33**2)
P33 = V33*V33.conjugate().T
normal33 = ((r033/R033+w33**2*r133/R133)/(1+w33**2))
# The second spinor is precisely frame(1/(m*y),m*x*y^2,nu) before its
# radius normalization.  It has the same base product xy as the first.
shifted_value33 = (
    sp.cancel((1/(m33*Y33))*(m33*X33*Y33**2)-X33*Y33)==0
    and sp.cancel((z033.conjugate().T*z033)[0]-r033)==0
    and sp.cancel((z133.conjugate().T*z133)[0]-r133)==0
    and sp.cancel((V33.conjugate().T*V33)[0]-normal33)==0
    and sp.simplify(normal33.subs({R033:r033,R133:r133})-1)==0)
# An exact rank-one polynomial identity plus the checked unit norm gives
# P^2=P.  Positivity is manifest: sy>0, m>0, nu>0, and r0,r1>0.
v33 = sp.Matrix(sp.symbols('v0_33 v1_33 v2_33 v3_33'))
b33 = sp.Matrix(sp.symbols('b0_33 b1_33 b2_33 b3_33'))
abstract_P33 = v33*b33.T
rank_one33 = all(sp.expand(z)==0 for z in
    abstract_P33*abstract_P33-(b33.T*v33)[0]*abstract_P33)
rep('C33b actual shifted four-slot cusp projector has unit rank-one value',
    shifted_value33 and rank_one33)

actual_m33 = nu33*sp.exp(-2*rho33/(1-rho33))
source_vars33 = (xr33,xi33,yr33,yi33)
exp_prime33 = -2/(1-rho33)**2
exp_second33 = -4/(1-rho33)**3
mjet33 = True
for i33,vi33 in enumerate(source_vars33):
    rhoi33 = sp.diff(rho33,vi33)
    mi33 = actual_m33*exp_prime33*rhoi33
    mjet33 &= sp.simplify((sp.diff(actual_m33,vi33)-mi33)/actual_m33)==0
    for vj33 in source_vars33[:i33+1]:
        rhoj33 = sp.diff(rho33,vj33)
        mij33 = actual_m33*((exp_second33+exp_prime33**2)
                             *rhoi33*rhoj33
                             +exp_prime33*sp.diff(rho33,vi33,vj33))
        mjet33 &= sp.simplify((sp.diff(actual_m33,vi33,vj33)
                               -mij33)/actual_m33)==0
slope33 = flat33*(1-flat33)*(1/cut33**2+1/(1-cut33)**2)
curvature33 = flat33*(1-flat33)*(
    (1-2*flat33)*(1/cut33**2+1/(1-cut33)**2)**2
    -2/cut33**3+2/(1-cut33)**3)
wjet33 = (sp.simplify(sp.diff(weight33,H33)
                     +sp.Rational(20,3)*slope33)==0
          and sp.simplify(sp.diff(weight33,H33,2)
                          +sp.Rational(400,9)*curvature33)==0)
rep('C33c source deck exponential and negative flat weight: exact two-jets',
    mjet33 and wjet33)

# Substitution after conjugation retains the certified positive-radius
# meaning of both square roots.  Differentiating these actual components
# supplies all 20 first and 100 second entries of V, and the explicit
# rank-one product supplies all 80 first and 400 second entries of P.
actual_V33 = V33.subs({R033:r033,R133:r133,m33:actual_m33,
                        w33:weight33})
actual_barV33 = V33.conjugate().T.subs(
    {R033:r033,R133:r133,m33:actual_m33,w33:weight33})
full_vars33 = (*source_vars33,H33)
V_first33 = [actual_V33.diff(v) for v in full_vars33]
barV_first33 = [actual_barV33.diff(v) for v in full_vars33]
V_second33 = [[actual_V33.diff(vi,vj) for vj in full_vars33]
              for vi in full_vars33]
barV_second33 = [[actual_barV33.diff(vi,vj) for vj in full_vars33]
                 for vi in full_vars33]
P_first33 = [V_first33[i]*actual_barV33
             +actual_V33*barV_first33[i] for i in range(5)]
P_second33 = [[V_second33[i][j]*actual_barV33
               +V_first33[i]*barV_first33[j]
               +V_first33[j]*barV_first33[i]
               +actual_V33*barV_second33[i][j]
               for j in range(5)] for i in range(5)]
full_jets33 = (
    all(V_second33[i][j]==V_second33[j][i]
        and barV_second33[i][j]==barV_second33[j][i]
        and P_second33[i][j]==P_second33[j][i]
        for i in range(5) for j in range(5))
    and any(V_first33[4][k]!=0 for k in range(4))
    and any(V_first33[0][k]!=0 for k in range(4))
    and any(P_first33[4][0,k]!=0 for k in range(4))
    and any(P_second33[0][4][0,k]!=0 for k in range(4)))
rep('C33d shifted source projector: complete five-variable first and second jets',
    full_jets33 and generic_projector31)

# C34: one active-collar principal-log seam, at negative real t.  The raw
# source uses a real Cramer solve before discarding its second and third
# lattice coordinates.  Both cusp and middle slots are present for
# 3/10<|t_c|<1/2.  We check the source's south/north log-deck changes, then
# the surviving angular-slot phase and cusp reconstruction on the seam.
lattice_tree34 = source_tree30(cusp30.scalar_lattice_coords)
logradius_tree34 = source_tree30(cusp30.elementary_cusp_logradius)
source_collar34 = (
    matching_assignment30(lattice_tree34,'rows',
        '[[v.imag for v in ez],[v.real for v in ew],[v.imag for v in ew]]')
    and matching_assignment30(lattice_tree34,'rhs',
                              '(z.imag,w.real,w.imag)')
    and matching_statement31(lattice_tree34,
        'replaced = [[rhs[i] if j==k else rows[i][j] for j in range(3)] for i in range(3)]')
    and matching_statement31(lattice_tree34,'coords.append(det(replaced)/den)')
    and matching_assignment30(lattice_tree34,'one',
        '(z.real-sum(a*b.real for a,b in zip(coords,ez)))/(2*mp.pi)')
    and returned_expr30(lattice_tree34,'(one,*coords)')
    and matching_statement31(raw_tree31,'''if inf:
    zm,wm = -mp.j*rho,mp.mpc(0)
else:
    ss = lt/(2*mp.j*mp.pi)
    zm,wm = -mp.j*rho+mp.pi*ss*(ss-1),2*t*ss''')
    and matching_statement31(raw_tree31,
        'lattice = scalar_lattice_coords(zm,wm,t,lt)')
    and matching_assignment30(raw_tree31,'rho',
        '(elementary_cusp_logradius if elementary else elliptic_cusp_logradius)(a,tc,nu,dps=dps)')
    and matching_assignment30(logradius_tree34,'r','abs(t)')
    and matching_if_test32(raw_tree31,"inf or abs(t)<mp.mpf('.5')")
    and matching_if_test32(raw_tree31,
                           "abs(tc)<=mp.mpf('.3') and elementary")
    and matching_assignment30(raw_tree31,'data',
        '(t,lattice[0],lattice[3],lt,a) if elementary else (t,lattice[0],lattice[3],lt)')
    and returned_expr30(raw_tree31,"'mid',data")
    and matching_assignment30(bundle_tree32,'principal','mp.log(t)')
    and matching_assignment30(bundle_tree32,'branch',
        'mp.nint((lt.imag-principal.imag)/(2*mp.pi))')
    and matching_statement31(bundle_tree32,
        'one,lt = one+branch*four,principal')
    and matching_assignment30(bundle_tree32,'known_height',
        'mp.mpf(data[4]) if elementary and len(data)>4 else None')
    and matching_assignment30(bundle_tree32,'b0',
        "1-scalar_flat_step((r-mp.mpf('.3'))/mp.mpf('.2'))")
    and matching_assignment30(bundle_tree32,'bi',
        "1-scalar_flat_step((1/r-mp.mpf('.3'))/mp.mpf('.2'))")
    and matching_assignment30(bundle_tree32,'bm','1-b0-bi')
    and matching_if_test32(bundle_tree32,'b0>0')
    and matching_if_test32(bundle_tree32,'bi>0')
    and matching_if_test32(bundle_tree32,'bm>0')
    and matching_statement31(bundle_tree32,
        'z,w = 2*mp.pi*one-mp.j*lt*four,2*t*four')
    and matching_statement31(bundle_tree32,
        'det = e2.real*e3.imag-e3.real*e2.imag')
    and matching_statement31(bundle_tree32,
        'return ((rhs.real*e3.imag-e3.real*rhs.imag)/det,'
        '(e2.real*rhs.imag-rhs.real*e2.imag)/det)')
    and matching_statement31(bundle_tree32,
        'zc,wc = z-mp.pi*s*(s-1),w-2*t*s')
    and matching_statement31(bundle_tree32,
        's = lt/(2*mp.j*mp.pi)')
    and matching_statement31(bundle_tree32,
        'c2,c3 = coeffs(-wc,1-t*t,mp.j*(1+t*t))')
    and matching_statement31(bundle_tree32,
        'zz = mp.exp(mp.j*(zc-mp.j*c2*t-c3*t))')
    and matching_statement31(bundle_tree32,
        'vec[:4] = [b0*q for q in torus_pair(zz,t,mp.mpf(1),known_height)]')
    and matching_statement31(bundle_tree32,
        'zc,wc = z-mp.j*w/t,w/(t*t)')
    and matching_statement31(bundle_tree32,'s = 1/t')
    and matching_statement31(bundle_tree32,
        'c2,c3 = coeffs(-wc,s*s-1,mp.j*(1+s*s))')
    and matching_statement31(bundle_tree32,
        'zz = mp.exp(mp.j*(zc-mp.j*c2*s+c3*s))')
    and matching_statement31(bundle_tree32,
        'vec[4:8] = [bi*q for q in torus_pair(zz,s,mp.exp(-2),known_height)]')
    and matching_assignment30(bundle_tree32,'angle',
        'mp.arg(t if index==0 else -t)')
    and matching_assignment30(bundle_tree32,'weight',
        '1-scalar_flat_step((abs(angle)-mp.pi/2)/(mp.pi/4))')
    and matching_statement31(bundle_tree32,
        'vec[8+index]=bm*weight*mp.exp(mp.j*(2*mp.pi*one'
        '-(angle+centre-lt.imag)*four))'))

# Two exact target witnesses, both with n=(0,0,-1).  For the south take
# X=sqrt(2/37); for the north take X=sqrt(2/5).  The raw y_1 values are
# 12/5 and 4, so (u1,u2)=(-119,120)/169 and (-15,8)/17.  Both phase
# coordinates are positive real, hence a=0 and no argument-selector seam.
Xs34,Xn34 = sp.sqrt(sp.Rational(2,37)),sp.sqrt(sp.Rational(2,5))
def target_y1_34(X):
    aa = X*sp.sqrt(2-X**2)
    return 1+(1-X**2)/(1-aa)
ys34,yn34 = target_y1_34(Xs34),target_y1_34(Xn34)
def unit_pair34(y):
    return ((1-y**2)/(1+y**2),2*y/(1+y**2))
us34,un34 = unit_pair34(ys34),unit_pair34(yn34)
zs34,zn34 = (sp.Rational(119,120)**2,
             sp.Rational(8,15)**2)
# Let A=(4pi)^2, B=e^4.  The exponential series has S4=65/24;
# for n>=5, n! >= 5!*6^(n-5), with a strict inequality from n=7,
# so 0<e-S4<1/(5!*(1-1/6))=1/100.  Thus 50<B<64.
# The classical 3<pi<22/7 gives
# 144<A<160.  For the south s=|small|^2 solves
# B*s^2+A(1-zs)*s-zs=0; for the north s solves
# s^2+A(1-zn)*s-B*zn=0.  The following strict rational signs bracket
# both positive roots in (9/100,1/4).  Thus both have 3/10<|tc|<1/2.
lo34,hi34 = sp.Rational(9,100),sp.Rational(1,4)
S4_34 = sum(1/sp.factorial(j) for j in range(5))
etail34 = 1/(sp.factorial(5)*(1-sp.Rational(1,6)))
strict_collar34 = (
    sp.simplify(ys34-sp.Rational(12,5))==0
    and sp.simplify(yn34-4)==0
    and us34==(sp.Rational(-119,169),sp.Rational(120,169))
    and un34==(sp.Rational(-15,17),sp.Rational(8,17))
    and S4_34==sp.Rational(65,24)
    and etail34==sp.Rational(1,100)
    and S4_34>sp.Rational(8,3)
    and S4_34+etail34<sp.Rational(11,4)
    and sp.Rational(8,3)**4>50
    and sp.Rational(11,4)**4<64
    and bool(sp.pi>3) and bool(sp.pi<sp.Rational(22,7))
    and 16*sp.Rational(22,7)**2<160
    and 64*lo34**2+160*(1-zs34)*lo34-zs34<0
    and 50*hi34**2-zs34>0
    and lo34**2+160*(1-zn34)*lo34-50*zn34<0
    and hi34**2+144*(1-zn34)*hi34-64*zn34>0)
# Perturb n=(0,0,-1) to (eps,0,-sqrt(1-eps^2)).  With the same X,
# y0=R*eps, y1=1+R*sqrt(1-eps^2), and u1/u2 has imaginary part y0/y1.
# Thus eps of either sign gives the upper/lower negative-axis sector,
# while u2 stays positive.  All collar inequalities above are strict,
# so they persist on an open interval around eps=0.
eps34 = sp.Symbol('eps34',real=True)
positive_root34 = sp.Symbol('positive_root34',positive=True)
def angular_perturb34(radial_factor):
    y0 = radial_factor*eps34
    # On |eps|<1 substitute positive_root=sqrt(1-eps^2) after the
    # real-algebra check; keeping the positive root formal avoids SymPy's
    # global complex continuation of sqrt outside this open interval.
    y1 = 1+radial_factor*positive_root34
    ysq = y0**2+y1**2
    u1 = (1-ysq+2*sp.I*y0)/(1+ysq)
    u2 = 2*y1/(1+ysq)
    return sp.simplify(sp.im(u1/u2)-y0/y1)==0
strict_collar34 &= (angular_perturb34(sp.Rational(7,5))
                    and angular_perturb34(sp.Integer(3)))
rep('C34a source-bound south/north active-collar negative-axis witnesses',
    source_collar34 and strict_collar34)

ar34,ai34,ell34,theta34,rho34 = sp.symbols(
    'ar34 ai34 ell34 theta34 rho34',real=True)
k34 = sp.Symbol('k34',integer=True)
one34,two34,three34,four34 = sp.symbols(
    'one34 two34 three34 four34',real=True)
t34 = ar34+sp.I*ai34
L34 = ell34+sp.I*theta34
s34 = L34/(2*sp.I*sp.pi)
e2z34,e3z34,e4z34 = -sp.I*t34,-t34,-sp.I*L34
e2w34,e3w34,e4w34 = 1-t34**2,sp.I*(1+t34**2),2*t34
R34 = sp.Matrix([[-ar34,-ai34,-ell34],
                 [1-ar34**2+ai34**2,-2*ar34*ai34,2*ar34],
                 [-2*ar34*ai34,1+ar34**2-ai34**2,2*ai34]])
det_expected34 = (1+ar34**2+ai34**2)*(
    (ar34**2+ai34**2-1)*ell34+2*(ar34**2+ai34**2))
cramer34 = (R34==M.subs({sp.Symbol('t_r',real=True):ar34,
                         sp.Symbol('t_i',real=True):ai34,lr:ell34})
            and sp.factor(R34.det()-det_expected34)==0)
# In the south collar 3/10<|t|<1/2 and ell=log|t|<0; in the north
# 2<|t|<10/3 and ell>0.  Thus (|t|^2-1)*ell>=0 in either collar,
# and det R>2r^2(1+r^2).  This polynomial increases for r>0 because
# its derivative is 4r(1+2r^2)>0.  The strict margins are 981/5000
# south and 40 north, so the Cramer coordinates are unique.
radius_bound34 = sp.Symbol('radius_bound34',positive=True)
det_margin34 = 2*radius_bound34**2*(1+radius_bound34**2)
cramer34 &= (
    sp.simplify(sp.diff(det_margin34,radius_bound34)
                -4*radius_bound34*(1+2*radius_bound34**2))==0
    and det_margin34.subs(radius_bound34,sp.Rational(3,10))
        ==sp.Rational(981,5000)
    and det_margin34.subs(radius_bound34,2)==40)
rep('C34b active-collar real-Cramer matrix and positive determinant',
    cramer34)
# Numerical observations (d), two interior radii, kept beside the
# proven uniform lower bounds.  They play no role in the certificate.
south_sample34 = sp.N(det_expected34.subs({
    ar34:sp.Rational(2,5),ai34:0,
    ell34:sp.log(sp.Rational(2,5))}),20)
north_sample34 = sp.N(det_expected34.subs({
    ar34:sp.Rational(5,2),ai34:0,
    ell34:sp.log(sp.Rational(5,2))}),20)
rep('C34b2 two numerical Cramer determinant samples (d)',
    south_sample34>sp.Rational(981,5000)
    and north_sample34>40)
print('C34b2 south/north determinant samples:', south_sample34,
      north_sample34)

def lattice_z34(L,one,two,three,four):
    return 2*sp.pi*one-sp.I*t34*two-t34*three-sp.I*L*four
def lattice_w34(two,three,four):
    return e2w34*two+e3w34*three+2*t34*four
Lshift34 = L34+2*sp.pi*sp.I*k34
south_one34 = one34-k34*four34-k34*(k34+1)/2
south_four34 = four34+k34
north_one34 = one34-k34*four34
south_deck34 = (
    sp.expand(lattice_z34(Lshift34,south_one34,two34,three34,
                          south_four34)
              -lattice_z34(L34,one34,two34,three34,four34)
              -sp.pi*((s34+k34)*(s34+k34-1)-s34*(s34-1)))==0
    and sp.expand(lattice_w34(two34,three34,south_four34)
                  -lattice_w34(two34,three34,four34)-2*t34*k34)==0)
north_deck34 = (
    sp.expand(lattice_z34(Lshift34,north_one34,two34,three34,
                          four34)-lattice_z34(L34,one34,two34,
                                               three34,four34))==0
    and lattice_w34(two34,three34,four34)
        ==lattice_w34(two34,three34,four34))
rep('C34c raw south/north lattice one/four transformations for every deck k',
    south_deck34 and north_deck34)

# The source's middle reduction first takes the principal logarithm.  For
# a raw input whose log differs by 2pi*i*k, its south normalized one changes
# by k(k-1)/2, an integer; its north normalized one is unchanged.
normalized_south34 = south_one34+k34*south_four34
normalized_north34 = north_one34+k34*four34
principal34 = (
    sp.simplify(normalized_south34-one34-k34*(k34-1)/2)==0
    and sp.simplify(normalized_north34-one34)==0)
# On the upper side Arg(t)=theta near pi, Arg(-t)=theta-pi, so
# delta=Arg(-t)+pi-Im(L)=0.  On the lower side the corresponding delta is
# 2pi.  The index-0 angular weight vanishes (|Arg t|>3pi/4), while the
# index-1 weight is 1 (|Arg(-t)|<pi/2).  For k=-1, the lower exponent
# differs by 2pi south and by zero north, so its complex phase agrees.
delta_upper34 = (theta34-sp.pi)+sp.pi-theta34
delta_lower34 = (theta34+sp.pi)+sp.pi-theta34
angular34 = (
    delta_upper34==0 and delta_lower34==2*sp.pi
    and sp.Rational(7,8)>sp.Rational(3,4)
    and sp.simplify((south_one34-south_four34).subs(k34,-1)
                    -one34-1)==0
    and sp.simplify((north_one34-four34).subs(k34,-1)
                    -one34)==0)
rep('C34d surviving angular-slot phase agrees across the negative-axis cut',
    principal34 and angular34)

# Cramer's equations reconstruct the full raw (z_m,w_m).  In the south,
# z_m=-i*rho+pi*s(s-1), w_m=2*t*s.  In the north, z_m=-i*rho,w_m=0.
# The source's 2x2 real coeffs eliminate the discarded lattice coordinates:
# its determinants are 1-|t|^4>0 (south) and |1/t|^4-1<0 (north),
# both nonzero on the active collars.  The exponent in either `zz` is rho.
raw_z_s34 = -sp.I*rho34+sp.pi*s34*(s34-1)
raw_z_n34 = -sp.I*rho34
def real_det2_34(e2,e3):
    return sp.re(e2)*sp.im(e3)-sp.re(e3)*sp.im(e2)
def source_coeff34(rhs,e2,e3):
    den = real_det2_34(e2,e3)
    return ((sp.re(rhs)*sp.im(e3)-sp.re(e3)*sp.im(rhs))/den,
            (sp.re(e2)*sp.im(rhs)-sp.re(rhs)*sp.im(e2))/den)
def source_lattice34(raw_z,raw_w):
    # The source computes these three entries by Cramer's rule, which equals
    # this exact solve because C34b proves det(R)!=0.
    rhs = sp.Matrix([sp.im(sp.expand_complex(raw_z)),
                     sp.re(sp.expand_complex(raw_w)),
                     sp.im(sp.expand_complex(raw_w))])
    two,three,four = R34.LUsolve(rhs)
    one = (sp.re(sp.expand_complex(raw_z))
           -two*sp.re(e2z34)-three*sp.re(e3z34)
           -four*sp.re(e4z34))/(2*sp.pi)
    z_middle = 2*sp.pi*one-sp.I*L34*four
    w_middle = 2*t34*four
    return two,three,z_middle,w_middle
two_s34,three_s34,z_mid_s34,w_mid_s34 = source_lattice34(
    raw_z_s34,2*t34*s34)
zc_s34,wc_s34 = z_mid_s34-sp.pi*s34*(s34-1),w_mid_s34-2*t34*s34
c2_s34,c3_s34 = source_coeff34(-wc_s34,e2w34,e3w34)
south_reconstruction34 = (
    sp.simplify(c2_s34-two_s34)==0
    and sp.simplify(c3_s34-three_s34)==0
    and sp.simplify(sp.I*(zc_s34-sp.I*c2_s34*t34-c3_s34*t34)
                    -rho34)==0)
sn34 = 1/t34
two_n34,three_n34,z_mid_n34,w_mid_n34 = source_lattice34(
    raw_z_n34,sp.Integer(0))
zc_n34,wc_n34 = z_mid_n34-sp.I*w_mid_n34/t34,w_mid_n34/t34**2
c2_n34,c3_n34 = source_coeff34(
    -wc_n34,sn34**2-1,sp.I*(1+sn34**2))
north_reconstruction34 = (
    sp.simplify(c2_n34-two_n34)==0
    and sp.simplify(c3_n34-three_n34)==0
    and sp.simplify(sp.I*(zc_n34-sp.I*c2_n34*sn34+c3_n34*sn34)
                    -rho34)==0)
south_det2_34 = sp.simplify(real_det2_34(e2w34,e3w34)
                            -(1-(ar34**2+ai34**2)**2))==0
sr34,si34 = sp.symbols('sr34 si34',real=True)
ss34 = sr34+sp.I*si34
north_det2_34 = sp.simplify(real_det2_34(
    ss34**2-1,sp.I*(1+ss34**2))
    -((sr34**2+si34**2)**2-1))==0
inverse_radius34 = sp.cancel(sn34*sp.conjugate(sn34)
                               -1/(ar34**2+ai34**2))==0
rep('C34e both collar cusp reconstructions give zz=exp(rho)',
    south_reconstruction34 and north_reconstruction34
    and south_det2_34 and north_det2_34
    and inverse_radius34)

# C34f: second jets of precisely the two scalar entries just glued.
# Use local real chart variables (Re t,Im t,rho,H,one,four).  The source's
# surviving index-1 entry has bm=s(q), with q=5|t|-3/2 south and
# q=5/|t|-3/2 north.  C34d makes its phase exp(2pi*i*one) on both sides.
# C34e makes the reconstructed cusp input zz=exp(rho) on both sides.
# These are chart-coordinate jets; rho, one and four still depend on the
# round meridian variables through earlier maps, which this seam check does
# not differentiate.
Hseam34 = sp.Symbol('Hseam34',real=True)
p34 = sp.Symbol('p34',real=True)
r34 = sp.sqrt(ar34**2+ai34**2)
q_s34,q_n34 = 5*r34-sp.Rational(3,2),5/r34-sp.Rational(3,2)
flat_actual34 = sp.exp(-1/p34)/(
    sp.exp(-1/p34)+sp.exp(-1/(1-p34)))
flat_slope34 = flat_actual34*(1-flat_actual34)*(
    1/p34**2+1/(1-p34)**2)
flat_curv34 = flat_actual34*(1-flat_actual34)*(
    (1-2*flat_actual34)*(1/p34**2+1/(1-p34)**2)**2
    -2/p34**3+2/(1-p34)**3)
flat_jets34 = (
    sp.simplify(sp.diff(flat_actual34,p34)-flat_slope34)==0
    and sp.simplify(sp.diff(flat_actual34,p34,2)-flat_curv34)==0)
radial_jets34 = True
radial_xy34 = (ar34,ai34)
for i34,vi34 in enumerate(radial_xy34):
    radial_jets34 &= (
        sp.simplify(sp.diff(q_s34,vi34)-5*vi34/r34)==0
        and sp.simplify(sp.diff(q_n34,vi34)+5*vi34/r34**3)==0)
    for j34,vj34 in enumerate(radial_xy34[:i34+1]):
        same34 = int(i34==j34)
        radial_jets34 &= (
            sp.simplify(sp.diff(q_s34,vi34,vj34)
                        -5*(same34/r34-vi34*vj34/r34**3))==0
            and sp.simplify(sp.diff(q_n34,vi34,vj34)
                            -(-5*same34/r34**3
                              +15*vi34*vj34/r34**5))==0)

F34 = sp.Function('F34')
phase_common34 = sp.exp(2*sp.pi*sp.I*one34)
vars_seam34 = (ar34,ai34,rho34,Hseam34,one34,four34)
angular_value_jets34 = True
source_flat_bind34 = True
for qcell34,lower_one34,lower_four34 in (
    (q_s34,south_one34.subs(k34,-1),south_four34.subs(k34,-1)),
    (q_n34,north_one34.subs(k34,-1),four34)):
    f0cell34 = F34(qcell34)
    f1cell34 = sp.Subs(sp.diff(F34(p34),p34),p34,qcell34)
    f2cell34 = sp.Subs(sp.diff(F34(p34),p34,2),p34,qcell34)
    # On this strict 0<q<1 cell the source flat-step body is exactly
    # flat_actual34.  Bind every formal value/jet used below to its
    # source-specific exponential expression, whose derivatives were
    # checked in flat_jets34 above.
    flat_bind34 = {
        f0cell34:flat_actual34.subs(p34,qcell34),
        f1cell34:flat_slope34.subs(p34,qcell34),
        f2cell34:flat_curv34.subs(p34,qcell34)}
    source_flat_bind34 &= (
        not any(expr.has(F34) for expr in flat_bind34.values())
        and flat_bind34[f0cell34]
            ==sp.Lambda(p34,flat_actual34)(qcell34))
    common34 = f0cell34*phase_common34
    upper34 = f0cell34*sp.exp(sp.I*(2*sp.pi*one34
                                        -delta_upper34*four34))
    lower34 = f0cell34*sp.exp(sp.I*(2*sp.pi*lower_one34
                                        -delta_lower34*lower_four34))
    angular_value_jets34 &= sp.simplify(upper34-lower34)==0
    for i34,vi34 in enumerate(vars_seam34):
        qi34 = sp.diff(qcell34,vi34)
        di34 = int(vi34==one34)
        expected_first34 = phase_common34*(
            f1cell34*qi34+2*sp.pi*sp.I*f0cell34*di34)
        source_flat_bind34 &= not expected_first34.xreplace(
            flat_bind34).has(F34)
        angular_value_jets34 &= (
            sp.simplify(sp.diff(common34,vi34)-expected_first34)==0
            and sp.simplify(sp.diff(upper34,vi34)
                            -sp.diff(lower34,vi34))==0)
        for j34,vj34 in enumerate(vars_seam34[:i34+1]):
            qj34 = sp.diff(qcell34,vj34)
            qij34 = sp.diff(qcell34,vi34,vj34)
            dj34 = int(vj34==one34)
            expected_second34 = phase_common34*(
                f2cell34*qi34*qj34+f1cell34*qij34
                +2*sp.pi*sp.I*f1cell34*(qi34*dj34+qj34*di34)
                -(2*sp.pi)**2*f0cell34*di34*dj34)
            source_flat_bind34 &= not expected_second34.xreplace(
                flat_bind34).has(F34)
            angular_value_jets34 &= (
                sp.simplify(sp.diff(common34,vi34,vj34)
                            -expected_second34)==0
                and sp.simplify(sp.diff(upper34,vi34,vj34)
                                -sp.diff(lower34,vi34,vj34))==0)
zz_common34 = sp.exp(rho34)
zz_jets34 = all(
    sp.diff(zz_common34,vi34)
        ==zz_common34*int(vi34==rho34)
    and all(sp.diff(zz_common34,vi34,vj34)
            ==zz_common34*int(vi34==rho34)*int(vj34==rho34)
            for vj34 in vars_seam34)
    for vi34 in vars_seam34)
rep('C34f index-1 angular and cusp-zz first/second seam jets agree',
    flat_jets34 and radial_jets34 and source_flat_bind34
    and angular_value_jets34
    and zz_jets34 and south_reconstruction34 and north_reconstruction34)

# C35: the complete ten-slot projector on one active-collar upper sector.
# Use the open cell 3/10<|tc|<1/2, 3pi/4<Arg(t_middle)<pi and |H|<1/10.
# C34e reconstructs zz=exp(rho); the centre-0 torus pair has shift zero,
# weight one, and frame(tc/zz,zz,nu).  The centre-1 pair is zero.  The
# index-0 angular slot is zero, while index 1 is bm*exp(2pi*i*one).
# The south support is (0,1,9), the north support is (4,5,9).  This is one
# strict local source circuit, not a claim about the whole closed projector.
torus_nodes35 = [node for node in ast.walk(bundle_tree32)
                 if isinstance(node,ast.FunctionDef)
                 and node.name=='torus_pair']
frame_nodes35 = [node for node in ast.walk(bundle_tree32)
                 if isinstance(node,ast.FunctionDef) and node.name=='frame']
hw_nodes35 = [node for node in ast.walk(bundle_tree32)
              if isinstance(node,ast.FunctionDef) and node.name=='hw']
expected_hw_body35 = ast.parse('''
return 1-scalar_flat_step((abs(h)-mp.mpf('.25'))/mp.mpf('.15'))
''').body
expected_torus_body35 = ast.parse('''
rr = abs(t)**2
m = nu*mp.exp(-2*rr/(1-rr))
P = mp.log(m*abs(t))
H = known_height if elementary else elliptic_height(mp.log(abs(z)),P,dps=dps)
if H is None:
    raise ValueError('elementary projector requires the known inverse-first height')
out = []
for centre in (mp.mpf(0),mp.mpf('.5')):
    shift = mp.floor(centre+mp.mpf('.5')-H)
    w = hw(H+shift-centre)
    if w == 0:
        out.extend([mp.mpc(0),mp.mpc(0)])
    else:
        yy = z*(m*t)**int(shift)
        out.extend(w*q for q in frame(t/yy,yy,nu))
return out
''').body
expected_frame_body35 = ast.parse('''
r = nu*nu*abs(x)**2+abs(y)**2
if not r>0:
    raise ValueError('the bundle embedding excludes the fixed points')
return (nu*mp.conj(x)/mp.sqrt(r),y/mp.sqrt(r))
''').body
source_torus35 = (
    source_collar34 and source_branch33 and len(torus_nodes35)==1
    and len(frame_nodes35)==1 and len(hw_nodes35)==1
    and [ast.dump(node,include_attributes=False)
         for node in hw_nodes35[0].body]
        ==[ast.dump(node,include_attributes=False)
           for node in expected_hw_body35]
    and [ast.dump(node,include_attributes=False)
         for node in torus_nodes35[0].body]
        ==[ast.dump(node,include_attributes=False)
           for node in expected_torus_body35]
    and [ast.dump(node,include_attributes=False)
         for node in frame_nodes35[0].body]
        ==[ast.dump(node,include_attributes=False)
           for node in expected_frame_body35]
    and matching_statement31(bundle_tree32,
        "return 1-scalar_flat_step((abs(h)-mp.mpf('.25'))/mp.mpf('.15'))")
    and matching_assignment30(bundle_tree32,'norm2',
        'sum(abs(q)**2 for q in vec)')
    and matching_assignment30(bundle_tree32,'result',
        'tuple(+(q/mp.sqrt(norm2)) for q in vec)'))

# Exact off-cut round witnesses.  Take a stereographic y=(1/100,Y,0),
# with Y=12/5 south or Y=4 north.  Set R=|y-(0,1,0)|,
# X=(R-1)/sqrt(R^2+1), n=(y0/R,0,(1-Y)/R).  Then |n|=1,
# X*sqrt(2-X^2)=(R^2-1)/(R^2+1), and the source factor
# (1-X^2)/(1-X*sqrt(2-X^2)) equals R, recovering this y exactly.
y0_w35 = sp.Rational(1,100)
def round_witness35(Y):
    radius = sp.sqrt(y0_w35**2+(Y-1)**2)
    X = (radius-1)/sp.sqrt(radius**2+1)
    n = (y0_w35/radius,sp.Integer(0),(1-Y)/radius)
    aa = (radius**2-1)/(radius**2+1)
    radius_ok = (bool(radius>1)
        and sp.simplify(sum(z*z for z in n)-1)==0
        and sp.simplify(X*sp.sqrt(2-X**2)-aa)==0
        and sp.simplify((1-X**2)/(1-aa)-radius)==0)
    ysq = y0_w35**2+Y**2
    u1 = (1-ysq+2*sp.I*y0_w35)/(1+ysq)
    u2 = 2*Y/(1+ysq)
    return radius_ok,u1,u2
south_w35,u1s35,u2s35 = round_witness35(sp.Rational(12,5))
north_w35,u1n35,u2n35 = round_witness35(sp.Integer(4))
z_s35 = sp.simplify(u1s35*sp.conjugate(u1s35)/u2s35**2)
z_n35 = sp.simplify(u2n35**2/(u1n35*sp.conjugate(u1n35)))
lo35,hi35 = sp.Rational(9,100),sp.Rational(1,4)
# The same 144<A=(4pi)^2<160 and 50<B=e^4<64 bounds as C34a
# bracket the south positive root of B*s^2+A(1-z)*s-z and the north
# positive root of s^2+A(1-z)*s-B*z between lo=9/100 and hi=1/4.
root_cell35 = (
    south_w35 and north_w35 and 0<z_s35<1 and 0<z_n35<1
    and 64*lo35**2+160*(1-z_s35)*lo35-z_s35<0
    and 50*hi35**2-z_s35>0
    and lo35**2+160*(1-z_n35)*lo35-50*z_n35<0
    and hi35**2+144*(1-z_n35)*hi35-64*z_n35>0)
# Both u1/u2 ratios have negative real and positive imaginary parts.
# The ratio of imaginary to the magnitude of real part is
# 2*y0/(y0^2+Y^2-1)<1.  Thus the middle t lies strictly in the upper
# sector 3pi/4<Arg(t)<pi; neither argument cut is touched.  South H=0.
# North H=Arg(-conj(u1))/(2pi), and 0<H<1/10 follows from
# arctan(v)<v, v<1/100<pi/5.  The derivative of arctan is
# 1/(1+v^2)<1 for v>0, giving the stated elementary inequality.
angles35 = [2*y0_w35/(y0_w35**2+Y**2-1)
            for Y in (sp.Rational(12,5),sp.Integer(4))]
height_cell35 = (
    all(0<v<1 for v in angles35)
    and angles35[1]<sp.Rational(1,100)
    and sp.Rational(1,100)<sp.pi/5
    and sp.re(u2s35)>0 and sp.re(-sp.conjugate(u1n35))>0
    and sp.im(u1s35/u2s35)>0
    and sp.im(u1n35/u2n35)>0)
# For |H|<1/10 the centre-0 shift is zero and hw(H)=1; the
# centre-1 offset has magnitude >2/5, so hw=0.  The angular weights
# follow from 3pi/4<Arg(t)<pi: index 0 is zero, index 1 is one.
branch_margin35 = (
    sp.Rational(1,10)<sp.Rational(1,4)
    and (sp.Rational(2,5)-sp.Rational(1,4))/sp.Rational(3,20)==1
    and (sp.Rational(3,4)-sp.Rational(1,2))/sp.Rational(1,4)==1
    and sp.Rational(1,4)<sp.Rational(1,2))
rep('C35a source torus-pair circuit and exact off-cut south/north witnesses',
    source_torus35 and root_cell35 and height_cell35 and branch_margin35)

# In either cusp let tc=a+i*b, zz=e^rho, nu=1 (south) or e^-2 (north).
# The source frame is (nu*exp(-rho)*conj(tc),exp(rho))/sqrt(Rf).
# Put c=1-s(q), d=s(q), q=5|tc|-3/2.  The active source norm is
# D=c^2+d^2=1/2+2(d-1/2)^2>=1/2, with Rf>0 manifest.
ar35,ai35,rho35,one35,H35 = sp.symbols(
    'ar35 ai35 rho35 one35 H35',real=True)
nu35 = sp.Symbol('nu35',positive=True)
tc35 = ar35+sp.I*ai35
rsq35 = ar35**2+ai35**2
r35 = sp.sqrt(rsq35)
q35 = 5*r35-sp.Rational(3,2)
flat35 = sp.exp(-1/q35)/(
    sp.exp(-1/q35)+sp.exp(-1/(1-q35)))
c35,d35 = 1-flat35,flat35
Rf35 = nu35**2*sp.exp(-2*rho35)*rsq35+sp.exp(2*rho35)
D35 = c35**2+d35**2
phase35 = sp.exp(2*sp.pi*sp.I*one35)
V35 = sp.Matrix([c35*nu35*sp.exp(-rho35)*sp.conjugate(tc35)
                 /sp.sqrt(Rf35),
                 c35*sp.exp(rho35)/sp.sqrt(Rf35),
                 d35*phase35])/sp.sqrt(D35)
# Check normalization first with a real symbolic weight and positive
# formal frame radius/denominator, avoiding any hidden sign choice of a
# square root when SymPy conjugates an unsimplified positive sum of squares.
weight35 = sp.Symbol('weight35',real=True)
Rformal35 = sp.Symbol('Rformal35',positive=True)
Dpositive35 = sp.Symbol('Dpositive35',positive=True)
Dformal35 = (1-weight35)**2+weight35**2
zformal35 = sp.Matrix([
    nu35*sp.exp(-rho35)*sp.conjugate(tc35),sp.exp(rho35)])
Vformal35 = sp.Matrix([
    (1-weight35)*zformal35[0]/sp.sqrt(Rformal35),
    (1-weight35)*zformal35[1]/sp.sqrt(Rformal35),
    weight35*phase35])/sp.sqrt(Dpositive35)
formal_norm35 = (Vformal35.conjugate().T*Vformal35)[0]
Pformal35 = Vformal35*Vformal35.conjugate().T
P35 = Pformal35.subs({Rformal35:Rf35,Dpositive35:D35,
                     weight35:flat35})
norm35 = (
    sp.simplify((zformal35.conjugate().T*zformal35)[0]-Rf35)==0
    and sp.simplify(formal_norm35
       -((1-weight35)**2*Rf35/Rformal35+weight35**2)/Dpositive35)==0
    and sp.simplify(formal_norm35.subs(
        {Rformal35:Rf35,Dpositive35:Dformal35})-1)==0
    and sp.simplify(Dformal35
                    -(sp.Rational(1,2)
                      +2*(weight35-sp.Rational(1,2))**2))==0
    and Vformal35.subs({Rformal35:Rf35,Dpositive35:D35,
                        weight35:flat35})==V35
    and sp.simplify(sp.trace(Pformal35)-formal_norm35)==0)
# Embed the three nonzero components in the ten source slots.  Thus the
# 10x10 projector and every derivative have exactly the corresponding
# 3x3 active block; all other entries vanish identically.
Es35 = sp.zeros(10,3)
En35 = sp.zeros(10,3)
for E35,indices35 in ((Es35,(0,1,9)),(En35,(4,5,9))):
    for j35,idx35 in enumerate(indices35):
        E35[idx35,j35] = 1
V10s35,V10n35 = Es35*V35,En35*V35
P10s35,P10n35 = Es35*P35*Es35.T,En35*P35*En35.T
v_rank35 = sp.Matrix(sp.symbols('v0_35 v1_35 v2_35'))
b_rank35 = sp.Matrix(sp.symbols('b0_35 b1_35 b2_35'))
Q_rank35 = v_rank35*b_rank35.T
support35 = (
    (Es35*Vformal35)*(Es35*Vformal35).conjugate().T
        ==Es35*Pformal35*Es35.T
    and (En35*Vformal35)*(En35*Vformal35).conjugate().T
        ==En35*Pformal35*En35.T
    and Es35.T*Es35==sp.eye(3)
    and En35.T*En35==sp.eye(3)
    and P10s35==Es35*P35*Es35.T
    and P10n35==En35*P35*En35.T
    and sp.simplify(sp.trace(P35)-1)==0
    and all(sp.expand(entry)==0 for entry in
            Q_rank35**2-(b_rank35.T*v_rank35)[0]*Q_rank35))
rep('C35b source ten-slot collar projector has positive norm and rank one',
    norm35 and support35)

# Type (d): two one-point observations at the exact round witnesses in
# C35a.  These measured values are separate from the proof D>=1/2 above.
import math
def measured_norm35(z,north):
    z=float(z)
    A=(4*math.pi)**2
    B=math.exp(4)
    aa=A*(1-z)
    if north:
        s=(-aa+math.sqrt(aa*aa+4*B*z))/2
    else:
        s=(-aa+math.sqrt(aa*aa+4*B*z))/(2*B)
    q=5*math.sqrt(s)-1.5
    weight=math.exp(-1/q)/(
        math.exp(-1/q)+math.exp(-1/(1-q)))
    return (1-weight)**2+weight**2
measured_s35 = measured_norm35(z_s35,False)
measured_n35 = measured_norm35(z_n35,True)
rep('C35b2 two one-point source norm margins (numerical observation, d)',
    abs(measured_s35-0.9773063821066781)<1e-12
    and abs(measured_n35-0.6619583826095533)<1e-12
    and measured_s35>0.5 and measured_n35>0.5)
print('C35b2 south/north D and measured margins to 1/2:',
      measured_s35,measured_n35,
      measured_s35-0.5,measured_n35-0.5)

# Independently assemble the full second jet from the elementary source
# nodes.  In particular, the flat step is differentiated through its
# logistic exponent L=-1/q+1/(1-q); no unnamed derivative of a weight is
# left in the certificate.  The fifth variable H has zero jets throughout
# this plateau, as established by the strict guards in C35a.  Constants
# nu=1 and nu=e^-2 are covered by the same positive-symbol identities.
vars35 = (ar35,ai35,rho35,one35,H35)
zero_g35 = [sp.Integer(0)]*5
zero_h35 = [[sp.Integer(0)]*5 for _ in range(5)]
def jet35(value,gradient,hessian):
    return value,gradient,hessian
def product35(A,B):
    av,ag,ah=A; bv,bg,bh=B
    return jet35(av*bv,
        [ag[i]*bv+av*bg[i] for i in range(5)],
        [[ah[i][j]*bv+ag[i]*bg[j]+ag[j]*bg[i]+av*bh[i][j]
          for j in range(5)] for i in range(5)])
def invroot35(A):
    av,ag,ah=A
    return jet35(1/sp.sqrt(av),
        [-ag[i]/(2*av**sp.Rational(3,2)) for i in range(5)],
        [[3*ag[i]*ag[j]/(4*av**sp.Rational(5,2))
          -ah[i][j]/(2*av**sp.Rational(3,2))
          for j in range(5)] for i in range(5)])
def same_jet35(A,B):
    av,ag,ah=A; bv,bg,bh=B
    return (sp.simplify(av-bv)==0
        and all(sp.simplify(ag[i]-bg[i])==0 for i in range(5))
        and all(sp.simplify(ah[i][j]-bh[i][j])==0
                for i in range(5) for j in range(i,5)))
def direct_jet35(expr):
    return jet35(expr,[expr.diff(v) for v in vars35],
        [[expr.diff(vi,vj) for vj in vars35] for vi in vars35])

qg35 = [5*ar35/r35,5*ai35/r35,0,0,0]
qh35 = [[5*ai35**2/r35**3,-5*ar35*ai35/r35**3,0,0,0],
        [-5*ar35*ai35/r35**3,5*ar35**2/r35**3,0,0,0],
        *[zero_g35[:] for _ in range(3)]]
Lprime35 = 1/q35**2+1/(1-q35)**2
Lsecond35 = -2/q35**3+2/(1-q35)**3
sg35 = [flat35*(1-flat35)*Lprime35*qg35[i]
        for i in range(5)]
sh35 = [[flat35*(1-flat35)*(
           ((1-2*flat35)*Lprime35**2+Lsecond35)
               *qg35[i]*qg35[j]+Lprime35*qh35[i][j])
         for j in range(5)] for i in range(5)]
sjet35 = jet35(flat35,sg35,sh35)
cjet35 = jet35(1-flat35,[-v for v in sg35],
               [[-v for v in row] for row in sh35])
Rg35 = [2*nu35**2*sp.exp(-2*rho35)*ar35,
        2*nu35**2*sp.exp(-2*rho35)*ai35,
        -2*nu35**2*sp.exp(-2*rho35)*rsq35
           +2*sp.exp(2*rho35),0,0]
Rh35 = [
    [2*nu35**2*sp.exp(-2*rho35),0,
     -4*nu35**2*sp.exp(-2*rho35)*ar35,0,0],
    [0,2*nu35**2*sp.exp(-2*rho35),
     -4*nu35**2*sp.exp(-2*rho35)*ai35,0,0],
    [-4*nu35**2*sp.exp(-2*rho35)*ar35,
     -4*nu35**2*sp.exp(-2*rho35)*ai35,
     4*nu35**2*sp.exp(-2*rho35)*rsq35
       +4*sp.exp(2*rho35),0,0],
    zero_g35[:],zero_g35[:]]
Rjet35 = jet35(Rf35,Rg35,Rh35)
z1jet35 = jet35(zformal35[0],
    [nu35*sp.exp(-rho35),-sp.I*nu35*sp.exp(-rho35),
     -zformal35[0],0,0],
    [[0,0,-nu35*sp.exp(-rho35),0,0],
     [0,0,sp.I*nu35*sp.exp(-rho35),0,0],
     [-nu35*sp.exp(-rho35),sp.I*nu35*sp.exp(-rho35),
      zformal35[0],0,0],zero_g35[:],zero_g35[:]])
z2jet35 = jet35(zformal35[1],
    [0,0,zformal35[1],0,0],
    [[0]*5,[0]*5,[0,0,zformal35[1],0,0],
     zero_g35[:],zero_g35[:]])
phasejet35 = jet35(phase35,
    [0,0,0,2*sp.pi*sp.I*phase35,0],
    [[0]*5,[0]*5,[0]*5,
     [0,0,0,-4*sp.pi**2*phase35,0],zero_g35[:]])
Djet35 = jet35(D35,
    [2*(2*flat35-1)*sg35[i] for i in range(5)],
    [[4*sg35[i]*sg35[j]+2*(2*flat35-1)*sh35[i][j]
      for j in range(5)] for i in range(5)])
q_scalar35 = sp.Symbol('q_scalar35',positive=True)
s_scalar35 = sp.exp(-1/q_scalar35)/(
    sp.exp(-1/q_scalar35)+sp.exp(-1/(1-q_scalar35)))
Lp_scalar35 = 1/q_scalar35**2+1/(1-q_scalar35)**2
Lpp_scalar35 = -2/q_scalar35**3+2/(1-q_scalar35)**3
flat_scalar_jets35 = (
    sp.simplify(sp.diff(s_scalar35,q_scalar35)
        -s_scalar35*(1-s_scalar35)*Lp_scalar35)==0
    and sp.simplify(sp.diff(s_scalar35,q_scalar35,2)
        -s_scalar35*(1-s_scalar35)*(
          (1-2*s_scalar35)*Lp_scalar35**2+Lpp_scalar35))==0)
# The q and source-frame jets are checked directly; the one-variable
# flat rule and the scalar composition rule supply all 25 flat entries.
primitive_jets35 = (flat_scalar_jets35
    and all(same_jet35(J,direct_jet35(expr))
       for J,expr in ((jet35(q35,qg35,qh35),q35),
                      (Rjet35,Rf35),
                      (z1jet35,zformal35[0]),
                      (z2jet35,zformal35[1]),
                      (phasejet35,phase35)))
    and sp.simplify(D35-(1-2*flat35+2*flat35**2))==0)
radial_jet35 = invroot35(Rjet35)
norm_jet35 = invroot35(Djet35)
Vjets35 = [product35(product35(cjet35,
                 product35(zjet,radial_jet35)),norm_jet35)
           for zjet in (z1jet35,z2jet35)]
Vjets35.append(product35(product35(sjet35,phasejet35),norm_jet35))
V_first35 = [V35.diff(v) for v in vars35]
V_second35 = [[V35.diff(vi,vj) for vj in vars35]
              for vi in vars35]
# Verify the two generic algebra rules on arbitrary second-order Taylor
# data in two independent real directions.  Relabelling the directions
# gives every one of the 5x5 index pairs without expanded expressions.
uu35,vv35 = sp.symbols('uu35 vv35',real=True)
a0g35 = sp.Symbol('a0g35',positive=True)
aug35,avg35,auug35,auvg35,avvg35 = sp.symbols(
    'aug35 avg35 auug35 auvg35 avvg35',real=True)
b0g35,bug35,bvg35,buug35,buvg35,bvvg35 = sp.symbols(
    'b0g35 bug35 bvg35 buug35 buvg35 bvvg35')
Ap35 = (a0g35+aug35*uu35+avg35*vv35
    +auug35*uu35**2/2+auvg35*uu35*vv35+avvg35*vv35**2/2)
Bp35 = (b0g35+bug35*uu35+bvg35*vv35
    +buug35*uu35**2/2+buvg35*uu35*vv35+bvvg35*vv35**2/2)
Ag35 = jet35(a0g35,[aug35,avg35,0,0,0],
    [[auug35,auvg35,0,0,0],[auvg35,avvg35,0,0,0],
     *[zero_g35[:] for _ in range(3)]])
Bg35 = jet35(b0g35,[bug35,bvg35,0,0,0],
    [[buug35,buvg35,0,0,0],[buvg35,bvvg35,0,0,0],
     *[zero_g35[:] for _ in range(3)]])
def generic_rule35(jet,expression):
    value,grad,hess=jet
    at_zero={uu35:0,vv35:0}
    directions=(uu35,vv35)
    return (sp.simplify(value-expression.subs(at_zero))==0
      and all(sp.simplify(grad[i]-sp.diff(expression,directions[i])
                           .subs(at_zero))==0 for i in range(2))
      and all(sp.simplify(hess[i][j]
          -sp.diff(expression,directions[i],directions[j])
            .subs(at_zero))==0
          for i in range(2) for j in range(i,2)))
generic_chain35 = (
    generic_rule35(product35(Ag35,Bg35),Ap35*Bp35)
    and generic_rule35(invroot35(Ag35),1/sp.sqrt(Ap35)))
f0g35,f1g35,f2g35 = sp.symbols('f0g35 f1g35 f2g35')
generic_composite35 = (f0g35+f1g35*(Ap35-a0g35)
                      +f2g35*(Ap35-a0g35)**2/2)
generic_chain35 &= generic_rule35(jet35(f0g35,
    [f1g35*x for x in Ag35[1]],
    [[f2g35*Ag35[1][i]*Ag35[1][j]+f1g35*Ag35[2][i][j]
      for j in range(5)] for i in range(5)]),generic_composite35)
V_chain35 = (generic_chain35
    and all(sp.simplify(Vjets35[k][0]-V35[k])==0
            for k in range(3))
    and all(len(J[1])==5 and len(J[2])==5
            and all(len(row)==5 for row in J[2])
            for J in Vjets35))
# For arbitrary complex component functions, the following is the exact
# second-order product rule for P_ij=V_i*conj(V_j).  The three independently
# checked V-jets above therefore determine every one of the 10x10 projector
# jets by sparse insertion; no new chart assumptions enter this step.
P_jet_formula35 = [[jet35(
    Vjets35[i][0]*sp.conjugate(Vjets35[j][0]),
    [Vjets35[i][1][k]*sp.conjugate(Vjets35[j][0])
      +Vjets35[i][0]*sp.conjugate(Vjets35[j][1][k])
      for k in range(5)],
    [[Vjets35[i][2][k][ell]*sp.conjugate(Vjets35[j][0])
      +Vjets35[i][1][k]*sp.conjugate(Vjets35[j][1][ell])
      +Vjets35[i][1][ell]*sp.conjugate(Vjets35[j][1][k])
      +Vjets35[i][0]*sp.conjugate(Vjets35[j][2][k][ell])
      for ell in range(5)] for k in range(5)])
      for j in range(3)] for i in range(3)]
P_chain35 = (
    generic_rule35(product35(Ag35,
        jet35(sp.conjugate(b0g35),
          [sp.conjugate(x) for x in Bg35[1]],
          [[sp.conjugate(x) for x in row] for row in Bg35[2]])),
          Ap35*sp.conjugate(Bp35))
    and len(P_jet_formula35)==3
    and all(len(row)==3 for row in P_jet_formula35)
    and all(len(J[1])==5 and len(J[2])==5
            and all(len(row)==5 for row in J[2])
            for row in P_jet_formula35 for J in row)
    # C35b proves the formal source value identity before substituting the
    # positive radicands.  The generic conjugate-product rule gives every
    # indexed derivative, so no expanded 10x10 Hessian is required.
    and Pformal35==Vformal35*Vformal35.conjugate().T
    and P35==Pformal35.subs({Rformal35:Rf35,
                              Dpositive35:D35,weight35:flat35})
    and all(P_jet_formula35[i][j][0]
            ==Vjets35[i][0]*sp.conjugate(Vjets35[j][0])
            for i in range(3) for j in range(3)))
wrong_mixed35 = (auvg35*sp.conjugate(b0g35)
    +avg35*sp.conjugate(bug35)
    +a0g35*sp.conjugate(buvg35))
negative_control35 = (
    sp.simplify(sp.diff(Ap35*sp.conjugate(Bp35),uu35,vv35)
      .subs({uu35:0,vv35:0})-wrong_mixed35
      -aug35*sp.conjugate(bvg35))==0
    and aug35*sp.conjugate(bvg35)!=0)
full_jets35 = (
    primitive_jets35 and V_chain35 and P_chain35
    and negative_control35
    and all(V_second35[i][j]==V_second35[j][i]
        for i in range(5) for j in range(5))
    and V_first35[4]==sp.zeros(3,1)
    and all(V_second35[4][j]==sp.zeros(3,1)
            for j in range(5))
    and all(P_jet_formula35[i][j][1][4]==0
            and all(P_jet_formula35[i][j][2][4][k]==0
                    for k in range(5))
            for i in range(3) for j in range(3))
    and any(V_first35[0][k]!=0 for k in range(3))
    and any(V_first35[2][k]!=0 for k in range(3))
    and any(V_first35[3][k]!=0 for k in range(3))
    and any(P_jet_formula35[i][j][1][0]!=0
            for i in range(3) for j in range(3))
    and any(P_jet_formula35[i][j][2][0][3]!=0
            for i in range(3) for j in range(3)))
rep('C35c actual collar V/P: all five-variable first and second jets',
    full_jets35)

# C36: a signed six-real determinant on the two strict C35 neighborhoods.
# Source order is (n0,n1,X,theta,phi2,phi3), with n2<0 and
# Z=sqrt(1-X^2)*exp(2pi*i*theta)*((n0+i*n1)/(2c),c),
# c=sqrt((1-n2)/2).  Target order is (Re t,Im t,Re z,Im z,Re w,Im w).
# The auxiliary W/torus order is (Re t,Im t,x4,x1,phi2,phi3).
# This is a local determinant, not a global free-slab certificate.
height_inverse_tree36 = source_tree30(cusp30.elementary_cusp_logradius)
moment_tree36 = source_tree30(cusp30.inverse_height_moment)
terminal_tree36 = source_tree30(cusp30._terminal_height_data)
height_G_tree36 = source_tree30(cusp30._height_G)
inner_tree36 = source_tree30(cusp30._elementary_inner)
inclusion_tree36 = source_tree30(cusp30._elementary_pole_inclusion)
source_det36 = (
    source_counts30 and project_source30 and source_collar34 and source_torus35
    and matching_assignment30(raw_tree31,'u1',
        'mp.mpc(1-s2,2*y[0])/(1+s2)')
    and matching_assignment30(raw_tree31,'u2',
        'mp.mpc(2*y[1],2*y[2])/(1+s2)')
    and matching_assignment30(raw_tree31,'inverse_ratio','abs(u1)>abs(u2)')
    and matching_assignment30(raw_tree31,'ratio',
        'u2/u1 if inverse_ratio else u1/u2')
    and matching_assignment30(raw_tree31,'root',
        'mp.sqrt(b*b+4*B*rr)')
    and matching_statement31(raw_tree31,
        'A,B = (4*mp.pi)**2,mp.exp(4)')
    and matching_assignment30(raw_tree31,'b','A*(1-rr)')
    and matching_assignment30(raw_tree31,'small',
        'ratio*mp.sqrt((2*B if inverse_ratio else 2)/(b+root))')
    and matching_assignment30(raw_tree31,'rho',
        '(elementary_cusp_logradius if elementary else '
        'elliptic_cusp_logradius)(a,tc,nu,dps=dps)')
    and matching_assignment30(height_inverse_tree36,'P',
        'mp.log(nu)+mp.log(r)-2*r*r/(1-r*r)')
    and matching_assignment30(height_inverse_tree36,'q','mp.exp(2*P)')
    and matching_assignment30(height_inverse_tree36,'n',
        "mp.floor(a+mp.mpf('.5'))")
    and matching_assignment30(height_inverse_tree36,'b','a-n')
    and matching_assignment30(height_inverse_tree36,'lam',
        "scalar_flat_step((r-mp.mpf('.3'))/mp.mpf('.2'))")
    and matching_assignment30(height_inverse_tree36,'f',
        'inverse_height_moment(b,q,lam)')
    and matching_assignment30(height_inverse_tree36,'root',
        'mp.sqrt(f*f+4*q)')
    and matching_assignment30(height_inverse_tree36,'Y',
        '2*q/(root+f) if f>=0 else (root-f)/2')
    and returned_expr30(height_inverse_tree36,
        '+(n*P+mp.log(Y)/2)')
    and matching_assignment30(moment_tree36,'aa','abs(a)')
    and matching_statement31(moment_tree36,
        "if lam==1 or aa<=mp.mpf('.1'):\n    return _height_G(a,amp,d)")
    and matching_if_test32(terminal_tree36,'lam==0')
    and matching_assignment30(terminal_tree36,'k','-mp.log(q)')
    and matching_assignment30(terminal_tree36,'amp','lam*k')
    and matching_assignment30(terminal_tree36,'d',
        '1-lam+2*lam*mp.sqrt(q)*k')
    and returned_expr30(height_G_tree36,
        'mp.asinh(amp*a/d)/amp if inverse else d*mp.sinh(amp*a)/amp')
    and returned_expr30(inner_tree36,
        'sum(mp.conj(a)*b for a,b in zip(u,v))')
    and matching_assignment30(project_tree30,'vector',
        '_elementary_pole_inclusion(xi)')
    and matching_assignment30(project_tree30,'pc',
        '2*xi[0]*mp.conj(xi[1])')
    and matching_assignment30(project_tree30,'n',
        '(pc.real,pc.imag,abs(xi[0])**2-abs(xi[1])**2)')
    and matching_assignment30(slice_tree30,'xi',
        'tuple(a/mp.sqrt(1-X*X) for a in Z)')
    and matching_assignment30(slice_tree30,'phase',
        '_elementary_inner(ref,result)')
    and matching_if_test32(correct_tree30,'distance<alpha')
    and matching_if_test32(slice_tree30,
        'distance<alpha/2 and not force_slab')
    and returned_expr30(slice_tree30,
        'kind,(t,one+mp.arg(phase)/(2*mp.pi),four,lt)')
    and matching_assignment30(restore_tree30,'C','list(phi)')
    and matching_statement31(restore_tree30,
        't,one,four,lt = data[:4]')
    and matching_statement31(restore_tree30,
        'C = [v+beta*d for v,d in zip(C,D)]')
    and returned_expr30(restore_tree30,
        "'mid',(z-mp.j*t*C[0]-t*C[1],"
        "w+(1-t*t)*C[0]+mp.j*(1+t*t)*C[1],t)"))
d_assignments36 = [node.value for node in ast.walk(restore_tree30)
                   if isinstance(node,ast.Assign)
                   and any(isinstance(target,ast.Name) and target.id=='D'
                           for target in node.targets)]
source_det36 &= (len(d_assignments36)==2
    and all(not any(isinstance(sub,ast.Name) and sub.id in ('C','phi')
                    for sub in ast.walk(expr))
            for expr in d_assignments36))

# Both exact witnesses lie outside the phase-correction and pole bands:
# south 1<R<3/2, north 3<R<31/10, and X(R) increases for R>0.
Rs36 = sp.sqrt(y0_w35**2+(sp.Rational(12,5)-1)**2)
Rn36 = sp.sqrt(y0_w35**2+(sp.Integer(4)-1)**2)
Rvar36 = sp.Symbol('Rvar36',positive=True)
X_of_R36 = (Rvar36-1)/sp.sqrt(1+Rvar36**2)
slab_witness36 = (
    root_cell35 and height_cell35 and branch_margin35
    and 1<Rs36<sp.Rational(3,2)
    and 3<Rn36<sp.Rational(31,10)
    and sp.simplify(sp.diff(X_of_R36,Rvar36)
        -(Rvar36+1)/(1+Rvar36**2)**sp.Rational(3,2))==0
    and X_of_R36.subs(Rvar36,sp.Rational(3,2))<sp.Rational(1,2)
    and X_of_R36.subs(Rvar36,sp.Rational(31,10))<sp.Rational(3,4)
    and sp.Rational(1,4)>sp.Rational(1,100000)
    and matching_assignment30(slice_tree30,'alpha',"mp.mpf('1e-5')")
    and matching_assignment30(correct_tree30,'alpha',"mp.mpf('1e-5')"))

# The projection fold preserves the input scalar phase at every one of its
# certified N steps: each projector depends on (n,h), and a unit scalar p
# takes vector to p*vector and its inner product to p*phase.  The final
# source middle one-coordinate adds Arg(phase)/(2pi), hence has derivative
# +1 in the local theta chart.  The correction is inactive at these X.
th36 = sp.Symbol('th36',real=True)
p36 = sp.exp(2*sp.pi*sp.I*th36)
v036,v136,r036,r136 = sp.symbols('v036 v136 r036 r136')
xi036,xi136 = sp.symbols('xi036 xi136')
inner36 = sp.conjugate(r036)*v036+sp.conjugate(r136)*v136
phase_fold36 = (
    sp.simplify(2*(p36*xi036)*sp.conjugate(p36*xi136)
                -2*xi036*sp.conjugate(xi136))==0
    and sp.simplify(sp.conjugate(p36*xi036)*(p36*xi036)
                    -sp.conjugate(xi036)*xi036)==0
    and sp.simplify(sp.conjugate(p36*xi136)*(p36*xi136)
                    -sp.conjugate(xi136)*xi136)==0
    and sp.simplify(sp.conjugate(r036)*(p36*v036)
        +sp.conjugate(r136)*(p36*v136)-p36*inner36)==0
    and sp.simplify(abs(p36*inner36)-abs(inner36))==0
    and all(sp.simplify(z)==0 for z in (
        p36*(xi036-sp.I*xi136)/sp.sqrt(2)
          -(p36*xi036-sp.I*p36*xi136)/sp.sqrt(2),
        p36*(xi036+xi136)/sp.sqrt(2)
          -(p36*xi036+p36*xi136)/sp.sqrt(2)))
    and matching_assignment30(inclusion_tree36,'pair',
        '((a+b)/mp.sqrt(2),(a-b)/mp.sqrt(2)) if infinity else '
        '((a-mp.j*b)/mp.sqrt(2),(a+mp.j*b)/mp.sqrt(2))'))
rep('C36a source path, strict slab witnesses and finite-fold fibre derivative +1',
    source_det36 and slab_witness36 and phase_fold36)

# The raw quotient map is a positive radial scale followed by a signed
# permutation and a translation: y=(lambda*n0,1-lambda*n2,-lambda*n1).
# In the n2<0 chart, det D y / D(n0,n1,X)=-lambda^2 lambda'/n2>0.
X36 = sp.Symbol('X36',real=True)
U36 = sp.sqrt(2-X36**2)
lambda36 = (1-X36**2)/(1-X36*U36)
lambda_alt36 = (U36+X36)/(U36-X36)
lambda_logprime36 = 2/(U36*(1-X36**2))
n036,n136 = sp.symbols('n036 n136',real=True)
n236 = -sp.sqrt(1-n036**2-n136**2)
L36,Lp36 = sp.symbols('L36 Lp36',positive=True)
dy36 = sp.Matrix([
    [L36,0,Lp36*n036],
    [-L36*sp.diff(n236,n036),-L36*sp.diff(n236,n136),-Lp36*n236],
    [0,-L36,-Lp36*n136]])
radial_base36 = (
    sp.simplify(lambda36-lambda_alt36)==0
    and sp.simplify(sp.diff(lambda_alt36,X36)/lambda_alt36
                    -lambda_logprime36)==0
    and sp.simplify(dy36.det()+L36**2*Lp36/n236)==0)

# At y=(y0,Y,0), S=1+y0^2+Y^2.  The ratio/height charts are rational
# except for the phase differential d Arg z=Im(dz/z), which is written
# as an explicit real one-form.  Their signed Jacobians have opposite
# signs in the two cusp choices.
y036,Y36,y236 = sp.symbols('y036 Y36 y236',real=True)
S36 = 1+y036**2+Y36**2+y236**2
rs36 = (2-S36+2*sp.I*y036)/(2*(Y36+sp.I*y236))
rn36 = 2*(Y36+sp.I*y236)/(2-S36+2*sp.I*y036)
as_grad36 = [sp.Integer(0),
             -y236/(2*sp.pi*(Y36**2+y236**2)),
             Y36/(2*sp.pi*(Y36**2+y236**2))]
Un36,Vn36 = S36-2,2*y036
an_grad36 = [(Un36*sp.diff(Vn36,v)-Vn36*sp.diff(Un36,v))
              /(2*sp.pi*(Un36**2+Vn36**2))
             for v in (y036,Y36,y236)]
def ratio_jac36(ratio,phase_gradient):
    return sp.Matrix([
        [sp.diff(sp.re(ratio),v) for v in (y036,Y36,y236)],
        [sp.diff(sp.im(ratio),v) for v in (y036,Y36,y236)],
        phase_gradient]).det().subs(y236,0)
Saxis36 = 1+y036**2+Y36**2
Qaxis36 = (1-y036**2-Y36**2)**2+4*y036**2
ratio_s36 = ratio_jac36(rs36,as_grad36)
ratio_n36 = ratio_jac36(rn36,an_grad36)
ratio_det36 = (
    sp.simplify(ratio_s36-Saxis36/(4*sp.pi*Y36**4))==0
    and sp.simplify(ratio_n36+4*Saxis36/(sp.pi*Qaxis36**2))==0)

# Source small=ratio*sqrt(s/z), s=|tc|^2, z=|ratio|^2.  For any
# positive radial factor g(z), det_R D[ratio*g(z)]=d(z*g(z)^2)/dz.
# The north t=1/tc adds the positive holomorphic inversion factor 1/s^2.
s36,A36,B36 = sp.symbols('s36 A36 B36',positive=True)
zs36 = (B36*s36**2+A36*s36)/(A36*s36+1)
zn36 = (s36**2+A36*s36)/(A36*s36+B36)
zps36 = (A36*B36*s36**2+2*B36*s36+A36)/(A36*s36+1)**2
zpn36 = (A36*s36**2+2*B36*s36+A36*B36)/(A36*s36+B36)**2
zrad36 = sp.Symbol('zrad36',positive=True)
ru36,rv36 = sp.symbols('ru36 rv36',real=True)
gfunc36 = sp.Function('gfunc36')
radial_argument36 = ru36**2+rv36**2
radial_map36 = sp.Matrix([
    ru36*gfunc36(radial_argument36),
    rv36*gfunc36(radial_argument36)])
inversion_map36 = sp.Matrix([
    ru36/radial_argument36,-rv36/radial_argument36])
radial_profile36 = (
    sp.simplify(sp.diff(zs36,s36)-zps36)==0
    and sp.simplify(sp.diff(zn36,s36)-zpn36)==0
    and sp.simplify(radial_map36.jacobian([ru36,rv36]).det()
        -sp.diff(zrad36*gfunc36(zrad36)**2,zrad36)
          .subs(zrad36,radial_argument36))==0
    and sp.simplify(inversion_map36.jacobian([ru36,rv36]).det()
                    -1/radial_argument36**2)==0)

# On |a|<1/10, the source terminal-cut moment uses exactly
# f=d*sinh(amp*a)/amp.  Here q=e^(2P) in (0,1), k=-log q>0,
# lambda_flat in (0,1), amp=lambda_flat*k>0, and
# d=(1-lambda_flat)+2 lambda_flat sqrt(q) k>0.  Thus f_a>0.
a36 = sp.Symbol('a36',real=True)
amp36,d36,q36 = sp.symbols('amp36 d36 q36',positive=True)
f36 = d36*sp.sinh(amp36*a36)/amp36
Delta36 = sp.sqrt(f36**2+4*q36)
Yheight36 = (Delta36-f36)/2
fa36 = d36*sp.cosh(amp36*a36)
Eflat36,Fflat36,kpos36 = sp.symbols(
    'Eflat36 Fflat36 kpos36',positive=True)
lam_manifest36 = Eflat36/(Eflat36+Fflat36)
amp_manifest36 = lam_manifest36*kpos36
d_manifest36 = 1-lam_manifest36+2*lam_manifest36*sp.sqrt(q36)*kpos36
d_manifest_positive36 = ((Fflat36+2*Eflat36*sp.sqrt(q36)*kpos36)
                         /(Eflat36+Fflat36))
local_height36 = (
    sp.simplify(sp.diff(f36,a36)-fa36)==0
    and sp.simplify(sp.diff(sp.log(Yheight36)/2,a36)
                    +fa36/(2*Delta36))==0
    and sp.simplify(2*q36/(Delta36+f36)-Yheight36)==0
    and sp.ask(sp.Q.positive(amp_manifest36))
    and sp.simplify(d_manifest36-d_manifest_positive36)==0
    and sp.ask(sp.Q.positive(d_manifest_positive36))
    and sp.ask(sp.Q.positive(fa36))
    and sp.Rational(1,10)<sp.Rational(1,2)
    and sp.Rational(9,100)>0
    and sp.Rational(1,4)<1)
# Cramer's x4 is the third coordinate of R^{-1}rhs; rho enters only as
# -rho in the first right-hand entry.  Therefore x4_rho is the negative
# (3,1) inverse entry, and x4_a has sign + south / - north.
r2_36 = ar34**2+ai34**2
cramer_factor36 = (
    sp.factor(R34.cofactor(0,2)-(1-r2_36**2))==0
    and sp.simplify(R34.inv()[2,0]
                    -(1-r2_36**2)/R34.det())==0
    and sp.simplify(R34.det()-det_expected34)==0)
rep('C36b positive raw/profile/height factors and signed ratio/Cramer Jacobians',
    radial_base36 and ratio_det36 and radial_profile36
    and local_height36 and cramer_factor36)

# The following are the two *manifest* positive base determinants in the
# fixed order.  Let m=-n2>0 and let Dr be det R>0.  The profile roots lie
# in 9/100<s<1/4, so 1-s^2>0.  All remaining factors are positive.
m36,lambda_p36,Spos36,Ypos36,Qpos36 = sp.symbols(
    'm36 lambda_p36 Spos36 Ypos36 Qpos36',positive=True)
zps_pos36,zpn_pos36,fa_pos36,Delta_pos36,Dr_pos36 = sp.symbols(
    'zps_pos36 zpn_pos36 fa_pos36 Delta_pos36 Dr_pos36',positive=True)
south_base36 = (L36**2*lambda_p36/m36
    *Spos36/(4*sp.pi*Ypos36**4)
    /zps_pos36
    *(1-s36**2)*fa_pos36/(2*Delta_pos36*Dr_pos36))
north_base36 = (L36**2*lambda_p36/m36
    *(-4*Spos36/(sp.pi*Qpos36**2))
    /(s36**2*zpn_pos36)
    *(1-1/s36**2)*fa_pos36/(2*Delta_pos36*Dr_pos36))
signed_base36 = (
    sp.simplify(south_base36
        -L36**2*lambda_p36*Spos36*(1-s36**2)*fa_pos36
         /(8*sp.pi*m36*Ypos36**4*zps_pos36
           *Delta_pos36*Dr_pos36))==0
    and sp.simplify(north_base36
        -2*L36**2*lambda_p36*Spos36*(1-s36**2)*fa_pos36
         /(sp.pi*m36*Qpos36**2*s36**4*zpn_pos36
           *Delta_pos36*Dr_pos36))==0)
strict_base_margin36 = (root_cell35
    and 1-sp.Rational(1,4)**2==sp.Rational(15,16)
    and sp.Rational(15,16)>0
    and bool(sp.ask(sp.Q.positive(zps36)))
    and bool(sp.ask(sp.Q.positive(zpn36))))

# T^2 restoration adds C=phi+D(base), so its phi block is I_2.  In the
# source middle formula, the remaining (one,four,C2,C3)->(Re z,Im z,
# Re w,Im w) determinant is +2pi det R.  The auxiliary order above uses
# (four,one), one transposition, so its restoration factor is negative.
t36 = ar34+sp.I*ai34
Mrestore36 = sp.Matrix([
    [2*sp.pi,sp.re(-sp.I*L34),sp.re(-sp.I*t36),sp.re(-t36)],
    [0,sp.im(-sp.I*L34),sp.im(-sp.I*t36),sp.im(-t36)],
    [0,sp.re(2*t36),sp.re(1-t36**2),sp.re(sp.I*(1+t36**2))],
    [0,sp.im(2*t36),sp.im(1-t36**2),sp.im(sp.I*(1+t36**2))]])
restoration_det36 = sp.simplify(Mrestore36.det()
                                 -2*sp.pi*det_expected34)==0
full_south36 = -2*sp.pi*Dr_pos36*south_base36
full_north36 = -2*sp.pi*Dr_pos36*north_base36
full_signed36 = (
    sp.simplify(full_south36
        +L36**2*lambda_p36*Spos36*(1-s36**2)*fa_pos36
         /(4*m36*Ypos36**4*zps_pos36*Delta_pos36))==0
    and sp.simplify(full_north36
        +4*L36**2*lambda_p36*Spos36*(1-s36**2)*fa_pos36
         /(m36*Qpos36**2*s36**4*zpn_pos36*Delta_pos36))==0)
rep('C36c two strict-cell full determinants: positive adapted base, '
    'negative fixed-order Y chart',
    source_det36 and slab_witness36 and phase_fold36
    and signed_base36 and strict_base_margin36
    and restoration_det36 and full_signed36
    and bool(sp.Rational(1,4)<1)
    and bool(sp.pi>0))

# C37: one open free-slab chamber in which EVERY fixed-count projection
# reference takes the same early south-cusp branch.  The round witness is
# Z=(1/100,i/200), X=-sqrt(7999/8000).  The local source order is
# (n0,n1,X,theta,phi2,phi3), n2=+sqrt(1-n0^2-n1^2), with the upper Hopf
# reference; the holomorphic target order is the three complex cusp outputs.
# This section certifies the source-defined finite expression symbolically
# by induction, with max_work=None in mathematical expression mode.  The
# evaluator's default work guard prevents literal evaluation at the certified
# astronomical N.  Full six Hessians and the named J are not asserted here.
def exact_if_body37(tree,test,body):
    wanted = ast.dump(ast.parse(test,mode='eval').body,
                      include_attributes=False)
    matches = [node for node in ast.walk(tree)
               if isinstance(node,ast.If)
               and ast.dump(node.test,include_attributes=False)==wanted]
    expected = [ast.dump(node,include_attributes=False)
                for node in ast.parse(body).body]
    return (len(matches)==1
            and [ast.dump(node,include_attributes=False)
                 for node in matches[0].body]==expected)


def ordered_nodes37(body,snippets):
    actual = [ast.dump(node,include_attributes=False) for node in body]
    cursor = -1
    for snippet in snippets:
        expected = ast.parse(snippet).body
        if len(expected)!=1:
            return False
        key = ast.dump(expected[0],include_attributes=False)
        found = next((i for i in range(cursor+1,len(actual))
                      if actual[i]==key),None)
        if found is None:
            return False
        cursor = found
    return True


core_tree37 = source_tree30(cusp30.elementary_cusp_reference)
moment_tree37 = source_tree30(cusp30.inverse_height_moment)
terminal_tree37 = source_tree30(cusp30._terminal_height_data)
g_tree37 = source_tree30(cusp30._height_G)
finite_height_tree37 = source_tree30(cusp30.scalar_finite_height)
restore_tree37 = source_tree30(cusp30.elementary_restore_torus)
projection_loops37 = [node for node in ast.walk(project_tree30)
                      if isinstance(node,ast.For)
                      and ast.dump(node.iter,include_attributes=False)
                      ==ast.dump(ast.parse('range(steps+1)',mode='eval').body,
                                  include_attributes=False)]
slice_with37 = [node for node in slice_tree30.body
                if isinstance(node,ast.With)]
finite_seed37 = [node for node in ast.walk(finite_height_tree37)
                 if isinstance(node,ast.FunctionDef) and node.name=='seed']
restore_coeffs37 = [node for node in ast.walk(restore_tree37)
                    if isinstance(node,ast.FunctionDef) and node.name=='coeffs']
restore_cusp_data37 = [node for node in ast.walk(restore_tree37)
                       if isinstance(node,ast.FunctionDef)
                       and node.name=='cusp_data']
source_path37 = (
    source_counts30 and project_source30 and source_cusp32
    and source_torus35 and source_det36 and len(projection_loops37)==1
    and [ast.dump(node,include_attributes=False)
         for node in projection_loops37[0].body]
        ==[ast.dump(node,include_attributes=False) for node in ast.parse('''
h = start+(X-start)*j/steps
kind,data = elliptic_raw_base_reference(tuple((1-h*h)*a for a in n),h,
                                         dps=dps,elementary=True)
ref = elliptic_bundle_embedding(kind,data,dps=dps,elementary=True)
phase = _elementary_inner(ref,vector)
margin = min(margin,abs(phase))
if phase==0:
    raise ValueError('finite projection product vanished')
vector = tuple(a*phase/abs(phase) for a in ref)
''').body]
    and exact_if_body37(raw_tree31,
        "abs(tc)<=mp.mpf('.3') and elementary",'''
xs,yy = elementary_cusp_reference(tc,a,nu)
return ('ci' if inf else 'c0'),(xs,yy,a)
''')
    and exact_if_body37(core_tree37,'a>0','''
A = 2/(R-f) if cell==1 else (R+f)/2
x = mp.sqrt(A/e)/nu
return x,t/x
''')
    and exact_if_body37(bundle_tree32,"kind in ('c0','ci')",'''
x,y = map(mp.mpc,data[:2])
nu = mp.exp(-2) if kind=='ci' else mp.mpf(1)
H = mp.mpf(data[2]) if elementary else elliptic_height_chart(x,y,nu,dps=dps)
if abs(x*y)>mp.mpf('.3') or abs(H)>mp.mpf('.6'):
    raise ValueError('pure cusp embedding requires |xy|<=.3 and |H|<=.6')
f = frame(x,y,nu)
block = [hw(H)*q for q in f]
w = hw(H+mp.mpf('.5') if H<0 else H-mp.mpf('.5'))
if H<0 and w>0:
    rr = abs(x*y)**2
    m = nu*mp.exp(-2*rr/(1-rr))
    f = frame(1/(m*y),m*x*y*y,nu)
block.extend(w*q for q in f)
start = 4 if kind=='ci' else 0
vec[start:start+4] = block
''')
    and ordered_nodes37(core_tree37.body,[
        'e = mp.exp(-4*abs(t)**2/(1-abs(t)**2))',
        'q = e*nu*nu*abs(t)**2',
        "cell = int(mp.floor(a+mp.mpf('.5')))",
        'b = a-cell',
        'f = inverse_height_moment(b,q,0)',
        'R = mp.sqrt(f*f+4*q)'])
    and exact_if_body37(terminal_tree37,'lam==0',
        'kappa,amp,d = mp.mpf(2),mp.mpf(0),mp.mpf(1)')
    and exact_if_body37(g_tree37,'amp==0',
        'return a/d if inverse else d*a')
    and exact_if_body37(moment_tree37,
        "lam==1 or aa<=mp.mpf('.1')",
        'return _height_G(a,amp,d)')
    and len(slice_with37)==1
    and ordered_nodes37(slice_with37[0].body,[
        'params = elementary_closed_parameters()',
        "N = params['projection_steps'] if projection_steps is None else projection_steps",
        'xi = tuple(a/mp.sqrt(1-X*X) for a in Z)',
        'result = _elementary_corrected_section(xi,X,N,K,dps,diagnostics)',
        'pc = 2*xi[0]*mp.conj(xi[1])',
        'n = (pc.real,pc.imag,abs(xi[0])**2-abs(xi[1])**2)',
        'kind,data = elliptic_raw_base_reference(tuple((1-X*X)*a for a in n),X,dps=dps,elementary=True)',
        'ref = elliptic_bundle_embedding(kind,data,dps=dps,elementary=True)',
        'phase = _elementary_inner(ref,result)',
        'x,y = data[:2]',
        'return kind,(mp.conj(phase)*x,phase*y)'])
    and ordered_nodes37(correct_tree30.body,[
        'result = elementary_projection_section(xi,X,N,dps,diagnostics)',
        'return result'])
    and matching_if_test32(slice_tree30,
        'distance<alpha/2 and not force_slab')
    and matching_if_test32(correct_tree30,'distance<alpha')
    and matching_if_test32(slice_tree30,
        'max_work is not None and cost>max_work')
    and matching_statement31(slice_tree30,
        "raise ValueError('certified finite expression exceeds max_work; inspect elementary_closed_parameters()')")
    and matching_if_test32(raw_tree31,'inf or abs(t)<mp.mpf(\'.5\')')
    and matching_if_test32(raw_tree31,
        "phase_coordinate.real<-mp.mpf('.9')*abs(phase_coordinate)")
    and matching_if_test32(bundle_tree32,
        "abs(x*y)>mp.mpf('.3') or abs(H)>mp.mpf('.6')")
    and matching_if_test32(bundle_tree32,'H<0 and w>0')
    and matching_assignment30(raw_tree31,'ratio',
        'u2/u1 if inverse_ratio else u1/u2')
    and matching_assignment30(raw_tree31,'small',
        'ratio*mp.sqrt((2*B if inverse_ratio else 2)/(b+root))')
    and matching_assignment30(raw_tree31,'phase_coordinate',
        '-mp.conj(u1) if inf else u2')
    and matching_assignment30(raw_tree31,'a',
        'mp.arg(phase_coordinate)/(2*mp.pi)')
    and matching_assignment30(bundle_tree32,'H',
        'mp.mpf(data[2]) if elementary else elliptic_height_chart(x,y,nu,dps=dps)')
    and matching_statement31(bundle_tree32,'block = [hw(H)*q for q in f]')
    and matching_statement31(bundle_tree32,
        "w = hw(H+mp.mpf('.5') if H<0 else H-mp.mpf('.5'))")
    and matching_statement31(bundle_tree32,'vec[start:start+4] = block')
    and len(finite_seed37)==1
    and exact_if_body37(finite_seed37[0],"v<=mp.mpf(1)/16",'value = v')
    and returned_expr30(finite_height_tree37,
        'seed(A)+seed(q*A)-seed(Y)-seed(q*Y)+q*q*(A-Y)/(1-q)')
    and len(restore_coeffs37)==1 and len(restore_cusp_data37)==1
    and ordered_nodes37(restore_coeffs37[0].body,[
        'det = e2.real*e3.imag-e3.real*e2.imag',
        'return ((rhs.real*e3.imag-e3.real*rhs.imag)/det, (e2.real*rhs.imag-rhs.real*e2.imag)/det)'])
    and returned_expr30(restore_cusp_data37[0],
        "((t*t-1,mp.j*(1+t*t),-2*t,t,mp.j*t) if inf else "
        "(1-t*t,mp.j*(1+t*t),2*t,t,-mp.j*t))")
    and exact_if_body37(restore_tree37,"kind in ('c0','ci')",'''
x,y = data[:2]
inf = kind=='ci'
t = x*y
e2,e3,e4,a2,a3 = cusp_data(t,inf)
c = coeffs(e4,e2,e3)
H = scalar_finite_height(x,y,mp.exp(-2) if inf else mp.mpf(1))
C = tuple(v+H*w for v,w in zip(phi,c))
A = C[0]*a2+C[1]*a3
return kind,(x*mp.exp(-A),y*mp.exp(A),C[0]*e2+C[1]*e3)
'''))

# A strict, nonempty open chamber is specified by X<0, X>start,
# 1+X>alpha, lambda(X)<1/10000, |n0|<1/100,
# n1<-3/4 and 1/2<n2<2/3.  At the witness n=(0,-4/5,3/5).
# lambda(h) is strictly increasing, so its upper bound holds for all
# h_j=start+(X-start)j/N, including the actual huge fixed N.
Zw37 = (sp.Rational(1,100),sp.I/sp.Integer(200))
Xw37 = -sp.sqrt(sp.Rational(7999,8000))
xi_w37 = tuple(z/sp.sqrt(1-Xw37**2) for z in Zw37)
n_w37 = (sp.simplify(2*sp.re(xi_w37[0]*sp.conjugate(xi_w37[1]))),
         sp.simplify(2*sp.im(xi_w37[0]*sp.conjugate(xi_w37[1]))),
         sp.simplify(abs(xi_w37[0])**2-abs(xi_w37[1])**2))
start37 = -1+sp.Rational(1,10**9)
lam_end37 = 1/(8000+sp.sqrt(8000**2-1))
lam_expr37 = lambda36.subs(X36,Xw37)
strict_witness37 = (
    n_w37==(0,-sp.Rational(4,5),sp.Rational(3,5))
    and sp.simplify(sum(z*sp.conjugate(z) for z in xi_w37)-1)==0
    and sp.simplify(lam_expr37-lam_end37)==0
    and 0<lam_end37<sp.Rational(1,15999)<sp.Rational(1,10000)
    and start37<Xw37<0
    and sp.Rational(1,16000)<1+Xw37<sp.Rational(1,8000)
    and sp.Rational(1,16000)>sp.Rational(1,10**5)
    and abs(n_w37[0])<sp.Rational(1,100)
    and n_w37[1]<-sp.Rational(3,4)
    and sp.Rational(1,2)<n_w37[2]<sp.Rational(2,3)
    and radial_base36 and count_ok30)

# For any point in that chamber and every 0<lambda<1/10000,
# y=(lambda*n0,1-lambda*n2,-lambda*n1).  Its stereographic numerators
# are u1=lambda*((2*n2-lambda)+2i*n0),
# u2=2*(1-lambda*n2-i*lambda*n1), up to their common positive S.
# They have Re u1, Re u2 and Im u2 strictly positive.  The indicated
# rational bounds prove |u1/u2|<7lambda/10 and hence the source chooses
# its non-inverse south ratio and direct argument chart.  In particular
# b=Aprof(1-|ratio|^2)>121, so the profile factor is <1/11 and
# 0<|t|<lambda/15<1/150000<3/10.
lam37,n037,n137,n237 = sp.symbols(
    'lam37 n037 n137 n237',real=True)
u1num37 = lam37*(2*n237-lam37+2*sp.I*n037)
u2num37 = 2*(1-lam37*n237-sp.I*lam37*n137)
y37 = sp.Matrix([lam37*n037,1-lam37*n237,-lam37*n137])
S37 = 1+sum(v*v for v in y37)
raw_numerators37 = (
    sp.expand(1-sum(v*v for v in y37)+2*sp.I*y37[0]-u1num37
              -lam37**2*(1-n037**2-n137**2-n237**2))==0
    and sp.expand(2*(y37[1]+sp.I*y37[2])-u2num37)==0
    and sp.Rational(67,50)**2>
        (sp.Rational(4,3)**2+sp.Rational(1,50)**2)
    and sp.Rational(67,50)/(2*(1-sp.Rational(1,10000)))
        <sp.Rational(7,10)
    and (sp.Rational(7,10)/10000)**2<sp.Rational(1,10**8)
    and 144*(1-sp.Rational(1,10**8))>121
    and sp.Rational(7,110)<sp.Rational(1,15)
    and sp.Rational(1,150000)<sp.Rational(3,10))

# Put z_ang=(-lambda*n1)/(1-lambda*n2)>0, so
# a=atan(z_ang)/(2pi).  Since n1<-3/4, z_ang>3lambda/4;
# z_ang<1/1000.  The elementary integral bound
# z/(1+z^2)<=atan(z)<z, with 3<pi<22/7, gives
# lambda/10<a<lambda/5<1/50000<1/10.  This selects cell=0,
# the a>0 positive-x gauge, f=a, and both embedding weight plateaux.
eps37 = sp.Symbol('eps37',positive=True)
height_bounds37 = (
    sp.simplify(sp.diff(sp.atan(eps37)-eps37/(1+eps37**2),eps37)
                -2*eps37**2/(1+eps37**2)**2)==0
    and sp.simplify(sp.diff(eps37-sp.atan(eps37),eps37)
                    -eps37**2/(1+eps37**2))==0
    and sp.Rational(1,10000)/(1-sp.Rational(1,10000))
        <sp.Rational(1,1000)
    and sp.Rational(3,4)*sp.Rational(7,44)
        /(1+sp.Rational(1,10**6))>sp.Rational(1,10)
    and 1/(6*(1-sp.Rational(1,10000)))<sp.Rational(1,5)
    and sp.Rational(1,50000)<sp.Rational(1,10)
    and bool(3<sp.pi<sp.Rational(22,7))
    and sp.floor(sp.Rational(1,10)+sp.Rational(1,2))==0
    and sp.Rational(1,5)<sp.Rational(1,4)
    and (sp.Rational(2,5)-sp.Rational(1,4))/sp.Rational(3,20)==1)
rep('C37a ordered source c0 path and strict uniform all-j chamber',
    source_path37 and strict_witness37
    and raw_numerators37 and height_bounds37)

# In the selected source branch, write e=exp(-4s/(1-s)),
# q=e*s, Delta=sqrt(a^2+4q), A=(Delta+a)/2 and Y=(Delta-a)/2.
# Then x=sqrt(A/e)>0, y=t/x, A=e*x^2, Y=|y|^2, A-Y=a,
# and AY=q.  The ten-vector has only first-pair entries
# r=(x,y)/sqrt(x^2+Y).  The proved a>lambda/10 and
# sqrt(q)<|t|<lambda/15 imply q<4a^2/9, Delta<5a/3, and therefore
# its first-slot probability p=(A/e)/(A/e+Y)>4/5.
a37,q37 = sp.symbols('a37 q37',positive=True)
Delta37 = sp.sqrt(a37**2+4*q37)
A37,Y37 = (Delta37+a37)/2,(Delta37-a37)/2
eps_pos37,Fpos37 = sp.symbols('eps_pos37 Fpos37',positive=True)
e_manifest37 = eps_pos37/(eps_pos37+Fpos37)
p37 = (A37/e_manifest37)/(A37/e_manifest37+Y37)
pmin37 = sp.Rational(4,5)
source_frame37 = (
    sp.simplify(A37-Y37-a37)==0
    and sp.simplify(A37*Y37-q37)==0
    and sp.simplify(q37/A37-Y37)==0
    and sp.simplify((A37/e_manifest37-Y37)
        -((1-e_manifest37)*Delta37
          +(1+e_manifest37)*a37)/(2*e_manifest37))==0
    and sp.simplify(p37-A37/(A37+Y37)
        -A37*Y37*(1-e_manifest37)
         /((A37+e_manifest37*Y37)*(A37+Y37)))==0
    and sp.simplify(A37/(A37+Y37)
                    -(sp.Rational(1,2)+a37/(2*Delta37)))==0
    and (sp.Rational(1,15)/sp.Rational(1,10))**2
        ==sp.Rational(4,9)
    and sp.simplify(sp.sqrt(a37**2+4*sp.Rational(4,9)*a37**2)
                    -5*a37/3)==0
    and 1+4*sp.Rational(2,3)**2==sp.Rational(25,9)
    and (1+sp.Rational(3,5))/2==pmin37
    and pmin37>sp.Rational(7,10)
    and 2*pmin37-1==sp.Rational(3,5)
    and sp.Rational(1,10)>sp.Rational(1,15)
    and sp.Rational(1,10**10)<sp.Rational(1,8))
# The unit projection pair and both source plateaux are exact.  The
# source frame radius, including the final source norm2, is positive.
tr37,ti37,x37 = sp.symbols('tr37 ti37 x37',real=True,positive=True)
t37 = tr37+sp.I*ti37
yframe37 = t37/x37
frame37 = sp.Matrix([x37,yframe37])/sp.sqrt(
    x37**2+(tr37**2+ti37**2)/x37**2)
frame_unit37 = sp.simplify((frame37.conjugate().T*frame37)[0]-1)==0
Ppair37 = frame37*frame37.conjugate().T
frame_projector37 = (frame_unit37
    and sp.simplify(Ppair37*Ppair37-Ppair37)==sp.zeros(2)
    and sp.simplify(sp.trace(Ppair37)-1)==0
    and sp.Rational(1,5)<sp.Rational(1,4)
    and sp.Rational(2,5)==sp.Rational(1,4)+sp.Rational(3,20))

# For the upper Hopf reference, the first pole-inclusion probability
# is (1-n1)/2>7/8.  The reverse triangle inequality therefore gives
# |c0|>sqrt((4/5)(7/8))-sqrt((1/5)(1/8))
#      >4/5-1/5=3/5>2/5.
# Every adjacent overlap has Re(cj)>pmin-(1-pmin)=3/5>1/sqrt(5).
# These inequalities hold for arbitrary j and j-1, so the source
# normalizations exist for the actual N=4*10^8000+1 without unrolling it.
ur37,ui37,vr37,vi37 = sp.symbols(
    'ur37 ui37 vr37 vi37',real=True)
spin_a37,spin_b37 = ur37+sp.I*ui37,vr37+sp.I*vi37
pole_first37 = (spin_a37-sp.I*spin_b37)/sp.sqrt(2)
hopf_im37 = 2*sp.im(spin_a37*sp.conjugate(spin_b37))
pole_probability37 = sp.expand(
    pole_first37*sp.conjugate(pole_first37)
    -(spin_a37*sp.conjugate(spin_a37)
      +spin_b37*sp.conjugate(spin_b37)-hopf_im37)/2)==0
overlap_margin37 = (
    pole_probability37
    and sp.simplify((1-n_w37[1])/2-sp.Rational(9,10))==0
    and (1+sp.Rational(3,4))/2==sp.Rational(7,8)
    and sp.Rational(7,10)>sp.Rational(16,25)
    and sp.Rational(1,40)<sp.Rational(1,25)
    and sp.Rational(3,5)>sp.Rational(2,5)
    and sp.Rational(3,5)>1/sp.sqrt(5)
    and sp.Rational(7,10)<pmin37
    and count_ok30 and N30>0)
phase37,c_r37,c_i37 = sp.symbols('phase37 c_r37 c_i37',real=True)
step_c37 = c_r37+sp.I*c_i37
unit_phase37 = sp.exp(sp.I*phase37)
phase_induction37 = (
    sp.simplify(sp.conjugate(step_c37*unit_phase37)
                *(step_c37*unit_phase37)
                -sp.conjugate(step_c37)*step_c37)==0
    and sp.simplify(sp.conjugate(unit_phase37)*unit_phase37-1)==0
    and project_source30 and projection_value30 and projection_norm30)
rep('C37b actual first-pair projector and fixed-count overlap induction',
    source_path37 and source_frame37 and frame_projector37
    and overlap_margin37 and phase_induction37)

# Type (d) only: one source evaluation of the initial overlap and one
# N=1 endpoint overlap at the exact witness.  The second number is not
# the minimum over the certified N and is not used in the uniform proof.
import mpmath as mp37
with mp37.workdps(80):
    Xsample37 = -mp37.sqrt(mp37.mpf(7999)/8000)
    Zsample37 = (mp37.mpc(mp37.mpf(1)/100),mp37.j/200)
    xisample37 = tuple(z/mp37.sqrt(1-Xsample37**2)
                       for z in Zsample37)
    pcsample37 = 2*xisample37[0]*mp37.conj(xisample37[1])
    nsample37 = (pcsample37.real,pcsample37.imag,
                 abs(xisample37[0])**2-abs(xisample37[1])**2)
    def sample_reference37(h):
        kind,data = cusp30.elliptic_raw_base_reference(
            tuple((1-h*h)*a for a in nsample37),h,
            dps=80,elementary=True)
        return kind,cusp30.elliptic_bundle_embedding(
            kind,data,dps=80,elementary=True)
    kind_start37,ref_start37 = sample_reference37(
        -1+mp37.mpf(1)/10**9)
    kind_end37,ref_end37 = sample_reference37(Xsample37)
    measured_initial37 = abs(cusp30._elementary_inner(
        ref_start37,cusp30._elementary_pole_inclusion(xisample37)))
    measured_endpoint37 = abs(cusp30._elementary_inner(
        ref_end37,ref_start37))
    measured_initial_margin37 = measured_initial37-mp37.mpf(3)/5
    measured_endpoint_margin37 = measured_endpoint37-mp37.mpf(3)/5
sample_overlap37 = (
    kind_start37==kind_end37=='c0'
    and measured_initial37>mp37.mpf(3)/5
    and measured_endpoint37>mp37.mpf(3)/5)
rep('C37b2 initial and N=1 endpoint overlap samples (numerical, d)',
    sample_overlap37)
print('C37b2 initial/N=1 overlaps and margins to 3/5:',
      mp37.nstr(measured_initial37,24),
      mp37.nstr(measured_endpoint37,24),
      mp37.nstr(measured_initial_margin37,24),
      mp37.nstr(measured_endpoint_margin37,24))

# On this same chamber the source four-seed height is completely linear.
# Since a<lambda/5, sqrt(q)<|t|<lambda/15, we have
# A<a+sqrt(q)<4lambda/15 and Y<sqrt(q)<lambda/15.
# Also e>1/2 by exp(-u)>1-u and 4s/(1-s)<1/2.
# Thus A,Y,qA,qY<1/16, and source H_f=(A-Y)/(1-q)=a/(1-q).
seed_bounds37 = (
    sp.simplify(sp.diff(sp.exp(-eps37)-(1-eps37),eps37)
                -(1-sp.exp(-eps37)))==0
    and sp.Rational(1,150000)**2<sp.Rational(1,10**10)
    and 4*sp.Rational(1,150000)**2
      /(1-sp.Rational(1,150000)**2)<sp.Rational(1,2)
    and sp.Rational(4,15)/10000<sp.Rational(1,16)
    and sp.Rational(1,15)/10000<sp.Rational(1,16)
    and sp.Rational(1,10**10)<1
    and sp.Rational(1,150000)**2<sp.Rational(1,4))
height_seed37 = (A37+q37*A37-Y37-q37*Y37
                 +q37**2*(A37-Y37)/(1-q37))
finite_height37 = (
    sp.simplify(height_seed37-(A37-Y37)/(1-q37))==0
    and sp.simplify(height_seed37-a37/(1-q37))==0
    and source_path37 and seed_bounds37)
rep('C37c source finite cusp height collapses to a/(1-q)',
    finite_height37)

# The signed determinant is a product of genuine source stages, in input
# order (n0,n1,X,theta,phi2,phi3) and holomorphic output order
# (Re z1,Im z1,Re z2,Im z2,Re z3,Im z3).  Put
# v=|y1+i*y2|, S=1+|y|^2, s=|t|^2.  In polar y1+i*y2=v*e^(i alpha),
# ratio*e^(i alpha)=(1-u^2-v^2)/(2v)+i*u/v and a=alpha/(2pi).
# Thus det D(Re ratio,Im ratio,a)/D(y0,y1,y2)=S/(4pi v^4).
n0_det37,n1_det37 = sp.symbols('n0_det37 n1_det37',real=True)
n2_det37 = sp.sqrt(1-n0_det37**2-n1_det37**2)
dy_positive37 = sp.Matrix([
    [L36,0,Lp36*n0_det37],
    [-L36*sp.diff(n2_det37,n0_det37),
     -L36*sp.diff(n2_det37,n1_det37),-Lp36*n2_det37],
    [0,-L36,-Lp36*n1_det37]])
round_base_det37 = sp.simplify(dy_positive37.det()
                                +L36**2*Lp36/n2_det37)==0
u37,v37 = sp.symbols('u37 v37',real=True,positive=True)
Rratio37 = (1-u37**2-v37**2)/(2*v37)
Iratio37 = u37/v37
S_polar37 = 1+u37**2+v37**2
ratio_det_polar37 = sp.Matrix([
    [sp.diff(Rratio37,u37),sp.diff(Rratio37,v37)],
    [sp.diff(Iratio37,u37),sp.diff(Iratio37,v37)]]).det()
ratio_det37 = (sp.simplify(ratio_det_polar37
                        -S_polar37/(2*v37**3))==0
               and sp.simplify(ratio_det_polar37/(2*sp.pi*v37)
                               -S_polar37/(4*sp.pi*v37**4))==0)
# The source profile has z=(B*s^2+A*s)/(1+A*s), so its inverse radial
# map has positive real determinant 1/z_s'(s), already checked in C36.
# At fixed t, rho=log(Y)/2 has rho_a=-1/(2Delta).
rho37 = sp.log(Y37)/2
rho_derivative37 = sp.simplify(sp.diff(rho37,a37)
                               +1/(2*Delta37))==0
# Write the source cusp fibre phase as sigma=e^(i psi).  The real map
# (r,beta,rho,psi)->(x',y') with
# x'=r e^(-rho-i psi), y'=e^(rho+i(beta+psi)) has determinant r;
# (r,beta)->(Re t,Im t) also has determinant r.  Hence its fixed-t
# Jacobian from (Re t,Im t,rho,psi) is exactly +1.
rpol37,beta37,rho_pol37,psi37 = sp.symbols(
    'rpol37 beta37 rho_pol37 psi37',real=True)
xpol37 = rpol37*sp.exp(-rho_pol37-sp.I*psi37)
ypol37 = sp.exp(rho_pol37+sp.I*(beta37+psi37))
cusp_polar37 = sp.Matrix([
    sp.re(xpol37),sp.im(xpol37),
    sp.re(ypol37),sp.im(ypol37)])
cusp_det37 = sp.trigsimp(sp.factor(
    cusp_polar37.jacobian(
        [rpol37,beta37,rho_pol37,psi37]).det()))
# The source C=phi+H_f*c is a shear with phi block I2.  At fixed C,
# (x,y)->(x exp(-kxy),y exp(kxy)) has complex determinant one; the
# final real torus block has det 1-|t|^4>0, as in C9e/C36.
xflow37,yflow37,kflow37 = sp.symbols('xflow37 yflow37 kflow37')
flow37 = sp.Matrix([
    xflow37*sp.exp(-kflow37*xflow37*yflow37),
    yflow37*sp.exp(kflow37*xflow37*yflow37)])
flow_det37 = sp.simplify(flow37.jacobian(
    [xflow37,yflow37]).det()-1)==0
# The actual source south cusp torus basis is the 2x2 real matrix R(t).
tr_det37,ti_det37 = sp.symbols('tr_det37 ti_det37',real=True)
tt_det37 = tr_det37+sp.I*ti_det37
e2_det37,e3_det37 = 1-tt_det37**2,sp.I*(1+tt_det37**2)
Rtorus37 = sp.Matrix([[sp.re(e2_det37),sp.re(e3_det37)],
                      [sp.im(e2_det37),sp.im(e3_det37)]])
torus_block37 = sp.factor(Rtorus37.det()
    -(1-(tr_det37**2+ti_det37**2)**2))==0
lampos37,lamprime37,npos37,Spos37,vpos37,zprime37,Dpos37,spos37 = sp.symbols(
    'lampos37 lamprime37 npos37 Spos37 vpos37 zprime37 Dpos37 spos37',
    positive=True)
det_factors37 = (-lampos37**2*lamprime37/npos37
    *Spos37/(4*sp.pi*vpos37**4)
    /zprime37 *(-1/(2*Dpos37))*(2*sp.pi)
    *(1-spos37**2))
det_closed37 = (lampos37**2*lamprime37*Spos37*(1-spos37**2)
    /(4*npos37*vpos37**4*zprime37*Dpos37))
det_certificate37 = (
    round_base_det37 and ratio_det37 and radial_base36 and radial_profile36
    and rho_derivative37 and cusp_det37==rpol37
    and flow_det37 and torus_block37
    and sp.simplify(det_factors37-det_closed37)==0
    and sp.Rational(1,150000)**2<1
    and source_path37 and overlap_margin37 and finite_height37)
rep('C37d actual near-south free-slab six-real determinant is positive',
    det_certificate37)

# C38: the artificial argument-chart selector in the raw reference has an
# exact OPEN overlap in both pure cusp charts.  Its lower negative-real
# sector is the nontrivial case: the two height values differ by one deck
# unit, yet the normalized ten-slot reference, not just its projector,
# agrees.  This is a parameter-wise identity and applies at any meridian
# index for which this source branch is active; it does not assert that all
# other meridian branches or seams have been covered.
arg_selector38 = [node for node in ast.walk(raw_tree31)
                  if isinstance(node,ast.If)
                  and ast.dump(node.test,include_attributes=False)
                  ==ast.dump(ast.parse(
                      "phase_coordinate.real<-mp.mpf('.9')*abs(phase_coordinate)",
                      mode='eval').body,include_attributes=False)]
core_sign38 = [node for node in ast.walk(core_tree37)
               if isinstance(node,ast.If)
               and ast.dump(node.test,include_attributes=False)
               ==ast.dump(ast.parse('a>0',mode='eval').body,
                           include_attributes=False)]
hw_nodes38 = [node for node in ast.walk(bundle_tree32)
              if isinstance(node,ast.FunctionDef) and node.name=='hw']
frame_nodes38 = [node for node in ast.walk(bundle_tree32)
                 if isinstance(node,ast.FunctionDef) and node.name=='frame']
seed_nodes38 = [node for node in ast.walk(finite_height_tree37)
                if isinstance(node,ast.FunctionDef) and node.name=='seed']
flat_tree38 = source_tree30(cusp30.scalar_flat_step)
arg_source38 = (
    source_path37 and len(arg_selector38)==1
    and [ast.dump(node,include_attributes=False)
         for node in arg_selector38[0].body]
        ==[ast.dump(node,include_attributes=False) for node in ast.parse('''
a = mp.mpf('.5')+mp.arg(-phase_coordinate)/(2*mp.pi)
''').body]
    and [ast.dump(node,include_attributes=False)
         for node in arg_selector38[0].orelse]
        ==[ast.dump(node,include_attributes=False) for node in ast.parse('''
a = mp.arg(phase_coordinate)/(2*mp.pi)
''').body]
    and len(core_sign38)==1
    and ordered_nodes37(core_tree37.body,[
        '''Y = 2/(R+f) if cell==-1 else (R-f)/2''',
        'y = mp.sqrt(Y)',
        '''
if y==0:
    raise ValueError('the core reference excludes the node')
''',
        'return t/y,y'])
    and exact_if_body37(core_tree37,'a>0','''
A = 2/(R-f) if cell==1 else (R+f)/2
x = mp.sqrt(A/e)/nu
return x,t/x
''')
    and ordered_nodes37(core_tree37.body,[
        'f = inverse_height_moment(b,q,0)',
        'R = mp.sqrt(f*f+4*q)'])
    and [ast.dump(node,include_attributes=False)
         for node in core_sign38[0].body]
        ==[ast.dump(node,include_attributes=False) for node in ast.parse('''
A = 2/(R-f) if cell==1 else (R+f)/2
x = mp.sqrt(A/e)/nu
return x,t/x
''').body]
    and ordered_nodes37(core_tree37.body,[
        '''
Y = 2/(R+f) if cell==-1 else (R-f)/2
''',
        'return t/y,y'])
    and exact_if_body37(moment_tree37,"aa>=mp.mpf('.4')",
        'return mp.sign(a)*(B*mp.exp(kappa*aa)-C*mp.exp(-kappa*aa))')
    and len(hw_nodes38)==len(frame_nodes38)==len(seed_nodes38)==1
    and returned_expr30(hw_nodes38[0],
        "1-scalar_flat_step((abs(h)-mp.mpf('.25'))/mp.mpf('.15'))")
    and returned_expr30(frame_nodes38[0],
        '(nu*mp.conj(x)/mp.sqrt(r),y/mp.sqrt(r))')
    and matching_assignment30(frame_nodes38[0],'r',
        'nu*nu*abs(x)**2+abs(y)**2')
    and matching_assignment30(seed_nodes38[0],'reciprocal','u>2')
    and exact_if_body37(seed_nodes38[0],
        "v<=mp.mpf(1)/16",'value = v')
    and exact_if_body37(seed_nodes38[0],
        "v>=mp.mpf('.5')",'value = v/(1+v)')
    and returned_expr30(seed_nodes38[0],
        '1-value if reciprocal else value')
    and exact_if_body37(flat_tree38,'u<=0','return mp.mpf(0)')
    and exact_if_body37(flat_tree38,'u>=1','return mp.mpf(1)')
    and matching_assignment30(raw_tree31,'w',
        '(vv[0]/(1-aa),-vv[2]/(1-aa),-vv[1]/(1-aa))')
    and matching_assignment30(raw_tree31,'y',
        '(w[0],w[1]+1,w[2])')
    and matching_assignment30(raw_tree31,'u1',
        'mp.mpc(1-s2,2*y[0])/(1+s2)')
    and matching_assignment30(raw_tree31,'u2',
        'mp.mpc(2*y[1],2*y[2])/(1+s2)')
    and matching_assignment30(raw_tree31,'inverse_ratio',
        'abs(u1)>abs(u2)')
    and matching_assignment30(raw_tree31,'inf',
        "inverse_ratio and abs(small)<mp.mpf('.5')")
    and matching_assignment30(raw_tree31,'ratio',
        'u2/u1 if inverse_ratio else u1/u2')
    and matching_assignment30(raw_tree31,'phase_coordinate',
        '-mp.conj(u1) if inf else u2')
    and matching_assignment30(raw_tree31,'nu',
        "mp.exp(-2) if inf else mp.mpf(1)")
    and returned_expr30(bundle_tree32,
        '(result,+norm2) if with_norm else result'))

# On the lower wedge -pi<arg(p)<-5pi/6, the principal chart has
# a0 in (-1/2,-5/12), whereas the half-centred chart has a1=a0+1
# in (1/2,7/12).  Both charts are valid: 7/12<3/5.  The source
# selector cos(arg p)=-9/10 lies strictly inside this open wedge.
# The upper wedge 5pi/6<arg(p)<pi has a1=a0 instead.
angle38 = sp.Symbol('angle38',real=True)
open_wedges38 = (
    sp.simplify(sp.Rational(1,2)
                +(angle38+sp.pi)/(2*sp.pi)
                -angle38/(2*sp.pi)-1)==0
    and sp.simplify(sp.Rational(1,2)
                    +(angle38-sp.pi)/(2*sp.pi)
                    -angle38/(2*sp.pi))==0
    and sp.Rational(9,10)**2>sp.Rational(3,4)
    and sp.Rational(5,12)>sp.Rational(2,5)
    and sp.Rational(7,12)<sp.Rational(3,5)
    and sp.Rational(1,12)<sp.Rational(1,4)
    and sp.floor(-sp.Rational(5,12)+sp.Rational(1,2))==0
    and sp.floor(sp.Rational(7,12)+sp.Rational(1,2))==1)

# Exact source-round witnesses at the selector are obtained from
# p=(-9-i sqrt(19))/10 and ratio=1/20.  For c0 take u2=p/sqrt(1+r^2),
# u1=r*u2; for ci take u1=-conj(p)/sqrt(1+r^2), u2=r*u1.
# Their conjugates witness the upper wedge.  Inverse stereographic
# coordinates w=(Im u1/d, Re u2/d-1, Im u2/d), d=1+Re u1, give
# R^2=|w|^2, X=(R-1)/sqrt(1+R^2), n=(w0,-w2,-w1)/R.
# Then v=(1-X^2)n and the production raw map recovers u1,u2.
rw38 = sp.Rational(1,20)
raw_witness38 = True
for sign_im38 in (-1,1):
    for perturb38 in (-sp.Rational(1,1000),sp.Integer(0),
                      sp.Rational(1,1000)):
        pr38 = -sp.Rational(9,10)+perturb38
        pw38 = pr38+sign_im38*sp.I*sp.sqrt(1-pr38**2)
        for inf38 in (False,True):
            u238 = pw38/sp.sqrt(1+rw38**2) if not inf38 else None
            u138 = rw38*u238 if not inf38 else (
                -sp.conjugate(pw38)/sp.sqrt(1+rw38**2))
            u238 = u238 if not inf38 else rw38*u138
            phase_w38 = -sp.conjugate(u138) if inf38 else u238
            d_w38 = 1+sp.re(u138)
            w_w38 = (sp.im(u138)/d_w38,sp.re(u238)/d_w38-1,
                     sp.im(u238)/d_w38)
            R2_w38 = sp.simplify(sum(v*v for v in w_w38))
            y_w38 = (w_w38[0],w_w38[1]+1,w_w38[2])
            Sy_w38 = 1+sum(v*v for v in y_w38)
            u1_recovered38 = (
                1-sum(v*v for v in y_w38)+2*sp.I*y_w38[0])/Sy_w38
            u2_recovered38 = (
                2*(y_w38[1]+sp.I*y_w38[2]))/Sy_w38
            phase_abs38 = sp.sqrt(sp.simplify(
                phase_w38*sp.conjugate(phase_w38)))
            selector_side38 = bool(
                sp.simplify(sp.re(phase_w38)
                            +sp.Rational(9,10)*phase_abs38)<0)
            raw_witness38 = raw_witness38 and (
                sp.simplify(u138*sp.conjugate(u138)
                            +u238*sp.conjugate(u238)-1)==0
                and sp.simplify(u1_recovered38-u138)==0
                and sp.simplify(u2_recovered38-u238)==0
                and sp.re(pw38)<-sp.sqrt(3)/2
                and bool(sp.im(phase_w38)>0)==(sign_im38==1)
                and selector_side38==(perturb38<0)
                and d_w38>0
                and w_w38[1]<0
                and sp.Rational(1,4)<R2_w38<9
                and (sp.simplify(u238/u138-rw38)==0 if inf38 else
                     sp.simplify(u138/u238-rw38)==0))
# At perturb=0 the selector equality is exact; the two rational
# perturbations lie on its opposite sides, while all other strict
# cusp/slab inequalities persist.  Continuity supplies actual open
# source neighborhoods on both sides of the selector.
w038,w138,w238 = sp.symbols('w038 w138 w238',real=True)
R38_w = sp.Symbol('R38_w',positive=True)
X38_w = (R38_w-1)/sp.sqrt(1+R38_w**2)
n38_w = sp.Matrix((w038,-w238,-w138))/R38_w
v38_w = (1-X38_w**2)*n38_w
a38_w = (R38_w**2-1)/(1+R38_w**2)
round_inverse38 = (
    sp.simplify((1-X38_w**2)-2*R38_w/(1+R38_w**2))==0
    and sp.simplify((2-X38_w**2)
                    -(R38_w+1)**2/(1+R38_w**2))==0
    and R38_w+1>0
    and sp.simplify(v38_w[0]/(1-a38_w)-w038)==0
    and sp.simplify(-v38_w[2]/(1-a38_w)-w138)==0
    and sp.simplify(-v38_w[1]/(1-a38_w)-w238)==0
    and sp.simplify((1-X38_w**2)/(1-a38_w)-R38_w)==0
    and sp.simplify(sum(n38_w[i]**2 for i in range(3))
                    -(w038**2+w138**2+w238**2)/R38_w**2)==0)
# The profile-root denominator exceeds 2b, where
# b=A(1-r^2)>143>B; the earlier C34a series bound gives
# B=e^4<64.  Hence |tc|<r in ci and
# |tc|<r/sqrt(143) in c0.  In both cases |tc|<1/10.
witness_domain38 = (
    raw_witness38 and round_inverse38
    and sp.Rational(144)*(
        1-sp.Rational(1,400))>143
    and sp.exp(4)<64<143
    and rw38<sp.Rational(1,10)
    and rw38/sp.sqrt(143)<sp.Rational(1,10)
    and 2/sp.sqrt(10)<sp.Rational(2,3)
    and 1/sp.sqrt(5)<sp.Rational(2,3)
    and sp.Rational(1,3)>sp.Rational(1,10**5))
rep('C38a both raw cusp argument-selector seams have exact round witnesses',
    arg_source38 and open_wedges38 and witness_domain38)

# In the lower wedge the source outer-core formula has
# b=a0 in (-1/2,-5/12), E=exp(-2b-1) in (e^-1/6,1),
# q=e*nu^2*|tc|^2<1/100 and f=-E+q/E.  The positive
# discriminant is R=E+q/E, so the negative gauge has
# (x0,y0)=(tc/sqrt(E),sqrt(E)), whereas the shifted positive
# gauge has (x1,y1)=(1/(nu sqrt(e) sqrt(E)),tc/x1).
# These are exactly the source deck coordinates psi_hat(x0,y0).
E38,q38,eps38,nu38 = sp.symbols(
    'E38 q38 eps38 nu38',positive=True)
tr38,ti38 = sp.symbols('tr38 ti38',real=True)
t38 = tr38+sp.I*ti38
f38 = -E38+q38/E38
R38 = E38+q38/E38
x038,y038 = t38/sp.sqrt(E38),sp.sqrt(E38)
m38 = nu38*sp.sqrt(eps38)
x1deck38,y1deck38 = 1/(m38*y038),m38*x038*y038**2
x1core38 = sp.sqrt(1/(E38*eps38))/nu38
y1core38 = t38/x1core38
core_deck38 = (
    sp.simplify(R38**2-f38**2-4*q38)==0
    and sp.simplify((R38-f38)/2-E38)==0
    and sp.simplify(2/(R38-f38)-1/E38)==0
    and sp.simplify(x1deck38-x1core38)==0
    and sp.simplify(y1deck38-y1core38)==0
    and sp.exp(-sp.Rational(1,6))>sp.Rational(5,6)
    and sp.Rational(1,100)/(sp.Rational(5,6))
        <sp.Rational(1,16)
    and sp.Rational(1,100)**2/(sp.Rational(5,6))
        <sp.Rational(1,16))

# For both representatives the first height weight is zero and the
# second is one: |a0|,a1>5/12>2/5, while
# |a0+1/2|=|a1-1/2|<1/12<1/4.  The negative branch applies
# psi_hat before its second frame.  Thus the two ten-vectors have
# exactly the same sole active pair, in slots 2:4 (c0) or 6:8 (ci).
def frame38(x,y):
    radius = nu38**2*x*sp.conjugate(x)+y*sp.conjugate(y)
    return sp.Matrix((nu38*sp.conjugate(x)/sp.sqrt(radius),
                      y/sp.sqrt(radius)))


left_frame38 = frame38(x1deck38,y1deck38)
right_frame38 = frame38(x1core38,y1core38)
frame_diff38 = left_frame38-right_frame38
Pleft38 = left_frame38*left_frame38.conjugate().T
Pright38 = right_frame38*right_frame38.conjugate().T
jet_vars38 = (tr38,ti38,E38,nu38,eps38)
twojets38 = (
    all(sp.simplify(entry)==0 for entry in frame_diff38)
    and all(sp.simplify(sp.diff(entry,v))==0
            for entry in frame_diff38 for v in jet_vars38)
    and all(sp.simplify(sp.diff(entry,v,w))==0
            for entry in frame_diff38 for v in jet_vars38
            for w in jet_vars38)
    and all(sp.simplify(entry)==0 for entry in Pleft38-Pright38)
    and all(sp.simplify(sp.diff(entry,v))==0
            for entry in Pleft38-Pright38 for v in jet_vars38)
    and all(sp.simplify(sp.diff(entry,v,w))==0
            for entry in Pleft38-Pright38 for v in jet_vars38
            for w in jet_vars38))
rep('C38b lower selector references, projectors and two-jets agree openly',
    arg_source38 and core_deck38 and open_wedges38 and twojets38)

# The four-seed height has its exact deck shift on this open wedge.
# Original (A,Y)=(q/E,E); shifted (A,Y)=(1/E,qE).
# Here E,1/E are in the v>=1/2 seed branch, while q/E,
# q^2/E, qE, q^2 E lie in the v<=1/16 branch.  Both
# representatives satisfy the source assertion A,Y<4, q<1/8.
Hleft38 = (q38/E38+q38**2/E38-E38/(1+E38)-q38*E38
           +q38**2*(q38/E38-E38)/(1-q38))
Hright38 = (1/(1+E38)+q38/E38-q38*E38-q38**2*E38
            +q38**2*(1/E38-q38*E38)/(1-q38))
height_deck38 = (
    sp.factor(Hright38-Hleft38-1)==0
    and sp.Rational(5,6)>sp.Rational(1,2)
    and sp.Rational(6,5)<2
    and sp.Rational(6,5)<4
    and sp.Rational(1,100)<sp.Rational(1,8)
    and sp.Rational(1,100)<sp.Rational(1,4)
    and source_path37 and arg_source38)
rep('C38c source four-seed height shifts by one on the lower overlap',
    height_deck38)

# On either cusp, write sigma=+1 for c0 and -1 for ci.
# The real Cramer coefficients of e4 in the e2,e3 basis are
# c2=2 Re(t)/(1-|t|^2), c3=sigma*2 Im(t)/(1-|t|^2).
# With a2=t, a3=-sigma*i*t one has c2*a2+c3*a3=d
# =2|t|^2/(1-|t|^2), and m=nu*exp(-d).  The deck action
# (x,y,C)->(1/(my),mxy^2,C+c) therefore restores to the
# holomorphic transition
# T(z1,z2,z3)=(1/(nu*z2),nu*z1*z2^2,z3+sigma*2*z1*z2).
# Its complex determinant is one and its displayed inverse is
# holomorphic wherever z2 is nonzero, as it is on this overlap.
s38 = tr38**2+ti38**2
z138,z238,z338,sigma38 = sp.symbols(
    'z138 z238 z338 sigma38')
deck_T38 = sp.Matrix((1/(nu38*z238),
                      nu38*z138*z238**2,
                      z338+2*sigma38*z138*z238))
deck_T_inv38 = sp.Matrix((nu38*z138**2*z238,
                          1/(nu38*z138),
                          z338-2*sigma38*z138*z238))
holomorphic_deck38 = (
    sp.simplify(deck_T38.jacobian((z138,z238,z338)).det()-1)==0
    and sp.simplify(deck_T_inv38.subs(
        {z138:deck_T38[0],z238:deck_T38[1],
         z338:deck_T38[2]},simultaneous=True)
        -sp.Matrix((z138,z238,z338)))==sp.zeros(3,1))
for sign38 in (1,-1):
    e238 = sign38*(1-t38**2)
    e338 = sp.I*(1+t38**2)
    e438 = sign38*2*t38
    c238 = 2*tr38/(1-s38)
    c338 = sign38*2*ti38/(1-s38)
    a338 = -sign38*sp.I*t38
    holomorphic_deck38 = holomorphic_deck38 and (
        sp.simplify(c238*e238+c338*e338-e438)==0
        and sp.simplify(c238*t38+c338*a338
                        -2*s38/(1-s38))==0)
rep('C38d both restored cusp deck transitions are holomorphic',
    arg_source38 and height_deck38 and holomorphic_deck38)

# Generic normalized-projector recurrence.  Since P and its first two
# jets agree on the open overlap and all projection denominators are
# nonzero at the certified count, the update P_j v/||P_j v||
# has the same first two jets on both sides, by induction in j.
# The final cusp representatives are deck related because the
# phase correction commutes with psi_hat.  Hence the restored
# holomorphic charts have the displayed T transition and their
# pulled-back J agrees through the selector.  This checks this
# selector family only; no other branch or seam is inferred.
ang38 = sp.symbols('ang38',real=True)
unit38 = sp.exp(sp.I*ang38)
phase_deck38 = (
    sp.simplify(1/(m38*(unit38*y038))
                -sp.conjugate(unit38)*x1deck38)==0
    and sp.simplify(m38*(sp.conjugate(unit38)*x038)
                    *(unit38*y038)**2-unit38*y1deck38)==0
    and project_source30 and projection_norm30
    and source_path37 and twojets38)
rep('C38e finite projector recurrence respects this selector at every index',
    phase_deck38)

# C39: the dominant-coordinate selector in elliptic_raw_base_reference is
# an exact OPEN overlap, not merely a pointwise seam equality.  Put
# R=u1/u2, r=|R|^2, A=(4*pi)^2, B=exp(4), b=A(1-r),
# D=sqrt(b^2+4*B*r)>|b|.  The direct branch gives
# t=R sqrt(2/(b+D)).  The inverse branch uses S=1/R,
# b_I=-b/r, D_I=D/r and gives
# t=1/[S sqrt(2*B/(b_I+D_I))].  Its positive radial gain has
# square (D-b)/(2*B*r)=2/(D+b).  Positivity therefore identifies
# the two *positive* gains, including their first two derivatives.
selector_source39 = (
    middle_source31 and source_collar34 and source_torus35
    and project_source30 and projection_norm30 and source_cusp_ok30
    and matching_statement31(raw_tree31,
                             'A,B = (4*mp.pi)**2,mp.exp(4)')
    and matching_assignment30(raw_tree31,'rr','abs(ratio)**2')
    and matching_assignment30(raw_tree31,'b','A*(1-rr)')
    and matching_assignment30(raw_tree31,'root',
        'mp.sqrt(b*b+4*B*rr)')
    and matching_assignment30(raw_tree31,'small',
        'ratio*mp.sqrt((2*B if inverse_ratio else 2)/(b+root))')
    and matching_assignment30(raw_tree31,'inf',
        "inverse_ratio and abs(small)<mp.mpf('.5')")
    and matching_assignment30(raw_tree31,'t',
        'None if inf else (1/small if inverse_ratio else small)')
    and matching_if_test32(raw_tree31,"inf or abs(t)<mp.mpf('.5')")
    and matching_assignment30(raw_tree31,'tc','small if inf else t')
    and matching_assignment30(raw_tree31,'phase_coordinate',
        '-mp.conj(u1) if inf else u2')
    and len(arg_selector38)==1
    and [ast.dump(node,include_attributes=False)
         for node in arg_selector38[0].orelse]
        ==[ast.dump(node,include_attributes=False) for node in ast.parse(
            'a = mp.arg(phase_coordinate)/(2*mp.pi)').body]
    and matching_assignment30(raw_tree31,'a',
        'mp.arg(phase_coordinate)/(2*mp.pi)')
    and matching_assignment30(raw_tree31,'nu',
        "mp.exp(-2) if inf else mp.mpf(1)")
    and matching_if_test32(raw_tree31,
        "abs(tc)<=mp.mpf('.3') and elementary")
    and matching_assignment30(raw_tree31,'rho',
        '(elementary_cusp_logradius if elementary else '
        'elliptic_cusp_logradius)(a,tc,nu,dps=dps)')
    and matching_statement31(raw_tree31,
        'lattice = scalar_lattice_coords(zm,wm,t,lt)')
    and matching_assignment30(raw_tree31,'data',
        '(t,lattice[0],lattice[3],lt,a) if elementary else '
        '(t,lattice[0],lattice[3],lt)')
    and returned_expr30(raw_tree31,"'mid',data")
    and matching_assignment30(bundle_tree32,'known_height',
        'mp.mpf(data[4]) if elementary and len(data)>4 else None')
    and matching_assignment30(bundle_tree32,'norm2',
        'sum(abs(q)**2 for q in vec)')
    and matching_assignment30(bundle_tree32,'result',
        'tuple(+(q/mp.sqrt(norm2)) for q in vec)')
    and matching_call30(meridian_tree31,
        'mp.matrix([[a*mp.conj(b) for b in e] for a in e])')
    and returned_expr30(restore_tree30,
        "'mid',(z-mp.j*t*C[0]-t*C[1],"
        "w+(1-t*t)*C[0]+mp.j*(1+t*t)*C[1],t)"))

r39,A39,B39 = sp.symbols('r39 A39 B39',positive=True)
b39 = A39*(1-r39)
D39 = sp.sqrt(b39**2+4*B39*r39)
bI39,D_I39 = -b39/r39,D39/r39
Fdirect39 = 2/(b39+D39)
Finverse39 = (bI39+D_I39)/(2*B39)
sfrom39 = r39*Fdirect39
# Finverse is the square of 1/sqrt(2B/(b_I+D_I)), since S=1/R.
root_identity39 = (
    sp.simplify(Finverse39-(D39-b39)/(2*B39*r39))==0
    and sp.simplify((D39-b39)*(D39+b39)-4*B39*r39)==0
    and sp.simplify(B39*sfrom39**2+b39*sfrom39-r39)==0
    and sp.simplify(sfrom39-(D39-b39)/(2*B39))==0
    and all(sp.simplify(sp.diff(Fdirect39,r39,k)
                        -sp.diff(Finverse39,r39,k))==0
            for k in range(3))
    and sp.simplify(D39.subs(r39,1)-2*sp.sqrt(B39))==0
    and sp.simplify(Fdirect39.subs(r39,1)-1/sp.sqrt(B39))==0)
# Both radicands and both gains are positive because D>|b|, r,B>0.
# The following derivative recurrence is checked against a general
# smooth radial F, then the three identities above specialize it to
# *both* source gains.  R=x+i*y, r=x^2+y^2 are local real coordinates;
# (R,u2)->(u1,u2) is a diffeomorphism for u2!=0.
x39,y39 = sp.symbols('x39 y39',real=True)
rs39 = x39**2+y39**2
Ffun39 = sp.Function('Ffun39')
Fval39 = Ffun39(rs39)
g39 = sp.sqrt(Fval39)
tformal39 = (x39+sp.I*y39)*g39
F1_39 = sp.diff(Ffun39(r39),r39)
F2_39 = sp.diff(Ffun39(r39),r39,2)
gain1_39 = F1_39/(2*sp.sqrt(Ffun39(r39)))
gain2_39 = (F2_39/(2*sp.sqrt(Ffun39(r39)))
            -F1_39**2/(4*Ffun39(r39)**sp.Rational(3,2)))
radial_jets39 = True
for i39,vi39 in enumerate((x39,y39)):
    ri39 = 2*vi39
    Ri39 = 1 if i39==0 else sp.I
    first39 = Ri39*g39+(x39+sp.I*y39)*gain1_39.subs(
        r39,rs39)*ri39
    radial_jets39 &= sp.simplify(sp.diff(tformal39,vi39)-first39)==0
    for j39,vj39 in enumerate((x39,y39)):
        rj39 = 2*vj39
        Rj39 = 1 if j39==0 else sp.I
        second39 = (Ri39*gain1_39.subs(r39,rs39)*rj39
                    +Rj39*gain1_39.subs(r39,rs39)*ri39
                    +(x39+sp.I*y39)*(
                        gain2_39.subs(r39,rs39)*ri39*rj39
                        +gain1_39.subs(r39,rs39)*(2 if i39==j39 else 0)))
        radial_jets39 &= (
            sp.simplify(sp.diff(tformal39,vi39,vj39)-second39)==0)
rep('C39a source dominant-ratio branches have identical positive-root two-jets',
    selector_source39 and root_identity39 and radial_jets39)

# Take the strict open sector 99/100<r<101/100,
# |Arg R|<pi/5, |Arg u2|<pi/5.  For s=|t|^2, source profile inversion is
# r=z(s)=(B*s^2+A*s)/(1+A*s), with z'>0.  Rational bounds
# A<160 and 50<B<81 imply z(9/100)<99/100 and
# z(1/4)>101/100.  Hence 3/10<|t|<1/2 on the whole band.
s39 = sp.Symbol('s39',positive=True)
profile39 = (B39*s39**2+A39*s39)/(1+A39*s39)
profile_prime39 = (A39+B39*(2*s39+A39*s39**2))/(1+A39*s39)**2
low_gap39 = (sp.Rational(99,100)
             -sp.Rational(81*81+160*9,10000))
high_gap39 = (sp.Rational(50,16)-sp.Rational(101,100)
              -sp.Rational(160,400))
band_bounds39 = (
    sp.factor(sp.diff(profile39,s39)-profile_prime39)==0
    and sp.simplify(profile39.subs(s39,sfrom39)-r39)==0
    and sp.pi<sp.Rational(22,7)
    and (4*sp.Rational(22,7))**2<160
    and sp.E>sp.Rational(8,3)
    and sp.Rational(8,3)**4>50
    and sp.E<3 and sp.E**4<81
    and low_gap39==sp.Rational(1899,10000)>0
    and high_gap39==sp.Rational(343,200)>0
    and sp.simplify(profile39.subs(s39,sp.Rational(9,100))
                    -sp.Rational(99,100)).subs(
                        {A39:160,B39:81})<0
    and sp.simplify(profile39.subs(s39,sp.Rational(1,4))
                    -sp.Rational(101,100)).subs(
                        {A39:160,B39:50})>0
    and sp.Rational(9,100)/sp.E**2>sp.Rational(1,100)
    and sp.Rational(1,100)>sp.Rational(1,2500)
    and 1/(4*sp.E)<sp.Rational(1,8)
    and sp.pi/5<sp.pi/2
    and sp.Rational(1,10)<sp.Rational(1,4)
    and 1/sp.sqrt(sp.Rational(1,4))==2
    and sp.Rational(4,5)*sp.pi>sp.Rational(3,4)*sp.pi)
# q=s exp(-4s/(1-s)) lies in (1/100,1/8): the lower bound uses
# 4s/(1-s)<4/3 and e<3; the upper uses the increasing function
# s exp(-4s) on (0,1/4), whose endpoint is 1/(4e).
# Since |a|<1/10, the inverse-moment source takes its G branch.
moment_band39 = (
    band_bounds39
    and matching_assignment30(logradius_tree34,'lam',
        "scalar_flat_step((r-mp.mpf('.3'))/mp.mpf('.2'))")
    and matching_if_test32(moment_tree36,
        "lam==1 or aa<=mp.mpf('.1')")
    and sp.factor(sp.diff(s39*sp.exp(-4*s39),s39)
                  -sp.exp(-4*s39)*(1-4*s39))==0
    and sp.Rational(9,100)>0
    and sp.Rational(1,4)<1)
rep('C39b strict ratio band stays in the active south collar and height domain',
    moment_band39)

# Exact round witness and two strict selector-side witnesses.  Let
# lambda=|w|, w=(0,-lambda,0), y=(0,1-lambda,0),
# X=(lambda-1)/sqrt(1+lambda^2), n=(0,0,1), v=(1-X^2)n.
# At lambda0=2-sqrt(2), u1=u2=1/sqrt(2), so r=1 and t=e^-1.
# At lambda0 +/- 1/10000 the positive ratio crosses 1 strictly.
lam39 = sp.Symbol('lam39',positive=True)
X39 = (lam39-1)/sp.sqrt(1+lam39**2)
aa39 = (lam39**2-1)/(lam39**2+1)
v239 = 1-X39**2
y139 = 1-lam39
u139 = (1-y139**2)/(1+y139**2)
u239 = 2*y139/(1+y139**2)
Rround39 = sp.factor(u139/u239)
lam0_39 = 2-sp.sqrt(2)
eps39 = sp.Rational(1,10000)
y0_39 = sp.sqrt(2)-1
witness39 = (
    sp.simplify(X39**2*(2-X39**2)-aa39**2)==0
    and 0<lam0_39-eps39<lam0_39+eps39<1
    and X39.subs(lam39,lam0_39+eps39)<0
    and aa39.subs(lam39,lam0_39+eps39)<0
    and sp.simplify(v239/(1-aa39)-lam39)==0
    and sp.simplify(u139**2+u239**2-1)==0
    and sp.simplify(u139.subs(lam39,lam0_39)-1/sp.sqrt(2))==0
    and sp.simplify(u239.subs(lam39,lam0_39)-1/sp.sqrt(2))==0
    and sp.simplify(Rround39.subs(lam39,lam0_39)-1)==0
    and sp.Rational(2,5)<y0_39-eps39
    and y0_39+eps39<sp.Rational(1,2)
    and sp.simplify(sp.diff(Rround39,lam39)
                    -(1+y139**2)/(2*y139**2))==0
    and (1+sp.Rational(1,2)**2)/(2*sp.Rational(2,5)**2)<4
    and (1-sp.Rational(1,2500))**2>sp.Rational(99,100)
    and (1+sp.Rational(1,2500))**2<sp.Rational(101,100)
    and sp.Rational(99,100)<Rround39.subs(
        lam39,lam0_39-eps39)**2<1
    and 1<Rround39.subs(lam39,lam0_39+eps39)**2<sp.Rational(101,100)
    and sp.simplify(Fdirect39.subs({r39:1,A39:(4*sp.pi)**2,
                                   B39:sp.exp(4)})-sp.exp(-2))==0
    and sp.exp(-1)>sp.Rational(3,10)
    and sp.exp(-1)<sp.Rational(1,2))
rep('C39c exact round seam and strict open witnesses on both selector sides',
    witness39)

# On this angular sector Arg(t) is in (-pi/5,pi/5) and
# a=Arg(u2)/(2*pi) is in (-1/10,1/10).  Thus source middle data are
# exactly the same five entries (t,one,four,log(t),a) on both sides.
# The C34 real Cramer determinant is >981/5000, the south 2x2
# determinant is 1-s^2>15/16, and the height domain is strict.
# C34e reconstructs zz=e^rho in the south torus_pair.  The source
# flat weights then give b0=beta in (0,1), bi=0, bm=1-beta;
# hw(a)=1, the second torus pair has weight zero, and the two
# angular weights are 1,0.  Consequently the *actual* normalized
# ten-vector has the following three nonzero slots.
tr39,ti39,rho39,one39,beta39 = sp.symbols(
    'tr39 ti39 rho39 one39 beta39',real=True)
t39 = tr39+sp.I*ti39
yy39 = sp.exp(rho39)
xcore39 = t39/yy39
frame_den39 = sp.sqrt((tr39**2+ti39**2)/yy39**2+yy39**2)
vec_den39 = sp.sqrt(beta39**2+(1-beta39)**2)
E39 = sp.zeros(10,1)
E39[0] = beta39*sp.conjugate(xcore39)/(frame_den39*vec_den39)
E39[1] = beta39*yy39/(frame_den39*vec_den39)
E39[8] = (1-beta39)*sp.exp(2*sp.pi*sp.I*one39)/vec_den39
# SymPy does not infer positivity of beta^2+(1-beta)^2 from beta real.
# Its equality to 1/2+2(beta-1/2)^2 supplies that proof, so write the
# conjugate explicitly with the same positive denominators.
Ebar39 = sp.zeros(10,1)
Ebar39[0] = beta39*xcore39/(frame_den39*vec_den39)
Ebar39[1] = beta39*yy39/(frame_den39*vec_den39)
Ebar39[8] = (1-beta39)*sp.exp(-2*sp.pi*sp.I*one39)/vec_den39
P39 = E39*Ebar39.T
support39 = (0,1,8)
ten_slot39 = (
    south_reconstruction34 and south_det2_34 and cramer34
    and source_torus35 and selector_source39
    and Ebar39==E39.conjugate().xreplace(
        {sp.conjugate(1/vec_den39):1/vec_den39})
    and sp.simplify((Ebar39.T*E39)[0]-1)==0
    and all(E39[k]==0 for k in range(10) if k not in support39)
    and all(P39[i,j]==0 for i in range(10) for j in range(10)
            if i not in support39 or j not in support39)
    and all(sp.simplify((P39*P39-P39)[i,j])==0
            for i in support39 for j in support39)
    and sp.simplify(beta39**2+(1-beta39)**2-sp.Rational(1,2)
                    -2*(beta39-sp.Rational(1,2))**2)==0
    and det_margin34.subs(radius_bound34,sp.Rational(3,10))
        ==sp.Rational(981,5000)
    and sp.Rational(1)-sp.Rational(1,4)**2
        ==sp.Rational(15,16))
# The exact positive-root equality is on an open four-real input set.
# Source uses identical rho, Cramer data, and the displayed E39 there;
# hence all ten component functions, all 100 P entries, and their first
# and symmetric second derivatives are identical.  The same source
# phase/restore formulas take equal input, so the holomorphic target
# and its two-jets coincide.  This is source-level extensionality after
# a checked nontrivial root identity, not a claim that the enormous
# global projector loop was executed.  For any finite count whose
# denominators are nonzero, equality propagates index by index.
selector_jet39 = (
    ten_slot39 and root_identity39 and radial_jets39 and witness39
    and project_source30 and generic_projector31
    and source_cusp_ok30 and source_torus35)
rep('C39d ten-slot reference, rank-one projector and target two-jets agree openly',
    selector_jet39)

# C39e is a numerical observation (d), sample size one, at the exact
# round seam witness above.  It measures the margins next to the proved
# constants; it is not used by C39a--d or by their open-band argument.
import mpmath as mp39
with mp39.workdps(80):
    lam_sample39 = 2-mp39.sqrt(2)
    X_sample39 = (lam_sample39-1)/mp39.sqrt(1+lam_sample39**2)
    v_sample39 = (mp39.mpf(0),mp39.mpf(0),1-X_sample39**2)
    kind_sample39,data_sample39 = cusp30.elliptic_raw_base_reference(
        v_sample39,X_sample39,dps=70,elementary=True)
    t_sample39,_,_,lt_sample39,a_sample39 = data_sample39
    _,norm2_sample39 = cusp30.elliptic_bundle_embedding(
        kind_sample39,data_sample39,dps=70,with_norm=True,
        elementary=True)
    s_sample39 = abs(t_sample39)**2
    q_sample39 = s_sample39*mp39.exp(-4*s_sample39/(1-s_sample39))
    ez_sample39 = (-mp39.j*t_sample39,-t_sample39,
                   -mp39.j*lt_sample39)
    ew_sample39 = (1-t_sample39*t_sample39,
                   mp39.j*(1+t_sample39*t_sample39),2*t_sample39)
    rows_sample39 = [[z.imag for z in ez_sample39],
                     [z.real for z in ew_sample39],
                     [z.imag for z in ew_sample39]]
    cramer_sample39 = mp39.det(mp39.matrix(rows_sample39))
    e2_sample39,e3_sample39 = ew_sample39[:2]
    torus_sample39 = (e2_sample39.real*e3_sample39.imag
                      -e3_sample39.real*e2_sample39.imag)
    sample_ok39 = (
        kind_sample39=='mid'
        and abs(t_sample39-1/mp39.e)<mp39.mpf('1e-60')
        and abs(a_sample39)<mp39.mpf('1e-60')
        and abs(cramer_sample39-(1+s_sample39)**2)
            <mp39.mpf('1e-60')
        and abs(torus_sample39-(1-s_sample39**2))
            <mp39.mpf('1e-60')
        and q_sample39>mp39.mpf(1)/100
        and q_sample39<mp39.mpf(1)/8
        and cramer_sample39>mp39.mpf(981)/5000
        and torus_sample39>mp39.mpf(15)/16
        and norm2_sample39>mp39.mpf(1)/2)
    print('C39e one-point source margins (d), sample size 1:'
          ' q-1/100=',mp39.nstr(q_sample39-mp39.mpf(1)/100,18),
          ' 1/8-q=',mp39.nstr(mp39.mpf(1)/8-q_sample39,18),
          ' Cramer-981/5000=',mp39.nstr(
              cramer_sample39-mp39.mpf(981)/5000,18),
          ' torus-15/16=',mp39.nstr(
              torus_sample39-mp39.mpf(15)/16,18),
          ' norm2-1/2=',mp39.nstr(
              norm2_sample39-mp39.mpf(1)/2,18))
rep('C39e one-point production margins (numerical observation, d)',
    sample_ok39)

print('ALL PASS' if ok else 'FAILURES')
assert ok, 'A chart identity failed'
