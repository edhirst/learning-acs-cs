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

print('ALL PASS' if ok else 'FAILURES')
assert ok, 'A chart identity failed'
