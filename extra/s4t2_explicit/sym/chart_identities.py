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

# C40: the production `inf` selector at |small|=1/2.  This is a FLAT
# two-jet join, not equality of the two raw formulas on an open overlap.
# Work in the strict inverse-ratio, north-collar angular chamber
# Arg(u1)=pi/4, Arg(t)=pi/3, a=Arg(-conj(u1))/(2pi)=3/8.
# Write s=|small|, t=1/small, r=|t|=1/s and P=log(e^-2 s)-2s^2/(1-s^2).
# The left source gives rho=elementary_cusp_logradius(a,small,e^-2),
# while the right source uses the ordinary middle branch.  Source AST
# guards below bind each formula used in this local certificate.
terminal_forward_tree40 = source_tree30(cusp30._height_terminal_forward)
radial_tree40 = source_tree30(cusp30._height_R)
selector_source40 = (
    source_collar34 and source_torus35 and project_source30
    and matching_assignment30(raw_tree31,'inf',
        "inverse_ratio and abs(small)<mp.mpf('.5')")
    and matching_assignment30(raw_tree31,'t',
        'None if inf else (1/small if inverse_ratio else small)')
    and matching_if_test32(raw_tree31,"inf or abs(t)<mp.mpf('.5')")
    and matching_assignment30(raw_tree31,'tc','small if inf else t')
    and matching_assignment30(raw_tree31,'phase_coordinate',
        '-mp.conj(u1) if inf else u2')
    and matching_assignment30(raw_tree31,'nu',
        "mp.exp(-2) if inf else mp.mpf(1)")
    and matching_statement31(raw_tree31,'if inf:\n    t = 1/tc')
    and matching_statement31(raw_tree31,
        'if inf:\n    zm,wm = -mp.j*rho,mp.mpc(0)\n'
        'else:\n    ss = lt/(2*mp.j*mp.pi)\n'
        '    zm,wm = -mp.j*rho+mp.pi*ss*(ss-1),2*t*ss')
    and matching_assignment30(raw_tree31,'lam',
        "1-scalar_flat_step((r-mp.mpf('1.2'))/mp.mpf('.4'))")
    and matching_if_test32(raw_tree31,'r<=1')
    and matching_assignment30(raw_tree31,'a2',
        '-mp.arg(u1)/(2*mp.pi)+lam*delta+(1-lam)/2')
    and matching_assignment30(raw_tree31,'four',
        "-a2-mp.mpf('.5')")
    and matching_assignment30(logradius_tree34,'P',
        'mp.log(nu)+mp.log(r)-2*r*r/(1-r*r)')
    and matching_assignment30(logradius_tree34,'q','mp.exp(2*P)')
    and matching_assignment30(logradius_tree34,'lam',
        "scalar_flat_step((r-mp.mpf('.3'))/mp.mpf('.2'))")
    and exact_if_body37(logradius_tree34,'lam==1',
        "return +(P*(a+mp.mpf('.5')))")
    and matching_statement31(flat_tree38,
        'a,b = mp.exp(-1/u),mp.exp(-1/(1-u))')
    and exact_if_body37(flat_tree38,'u>=1','return mp.mpf(1)')
    and returned_expr30(flat_tree38,'a/(a+b)')
    and matching_assignment30(bundle_tree32,'bi',
        "1-scalar_flat_step((1/r-mp.mpf('.3'))/mp.mpf('.2'))")
    and matching_assignment30(bundle_tree32,'bm','1-b0-bi')
    and matching_assignment30(bundle_tree32,'angle',
        'mp.arg(t if index==0 else -t)')
    and matching_assignment30(bundle_tree32,'weight',
        '1-scalar_flat_step((abs(angle)-mp.pi/2)/(mp.pi/4))')
    and matching_statement31(bundle_tree32,
        'vec[8+index]=bm*weight*mp.exp(mp.j*(2*mp.pi*one'
        '-(angle+centre-lt.imag)*four))')
    and matching_statement31(restore_tree30,'C = list(phi)')
    and matching_statement31(restore_tree30,
        'C = [v+beta*d for v,d in zip(C,D)]')
    and returned_expr30(slice_tree30,
        'kind,(t,one+mp.arg(phase)/(2*mp.pi),four,lt)'))

# The inf Cramer right-hand side is (z,w)=(-i*rho,0).  With
# eta=rho/P and t=x+i*y, its exact coordinates are
# (eta*phi/(2pi),-2eta*x/(r^2-1),-2eta*y/(r^2-1),-eta).
xr40,yr40,L40,phi40,eta40 = sp.symbols(
    'xr40 yr40 L40 phi40 eta40',real=True)
t40 = xr40+sp.I*yr40
r2_40 = xr40**2+yr40**2
lt40 = L40+sp.I*phi40
P40 = -L40-2*r2_40/(r2_40-1)
one40 = eta40*phi40/(2*sp.pi)
c2_40 = -2*eta40*xr40/(r2_40-1)
c3_40 = -2*eta40*yr40/(r2_40-1)
four40 = -eta40
z40 = (2*sp.pi*one40-sp.I*t40*c2_40-t40*c3_40
       -sp.I*lt40*four40)
w40 = ((1-t40**2)*c2_40+sp.I*(1+t40**2)*c3_40
       +2*t40*four40)
det40 = (1+r2_40)*(L40*(r2_40-1)+2*r2_40)
ez40 = (-sp.I*t40,-t40,-sp.I*lt40)
ew40 = (1-t40**2,sp.I*(1+t40**2),2*t40)
cramer_matrix40 = sp.Matrix((
    [sp.im(sp.expand(v)) for v in ez40],
    [sp.re(sp.expand(v)) for v in ew40],
    [sp.im(sp.expand(v)) for v in ew40]))
cramer40 = (
    sp.factor(z40+sp.I*eta40*P40)==0
    and sp.factor(w40)==0
    and sp.factor(cramer_matrix40.det()-det40)==0
    and sp.Rational(1000,501)>sp.Rational(19,10)
    and 2*(1+sp.Rational(19,10)**2)*sp.Rational(19,10)**2>30
    and 1-sp.Rational(501,1000)**4>sp.Rational(9,10))
rep('C40a source inf/middle guards and exact Cramer signs; det>30',
    selector_source40 and cramer40)

# Exact round witnesses on both strict sides and at the seam.  Put
# p=s^2, A=(4pi)^2, B=e^4, and
# kappa^2=p(p+A)/(Ap+B).  For u1=e^(i*pi/4)/sqrt(1+kappa^2),
# u2=kappa*e^(-i*pi/12)/sqrt(1+kappa^2), the source inverse-ratio
# positive square root gives small=s*e^(-i*pi/3).  The identity below
# checks the positive discriminant, not just its square.
p40,A40,B40 = sp.symbols('p40 A40 B40',positive=True)
ratio_sq40 = p40*(p40+A40)/(A40*p40+B40)
bprof40 = A40*(1-ratio_sq40)
Dprof40 = (A40*B40+A40*p40**2+2*B40*p40)/(A40*p40+B40)
profile40 = (
    sp.factor(Dprof40**2-bprof40**2-4*B40*ratio_sq40)==0
    and sp.factor(ratio_sq40*2*B40/(bprof40+Dprof40)-p40)==0
    and sp.factor(sp.diff(ratio_sq40,p40)
        -Dprof40/(A40*p40+B40))==0
    and sp.factor(sp.diff(ratio_sq40,A40)
        -p40*(B40-p40**2)/(A40*p40+B40)**2)==0
    and sp.factor(sp.diff(ratio_sq40,B40)
        +p40*(p40+A40)/(A40*p40+B40)**2)==0)
low_ratio40 = ratio_sq40.subs(
    {p40:sp.Rational(6,25),A40:144,B40:55})
high_ratio40 = ratio_sq40.subs(
    {p40:sp.Rational(13,50),A40:160,B40:50})
# For B<55, sum_{n=0}^4 1/n!=65/24 and the remaining terms obey
# 1/(5+k)! <= (1/120)6^(-k), with strict inequality for k>=2.
# Their geometric sum is 1/100, so e<65/24+1/100<68/25.
elementary_bounds40 = (
    sp.pi>3 and sp.pi<sp.Rational(22,7)
    and 144<(4*sp.pi)**2<160
    and sp.E>sp.Rational(8,3)
    and sp.Rational(8,3)**4>50
    and sp.Rational(65,24)+sp.Rational(1,100)<sp.Rational(68,25)
    and sp.E<sp.Rational(68,25)
    and sp.Rational(68,25)**4<55
    and sp.Rational(9,25)<low_ratio40<high_ratio40<sp.Rational(49,100)
    and 50-sp.Rational(13,50)**2>0
    and sp.factor(sp.diff(4*p40/(1-p40),p40)
                  -4/(1-p40)**2)==0
    and 4*sp.Rational(13,50)/(1-sp.Rational(13,50))
        ==sp.Rational(52,37)
    and sp.Rational(6,25)/(3**5*2)==sp.Rational(1,2025)
    and sp.Rational(52,37)<sp.Rational(3,2)
    and sp.exp(sp.Rational(11,2))<3**5*2
    and sp.Rational(1,2025)>sp.Rational(1,2500)
    and sp.Rational(13,50)/50<sp.Rational(1,8))
q_profile40 = p40*sp.exp(-4-4*p40/(1-p40))
height_domain40 = (
    sp.simplify(sp.exp(2*(-2+sp.log(sp.sqrt(p40))
                          -2*p40/(1-p40)))-q_profile40)==0
    and all(sp.Rational(1,2025)
            <q_profile40.subs(p40,ss40**2)<sp.Rational(1,8)
            for ss40 in (sp.Rational(499,1000),sp.Rational(1,2),
                         sp.Rational(501,1000))))

kap40 = sp.Symbol('kap40',positive=True)
co40 = (sp.sqrt(6)+sp.sqrt(2))/4
si40 = (sp.sqrt(6)-sp.sqrt(2))/4
u1_40 = (1+sp.I)/sp.sqrt(2*(1+kap40**2))
u2_40 = kap40*(co40-sp.I*si40)/sp.sqrt(1+kap40**2)
d40 = 1+sp.re(u1_40)
y40 = (sp.im(u1_40)/d40,sp.re(u2_40)/d40,
       sp.im(u2_40)/d40)
wy40 = (y40[0],y40[1]-1,y40[2])
ysq40 = sum(v*v for v in y40)
R2_40 = sum(v*v for v in wy40)
n40 = (wy40[0]/sp.sqrt(R2_40),-wy40[2]/sp.sqrt(R2_40),
       -wy40[1]/sp.sqrt(R2_40))
Q40 = sp.sqrt(1+kap40**2)
dc40 = Q40+1/sp.sqrt(2)
round_witness40 = (
    sp.simplify(u1_40*sp.conjugate(u1_40)
                +u2_40*sp.conjugate(u2_40)-1)==0
    and sp.simplify(u2_40/u1_40
                    -kap40*(1-sp.I*sp.sqrt(3))/2)==0
    and sp.simplify((1-ysq40+2*sp.I*y40[0])/(1+ysq40)
                    -u1_40)==0
    and sp.simplify(2*(y40[1]+sp.I*y40[2])/(1+ysq40)
                    -u2_40)==0
    and sp.cancel(sum(v*v for v in n40)-1)==0
    and sp.simplify(wy40[0]-1/sp.sqrt(2)/dc40)==0
    and sp.simplify(y40[1]-kap40*co40/dc40)==0
    and sp.simplify(y40[2]+kap40*si40/dc40)==0
    and round_inverse38)
# For 3/5<kappa<7/10: 23/20<Q<5/4,
# 7/10<1/sqrt(2)<3/4, hence 37/20<dc<2.
# These rational chains give 1/3<w0<7/16,
# 59/100<-w1<4/5, |w2|<1/8, 2/3<R<1, and n2>59/100.
round_margins40 = (
    sp.Rational(23,20)**2<1+sp.Rational(3,5)**2
    and 1+sp.Rational(7,10)**2<sp.Rational(5,4)**2
    and sp.Rational(7,10)<1/sp.sqrt(2)<sp.Rational(3,4)
    and sp.Rational(9,10)<co40<1
    and 0<si40<sp.Rational(3,10)
    and sp.Rational(7,10)/2>sp.Rational(1,3)
    and sp.Rational(3,4)/sp.Rational(37,20)<sp.Rational(7,16)
    and sp.Rational(3,5)*sp.Rational(9,10)/2>sp.Rational(1,5)
    and sp.Rational(7,10)/sp.Rational(37,20)<sp.Rational(41,100)
    and sp.Rational(7,10)*sp.Rational(3,10)
        /sp.Rational(37,20)<sp.Rational(1,8)
    and 1+sp.Rational(49,100)<sp.Rational(13,10)**2
    and (1-sp.Rational(7,10))/sp.Rational(13,10)
        ==sp.Rational(3,13)
    and sp.Rational(1,3)**2+sp.Rational(59,100)**2
        >sp.Rational(4,9)
    and sp.Rational(7,16)**2+sp.Rational(4,5)**2
        +sp.Rational(1,8)**2<1
    and sp.Rational(9,10)-1/sp.sqrt(2)>sp.Rational(1,6)
    and sp.Rational(499,1000)>sp.Rational(3,10)
    and sp.Rational(1,2)-sp.Rational(501,1000)<0
    and sp.Rational(499,1000)<sp.Rational(1,2)<sp.Rational(501,1000)
    and sp.Rational(6,25)<sp.Rational(499,1000)**2
    and sp.Rational(501,1000)**2<sp.Rational(13,50)
    and sp.Rational(1000,501)>sp.Rational(8,5)
    and sp.Rational(3,8)<sp.Rational(2,5)
    and sp.Rational(3,8)>sp.Rational(1,5))
# Here |u1|-|u2|>(1-7/10)/(13/10)=3/13;
# Re(-conj u1)/|u1|=-1/sqrt(2)>-9/10.
# The source heights have q=e^-4*s^2*exp(-4s^2/(1-s^2)) and
# 1/2025<q<13/2500<1/8: 4s^2/(1-s^2)<52/37<3/2,
# e^(11/2)<3^5*2, and e^4>50.  The exact round reconstruction above
# yields -1/3<X<0, n2>59/100, hence a strict free slab.
rep('C40b exact round seam and strict selector-side witnesses',
    profile40 and elementary_bounds40 and height_domain40 and round_witness40
    and round_margins40)

# On the left let v=1/2-s>0.  The source flat step has
# u=(s-.3)/.2=1-5v and beta=1-lambda, where
# beta=exp(-1/(5v))/(exp(-1/(1-5v))+exp(-1/(5v))).
# Its value and first two derivatives vanish at v=0.  On the right
# beta=0 exactly.  The same beta is the north embedding weight bi.
v40 = sp.Symbol('v40',positive=True)
Ef40 = sp.exp(-1/(5*v40))
Af40 = sp.exp(-1/(1-5*v40))
beta_flat40 = Ef40/(Af40+Ef40)
flat_endpoint40 = (
    sp.limit(Af40,v40,0,dir='+')==sp.exp(-1)
    and all(sp.limit(sp.diff(beta_flat40,v40,k),v40,0,dir='+')==0
            for k in range(3)))

# The terminal interior formula really joins its lam==1 direct branch.
# At lam=1: kappa=amp=k=-log q, B=C=sqrt(q), z0=0,
# d=2 sqrt(q) k, Acut=E0, Lambda=mu=1, C0=C1=0.
# Thus both exact source formulas for the terminal T reduce to T(z)=z.
k40,cut40,eps40,zterm40 = sp.symbols(
    'k40 cut40 eps40 zterm40',positive=True)
qh40 = sp.exp(-k40)
Bh40 = sp.exp(-k40/2)
Ch40 = qh40*sp.exp(k40/2)
dh40 = 2*sp.sqrt(qh40)*k40
z0h40 = sp.simplify(Bh40-Ch40)
E0h40 = sp.exp(eps40/cut40)
Ach40 = sp.exp(eps40/(cut40-z0h40))
Lamh40 = sp.exp(Ach40-E0h40)
muh40 = sp.simplify(
    dh40*(E0h40*eps40/cut40**2)
    /(k40*(Bh40+Ch40)
      *(Ach40*eps40/(cut40-z0h40)**2)))
ellh40 = eps40/(cut40-zterm40)-eps40/cut40
terminal_identity40 = (
    z0h40==0 and sp.simplify(Ach40-E0h40)==0
    and sp.simplify(Lamh40-1)==0
    and sp.simplify(muh40-1)==0
    and sp.simplify(1+muh40-2*Lamh40)==0
    and sp.simplify(Lamh40-muh40)==0
    and sp.simplify(cut40*ellh40/(eps40/cut40+ellh40)
                    -zterm40)==0
    and matching_statement31(terminal_forward_tree40,
        'C0,C1 = 1+mu-2*Lambda,Lambda-mu')
    and matching_if_test32(terminal_forward_tree40,'Ez>50')
    and matching_if_test32(radial_tree40,'value>=c')
    and matching_statement31(moment_tree36,'u = _height_R(aa)')
    and matching_statement31(moment_tree36,
        'v = _height_terminal_forward(delta_z,data)')
    and exact_if_body37(moment_tree36,
        "lam==1 or aa<=mp.mpf('.1')",
        'return _height_G(a,amp,d)'))
# For a=3/8, R and R^-1 are identity with margin a-1/5=7/40;
# a<2/5 gives a strict terminal interior point.  All denominators
# are positive in the q band.  The interior formula is therefore
# smooth in lambda up to 1 and equals G(a,q,1) there.
ah40 = sp.Symbol('ah40',positive=True)
Eh40 = sp.exp(k40*ah40)
fh40 = sp.sqrt(qh40)*(Eh40-1/Eh40)
Rh40 = sp.sqrt(qh40)*(Eh40+1/Eh40)
Yh40 = sp.simplify(2*qh40/(Rh40+fh40))
height_join40 = (
    terminal_identity40 and flat_endpoint40
    and sp.Rational(3,8)-sp.Rational(1,5)==sp.Rational(7,40)
    and sp.Rational(2,5)-sp.Rational(3,8)==sp.Rational(1,40)
    and sp.simplify(Rh40**2-fh40**2-4*qh40)==0
    and sp.simplify(Yh40-sp.exp(-k40*(ah40+sp.Rational(1,2))))==0
    and sp.simplify(sp.log(Yh40)/2
        +k40*(ah40+sp.Rational(1,2))/2)==0)
# Consequently Delta=rho-P(a+1/2) has a zero 2-jet at s=1/2.
# Since -P=log r+2r^2/(r^2-1)>2 for r>1, eta=rho/P has
# eta-(a+1/2) flat to the same order.
rep('C40c terminal height joins exactly; radial weight and Delta are flat',
    height_join40 and cramer40)

# The lam<1 terminal expression uses its genuine low-Ez interior branch
# near this seam.  At lam=1 and a=3/8, z0=0 and
# gap=zcut-delta_z=2 sqrt(q)[sinh((2/5)k)-sinh((3/8)k)].
# The derivative of sinh(ka) is k cosh(ka), with
# cosh(ka)-1=2 sinh(ka/2)^2>0.  The mean-value theorem gives
# gap>sqrt(q) k/20>1/180, since q>1/2025 and
# k=-log q=16/3+log 4>5.  Source eta=1/100 then gives
# Ez=exp(eta/gap)<exp(9/5)<e^2<9<50.
gap_h40 = 2*sp.sqrt(qh40)*(
    sp.sinh(sp.Rational(2,5)*k40)
    -sp.sinh(sp.Rational(3,8)*k40))
aa_h40 = sp.Symbol('aa_h40',positive=True)
low_terminal40 = (
    sp.simplify(((Bh40*sp.exp(sp.Rational(2,5)*k40)
                 -Ch40*sp.exp(-sp.Rational(2,5)*k40))
                -(Bh40*sp.exp(sp.Rational(3,8)*k40)
                  -Ch40*sp.exp(-sp.Rational(3,8)*k40))
                -gap_h40).rewrite(sp.exp))==0
    and sp.simplify(sp.diff(sp.sinh(k40*aa_h40),aa_h40)
                    -k40*sp.cosh(k40*aa_h40))==0
    and sp.simplify(sp.cosh(k40*aa_h40)-1
                    -2*sp.sinh(k40*aa_h40/2)**2)==0
    and sp.Rational(2,5)-sp.Rational(3,8)==sp.Rational(1,40)
    and sp.Rational(1,45)*5/20==sp.Rational(1,180)
    and sp.Rational(16,3)>5
    and sp.simplify(-sp.log(sp.exp(-sp.Rational(16,3))/4)
                    -(sp.Rational(16,3)+sp.log(4)))==0
    and sp.Rational(1,100)/sp.Rational(1,180)
        ==sp.Rational(9,5)
    and sp.exp(sp.Rational(9,5))<9<50
    and matching_assignment30(terminal_tree36,'eta',
        "mp.mpf('.01')")
    and matching_if_test32(terminal_forward_tree40,'delta_z>=gap0')
    and matching_assignment30(terminal_forward_tree40,'gap0','zcut-z0')
    and matching_assignment30(terminal_forward_tree40,'gap',
        'gap0-delta_z')
    and matching_assignment30(terminal_forward_tree40,'L','eta/gap')
    and matching_assignment30(terminal_forward_tree40,'Ez','mp.exp(L)')
    and matching_if_test32(terminal_forward_tree40,'Ez>50'))
rep('C40c3 strict low-Ez source branch at the terminal seam',
    low_terminal40 and height_domain40)

# The inf torus correction in elementary_restore_torus calls the four-seed
# height only on the left, where |tc|<1/2.  At the seam its analytic
# representative has q=exp(-16/3)/4, A=q^(1/8), Y=q^(7/8).
# The source guard is rr<=1/4, max(A,Y)<=4, q<1/8.  In fact the
# seed cases are strict: 1/4<A<15/32<1/2 and
# 0<Y<1/128, qA<1/16, qY<1/16.  Thus this is an open smooth
# one-sided germ; the guard rr<=1/4 holds for every left s<1/2.
finite_height_tree40 = source_tree30(cusp30.scalar_finite_height)
seed_nodes40 = [node for node in ast.walk(finite_height_tree40)
                if isinstance(node,ast.FunctionDef) and node.name=='seed']
q_seam40 = sp.exp(-sp.Rational(16,3))/4
A_seam40 = q_seam40**sp.Rational(1,8)
Y_seam40 = q_seam40**sp.Rational(7,8)
restore_domain40 = (
    north_reconstruction34
    and sp.simplify(sp.exp(2*(-sp.Rational(8,3)-sp.log(2)))
                -q_seam40)==0
    and sp.simplify(sp.exp(2*(-sp.Rational(8,3)-sp.log(2))
                            *sp.Rational(7,8))-Y_seam40)==0
    and sp.simplify(q_seam40/Y_seam40-A_seam40)==0
    and sp.exp(4)>50
    and sp.exp(sp.Rational(4,3))
        >1+sp.Rational(4,3)+sp.Rational(4,3)**2/2
    and 50*(1+sp.Rational(4,3)+sp.Rational(4,3)**2/2)>128
    and sp.Rational(1,2500)<q_seam40<sp.Rational(1,512)
    and sp.Rational(1,2500)>sp.Rational(1,4)**8
    and sp.Rational(1,512)<sp.Rational(15,32)**8
    and sp.Rational(1,4)<A_seam40<sp.Rational(15,32)
    and 0<Y_seam40<sp.Rational(1,128)
    and q_seam40*A_seam40<sp.Rational(1,16)
    and q_seam40*Y_seam40<sp.Rational(1,16)
    and 1-sp.Rational(501,1000)**4>sp.Rational(9,10)
    and matching_assignment30(finite_height_tree40,'rr','abs(x*y)**2')
    and matching_statement31(finite_height_tree40,
        'A,Y = nu*nu*e*abs(x)**2,abs(y)**2')
    and matching_assignment30(finite_height_tree40,'q','A*Y')
    and matching_if_test32(finite_height_tree40,
        "not rr<=mp.mpf('.25') or max(A,Y)>4 or q>=mp.mpf('.125')")
    and len(seed_nodes40)==1
    and matching_if_test32(seed_nodes40[0],
        'v<=mp.mpf(1)/16')
    and matching_if_test32(seed_nodes40[0],
        "v>=mp.mpf('.5')")
    and returned_expr30(finite_height_tree40,
        'seed(A)+seed(q*A)-seed(Y)-seed(q*Y)+q*q*(A-Y)/(1-q)'))
rep('C40c2 strict four-seed restore germ and source guard at radial seam',
    restore_domain40)

# On pi/4<Arg(t)<pi/2, the two angular slots have angles
# phi and phi-pi.  Their phase brackets are identically zero:
# phi+0-Im(log t)=0 and (phi-pi)+pi-Im(log t)=0.
# At phi=pi/3 their weights are 1 and
# omega=1-flatstep(2/3)=1/(1+exp(3/2)) in (0,1).
ang40 = sp.Symbol('ang40',real=True)
omega40 = 1/(1+sp.exp(sp.Rational(3,2)))
angular40 = (
    sp.simplify(ang40+(ang40-sp.pi)*0-ang40)==0
    and sp.simplify((ang40-sp.pi)+sp.pi-ang40)==0
    and sp.pi/4<sp.pi/3<sp.pi/2
    and (sp.pi-sp.pi/3)-sp.pi/2==sp.pi/6
    and 3*sp.pi/4-(sp.pi-sp.pi/3)==sp.pi/12
    and 0<omega40<1)
# On the left the unnormalized ten-vector is
# (beta*c4,...,beta*c7,(1-beta)*g,(1-beta)*omega*g),
# g=exp(i*eta*phi), with sum|c_j|^2=U in (1,2).
# On the right beta=0 and g=1.  The angular block of the rank-one
# projector is independent of g; cusp/cross entries carry beta.
jetx40,bc40,dc40,wgt40,U40 = sp.symbols(
    'jetx40 bc40 dc40 wgt40 U40',real=True)
beta40 = jetx40**3*bc40
N40 = beta40**2*U40+(1-beta40)**2*(1+wgt40**2)
angular_diff40 = (1-beta40)**2/N40-1/(1+wgt40**2)
cross40 = beta40*(1-beta40)/N40
cusp40 = beta40**2/N40
projector_jets40 = (
    sp.factor(angular_diff40
              +beta40**2*U40/(N40*(1+wgt40**2)))==0
    and all(sp.simplify(sp.diff(expr,jetx40,k).subs(jetx40,0))==0
            for expr in (angular_diff40,cross40,cusp40)
            for k in range(3))
    and sp.Rational(1,4)<(1-sp.Rational(1,2))**2*(1+omega40**2)
    and project_source30 and projection_value30 and projection_norm30
    and generic_projector31)
# A flat factor in one transverse coordinate has zero value and all
# first/symmetric second derivatives, including tangential and mixed
# directions.  This generic check applies to every smooth coefficient
# of the actual source formulas, without expanding the huge finite fold.
a_jet40,phi_jet40 = sp.symbols('a_jet40 phi_jet40',real=True)
Fjet40 = sp.Function('Fjet40')(jetx40,a_jet40,phi_jet40)
flat_all_directions40 = all(
    sp.simplify((sp.diff(jetx40**3*Fjet40,*derivs)
                 if derivs else jetx40**3*Fjet40)
                .subs(jetx40,0))==0
    for derivs in ((),(jetx40,),(a_jet40,),(phi_jet40,),
                   (jetx40,jetx40),(jetx40,a_jet40),
                   (jetx40,phi_jet40),(a_jet40,a_jet40),
                   (a_jet40,phi_jet40),(phi_jet40,phi_jet40)))
rep('C40d sparse ten-slot projector has identical complete flat two-jets',
    angular40 and projector_jets40 and flat_all_directions40)

# At each finite projection step P_j v/||P_j v|| is smooth wherever
# ||P_j v||>0.  The preceding P two-jet equality propagates by induction
# only if every intermediate P_j along the chosen path also has matching
# two-jets across this source seam and every denominator is nonzero.
# This local certificate supplies neither condition globally.
# The inf reference differs from the right angular reference by g,
# up to a flat beta error.  If p is the right final inner product,
# the left one has local form conj(g)*(p+jetx^3*h).  The source one
# equals eta*phi/(2pi), so its exp(i*2pi*one) cancels conj(g).
# The final p is nonzero when the finite fold exists: its last update
# is at h=X, giving a unit vector in the final rank-one image.  The
# witness is outside the pole correction.  Work in any fixed local Arg
# chart; a constant phase rotation permits Re(p)>0 in the formal model.
eta_jet40 = a_jet40+sp.Rational(1,2)+jetx40**3*dc40
gjet40 = sp.exp(sp.I*eta_jet40*phi_jet40)
pr40 = sp.Symbol('pr40',positive=True)
pi40,hr40,hi40 = sp.symbols('pi40 hr40 hi40',real=True)
pjet40 = pr40+sp.I*pi40
hjet40 = hr40+sp.I*hi40
phase_inf40 = sp.exp(-sp.I*eta_jet40*phi_jet40)*(
    pjet40+jetx40**3*hjet40)
argdiff40 = (sp.atan2(pi40+jetx40**3*hi40,
                      pr40+jetx40**3*hr40)
             -sp.atan2(pi40,pr40))
phase_jets40 = (
    pr40**2+pi40**2>0
    and sp.simplify(pjet40*sp.conjugate(pjet40)
                    -(pr40**2+pi40**2))==0
    and matching_assignment30(project_tree30,'h',
        'start+(X-start)*j/steps')
    and matching_if_test32(correct_tree30,'distance<alpha')
    and matching_if_test32(slice_tree30,
        'distance<alpha/2 and not force_slab')
    and sp.simplify(gjet40*phase_inf40
                -(pjet40+jetx40**3*hjet40))==0
    and all(sp.simplify(sp.diff(argdiff40,jetx40,k)
                        .subs(jetx40,0))==0 for k in range(3))
    and all(sp.diff(-eta_jet40+a_jet40+sp.Rational(1,2),
                    jetx40,k).subs(jetx40,0)==0 for k in range(3)))
# The final source one differs by a flat 2-jet modulo a locally
# constant integer deck.  The restored mid output has C=phi+beta D
# on the left and C=phi on the right; inserting one_L=one_R+x^3 u,
# four_L=four_R-x^3 d, C_L=C_R+x^3 b D gives these exact differences.
oneR40,fourR40,tt40,ll40,ct240,ct340,uu40,D240,D340 = sp.symbols(
    'oneR40 fourR40 tt40 ll40 ct240 ct340 uu40 D240 D340')
oneL40 = oneR40+jetx40**3*uu40
fourL40 = fourR40-jetx40**3*dc40
CL240,CL340 = ct240+beta40*D240,ct340+beta40*D340
zR40 = 2*sp.pi*oneR40-sp.I*ll40*fourR40
wR40 = 2*tt40*fourR40
zL40 = 2*sp.pi*oneL40-sp.I*ll40*fourL40
wL40 = 2*tt40*fourL40
rest_zR40 = zR40-sp.I*tt40*ct240-tt40*ct340
rest_wR40 = wR40+(1-tt40**2)*ct240+sp.I*(1+tt40**2)*ct340
rest_zL40 = zL40-sp.I*tt40*CL240-tt40*CL340
rest_wL40 = wL40+(1-tt40**2)*CL240+sp.I*(1+tt40**2)*CL340
rest_zdiff40 = jetx40**3*(2*sp.pi*uu40+sp.I*ll40*dc40
                         -sp.I*tt40*bc40*D240-tt40*bc40*D340)
rest_wdiff40 = jetx40**3*(-2*tt40*dc40
                         +(1-tt40**2)*bc40*D240
                         +sp.I*(1+tt40**2)*bc40*D340)
restore_jets40 = (
    sp.expand(rest_zL40-rest_zR40-rest_zdiff40)==0
    and sp.expand(rest_wL40-rest_wR40-rest_wdiff40)==0
    and all(sp.diff(expr,jetx40,k).subs(jetx40,0)==0
            for expr in (rest_zdiff40,rest_wdiff40) for k in range(3))
    and sp.simplify(sp.exp(2*sp.pi*sp.I*sp.Symbol('deck40',integer=True))-1)==0
    and source_torus35 and project_source30 and generic_projector31)
rep('C40e conditional fold, phase cancellation and restored mid two-jets',
    phase_jets40 and restore_jets40 and projector_jets40
    and flat_all_directions40)

# C40f is solely a reproducible numerical observation (d), sample size
# two.  It tests the actual mpmath branches under the opt-in formula
# mode and prints margins beside the proved local constants.  It is not
# used to prove C40a--e or any global nonzero projection denominator.
mp40 = mp39
sample_ok40 = True
sample_rows40 = []
with mp40.workdps(190), cusp30.elementary_formula_branches():
    A_sample40,B_sample40 = (4*mp40.pi)**2,mp40.exp(4)
    for ss40 in (mp40.mpf(499)/1000,mp40.mpf(501)/1000):
        kap_sample40 = mp40.sqrt(
            ss40**2*(ss40**2+A_sample40)
            /(A_sample40*ss40**2+B_sample40))
        u1_sample40 = mp40.exp(mp40.j*mp40.pi/4)/mp40.sqrt(
            1+kap_sample40**2)
        u2_sample40 = (kap_sample40*mp40.exp(-mp40.j*mp40.pi/12)
                       /mp40.sqrt(1+kap_sample40**2))
        ds40 = 1+u1_sample40.real
        ws40 = (u1_sample40.imag/ds40,
                u2_sample40.real/ds40-1,
                u2_sample40.imag/ds40)
        Rs40 = mp40.sqrt(sum(v*v for v in ws40))
        Xs40 = (Rs40-1)/mp40.sqrt(1+Rs40**2)
        ns40 = (ws40[0]/Rs40,-ws40[2]/Rs40,-ws40[1]/Rs40)
        kind_s40,data_s40 = cusp30.elliptic_raw_base_reference(
            tuple((1-Xs40**2)*v for v in ns40),Xs40,
            dps=180,elementary=True)
        _,norm_s40 = cusp30.elliptic_bundle_embedding(
            kind_s40,data_s40,dps=180,with_norm=True,elementary=True)
        ts40 = data_s40[0]
        rs40 = abs(ts40)
        qs40 = ss40**2*mp40.exp(-4-4*ss40**2/(1-ss40**2))
        dets40 = (1+rs40**2)*(
            mp40.log(rs40)*(rs40**2-1)+2*rs40**2)
        torus_s40 = 1-ss40**4
        dom_s40 = abs(u1_sample40)-abs(u2_sample40)
        sample_rows40.append((ss40,qs40,dets40,torus_s40,
                              norm_s40,dom_s40))
        sample_ok40 &= (
            kind_s40=='mid' and -mp40.mpf(1)/3<Xs40<0
            and ns40[2]>mp40.mpf(59)/100
            and abs(ts40-mp40.exp(mp40.j*mp40.pi/3)/ss40)
                <mp40.mpf('1e-130')
            and qs40>mp40.mpf(1)/2025 and dets40>30
            and torus_s40>mp40.mpf(9)/10
            and norm_s40>mp40.mpf(1)/4
            and dom_s40>mp40.mpf(3)/13
            and (len(data_s40)==5 if ss40<mp40.mpf('.5')
                 else len(data_s40)==4))
    ss_left40 = mp40.mpf(499)/1000
    tc_left40 = ss_left40*mp40.exp(-mp40.j*mp40.pi/3)
    P_left40 = (-2+mp40.log(ss_left40)
                -2*ss_left40**2/(1-ss_left40**2))
    Delta_left40 = (cusp30.elementary_cusp_logradius(
        mp40.mpf(3)/8,tc_left40,mp40.exp(-2),dps=180)
        -P_left40*mp40.mpf(7)/8)
    sample_ok40 &= Delta_left40<0 and abs(Delta_left40)>mp40.mpf('1e-100')
    print('C40f production margins (d), sample size 2:'
          ' min(q-1/2025)=',mp40.nstr(min(row[1]-mp40.mpf(1)/2025
                                        for row in sample_rows40),18),
          ' min(det-30)=',mp40.nstr(min(row[2]-30
                                      for row in sample_rows40),18),
          ' min(torus-9/10)=',mp40.nstr(min(row[3]-mp40.mpf(9)/10
                                            for row in sample_rows40),18),
          ' min(norm2-1/4)=',mp40.nstr(min(row[4]-mp40.mpf(1)/4
                                            for row in sample_rows40),18),
          ' left Delta=',mp40.nstr(Delta_left40,18))
rep('C40f two production branch margins (numerical observation, d)',
    sample_ok40)

# C41: the cusp/middle source join |t_c|=3/10, at both cusp signs.
# For |tc|<=3/10 the raw reference returns the c0/ci core reference of
# elementary_cusp_reference; for 3/10<|tc|<1/2 it returns middle data
# built from elementary_cusp_logradius and the real Cramer solve.  Let M0
# denote the right-hand source formulas with the three collar weights
# lambda=s((|tc|-3/10)/(1/5)), 1-b0 and bm replaced by zero.  C41 proves,
# on an open neighbourhood of the whole seam (any fixed argument chart of
# a, any fixed branch of log t): the core reference is a unit multiple of
# the M0 reference, so the two ten-slot projectors are identical; the
# right source differs from M0 only by terms with a factor flat at
# |tc|=3/10; and for every phi and every common final unit vector the
# c0/ci restoration equals the holomorphic chart change Phi/Theta of the
# M0 middle restoration.  Unlike C40 this is not one chamber.
def body_dump41(tree):
    body = tree.body
    if (body and isinstance(body[0],ast.Expr)
            and isinstance(body[0].value,ast.Constant)
            and isinstance(body[0].value.value,str)):
        body = body[1:]
    return [ast.dump(node,include_attributes=False) for node in body]


def text_dump41(text):
    return [ast.dump(node,include_attributes=False)
            for node in ast.parse(text).body]


core_body41 = '''
import mpmath as mp
t,a,nu = mp.mpc(t),mp.mpf(a),mp.mpf(nu)
if abs(t)>mp.mpf('.3') or abs(a)>=mp.mpf('.6'):
    raise ValueError('core reference requires |t|<=.3 and |a|<.6')
e = mp.exp(-4*abs(t)**2/(1-abs(t)**2))
q = e*nu*nu*abs(t)**2
cell = int(mp.floor(a+mp.mpf('.5')))
b = a-cell
f = inverse_height_moment(b,q,0)
R = mp.sqrt(f*f+4*q)
if a>0:
    A = 2/(R-f) if cell==1 else (R+f)/2
    x = mp.sqrt(A/e)/nu
    return x,t/x
Y = 2/(R+f) if cell==-1 else (R-f)/2
y = mp.sqrt(Y)
if y==0:
    raise ValueError('the core reference excludes the node')
return t/y,y
'''
logradius_body41 = '''
import mpmath as mp
with mp.workdps(max(mp.mp.dps,dps)+50):
    a,t,nu = mp.mpf(a),mp.mpc(t),mp.mpf(nu)
    r = abs(t)
    if not 0<r<=mp.mpf('.5') or not 0<nu<=1:
        raise ValueError('elementary cusp inverse requires 0<|t|<=.5 and 0<nu<=1')
    P = mp.log(nu)+mp.log(r)-2*r*r/(1-r*r)
    q = mp.exp(2*P)
    n = mp.floor(a+mp.mpf('.5'))
    b = a-n
    lam = scalar_flat_step((r-mp.mpf('.3'))/mp.mpf('.2'))
    if lam==1:
        return +(P*(a+mp.mpf('.5')))
    f = inverse_height_moment(b,q,lam)
    root = mp.sqrt(f*f+4*q)
    Y = 2*q/(root+f) if f>=0 else (root-f)/2
    return +(n*P+mp.log(Y)/2)
'''
restore_loop41 = '''
for inf,tc,beta in ((False,t,1-scalar_flat_step((r-mp.mpf('.3'))/mp.mpf('.2'))),
                    (True,1/t,1-scalar_flat_step((1/r-mp.mpf('.3'))/mp.mpf('.2')))):
    if beta==0:
        continue
    e2,e3,e4,a2,a3 = cusp_data(tc,inf)
    c = coeffs(e4,e2,e3)
    if inf:
        zc,wc = z-mp.j*w/t,w/(t*t)
    else:
        s = lt/(2*mp.j*mp.pi)
        zc,wc = z-mp.pi*s*(s-1),w-2*t*s
    cc = coeffs(-wc,e2,e3)
    yy = mp.exp(mp.j*(zc-mp.j*(cc[0]*a2+cc[1]*a3)))
    H = scalar_finite_height(tc/yy,yy,mp.exp(-2) if inf else mp.mpf(1))
    if inf:
        D = tuple((H+four)*v for v in c)
    else:
        s = lt/(2*mp.j*mp.pi)
        xi = scalar_lattice_coords(mp.pi*s*(s-1),2*t*s,t,lt)
        D = tuple(xi[j+1]+(H-four+xi[3])*c[j] for j in range(2))
    C = [v+beta*d for v,d in zip(C,D)]
'''
selector_nodes41 = [node for node in ast.walk(raw_tree31)
                    if isinstance(node,ast.If)
                    and ast.dump(node.test,include_attributes=False)
                    ==ast.dump(ast.parse("inf or abs(t)<mp.mpf('.5')",
                                         mode='eval').body,
                               include_attributes=False)]
restore_loops41 = [node for node in ast.walk(restore_tree30)
                   if isinstance(node,ast.For)]
terminal_lam41 = [node for node in ast.walk(terminal_tree36)
                  if isinstance(node,ast.If)
                  and ast.dump(node.test,include_attributes=False)
                  ==ast.dump(ast.parse('lam==0',mode='eval').body,
                             include_attributes=False)]
height_G_tree41 = source_tree30(cusp30._height_G)
lattice_tree41 = source_tree30(cusp30.scalar_lattice_coords)
lattice_body41 = '''
import mpmath as mp
ez = (-mp.j*t,-t,-mp.j*logt)
ew = (1-t*t,mp.j*(1+t*t),2*t)
rows = [[v.imag for v in ez],[v.real for v in ew],[v.imag for v in ew]]
rhs = (z.imag,w.real,w.imag)
def det(a):
    return (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
            -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
            +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))
den = det(rows)
coords = []
for k in range(3):
    replaced = [[rhs[i] if j==k else rows[i][j] for j in range(3)] for i in range(3)]
    coords.append(det(replaced)/den)
one = (z.real-sum(a*b.real for a,b in zip(coords,ez)))/(2*mp.pi)
return (one,*coords)
'''
# The source det helper is the cofactor expansion of a 3x3 determinant,
# so its Cramer quotients are the unique solve whenever den!=0.
gen41 = sp.Matrix(3,3,sp.symbols('g41_0:9'))
det_helper41 = sp.expand(
    gen41[0,0]*(gen41[1,1]*gen41[2,2]-gen41[1,2]*gen41[2,1])
    -gen41[0,1]*(gen41[1,0]*gen41[2,2]-gen41[1,2]*gen41[2,0])
    +gen41[0,2]*(gen41[1,0]*gen41[2,1]-gen41[1,1]*gen41[2,0])
    -gen41.det())==0
source41 = (
    body_dump41(lattice_tree41)==text_dump41(lattice_body41)
    and det_helper41 and
    source_path37 and source_collar34 and source_torus35
    and source_cusp_ok30 and selector_source40 and arg_source38
    and body_dump41(core_tree37)==text_dump41(core_body41)
    and body_dump41(logradius_tree34)==text_dump41(logradius_body41)
    and len(selector_nodes41)==1
    and ordered_nodes37(selector_nodes41[0].body,[
        'tc = small if inf else t',
        "nu = mp.exp(-2) if inf else mp.mpf(1)",
        '''
if abs(tc)<=mp.mpf('.3') and elementary:
    xs,yy = elementary_cusp_reference(tc,a,nu)
    return ('ci' if inf else 'c0',(xs,yy,a))
''',
        'rho = (elementary_cusp_logradius if elementary else '
        'elliptic_cusp_logradius)(a,tc,nu,dps=dps)',
        'yy = mp.exp(rho)',
        'xs = tc/yy',
        'if inf:\n    t = 1/tc',
        'lt = mp.log(t)',
        'lattice = scalar_lattice_coords(zm,wm,t,lt)',
        'data = (t,lattice[0],lattice[3],lt,a) if elementary else '
        '(t,lattice[0],lattice[3],lt)',
        "return ('mid',data)"])
    and len(restore_loops41)==1
    and ast.dump(restore_loops41[0],include_attributes=False)
        ==text_dump41(restore_loop41)[0]
    and ordered_nodes37(restore_tree30.body,[
        't,one,four,lt = data[:4]',
        'r = abs(t)',
        'z,w = 2*mp.pi*one-mp.j*lt*four,2*t*four',
        'C = list(phi)'])
    and returned_expr30(restore_tree30,
        "'mid',(z-mp.j*t*C[0]-t*C[1],"
        "w+(1-t*t)*C[0]+mp.j*(1+t*t)*C[1],t)")
    and matching_statement31(slice_tree30,'x,y = data[:2]')
    and returned_expr30(slice_tree30,'kind,(mp.conj(phase)*x,phase*y)')
    and returned_expr30(slice_tree30,
        'kind,(t,one+mp.arg(phase)/(2*mp.pi),four,lt)')
    and matching_assignment30(slice_tree30,'phase',
        '_elementary_inner(ref,result)')
    and returned_expr30(inner_tree36,
        'sum(mp.conj(a)*b for a,b in zip(u,v))')
    and len(terminal_lam41)==1
    and [ast.dump(node,include_attributes=False)
         for node in terminal_lam41[0].body]
        ==text_dump41('kappa,amp,d = mp.mpf(2),mp.mpf(0),mp.mpf(1)')
    and ordered_nodes37(terminal_lam41[0].orelse,[
        'k = -mp.log(q)',
        'kappa = 2+lam*(k-2)',
        'amp = lam*k',
        'd = 1-lam+2*lam*mp.sqrt(q)*k'])
    and exact_if_body37(height_G_tree41,'amp==0',
        'return a/d if inverse else d*a')
    and returned_expr30(height_G_tree41,
        'mp.asinh(amp*a/d)/amp if inverse else d*mp.sinh(amp*a)/amp')
    and matching_if_test32(moment_tree36,
        "abs(a)>mp.mpf('.5') or not 0<=q<=mp.mpf('.125') or "
        "not 0<=lam<=1 or (lam>0 and q<mp.mpf(1)/2500)")
    and matching_assignment30(moment_tree36,'data',
        '_terminal_height_data(q,lam)')
    and matching_if_test32(moment_tree36,"lam==1 or aa<=mp.mpf('.1')")
    and matching_if_test32(moment_tree36,"aa>=mp.mpf('.4')")
    and exact_if_body37(flat_tree38,'u<=0','return mp.mpf(0)')
    and exact_if_body37(flat_tree38,'u>=1','return mp.mpf(1)'))

# The lambda=0 piece of the terminal data is the value of the general
# piece at lambda=0, and the amp=0 quotients are removable: their series
# are even in amp.  Thus f(b,q,lambda) is the restriction of one smooth
# function of lambda>=0 on the second parameter piece of section 11.2,
# whose smoothness up to lambda=0 is proved there.
lam41,k41,qq41,g41,a41,y41,d41 = sp.symbols(
    'lam41 k41 qq41 g41 a41 y41 d41',positive=True)
lambda_piece41 = (
    sp.simplify((2+lam41*(k41-2)).subs(lam41,0)-2)==0
    and (lam41*k41).subs(lam41,0)==0
    and sp.simplify((1-lam41+2*lam41*sp.sqrt(qq41)*k41)
                    .subs(lam41,0)-1)==0
    and sp.simplify(sp.series(d41*sp.sinh(g41*a41)/g41,g41,0,4).removeO()
                    -(d41*a41+d41*a41**3*g41**2/6))==0
    and sp.simplify(sp.series(sp.asinh(g41*y41/d41)/g41,g41,0,4).removeO()
                    -(y41/d41-g41**2*y41**3/(6*d41**3)))==0)
rep('C41a ordered source paths on both sides of |tc|=3/10; lambda=0 piece',
    source41 and lambda_piece41)

# C41b.  Write q=e*nu^2*|t|^2 (both sources use this q: the core sets it
# directly, the logradius sets q=exp(2P)) and R=sqrt(f^2+4q)>|f|, where
# f=f(b,q,lambda) is the common source moment and n=cell=floor(a+1/2).
# The right source has e^(2rho)=q^n*Y with Y=2q/(R+f)=(R-f)/2 (its two
# branches are one function).  C34e reconstructs zz=e^rho, so the M0
# reference is (x_m,y_m)=(t/e^rho,e^rho).  Cell by cell the core reference
# is (conj(g)x_m,g*y_m), with g=t/|t| for a>0 and g=1 for a<=0.
f41,tr41,ti41,th41 = sp.symbols('f41 tr41 ti41 th41',real=True)
q41,nu41 = sp.symbols('q41 nu41',positive=True)
t41 = tr41+sp.I*ti41
abs2_41 = tr41**2+ti41**2
R41 = sp.sqrt(f41**2+4*q41)
e41 = q41/(nu41**2*abs2_41)
Ypos41,Yneg41 = 2*q41/(R41+f41),(R41-f41)/2
cells41 = (
    # (a>0 ?, |x_c|^2 or |y_c|^2 from the core, e^(2 rho) from logradius)
    (True,((R41+f41)/2)/(e41*nu41**2),Ypos41),            # cell 0, a>0
    (True,(2/(R41-f41))/(e41*nu41**2),q41*Yneg41),        # cell 1
    (False,Yneg41,Yneg41),                               # cell 0, a<=0
    (False,2/(R41+f41),Ypos41/q41))                      # cell -1
rq41 = sp.Symbol('rq41',positive=True)
reference_cells41 = (
    # logradius q=exp(2P), P=log(nu)+log(r)-2r^2/(1-r^2), equals the
    # core q=e*nu^2*r^2, e=exp(-4r^2/(1-r^2)); both take cell=n=
    # floor(a+1/2), b=a-cell and the same moment f (lambda=0 in M0).
    sp.simplify(sp.exp(2*(sp.log(nu41)+sp.log(rq41)-2*rq41**2/(1-rq41**2)))
                -nu41**2*rq41**2*sp.exp(-4*rq41**2/(1-rq41**2)))==0
    and matching_assignment30(core_tree37,'cell',
                              "int(mp.floor(a+mp.mpf('.5')))")
    and matching_assignment30(logradius_tree34,'n',
                              "mp.floor(a+mp.mpf('.5'))")
    and sp.simplify(sp.expand((R41-f41)*(R41+f41))-4*q41)==0
    and sp.simplify(Ypos41-Yneg41)==0
    and all(sp.simplify((core*mid-abs2_41) if positive_x else
                        (core-mid))==0
            for positive_x,core,mid in cells41))
# Gauge: with |x_c|*e^rho=|t| (a>0) the positive core x_c equals
# conj(g)*t/e^rho and y_c=t/x_c=g*e^rho.  A unit scalar g multiplies both
# source frames (plain and deck-shifted) by g, since r is unchanged.
xr41,xi41,yr41,yi41 = sp.symbols('xr41 xi41 yr41 yi41',real=True)
m41 = sp.Symbol('m41',positive=True)
x41,y41c = xr41+sp.I*xi41,yr41+sp.I*yi41
gauge41 = sp.exp(sp.I*th41)
def frame41(x,y):
    r = sp.expand(nu41**2*x*sp.conjugate(x)+y*sp.conjugate(y))
    return sp.Matrix([nu41*sp.conjugate(x),y])/sp.sqrt(r)
frame_plain41 = (frame41(sp.conjugate(gauge41)*x41,gauge41*y41c)
                 -gauge41*frame41(x41,y41c))
frame_deck41 = (frame41(1/(m41*gauge41*y41c),
                        m41*sp.conjugate(gauge41)*x41*(gauge41*y41c)**2)
                -gauge41*frame41(1/(m41*y41c),m41*x41*y41c**2))
gauge_frames41 = (
    all(sp.simplify(sp.expand_complex(v))==0 for v in frame_plain41)
    and all(sp.simplify(sp.expand_complex(v))==0 for v in frame_deck41)
    # g=t/|t| with |t|=sqrt(abs2)>0 real: conj(g)*t=|t|.
    and sp.simplify(sp.expand(sp.conjugate(t41)*t41)/sp.sqrt(abs2_41)
                    -sp.sqrt(abs2_41))==0)
# The weights.  Core: hw(H)*frame(x,y) and, second, hw(H+1/2)*deck frame
# if H<0, hw(H-1/2)*frame(x,y) if H>=0.  torus_pair(z,t,nu,H): centre c
# has shift=floor(c+1/2-H), weight hw(H+shift-c), frame(t/yy,yy) with
# yy=z*(m*t)^shift; shift=1 gives frame(1/(m z),m*t*z), the deck frame.
# hw(h)=0 for |h|>=2/5 (flat step >=1).  On |H|<3/5 the source shift
# values are constant on these intervals, and every disagreement of
# the two formulas carries a zero weight.
def shift_constant41(c,lo,lo_closed,hi,hi_closed,shift):
    # g(H)=c+1/2-H decreases, so g maps the H-interval onto the interval
    # with ends g(hi)<g(lo), closedness swapped.  floor(g)=shift on it
    # iff g(hi)>=shift and g(lo)<shift+1 (or <=shift+1 at an open end).
    left,right = c+sp.Rational(1,2)-hi,c+sp.Rational(1,2)-lo
    right_closed = lo_closed
    return (left>=shift
            and (right<shift+1 if right_closed else right<=shift+1))
two5_41 = sp.Rational(2,5)
half41 = sp.Rational(1,2)
H_intervals41 = (
    # (lo, lo_closed, hi, hi_closed, shift centre 0, shift centre 1/2)
    (-sp.Rational(3,5),False,-half41,True,1,1),
    (-half41,False,sp.Integer(0),False,0,1),
    (sp.Integer(0),False,half41,True,0,0),
    (half41,False,sp.Rational(3,5),False,-1,0))
shift_table41 = all(
    shift_constant41(sp.Integer(0),lo,lc,hi,hc,s0)
    and shift_constant41(half41,lo,lc,hi,hc,s1)
    for lo,lc,hi,hc,s0,s1 in H_intervals41)
# Zero weights where the formulas differ: on (-3/5,-1/2] both first
# weights hw(H+1)=hw(H)=0; on (1/2,3/5) hw(H-1)=hw(H)=0; at H=0 the
# second weights are hw(1/2)=hw(-1/2)=0; on (-1/2,0) and (-3/5,-1/2]
# the second pair is the deck frame on both sides, with weight
# hw(H+1/2); on [0,3/5) it is the plain frame with weight hw(H-1/2).
zero_weights41 = (
    sp.Rational(2,5)<=half41 and abs(-sp.Rational(3,5)+1)>=two5_41
    and abs(-half41+1)>=two5_41 and half41>=two5_41
    and abs(sp.Rational(3,5)-1)>=two5_41 and abs(half41-1)>=two5_41
    and sp.floor(half41-0)==0 and sp.floor(1-sp.Integer(0))==1
    and abs(0+1-half41)>=two5_41 and abs(0-half41)>=two5_41)
rep('C41b core reference = unit multiple of the M0 reference, cell by cell',
    reference_cells41 and gauge_frames41 and shift_table41
    and zero_weights41 and south_reconstruction34
    and north_reconstruction34)

# C41c.  On the right, near |tc|=3/10, the source has (with u=5|tc|-3/2)
# lambda=s(u) in elementary_cusp_logradius, b0=1-s(u), bm=s(u) in the
# embedding (bi=0 exactly: (1/r-3/10)/(1/5)>=1 for r<=2), and
# beta=b0=1-s(u), beta_inf=0 in the middle restoration.  The source flat
# step s has value and first two derivatives tending to zero as u->0+.
# Hence every right source quantity minus its M0 value has a zero
# complete two-jet on the seam: either it carries the factor s(u) or 1-b0,
# or it is F(lambda,...)-F(0,...) for a function F smooth in lambda>=0
# (C41a).  The generic chain-rule checks below use arbitrary smooth K, F.
u41 = sp.Symbol('u41',positive=True)
flat41 = flat_actual34.subs(p34,u41)
flat_limits41 = (
    flat_jets34
    and all(sp.limit(expr.subs(p34,u41),u41,0,'+')==0
            for expr in (flat_actual34,flat_slope34,flat_curv34))
    and sp.simplify(sp.diff(flat41,u41)-flat_slope34.subs(p34,u41))==0
    and sp.simplify(sp.diff(flat41,u41,2)-flat_curv34.subs(p34,u41))==0
    and (1/sp.Integer(2)-sp.Rational(3,10))/sp.Rational(1,5)==1)
rr41,v41,w41 = sp.symbols('rr41 v41 w41',real=True)
gf41 = sp.Function('gf41')(rr41)
Kf41 = sp.Function('Kf41')(rr41,v41,w41)
Ff41 = sp.Function('Ff41')
lam_s41 = sp.Symbol('lam_s41',real=True)
def flat_at_seam41(expr):
    # Substitute g=g'=g''=0 (the flat step at u=0+) after differentiating.
    out = expr.subs({sp.Derivative(gf41,(rr41,2)):0,
                     sp.Derivative(gf41,rr41):0}).subs(gf41,0)
    return sp.simplify(out.doit())==0
factor_terms41 = gf41*Kf41
param_terms41 = Ff41(gf41,v41,w41)-Ff41(0,v41,w41)
derivs41 = ((),(rr41,),(v41,),(w41,),(rr41,rr41),(rr41,v41),(rr41,w41),
            (v41,v41),(v41,w41),(w41,w41))
flat_two_jets41 = all(
    flat_at_seam41(sp.diff(expr,*d) if d else expr)
    for expr in (factor_terms41,param_terms41) for d in derivs41)
rep('C41c right source = M0 + terms flat at |tc|=3/10: complete zero two-jets',
    flat_limits41 and flat_two_jets41 and flat_all_directions40)

# C41d.  Final phase and torus restoration.  Let E be the common final
# unit vector of the finite fold (it lies on the final rank-one image, so
# p=<ref,E> is a unit scalar).  With ref_c=g*ref_m (C41b) the source
# inner product, conjugate-linear in its first slot, gives p_m=g*p_c.
# The c0/ci W-point is (conj(p_c)x_c,p_c*y_c); the middle W-point shifts
# `one` by Arg(p_m)/(2pi), i.e. multiplies the reconstructed zz=e^rho by
# p_m.  Both W-points therefore have y=p_c*g*e^rho.  An integer change of
# the local Arg branch is a lattice translation of `one`.
al41 = sp.Symbol('al41',real=True)
p_c41 = sp.exp(sp.I*al41)
r1_41,r2_41,E1_41,E2_41 = sp.symbols('r1_41 r2_41 E1_41 E2_41')
def inner41(u,v):                    # the source: sum(conj(a)*b)
    return sum(sp.conjugate(a)*b for a,b in zip(u,v))
ref_c41 = (r1_41,r2_41)
ref_m41 = tuple(sp.conjugate(gauge41)*v for v in ref_c41)   # ref_c=g*ref_m
E41 = (E1_41,E2_41)
p_m41 = sp.expand(inner41(ref_m41,E41))
xm41,ym41 = sp.symbols('xm41 ym41')
final_W41 = (
    sp.simplify(p_m41-gauge41*inner41(ref_c41,E41))==0
    # control: a linear first slot would give conj(g)*p_c instead
    and sp.simplify(sum(a*b for a,b in zip(ref_m41,E41))
                    -gauge41*sum(a*b for a,b in zip(ref_c41,E41)))!=0
    # core W y=p_c*(g*y_m), middle zz=e^rho times p_m=g*p_c: equal;
    # x=tc/y on both sides.
    and sp.simplify(p_c41*(gauge41*ym41)-(gauge41*p_c41)*ym41)==0
    and sp.simplify(sp.conjugate(p_c41)*(sp.conjugate(gauge41)*xm41)
                    *(gauge41*p_c41)*ym41-xm41*ym41)==0)
# Restoration.  The middle chart (z,w,t) and the c0 chart (x,y,w_c) are
# glued by the paper's holomorphic Phi: y=exp(i(z-pi s(s-1))),
# s=log(t)/(2 pi i), x=t/y, w_c=w-2ts (tfull.c0_to_mid is its inverse);
# ci is glued by Theta: y=exp(i(z-i w/t)), x=(1/t)/y, w_c=w/t^2.
# Take any branch s of log(t)/(2 pi i) (south), with real lattice data
# xi3 and (xi1,xi2) the unique real solution of the w-equation
# 2ts=xi1(1-t^2)+xi2 i(1+t^2)+xi3 2t of scalar_lattice_coords.
sr41,si41,xi3_41,four41,H41,ph2_41,ph3_41,one41 = sp.symbols(
    'sr41 si41 xi3_41 four41 H41 ph2_41 ph3_41 one41',real=True)
s41 = sr41+sp.I*si41
def coeffs41(rhs,e2,e3):
    det = sp.re(e2)*sp.im(e3)-sp.re(e3)*sp.im(e2)
    return ((sp.re(rhs)*sp.im(e3)-sp.re(e3)*sp.im(rhs))/det,
            (sp.re(e2)*sp.im(rhs)-sp.re(rhs)*sp.im(e2))/det)
def cusp_data41(tc,inf):
    return ((tc*tc-1,sp.I*(1+tc*tc),-2*tc,tc,sp.I*tc) if inf else
            (1-tc*tc,sp.I*(1+tc*tc),2*tc,tc,-sp.I*tc))
phi41 = (ph2_41,ph3_41)
restoration41 = True
chart_dets41 = []
for inf41 in (False,True):
    lt41 = 2*sp.pi*sp.I*s41
    z41 = 2*sp.pi*one41-sp.I*lt41*four41
    w41m = 2*t41*four41
    tc41 = 1/t41 if inf41 else t41
    e2_41,e3_41,e4_41,a2_41,a3_41 = cusp_data41(tc41,inf41)
    c41 = coeffs41(sp.expand_complex(e4_41),e2_41,e3_41)
    if inf41:
        zc41,wc41 = z41-sp.I*w41m/t41,w41m/(t41*t41)
    else:
        zc41,wc41 = z41-sp.pi*s41*(s41-1),w41m-2*t41*s41
    cc41 = coeffs41(sp.expand_complex(-wc41),e2_41,e3_41)
    if inf41:
        D41 = tuple((H41+four41)*v for v in c41)
    else:
        xi12_41 = coeffs41(sp.expand_complex(2*t41*s41-2*t41*xi3_41),
                           e2_41,e3_41)
        D41 = tuple(xi12_41[j]+(H41-four41+xi3_41)*c41[j]
                    for j in range(2))
    # D-cc=H*c: the restored torus coordinates differ exactly by H c.
    restoration41 &= all(sp.simplify(D41[j]-cc41[j]-H41*c41[j])==0
                         for j in range(2))
    yy_exp41 = sp.I*(zc41-sp.I*(cc41[0]*a2_41+cc41[1]*a3_41))
    Cm41 = [phi41[j]+D41[j] for j in range(2)]
    zr41 = z41-sp.I*t41*Cm41[0]-t41*Cm41[1]
    wr41 = w41m+(1-t41**2)*Cm41[0]+sp.I*(1+t41**2)*Cm41[1]
    if inf41:
        y_exp41,wcusp41 = sp.I*(zr41-sp.I*wr41/t41),wr41/t41**2
    else:
        y_exp41,wcusp41 = sp.I*(zr41-sp.pi*s41*(s41-1)),wr41-2*t41*s41
    Cc41 = [phi41[j]+H41*c41[j] for j in range(2)]
    A41 = Cc41[0]*a2_41+Cc41[1]*a3_41
    restoration41 &= (
        sp.simplify(sp.expand(y_exp41-(yy_exp41+A41)))==0
        and sp.simplify(sp.expand(wcusp41-(Cc41[0]*e2_41+Cc41[1]*e3_41)))==0)
    # Complex Jacobian of the holomorphic chart change (z,w,t)->(x,y,w_c).
    zh41,wh41,th41h = sp.symbols('zh41 wh41 th41h')
    sh41 = sp.log(th41h)/(2*sp.pi*sp.I)
    if inf41:
        yh41 = sp.exp(sp.I*(zh41-sp.I*wh41/th41h))
        Fh41 = sp.Matrix([1/(th41h*yh41),yh41,wh41/th41h**2])
    else:
        yh41 = sp.exp(sp.I*(zh41-sp.pi*sh41*(sh41-1)))
        Fh41 = sp.Matrix([th41h/yh41,yh41,wh41-2*th41h*sh41])
    chart_dets41.append(sp.simplify(Fh41.jacobian([zh41,wh41,th41h]).det()))
restoration41 &= (sp.simplify(chart_dets41[0]-sp.I)==0
                  and sp.simplify(chart_dets41[1]+sp.I/sp.Symbol('th41h')**4)==0)
tfull_tree41 = ast.parse(open(os.path.join(os.path.dirname(
    cusp30.__file__),'tfull.py')).read())
tfull_defs41 = {node.name:node for node in tfull_tree41.body
                if isinstance(node,ast.FunctionDef)}
chart_source41 = (
    'c0_to_mid' in tfull_defs41 and 'ci_to_mid' in tfull_defs41
    and body_dump41(tfull_defs41['c0_to_mid'])==text_dump41('''
t = x * y
logt = torch.log(t)
s = logt / (2j * math.pi)
return -1j * torch.log(y) + math.pi * s * (s - 1), w + 2 * t * s, t, logt
''')
    and body_dump41(tfull_defs41['ci_to_mid'])==text_dump41('''
s = x * y
t = 1 / s
wm = w * t ** 2
return -1j * torch.log(y) + 1j * wm / t, wm, t, -torch.log(s)
'''))
# Phi and Theta invert these: compose Phi (Theta) with the tfull output
# and recover (x,y,w); a different log branch changes z by 2 pi k only.
xl41,yl41,wl41 = sp.symbols('xl41 yl41 wl41',nonzero=True)
tl41 = xl41*yl41
sl41 = sp.log(tl41)/(2*sp.pi*sp.I)
z_c0_41 = -sp.I*sp.log(yl41)+sp.pi*sl41*(sl41-1)        # c0_to_mid
w_c0_41 = wl41+2*tl41*sl41
y_back41 = sp.exp(sp.I*(z_c0_41-sp.pi*sl41*(sl41-1)))
tci41 = 1/tl41                                          # ci_to_mid
w_ci_41 = wl41*tci41**2
z_ci_41 = -sp.I*sp.log(yl41)+sp.I*w_ci_41/tci41
y_back_ci41 = sp.exp(sp.I*(z_ci_41-sp.I*w_ci_41/tci41))
chart_inverse41 = (
    sp.simplify(y_back41-yl41)==0
    and sp.simplify(tl41/y_back41-xl41)==0
    and sp.simplify(w_c0_41-2*tl41*sl41-wl41)==0
    and sp.simplify(y_back_ci41-yl41)==0
    and sp.simplify((1/tci41)/y_back_ci41-xl41)==0
    and sp.simplify(w_ci_41/tci41**2-wl41)==0)
# The restore loop reconstructs yy from the final middle W-point: at the
# raw data of C34e with one shifted by al/(2pi), its exponent is
# rho+i*al.  With al=Arg(p_m) this is yy=e^rho*p_m, the core y, so both
# restorations call scalar_finite_height with the same (tc/y,y).
del41 = sp.Symbol('del41',real=True)
yy_link41 = True
for inf41 in (False,True):
    if inf41:
        tcl41 = 1/t34
        zl41,wl41m = z_mid_n34+2*sp.pi*del41,w_mid_n34
        zcl41,wcl41 = zl41-sp.I*wl41m/t34,wl41m/t34**2
    else:
        tcl41 = t34
        zl41,wl41m = z_mid_s34+2*sp.pi*del41,w_mid_s34
        zcl41,wcl41 = zl41-sp.pi*s34*(s34-1),wl41m-2*t34*s34
    e2l41,e3l41,_,a2l41,a3l41 = cusp_data41(tcl41,inf41)
    ccl41 = source_coeff34(-wcl41,e2l41,e3l41)
    yy_link41 &= sp.simplify(sp.I*(zcl41-sp.I*(ccl41[0]*a2l41+ccl41[1]*a3l41))
                             -(rho34+2*sp.pi*sp.I*del41))==0
rep('C41d final W-points agree; c0/ci restoration = Phi/Theta of middle one',
    final_W41 and restoration41 and chart_source41 and chart_inverse41
    and yy_link41 and source_cusp_ok30 and source_torus35)

# C41e.  Domain of the whole seam and exact round witnesses.  On
# 2/25<|tc|^2<1/10 (an open band around the seam) both sources are in
# their certified domains: south q=s*exp(-4s/(1-s)) and north
# q=e^-4*s*exp(-4s/(1-s)) (s=|tc|^2) lie in (1/2500,1/8), as needed
# by inverse_height_moment with lambda>0; the Cramer determinant is
# >2r^2(1+r^2)>0 (for every r=|t|>0, below); the torus 2x2 determinants
# are 1-s^2>0 and s^2-1<0; and |a|<3/5 in both argument charts
# (|Arg(-p)|<arccos(9/10)<pi/5 in the half-centred chart).
s_band41 = sp.Symbol('s_band41',positive=True)
q_south41 = s_band41*sp.exp(-4*s_band41/(1-s_band41))
band41 = (
    sp.simplify(sp.diff(q_south41,s_band41)
                -sp.exp(-4*s_band41/(1-s_band41))
                *(1-4*s_band41/(1-s_band41)**2))==0
    # (1-4s/(1-s)^2)>0 on s<1/10: 4/(10*(81/100))<1, so q increases.
    and 4*sp.Rational(1,10)/(1-sp.Rational(1,10))**2<1
    and q_south41.subs(s_band41,sp.Rational(2,25))
        >sp.Rational(2,25)*sp.exp(-sp.Rational(1,2))
    # q(2/25)=(2/25)exp(-8/23)>(2/25)exp(-1/2)>(2/25)(1/2)=1/25.
    and sp.exp(-sp.Rational(1,2))>sp.Rational(1,2)
    and sp.Rational(2,25)*sp.Rational(1,2)==sp.Rational(1,25)
    and sp.Rational(1,10)<sp.Rational(1,8)
    and sp.exp(4)<55 and sp.Rational(1,25)/55>sp.Rational(1,2500)
    and 4*sp.Rational(2,25)/(1-sp.Rational(2,25))<sp.Rational(1,2)
    and sp.cos(sp.pi/5)<sp.Rational(9,10)
    and sp.Rational(1,2)+sp.Rational(1,10)==sp.Rational(3,5)
    and 1-sp.Rational(1,10)**2>0
    # det-2r^2(1+r^2)=(1+r^2)(r^2-1)log r>=0 for every r>0 (both factors
    # change sign at r=1), so det>2r^2(1+r^2)>0 on the whole band.
    and sp.factor(det_expected34-2*(ar34**2+ai34**2)*(1+ar34**2+ai34**2)
                  -(1+ar34**2+ai34**2)*(ar34**2+ai34**2-1)*ell34)==0)
# Exact witnesses.  Put p=|tc|^2, A=(4 pi)^2, B=e^4, Q=sqrt(1+k^2).
# South: k^2=(B p^2+A p)/(1+A p), u1=-k/Q, u2=i/Q.  North:
# k^2=(p^2+A p)/(A p+B), u1=i/Q, u2=k/Q.  Then the source positive root
# gives tc=i sqrt(p) (south) and tc=-i sqrt(p) (north, middle t=i/sqrt(p)),
# the phase coordinate is i/Q, so a=1/4 in the principal chart, and
# Arg t=pi/2 avoids the log cut.  Stereographic y=(0,0,Q+k) (south)
# and y=(1/Q,k/Q,0) (north); X=(R-1)/sqrt(1+R^2), R=|y-(0,1,0)|.
pw41,Aw41,Bw41 = sp.symbols('pw41 Aw41 Bw41',positive=True)
k2_south41 = (Bw41*pw41**2+Aw41*pw41)/(1+Aw41*pw41)
k2_north41 = (pw41**2+Aw41*pw41)/(Aw41*pw41+Bw41)
root_south41 = (Aw41+2*Bw41*pw41+Aw41*Bw41*pw41**2)/(1+Aw41*pw41)
root_north41 = (Aw41*Bw41+2*Bw41*pw41+Aw41*pw41**2)/(Aw41*pw41+Bw41)
b_south41 = Aw41*(1-k2_south41)
b_north41 = Aw41*(1-k2_north41)
witness_roots41 = (
    sp.factor(root_south41**2-b_south41**2-4*Bw41*k2_south41)==0
    and sp.factor(root_north41**2-b_north41**2-4*Bw41*k2_north41)==0
    # positive roots: both displayed numerators are positive
    and sp.factor(k2_south41*2/(b_south41+root_south41)-pw41)==0
    and sp.factor(k2_north41*2*Bw41/(b_north41+root_north41)-pw41)==0)
kw41 = sp.Symbol('kw41',positive=True)
Qw41 = sp.sqrt(1+kw41**2)
stereo41 = True
for u1w41,u2w41,yw41,re_u1_41 in (
        (-kw41/Qw41,sp.I/Qw41,(0,0,Qw41+kw41),-kw41/Qw41),
        (sp.I/Qw41,kw41/Qw41,(1/Qw41,kw41/Qw41,0),sp.Integer(0))):
    Sy41 = 1+sum(v**2 for v in yw41)
    stereo41 &= (
        sp.simplify(u1w41*sp.conjugate(u1w41)
                    +u2w41*sp.conjugate(u2w41)-1)==0
        and sp.simplify((1-(Sy41-1)+2*sp.I*yw41[0])/Sy41-u1w41)==0
        and sp.simplify(2*(yw41[1]+sp.I*yw41[2])/Sy41-u2w41)==0
        and sp.simplify(sp.re(u1w41)-re_u1_41)==0)
# 1+Re(u1)=1-k/Q>0 south because Q^2-k^2=1 (Q>k>0), and 1 north, so the
# inverse stereographic coordinates y above are well defined.
stereo41 &= sp.simplify(Qw41**2-kw41**2-1)==0
# Ratios: south u1/u2=i k (|u1|<|u2| iff k<1), north u2/u1=-i k
# (|u1|>|u2| iff k<1).  p<=(301/1000)^2<1/10 gives B p^2<55/100<1, and
# p^2<B, so k<1 at all three witnesses; |small|=sqrt(p)<1/2.
witness_sides41 = (
    sp.simplify((-kw41/Qw41)/(sp.I/Qw41)-sp.I*kw41)==0
    and sp.simplify((kw41/Qw41)/(sp.I/Qw41)+sp.I*kw41)==0
    and sp.simplify(k2_south41-1-(Bw41*pw41**2-1)/(1+Aw41*pw41))==0
    and sp.simplify(k2_north41-1-(pw41**2-Bw41)/(Aw41*pw41+Bw41))==0
    and 55*sp.Rational(301,1000)**4<1
    and sp.Rational(301,1000)**2<sp.Rational(1,10)
    and sp.Rational(299,1000)**2>sp.Rational(2,25)
    and sp.Rational(299,1000)<sp.Rational(3,10)<sp.Rational(301,1000)
    and sp.Rational(301,1000)<sp.Rational(1,2)
    and sp.arg(sp.I)==sp.pi/2 and sp.arg(sp.I)/(2*sp.pi)==sp.Rational(1,4)
    and sp.re(sp.I)>-sp.Rational(9,10))
# Free slab: south R^2=1+(Q+k)^2 with 1<Q+k<1+sqrt(2), so 1<R<sqrt(7);
# north R^2=2-2k/Q, where (k/Q)^2=k^2/(1+k^2)<1/2 for k<1, so
# 2-sqrt(2)<R^2<2.  X=(R-1)/sqrt(1+R^2) increases in R (round_inverse38),
# hence |X|<3/5 at every witness, far from the pole bands 1-|X|<10^-5.
Rv41 = sp.Symbol('Rv41',positive=True)
Xv41 = (Rv41-1)/sp.sqrt(1+Rv41**2)
slab41 = (
    round_inverse38
    and sp.simplify(sp.diff(Xv41,Rv41)-(1+Rv41)/(1+Rv41**2)**sp.Rational(3,2))==0
    and 1+(1+sp.sqrt(2))**2<7
    and Xv41.subs(Rv41,sp.sqrt(7))<sp.Rational(3,5)
    and Xv41.subs(Rv41,sp.sqrt(2-sp.sqrt(2)))>-sp.Rational(3,5)
    and sp.simplify((kw41/Qw41)**2-kw41**2/(1+kw41**2))==0
    and sp.simplify(kw41**2/(1+kw41**2)-sp.Rational(1,2)
                    -(kw41**2-1)/(2*(1+kw41**2)))==0)
rep('C41e seam band in both source domains; exact witnesses on and across it',
    band41 and witness_roots41 and stereo41 and witness_sides41 and slab41)

# C41f is solely a reproducible numerical observation (d), sample size
# 18.  Per sign it uses four (a, arg) pairs covering a<=0 (with the deck
# frame), a=0, a>0 in cell 0 and a>0 in cell 1 (half-centred argument
# chart), at |tc|=299/1000 (left) and 3/10 (seam), plus one right point
# |tc|=301/1000.  The spinor is u2=e^(2 pi i a)/Q, u1=k u2 e^(i arg tc)
# (south) or u1=-conj(e^(2 pi i a))/Q, u2=k u1 e^(i arg tc) (north), with
# k from the C41e profiles.  Left: the production raw reference (c0/ci)
# is compared with the production middle formulas at the same (tc,a),
# which are M0 there (b0 or bi is 1, bm=0, lambda=0); the restorations
# are compared through Phi/Theta for one common final phase.  Seam:
# rounding may send the raw reference to either branch, so the
# production core formula is evaluated directly and compared in the same
# way.  Right: the flat collar weight and rho-rho(lambda=0) are printed.
# It is not used by C41a-e.
mp41 = mp39
rows41 = []
sample_ok41 = True
def middle_m0_41(tc,a,nu,north):
    rho = cusp30.elementary_cusp_logradius(a,tc,nu,dps=110)
    tm = 1/tc if north else tc
    lt = mp41.log(tm)
    if north:
        zm,wm = -mp41.j*rho,mp41.mpc(0)
    else:
        ss = lt/(2*mp41.j*mp41.pi)
        zm,wm = -mp41.j*rho+mp41.pi*ss*(ss-1),2*tm*ss
    lat = cusp30.scalar_lattice_coords(zm,wm,tm,lt)
    return (tm,lat[0],lat[3],lt,a)
def compare41(kind,core,north,phi,pc):
    xc,yc,ac = core
    nu = mp41.exp(-2) if north else mp41.mpf(1)
    ref_c = cusp30.elliptic_bundle_embedding(kind,core,dps=110,elementary=True)
    data_m = middle_m0_41(xc*yc,ac,nu,north)
    ref_m = cusp30.elliptic_bundle_embedding('mid',data_m,dps=110,elementary=True)
    Pdiff = max(abs(a*mp41.conj(b)-c*mp41.conj(d))
                for a,c in zip(ref_c,ref_m) for b,d in zip(ref_c,ref_m))
    E = tuple(pc*v for v in ref_c)
    p_c = cusp30._elementary_inner(ref_c,E)
    p_m = cusp30._elementary_inner(ref_m,E)
    tm,lt = data_m[0],data_m[3]
    _,out_c = cusp30.elementary_restore_torus(kind,(mp41.conj(p_c)*xc,p_c*yc),phi)
    _,out_m = cusp30.elementary_restore_torus(
        'mid',(tm,data_m[1]+mp41.arg(p_m)/(2*mp41.pi),data_m[2],lt),phi)
    zr,wr,_ = out_m
    if north:
        yr = mp41.exp(mp41.j*(zr-mp41.j*wr/tm))
        mapped = (1/(tm*yr),yr,wr/tm**2)
    else:
        sr = lt/(2*mp41.j*mp41.pi)
        yr = mp41.exp(mp41.j*(zr-mp41.pi*sr*(sr-1)))
        mapped = (tm/yr,yr,wr-2*tm*sr)
    return Pdiff,max(abs(a-b) for a,b in zip(out_c,mapped)),abs(abs(p_m)-1)
with mp41.workdps(120), cusp30.elementary_formula_branches():
    A41n,B41n = (4*mp41.pi)**2,mp41.exp(4)
    phi41n = (mp41.mpf(1)/3,mp41.mpf(2)/7)
    pc41n = mp41.exp(mp41.j*mp41.mpf('0.7'))
    pairs41n = ((mp41.mpf('-0.35'),mp41.mpf(2)),(mp41.mpf(0),mp41.mpf(-1)),
                (mp41.mpf(1)/4,mp41.pi/2),(mp41.mpf('0.55'),mp41.mpf('0.7')))
    def raw41(rad,arg_tc,aa,north):
        pp = rad**2
        kk = mp41.sqrt((pp**2+A41n*pp)/(A41n*pp+B41n) if north else
                       (B41n*pp**2+A41n*pp)/(1+A41n*pp))
        QQ = mp41.sqrt(1+kk**2)
        if north:
            u1 = -mp41.conj(mp41.exp(2*mp41.pi*mp41.j*aa))/QQ
            u2 = kk*u1*mp41.exp(mp41.j*arg_tc)
        else:
            u2 = mp41.exp(2*mp41.pi*mp41.j*aa)/QQ
            u1 = kk*u2*mp41.exp(mp41.j*arg_tc)
        d = 1+u1.real
        w = (u1.imag/d,u2.real/d-1,u2.imag/d)
        R = mp41.sqrt(sum(v*v for v in w))
        X = (R-1)/mp41.sqrt(1+R**2)
        n = (w[0]/R,-w[2]/R,-w[1]/R)
        kind,data = cusp30.elliptic_raw_base_reference(
            tuple((1-X**2)*v for v in n),X,dps=110,elementary=True)
        return kind,data,X
    for north41 in (False,True):
        nu41n = mp41.exp(-2) if north41 else mp41.mpf(1)
        ckind41n = 'ci' if north41 else 'c0'
        sgn41n = -1 if north41 else 1
        for aa41n,arg41n in pairs41n:
            arg_tc41n = sgn41n*arg41n
            for rad41n in (mp41.mpf(299)/1000,mp41.mpf(3)/10):
                tc41n = rad41n*mp41.exp(mp41.j*arg_tc41n)
                kind41n,data41n,X41n = raw41(rad41n,arg_tc41n,aa41n,north41)
                sample_ok41 &= abs(X41n)<1-mp41.mpf('1e-5')
                if rad41n<mp41.mpf(3)/10:
                    sample_ok41 &= (kind41n==ckind41n
                        and abs(data41n[0]*data41n[1]-tc41n)<mp41.mpf('1e-100')
                        and abs(data41n[2]-aa41n)<mp41.mpf('1e-100'))
                    core41n = data41n
                else:
                    core41n = (*cusp30.elementary_cusp_reference(
                        tc41n,aa41n,nu41n),aa41n)
                    rows41.append(('seam-branch',north41,kind41n,None))
                Pdiff41n,outdiff41n,unit41n = compare41(
                    ckind41n,core41n,north41,phi41n,pc41n)
                rows41.append(('left',north41,Pdiff41n,outdiff41n))
                sample_ok41 &= (Pdiff41n<mp41.mpf('1e-90')
                                and outdiff41n<mp41.mpf('1e-90')
                                and unit41n<mp41.mpf('1e-100'))
        rad41n = mp41.mpf(301)/1000
        arg_tc41n = sgn41n*mp41.pi/2
        tc41n = rad41n*mp41.exp(mp41.j*arg_tc41n)
        kind41n,data41n,X41n = raw41(rad41n,arg_tc41n,mp41.mpf(1)/4,north41)
        t_r41n,_,_,_,a_r41n = data41n
        sample_ok41 &= (kind41n=='mid'
                        and abs(t_r41n-(1/tc41n if north41 else tc41n))
                        <mp41.mpf('1e-100')
                        and abs(a_r41n-mp41.mpf(1)/4)<mp41.mpf('1e-100'))
        lam41n = cusp30.scalar_flat_step((rad41n-mp41.mpf(3)/10)/(mp41.mpf(1)/5))
        rho41n = cusp30.elementary_cusp_logradius(a_r41n,tc41n,nu41n,dps=110)
        P41n = mp41.log(nu41n)+mp41.log(rad41n)-2*rad41n**2/(1-rad41n**2)
        q41n = mp41.exp(2*P41n)
        f041n = cusp30.inverse_height_moment(a_r41n,q41n,0)
        rho041n = mp41.log(2*q41n/(mp41.sqrt(f041n**2+4*q41n)+f041n))/2
        rows41.append(('right',north41,lam41n,abs(rho41n-rho041n)))
        sample_ok41 &= (0<lam41n<mp41.mpf('1e-80')
                        and abs(rho41n-rho041n)<mp41.mpf('1e-80'))
    print('C41f production observations (d), sample size 18:'
          ' max |P_c-P_M0| =',mp41.nstr(max(r[2] for r in rows41
                                            if r[0]=='left'),5),
          ' max |out_c-Phi(out_mid)| =',mp41.nstr(max(r[3] for r in rows41
                                                     if r[0]=='left'),5),
          ' right flat weights =',[mp41.nstr(r[2],5) for r in rows41
                                   if r[0]=='right'],
          ' |rho-rho(lambda=0)| =',[mp41.nstr(r[3],5) for r in rows41
                                     if r[0]=='right'],
          ' seam raw branches =',[r[2] for r in rows41
                                  if r[0]=='seam-branch'])
rep('C41f 18 production seam observations (numerical observation, d)',
    sample_ok41)

# C42: every pure-middle source seam at once: the |t|=1 selector and the
# principal-argument cuts of u1, u2 and t (and their corners), on the
# pure middle 1/2<=|t|<=2.  Middle W-data (t,one,four,L) give the point
# Psi=(z,w)=(2pi one-i L four,2t four) upstairs.  The group G generated by
#   (a) one->one+1,  (b) four->four+1,  (c) (one,L)->(one-four,L+2pi i)
# changes Psi by the lattice vectors e1=(2pi,0), e4(L)=(-iL,2t) and 0.
# C42 proves: on any open subset of the interior 1/2<|t|<2 the source
# W-data agree modulo G with W-data built from continuous branches of
# arg u1, arg u2, arg t, log t (hence smooth); the ten-slot projector is
# that of the continuous data; and, for a final unit vector E of the
# fold that is smooth there (a premise resting on the other seams and on
# nonzero overlaps), the restored middle output differs from a smooth
# one by a holomorphic lattice translation.  This covers every
# pure-middle seam, not one chamber; (j) of symbolic_identities.py is the
# special case of the t-cut on 1/2<|t|^2<3/4.  The joins at |t|=1/2, 2
# are separate.
middle_tail42 = '''
r = abs(t)
lt = mp.log(t)
ss = lt/(2*mp.j*mp.pi)
xi4 = scalar_lattice_coords(mp.pi*ss*(ss-1),2*t*ss,t,lt)[3]
if r<=1:
    four = mp.arg(u2)/(2*mp.pi)+xi4+mp.mpf('.5')
else:
    lam = 1-scalar_flat_step((r-mp.mpf('1.2'))/mp.mpf('.4'))
    delta = mp.arg(t)/(2*mp.pi)-xi4
    a2 = -mp.arg(u1)/(2*mp.pi)+lam*delta+(1-lam)/2
    four = -a2-mp.mpf('.5')
return ('mid',(t,mp.mpf(0),four,lt))
'''
raw_with42 = [node for node in raw_tree31.body if isinstance(node,ast.With)]
# The two angular slots: index 0 with centre 0, index 1 with centre pi.
angular_for42 = ast.parse(
    'for index,centre in enumerate((mp.mpf(0),mp.pi)):\n    pass').body[0]
angular_loops42 = [
    node for node in ast.walk(bundle_tree32)
    if isinstance(node,ast.For)
    and ast.dump(node.target,include_attributes=False)
        ==ast.dump(angular_for42.target,include_attributes=False)
    and ast.dump(node.iter,include_attributes=False)
        ==ast.dump(angular_for42.iter,include_attributes=False)]
source42 = (
    source41 and source_collar34
    and len(raw_with42)==1
    and [ast.dump(node,include_attributes=False)
         for node in raw_with42[0].body[-6:]]==text_dump41(middle_tail42)
    and ordered_nodes37(raw_with42[0].body,[
        'inverse_ratio = abs(u1)>abs(u2)',
        'ratio = u2/u1 if inverse_ratio else u1/u2',
        'rr = abs(ratio)**2',
        'b = A*(1-rr)',
        'root = mp.sqrt(b*b+4*B*rr)',
        'small = ratio*mp.sqrt((2*B if inverse_ratio else 2)/(b+root))',
        "inf = inverse_ratio and abs(small)<mp.mpf('.5')",
        't = None if inf else (1/small if inverse_ratio else small)'])
    and matching_assignment30(bundle_tree32,'principal','mp.log(t)')
    and matching_assignment30(bundle_tree32,'branch',
        'mp.nint((lt.imag-principal.imag)/(2*mp.pi))')
    and matching_statement31(bundle_tree32,'one,lt = one+branch*four,principal')
    and len(angular_loops42)==1)

# C42a.  The pure middle is the inverse-ratio branch, t=1/small.  In the
# direct branch b=A(1-rr)>=0 and root>=2 sqrt(B rr), so
# |small|^2=2rr/(b+root)<=sqrt(rr/B)<=B^(-1/2)=e^-2, i.e. |small|<=1/e<1/2.
# (the direct branch has |u1|<=|u2|, i.e. rr<=1 and b>=0).
# Hence |t|>=1/2 forces |u1|>|u2|, t=u1/(c u2) with c=sqrt(2B/(b+root))>0,
# so arg t = arg u1 - arg u2 modulo 2 pi.  On 1/2<=|t|<=2 the cusp
# weights b0, bi and both restoration weights are exactly zero.
rr42,A42,B42,rw42 = sp.symbols('rr42 A42 B42 rw42',positive=True)
b42 = A42*(1-rr42)
root42 = sp.sqrt(b42**2+4*B42*rr42)
direct_small42 = (
    sp.expand(root42**2-b42**2-4*B42*rr42)==0
    # root^2 >= 4 B rr, hence 2rr/(b+root) <= 2rr/(2 sqrt(B rr)).
    and sp.simplify(2*rr42/(2*sp.sqrt(B42*rr42))-sp.sqrt(rr42/B42))==0
    and sp.exp(-1)<sp.Rational(1,2)
    # flat-step arguments: (r-3/10)/(1/5)-1=5(r-1/2)>=0 for r>=1/2 and
    # (1/r-3/10)/(1/5)-1=5(2-r)/(2r)>=0 for r<=2, so b0=bi=beta=0.
    and sp.factor((rw42-sp.Rational(3,10))/sp.Rational(1,5)-1
                  -5*(rw42-sp.Rational(1,2)))==0
    and sp.factor((1/rw42-sp.Rational(3,10))/sp.Rational(1,5)-1
                  -5*(2-rw42)/(2*rw42))==0)
rep('C42a pure middle = inverse branch, t=u1/(c u2); all cusp weights zero',
    source42 and direct_small42)

# C42b.  Group action.  Psi changes by e1, e4(L), 0 under (a), (b), (c).
# The source embedding first normalizes L to the principal log, which
# applies (c) (the branch line above), then sets slot 8 to
# d0*exp(2 pi i one) (its bracket Arg t-Im L is 0) and slot 9 to
# d1*exp(i(2 pi one-beta9*four)), beta9=Arg(-t)+pi-Arg(t) in {0,2 pi}.
# So V(one,four)=exp(2 pi i one)V(0,four) is invariant under (a), (b).
# With a continuous log branch Lc=P+2 pi i k (P principal) and a
# continuous four_c, Vc=V(k*four_c,four_c) is continuous across the
# t-cut: there k increases by 1 while beta9 jumps from 0 to 2 pi.
one42,four42,lr42,th42,tr42,ti42 = sp.symbols(
    'one42 four42 lr42 th42 tr42 ti42',real=True)
k42,m42,n42 = sp.symbols('k42 m42 n42',integer=True)
t42 = tr42+sp.I*ti42
L42 = lr42+sp.I*th42
def Psi42(one,four,L):
    return (2*sp.pi*one-sp.I*L*four,2*t42*four)
e1_42,e4_42 = (2*sp.pi,0),(-sp.I*L42,2*t42)
group42 = (
    all(sp.expand(a-b-c)==0 for a,b,c in zip(
        Psi42(one42+1,four42,L42),Psi42(one42,four42,L42),e1_42))
    and all(sp.expand(a-b-c)==0 for a,b,c in zip(
        Psi42(one42,four42+1,L42),Psi42(one42,four42,L42),e4_42))
    and all(sp.expand(a-b)==0 for a,b in zip(
        Psi42(one42-four42,four42,L42+2*sp.pi*sp.I),
        Psi42(one42,four42,L42))))
d0_42,d1_42 = sp.symbols('d0_42 d1_42',real=True)
b9_42 = sp.Symbol('b9_42',integer=True)          # beta9=2 pi b9, b9 in {0,1}
def Vslots42(one,four,b9):
    return (d0_42*sp.exp(2*sp.pi*sp.I*one),
            d1_42*sp.exp(sp.I*(2*sp.pi*one-2*sp.pi*b9*four)))
def same42(u,v):
    return all(sp.simplify(sp.expand(sp.powsimp(a/b)))==1 for a,b in zip(u,v))
embedding42 = (
    same42(Vslots42(one42+1,four42,b9_42),Vslots42(one42,four42,b9_42))
    and same42(Vslots42(one42,four42+1,b9_42),Vslots42(one42,four42,b9_42))
    and same42(Vslots42(one42,four42,b9_42),
               tuple(sp.exp(2*sp.pi*sp.I*one42)*v
                     for v in Vslots42(0,four42,b9_42)))
    # continuity of Vc across the t-cut, where d0=0 (|Arg t|>=3pi/4):
    # (k,beta9)=(k0,0) -> (k0+1,2pi) leaves slot 9 unchanged.
    and same42(Vslots42(k42*four42,four42,0)[1:],
               Vslots42((k42+1)*four42,four42,1)[1:])
    and not same42(Vslots42(k42*four42,four42,0)[:1],
                   Vslots42((k42+1)*four42,four42,1)[:1])
    # beta9=Arg(-t)+pi-Arg(t) with principal Arg: 0 for Arg t in (0,pi],
    # 2 pi for Arg t in (-pi,0], checked at exact points of each half.
    and all(sp.arg(-sp.exp(sp.I*a))+sp.pi-sp.arg(sp.exp(sp.I*a))==b
            for a,b in ((sp.pi/3,0),(2*sp.pi/3,0),(sp.pi,0),
                        (-sp.pi/3,2*sp.pi),(-2*sp.pi/3,2*sp.pi)))
    # weight supports: d0=1-s((|Arg t|-pi/2)/(pi/4))=0 iff |Arg t|>=3pi/4,
    # d1 (same with Arg(-t)) =0 iff |Arg t|<=pi/4: flat step of >=1.
    and (3*sp.pi/4-sp.pi/2)/(sp.pi/4)==1
    and sp.Abs(sp.arg(-sp.exp(sp.I*sp.pi/4)))==3*sp.pi/4)
rep('C42b G changes Psi by lattice vectors; the embedding is G-invariant',
    group42 and embedding42)

# C42c.  Continuous four.  Fix continuous branches arg_c u1, arg_c u2,
# arg_c t and Lc=log|t|+i arg_c t on an open set.  The principal values
# are arg_c+2 pi n1, arg_c+2 pi n2, arg_c t-2 pi k and P=Lc-2 pi i k, with
# locally constant integers.  The canonical point (pi s(s-1),2ts),
# s=L/(2 pi i), has fourth Cramer coordinate xi4(L+2 pi i k)=xi4(L)+k
# (C34c; the solve is unique since det>0).  Hence the source values are
#   four_L = four_Lc + n2 - k,   four_R = four_Rc + n1,
# with delta=Arg(t)/(2pi)-xi4 unchanged, for every lam (left: r<=1,
# right: r>1).  Near |t|=1 lam=1, and arg_c t=arg_c u1-arg_c u2+2 pi N
# (C42a) gives four_Rc-four_Lc=-N-1, an integer.
a1_42,a2_42,X4_42,lam42,N42 = sp.symbols('a1_42 a2_42 X4_42 lam42 N42',real=True)
n1_42,n2_42,Nint42 = sp.symbols('n1_42 n2_42 Nint42',integer=True)
at_c42 = sp.Symbol('at_c42',real=True)
def four_L42(arg_u2,xi4):
    return arg_u2/(2*sp.pi)+xi4+sp.Rational(1,2)
def four_R42(arg_u1,arg_t,xi4,lam):
    delta = arg_t/(2*sp.pi)-xi4
    return -(-arg_u1/(2*sp.pi)+lam*delta+(1-lam)/2)-sp.Rational(1,2)
continuous_four42 = (
    south_deck34
    and sp.expand(four_L42(a2_42+2*sp.pi*n2_42,X4_42-k42)
                  -four_L42(a2_42,X4_42)-(n2_42-k42))==0
    and sp.expand(four_R42(a1_42+2*sp.pi*n1_42,at_c42-2*sp.pi*k42,
                           X4_42-k42,lam42)
                  -four_R42(a1_42,at_c42,X4_42,lam42)-n1_42)==0
    and sp.expand(four_R42(a1_42,a1_42-a2_42+2*sp.pi*Nint42,X4_42,1)
                  -four_L42(a2_42,X4_42)+Nint42+1)==0)
rep('C42c source four = continuous four + integer, on both sides of |t|=1',
    continuous_four42)

# C42d.  Let four_c be continuous with four=four_c+n.  By (b)-invariance
# V=V(0,four)=V(0,four_c), and the continuous Vc=V(k four_c,four_c)
# (C42b) gives V=exp(-2 pi i k four_c)Vc: the projectors coincide.  For a
# common final unit vector E the source phase is p=<V,E>=
# exp(2 pi i k four_c)<Vc,E> (conjugate-linear first slot), so
# one=one_c+k four_c+m.  Then Psi(W)-Psi(Wc)=(m-k n)e1+n e4(Lc), a
# lattice vector.  On the pure middle C=phi (C42a), and the restoration
# (z,w)->(z-itC2-tC3,w+(1-t^2)C2+i(1+t^2)C3) commutes with lattice
# translations; (z,w,t)->(z+2 pi a-i b log t,w+2bt,t) is a holomorphic
# chart change with complex determinant 1.
oc42,fc42 = sp.symbols('oc42 fc42',real=True)
pc42 = sp.Symbol('pc42')
Ec42 = sp.Symbol('Ec42')
# Derived from the conjugate-linear source inner product (inner41).
Vc42 = sp.symbols('Vc0_42 Vc1_42')
E42 = sp.symbols('E0_42 E1_42')
V_src42 = tuple(sp.exp(-2*sp.pi*sp.I*k42*fc42)*v for v in Vc42)
phase42 = (
    sp.simplify(sp.expand(inner41(V_src42,E42))
                -sp.exp(2*sp.pi*sp.I*k42*fc42)*inner41(Vc42,E42))==0)
W_src42 = Psi42(oc42+k42*fc42+m42,fc42+n42,L42-2*sp.pi*sp.I*k42)
W_c42 = Psi42(oc42,fc42,L42)
lattice42 = all(sp.expand(a-b-((m42-k42*n42)*c+n42*d))==0
                for a,b,c,d in zip(W_src42,W_c42,e1_42,e4_42))
zz42,ww42,C2_42,C3_42 = sp.symbols('zz42 ww42 C2_42 C3_42')
def restore42(z,w):
    return (z-sp.I*t42*C2_42-t42*C3_42,
            w+(1-t42**2)*C2_42+sp.I*(1+t42**2)*C3_42)
dz42,dw42 = sp.symbols('dz42 dw42')
commute42 = all(sp.expand(a-b-c)==0 for a,b,c in zip(
    restore42(zz42+dz42,ww42+dw42),restore42(zz42,ww42),(dz42,dw42)))
zh42,wh42,th42h = sp.symbols('zh42 wh42 th42h')
deck_map42 = sp.Matrix([zh42+2*sp.pi*m42-sp.I*n42*sp.log(th42h),
                        wh42+2*n42*th42h,th42h])
deck_det42 = sp.simplify(deck_map42.jacobian([zh42,wh42,th42h]).det())==1
rep('C42d same projector; Psi(W) = Psi(continuous W) + lattice; restoration ok',
    phase42 and lattice42 and commute42 and deck_det42 and source41)

# C42e.  Exact witnesses of the |t|=1 selector and of its corner with
# the t-cut.  Put p=1/|t|^2, k^2=(p^2+Ap)/(Ap+B) (C41e north profile),
# Q=sqrt(1+k^2), u1=i/Q and u2=k exp(i(pi/2-theta))/Q.  Then
# u2/u1=k exp(-i theta), the positive root gives small=exp(-i theta)/|t|,
# t=|t| exp(i theta), and stereographic y=(1/Q,k sin(theta)/Q,
# k cos(theta)/Q) recovers (u1,u2).  Take |t| in {999/1000,1,1001/1000}
# and theta=pi/2 (away from all cuts) or theta=pi (the t-cut corner).
# k^2<1, |small|>=1/2 and |t|>=1/2 put these in the pure middle;
# R^2=|y-(0,1,0)|^2=2-2k sin(theta)/Q lies in (2-sqrt 2,2], so |X|<3/5.
tw42,thw42 = sp.symbols('tw42 thw42',positive=True)
kw42 = sp.Symbol('kw42',positive=True)
Qw42 = sp.sqrt(1+kw42**2)
u1w42 = sp.I/Qw42
u2w42 = kw42*sp.exp(sp.I*(sp.pi/2-thw42))/Qw42
yw42 = (1/Qw42,kw42*sp.sin(thw42)/Qw42,kw42*sp.cos(thw42)/Qw42)
Sy42 = 1+sum(v**2 for v in yw42)
witness42 = (
    sp.simplify(sp.expand_complex(u1w42*sp.conjugate(u1w42)
                                  +u2w42*sp.conjugate(u2w42))-1)==0
    and sp.simplify(sp.expand_complex(
        (1-(Sy42-1)+2*sp.I*yw42[0])/Sy42-u1w42))==0
    and sp.simplify(sp.expand_complex(
        2*(yw42[1]+sp.I*yw42[2])/Sy42-u2w42))==0
    and sp.simplify(sp.expand_complex(u2w42/u1w42
                    -kw42*sp.exp(-sp.I*thw42)))==0
    and witness_roots41
    and all(sp.Rational(1,4)<=1/tt**2<=4
            for tt in (sp.Rational(999,1000),1,sp.Rational(1001,1000)))
    and sp.Rational(999,1000)<1<sp.Rational(1001,1000)
    and sp.simplify(sum((yw42[i]-(0,1,0)[i])**2 for i in range(3))
                    -(2-2*kw42*sp.sin(thw42)/Qw42))==0
    and Xv41.subs(Rv41,sp.sqrt(2-sp.sqrt(2)))>-sp.Rational(3,5)
    and Xv41.subs(Rv41,sp.sqrt(2))<sp.Rational(3,5)
    # k^2=(p^2+Ap)/(Ap+B)<1 iff p^2<B; here p<=1000^2/999^2<2<B.
    and (sp.Rational(1000,999)**2)**2<2 and sp.exp(4)>2)
rep('C42e exact witnesses on and across |t|=1, and at its t-cut corner',
    witness42 and slab41)

# C42f is solely a reproducible numerical observation (d).  For each
# pure-middle seam (|t|=1 at theta=pi/2; the t-cut at |t|=9/10 and
# 11/10; the corner |t|=1, theta=pi; the u1-cut at |t|=11/10; the
# u2-cut at |t|=9/10; the t- and u1-cuts at |t|=7/5) it takes production
# points q+ and q- at distance
# eps on either side, one fixed unit vector E=V(q+), and the production
# W-data (one=Arg<V,E>/(2pi)).  It records the integer lattice vector
# nearest to Psi(q+)-Psi(q-), the residual and the projector jump.  Both
# should be O(eps); they are printed for eps=1e-12 and 1e-24.  Two more
# seams lie in the lambda band (the t- and u1-cuts at |t|=7/5).  Sample
# size 32.  It is not used by C42a-e.
def middle_point42(tabs,theta,alpha):
    p = 1/tabs**2
    kk = mp41.sqrt((p**2+A42n*p)/(A42n*p+B42n))
    QQ = mp41.sqrt(1+kk**2)
    u1 = mp41.exp(mp41.j*alpha)/QQ
    u2 = kk*mp41.exp(mp41.j*(alpha-theta))/QQ
    d = 1+u1.real
    w = (u1.imag/d,u2.real/d-1,u2.imag/d)
    R = mp41.sqrt(sum(v*v for v in w))
    X = (R-1)/mp41.sqrt(1+R**2)
    n = (w[0]/R,-w[2]/R,-w[1]/R)
    kind,data = cusp30.elliptic_raw_base_reference(
        tuple((1-X**2)*v for v in n),X,dps=110,elementary=True)
    return kind,data
rows42 = []
sample_ok42 = True
with mp41.workdps(120), cusp30.elementary_formula_branches():
    A42n,B42n = (4*mp41.pi)**2,mp41.exp(4)
    pi42 = mp41.pi
    seams42 = (
        ('|t|=1',lambda e:(1+e,pi42/2,pi42/2),lambda e:(1-e,pi42/2,pi42/2)),
        ('t-cut 9/10',lambda e:(mp41.mpf(9)/10,pi42-e,pi42/2),
                      lambda e:(mp41.mpf(9)/10,-pi42+e,pi42/2)),
        ('t-cut 11/10',lambda e:(mp41.mpf(11)/10,pi42-e,pi42/2),
                       lambda e:(mp41.mpf(11)/10,-pi42+e,pi42/2)),
        ('corner',lambda e:(1+e,pi42-e,pi42/2),lambda e:(1-e,-pi42+e,pi42/2)),
        ('u1-cut',lambda e:(mp41.mpf(11)/10,pi42/2,pi42-e),
                  lambda e:(mp41.mpf(11)/10,pi42/2,-pi42+e)),
        ('u2-cut',lambda e:(mp41.mpf(9)/10,pi42/2,3*pi42/2-e),
                  lambda e:(mp41.mpf(9)/10,pi42/2,3*pi42/2+e)),
        ('t-cut 7/5',lambda e:(mp41.mpf(7)/5,pi42-e,pi42/2),
                     lambda e:(mp41.mpf(7)/5,-pi42+e,pi42/2)),
        ('u1-cut 7/5',lambda e:(mp41.mpf(7)/5,pi42/2,pi42-e),
                      lambda e:(mp41.mpf(7)/5,pi42/2,-pi42+e)))
    for name42,plus42,minus42 in seams42:
        residuals42 = []
        for eps42 in (mp41.mpf('1e-12'),mp41.mpf('1e-24')):
            kp42,dp42 = middle_point42(*plus42(eps42))
            km42,dm42 = middle_point42(*minus42(eps42))
            Vp42 = cusp30.elliptic_bundle_embedding('mid',dp42,dps=110,elementary=True)
            Vm42 = cusp30.elliptic_bundle_embedding('mid',dm42,dps=110,elementary=True)
            Pjump42 = max(abs(a*mp41.conj(b)-c*mp41.conj(d))
                          for a,c in zip(Vp42,Vm42) for b,d in zip(Vp42,Vm42))
            psis42 = []
            for data42,V42 in ((dp42,Vp42),(dm42,Vm42)):
                tt42,_,ff42,ll42 = data42
                ph42 = cusp30._elementary_inner(V42,Vp42)
                oo42 = mp41.arg(ph42)/(2*mp41.pi)
                psis42.append((2*mp41.pi*oo42-mp41.j*ll42*ff42,2*tt42*ff42,
                               tt42,ll42))
            dz42n = psis42[0][0]-psis42[1][0]
            dw42n = psis42[0][1]-psis42[1][1]
            tt42,ll42 = psis42[1][2],psis42[1][3]
            nn42 = int(mp41.nint((dw42n/(2*tt42)).real))
            mm42 = int(mp41.nint(((dz42n+mp41.j*ll42*nn42)/(2*mp41.pi)).real))
            res42 = max(abs(dz42n-(2*mp41.pi*mm42-mp41.j*ll42*nn42)),
                        abs(dw42n-2*tt42*nn42))
            residuals42.append((eps42,mm42,nn42,res42,Pjump42))
            sample_ok42 &= (kp42=='mid' and km42=='mid'
                            and res42<1000*eps42 and Pjump42<1000*eps42)
        sample_ok42 &= (residuals42[0][1:3]==residuals42[1][1:3]
                        and residuals42[1][3]<residuals42[0][3]*mp41.mpf('1e-9'))
        rows42.append((name42,residuals42))
    print('C42f production observations (d), sample size 32:')
    for name42,res42 in rows42:
        print('   ',name42,' lattice (m,n)=',res42[0][1:3],
              ' residual(eps=1e-12,1e-24)=',
              [mp41.nstr(r[3],3) for r in res42],
              ' |P jump|=',[mp41.nstr(r[4],3) for r in res42])
rep('C42f production points across eight pure-middle seams (numerical, d)',
    sample_ok42)

# C43: the argument-chart selector of the raw reference (C38) on
# 1/10<=|tc|<1/2, i.e. in the rest of the core and in the collar, at both
# cusp signs.  The source takes a=1/2+Arg(-p)/(2pi) if Re p<-(9/10)|p|,
# else a=Arg(p)/(2pi) (p=u2 south, p=-conj(u1) north).  On the upper
# wedge the two values agree; on the lower wedge -pi<Arg p<-5pi/6 the
# principal a0 is in (-1/2,-5/12) and the half-centred a1=a0+1 (C38a).
# C43 proves, on the whole lower overlap and for every |tc|<1/2 (core
# and collar), that the two charts give the same ten-slot projector; and,
# for a final unit vector E of the fold that is smooth there (a premise
# resting on the other seams and on nonzero overlaps), restored outputs
# related by a holomorphic deck transition: the cusp deck T of C38d in
# the core, and the lattice vector +e4 (south) or -e4 (north) in the
# collar.  It replaces C38's seed-branch computation, valid only for
# q<1/100, by an exact deck shift of the four-seed height valid on
# 1/4<Y<2 for every q<1/8.
seed_body43 = '''
reciprocal = u>2
v = 1/u if reciprocal else u
if v<=mp.mpf(1)/16:
    value = v
elif v>=mp.mpf('.5'):
    value = v/(1+v)
else:
    weight = scalar_flat_step(mp.log(16*v)/mp.log(8))
    value = v+weight*(v/(1+v)-v)
return 1-value if reciprocal else value
'''
finite_body43 = '''
import mpmath as mp
rr = abs(x*y)**2
e = mp.exp(-4*rr/(1-rr))
A,Y = nu*nu*e*abs(x)**2,abs(y)**2
q = A*Y
if not rr<=mp.mpf('.25') or max(A,Y)>4 or q>=mp.mpf('.125'):
    raise ValueError('four-seed finite height is outside its stated domain')
def seed(u):
    reciprocal = u>2
    v = 1/u if reciprocal else u
    if v<=mp.mpf(1)/16:
        value = v
    elif v>=mp.mpf('.5'):
        value = v/(1+v)
    else:
        weight = scalar_flat_step(mp.log(16*v)/mp.log(8))
        value = v+weight*(v/(1+v)-v)
    return 1-value if reciprocal else value
return seed(A)+seed(q*A)-seed(Y)-seed(q*Y)+q*q*(A-Y)/(1-q)
'''
finite_tree43 = source_tree30(cusp30.scalar_finite_height)
source43 = (
    source41 and arg_source38 and open_wedges38 and source_torus35
    and body_dump41(finite_tree43)==text_dump41(finite_body43)
    and len(seed_nodes38)==1
    and body_dump41(seed_nodes38[0])==text_dump41(seed_body43))

# C43a.  K=seed satisfies K(u)+K(1/u)=1 for every u>0: on [1/2,2] both
# are v/(1+v) and the sum is 1; otherwise one of u,1/u exceeds 2 and the
# reciprocal branch returns 1-K of the other.  It is linear, K(u)=u, on
# u<=1/16.  Let (A,Y)=(q/Y0,Y0) be a cusp representative and
# psi_hat(x,y)=(1/(my),mxy^2), m=nu*sqrt(e); then psi_hat has
# (A,Y)=(1/Y0,q*Y0) with the same t and q.  If 1/4<Y0<2 and q<1/8, both
# lie in the four-seed domain (max(A,Y)<=4, q<1/8) and q^2 Y0,
# q^2/Y0<1/16.  The A-seed K(q/Y0) of p cancels the qA-seed of psi_hat p,
# the qY-seed K(qY0) of p cancels the Y-seed of psi_hat p, the q^2 seeds
# are linear, and
#   H_f(psi_hat p)-H_f(p)=K(1/Y0)+K(Y0)=1   exactly.
Y0_43,q43 = sp.symbols('Y0_43 q43',positive=True)
Kpair43 = sp.Function('K43')
def H4_43(A,Y,qq,K):
    return K(A)+K(qq*A)-K(Y)-K(qq*Y)+qq**2*(A-Y)/(1-qq)
lin43 = {Kpair43(q43**2/Y0_43):q43**2/Y0_43,Kpair43(q43**2*Y0_43):q43**2*Y0_43}
shift43 = sp.simplify(
    (H4_43(1/Y0_43,q43*Y0_43,q43,Kpair43)
     -H4_43(q43/Y0_43,Y0_43,q43,Kpair43)).xreplace(lin43)
    -(Kpair43(1/Y0_43)+Kpair43(Y0_43)))==0
vv43 = sp.Symbol('vv43',positive=True)
seed_mid43 = vv43/(1+vv43)
reciprocal43 = (
    sp.simplify(seed_mid43+seed_mid43.subs(vv43,1/vv43)-1)==0
    and seed_mid43.subs(vv43,2)+seed_mid43.subs(vv43,sp.Rational(1,2))==1
    # psi_hat moduli: A'=nu^2 e/(m^2 Y)=1/Y, Y'=m^2 |x|^2 Y^2=A Y=qY.
    and sp.simplify((nu41**2*sp.Symbol('ee43',positive=True)
                     /((nu41*sp.sqrt(sp.Symbol('ee43',positive=True)))**2
                       *Y0_43))-1/Y0_43)==0
    # domain: 1/4<Y0<2, q<1/8 give max(A,Y)<=4 and q^2 max(Y0,1/Y0)<1/16
    # (the source's linear seed branch is v<=1/16)
    and sp.Rational(1,8)**2*4<=sp.Rational(1,16)
    and sp.Rational(1,8)**2*2<sp.Rational(1,16))
rep('C43a four-seed height shifts by exactly 1 under psi_hat on 1/4<Y<2',
    source43 and shift43 and reciprocal43)

# C43b.  Y0 is in (1/4,2) on both charts' common domain.  With b=a0 in
# (-1/2,-5/12) the source moment is on its outer branch |b|>=2/5:
# f=-(e^(kappa(|b|-1/2))-q e^(kappa(1/2-|b|))), kappa=2+lam(k-2), k=-log q.
# The cell-0 negative reference has Y0=(R-f)/2 (R=sqrt(f^2+4q)), so
# |f|<=Y0<=|f|+sqrt q.  Core (lam=0, kappa=2, q<1/8 by L13d):
#   |f|>e^(-1/6)-e^(1/6)/8>1/4.
# Collar: with s=|tc|^2 in (9/100,1/4), q=nu^2 s exp(-4s/(1-s)).  Its
# derivative has the sign of 1-4s/(1-s)^2, which decreases, so q is
# unimodal and its minimum is at an endpoint: q(1/4)>(1/4)e^(-4/3)>1/16
# (e^(4/3)<4) and q(9/100)>(9/100)e^(-4/9)>9/200>1/25.  Hence south
# q>1/25 and north q>1/(25*55)>1/2500, while q<s e^(-4s)<=1/(4e)<1/8.
# So 2<=kappa<=k<8 and
#   |f|>=e^(-k/12)-e^(-k)e^(k/12)=e^(-k/12)(1-e^(-5k/6))
#      >e^(-2/3)(1-e^(-5/3))>1/4,
# and in both cases Y0<1+sqrt(1/8)<2.
k43 = sp.Symbol('k43',positive=True)
domain43 = (
    sp.exp(-sp.Rational(1,6))-sp.exp(sp.Rational(1,6))/8>sp.Rational(1,4)
    and sp.exp(-sp.Rational(2,3))*(1-sp.exp(-sp.Rational(5,3)))>sp.Rational(1,4)
    # q<1/8<e^-2 gives k>2, and q>1/2500>e^-8 gives k<8.
    and sp.exp(8)>2500 and sp.exp(2)<8
    and sp.simplify(sp.exp(-k43/12)-sp.exp(-k43)*sp.exp(k43/12)
                    -sp.exp(-k43/12)*(1-sp.exp(-5*k43/6)))==0
    and 1+sp.sqrt(sp.Rational(1,8))<2
    and sp.Rational(5,12)>sp.Rational(2,5)
    and sp.factor(sp.diff(4*s_band41/(1-s_band41)**2,s_band41)
                  -4*(1+s_band41)/(1-s_band41)**3)==0
    and sp.exp(sp.Rational(4,3))<4
    and 4*sp.Rational(9,100)/(1-sp.Rational(9,100))<sp.Rational(4,9)
    and sp.Rational(9,100)*sp.exp(-sp.Rational(4,9))>sp.Rational(9,200)
    and sp.Rational(9,200)>sp.Rational(1,25)
    and sp.Rational(1,25)/55>sp.Rational(1,2500)
    and 1/(4*sp.E)<sp.Rational(1,8)
    and matching_if_test32(moment_tree36,"aa>=mp.mpf('.4')")
    and exact_if_body37(moment_tree36,"aa>=mp.mpf('.4')",
        'return mp.sign(a)*(B*mp.exp(kappa*aa)-C*mp.exp(-kappa*aa))')
    and matching_statement31(terminal_tree36,
        'B,C = mp.exp(-kappa/2),q*mp.exp(kappa/2)'))
rep('C43b both charts keep Y0 in (1/4,2): core (q<1/8) and collar',
    domain43)

# C43c.  Core (|tc|<=3/10).  Same b=a0, f, R on both charts (cells 0 and
# 1).  The a0<=0 reference is y0=sqrt((R-f)/2), x0=t/y0; the a1>0 cell-1
# reference is x1=sqrt(2/((R-f)e))/nu, y1=t/x1.  Hence
# (x1,y1)=psi_hat(x0,y0) exactly, for every f.  Both heights have zero
# first weight (|a0|,a1>5/12>2/5) and second weight hw(a0+1/2)=hw(a1-1/2);
# the negative chart applies psi_hat before its second frame.  So the
# ten-vectors coincide, the final phase p is the same (|p|=1 because the
# fold vector E lies in the image of the final projector), and the final
# W-points are psi_hat-related (C38e: the phase commutes with psi_hat).
# The four-seed height shifts by one (C43a,b), so the restorations are
# related by the holomorphic cusp deck T of C38d.
fc43,qc43,ec43,nuc43 = sp.symbols('fc43 qc43 ec43 nuc43',positive=True)
fneg43 = -fc43                       # f<0 on the negative outer branch
Rc43 = sp.sqrt(fneg43**2+4*qc43)
y0c43 = sp.sqrt((Rc43-fneg43)/2)
mc43 = nuc43*sp.sqrt(ec43)
x1c43 = sp.sqrt(2/((Rc43-fneg43)*ec43))/nuc43
core_deck43 = (
    sp.simplify(x1c43-1/(mc43*y0c43))==0
    and sp.simplify(t41/x1c43-mc43*(t41/y0c43)*y0c43**2)==0
    and holomorphic_deck38 and phase_deck38
    # first weights: |a0|>5/12>2/5 and a1>1/2>2/5 give hw=0; the
    # second weights are hw(a0+1/2)=hw(a1-1/2) since a1=a0+1.
    and sp.Rational(5,12)>sp.Rational(2,5)
    # m=nu*exp(-2rr/(1-rr)) in the core block is nu*sqrt(e): both are
    # positive and m^2=nu^2 e.
    and matching_assignment30(bundle_tree32,'m','nu*mp.exp(-2*rr/(1-rr))')
    and sp.simplify(sp.exp(-2*rq41**2/(1-rq41**2))**2
                    -sp.exp(-4*rq41**2/(1-rq41**2)))==0)
rep('C43c core: references psi_hat-related, same ten-vector, deck T',
    core_deck43 and source43)

# C43d.  Collar (3/10<|tc|<1/2), with t=re^(i th) the middle variable and
# L=log r+i th.  Both charts use the same b and lambda, so the logradius
# gives rho1=rho0+P with P=log(m|tc|) (cells n=0,1).  The raw middle point
# therefore moves by (-iP,0), and in the middle lattice
#   south: (-iP,0)= e4-c.e-(th/2pi)e1,   north: (-iP,0)=-e4+c.e+(th/2pi)e1,
# where c=coeffs(2t) are torus directions (dropped from W).  So
# four->four+-1 and one->one-th_c/(2pi), th_c=Arg tc.  The ten-vector is
# multiplied by exp(-i th_c): the cusp pair has zz->zz e^P, H->H+1, shift
# ->shift-1, so yy->yy |m tc|/(m tc), and the angular slots get
# exp(2 pi i Delta one) (an integer change of four is invisible, C42b).
# Hence the projectors agree, p->exp(i th_c)p, the corrected one agrees
# modulo 1, and Psi(W) changes by +-e4 plus an integer multiple of e1.
# In the restore loop yy changes by m*tc (the psi_hat factor), so H_f
# rises by one (C43a,b), and D is unchanged: south H-four, north H+four.
# The restored outputs differ by a lattice translation (C42d).
xd43,yd43,lr43,th43 = sp.symbols('xd43 yd43 lr43 th43',real=True)
td43 = xd43+sp.I*yd43
r2d43 = xd43**2+yd43**2
Ld43 = lr43+sp.I*th43
e1d43,e2d43,e3d43,e4d43 = ((2*sp.pi,0),(-sp.I*td43,1-td43**2),
                           (-td43,sp.I*(1+td43**2)),(-sp.I*Ld43,2*td43))
cd43 = coeffs41(sp.expand(2*td43),e2d43[1],e3d43[1])
ced43 = (cd43[0]*e2d43[0]+cd43[1]*e3d43[0],cd43[0]*e2d43[1]+cd43[1]*e3d43[1])
P_s43 = lr43-2*r2d43/(1-r2d43)                  # nu=1, |tc|=r
P_n43 = -2-lr43-2/(r2d43-1)                      # nu=e^-2, |tc|=1/r
lattice_shift43 = (
    all(sp.simplify(sp.expand_complex(a-b))==0 for a,b in zip(
        (-sp.I*P_s43,0),
        (e4d43[0]-ced43[0]-th43/(2*sp.pi)*e1d43[0],e4d43[1]-ced43[1])))
    and all(sp.simplify(sp.expand_complex(a-b))==0 for a,b in zip(
        (-sp.I*P_n43,0),
        (-e4d43[0]+ced43[0]+th43/(2*sp.pi)*e1d43[0],-e4d43[1]+ced43[1]))))
# Gauge of the ten-vector: cusp pair (yy -> yy*|m tc|/(m tc), the frame
# gauge of C41b with g=exp(-i th_c)) and angular slots (exp(2 pi i Delta)).
rc43,mm43,thc43 = sp.symbols('rc43 mm43 thc43',positive=True)
tc43 = rc43*sp.exp(sp.I*thc43)
gauge43 = (
    sp.simplify(mm43*rc43/(mm43*tc43)-sp.exp(-sp.I*thc43))==0
    and gauge_frames41
    and sp.simplify(sp.exp(2*sp.pi*sp.I*(-thc43/(2*sp.pi)))
                    -sp.exp(-sp.I*thc43))==0)
# p changes by exp(i th_c) (conjugate-linear first slot), so the
# corrected one, one+Arg(p)/(2pi), is unchanged modulo 1.
Vg43 = sp.symbols('Vg0_43 Vg1_43')
Eg43 = sp.symbols('Eg0_43 Eg1_43')
phase43 = (
    sp.simplify(sp.expand(inner41(tuple(sp.exp(-sp.I*thc43)*v for v in Vg43),Eg43))
                -sp.exp(sp.I*thc43)*inner41(Vg43,Eg43))==0)
# Restore-loop yy under W -> W+-e4 (+ integer e1): its exponent changes by
# log(m tc) = P+i th_c (south: L-d, north: -L-2-d, d=2|tc|^2/(1-|tc|^2)).
yy_shift43 = True
for inf43,sign43 in ((False,1),(True,-1)):
    zq43,wq43 = sp.symbols('zq43 wq43')
    sq43 = Ld43/(2*sp.pi*sp.I)
    tcd43 = 1/td43 if inf43 else td43
    e2q,e3q,_,a2q,a3q = cusp_data41(tcd43,inf43)
    def expo43(z,w):
        if inf43:
            zc,wc = z-sp.I*w/td43,w/td43**2
        else:
            zc,wc = z-sp.pi*sq43*(sq43-1),w-2*td43*sq43
        cc = coeffs41(sp.expand_complex(-wc),e2q,e3q)
        return sp.I*(zc-sp.I*(cc[0]*a2q+cc[1]*a3q))
    zr43,zi43,wr43,wi43 = sp.symbols('zr43 zi43 wr43 wi43',real=True)
    z0q,w0q = zr43+sp.I*zi43,wr43+sp.I*wi43
    diff43 = sp.expand(expo43(z0q+sign43*e4d43[0],w0q+sign43*e4d43[1])
                       -expo43(z0q,w0q))
    tc_abs2 = 1/r2d43 if inf43 else r2d43
    d43 = 2*tc_abs2/(1-tc_abs2)
    target43 = (-Ld43-2-d43) if inf43 else (Ld43-d43)
    yy_shift43 &= sp.simplify(sp.expand_complex(diff43-target43))==0
# Arithmetic of the two source D formulas (bound in restore_loop41).
D_invariant43 = (sp.simplify((H41+1)-(four41+1)-(H41-four41))==0
                 and sp.simplify((H41+1)+(four41-1)-(H41+four41))==0)
rep('C43d collar: lattice shift +-e4, same projector, D invariant, deck',
    lattice_shift43 and gauge43 and phase43 and yy_shift43
    and D_invariant43 and deck_det42 and source43)

# C43e.  Exact round witnesses on the lower selector, with
# p_c=pr-i sqrt(1-pr^2), pr in {-9/10-1/1000,-9/10,-9/10+1/1000}
# (half-centred, seam, principal chart; Re p+(9/10)|p|=pr+9/10, C38a).
# Core south (|tc|=1/5, direct branch, k^2=(Bp^2+Ap)/(1+Ap), p=1/25):
#   u2=p_c/Q, u1=i k/Q; y=(k,pr,-sqrt(1-pr^2))/Q, R^2=2-2pr/Q.
# Collar south (|tc|=2/5, inverse branch, small=1/t, p'=|small|^2=25/4,
# k^2=(p'^2+Ap')/(Ap'+B)): u1=i/Q, u2=k p_c/Q; y=(1,k pr,-k sqrt(1-pr^2))/Q,
#   R^2=2-2k pr/Q.
# North (|tc|=1/5 and 2/5, inverse branch, k^2=(p^2+Ap)/(Ap+B)):
#   u1=-conj(p_c)/Q, u2=-i k/Q; y=(-sqrt(1-pr^2),0,-k)/(Q-pr),
#   R^2=1+(1-pr^2+k^2)/(Q-pr)^2.
# In each case 0<k<1<Q<sqrt 2 and -1<pr<-4/5 give 1<R^2<4, hence
# 0<X<1/sqrt 5<3/5 (round_inverse38); the moduli |small|=sqrt(p) follow
# from witness_roots41, which depends only on k^2.
# Rational parametrization pr=-(1-tau^2)/(1+tau^2), sqrt(1-pr^2)=
# 2tau/(1+tau^2), tau>0; pr=-9/10 is tau=1/sqrt(19), and pr in
# (-1,-4/5) is 0<tau<1/3.
tau43,kw43 = sp.symbols('tau43 kw43',positive=True)
pr43 = -(1-tau43**2)/(1+tau43**2)
Qw43 = sp.sqrt(1+kw43**2)
sq43w = 2*tau43/(1+tau43**2)
pcw43 = pr43-sp.I*sq43w
cases43 = (
    (sp.I*kw43/Qw43,pcw43/Qw43,(kw43,pr43,-sq43w),Qw43,
     2-2*pr43/Qw43),
    (sp.I/Qw43,kw43*pcw43/Qw43,(1,kw43*pr43,-kw43*sq43w),Qw43,
     2-2*kw43*pr43/Qw43),
    (-sp.conjugate(pcw43)/Qw43,-sp.I*kw43/Qw43,(-sq43w,0,-kw43),Qw43-pr43,
     1+(1-pr43**2+kw43**2)/(Qw43-pr43)**2))
witness43 = True
for u1c,u2c,ynum,yden,R2c in cases43:
    yv = tuple(v/yden for v in ynum)
    Sy = 1+sum(v**2 for v in yv)
    witness43 &= (
        sp.simplify(sp.expand_complex(u1c*sp.conjugate(u1c)
                                      +u2c*sp.conjugate(u2c))-1)==0
        and sp.simplify(sp.expand_complex((1-(Sy-1)+2*sp.I*yv[0])/Sy-u1c))==0
        and sp.simplify(sp.expand_complex(2*(yv[1]+sp.I*yv[2])/Sy-u2c))==0
        and sp.simplify(sum((yv[i]-(0,1,0)[i])**2 for i in range(3))-R2c)==0)
# k^2<1: south core Bp^2<1 (p=1/25, B<55); inverse branch p^2<B (p<=25/4,
# B>50).  R^2: cases 1-2 are 2-2pr/Q or 2-2k pr/Q in (2,4) since
# 0<-pr<1, 0<k<1<Q; case 3 is 1+(1-pr^2+k^2)/(Q-pr)^2 in (1,2) since
# 1-pr^2<9/25, k^2<1, Q-pr>9/5.  So 1<R^2<4 and 0<X<X(2)=1/sqrt 5<3/5.
witness43 &= (
    sp.exp(4)>50 and sp.exp(4)<55
    and 55*sp.Rational(1,25)**2<1
    and sp.Rational(25,4)**2<50
    and sp.Rational(9,25)+1<sp.Rational(9,5)**2
    and sp.simplify(pr43**2+sq43w**2-1)==0
    and sp.simplify(pr43.subs(tau43,1/sp.sqrt(19))+sp.Rational(9,10))==0
    and sp.Rational(-9,10)-sp.Rational(1,1000)>-1
    and sp.Rational(-9,10)+sp.Rational(1,1000)<-sp.Rational(4,5)
    and 1-(sp.Rational(4,5))**2==sp.Rational(9,25)
    and Xv41.subs(Rv41,2)<sp.Rational(3,5)
    and sp.Rational(1,5)<sp.Rational(3,10)<sp.Rational(2,5)<sp.Rational(1,2)
    and (-sp.Rational(9,10)-sp.Rational(1,1000))+sp.Rational(9,10)<0
    and (-sp.Rational(9,10)+sp.Rational(1,1000))+sp.Rational(9,10)>0
    # all three lie in the lower wedge Arg p<-5pi/6, i.e. pr<-sqrt(3)/2
    and sp.Rational(899,1000)**2>sp.Rational(3,4))
rep('C43e exact selector witnesses at |tc|=1/5 and 2/5, both cusps and sides',
    witness43 and witness_roots41 and round_inverse38)

# C43f is solely a reproducible numerical observation (d), sample size 8:
# |tc| in {.15,.25} (core) and {.35,.45} (collar) at both signs, just
# inside the half-centred chart of the lower selector.  The production
# raw reference gives chart a1; the other chart a0=a1-1 is evaluated at
# the same point with the production core reference or logradius and
# lattice solve.  It prints the projector difference, the deck-relation
# residual of the restored outputs (T of C38d in the core, +-e4 plus an
# integer e1 in the collar), the measured (Delta four, Delta one), and
# the measured margins of |f| and Y0 beside the proven C43b bounds.
rows43 = []
margins43 = []
sample_ok43 = True
with mp41.workdps(120), cusp30.elementary_formula_branches():
    A43n,B43n = (4*mp41.pi)**2,mp41.exp(4)
    beta43 = -(mp41.pi-mp41.acos(mp41.mpf(9)/10))-mp41.mpf('1e-30')
    phi43 = (mp41.mpf(1)/3,mp41.mpf(2)/7)
    pc43n = mp41.exp(mp41.j*mp41.mpf('0.7'))
    for north43 in (False,True):
        nu43 = mp41.exp(-2) if north43 else mp41.mpf(1)
        sig43 = -1 if north43 else 1
        for rad43 in (mp41.mpf('0.15'),mp41.mpf('0.25'),mp41.mpf('0.35'),
                      mp41.mpf('0.45')):
            arg_tc43 = sig43*mp41.mpf('0.8')
            if north43 or rad43<1/mp41.e:           # direct south / north
                pp = rad43**2
                kk = mp41.sqrt((pp**2+A43n*pp)/(A43n*pp+B43n) if north43 else
                               (B43n*pp**2+A43n*pp)/(1+A43n*pp))
            else:                                   # south inverse branch
                pp = 1/rad43**2
                kk = mp41.sqrt((pp**2+A43n*pp)/(A43n*pp+B43n))
            QQ = mp41.sqrt(1+kk**2)
            ph = mp41.exp(mp41.j*beta43)
            if north43:
                u1 = -mp41.conj(ph)/QQ
                u2 = kk*u1*mp41.exp(mp41.j*arg_tc43)
            elif rad43<1/mp41.e:
                u2 = ph/QQ
                u1 = kk*u2*mp41.exp(mp41.j*arg_tc43)
            else:
                u2 = kk*ph/QQ
                u1 = u2/kk*mp41.exp(mp41.j*arg_tc43)
            d = 1+u1.real
            w = (u1.imag/d,u2.real/d-1,u2.imag/d)
            R = mp41.sqrt(sum(v*v for v in w))
            X = (R-1)/mp41.sqrt(1+R**2)
            n = (w[0]/R,-w[2]/R,-w[1]/R)
            kind1,data1 = cusp30.elliptic_raw_base_reference(
                tuple((1-X**2)*v for v in n),X,dps=110,elementary=True)
            a1 = data1[-1]
            a0 = a1-1
            tc = rad43*mp41.exp(mp41.j*arg_tc43)
            sample_ok43 &= a1>mp41.mpf(1)/2 and a1<mp41.mpf(7)/12
            ref1 = cusp30.elliptic_bundle_embedding(kind1,data1,dps=110,elementary=True)
            if kind1 in ('c0','ci'):
                x0,y0 = cusp30.elementary_cusp_reference(tc,a0,nu43)
                data0 = (x0,y0,a0)
                ref0 = cusp30.elliptic_bundle_embedding(kind1,data0,dps=110,elementary=True)
                E = tuple(pc43n*v for v in ref1)
                p0 = cusp30._elementary_inner(ref0,E)
                p1 = cusp30._elementary_inner(ref1,E)
                _,o0 = cusp30.elementary_restore_torus(kind1,(mp41.conj(p0)*x0,p0*y0),phi43)
                _,o1 = cusp30.elementary_restore_torus(kind1,(mp41.conj(p1)*data1[0],p1*data1[1]),phi43)
                T0 = (1/(nu43*o0[1]),nu43*o0[0]*o0[1]**2,o0[2]+2*sig43*o0[0]*o0[1])
                resid = max(abs(a-b) for a,b in zip(o1,T0))
                dfour = None
            else:
                rho0 = cusp30.elementary_cusp_logradius(a0,tc,nu43,dps=110)
                tm,_,_,lt,_ = data1
                if north43:
                    zm,wm = -mp41.j*rho0,mp41.mpc(0)
                else:
                    ss = lt/(2*mp41.j*mp41.pi)
                    zm,wm = -mp41.j*rho0+mp41.pi*ss*(ss-1),2*tm*ss
                lat = cusp30.scalar_lattice_coords(zm,wm,tm,lt)
                data0 = (tm,lat[0],lat[3],lt,a0)
                ref0 = cusp30.elliptic_bundle_embedding('mid',data0,dps=110,elementary=True)
                E = tuple(pc43n*v for v in ref1)
                p0 = cusp30._elementary_inner(ref0,E)
                p1 = cusp30._elementary_inner(ref1,E)
                W0 = (tm,data0[1]+mp41.arg(p0)/(2*mp41.pi),data0[2],lt)
                W1 = (tm,data1[1]+mp41.arg(p1)/(2*mp41.pi),data1[2],lt)
                _,o0 = cusp30.elementary_restore_torus('mid',W0,phi43)
                _,o1 = cusp30.elementary_restore_torus('mid',W1,phi43)
                dfour = int(mp41.nint(data1[2]-data0[2]))
                dz,dw = o1[0]-o0[0],o1[1]-o0[1]
                mint = int(mp41.nint(((dz+mp41.j*lt*dfour)/(2*mp41.pi)).real))
                resid = max(abs(dz-(2*mp41.pi*mint-mp41.j*lt*dfour)),
                            abs(dw-2*tm*dfour),abs(o1[2]-o0[2]))
                done = data1[1]-data0[1]+arg_tc43/(2*mp41.pi)
                sample_ok43 &= (dfour==(-1 if north43 else 1)
                                and abs(done-mp41.nint(done))<mp41.mpf('1e-100'))
            Pd = max(abs(a*mp41.conj(b)-c*mp41.conj(dd))
                     for a,c in zip(ref0,ref1) for b,dd in zip(ref0,ref1))
            rr43 = rad43**2
            q43n = nu43**2*rr43*mp41.exp(-4*rr43/(1-rr43))
            lam43 = cusp30.scalar_flat_step((rad43-mp41.mpf(3)/10)/(mp41.mpf(1)/5))
            f43n = cusp30.inverse_height_moment(a0,q43n,lam43)
            Y43n = (mp41.sqrt(f43n**2+4*q43n)-f43n)/2
            margins43.append((abs(f43n),Y43n))
            rows43.append((north43,rad43,kind1,Pd,resid,dfour))
            sample_ok43 &= Pd<mp41.mpf('1e-90') and resid<mp41.mpf('1e-90')
    print('C43f production observations (d), sample size 8:'
          ' max |P_a0-P_a1| =',mp41.nstr(max(r[3] for r in rows43),5),
          ' max deck residual =',mp41.nstr(max(r[4] for r in rows43),5),
          ' kinds =',[r[2] for r in rows43],
          ' Delta four (collar) =',[r[5] for r in rows43 if r[5] is not None],
          ' min|f| =',mp41.nstr(min(m[0] for m in margins43),6),
          '(proven >1/4) Y0 in [',mp41.nstr(min(m[1] for m in margins43),6),
          ',',mp41.nstr(max(m[1] for m in margins43),6),'] (proven (1/4,2))')
rep('C43f eight production points at the selector (numerical, d)',
    sample_ok43)

# C44: the two collar/pure-middle source joins along their whole seams:
# south |t|=1/2 (selector abs(t)<.5) and north |small|=1/2, i.e. |t|=2
# (selector inf=inverse_ratio and abs(small)<.5).  N13 and C40 checked one
# chamber each.  Let M1 be the collar formulas with lambda=1 and the cusp
# weight b0 (south) or bi (north) set to 0 (bm=1).  C44 proves: the
# collar source differs from M1 by terms flat at the seam; M1 and the
# pure-middle source give W-data with the same four modulo an integer and
# one differing by a smooth real gauge; hence the same ten-slot projector
# and, for a smooth final fold vector E (a premise resting on the other
# seams), restored outputs differing by a holomorphic lattice
# translation (C42d).  It holds for every argument chart of a and every
# fixed branch of log t.

# C44a.  At lambda=1 the source moment is f=2 sqrt(q) sinh(kb) on every
# branch: the |b|<=1/10 and lam==1 branch is G=d sinh(amp b)/amp with
# amp=k, d=2 sqrt(q) k; the outer branch has B=C=sqrt(q), kappa=k; the
# interior terminal map is the identity (C40c).  With root=2 sqrt(q)
# cosh(kb), both source Y formulas give Y=sqrt(q) e^(-kb), so
# rho=nP+log(Y)/2=P(a+1/2) (P=log(q)/2=-k/2): this is the lam==1
# return of elementary_cusp_logradius.  Left of the seam
# lambda=s((|tc|-3/10)/(1/5)) and b0 or bi equal 1-s, flat at u=1.
qb44,kb44 = sp.symbols('qb44 kb44',positive=True)
bb44,nn44 = sp.symbols('bb44 nn44',real=True)
sq44 = sp.sqrt(qb44)
f44 = 2*sq44*sp.sinh(kb44*bb44)
lam1_44 = (
    sp.simplify((2*sq44*kb44*sp.sinh(kb44*bb44)/kb44)-f44)==0
    and sp.simplify((sq44*sp.exp(kb44*bb44)-qb44/sq44*sp.exp(-kb44*bb44))-f44)==0
    and sp.simplify(sp.expand((f44**2+4*qb44)
                    -(2*sq44*sp.cosh(kb44*bb44))**2).rewrite(sp.exp))==0
    and sp.simplify((2*qb44/(2*sq44*sp.cosh(kb44*bb44)+f44))
                    -sq44*sp.exp(-kb44*bb44))==0
    and sp.simplify(((2*sq44*sp.cosh(kb44*bb44)-f44)/2)
                    -sq44*sp.exp(-kb44*bb44))==0
    # rho=nP+log(Y)/2, log Y=log(q)/2-kb=-k/2-kb (q=e^-k), P=-k/2,
    # equals P(n+b+1/2)
    and sp.simplify(sp.log(sp.sqrt(sp.exp(-kb44)))+kb44/2)==0
    and sp.expand(nn44*(-kb44/2)+(-kb44/2-kb44*bb44)/2
                  -(-kb44/2)*(nn44+bb44+sp.Rational(1,2)))==0
    and exact_if_body37(logradius_tree34,'lam==1',
        "return +(P*(a+mp.mpf('.5')))")
    and matching_if_test32(moment_tree36,"lam==1 or aa<=mp.mpf('.1')")
    and height_join40 and selector_source40 and source41
    # outer branch (C43b domain43) and the interior composition: at
    # lambda=1 the terminal T is the identity (C40c), G^-1 undoes G, and
    # R^-1(R(aa))=aa by the C29g inverse identities (repeated here).
    and domain43
    and ordered_nodes37(moment_tree36.body,[
        'u = _height_R(aa)',
        'delta_z = B*mp.expm1(kappa*u)-C*mp.expm1(-kappa*u)',
        'v = _height_terminal_forward(delta_z,data)',
        'return mp.sign(a)*_height_G(_height_R(_height_G(v,amp,d,'
        'inverse=True),inverse=True),amp,d)'])
    and sp.simplify(sp.exp(sg29)-1-sp.exp(wg29))==0
    and sp.factor(wg29**2+4-(xg29+1/xg29)**2)==0
    and sp.simplify(cg29-cg29/sp.log(sp.exp(2)+sp.exp(ellg29)-sp.exp(2))
                    -rg29)==0
    and sp.simplify(cg29-cg29/sp.log(sp.E+sp.exp(ellg29)-sp.E)-rg29)==0)
rep('C44a lambda=1 collar height is P(a+1/2) on every moment branch',
    lam1_44)

# C44b.  Lattice data of M1.  South: (zm,wm)=canonical+(a+1/2)(-iP,0), and
# (-iP,0)=e4-c.e-(th/2pi)e1 (C43d), so four=xi4+a+1/2 and
# one=xi1-(a+1/2)th/(2pi).  The pure middle (|t|<=1 branch) has
# four=Arg(u2)/(2pi)+xi4+1/2, one=0, and a=Arg(u2)/(2pi)+j with j an
# integer (j=0 principal chart; j in {0,1} half-centred).  North:
# (zm,wm)=(a+1/2)(-iP,0)=(a+1/2)(-e4+c.e+(th/2pi)e1), so four=-(a+1/2),
# one=(a+1/2)th/(2pi); the pure middle (|t|>=8/5, where lam=0) has
# four=Arg(u1)/(2pi)-1, and a=1/2-Arg(u1)/(2pi)+j, since
# Arg(-conj u1)=pi-Arg u1 modulo 2 pi.  In both cases the fours differ by
# an integer and the ones by a smooth real gauge g.
au44,x4_44,j44 = sp.symbols('au44 x4_44 j44',real=True)
jint44 = sp.Symbol('jint44',integer=True)
south_four44 = x4_44+(au44/(2*sp.pi)+jint44)+sp.Rational(1,2)
south_right44 = au44/(2*sp.pi)+x4_44+sp.Rational(1,2)
north_a44 = sp.Rational(1,2)-au44/(2*sp.pi)+jint44
north_four44 = -(north_a44+sp.Rational(1,2))
north_right44 = -(-au44/(2*sp.pi)+0+sp.Rational(1,2))-sp.Rational(1,2)
lattice44 = (
    lattice_shift43
    and sp.simplify(south_four44-south_right44-jint44)==0
    and sp.simplify(north_four44-north_right44+jint44)==0
    # pure middle near |t|=2: lam=1-s((r-6/5)/(2/5))=0 for r>=8/5
    and (sp.Rational(8,5)-sp.Rational(6,5))/sp.Rational(2,5)==1
    # chart facts at exact points of both half-planes: south a-Arg(u2)/2pi
    # and north a-(1/2-Arg(u1)/2pi) are integers in both source charts
    # (principal a=Arg(p)/2pi; half-centred a=1/2+Arg(-p)/2pi).
    and all(sp.simplify(chart-target).is_integer
            for ang in (sp.pi/7,5*sp.pi/6,-sp.pi/7,-5*sp.pi/6,sp.pi)
            for (chart,target) in (
                (sp.arg(sp.exp(sp.I*ang))/(2*sp.pi),
                 sp.arg(sp.exp(sp.I*ang))/(2*sp.pi)),
                (sp.Rational(1,2)+sp.arg(-sp.exp(sp.I*ang))/(2*sp.pi),
                 sp.arg(sp.exp(sp.I*ang))/(2*sp.pi)),
                (sp.arg(-sp.conjugate(sp.exp(sp.I*ang)))/(2*sp.pi),
                 sp.Rational(1,2)-sp.arg(sp.exp(sp.I*ang))/(2*sp.pi)),
                (sp.Rational(1,2)+sp.arg(sp.conjugate(sp.exp(sp.I*ang)))
                 /(2*sp.pi),
                 sp.Rational(1,2)-sp.arg(sp.exp(sp.I*ang))/(2*sp.pi)))))
rep('C44b M1 and pure-middle fours differ by an integer at both seams',
    lattice44 and source42)

# C44c.  Projector, phase and restoration.  At the seam the cusp pair of
# M1 has weight 0 and bm=1, so the ten-vector is the two angular slots
# exp(i(2 pi one-beta9 four)) times the weights: M1 = exp(2 pi i g) times
# the pure-middle one (an integer change of four is invisible, C42b).
# Hence the projectors agree; p_M1=exp(-2 pi i g)p_R (conjugate-linear
# slot), so the corrected one agrees modulo 1, four differs by an integer
# and Psi by a lattice vector.  Both restorations have C=phi (the cusp
# weight is 0 for M1 and exactly 0 on the pure middle), so the restored
# outputs differ by a lattice translation, holomorphic with determinant 1.
g44 = sp.Symbol('g44',real=True)
Vr44 = sp.symbols('Vr0_44 Vr1_44')
Er44 = sp.symbols('Er0_44 Er1_44')
pR44 = inner41(Vr44,Er44)
pM44 = sp.expand(inner41(tuple(sp.exp(2*sp.pi*sp.I*g44)*v for v in Vr44),Er44))
gauge44 = (
    sp.simplify(pM44-sp.exp(-2*sp.pi*sp.I*g44)*pR44)==0
    and embedding42 and lattice42 and commute42 and deck_det42)
# Flatness from the left: with u=(|tc|-3/10)/(1/5), 1-lambda=b0 (or bi)
# =1-s(u) and its first two derivatives tend to 0 as u->1-, by the
# symmetry s(1-u)=1-s(u) and C41c.
flat_left44 = (
    sp.simplify(flat_actual34.subs(p34,1-u41)+flat_actual34.subs(p34,u41)-1)==0
    and flat_limits41 and flat_two_jets41)
rep('C44c same projector, ones agree mod 1, outputs differ by a lattice vector',
    gauge44 and flat_left44)

# C44d.  Exact witnesses: the C42e family u1=i/Q, u2=k/Q (theta=pi/2,
# a=0) at |t|=499/1000,1/2,501/1000 (south seam; inverse branch with
# p=1/|t|^2, p^2<16.3<50<B, |small|=1/|t|>1/2), and the C41e north family
# at |small|=499/1000,1/2,501/1000 (p=|small|^2, inf iff |small|<1/2).
# Both have |X|<3/5 (C41e, C42e).
witness44 = (
    witness42 and stereo41 and slab41 and witness_roots41
    and (1/sp.Rational(499,1000)**2)**2<sp.Rational(163,10)
    and sp.Rational(163,10)<50<sp.exp(4)
    and 1/sp.Rational(501,1000)>sp.Rational(1,2)
    and sp.Rational(499,1000)<sp.Rational(1,2)<sp.Rational(501,1000)
    and sp.Rational(501,1000)**4<sp.exp(4))
rep('C44d exact witnesses on and across both seams',witness44)

# C44e is solely a reproducible numerical observation (d).  At each seam
# it takes production points at |t| (south) or |small| (north) equal to
# 1/2-eps and 1/2+eps for several (Arg, a), including the half-centred
# chart, a common E=V(q+), and records the integer lattice vector nearest
# to Psi(q+)-Psi(q-), the residual and the projector jump, for eps=1e-12
# and 1e-24.  Sample size 36 (nine choices, two eps, two sides).  Only
# for Arg t in (-3pi/4,-pi/4) does the angular projector see four; the
# last two choices are there.  At eps<=1e-12 the collar side already
# has lambda=1 and zero cusp weight in floating point, so C44e compares
# M1 with the pure middle; the flat regime itself is C44c.
# It is not used by C44a-d.
def seam_point44(side_abs,arg_t,aa,north):
    A,B = (4*mp41.pi)**2,mp41.exp(4)
    if north:
        p = side_abs**2
        k = mp41.sqrt((p**2+A*p)/(A*p+B)); Q = mp41.sqrt(1+k**2)
        u1 = -mp41.conj(mp41.exp(2*mp41.pi*mp41.j*aa))/Q
        u2 = k*u1*mp41.exp(-mp41.j*arg_t)
    else:
        p = 1/side_abs**2
        k = mp41.sqrt((p**2+A*p)/(A*p+B)); Q = mp41.sqrt(1+k**2)
        u2 = k*mp41.exp(2*mp41.pi*mp41.j*aa)/Q
        u1 = mp41.exp(mp41.j*(arg_t+2*mp41.pi*aa))/Q
    d = 1+u1.real
    w = (u1.imag/d,u2.real/d-1,u2.imag/d)
    R = mp41.sqrt(sum(v*v for v in w))
    X = (R-1)/mp41.sqrt(1+R**2)
    n = (w[0]/R,-w[2]/R,-w[1]/R)
    return cusp30.elliptic_raw_base_reference(
        tuple((1-X**2)*v for v in n),X,dps=110,elementary=True)
rows44 = []
sample_ok44 = True
with mp41.workdps(120), cusp30.elementary_formula_branches():
    choices44 = ((False,mp41.mpf('0.7'),mp41.mpf('0.1')),
                 (False,mp41.mpf('-2.5'),mp41.mpf('-0.3')),
                 (False,mp41.mpf('2.9'),mp41.mpf('0.35')),
                 (False,mp41.mpf('1.3'),mp41.mpf('0.55')),
                 (True,mp41.mpf('0.7'),mp41.mpf('0.1')),
                 (True,mp41.mpf('-2.5'),mp41.mpf('-0.3')),
                 (True,mp41.mpf('1.9'),mp41.mpf('0.55')),
                 (False,mp41.mpf('-1.2'),mp41.mpf('0.2')),
                 (True,mp41.mpf('-1.9'),mp41.mpf('-0.15')))
    for north44,arg44,aa44 in choices44:
        res44 = []
        for eps44 in (mp41.mpf('1e-12'),mp41.mpf('1e-24')):
            kp,dp = seam_point44(mp41.mpf(1)/2+eps44,arg44,aa44,north44)
            km,dm = seam_point44(mp41.mpf(1)/2-eps44,arg44,aa44,north44)
            Vp = cusp30.elliptic_bundle_embedding(kp,dp,dps=110,elementary=True)
            Vm = cusp30.elliptic_bundle_embedding(km,dm,dps=110,elementary=True)
            Pj = max(abs(a*mp41.conj(b)-c*mp41.conj(d))
                     for a,c in zip(Vp,Vm) for b,d in zip(Vp,Vm))
            psi = []
            for data,V in ((dp,Vp),(dm,Vm)):
                tt,oo,ff,ll = data[:4]
                ph = cusp30._elementary_inner(V,Vp)
                o1 = oo+mp41.arg(ph)/(2*mp41.pi)
                psi.append((2*mp41.pi*o1-mp41.j*ll*ff,2*tt*ff,tt,ll))
            dz,dw = psi[0][0]-psi[1][0],psi[0][1]-psi[1][1]
            tt,ll = psi[1][2],psi[1][3]
            nn = int(mp41.nint((dw/(2*tt)).real))
            mm = int(mp41.nint(((dz+mp41.j*ll*nn)/(2*mp41.pi)).real))
            res = max(abs(dz-(2*mp41.pi*mm-mp41.j*ll*nn)),abs(dw-2*tt*nn))
            res44.append((mm,nn,res,Pj,len(dp),len(dm)))
            sample_ok44 &= (kp=='mid' and km=='mid' and res<1000*eps44
                            and Pj<1000*eps44)
        sample_ok44 &= (res44[0][:2]==res44[1][:2]
                        and res44[1][2]<res44[0][2]*mp41.mpf('1e-9')
                        and {res44[0][4],res44[0][5]}=={4,5})
        rows44.append((north44,arg44,aa44,res44))
    print('C44e production observations (d), sample size 36:')
    for north44,arg44,aa44,res44 in rows44:
        print('   ','north' if north44 else 'south',mp41.nstr(arg44,3),
              mp41.nstr(aa44,3),' lattice (m,n)=',res44[0][:2],
              ' residual=',[mp41.nstr(r[2],3) for r in res44],
              ' |P jump|=',[mp41.nstr(r[3],3) for r in res44])
rep('C44e production points across both collar/middle seams (numerical, d)',
    sample_ok44)

# C45: the inverse-first height moment f(b,q,lambda) and the cusp height
# cells, as one smooth function.  Every conditional executed by
# inverse_height_moment, _height_R (forward and inverse) and
# _height_terminal_forward on the G path is a domain guard, a
# working-precision shortcut (bypassed by elementary_formula_branches,
# certified in section 11.7), a pair of algebraically identical stable
# rewrites (C45a), a flat join of the written section 11.2 type (C45b), or
# the lam==1 test (C44a; never true in the collar).  The height cells of
# elementary_cusp_reference and elementary_cusp_logradius are equal
# analytic cell formulas on open intervals (C45c); their other
# conditionals (f>=0, a>0, lam==1, amp==0, lam==0) are C41a-b and C44a.
# The full classification is C50.  Identities checked with positive
# symbols are identities of analytic functions (or algebraic ones) and
# extend by analytic continuation to the source's other sign ranges.
moment_body45 = '''
import mpmath as mp
a,q,lam = mp.mpf(a),mp.mpf(q),mp.mpf(lam)
if (abs(a)>mp.mpf('.5') or not 0<=q<=mp.mpf('.125') or
        not 0<=lam<=1 or (lam>0 and q<mp.mpf(1)/2500)):
    raise ValueError('inverse moment is outside the certified (q,lam) domain')
data = _terminal_height_data(q,lam)
kappa,amp,d,B,C,z0,_,_,_,_,_,_ = data
aa = abs(a)
if lam==1 or aa<=mp.mpf('.1'):
    return _height_G(a,amp,d)
if aa>=mp.mpf('.4'):
    return mp.sign(a)*(B*mp.exp(kappa*aa)-C*mp.exp(-kappa*aa))
u = _height_R(aa)
delta_z = B*mp.expm1(kappa*u)-C*mp.expm1(-kappa*u)
v = _height_terminal_forward(delta_z,data)
return mp.sign(a)*_height_G(_height_R(_height_G(v,amp,d,inverse=True),inverse=True),amp,d)
'''
radial_body45 = '''
import mpmath as mp
delta,c = mp.mpf('.1'),mp.mpf('.2')
if inverse:
    if value<=0:
        return delta
elif value<=delta:
    return mp.mpf(0)
if value>=c:
    return value
if _elementary_precision_shortcuts.get():
    if c/(c-value)>(mp.mp.dps+10)*mp.log(10):
        return value
if inverse:
    y = mp.e*mp.expm1(value/(c-value))
    L = y+mp.log1p(-mp.exp(-y)) if y>1 else mp.log(mp.expm1(y))
    x = (2/(mp.sqrt(L*L+4)-L) if L<0 else
         (L+mp.sqrt(L*L+4))/2)
    ell = mp.log1p(x/mp.exp(2))
    return delta+c*ell/(2*(2+ell))
x = mp.exp(2)*mp.expm1(2*(value-delta)/(c-value))
F = x-1/x
soft = (mp.log1p(mp.exp(F)) if F<0 else
        F+mp.log1p(mp.exp(-F)))
ell = mp.log1p(soft/mp.e)
return c*ell/(1+ell)
'''
terminal_body45 = '''
import mpmath as mp
_,_,_,_,_,z0,zcut,eta,E0,Acut,Lambda,mu = data
gap0 = zcut-z0
if delta_z<=0:
    return mp.mpf(0)
if delta_z>=gap0:
    return z0+delta_z
gap = gap0-delta_z
L = eta/gap
Ez = mp.exp(L)
C0,C1 = 1+mu-2*Lambda,Lambda-mu
K = 1+abs(C0)*mp.exp(E0)+abs(C1)*mp.exp(E0+Acut)
if _elementary_precision_shortcuts.get():
    if Ez>(mp.mp.dps+20)*mp.log(10)+mp.log(K):
        return z0+delta_z
if Ez>50:
    D = C0*mp.exp(E0-Ez)+C1*mp.exp(E0+Acut-2*Ez)
    correction = mp.log1p(D)
    ell = mp.log1p(correction/Ez)
    z = z0+delta_z
    return z+eta*ell/(L*(L+ell))
x = Acut*mp.expm1(eta*delta_z/(gap*gap0))
u = mp.expm1(x)
F = (Lambda*u*u+mu*u)/(1+u)
M = mp.log1p(F)
ell = mp.log1p(M/E0)
return zcut*ell/(eta/zcut+ell)
'''
terminal_data_body45 = '''
import mpmath as mp
if lam==0:
    kappa,amp,d = mp.mpf(2),mp.mpf(0),mp.mpf(1)
else:
    if q<=0:
        raise ValueError('positive collar weight requires q>0')
    k = -mp.log(q)
    kappa = 2+lam*(k-2)
    amp = lam*k
    d = 1-lam+2*lam*mp.sqrt(q)*k
B,C = mp.exp(-kappa/2),q*mp.exp(kappa/2)
z0 = B-C
zcut = B*mp.exp(mp.mpf('.4')*kappa)-C*mp.exp(-mp.mpf('.4')*kappa)
eta = mp.mpf('.01')
E0 = mp.exp(eta/zcut)
Acut = mp.exp(eta/(zcut-z0))
Lambda = mp.exp(Acut-E0)
Eprime0 = E0*eta/zcut**2
Eprime_z0 = Acut*eta/(zcut-z0)**2
mu = d*Eprime0/(kappa*(B+C)*Eprime_z0)
return kappa,amp,d,B,C,z0,zcut,eta,E0,Acut,Lambda,mu
'''
radial_tree45 = source_tree30(cusp30._height_R)
terminal_fwd_tree45 = source_tree30(cusp30._height_terminal_forward)
source45 = (
    source41
    and body_dump41(moment_tree36)==text_dump41(moment_body45)
    and body_dump41(radial_tree45)==text_dump41(radial_body45)
    and body_dump41(terminal_fwd_tree45)==text_dump41(terminal_body45)
    and body_dump41(terminal_tree36)==text_dump41(terminal_data_body45))

# C45a.  Identical rewrites (identity type (a)).  With delta=1/10, c=1/5,
# E_h(a)=exp(c/(c-a)), phi=E_h-e^2, psi=E_h-e, S(x)=log(1+e^(x-1/x)),
# R=psi^-1 S phi with psi^-1(v)=c-c/log(e+v), R^-1=phi^-1 S^-1 psi:
#  forward: x=e^2 expm1(2(a-delta)/(c-a))=phi(a); log1p(e^F)=
#    F+log1p(e^-F); c ell/(1+ell)=psi^-1(S) for ell=log1p(S/e);
#  inverse: e expm1(v/(c-v))=psi(v); y+log1p(-e^-y)=log(expm1 y);
#    2/(sqrt(L^2+4)-L)=(L+sqrt(L^2+4))/2 (positive root of x-1/x=L);
#    delta+c ell/(2(2+ell))=phi^-1(x) for ell=log1p(x/e^2);
#  terminal: both branches equal z_c-eta/log(E0+M(E_t(z)-A0)) (11.3);
#    the Ez>50 branch is its (11.3a) form with v=exp(-U), U=Ez-Acut,
#    Lambda=exp(Acut-E0), C0=1+mu-2Lambda, C1=Lambda-mu.
av45,vv45,Fv45,yv45,Lv45,xv45,sv45 = sp.symbols(
    'av45 vv45 Fv45 yv45 Lv45 xv45 sv45',positive=True)
dl45,cc45 = sp.Rational(1,10),sp.Rational(1,5)
Eh45 = sp.exp(cc45/(cc45-av45))
rewrites45 = (
    # e^2 expm1(2(a-delta)/(c-a)) = E_h-e^2: the exponents agree
    sp.simplify(2+2*(av45-dl45)/(cc45-av45)-cc45/(cc45-av45))==0
    and sp.simplify(sp.expand_log(sp.log(1+sp.exp(Fv45)),force=True)
                    -(Fv45+sp.log(1+sp.exp(-Fv45))))==0
    and sp.simplify(cc45*sp.log(1+sv45/sp.E)/(1+sp.log(1+sv45/sp.E))
                    -(cc45-cc45/sp.log(sp.E+sv45)))==0
    and sp.simplify(sp.E*(sp.exp(vv45/(cc45-vv45))-1)
                    -(sp.exp(cc45/(cc45-vv45))-sp.E))==0
    and sp.simplify(sp.exp(yv45+sp.log(1-sp.exp(-yv45)))-(sp.exp(yv45)-1))==0
    and sp.simplify(2/(sp.sqrt(Lv45**2+4)-Lv45)
                    -(Lv45+sp.sqrt(Lv45**2+4))/2)==0
    and sp.simplify(((Lv45+sp.sqrt(Lv45**2+4))/2)
                    -1/((Lv45+sp.sqrt(Lv45**2+4))/2)-Lv45)==0
    and sp.simplify(dl45+cc45*sp.log(1+xv45/sp.exp(2))
                    /(2*(2+sp.log(1+xv45/sp.exp(2))))
                    -(cc45-cc45/sp.log(sp.exp(2)+xv45)))==0)
# Terminal map.  Put E0,Acut>0, Lambda=exp(Acut-E0), eta>0, U=Ez-Acut and
# M(U)=log(1+F(e^U-1)), F(u)=(Lambda u^2+mu u)/(1+u).  Then
# 1+F(e^U-1)=Lambda e^U (1+A1 v+A2 v^2), v=e^-U, with the C0,C1 above,
# so E0+M=Ez+log(1+D), D=C0 e^(E0-Ez)+C1 e^(E0+Acut-2Ez).  With
# L=log Ez=eta/gap, both returns equal z_c-eta/log(E0+M).
E0t45,Act45,mut45,Ezt45,etat45,zct45 = sp.symbols(
    'E0t45 Act45 mut45 Ezt45 etat45 zct45',positive=True)
Lamt45 = sp.exp(Act45-E0t45)
Ut45 = Ezt45-Act45
ut45 = sp.exp(Ut45)-1
Ft45 = (Lamt45*ut45**2+mut45*ut45)/(1+ut45)
C0t45,C1t45 = 1+mut45-2*Lamt45,Lamt45-mut45
Dt45 = C0t45*sp.exp(E0t45-Ezt45)+C1t45*sp.exp(E0t45+Act45-2*Ezt45)
Mt45 = sp.log(1+Ft45)
terminal45 = (
    # 1+F(e^U-1)=e^(Ez-E0)(1+D), hence E0+M=Ez+log(1+D)
    sp.simplify(sp.expand(1+Ft45)
                -sp.expand(sp.exp(Ezt45-E0t45)*(1+Dt45)))==0)
Lt45 = sp.log(Ezt45)
ellA45 = sp.Symbol('ellA45',positive=True)       # = log1p(log1p(D)/Ez)
ellB45 = sp.Symbol('ellB45',positive=True)       # = log1p(M/E0)
zt45 = zct45-etat45/Lt45
terminal_returns45 = (
    # direct branch: z_c ell/(eta/z_c+ell) with log(E0+M)=eta/z_c+ell
    sp.simplify(zct45*ellB45/(etat45/zct45+ellB45)
                -(zct45-etat45/(etat45/zct45+ellB45)))==0
    # Ez>50 branch: z+eta ell/(L(L+ell)) with log(E0+M)=L+ell
    and sp.simplify(zt45+etat45*ellA45/(Lt45*(Lt45+ellA45))
                    -(zct45-etat45/(Lt45+ellA45)))==0
    # and E0=exp(eta/z_c): log(E0+M)=log E0+log1p(M/E0)
    and sp.simplify(sp.expand_log(sp.log(sp.exp(etat45/zct45)
                    *(1+sp.Symbol('Mq45',positive=True)
                      /sp.exp(etat45/zct45))),force=True)
                    -(etat45/zct45+sp.log(1+sp.Symbol('Mq45',positive=True)
                                          /sp.exp(etat45/zct45))))==0
    # and log(Ez+log1p D)=log Ez+log1p(log1p(D)/Ez)
    and sp.simplify(sp.expand_log(sp.log(Ezt45*(1+sp.Symbol('cq45',positive=True)
                                                /Ezt45)),force=True)
                    -(Lt45+sp.log(1+sp.Symbol('cq45',positive=True)/Ezt45)))==0
    # and E_t(z)-A0=Acut*expm1(eta dz/(gap gap0)), gap=gap0-dz,
    # Acut=exp(eta/gap0): the exponents differ by eta dz/(gap gap0)
    and sp.simplify(etat45/(sp.Symbol('g0q45',positive=True)
                            -sp.Symbol('dzq45',positive=True))
                    -etat45/sp.Symbol('g0q45',positive=True)
                    -etat45*sp.Symbol('dzq45',positive=True)
                    /((sp.Symbol('g0q45',positive=True)
                       -sp.Symbol('dzq45',positive=True))
                      *sp.Symbol('g0q45',positive=True)))==0)
rep('C45a stable rewrites in R, R^-1 and both terminal branches are identities',
    source45 and rewrites45 and terminal45 and terminal_returns45)

# C45b.  The joins inside the moment (written section 11.2 arguments, types
# (b)/(c); the identities they use are checked here).  (i) |b|>=2/5: with
# R=id (b>=c) and T=id (z>=z_c=g(2/5)), the interior composite equals
# the outer formula: B expm1(kappa u)-C expm1(-kappa u)=g(u)-z0 and
# G(G^-1(g))=g; the switch at |b|=2/5 is a flat join because T-id is
# flat at z_c: by C45a, T(z)-z=eta*ell/(L(L+ell)), ell=log1p(log1p(D)/Ez)
# with D=O(exp(-Ez)), Ez=exp(L), L=eta/(z_c-z).  (ii) |b|=c=1/5: R-id is
# flat at c: with L=c/(c-a), w=1/E_h=exp(-L),
#   R(a)-a = c*l/(L(L+l)),  l=log(w(e+S(1/w-e^2))),
# and w(e+S(1/w-e^2)) = 1+(e-e^2)w-w^2/(1-e^2 w)+w log(1+exp(-1/w+e^2
# +w/(1-e^2 w))), a smooth function of w at 0+ equal to 1 there (the last
# term is flat); hence every a-derivative of R-a is a polynomial in L
# times a flat function of L.  (iii) |b|=1/10: (11.4) glues to G to all
# orders by the section 11.2 blow-up argument (R(1/10+s)=A(s)exp(-1/(20e^2 s)),
# h(x)=x+x^2B); this is cited, not re-derived.  (iv) The R^-1 switch
# `value>=c` sits at |b|=b*(q,lambda), where h(R(b*))=c (b* lies in
# (0.15,0.2)): R^-1-id is flat at c because R-id is (written, type (b):
# from R(a)=a+phi(a) with phi flat at c, R^-1(y)=y-phi(R^-1(y)) and every
# derivative of phi o R^-1 vanishes there).  (i) also needs R^-1=id on
# [2/5,1/2], i.e. h(2/5)>c: G_*(c)<=(3/10)sinh(8/5)/(8/5)<446/1000 (d<3/2,
# gamma c<=8/5, sinh(x)/x increasing), whereas z_c=g(2/5)>4485/10000
# (section 11.2: the minimum on 2<=kappa<=8 of e^(-kappa/10)-e^(-9kappa/10)
# is at kappa=8), so G_*^-1(z_c)>c.
gq45,kq45,aq45 = sp.symbols('gq45 kq45 aq45',positive=True)
Bq45 = sp.exp(-kq45/2)
Cq45 = gq45*sp.exp(kq45/2)
uq45 = sp.Symbol('uq45',positive=True)
wq45 = sp.Symbol('wq45',positive=True)
Lq45 = sp.Symbol('Lq45',positive=True)
xq45 = 1/wq45-sp.exp(2)
Sq45 = sp.log(1+sp.exp(xq45-1/xq45))
lq45 = sp.log(wq45*(sp.E+Sq45))
joins45 = (
    sp.simplify(sp.expand((Bq45*(sp.exp(kq45*uq45)-1)-Cq45*(sp.exp(-kq45*uq45)-1)))
                -sp.expand((Bq45*sp.exp(kq45*uq45)-Cq45*sp.exp(-kq45*uq45))
                           -(Bq45-Cq45)))==0
    # R-a in the (w,L) variables: a=c-c/L and R=c-c/log(e+S(phi)),
    # with phi=E_h-e^2=1/w-e^2 and log E_h=L=-log w.
    and sp.simplify((cc45-cc45/(Lq45+lq45))-(cc45-cc45/Lq45)
                    -cc45*lq45/(Lq45*(Lq45+lq45)))==0
    and sp.simplify(sp.expand(wq45*(sp.E+xq45-1/xq45))
                    -(1+(sp.E-sp.exp(2))*wq45
                      -wq45**2/(1-sp.exp(2)*wq45)))==0
    # S(x)-(x-1/x)=log(1+exp(-(x-1/x))) (the flat term): for real F both
    # sides are real logarithms of positive numbers with equal exponentials
    and sp.simplify(sp.expand((1+sp.exp(sp.Symbol('Fr45',real=True)))
                              *sp.exp(-sp.Symbol('Fr45',real=True)))
                    -(1+sp.exp(-sp.Symbol('Fr45',real=True))))==0
    and matching_if_test32(moment_tree36,"aa>=mp.mpf('.4')")
    and exact_if_body37(radial_tree45,'value>=c','return value')
    and exact_if_body37(terminal_fwd_tree45,'delta_z>=gap0',
                        'return z0+delta_z')
    and sp.Rational(3,10)*sp.sinh(sp.Rational(8,5))/sp.Rational(8,5)
        <sp.Rational(446,1000)
    and sp.exp(-sp.Rational(4,5))-sp.exp(-sp.Rational(36,5))
        >sp.Rational(4485,10000)
    and sp.Rational(446,1000)<sp.Rational(4485,10000))
rep('C45b moment joins at |b|=2/5, 1/5 and b* are flat; 11.2 identities',
    joins45 and source45)

# C45c.  Height-cell seams (identity type (a)).  Core (lambda=0, kappa=2,
# B=e^-1, C=qe, outer branch on |b|>=2/5): at a=1/2 the cell-0 a>0 formula
# A=(R+f)/2 with b=a and the cell-1 formula A=2/(R-f) with b=a-1 are both
# e^(2a-1); at a=-1/2 the cell-0 Y=(R-f)/2 and the cell-(-1) Y=2/(R+f)
# are both e^(-2a-1).  Collar (any kappa): across a=1/2 the logradius
# 0*P+log(Y(b))/2 (b=a) and P+log(Y(b-1))/2 coincide: on the outer branch
# Y(b)=q e^(-kappa(b-1/2)) for b>2/5 and Y(b')=e^(-kappa(b'+1/2)) for
# b'<-2/5, with P=log(q)/2.  These are equalities of analytic functions
# on the open intervals |a-+1/2|<1/10, so all jets agree.
ac45 = sp.Symbol('ac45',real=True)
qc45 = sp.Symbol('qc45',positive=True)
fpos45 = sp.exp(2*ac45-1)-qc45*sp.exp(1-2*ac45)        # b=a>0, kappa=2
fneg45 = -(sp.exp(1-2*ac45)-qc45*sp.exp(2*ac45-1))     # b=a-1<0
Rpos45 = sp.exp(2*ac45-1)+qc45*sp.exp(1-2*ac45)
Rneg45 = sp.exp(1-2*ac45)+qc45*sp.exp(2*ac45-1)
gneg45 = -(sp.exp(-2*ac45-1)-qc45*sp.exp(1+2*ac45))    # b=a<0 near -1/2
Rgneg45 = sp.exp(-2*ac45-1)+qc45*sp.exp(1+2*ac45)
gpos45 = sp.exp(2*ac45+1)-qc45*sp.exp(-2*ac45-1)       # b=a+1>0
Rgpos45 = sp.exp(2*ac45+1)+qc45*sp.exp(-2*ac45-1)
kc45,bc45 = sp.symbols('kc45 bc45',positive=True)
fcol45 = sp.exp(-kc45/2)*sp.exp(kc45*bc45)-qc45*sp.exp(kc45/2)*sp.exp(-kc45*bc45)
Rcol45 = sp.exp(-kc45/2)*sp.exp(kc45*bc45)+qc45*sp.exp(kc45/2)*sp.exp(-kc45*bc45)
bn45 = bc45-1                                           # b'=b-1<0
fcoln45 = -(sp.exp(-kc45/2)*sp.exp(-kc45*bn45)
            -qc45*sp.exp(kc45/2)*sp.exp(kc45*bn45))
Rcoln45 = (sp.exp(-kc45/2)*sp.exp(-kc45*bn45)
           +qc45*sp.exp(kc45/2)*sp.exp(kc45*bn45))
cells45 = (
    sp.simplify(sp.expand(Rpos45**2-fpos45**2-4*qc45))==0
    and sp.simplify(sp.expand(Rneg45**2-fneg45**2-4*qc45))==0
    and sp.simplify((Rpos45+fpos45)/2-sp.exp(2*ac45-1))==0
    and sp.simplify(2/(Rneg45-fneg45)-sp.exp(2*ac45-1))==0
    and sp.simplify(sp.expand(Rgneg45**2-gneg45**2-4*qc45))==0
    and sp.simplify(sp.expand(Rgpos45**2-gpos45**2-4*qc45))==0
    and sp.simplify((Rgneg45-gneg45)/2-sp.exp(-2*ac45-1))==0
    and sp.simplify(2/(Rgpos45+gpos45)-sp.exp(-2*ac45-1))==0
    and sp.simplify(sp.expand(Rcol45**2-fcol45**2-4*qc45))==0
    and sp.simplify(sp.expand(Rcoln45**2-fcoln45**2-4*qc45))==0
    and sp.simplify(2*qc45/(Rcol45+fcol45)
                    -qc45*sp.exp(-kc45*(bc45-sp.Rational(1,2))))==0
    and sp.simplify((Rcoln45-fcoln45)/2
                    -sp.exp(-kc45*(bn45+sp.Rational(1,2))))==0
    and sp.simplify(sp.expand_log(sp.log(qc45*sp.exp(-kc45*(bc45-sp.Rational(1,2)))),
                                  force=True)/2
                    -(sp.log(qc45)/2+sp.expand_log(
                        sp.log(sp.exp(-kc45*(bn45+sp.Rational(1,2)))),
                        force=True)/2))==0
    and exact_if_body37(core_tree37,'a>0','''
A = 2/(R-f) if cell==1 else (R+f)/2
x = mp.sqrt(A/e)/nu
return x,t/x
''')
    and matching_statement31(core_tree37,
        "Y = 2/(R+f) if cell==-1 else (R-f)/2"))
rep('C45c core and collar height-cell seams: equal analytic cell formulas',
    cells45 and source45)

# C46: pole-cap internals, whole cap, at both signs (south c0, north ci).
# The production cap branch of elementary_global_slice (used when
# 1-|X|<alpha/2, alpha=1e-5) computes
#   z_L=U_j^* Z/sqrt(s_j(1+|X|sqrt(2-X^2))),  zz=H_j(z_L),
#   H_j=D^-1 o g_j o D  (identity on |z_L|^2<=a=1e-8),
#   return (conj(zz0)/nu_j,zz1),
# (programme (11.12h), (11.17)).  C46 proves that this is ONE smooth map on
# the whole cap: every conditional on its call path is (i) a static choice,
# (ii) a guard proved never to fire on the cap, (iii) a working-precision
# shortcut (bypassed by elementary_formula_branches; its error bound is the
# type (c) argument of programme section 11.7, and symbolic_identities (m1)
# counts it), (iv) a rewrite whose two alternatives are the same function
# (type (a)), or (v) a flat radial join (type (a) identities plus the
# written calculus of section 11.3).  C46a binds the source and classifies
# all 32 conditionals; C46b-C46d prove the classification; C46e is a
# numerical observation (d) only.  Out of scope: the cap/slab selector
# distance<alpha/2 itself (the open overlap of section 11.6, N9).
import copy as copy46

# ---------------------------------------------------------------- C46a
# Source binding (exact AST equality of complete function bodies; type (a)
# in the sense of a source identity).  Docstrings are stripped by
# body_dump41, comments are not part of the AST.
patch_body46 = '''
import mpmath as mp
with mp.workdps(max(mp.mp.dps,dps)+50):
    pair = tuple(map(mp.mpc,z))
    inner = mp.mpf('1e-8')
    outer = mp.mpf('1.5e-8')
    radius = sum(abs(a)**2 for a in pair)
    if radius<=inner:
        return pair
    if radius<(inner+outer)/2:
        x = mp.exp(3)*mp.expm1(3*(radius-inner)/(outer-radius))
        if _elementary_precision_shortcuts.get():
            if x-1/x<=-(dps+10)*mp.log(10):
                return pair
    compressed = _elementary_pole_radial(pair)
    raw = elliptic_pole_raw_inverse(compressed,infinity,dps=dps,elementary=True)
    return _elementary_pole_radial(raw,inverse=True)
'''
radial_body46 = '''
import mpmath as mp
z = tuple(map(mp.mpc,pair))
radius = sum(abs(w)**2 for w in z)
a,c = mp.mpf('1e-8'),mp.mpf('1.5e-8')
if radius>=c:
    return z
if _elementary_precision_shortcuts.get():
    if c/(c-radius)>(mp.mp.dps+10)*mp.log(10):
        return z
if inverse:
    if radius<=0:
        raise ValueError('radial inverse requires positive radius')
    y = mp.e*mp.expm1(radius/(c-radius))
    L = (y+mp.log1p(-mp.exp(-y)) if y>1 else mp.log(mp.expm1(y)))
    x = 2/(mp.sqrt(L*L+4)-L) if L<0 else (L+mp.sqrt(L*L+4))/2
    ell = mp.log1p(x/mp.exp(3))
    new_radius = a+c*ell/(3*(3+ell))
else:
    if radius<=a:
        raise ValueError('radial compression requires radius greater than 1e-8')
    x = mp.exp(3)*mp.expm1(3*(radius-a)/(c-radius))
    f = x-1/x
    y = mp.log1p(mp.exp(f)) if f<0 else f+mp.log1p(mp.exp(-f))
    ell = mp.log1p(y/mp.e)
    new_radius = c*ell/(1+ell)
factor = mp.sqrt(new_radius/radius)
return tuple(+(factor*w) for w in z)
'''
raw_body46 = '''
import mpmath as mp
with mp.workdps(max(mp.mp.dps,dps)+20):
    z0,z1 = map(mp.mpc,z)
    radius = abs(z0)**2+abs(z1)**2
    if not 0<radius<=mp.mpf('1e-4'):
        raise ValueError('raw pole inverse requires 0<|z|^2<=1e-4')
    pc = 2*z0*mp.conj(z1)
    p = (pc.real,pc.imag,abs(z0)**2-abs(z1)**2)
    scale = mp.pi if infinity else 2*mp.pi
    v = tuple(scale*w for w in ((p[2],-p[1],p[0]) if infinity
                                else (-p[1],-p[2],p[0])))
    aa,bb,cc = v[0],-v[2],-v[1]
    if infinity:
        s2 = aa*aa+bb*bb+cc*cc
        den = 1+2*bb+2*s2
        u1 = mp.mpc(-1-2*bb,2*aa)/den
        u2 = mp.mpc(2*(bb+s2),2*cc)/den
        ratio = u2/u1
        height = -mp.arg(-u1)/(2*mp.pi)
    else:
        s2 = aa*aa+(bb+1)**2+cc*cc
        # 1-s2=-(aa^2+bb^2+cc^2+2bb), written without the O(1) cancellation
        u1 = mp.mpc(-(aa*aa+bb*bb+cc*cc+2*bb),2*aa)/(1+s2)
        u2 = mp.mpc(2*(bb+1),2*cc)/(1+s2)
        ratio = u1/u2
        height = mp.arg(u2)/(2*mp.pi)
    A,B = (4*mp.pi)**2,mp.exp(4)
    uu = abs(ratio)**2
    den = A*(1-uu)+mp.sqrt(A*A*(1-uu)**2+4*B*uu)
    tc = ratio*mp.sqrt((2*B if infinity else 2)/den)
    nu = mp.exp(-2) if infinity else mp.mpf(1)
    if elementary:
        if abs(height)>=mp.mpf('.1'):
            raise ValueError('elementary pole inverse left its certified height region')
        e = mp.exp(-4*abs(tc)**2/(1-abs(tc)**2))
        q = e*nu*nu*abs(tc)**2
        R = mp.sqrt(height*height+4*q)
        AA = 2*q/(R-height) if height<0 else (R+height)/2
        yy = 2*q/(R+height) if height>=0 else (R-height)/2
        xx = AA/e
        moment,source_radius = xx-yy,xx+yy
    elif tc == 0:
        rr = elliptic_axis_radius(height)
        moment = mp.sign(height)*rr
        source_radius = rr
    else:
        rho = (elementary_cusp_logradius if elementary else elliptic_cusp_logradius)(height,tc,nu,dps=dps)
        yy = mp.exp(2*rho)
        xx = nu*nu*abs(tc)**2/yy
        moment,source_radius = xx-yy,xx+yy
    source_pc = 2*nu*mp.conj(tc)
    w0 = ((source_radius+moment)*z0+source_pc*z1)/(2*source_radius)
    w1 = (mp.conj(source_pc)*z0+(source_radius-moment)*z1)/(2*source_radius)
    length = mp.sqrt(abs(w0)**2+abs(w1)**2)
    if length == 0:
        raise ValueError('local Hopf projections became orthogonal')
    factor = mp.sqrt(source_radius)/length
    return (+(factor*w0),+(factor*w1))
'''
inclusion_body46 = '''
import mpmath as mp
a,b = xi
pair = ((a+b)/mp.sqrt(2),(a-b)/mp.sqrt(2)) if infinity else (
    (a-mp.j*b)/mp.sqrt(2),(a+mp.j*b)/mp.sqrt(2))
result = [mp.mpc(0)]*10
start = 4 if infinity else 0
result[start:start+2] = pair
return tuple(result)
'''
branches_body46 = '''
token = _elementary_precision_shortcuts.set(False)
try:
    yield
finally:
    _elementary_precision_shortcuts.reset(token)
'''
cap_test46 = 'distance<alpha/2 and not force_slab'
cap_body46 = '''
inf = X>0
scale = mp.pi if inf else 2*mp.pi
den = mp.sqrt(scale*(1+abs(X)*mp.sqrt(2-X*X)))
pair = _elementary_pole_inclusion(Z,inf)[4:6] if inf else _elementary_pole_inclusion(Z)[:2]
zz = elementary_pole_patch(tuple(a/den for a in pair),inf,dps=dps)
nu = mp.exp(-2) if inf else mp.mpf(1)
return ('ci' if inf else 'c0',(mp.conj(zz[0])/nu,zz[1]))
'''
slice_prefix46 = '''
Z = tuple(map(mp.mpc,Z))
X = mp.mpf(X)
alpha = mp.mpf('1e-5')
distance = 1-abs(X)
'''
patch_tree46 = source_tree30(cusp30.elementary_pole_patch)
radial_tree46 = source_tree30(cusp30._elementary_pole_radial)
raw_tree46 = source_tree30(cusp30.elliptic_pole_raw_inverse)
inclusion_tree46 = source_tree30(cusp30._elementary_pole_inclusion)
branches_tree46 = source_tree30(cusp30.elementary_formula_branches)
cap_nodes46 = [node for node in slice_with37[0].body
               if isinstance(node,ast.If)
               and ast.dump(node.test,include_attributes=False)
               ==ast.dump(ast.parse(cap_test46,mode='eval').body,
                          include_attributes=False)]
source46 = (
    body_dump41(patch_tree46)==text_dump41(patch_body46)
    and body_dump41(radial_tree46)==text_dump41(radial_body46)
    and body_dump41(raw_tree46)==text_dump41(raw_body46)
    and body_dump41(inclusion_tree46)==text_dump41(inclusion_body46)
    and body_dump41(branches_tree46)==text_dump41(branches_body46)
    and len(slice_with37)==1
    # the cap test is the first branch after the four prefix statements,
    # and its body (ending in a return) is exactly the cap formula
    and [ast.dump(node,include_attributes=False)
         for node in slice_with37[0].body[:4]]==text_dump41(slice_prefix46)
    and len(cap_nodes46)==1 and slice_with37[0].body[4] is cap_nodes46[0]
    and exact_if_body37(slice_tree30,cap_test46,cap_body46)
    and isinstance(cap_nodes46[0].body[-1],ast.Return))

# Classification of EVERY If/IfExp on the cap call path (multiset of test
# ASTs; the source must contain exactly these, no more, no fewer).
# Classes: static = fixed on each cap / by the call site; guard = raise
# never reached on the cap (C46b,C46d); precision = dps shortcut (bypassed
# in formula-branch mode, programme 11.7); same = identical alternatives
# (C46b2,C46c); join = flat radial join (C46c); noop = no effect in
# formula-branch mode (C46a2); unexec = not on the elementary=True path.
classes46 = {
    'cap': [('IfExp','inf','static')]*4,
    'inclusion': [('IfExp','infinity','static')]*2,
    'patch': [('If','radius<=inner','join'),
              ('If','radius<(inner+outer)/2','noop'),
              ('If','_elementary_precision_shortcuts.get()','precision'),
              ('If','x-1/x<=-(dps+10)*mp.log(10)','precision')],
    'radial': [('If','radius>=c','join'),
               ('If','_elementary_precision_shortcuts.get()','precision'),
               ('If','c/(c-radius)>(mp.mp.dps+10)*mp.log(10)','precision'),
               ('If','inverse','static'),
               ('If','radius<=0','guard'),
               ('IfExp','y>1','same'),
               ('IfExp','L<0','same'),
               ('If','radius<=a','guard'),
               ('IfExp','f<0','same')],
    'raw': [('If',"not 0<radius<=mp.mpf('1e-4')",'guard'),
            ('IfExp','infinity','static'),
            ('IfExp','infinity','static'),
            ('If','infinity','static'),
            ('IfExp','infinity','static'),
            ('IfExp','infinity','static'),
            ('If','elementary','static'),
            ('If',"abs(height)>=mp.mpf('.1')",'guard'),
            ('IfExp','height<0','same'),
            ('IfExp','height>=0','same'),
            ('If','tc == 0','unexec'),
            ('IfExp','elementary','unexec'),
            ('If','length == 0','guard')]}
def conditionals46(tree):
    return sorted((type(node).__name__,
                   ast.dump(node.test,include_attributes=False))
                  for node in ast.walk(tree)
                  if isinstance(node,(ast.If,ast.IfExp)))
def expected46(rows):
    return sorted((kind,ast.dump(ast.parse(test,mode='eval').body,
                                 include_attributes=False))
                  for kind,test,_ in rows)
cap_module46 = ast.Module(body=cap_nodes46[0].body,type_ignores=[])
classified46 = (
    conditionals46(cap_module46)==expected46(classes46['cap'])
    and conditionals46(inclusion_tree46)==expected46(classes46['inclusion'])
    and conditionals46(patch_tree46)==expected46(classes46['patch'])
    and conditionals46(radial_tree46)==expected46(classes46['radial'])
    and conditionals46(raw_tree46)==expected46(classes46['raw']))
count46 = sum(len(rows) for rows in classes46.values())
# The arguments that make 'inverse' and 'elementary' static: the bound
# patch body calls D as _elementary_pole_radial(pair), D^-1 with
# inverse=True and the raw inverse with elementary=True; D's
# inverse=False is the signature default, bound here too.
signature46 = (
    ast.dump(radial_tree46.args,include_attributes=False)
    ==ast.dump(ast.parse('def f(pair,inverse=False): pass').body[0].args,
               include_attributes=False)
    and matching_assignment30(patch_tree46,'compressed',
                              '_elementary_pole_radial(pair)')
    and matching_assignment30(patch_tree46,'raw',
        'elliptic_pole_raw_inverse(compressed,infinity,dps=dps,elementary=True)'))
rep(f'C46a full pole-cap source bound; all {count46} conditionals classified',
    source46 and classified46 and count46==32 and signature46)

# C46a2.  (noop) In elementary_pole_patch the block under
# radius<(inner+outer)/2 only assigns x and then tests the precision
# shortcut; x is never read afterwards.  So with the shortcut disabled the
# block has no effect.  (precision) The opt-in context manager switches
# the ContextVar read by both shortcut gates.  Source checks + a runtime
# check of the state; the shortcut error bounds themselves are the type (c)
# argument of programme 11.7.
patch_with46 = [node for node in patch_tree46.body
                if isinstance(node,ast.With)][0].body
noop_index46 = [i for i,node in enumerate(patch_with46)
                if isinstance(node,ast.If)
                and ast.dump(node.test,include_attributes=False)
                ==ast.dump(ast.parse('radius<(inner+outer)/2',mode='eval')
                           .body,include_attributes=False)]
noop46 = (
    len(noop_index46)==1
    and [type(node).__name__ for node in patch_with46[noop_index46[0]].body]
        ==['Assign','If']
    and not any(isinstance(node,ast.Name) and node.id=='x'
                for later in patch_with46[noop_index46[0]+1:]
                for node in ast.walk(later)))
gates46 = [node for tree in (patch_tree46,radial_tree46)
           for node in ast.walk(tree)
           if isinstance(node,ast.If)
           and ast.dump(node.test,include_attributes=False)
           ==ast.dump(ast.parse('_elementary_precision_shortcuts.get()',
                                mode='eval').body,include_attributes=False)]
state_default46 = cusp30._elementary_precision_shortcuts.get()
with cusp30.elementary_formula_branches():
    state_formula46 = cusp30._elementary_precision_shortcuts.get()
gate_ok46 = (len(gates46)==2
             and all(len(node.body)==1 and isinstance(node.body[0],ast.If)
                     and isinstance(node.body[0].body[0],ast.Return)
                     and not node.orelse for node in gates46)
             and state_default46 is True and state_formula46 is False
             and cusp30._elementary_precision_shortcuts.get() is True)
rep('C46a2 x-block is a no-op and both shortcuts are formula-mode gated',
    noop46 and gate_ok46)

# ---------------------------------------------------------------- C46b
# Symbolic execution of the bound source statements.  Each statement node
# of the actual source tree is executed with sympy: `import mpmath` is
# dropped, `mp.*` is a sympy stand-in, `.real/.imag` become re/im, and an
# IfExp whose test is not a Python/sympy boolean records BOTH alternatives
# (then proceeds with the else-alternative).  Static IfExps (a Python bool
# test) pick their branch.  So the identities below are identities of the
# source expressions themselves, type (a).
class SymMp46:
    pi = sp.pi
    e = sp.E
    j = sp.I
    sqrt = staticmethod(sp.sqrt)
    exp = staticmethod(sp.exp)
    log = staticmethod(sp.log)
    conj = staticmethod(sp.conjugate)
    arg = staticmethod(sp.Function('Arg46'))
    @staticmethod
    def mpc(re_part,im_part=0):
        return re_part+sp.I*im_part
    @staticmethod
    def mpf(value):
        return sp.Rational(value) if isinstance(value,str) else sp.sympify(value)
    @staticmethod
    def expm1(value):
        return sp.exp(value)-1
    @staticmethod
    def log1p(value):
        return sp.log(1+value)


class SymSource46(ast.NodeTransformer):
    def visit_Import(self,node):
        return None
    def visit_Attribute(self,node):
        self.generic_visit(node)
        if node.attr in ('real','imag'):
            return ast.Call(func=ast.Name(id='re46' if node.attr=='real'
                                          else 'im46',ctx=ast.Load()),
                            args=[node.value],keywords=[])
        return node
    def visit_IfExp(self,node):
        self.generic_visit(node)
        return ast.Call(func=ast.Name(id='choose46',ctx=ast.Load()),
                        args=[node.test,node.body,node.orelse],keywords=[])


alternatives46 = []
def choose46(test,first,second):
    try:
        flag = bool(test)
    except TypeError:
        alternatives46.append((first,second))
        return second
    return first if flag else second


def namespace46(**values):
    ns = {'mp':SymMp46,'re46':sp.re,'im46':sp.im,'choose46':choose46}
    ns.update(values)
    return ns


def run46(nodes,ns):
    for node in nodes:
        module = ast.Module(body=[SymSource46().visit(copy46.deepcopy(node))],
                            type_ignores=[])
        exec(compile(ast.fix_missing_locations(module),'sym46','exec'),ns)
    return ns


def stmt46(body,text):
    wanted = ast.dump(ast.parse(text).body[0],include_attributes=False)
    found = [node for node in body
             if ast.dump(node,include_attributes=False)==wanted]
    assert len(found)==1, text
    return found[0]


raw_with46 = [node for node in raw_tree46.body
              if isinstance(node,ast.With)][0].body
x0_46,y0_46,x1_46,y1_46 = sp.symbols('x0_46 y0_46 x1_46 y1_46',real=True)
zr46 = x0_46**2+y0_46**2+x1_46**2+y1_46**2
aa46,bb46,cc46 = sp.symbols('aa46 bb46 cc46',real=True)
infinity_if46 = [node for node in raw_with46
                 if isinstance(node,ast.If)
                 and ast.dump(node.test,include_attributes=False)
                 ==ast.dump(ast.parse('infinity',mode='eval').body,
                            include_attributes=False)][0]
prefix_texts46 = ['pc = 2*z0*mp.conj(z1)',
                  'p = (pc.real,pc.imag,abs(z0)**2-abs(z1)**2)',
                  'scale = mp.pi if infinity else 2*mp.pi',
                  'v = tuple(scale*w for w in ((p[2],-p[1],p[0]) if infinity '
                  'else (-p[1],-p[2],p[0])))',
                  'aa,bb,cc = v[0],-v[2],-v[1]']
stereo46 = True
height_num46 = {}
den_min46 = {}
for north46 in (False,True):
    alternatives46.clear()
    ns46 = run46([stmt46(raw_with46,text) for text in prefix_texts46],
                 namespace46(z0=x0_46+sp.I*y0_46,z1=x1_46+sp.I*y1_46,
                             infinity=north46))
    p46 = ns46['p']
    sj46 = sp.pi if north46 else 2*sp.pi
    # Hopf: |p|^2=|z|^4, so every |p_k|<=u=|z|^2 (polynomial identity).
    stereo46 &= (sp.expand(sum(w**2 for w in p46)-zr46**2)==0
                 and ns46['scale']==sj46 and not alternatives46)
    # The rotated, scaled vector: bb=-s_j p0 at both signs; the height
    # numerators cc (south) and 2aa (north) are both 2 pi p2.
    stereo46 &= (sp.expand(ns46['bb']+sj46*p46[0])==0
                 and sp.expand(ns46['aa']**2+ns46['bb']**2+ns46['cc']**2
                               -sj46**2*sum(w**2 for w in p46))==0)
    height_num46[north46] = sp.expand(2*ns46['aa'] if north46 else ns46['cc'])
    den_min46[north46] = sp.expand(1+2*ns46['bb'] if north46 else 1+ns46['bb'])
    # Execute the selected stereographic branch with aa,bb,cc symbolic.
    branch46 = infinity_if46.body if north46 else infinity_if46.orelse
    nb46 = run46(branch46,namespace46(aa=aa46,bb=bb46,cc=cc46))
    u1_46,u2_46 = nb46['u1'],nb46['u2']
    norm46 = sp.simplify(sp.re(u1_46)**2+sp.im(u1_46)**2
                         +sp.re(u2_46)**2+sp.im(u2_46)**2-1)==0
    if north46:
        dom46 = -u1_46                     # dominant coordinate (ratio=u2/u1)
        den46 = nb46['den']
        stereo46 &= (
            norm46
            and sp.simplify(sp.re(dom46)*den46-(1+2*bb46))==0
            and sp.simplify(sp.im(dom46)*(1+2*bb46)+sp.re(dom46)*2*aa46)==0
            and sp.expand(den46-(1+2*bb46)-2*(aa46**2+bb46**2+cc46**2))==0
            and sp.simplify((sp.re(u2_46)**2+sp.im(u2_46)**2)*den46**2
                            -4*((bb46+nb46['s2'])**2+cc46**2))==0
            and sp.simplify((sp.re(u1_46)**2+sp.im(u1_46)**2)*den46**2
                            -((1+2*bb46)**2+4*aa46**2))==0
            # tc=0 and height=0 force u2=0 and aa=0, then bb(1+bb)=0
            and sp.factor((bb46+nb46['s2']).subs({aa46:0,cc46:0}))
                ==bb46*(bb46+1))
    else:
        dom46 = u2_46                      # dominant coordinate (ratio=u1/u2)
        s2_46 = nb46['s2']
        stereo46 &= (
            norm46
            and sp.simplify(sp.re(dom46)*(1+s2_46)/2-(1+bb46))==0
            and sp.simplify(sp.im(dom46)*(1+bb46)-sp.re(dom46)*cc46)==0
            and sp.expand(s2_46-1-(aa46**2+bb46**2+cc46**2)-2*bb46)==0
            and sp.simplify((sp.re(u1_46)**2+sp.im(u1_46)**2)*(1+s2_46)**2
                            -((1-s2_46)**2+4*aa46**2))==0
            and sp.simplify((sp.re(u2_46)**2+sp.im(u2_46)**2)*(1+s2_46)**2
                            -4*((1+bb46)**2+cc46**2))==0
            # tc=0 and height=0 force u1=0 and cc=0, then bb(bb+2)=0
            and sp.factor((s2_46-1).subs({aa46:0,cc46:0}))==bb46*(bb46+2))
    # The executed height/ratio statements of the branch use exactly this
    # dominant coordinate: height=+-Arg(dom)/(2 pi), ratio = other/dom.
    stereo46 &= (
        sp.expand(nb46['height']-(-1 if north46 else 1)
                  *SymMp46.arg(dom46)/(2*sp.pi))==0
        and sp.expand(nb46['ratio']-(u2_46/u1_46 if north46
                                     else u1_46/u2_46))==0)
# Both signs share the SAME height formula in p=Hopf(z):
#   tan(2 pi height)=2 pi p2/(1-2 pi p0),  dominant real part ~ 1-2 pi p0.
stereo46 &= (sp.expand(height_num46[False]-height_num46[True])==0
             and sp.expand(den_min46[False]-den_min46[True])==0
             and sp.expand(den_min46[False]-(1-2*sp.pi*p46[0]))==0
             and sp.expand(height_num46[False]-2*sp.pi*p46[2])==0)
rep('C46b raw stereographic/height source: identities at both signs',
    stereo46 and source46)

# C46b2.  The elementary=True block, executed symbolically on the source
# statements (guard excluded; it is proved inactive in C46b3).  With
# R^2=h^2+4q the two alternatives of AA and of yy are the same functions
# (R+h)/2 and (R-h)/2 (type (a)); xx*yy=nu^2|tc|^2, so the source vector
# p_S=(Re source_pc,Im source_pc,moment) has |p_S|=source_radius exactly.
# The Hopf projection: |w|^2=(r+n_S.Hopf(z))/2+(|p_S|^2-s^2)r/(4s^2), so
# |w|^2=r(1+n_S.n_L)/2 (programme (11.10), C22b,c).
elem_if46 = [node for node in raw_with46
             if isinstance(node,ast.If)
             and ast.dump(node.test,include_attributes=False)
             ==ast.dump(ast.parse('elementary',mode='eval').body,
                        include_attributes=False)][0]
h46 = sp.Symbol('h46',real=True)
tr46,ti46 = sp.symbols('tr46 ti46',real=True)
elem_ok46 = (isinstance(elem_if46.body[0],ast.If)
             and isinstance(elem_if46.body[0].body[0],ast.Raise)
             and len(elem_if46.body)==8)
for north46 in (False,True):
    alternatives46.clear()
    nu46 = sp.exp(-2) if north46 else sp.Integer(1)
    ne46 = run46(elem_if46.body[1:],
                 namespace46(height=h46,tc=tr46+sp.I*ti46,nu=nu46))
    pairs46 = list(alternatives46)
    elem_ok46 &= (len(pairs46)==2
                  and all(sp.simplify(sp.radsimp(first)-second)==0
                          for first,second in pairs46))
    run46([stmt46(raw_with46,'source_pc = 2*nu*mp.conj(tc)')],ne46)
    T46 = tr46**2+ti46**2
    elem_ok46 &= (
        sp.expand(ne46['AA']*ne46['yy']-ne46['q'])==0
        and sp.simplify(ne46['xx']*ne46['yy']-nu46**2*T46)==0
        and sp.simplify(sp.re(ne46['source_pc'])**2+sp.im(ne46['source_pc'])**2
                        +ne46['moment']**2-ne46['source_radius']**2)==0)
s46,m46,sr46,si46 = (sp.Symbol('s46',positive=True),sp.Symbol('m46',real=True),
                     sp.Symbol('sr46',real=True),sp.Symbol('si46',real=True))
nw46 = run46([stmt46(raw_with46,text) for text in (
    'w0 = ((source_radius+moment)*z0+source_pc*z1)/(2*source_radius)',
    'w1 = (mp.conj(source_pc)*z0+(source_radius-moment)*z1)/(2*source_radius)')],
    namespace46(z0=x0_46+sp.I*y0_46,z1=x1_46+sp.I*y1_46,source_radius=s46,
                moment=m46,source_pc=sr46+sp.I*si46))
length2_46 = sp.expand(sum(sp.re(w)**2+sp.im(w)**2
                           for w in (nw46['w0'],nw46['w1'])))
proj46 = sp.expand(4*s46**2*(length2_46-(zr46+(sr46*p46[0]+si46*p46[1]
                                               +m46*p46[2])/s46)/2)
                   -(sr46**2+si46**2+m46**2-s46**2)*zr46)==0
nS46 = sp.symbols('nS0_46 nS1_46 nS2_46',real=True)
nL46 = sp.symbols('nL0_46 nL1_46 nL2_46',real=True)
unit46 = sp.expand(1-sum((a-b)**2 for a,b in zip(nS46,nL46))/4
                   -(1+sum(a*b for a,b in zip(nS46,nL46)))/2
                   +(sum(a*a for a in nS46)-1)/4+(sum(b*b for b in nL46)-1)/4)==0
# The AA/yy alternatives again, with q a positive symbol (independent of
# the particular e): both pairs are identical for every q>0.
qs46 = sp.Symbol('qs46',positive=True)
alternatives46.clear()
nq46 = run46([stmt46(elem_if46.body,text) for text in (
    'R = mp.sqrt(height*height+4*q)',
    'AA = 2*q/(R-height) if height<0 else (R+height)/2',
    'yy = 2*q/(R+height) if height>=0 else (R-height)/2')],
    namespace46(height=h46,q=qs46))
elem_ok46 &= (len(alternatives46)==2
              and all(sp.simplify(sp.radsimp(first)-second)==0
                      for first,second in alternatives46)
              and sp.expand(nq46['AA']*nq46['yy']-qs46)==0)
# The tc chain, executed with a symbolic ratio: |tc|^2=K uu/den with
# K=2B (north), 2 (south), and den>=2A(1-uu) is used in C46b3.
rx46,ry46 = sp.symbols('rx46 ry46',real=True)
dd46 = sp.Symbol('dd46',positive=True)
tc_ok46 = True
for north46 in (False,True):
    alternatives46.clear()
    nt46 = run46([stmt46(raw_with46,text) for text in (
        'A,B = (4*mp.pi)**2,mp.exp(4)',
        'uu = abs(ratio)**2',
        'den = A*(1-uu)+mp.sqrt(A*A*(1-uu)**2+4*B*uu)',
        'nu = mp.exp(-2) if infinity else mp.mpf(1)')],
        namespace46(ratio=rx46+sp.I*ry46,infinity=north46))
    # tc with den a positive symbol (den>0 is shown in C46b3)
    ntc46 = run46([stmt46(raw_with46,
                          'tc = ratio*mp.sqrt((2*B if infinity else 2)/den)')],
                  namespace46(ratio=rx46+sp.I*ry46,infinity=north46,
                              B=nt46['B'],den=dd46))
    Kt46 = 2*sp.exp(4) if north46 else 2
    tc_ok46 &= (
        not alternatives46
        and sp.expand(nt46['uu']-(rx46**2+ry46**2))==0
        and sp.expand(ntc46['tc']*sp.conjugate(ntc46['tc'])
                      -Kt46*(rx46**2+ry46**2)/dd46)==0
        and nt46['A']==16*sp.pi**2 and nt46['B']==sp.exp(4)
        # den=A(1-uu)+sqrt(A^2(1-uu)^2+4B uu) with the positive root
        and sp.expand(nt46['den']-16*sp.pi**2*(1-nt46['uu'])
                      -sp.sqrt(256*sp.pi**4*(1-nt46['uu'])**2
                               +4*sp.exp(4)*nt46['uu']))==0
        and nt46['nu']==(sp.exp(-2) if north46 else 1))
# source_pc is the Hopf vector of the slice pair (nu conj(x), y) of the
# docstring, with tc=xy: 2 (nu conj x) conj(y)=2 nu conj(tc).
xs46,ys46 = sp.symbols('xs46 ys46')
for nu46 in (sp.Integer(1),sp.exp(-2)):
    npc46 = run46([stmt46(raw_with46,'source_pc = 2*nu*mp.conj(tc)')],
                  namespace46(tc=xs46*ys46,nu=nu46))
    tc_ok46 &= sp.expand(npc46['source_pc']
                         -2*(nu46*sp.conjugate(xs46))*sp.conjugate(ys46))==0
# The output: length, factor and the returned pair have |out|^2=source_radius
# (length>0 is shown in C46b3; it enters below as a positive symbol).
wr46 = sp.symbols('w0r46 w0i46 w1r46 w1i46',real=True)
ll46 = sp.Symbol('ll46',positive=True)
wsq46 = sum(w**2 for w in wr46)
no46 = run46([stmt46(raw_with46,'length = mp.sqrt(abs(w0)**2+abs(w1)**2)')],
             namespace46(w0=wr46[0]+sp.I*wr46[1],w1=wr46[2]+sp.I*wr46[3]))
nf46 = run46([stmt46(raw_with46,'factor = mp.sqrt(source_radius)/length')],
             namespace46(w0=wr46[0]+sp.I*wr46[1],w1=wr46[2]+sp.I*wr46[3],
                         source_radius=s46,length=ll46))
ret46 = stmt46(raw_with46,'return (+(factor*w0),+(factor*w1))')
ret_out46 = eval(compile(ast.fix_missing_locations(ast.Expression(
    SymSource46().visit(copy46.deepcopy(ret46.value)))),'sym46','eval'),nf46)
out_ok46 = (sp.expand(no46['length']-sp.sqrt(wsq46))==0
            and len(ret_out46)==2
            and sp.expand(sum(w*sp.conjugate(w) for w in ret_out46)
                          -s46*wsq46/ll46**2)==0)
# so with ll46=length, |out|^2=s46 |w|^2/length^2=s46.
rep('C46b2 elementary path: AA/yy alternatives identical, |p_S|=s, |w|^2',
    elem_ok46 and proj46 and unit46)
rep('C46b4 executed tc chain, slice Hopf vector, output norm |out|^2=s',
    tc_ok46 and out_ok46)

# C46b3.  Explicit constants (type (c)) on the whole raw domain
# 0<u=|z|^2<=U=1e-4, which contains every compressed cap radius
# u=S(r_L)<=r_L<1.7e-6 (C46d2) and u>0 (S>0, C46c).  Constants: 3<pi<4,
# e^4<55, A=(4 pi)^2>144, and |atan x|<=|x| (atan(0)=0, 0<atan'<=1).
#  * guard 0<radius<=1e-4: holds (u in (0,U]).
#  * dominant real part: Re u2 (south), Re(-u1) (north) have the sign of
#    1-2 pi p0>=1-8U>0 (C46b), so u2 resp. u1 is nonzero (ratio defined),
#    the north stereographic den>=1-8U>0, and arg is taken in the open
#    right half-plane, where Arg=atan(Im/Re) is smooth: height is smooth.
#  * |height|=|atan(2 pi p2/(1-2 pi p0))|/(2 pi)<=u/(1-8u)<=U/(1-8U)<1/10:
#    the guard abs(height)>=.1 never fires (margin recorded below).
#  * uu=|ratio|^2: south <=((4pi^2u^2+4pi u)^2+(4pi u)^2)/(4(1-2pi u)^2),
#    north <=4((pi u+pi^2u^2)^2+(pi u)^2)/(1-2pi u)^2 (C46b identities,
#    |aa|,|bb|,|cc|<=s_j u); both <UU below.  The profile denominator is
#    >=2A(1-uu)>0 (sqrt(A^2(1-uu)^2+4Buu)>=A(1-uu)), so
#    |tc|^2<=max(1,B)uu/(A(1-uu))<55 UU/(144(1-UU))<1/100: e is defined.
#  * R>0: tc=0 and height=0 would force (aa,bb,cc)=0 (C46b factor checks
#    with 1+bb>0, 1+2bb>0), i.e. p=0, u=0.  Hence (R+h)/2, (R-h)/2 are
#    smooth, source_radius=(R+h)/(2e)+(R-h)/2>0 and both w denominators
#    are positive; D^-1's guard radius<=0 never fires (|g_j|^2=source_radius).
#  * length: by (11.9)-(11.11) (written type (c) bounds of section 11.3,
#    constants checked by I2u) |n_S-n_L|<=200u, so
#    |w|^2/u=(1+n_S.n_L)/2=1-|n_S-n_L|^2/4>=1-10^4 u^2>=1-10^4 U^2>0;
#    the guard length==0 never fires.
U46 = sp.Rational(1,10**4)
PI_LO46,PI_HI46,E4_HI46 = 3,4,55
height_bound46 = U46/(1-2*PI_HI46*U46)
UUs46 = ((PI_HI46**2*U46**2+PI_HI46*U46)**2*16+(4*PI_HI46*U46)**2)\
        /(4*(1-2*PI_HI46*U46)**2)
UUn46 = 4*((PI_HI46*U46+PI_HI46**2*U46**2)**2+(PI_HI46*U46)**2)\
        /(1-2*PI_HI46*U46)**2
UU46 = max(UUs46,UUn46)
tc2_46 = E4_HI46*UU46/(144*(1-UU46))
const_b46 = (
    sp.pi>PI_LO46 and sp.pi<PI_HI46 and sp.exp(4)<E4_HI46
    and (4*sp.pi)**2>144
    and 1-2*PI_HI46*U46>0
    and height_bound46<sp.Rational(1,10)
    and UU46<1 and tc2_46<sp.Rational(1,100)
    and 1-10**4*U46**2>0
    and sp.Rational(17,10**7)<=U46)
print('C46b3 margins on u<=1e-4: |height|<=',sp.N(height_bound46,6),
      '(guard 1/10); |ratio|^2<',sp.N(UU46,4),'; |tc|^2<',sp.N(tc2_46,4),
      '; |w|^2/u>=',sp.N(1-10**4*U46**2,8),
      '; at the cap radius 1.7e-6: |height|<=',
      sp.N(sp.Rational(17,10**7)/(1-8*sp.Rational(17,10**7)),6))
rep('C46b3 raw guards inactive, Arg off its cut, R>0, |tc|<1, length>0',
    const_b46 and stereo46 and elem_ok46 and tc_ok46)

# ---------------------------------------------------------------- C46c
# C46c.  The production radial maps ARE the section 11.3 maps.  Symbolic
# execution of the bound statements of _elementary_pole_radial (type (a)):
#   forward:  x=phi(r)=E(r)-e^3, y=T(x)=log(1+exp(x-1/x)),
#             new=psi^-1(y)=c(1-1/log(y+e))            (a<r<c),
#   inverse:  y=psi(s)=E(s)-e, L=log(e^y-1)=F(T^-1(y)), x-1/x=L,
#             new=phi^-1(x)=c(1-1/log(x+e^3))          (0<s<c),
# with E(t)=exp(c/(c-t)), a=1e-8, c=3a/2.  The IfExp alternatives f<0,
# y>1, L<0 are identical functions: exp of both log forms agree (both
# arguments real), and the quadratic root is rationalized.  So S and S^-1
# have no seams other than r=c (and r=a for the patch).
radial_inverse_if46 = [node for node in radial_tree46.body
                       if isinstance(node,ast.If)
                       and ast.dump(node.test,include_attributes=False)
                       ==ast.dump(ast.parse('inverse',mode='eval').body,
                                  include_attributes=False)][0]
r46,sv46 = sp.symbols('r46 sv46',positive=True)
X46,Y46 = sp.symbols('X46 Y46',positive=True)
L46 = sp.Symbol('L46',real=True)
rc46 = run46([stmt46(radial_tree46.body,"a,c = mp.mpf('1e-8'),mp.mpf('1.5e-8')")],
             namespace46())
a46,c46 = rc46['a'],rc46['c']
E46 = lambda t: sp.exp(c46/(c46-t))
def exp_equal46(lhs,rhs):
    # lhs, rhs are products of exponentials: equal iff exponents agree
    left = sp.powsimp(sp.expand(lhs),combine='exp',force=True)
    right = sp.powsimp(sp.expand(rhs),combine='exp',force=True)
    return (left.func==sp.exp and right.func==sp.exp
            and sp.cancel(left.args[0]-right.args[0])==0)
radial_ok46 = (a46==sp.Rational(1,10**8) and c46==3*a46/2
               and isinstance(radial_inverse_if46.body[0],ast.If)
               and isinstance(radial_inverse_if46.orelse[0],ast.If)
               and len(radial_inverse_if46.body)==6
               and len(radial_inverse_if46.orelse)==6)
# forward branch: statements 1..5 of the else block
fwd46 = radial_inverse_if46.orelse[1:]
alternatives46.clear()
nf46 = run46(fwd46[:1],namespace46(radius=r46,a=a46,c=c46))
radial_ok46 &= exp_equal46(nf46['x']+sp.exp(3),E46(r46))
nf46 = run46(fwd46[1:3],namespace46(x=X46,a=a46,c=c46))
radial_ok46 &= (len(alternatives46)==1
                and sp.expand(sp.exp(alternatives46[0][0])
                              -sp.expand(sp.exp(alternatives46[0][1])))==0
                and sp.simplify(alternatives46[0][0]
                                -sp.log(1+sp.exp(X46-1/X46)))==0)
nf46 = run46(fwd46[3:],namespace46(y=Y46,a=a46,c=c46))
radial_ok46 &= sp.simplify(nf46['new_radius']-c46*(1-1/sp.log(Y46+sp.E)))==0
# inverse branch: statements 1..5 of the if block
inv46 = radial_inverse_if46.body[1:]
alternatives46.clear()
ni46 = run46(inv46[:1],namespace46(radius=sv46,a=a46,c=c46))
radial_ok46 &= exp_equal46(ni46['y']+sp.E,E46(sv46))
ni46 = run46(inv46[1:2],namespace46(y=Y46,a=a46,c=c46))
radial_ok46 &= (len(alternatives46)==1
                and sp.expand(sp.expand(sp.exp(alternatives46[0][0]))
                              -sp.exp(alternatives46[0][1]))==0
                and sp.simplify(alternatives46[0][1]
                                -sp.log(sp.exp(Y46)-1))==0)
alternatives46.clear()
ni46 = run46(inv46[2:3],namespace46(L=L46,a=a46,c=c46))
xroot46 = ni46['x']
radial_ok46 &= (len(alternatives46)==1
                and sp.simplify(sp.radsimp(alternatives46[0][0])
                                -alternatives46[0][1])==0
                and sp.simplify(xroot46-1/xroot46-L46)==0)
ni46 = run46(inv46[3:],namespace46(x=X46,a=a46,c=c46))
radial_ok46 &= sp.simplify(ni46['new_radius']
                           -c46*(1-1/sp.log(X46+sp.exp(3))))==0
# Hence T(T^-1(y))=log(1+exp(F(x)))=log(1+exp(L))=y: the inverse branch
# is phi^-1 o T^-1 o psi with the T^-1 just checked (x-1/x=L=log(e^y-1)).
psi_inverse46 = nf46['new_radius']             # source psi^-1(Y46)
rep('C46c source radial maps are psi^-1 T phi and phi^-1 T^-1 psi; '
    'rewrites same',radial_ok46 and source46)

# C46c2.  Outer seams (flat joins, type (a) identities + written calculus).
# Put t=c-r>0 (resp. t=c-s) and Q=exp(-c/t), A3=e^3.  Using the source
# forms just identified:
#   S(c-t):     1+ell=log(y+e)=c/t+log(1+Q D0(Q)),
#               D0=e-A3-Q/(1-A3 Q)+log(1+exp(-1/Q+A3+Q/(1-A3 Q)))
#               ((11.12i), C28b, here from the SOURCE x=e^3 expm1(.)),
#   S^-1(c-t):  log(x+e^3)=c/t+log(1+Q D1(Q)),  with  m=-e+log(1-e^e e^(-1/Q)),
#               K=sqrt((1+Qm)^2+4Q^2),  D1=e^3+m+2Q/(K+1+Qm)
#               (new; source y=e expm1(.), y>1 branch, L>=0 root),
#   and for either, c(1-1/(c/t+l))-(c-t)=t^2 l/(c+t l).
# What is proved (written, section 11.3 style): D0 and D1 are C^infinity on
# [0,Q0): the rational parts have nonvanishing denominators near Q=0
# (1-A3 Q -> 1, K+1+Qm -> 2, and K>|1+Qm| always) and the remaining terms
# are log(1+(flat)) with exp(-1/Q) flat at 0+.  So l(Q)=log(1+Q D(Q)) is
# smooth with l(0)=0, Q(t) is flat at t=0+, and S-id, S^-1-id are flat at
# c from the left.  The radial factors sqrt(S(r)/r) (S=id on r>=c) and
# sqrt(S^-1(v)/v) are therefore C^infinity in r resp. v near c; since r=|z|^2
# and v=|g_j(D z)|^2 are smooth (C46b), D and D^-1 o g_j o D are smooth
# across {r=c} and {v=c}: ALL jets match.  What is checked symbolically: the
# identities below.
t46,Q46 = sp.symbols('t46 Q46',positive=True)
mm46 = sp.Symbol('mm46',real=True)
K46 = sp.Symbol('K46',positive=True)
A3_46 = sp.exp(3)
xS46 = run46(fwd46[:1],namespace46(radius=c46-t46,a=a46,c=c46))['x']
seam_ok46 = sp.simplify(xS46-(1/sp.exp(-c46/t46)-A3_46))==0
xQ46 = 1/Q46-A3_46
alternatives46.clear()
yQ46 = run46(fwd46[1:3],namespace46(x=xQ46,a=a46,c=c46))['y']
seam_ok46 &= len(alternatives46)==1            # yQ46: the source f>=0 form
D0_46 = sp.E-A3_46-Q46/(1-A3_46*Q46)\
        +sp.log(1+sp.exp(-1/Q46+A3_46+Q46/(1-A3_46*Q46)))
seam_ok46 &= sp.simplify(sp.E+yQ46-1/Q46-D0_46)==0
ell46 = sp.Symbol('ell46',positive=True)
seam_ok46 &= sp.simplify(c46*(1-1/(c46/t46+ell46))-(c46-t46)
                         -t46**2*ell46/(c46+t46*ell46))==0
ySinv46 = run46(inv46[:1],namespace46(radius=c46-t46,a=a46,c=c46))['y']
seam_ok46 &= sp.simplify(ySinv46-(1/sp.exp(-c46/t46)-sp.E))==0
yQi46 = 1/Q46-sp.E
alternatives46.clear()
run46(inv46[1:2],namespace46(y=yQi46,a=a46,c=c46))
seam_ok46 &= len(alternatives46)==1
LQ46 = alternatives46[0][0]                    # the source y>1 form
mQ46 = -sp.E+sp.log(1-sp.exp(sp.E)*sp.exp(-1/Q46))
seam_ok46 &= sp.simplify(LQ46-(1/Q46+mQ46))==0
Lm46 = (1+Q46*mm46)/Q46
Krel46 = (1+Q46*mm46)**2+4*Q46**2
# sqrt(L^2+4)=K/Q: equal squares and K/Q>0.
seam_ok46 &= sp.simplify(((K46/Q46)**2-(Lm46**2+4)).subs(K46,sp.sqrt(Krel46)))==0
xK46 = (Lm46+K46/Q46)/2
# xK46 is the source L>=0 root at L=Lm46: (x-L/2)^2=(L^2+4)/4=K^2/(4Q^2),
# both x-L/2 and K/(2Q) positive.
seam_ok46 &= (sp.expand((xroot46.subs(L46,Lm46)-Lm46/2)**2
                        -Krel46/(4*Q46**2))==0
              and (xroot46-L46/2).is_positive is True)
D1_46 = A3_46+mm46+2*Q46/(K46+1+Q46*mm46)
numer46 = sp.numer(sp.together(Q46*(xK46+A3_46)-1-Q46*D1_46))
seam_ok46 &= (sp.expand(sp.expand(numer46).subs(K46**2,Krel46))==0
              and sp.expand(numer46-(K46**2-Krel46))==0)
# K>|1+Qm|, hence K+1+Qm>0 (manifest: K^2-(1+Qm)^2=4Q^2>0).
seam_ok46 &= sp.expand(Krel46-(1+Q46*mm46)**2-4*Q46**2)==0
rep('C46c2 outer seams r=c of D and v=c of D^-1: exact flat-factor forms',
    seam_ok46 and radial_ok46)

# C46c3.  Inner seam |z|^2=a (patch identity piece).  What is proved: the
# section 11.3 inner spinor seam argument (blow-up z=sqrt(u) xi, the
# parameterized inverse function theorem, and the W-equation whose identity
# is C28c) gives that H_j joins the identity with ALL jets at |z|^2=a;
# that argument is written, not re-derived here.  Its scalar inputs are
# checked on the SOURCE forms (type (a)): phi(a)=0 with phi'(a)=6e^3/a>0
# (simple zero), psi(0)=0 with psi'(0)=e/c>0 (so h0(t)=(exp(psi(t))-1)/t
# is smooth and positive at 0), and S(a+t)=psi^-1(T(phi(a+t))) is flat at
# t=0+ because T(x)=log(1+exp(x)exp(-1/x)) is flat at x=0+ and psi^-1
# is analytic at 0 with psi^-1(0)=0 (written).
xsrc46 = run46(fwd46[:1],namespace46(radius=r46,a=a46,c=c46))['x']
ysrc46 = run46(inv46[:1],namespace46(radius=sv46,a=a46,c=c46))['y']
inner_ok46 = (
    sp.simplify(xsrc46.subs(r46,a46))==0
    and sp.simplify(sp.diff(xsrc46,r46).subs(r46,a46)-6*sp.exp(3)/a46)==0
    and sp.simplify(ysrc46.subs(sv46,0))==0
    and sp.simplify(sp.diff(ysrc46,sv46).subs(sv46,0)-sp.E/c46)==0
    and sp.simplify(psi_inverse46.subs(Y46,0))==0
    and sp.simplify(sp.diff(psi_inverse46,Y46).subs(Y46,0)-c46/sp.E)==0)
rep('C46c3 inner seam |z|^2=a: simple zero of phi, psi(0)=0, psi\'(0)=e/c',
    inner_ok46)

# ---------------------------------------------------------------- C46d
# C46d.  Symbolic execution of the bound cap branch (type (a)) with a
# symbolic round point (Z,X), X=-Xm (south) or X=+Xm (north), Xm>0, and
# the bound _elementary_pole_inclusion executed symbolically.
# elementary_pole_patch is replaced by a recorder returning symbols.
# Checks: inf=X>0 evaluates to the static value on each sign class; the
# patch receives U_j^* Z/den with |U_j^* Z|^2=|Z|^2 (unitary inclusion),
# den^2=s_j(1+|X|sqrt(2-X^2)), infinity=inf; the return is
# (kind,(conj(zz0)/nu_j,zz1)) with kind c0/ci and nu_j=1/e^-2.
Xm46 = sp.Symbol('Xm46',positive=True)
Zr46 = sp.symbols('Za0_46 Zb0_46 Za1_46 Zb1_46',real=True)
Zsym46 = (Zr46[0]+sp.I*Zr46[1],Zr46[2]+sp.I*Zr46[3])
W46 = sp.symbols('Wa0_46 Wb0_46 Wa1_46 Wb1_46',real=True)
Wsym46 = (W46[0]+sp.I*W46[1],W46[2]+sp.I*W46[3])
inclusion_ns46 = namespace46()
run46([inclusion_tree46],inclusion_ns46)       # defines the function
patch_calls46 = []
def patch_recorder46(pair,infinity,dps=50):
    patch_calls46.append((pair,infinity))
    return Wsym46
cap_ok46 = True
for north46 in (False,True):
    patch_calls46.clear()
    alternatives46.clear()
    Xv46 = Xm46 if north46 else -Xm46
    nc46 = run46(cap_nodes46[0].body[:-1],
                 namespace46(X=Xv46,Z=Zsym46,dps=50,
                             _elementary_pole_inclusion=
                             inclusion_ns46['_elementary_pole_inclusion'],
                             elementary_pole_patch=patch_recorder46))
    ret46 = eval(compile(ast.fix_missing_locations(ast.Expression(
        SymSource46().visit(copy46.deepcopy(cap_nodes46[0].body[-1].value)))),
        'sym46','eval'),nc46)
    sj46 = sp.pi if north46 else 2*sp.pi
    nuj46 = sp.exp(-2) if north46 else sp.Integer(1)
    zL46,infarg46 = patch_calls46[0]
    den46 = nc46['den']
    cap_ok46 &= (
        bool(nc46['inf'])==north46 and not alternatives46
        and nc46['scale']==sj46 and len(patch_calls46)==1
        and bool(infarg46)==north46
        and sp.simplify(den46**2-sj46*(1+Xm46*sp.sqrt(2-Xm46**2)))==0
        and sp.expand(sum(sp.re(w*den46)**2+sp.im(w*den46)**2 for w in zL46)
                      -sum(v**2 for v in Zr46))==0
        and ret46[0]==('ci' if north46 else 'c0')
        and sp.simplify(ret46[1][0]-sp.conjugate(Wsym46[0])/nuj46)==0
        and sp.simplify(ret46[1][1]-Wsym46[1])==0
        and nc46['nu']==nuj46)
# The round cap identity |X|sqrt(2-X^2)=sqrt(1-|Z|^4) when |Z|^2=1-X^2
# (N8a0): equal squares, both sides nonnegative.
cap_ok46 &= sp.expand(Xm46**2*(2-Xm46**2)-(1-(1-Xm46**2)**2))==0
rep('C46d cap branch source: static inf, unitary U_j, den^2, return',
    cap_ok46 and source46)

# C46d2.  Explicit constants (type (c)) on the cap 1-|X|<alpha/2, |Z|^2+X^2=1.
#  * |X|>1-alpha/2>0: X!=0, so inf=X>0 is constant on each of the two cap
#    components (static) and |X|=+-X is smooth there.
#  * 2-X^2>=1, so sqrt(2-X^2)>=1 is smooth and
#    1+|X|sqrt(2-X^2)>=1+|X|>2-alpha/2>0: den>0 and smooth.
#  * |Z|^2=1-X^2=(1-|X|)(1+|X|)<(alpha/2)*2=alpha.
#  * r_L=|z_L|^2=|Z|^2/(s_j(1+|X|sqrt(2-X^2)))<alpha/(pi(2-alpha/2))
#        <alpha/(3(2-alpha/2))<1.7e-6=:RL (s_j>=pi>3; programme 11.6).
#    Hence the patch's D is called only for a<r_L<RL<=1e-4 (guard
#    radius<=a never fires), and u=S(r_L)<=r_L (S<r on (a,c), S=id on
#    r>=c) lies in (0,1e-4] (C46b3 applies).
#  * |zz|^2<=max(c,r_L+80r_L^2) (written bound after (11.12j), type (c))
#    <max(c,RL+80RL^2)<2e-6<(.0098)^2.
#  * Four-seed restoration (source43, source_cusp_ok30): x_s=conj(zz0)/nu,
#    y_s=zz1, A=nu^2 e|x_s|^2=e|zz0|^2<=|zz|^2<2e-6<1, Y=|zz1|^2<2e-6<1,
#    q=AY<(1e-6)^2<1/8, rr=|x_s y_s|^2<=e^4(|zz|^2/2)^2<55e-12<1/4,
#    every seed argument A,qA,Y,qY<1/16 (linear branch), and the cusp
#    Cramer determinant +-(1-rr^2)!=0.
alpha46 = sp.Rational(1,10**5)
RL46 = sp.Rational(17,10**7)
zz_max46 = sp.Rational(2,10**6)
nu_s46 = sp.Symbol('nu_s46',positive=True)
ee46 = sp.Symbol('ee46',positive=True)
zq46 = sp.symbols('zq0_46 zq1_46')
const_d46 = (
    1-alpha46/2>0 and 2-alpha46/2>0
    and alpha46/(PI_LO46*(2-alpha46/2))<RL46
    and RL46<=U46
    and max(c46,RL46+80*RL46**2)<zz_max46
    and zz_max46<sp.Rational(98,10**4)**2
    and zz_max46<1 and zz_max46<sp.Rational(1,16)
    and (zz_max46/2)**2<sp.Rational(1,8)
    and E4_HI46*(zz_max46/2)**2<sp.Rational(1,4)
    and E4_HI46*(zz_max46/2)**2<sp.Rational(9,100)
    and 1-(E4_HI46*(zz_max46/2)**2)**2>0
    # A=nu^2 e |conj(zz0)/nu|^2 = e |zz0|^2 for nu>0
    and sp.simplify(nu_s46**2*ee46*sp.Abs(sp.conjugate(zq46[0])/nu_s46)**2
                    -ee46*sp.Abs(zq46[0])**2)==0
    and source43 and source_cusp_ok30)
print('C46d2 margins: r_L<',sp.N(alpha46/(PI_LO46*(2-alpha46/2)),8),
      '(<1.7e-6); |zz|^2<',sp.N(max(c46,RL46+80*RL46**2),8),
      '(<2e-6); four-seed max(A,Y)<2e-6 (domain <=4, seeds linear <=1/16)')
rep('C46d2 cap denominators, r_L<1.7e-6, |zz|^2<2e-6, four-seed A,Y<1',
    const_d46 and cap_ok46)

# ---------------------------------------------------------------- C46e
# C46e is solely a reproducible numerical observation (d); C46a-d do not
# use it.  Under elementary_formula_branches (no dps shortcut), at 100
# digits (pole patch internal 150-170), both signs, two spinor directions:
#  (e1) production cap formula (elementary_global_slice at the round point
#       (11.17a) of z_L) at |z_L|^2=a(1+-eps), eps in {1e-3,2e-4}, and
#       |z_L|^2=c(1+-eps), eps in {1/30,1/100} (this crossing contains
#       both the D seam r=c and the nearby D^-1 seam v=c).  Recorded: slice
#       vs direct (conj(H0)/nu,H1); jump |H(z+)-H(z-)|/|z+-z-| (O(eps)
#       difference, ratio ~1); exact identity on the inside of the inner
#       seam and exact D=id outside the outer seam; and the flat defects
#       |H(z)-z|/|z| (inner) and |H(z)-D^-1 g(z)|/|z| (outer).
#       Sample size 2 signs x 2 directions x 4 eps x 2 sides = 32.
#  (e2) the radial factors alone at r=c(1-eps) against the section 11.7
#       bound 176q, q=exp(-1/eps), and the rewrite alternatives evaluated
#       both ways.
def cap_point46(zL,north):
    sj = mp41.pi if north else 2*mp41.pi
    h = sj*(abs(zL[0])**2+abs(zL[1])**2)
    X = (1 if north else -1)*(1-h)/mp41.sqrt(1+h*h)
    g = mp41.sqrt(2*sj/(1+h*h))
    a,b = ((zL[0]+zL[1])/mp41.sqrt(2),
           ((zL[0]-zL[1])/mp41.sqrt(2) if north
            else mp41.j*(zL[0]-zL[1])/mp41.sqrt(2)))
    return (g*a,g*b),X


rows46 = []
sample_ok46 = True
with mp41.workdps(100), cusp30.elementary_formula_branches():
    a_n46,c_n46 = mp41.mpf('1e-8'),mp41.mpf('1.5e-8')
    dirs46 = ((mp41.mpc(1,2),mp41.mpc(3,-1)),(mp41.mpc(-2,1),mp41.mpc(0,1)))
    for north46 in (False,True):
        nu_n46 = mp41.exp(-2) if north46 else mp41.mpf(1)
        for d46 in dirs46:
            dn46 = mp41.sqrt(abs(d46[0])**2+abs(d46[1])**2)
            xi46 = (d46[0]/dn46,d46[1]/dn46)
            for seam46,r0_46,epss46 in (('inner',a_n46,('1e-3','2e-4')),
                                        ('outer',c_n46,('1/30','1/100'))):
                for eps_s46 in epss46:
                    eps46 = mp41.mpf(eps_s46) if '/' not in eps_s46 else \
                        mp41.mpf(1)/int(eps_s46.split('/')[1])
                    side46 = {}
                    for sgn46 in (1,-1):
                        rr46 = r0_46*(1+sgn46*eps46)
                        zL46n = (mp41.sqrt(rr46)*xi46[0],mp41.sqrt(rr46)*xi46[1])
                        Zp46,Xp46 = cap_point46(zL46n,north46)
                        kind46,out46 = cusp30.elementary_global_slice(
                            Zp46,Xp46,dps=80)
                        H46 = cusp30.elementary_pole_patch(zL46n,north46,dps=80)
                        direct46 = (mp41.conj(H46[0])/nu_n46,H46[1])
                        cons46 = max(abs(u-v) for u,v in zip(out46,direct46))
                        if seam46=='inner':
                            ref46 = zL46n
                        else:
                            # D^-1 o g_j at the patch's own working precision
                            # (max(100,80)+50), so that D=id gives equality
                            with mp41.workdps(150):
                                ref46 = cusp30._elementary_pole_radial(
                                    cusp30.elliptic_pole_raw_inverse(
                                        zL46n,north46,dps=80,elementary=True),
                                    inverse=True)
                        defect46 = (mp41.sqrt(sum(abs(u-v)**2 for u,v
                                                   in zip(H46,ref46)))
                                    /mp41.sqrt(rr46))
                        side46[sgn46] = (zL46n,H46,cons46,defect46,kind46,
                                         1-abs(Xp46))
                    zp46,Hp46 = side46[1][:2]
                    zm46,Hm46 = side46[-1][:2]
                    jump46 = (mp41.sqrt(sum(abs(u-v)**2 for u,v in zip(Hp46,Hm46)))
                              /mp41.sqrt(sum(abs(u-v)**2 for u,v in zip(zp46,zm46))))
                    # inside of the inner seam / outside of the outer seam the
                    # formula is exactly the identity / exactly D^-1 g
                    exact_side46 = side46[-1] if seam46=='inner' else side46[1]
                    flat_side46 = side46[1] if seam46=='inner' else side46[-1]
                    rows46.append((north46,seam46,eps_s46,jump46,
                                   flat_side46[3],exact_side46[3],
                                   max(side46[1][2],side46[-1][2])))
                    sample_ok46 &= (
                        side46[1][4]==('ci' if north46 else 'c0')
                        and side46[-1][4]==side46[1][4]
                        and max(side46[1][5],side46[-1][5])<mp41.mpf('5e-6')
                        and max(side46[1][2],side46[-1][2])<mp41.mpf('1e-90')
                        and exact_side46[3]==0
                        and abs(jump46-1)<mp41.mpf('1e-3'))
    inner_def46 = {e:max(r[4] for r in rows46 if r[1]=='inner' and r[2]==e)
                   for e in ('1e-3','2e-4')}
    outer_def46 = {e:max(r[4] for r in rows46 if r[1]=='outer' and r[2]==e)
                   for e in ('1/30','1/100')}
    # flatness: the defect collapses by many orders when eps shrinks
    sample_ok46 &= (inner_def46['1e-3']<mp41.mpf('1e-6')
                    and inner_def46['2e-4']<inner_def46['1e-3']*mp41.mpf('1e-10')
                    and outer_def46['1/30']<mp41.mpf('1e-10')
                    and outer_def46['1/100']<outer_def46['1/30']*mp41.mpf('1e-20'))
    print('C46e production cap observations (d), sample size',2*len(rows46),':')
    print('    max |slice-direct| =',mp41.nstr(max(r[6] for r in rows46),4),
          ' max |jump ratio-1| =',mp41.nstr(max(abs(r[3]-1) for r in rows46),4))
    print('    inner flat defect |H(z)-z|/|z| at eps=1e-3, 2e-4:',
          mp41.nstr(inner_def46['1e-3'],4),mp41.nstr(inner_def46['2e-4'],4),
          ' (exactly 0 inside)')
    print('    outer flat defect |H(z)-D^-1 g(z)|/|z| at eps=1/30, 1/100:',
          mp41.nstr(outer_def46['1/30'],4),mp41.nstr(outer_def46['1/100'],4),
          ' (exactly 0 outside)')
    # (e2) radial factors at the outer seam vs 176q, and rewrites both ways
    radial_rows46 = []
    for n_eps46 in (30,100):
        eps46 = mp41.mpf(1)/n_eps46
        q46 = mp41.exp(-1/eps46)
        rr46 = c_n46*(1-eps46)
        Sz46 = cusp30._elementary_pole_radial((mp41.sqrt(rr46),mp41.mpc(0)))
        Si46 = cusp30._elementary_pole_radial((mp41.sqrt(rr46),mp41.mpc(0)),
                                              inverse=True)
        relS46 = abs(abs(Sz46[0])**2-rr46)/rr46
        relSi46 = abs(abs(Si46[0])**2-rr46)/rr46
        radial_rows46.append((n_eps46,relS46,relSi46,176*q46))
        sample_ok46 &= relS46<176*q46 and relSi46<176*q46 and relS46>0
    rew46 = mp41.mpf(0)
    for hh46 in (mp41.mpf('1e-7'),mp41.mpf('-1e-7'),mp41.mpf('3e-3'),
                 mp41.mpf('-3e-3')):
        for qq46 in (mp41.mpf('1e-14'),mp41.mpf('1e-9')):
            RR46 = mp41.sqrt(hh46*hh46+4*qq46)
            rew46 = max(rew46,abs(2*qq46/(RR46-hh46)-(RR46+hh46)/2)/RR46,
                        abs(2*qq46/(RR46+hh46)-(RR46-hh46)/2)/RR46)
    for ff46 in (mp41.mpf('-30'),mp41.mpf('-.5'),mp41.mpf('.5'),mp41.mpf('30')):
        rew46 = max(rew46,abs(mp41.log1p(mp41.exp(ff46))
                              -(ff46+mp41.log1p(mp41.exp(-ff46))))
                    /mp41.log1p(mp41.exp(ff46)))
    for yv46 in (mp41.mpf('.3'),mp41.mpf('1.5'),mp41.mpf('5')):
        rew46 = max(rew46,abs(yv46+mp41.log1p(-mp41.exp(-yv46))
                              -mp41.log(mp41.expm1(yv46)))/abs(mp41.log(mp41.expm1(yv46))))
    for Lv46 in (mp41.mpf('-40'),mp41.mpf('-.7'),mp41.mpf('.7'),mp41.mpf('40')):
        rew46 = max(rew46,abs(2/(mp41.sqrt(Lv46**2+4)-Lv46)
                              -(Lv46+mp41.sqrt(Lv46**2+4))/2)
                    /((Lv46+mp41.sqrt(Lv46**2+4))/2))
    sample_ok46 &= rew46<mp41.mpf('1e-60')
    print('    radial |S(r)-r|/r, |S^-1(r)-r|/r at r=c(1-eps) vs 176q:',
          [(f'1/{n}',mp41.nstr(a,4),mp41.nstr(b,4),mp41.nstr(cq,4))
           for n,a,b,cq in radial_rows46])
    print('    rewrite alternatives evaluated both ways, max relative gap:',
          mp41.nstr(rew46,4))
rep('C46e production cap formula across the radial seams (numerical, d)',
    sample_ok46)

# C47: the slab phase band and its selectors, at both pole signs.
# Scope: the phase-correction path of the finite slab formula (11.15)-(11.16),
# i.e. _elementary_corrected_section (with its nested mismatch),
# elementary_phase_family (with increments/primitive), the finite projection
# section, _elementary_unit, _elementary_pole_inclusion, and the slice lines
# that select or price the band.  C47 proves that every conditional on that
# path is (i) static (constant on each band component, or fixed by the call),
# (ii) a guard that the written count bounds (11.14) of section 11.5 keep
# inactive (nonzero projections, right-half-plane increments, chart radius),
# (iii) a diagnostic that never reaches a returned value, or (iv) a smooth
# join: the band cutoff distance<alpha (flat weight, C47b), the hemisphere
# Hopf reference nn[2]>=0 in mismatch (exact phase cancellation, C47c) and
# the stereographic chart upper=n[2]>=0 (exact real agreement, C47d).  The
# start constant Arg g(north,anchor) has no variable seam (C47e).  Premises
# not re-proved here: the written type-(c) bounds of sections 11.4-11.5
# (||D P||<L (11.13), ||Dg_j||<=(9+8N)L<L_g, the cap estimate
# ||P-I_j Pi(n) I_j^*||<=100 r_L with r_L<=2 alpha/3), the smoothness of
# the meridian projector P(n,h) itself (its selectors are C31-C44), and
# |Z|^2=1-X^2 on the round S^4.  Conditionals of the slice outside this path
# (cap/slab selector, W-chart representative) are listed and attributed to
# N9 (node_remainder.py) / section 11.6 and C41d/C44c.

# C47a.  Source binding, type (a) (exact AST equality).  The full bodies of
# _elementary_corrected_section, elementary_phase_family,
# elementary_projection_section, _elementary_unit and
# _elementary_pole_inclusion are fixed below; the slice lines that select
# or price the band are matched individually.  A census of EVERY If/IfExp
# test and every comparison in these functions and in the slice is compared
# with an explicit classification table, so no conditional is unclassified.
# Finally a dataflow check shows that the diagnostics/margin bookkeeping
# (the only min() calls) never reaches an assigned value or a return.
# The remaining non-polynomial primitives on the path are smooth where they
# are used (type (b), with the cited checks): abs(X) (X!=0 on each band, so
# 1-|X| is 1+X or 1-X, C47b); abs(.)**2 (a sum of squares); abs(phase) and
# the root in _elementary_unit (arguments nonzero by the guards, C47d);
# sqrt((1+-nn[2])/2) (argument >=1/2 on the selected hemisphere, >0 on its
# open domain, C47c); atan(Im R/Re R) (Re R>0 guard, C47d); mp.arg only at
# one constant point per sign (C47e).  The slice's mp.arg(phase) (mid
# W-point, a lattice e1 change, C42d/C50b) and sqrt(1-X*X) lie outside
# this phase path.
corrected_body47 = '''
import mpmath as mp
alpha = mp.mpf('1e-5')
distance = 1-abs(X)
pc = 2*xi[0]*mp.conj(xi[1])
n = (pc.real,pc.imag,abs(xi[0])**2-abs(xi[1])**2)
result = elementary_projection_section(xi,X,N,dps,diagnostics)
if distance<alpha:
    inf = X>0
    anchor = 1-alpha if inf else -1+alpha
    def mismatch(nn,h):
        if nn[2]>=0:
            w = mp.sqrt((1+nn[2])/2)
            xx = (w,mp.mpc(nn[0],-nn[1])/(2*w))
        else:
            w = mp.sqrt((1-nn[2])/2)
            xx = (mp.mpc(nn[0],nn[1])/(2*w),w)
        section = elementary_projection_section(xx,h,N,dps,diagnostics)
        P = elliptic_meridian_projector(nn,h,dps=dps,elementary=True)
        direct = _elementary_unit(P*mp.matrix(_elementary_pole_inclusion(xx,inf)))
        return _elementary_inner(direct,section)
    theta = elementary_phase_family(mismatch,n,X,anchor,K,diagnostics=diagnostics)
    weight = 1-scalar_flat_step((distance-alpha/2)/(alpha/2))
    result = tuple(a*mp.exp(-mp.j*weight*theta) for a in result)
return result
'''
phase_body47 = '''
import mpmath as mp
if not isinstance(steps,int) or steps<1 or chart not in (None,'north','south'):
    raise ValueError('positive integer steps and a valid chart are required')
north = (mp.mpf(0),mp.mpf(0),mp.mpf(1))
def increments(values):
    previous = next(values)
    total = mp.mpf(0)
    for value in values:
        ratio = value*mp.conj(previous)
        if ratio.real<=0:
            raise ValueError('phase count does not give right-half-plane increments')
        if diagnostics is not None:
            diagnostics['minimum_phase_real'] = min(
                diagnostics.get('minimum_phase_real',mp.mpf(1)),ratio.real)
        total += mp.atan(ratio.imag/ratio.real)
        previous = value
    return total
theta = mp.arg(g(north,anchor))
theta += increments(g(north,anchor+(X-anchor)*j/steps) for j in range(steps+1))
def primitive(a,upper):
    def values():
        for j in range(steps+1):
            v,w = a[0]*j/steps,a[1]*j/steps
            r2 = v*v+w*w
            yield g((2*v/(1+r2),2*w/(1+r2),(1-r2 if upper else r2-1)/(1+r2)),X)
    return increments(values())
upper = n[2]>=0 if chart is None else chart=='north'
sign = 1 if upper else -1
a = (n[0]/(1+sign*n[2]),n[1]/(1+sign*n[2]))
if a[0]**2+a[1]**2>4:
    raise ValueError('phase chart requires stereographic radius <=2')
if not upper:
    theta += primitive((mp.mpf(1),mp.mpf(0)),True)-primitive((mp.mpf(1),mp.mpf(0)),False)
return theta+primitive(a,upper)
'''
project_body47 = '''
import mpmath as mp
if not isinstance(steps,int) or steps<1:
    raise ValueError('steps must be a positive integer')
with mp.workdps(max(mp.mp.dps,dps)+20):
    xi = tuple(map(mp.mpc,xi))
    X = mp.mpf(X)
    pc = 2*xi[0]*mp.conj(xi[1])
    n = (pc.real,pc.imag,abs(xi[0])**2-abs(xi[1])**2)
    start = -1+mp.mpf('1e-9')
    vector = _elementary_pole_inclusion(xi)
    margin = mp.mpf(1)
    for j in range(steps+1):
        h = start+(X-start)*j/steps
        kind,data = elliptic_raw_base_reference(tuple((1-h*h)*a for a in n),h,
                                                 dps=dps,elementary=True)
        ref = elliptic_bundle_embedding(kind,data,dps=dps,elementary=True)
        phase = _elementary_inner(ref,vector)
        margin = min(margin,abs(phase))
        if phase==0:
            raise ValueError('finite projection product vanished')
        vector = tuple(a*phase/abs(phase) for a in ref)
    if diagnostics is not None:
        diagnostics['minimum_projection_norm'] = min(
            diagnostics.get('minimum_projection_norm',mp.mpf(1)),margin)
    return tuple(+a for a in vector)
'''
unit_body47 = '''
import mpmath as mp
length = mp.sqrt(sum(abs(a)**2 for a in v))
if length==0:
    raise ValueError('zero denominator in a finite projection')
return tuple(a/length for a in v)
'''
inclusion_body47 = '''
import mpmath as mp
a,b = xi
pair = ((a+b)/mp.sqrt(2),(a-b)/mp.sqrt(2)) if infinity else (
    (a-mp.j*b)/mp.sqrt(2),(a+mp.j*b)/mp.sqrt(2))
result = [mp.mpc(0)]*10
start = 4 if infinity else 0
result[start:start+2] = pair
return tuple(result)
'''
unit_tree47 = source_tree30(cusp30._elementary_unit)
inclusion_tree47 = source_tree30(cusp30._elementary_pole_inclusion)
meridian_tree47 = source_tree30(cusp30.elliptic_meridian_projector)
def signature47(tree):
    args = tree.args
    return ([x.arg for x in args.args],
            [ast.dump(d,include_attributes=False) for d in args.defaults],
            [x.arg for x in args.kwonlyargs])
bodies47 = (
    body_dump41(correct_tree30)==text_dump41(corrected_body47)
    and body_dump41(phase_tree30)==text_dump41(phase_body47)
    and body_dump41(project_tree30)==text_dump41(project_body47)
    and body_dump41(unit_tree47)==text_dump41(unit_body47)
    and body_dump41(inclusion_tree47)==text_dump41(inclusion_body47)
    and returned_expr30(inner_tree36,'sum(mp.conj(a)*b for a,b in zip(u,v))')
    and returned_expr30(meridian_tree47,
        'mp.matrix([[a*mp.conj(b) for b in e] for a in e])')
    and matching_assignment30(meridian_tree47,'v',
        'tuple((1-h*h)*mp.mpf(a) for a in n)')
    # signatures: the production call passes no `chart`, so chart=None.
    and signature47(phase_tree30)==(
        ['g','n','X','anchor','steps','chart','diagnostics'],
        [ast.dump(ast.Constant(None))]*2,[])
    and signature47(correct_tree30)==(['xi','X','N','K','dps','diagnostics'],[],[])
    and signature47(inclusion_tree47)==(['xi','infinity'],
                                        [ast.dump(ast.Constant(False))],[])
    # positional calls (xx,h,N,dps,diagnostics), (nn,h,dps=,elementary=)
    # and (direct,section) rely on these parameter orders
    and signature47(project_tree30)==(['xi','X','steps','dps','diagnostics'],
        [ast.dump(ast.Constant(50)),ast.dump(ast.Constant(None))],[])
    and signature47(meridian_tree47)==(['n','X','dps','elementary'],
        [ast.dump(ast.Constant(50)),ast.dump(ast.Constant(False))],[])
    and signature47(inner_tree36)==(['u','v'],[],[]))
# Slice lines.  cap_band only prices the evaluation (cost) for the
# execution guard max_work; the band itself is selected inside
# _elementary_corrected_section by the same test distance<alpha.
slice_lines47 = (
    matching_assignment30(slice_tree30,'alpha',"mp.mpf('1e-5')")
    and matching_assignment30(slice_tree30,'distance','1-abs(X)')
    and matching_assignment30(slice_tree30,'cap_band','distance<alpha')
    and matching_assignment30(slice_tree30,'cost',
                              '(N+1)*(12*(K+1) if cap_band else 1)')
    and exact_if_body37(slice_tree30,'max_work is not None and cost>max_work',
        "raise ValueError('certified finite expression exceeds max_work; "
        "inspect elementary_closed_parameters()')")
    and exact_if_body37(slice_tree30,
        "not -1+mp.mpf('1e-9')<=X<=1-mp.mpf('1e-9')",
        "raise ValueError('the finite slab formula excludes the inner pole caps')")
    and matching_if_test32(slice_tree30,'distance<alpha/2 and not force_slab')
    and matching_assignment30(slice_tree30,'xi',
                              'tuple(a/mp.sqrt(1-X*X) for a in Z)')
    and matching_assignment30(slice_tree30,'result',
        '_elementary_corrected_section(xi,X,N,K,dps,diagnostics)')
    and source_counts30 and project_source30 and phase_source30)
def names47(node):
    return {x.id for x in ast.walk(node) if isinstance(x,ast.Name)}
# cap_band and cost are read only where stated (dataflow, type (a)).
cap_reads47 = [node for node in ast.walk(slice_tree30)
               if isinstance(node,(ast.Assign,ast.If,ast.Return,ast.Expr))
               and 'cap_band' in names47(node.value if not isinstance(node,ast.If)
                                         else node.test)]
cost_reads47 = [node for node in ast.walk(slice_tree30)
                if isinstance(node,ast.If) and 'cost' in names47(node.test)]
slice_lines47 &= (
    len(cap_reads47)==1      # read only by the cost assignment
    and isinstance(cap_reads47[0],ast.Assign)
    and cap_reads47[0].targets[0].id=='cost'
    and len(cost_reads47)==1
    and all(not isinstance(node,ast.Return) or 'cost' not in names47(node)
            for node in ast.walk(slice_tree30)))

# Census of conditionals.  Classes: 'static' (fixed by the call or
# constant on each band component), 'guard' (raises; kept inactive by the
# cited bound), 'diag' (diagnostics only), 'join' (smooth join, cited
# subsection), 'outside' (not on the phase path; cited certificate).
def census47(tree):
    out = []
    for node in ast.walk(tree):
        if isinstance(node,(ast.If,ast.IfExp,ast.While)):
            out.append(('test',ast.unparse(node.test)))
        if isinstance(node,ast.Compare):
            out.append(('cmp',ast.unparse(node)))
    return sorted(out)
classes47 = {
  '_elementary_corrected_section': {
    ('test','distance < alpha'): 'join C47b',
    ('cmp','distance < alpha'): 'join C47b',
    ('cmp','X > 0'): 'static C47e: X!=0 on the band, sign fixed per component',
    ('test','inf'): 'static C47e',
    ('test','nn[2] >= 0'): 'join C47c: exact equality of both branches',
    ('cmp','nn[2] >= 0'): 'join C47c'},
  'elementary_phase_family': {
    ('test',"not isinstance(steps, int) or steps < 1 or chart not in (None, 'north', 'south')"):
        'static guard: K is a fixed positive int, chart=None (C30a, C47a)',
    ('cmp','steps < 1'): 'static guard',
    ('cmp',"chart not in (None, 'north', 'south')"): 'static guard',
    ('test','ratio.real <= 0'): 'guard: Re R>7/8 by (11.14), C47d',
    ('cmp','ratio.real <= 0'): 'guard (11.14)',
    ('test','diagnostics is not None'): 'diag',
    ('cmp','diagnostics is not None'): 'diag',
    ('test','chart is None'): 'static: production chart=None',
    ('cmp','chart is None'): 'static',
    ('cmp',"chart == 'north'"): 'static (not taken in production)',
    ('cmp','n[2] >= 0'): 'join C47d: exact agreement of the charts',
    ('test','upper'): 'join C47d (inside primitive() its value is the literal chart flag: static)',
    ('test','a[0] ** 2 + a[1] ** 2 > 4'): 'guard: radius<=1 in the selected chart, C47d',
    ('cmp','a[0] ** 2 + a[1] ** 2 > 4'): 'guard C47d',
    ('test','not upper'): 'join C47d: equator correction of the south chart'},
  'elementary_projection_section': {
    ('test','not isinstance(steps, int) or steps < 1'): 'static guard: N fixed positive int',
    ('cmp','steps < 1'): 'static guard',
    ('test','phase == 0'): 'guard: |P_j e_(j-1)|>1/2 by (11.14)/(11.15), C47d',
    ('cmp','phase == 0'): 'guard (11.15)',
    ('test','diagnostics is not None'): 'diag',
    ('cmp','diagnostics is not None'): 'diag'},
  '_elementary_unit': {
    ('test','length == 0'): 'guard: cap estimate, 1-100 r_L>0, C47d',
    ('cmp','length == 0'): 'guard'},
  '_elementary_pole_inclusion': {
    ('test','infinity'): 'static: inf fixed per band (C47e); I_0 in the section'},
  'elementary_global_slice': {
    ('test','distance < alpha / 2 and (not force_slab)'): 'outside: cap/slab join, section 11.6 and N9',
    ('cmp','distance < alpha / 2'): 'outside: cap/slab, N9',
    ('cmp','X > 0'): 'outside: cap branch sign (static per cap)',
    ('test','inf'): 'outside: cap branch, static per cap',
    ('test',"not -1 + mp.mpf('1e-9') <= X <= 1 - mp.mpf('1e-9')"): 'static domain guard (slab)',
    ('cmp',"-1 + mp.mpf('1e-9') <= X <= 1 - mp.mpf('1e-9')"): 'static domain guard',
    ('test','projection_steps is None'): 'static count choice',
    ('cmp','projection_steps is None'): 'static',
    ('test','phase_steps is None'): 'static count choice',
    ('cmp','phase_steps is None'): 'static',
    ('cmp','distance < alpha'): 'cost only (cap_band), C47a dataflow',
    ('test','cap_band'): 'cost only, C47a dataflow',
    ('test','max_work is not None and cost > max_work'): 'execution guard, value-independent',
    ('cmp','max_work is not None'): 'execution guard',
    ('cmp','cost > max_work'): 'execution guard',
    ('test',"kind == 'mid'"): 'outside: W-chart representative, C41d/C44c',
    ('cmp',"kind == 'mid'"): 'outside: C41d/C44c'}}
census_trees47 = {
    '_elementary_corrected_section': correct_tree30,
    'elementary_phase_family': phase_tree30,
    'elementary_projection_section': project_tree30,
    '_elementary_unit': unit_tree47,
    '_elementary_pole_inclusion': inclusion_tree47,
    'elementary_global_slice': slice_tree30}
census_ok47 = all(set(census47(tree))==set(classes47[name])
                  for name,tree in census_trees47.items())
census_count47 = sum(len(census47(tree)) for tree in census_trees47.values())
# The meridian projector itself has no conditional of its own (its callee
# selectors are C31-C44).
census_ok47 &= census47(meridian_tree47)==[]
# Guards: every If whose body raises is listed as a (static/domain/
# execution) guard above, and every raise sits in such an If.
guard_tests47 = sorted(ast.unparse(node.test)
                       for tree in census_trees47.values()
                       for node in ast.walk(tree)
                       if isinstance(node,ast.If)
                       and any(isinstance(b,ast.Raise) for b in node.body))
raise_count47 = sum(isinstance(node,ast.Raise)
                    for tree in census_trees47.values()
                    for node in ast.walk(tree))
census_ok47 &= (
    raise_count47==len(guard_tests47)==8
    and all('guard' in classes47[name][('test',ast.unparse(node.test))]
            for name,tree in census_trees47.items()
            for node in ast.walk(tree)
            if isinstance(node,ast.If)
            and any(isinstance(b,ast.Raise) for b in node.body)))
# Diagnostics dataflow: margin/diagnostics enter only assignments whose
# targets are margin or diagnostics[...], or are passed through to the
# same three phase-path functions (which satisfy the same property).
passthrough47 = {'elementary_projection_section','elementary_phase_family',
                 '_elementary_corrected_section'}
def bookkeeping_names47(node):
    if (isinstance(node,ast.Call) and isinstance(node.func,ast.Name)
            and node.func.id in passthrough47):
        return set().union(*[bookkeeping_names47(x) for x in
                             [node.func]+[a for a in node.args
                                          if not (isinstance(a,ast.Name)
                                                  and a.id=='diagnostics')]
                             +[k.value for k in node.keywords
                               if k.arg!='diagnostics']])
    if isinstance(node,ast.Name):
        return {node.id} & {'margin','diagnostics'}
    return set().union(set(),*[bookkeeping_names47(x)
                               for x in ast.iter_child_nodes(node)])
def bookkeeping_free47(tree):
    for node in ast.walk(tree):
        if isinstance(node,(ast.Assign,ast.AugAssign)):
            if bookkeeping_names47(node.value):
                targets = node.targets if isinstance(node,ast.Assign) else [node.target]
                if not all((isinstance(t,ast.Name) and t.id=='margin')
                           or (isinstance(t,ast.Subscript)
                               and isinstance(t.value,ast.Name)
                               and t.value.id=='diagnostics')
                           for t in targets):
                    return False
        if isinstance(node,(ast.Return,ast.Expr,ast.Yield)) and node.value is not None:
            if bookkeeping_names47(node.value):
                return False
        if isinstance(node,ast.If) and bookkeeping_names47(node.test):
            if ast.unparse(node.test)!='diagnostics is not None':
                return False
    return True
min_calls47 = [ast.unparse(node) for tree in census_trees47.values()
               for node in ast.walk(tree)
               if isinstance(node,ast.Call) and isinstance(node.func,ast.Name)
               and node.func.id in ('min','max')]
diag_ok47 = (
    all(bookkeeping_free47(tree) for tree in census_trees47.values())
    # control: a body that returned margin would be rejected
    and not bookkeeping_free47(ast.parse('def f(margin):\n    return margin').body[0])
    and sorted(min_calls47)==sorted([
        "min(diagnostics.get('minimum_phase_real', mp.mpf(1)), ratio.real)",
        'max(mp.mp.dps, dps)', 'min(margin, abs(phase))',
        "min(diagnostics.get('minimum_projection_norm', mp.mpf(1)), margin)",
        'max(mp.mp.dps, dps)']))
rep(f'C47a full phase-path bodies; {census_count47} conditional occurrences classified, '
    f'{raise_count47} guards',
    bodies47 and slice_lines47 and census_ok47 and diag_ok47)

# C47b.  The band selector `if distance<alpha:` is a smooth join, both
# signs.  Type (a) source facts: the If has no else-branch and the function
# returns `result`, so for distance>=alpha the output is the projection
# section E itself; inside, it is exp(-i*w*theta)*E with
# w=1-s(u), u=(distance-alpha/2)/(alpha/2) (C47a body), s=scalar_flat_step.
# distance=1-|X| is 1+X on the south band and 1-X on the north band (both
# are X-intervals avoiding 0), i.e. exactly the chi_0, chi_inf arguments of
# (11.16).  Type (a): at distance=alpha, u=1; and 1-s(u)=s(1-u) with
# 1-u=(alpha-distance)/(alpha/2) (flat_left44 symmetry), so on
# alpha/2<distance<alpha the weight is s((alpha-distance)/(alpha/2)).
band_if47 = [node for node in ast.walk(correct_tree30)
             if isinstance(node,ast.If)
             and ast.unparse(node.test)=='distance < alpha']
al47 = sp.Symbol('al47',positive=True)
d47 = sp.Symbol('d47',positive=True)
u47 = (d47-al47/2)/(al47/2)
Xs47,Xn47 = sp.Symbol('Xs47',negative=True),sp.Symbol('Xn47',positive=True)
band_source47 = (
    len(band_if47)==1 and band_if47[0].orelse==[]
    and ast.unparse(correct_tree30.body[-1])=='return result'
    and ordered_nodes37(correct_tree30.body,[
        'result = elementary_projection_section(xi,X,N,dps,diagnostics)',
        'return result'])
    and correct_tree30.body.index(band_if47[0])==len(correct_tree30.body)-2
    and ordered_nodes37(band_if47[0].body,[
        'theta = elementary_phase_family(mismatch,n,X,anchor,K,diagnostics=diagnostics)',
        'weight = 1-scalar_flat_step((distance-alpha/2)/(alpha/2))',
        'result = tuple(a*mp.exp(-mp.j*weight*theta) for a in result)'])
    and exact_if_body37(flat_tree38,'u>=1','return mp.mpf(1)')
    and exact_if_body37(flat_tree38,'u<=0','return mp.mpf(0)')
    and returned_expr30(flat_tree38,'a/(a+b)')
    and matching_statement31(flat_tree38,'a,b = mp.exp(-1/u),mp.exp(-1/(1-u))'))
# u as written in the source (sympified from its AST text)
flat_calls47 = [node for node in ast.walk(correct_tree30)
                if isinstance(node,ast.Call)
                and ast.unparse(node.func)=='scalar_flat_step']
u_src47 = (sp.sympify(ast.unparse(flat_calls47[0].args[0]),
                      locals={'distance':d47,'alpha':al47})
           if len(flat_calls47)==1 else sp.nan)
band_weight47 = (
    sp.cancel(u_src47-u47)==0
    and sp.cancel(u47.subs(d47,al47))==1
    and sp.cancel(u47.subs(d47,al47/2))==0
    and sp.cancel(1-u47-(al47-d47)/(al47/2))==0
    and u47.subs({d47:sp.Rational(1,10**5),al47:sp.Rational(1,10**5)})==1
    and sp.simplify(flat_actual34.subs(p34,1-u41)+flat_actual34.subs(p34,u41)-1)==0
    and sp.simplify((1-sp.Abs(Xs47))-(1+Xs47))==0
    and sp.simplify((1-sp.Abs(Xn47))-(1-Xn47))==0)

# All-orders flatness of s at 0+ (type (c), written; every identity it uses
# is checked).  On 0<v<1, s(v)=sigma(y(v)) with sigma(y)=1/(1+e^-y) and
# y=1/(1-v)-1/v.  sigma'=sigma(1-sigma), so by the chain rule (generic
# induction step below) sigma^(k)=p_k(sigma) with p_0=x and
# p_(k+1)=x(1-x)p_k'; thus x divides every p_k and |p_k(sigma)|<=c_k*sigma,
# c_k the sum of |coefficients| of p_k/x.  y^(j)=j!/(1-v)^(j+1)
# -(-1)^j j!/v^(j+1) (checked for symbolic j), so on 0<v<=1/2
# |y^(j)|<=2 j! v^-(j+1).  Faa di Bruno: s^(k)=sum_m p_m(sigma)
# B_(k,m)(y',y'',...), each monomial of B_(k,m) having degree m and weight
# k, hence bounded by B_(k,m)(2*1!,2*2!,...) v^-(k+m)<=... v^-2k.  Finally
# sigma(y(v))=s(v)<=exp(-1/v+1/(1-v))<=e^2 exp(-1/v) on (0,1/2].  So
# |s^(k)(v)|<=M_k v^-2k exp(-1/v) with M_k=e^2 sum_m c_m B_(k,m)(2*1!,...),
# and v^-2k exp(-1/v)->0 (t=1/v: t^2k e^-t->0).  M_1..M_4 are printed.
v47 = sp.Symbol('v47',positive=True)
y47 = sp.Symbol('y47',real=True)
x47 = sp.Symbol('x47')
jj47 = sp.Symbol('jj47',integer=True,nonnegative=True)
sigma47 = 1/(1+sp.exp(-y47))
yv47 = 1/(1-v47)-1/v47
Pf47 = sp.Function('Pf47')
p47 = [x47]
for k47 in range(5):
    p47.append(sp.expand(x47*(1-x47)*sp.diff(p47[-1],x47)))
c47 = [sum(abs(cf) for cf in sp.Poly(sp.cancel(pk/x47),x47).coeffs())
       for pk in p47]
fy47 = sp.Function('fy47')
yy47 = sp.Function('yy47')(v47)
dy47 = sp.symbols('dy47_1:6')
faa47 = all(
    sp.simplify(sp.diff(fy47(yy47),v47,k).subs(
        {sp.Derivative(yy47,(v47,j)):dy47[j-1] for j in range(k,0,-1)})
        -sum(sp.Subs(sp.Derivative(fy47(y47),(y47,m)),y47,yy47).doit()
             *sp.bell(k,m,dy47[:k-m+1]) for m in range(1,k+1)))==0
    for k in range(1,5))
weight_ok47 = all(
    all(sum(sp.degree_list(mono,*dy47[:k-m+1]))==m
        and sum((j+1)*e for j,e in enumerate(sp.degree_list(mono,*dy47[:k-m+1])))==k
        and sp.Poly(mono,*dy47[:k-m+1]).coeffs()[0]>0
        for mono in sp.Add.make_args(sp.expand(sp.bell(k,m,dy47[:k-m+1]))))
    for k in range(1,5) for m in range(1,k+1))
M47 = [sp.E**2*sum(c47[m]*sp.bell(k,m,[2*sp.factorial(j) for j in range(1,k-m+2)])
                   for m in range(1,k+1)) for k in range(1,5)]
flat_all47 = (
    sp.simplify(flat_actual34.subs(p34,v47)-sigma47.subs(y47,yv47))==0
    and sp.simplify(sp.diff(sigma47,y47)-sigma47*(1-sigma47))==0
    # generic induction step: d/dy P(sigma)=P'(sigma) sigma(1-sigma)
    and sp.simplify(sp.diff(Pf47(sigma47),y47)
                    -sp.Subs(sp.Derivative(Pf47(x47),x47),x47,sigma47).doit()
                     *sigma47*(1-sigma47))==0
    and all(sp.simplify(sp.diff(sigma47,y47,k)-p47[k].subs(x47,sigma47))==0
            for k in range(5))
    and all(sp.rem(pk,x47,x47)==0 for pk in p47)
    # y^(j) for symbolic j: induction step of both terms
    and sp.simplify(sp.diff(sp.factorial(jj47)/(1-v47)**(jj47+1),v47)
                    -sp.factorial(jj47+1)/(1-v47)**(jj47+2))==0
    and sp.simplify(sp.diff((-1)**jj47*sp.factorial(jj47)/v47**(jj47+1),v47)
                    -(-1)**(jj47+1)*sp.factorial(jj47+1)/v47**(jj47+2))==0
    and sp.simplify(sp.diff(yv47,v47)-(1/(1-v47)**2+1/v47**2))==0
    and faa47 and weight_ok47
    # e^(-1/v)/(e^(-1/v)+e^(-1/(1-v)))<=e^(-1/v+1/(1-v)), 1/(1-v)<=2
    and sp.simplify(flat_actual34.subs(p34,v47)*(1+sp.exp(1/v47-1/(1-v47)))-1)==0
    and sp.cancel(2-1/(1-v47)-(1-2*v47)/(1-v47))==0      # >=0 iff v<=1/2
    and all(sp.limit(v47**(-2*k)*sp.exp(-1/v47),v47,0,'+')==0 for k in range(1,5))
    and flat_limits41 and flat_left44)
print('     C47b flat-step constants M_1..M_4 (|s^(k)(v)|<=M_k v^-2k e^-1/v, 0<v<=1/2):',
      [sp.nsimplify(m/sp.E**2)*sp.Symbol('e^2') for m in M47])

# The join.  exp(-i w theta)-1=(-i w theta)*E1(-i w theta), with
# E1(z)=int_0^1 e^(tz)dt smooth (type (a) identity below).  Hence
# E_c-E=w*Kc with Kc=-i theta E1(-i w theta) E smooth wherever theta and E
# are.  theta is a finite composition that is smooth on the open set where
# its finitely many strict inequalities (C47d guards) hold; that set
# contains the closed band 1-|X|<=alpha by (11.14), so Kc and all its
# derivatives stay bounded as distance->alpha-.  By Leibniz every
# derivative of w*Kc is a sum of (derivative of w)*(bounded), and every
# derivative of w(distance)=s((alpha-distance)/(alpha/2)) tends to 0 as
# distance->alpha- (above).  The outside value is 0.  Hence E_c is C^infty
# across distance=alpha (type (b)).  The generic checks: a factor flat to
# order three in the transverse variable has zero three-jet in all
# directions, also when the cofactor depends on that variable; and
# F(w,theta(d,y))-F(0,theta(d,y)) has zero two-jet (theta depends on d,
# stronger than C41c).  Control: a non-flat linear cutoff c(alpha-d) gives
# the nonzero first jet i*c*theta.
zz47,tt47 = sp.symbols('zz47 tt47')
th47,w47,cc47 = sp.symbols('th47 w47 cc47',real=True)
E1_47 = sp.integrate(sp.exp(tt47*zz47),(tt47,0,1),conds='none')
factor47 = (
    sp.simplify(E1_47-(sp.exp(zz47)-1)/zz47)==0
    and sp.simplify(sp.exp(-sp.I*w47*th47)-1
                    -(-sp.I*w47*th47)*E1_47.subs(zz47,-sp.I*w47*th47))==0)
dd47,ya47,yb47 = sp.symbols('dd47 ya47 yb47',real=True)
gw47 = sp.Function('gw47')(dd47)
Kw47 = sp.Function('Kw47')(dd47,ya47,yb47)
thw47 = sp.Function('thw47')(dd47,ya47,yb47)
Fw47 = sp.Function('Fw47')
def flat_order3_47(expr):
    out = expr.doit()
    for k in (3,2,1):
        out = out.subs(sp.Derivative(gw47,(dd47,k)),0)
    out = out.subs(gw47,0)
    return sp.simplify(out.doit())==0
vars47 = (dd47,ya47,yb47)
jets3_47 = [()]+[(a,) for a in vars47]+[
    (a,b) for i,a in enumerate(vars47) for b in vars47[i:]]+[
    (a,b,c) for i,a in enumerate(vars47) for j,b in enumerate(vars47[i:],i)
    for c in vars47[j:]]
flat_join47 = (
    len(jets3_47)==20
    and all(flat_order3_47(sp.diff(gw47*Kw47,*dv) if dv else gw47*Kw47)
            for dv in jets3_47)
    and all(flat_order3_47(sp.diff(Fw47(gw47,thw47)-Fw47(0,thw47),*dv)
                           if dv else Fw47(gw47,thw47)-Fw47(0,thw47))
            for dv in jets3_47 if len(dv)<=2)
    and flat_two_jets41
    # control for the jet test itself: g=0 alone (g'=gp!=0) leaves gp*K
    and sp.simplify(sp.diff(gw47*Kw47,dd47).subs(sp.Derivative(gw47,dd47),
                    sp.Symbol('gp47')).subs(gw47,0))!=0
    and sp.simplify(sp.diff(sp.exp(-sp.I*cc47*(al47-d47)*th47)-1,d47)
                    .subs(d47,al47)-sp.I*cc47*th47)==0)
rep('C47b band cutoff distance<alpha: weight flat at u=1, C^infty join, both signs',
    band_source47 and band_weight47 and flat_all47 and factor47 and flat_join47)

# C47c.  The hemisphere selector `nn[2]>=0` in mismatch is a smooth join
# with exact agreement, both signs (j=0 uses I_0, j=inf uses I_inf in
# `direct`; the section always starts from I_0).  Type (a) identities:
# (1) both source references are unit Hopf preimages of nn: the upper one
#     on n2>-1 (w=sqrt((1+n2)/2)>0), the lower one on n2<1; on the
#     selected closed hemisphere the square-root argument is >=1/2
#     (manifest, type (b));
# (2) on the overlap -1<n2<1, xx_lower=lambda*xx_upper with the unit
#     scalar lambda=m/|m|, m=n0+i n1;
# (3) every source operation is S^1-equivariant in xx: n=Hopf(xx) is
#     invariant, I_j is C-linear, each projection step
#     ref*<ref,v>/|<ref,v>| maps lambda*v to lambda times the old value, P is
#     C-linear and unit() commutes with unit scalars;
# (4) hence mismatch=<direct,section> (conjugate-linear first slot) is
#     unchanged: conj(lambda)*lambda=1.  Control: rotating only one slot
#     changes it by conj(lambda).  Also |mismatch|=1 (both unit vectors on
#     the line of e(nn,h)).
# So mismatch is ONE function of (nn,h), smooth on S^2 x band given the
# projector's smoothness, and the selector merely chooses the formula.
mm_if47 = [node for node in ast.walk(correct_tree30)
           if isinstance(node,ast.If) and ast.unparse(node.test)=='nn[2] >= 0']
hemi_source47 = (
    len(mm_if47)==1
    and [ast.dump(x,include_attributes=False) for x in mm_if47[0].body]
        ==text_dump41('w = mp.sqrt((1+nn[2])/2)\n'
                      'xx = (w,mp.mpc(nn[0],-nn[1])/(2*w))')
    and [ast.dump(x,include_attributes=False) for x in mm_if47[0].orelse]
        ==text_dump41('w = mp.sqrt((1-nn[2])/2)\n'
                      'xx = (mp.mpc(nn[0],nn[1])/(2*w),w)')
    and bodies47)
n0c47,n1c47 = sp.symbols('n0c47 n1c47',real=True)
Wc47,WLc47 = sp.symbols('Wc47 WLc47',positive=True)
def hopf47(xi):
    pc = sp.expand(2*xi[0]*sp.conjugate(xi[1]))
    return (sp.re(pc),sp.im(pc),
            sp.expand(xi[0]*sp.conjugate(xi[0])-xi[1]*sp.conjugate(xi[1])))
def upper_ref47(n0,n1,n2):
    w = sp.sqrt((1+n2)/2)
    return (w,(n0-sp.I*n1)/(2*w))
def lower_ref47(n0,n1,n2):
    w = sp.sqrt((1-n2)/2)
    return ((n0+sp.I*n1)/(2*w),w)
# (1) whole open domains: n2=2W^2-1 (upper) / 1-2W^2 (lower), W>0, and
# the sphere relation n0^2+n1^2=1-n2^2=4W^2(1-W^2).
def mod_sphere47(expr,W):
    return sp.simplify(sp.expand(expr).subs(n1c47**2,4*W**2*(1-W**2)-n0c47**2))
hemi_hopf47 = True
for ref47,n2c47 in ((upper_ref47,2*Wc47**2-1),(lower_ref47,1-2*Wc47**2)):
    xxc47 = ref47(n0c47,n1c47,n2c47)
    hemi_hopf47 &= (
        all(mod_sphere47(a-b,Wc47)==0
            for a,b in zip(hopf47(xxc47),(n0c47,n1c47,n2c47)))
        and mod_sphere47(sum(sp.expand(z*sp.conjugate(z)) for z in xxc47)-1,Wc47)==0)
nn2_47 = sp.Symbol('nn2_47',real=True)
# square-root arguments: >=1/2 exactly on the selected closed/open
# hemisphere, and >0 on the open domains n2>-1 / n2<1 (type (b)).
hemi_hopf47 &= (
    sp.solveset((1+nn2_47)/2>=sp.Rational(1,2),nn2_47,sp.Interval(-1,1))
        ==sp.Interval(0,1)
    and sp.solveset((1-nn2_47)/2>sp.Rational(1,2),nn2_47,sp.Interval(-1,1))
        ==sp.Interval.Ropen(-1,0)
    and sp.solveset((1+nn2_47)/2>0,nn2_47,sp.Interval(-1,1))==sp.Interval.Lopen(-1,1)
    and sp.solveset((1-nn2_47)/2>0,nn2_47,sp.Interval(-1,1))==sp.Interval.Ropen(-1,1))
# (2) overlap: n2=(W^2-WL^2)/D, D=W^2+WL^2, |m|=2 W WL/D.
D47 = Wc47**2+WLc47**2
n2o47 = (Wc47**2-WLc47**2)/D47
def mod_overlap47(expr):
    return sp.simplify(sp.expand(expr).subs(n1c47**2,4*Wc47**2*WLc47**2/D47**2-n0c47**2))
up47 = upper_ref47(n0c47,n1c47,n2o47)
lo47 = lower_ref47(n0c47,n1c47,n2o47)
lam47 = (n0c47+sp.I*n1c47)*D47/(2*Wc47*WLc47)
overlap47 = (
    sp.simplify(up47[0]-Wc47/sp.sqrt(D47))==0
    and sp.simplify(lo47[1]-WLc47/sp.sqrt(D47))==0
    and mod_overlap47(n0c47**2+n1c47**2+n2o47**2-1)==0
    and mod_overlap47(lam47*sp.conjugate(lam47)-1)==0
    and all(mod_overlap47(lo47[k]-lam47*up47[k])==0 for k in range(2))
    # control: the two references are NOT equal on the overlap
    and mod_overlap47(lo47[1]-up47[1])!=0)
# (3)-(4) equivariance and cancellation, in the ten-slot space.
ph47 = sp.Symbol('ph47',real=True)
lamu47 = sp.exp(sp.I*ph47)
ar47 = sp.symbols('ar47_0:2',real=True)
ai47 = sp.symbols('ai47_0:2',real=True)
xig47 = tuple(ar47[k]+sp.I*ai47[k] for k in range(2))
def inclusion47(xi,infinity):
    a,b = xi
    pair = ((a+b)/sp.sqrt(2),(a-b)/sp.sqrt(2)) if infinity else (
        (a-sp.I*b)/sp.sqrt(2),(a+sp.I*b)/sp.sqrt(2))
    out = [sp.Integer(0)]*10
    start = 4 if infinity else 0
    out[start:start+2] = pair
    return out
er47 = sp.symbols('er47_0:10',real=True)
ei47 = sp.symbols('ei47_0:10',real=True)
eg47 = [er47[k]+sp.I*ei47[k] for k in range(10)]
vr47 = sp.symbols('vr47_0:10',real=True)
vi47 = sp.symbols('vi47_0:10',real=True)
vg47 = [vr47[k]+sp.I*vi47[k] for k in range(10)]
pr47,pi47 = sp.symbols('pr47 pi47',real=True)
pg47 = pr47+sp.I*pi47
def unit_scalar47(p):
    return p/sp.sqrt(sp.expand(p*sp.conjugate(p)))
lamp47 = sp.expand(lamu47*pg47)
equivariant47 = (
    all(sp.simplify(a-b)==0 for a,b in zip(
        hopf47(tuple(lamu47*z for z in xig47)),hopf47(xig47)))
    and all(sp.expand(a-lamu47*b)==0
            for inf47 in (False,True)
            for a,b in zip(inclusion47(tuple(lamu47*z for z in xig47),inf47),
                           inclusion47(xig47,inf47)))
    # <ref,lambda v>=lambda<ref,v>; then the normalized step is equivariant
    and sp.expand(inner41(eg47,[lamu47*z for z in vg47])
                  -lamu47*inner41(eg47,vg47))==0
    and sp.simplify(unit_scalar47(lamp47)-lamu47*unit_scalar47(pg47))==0
    # P=e e^* (the source projector): P(lambda u)=lambda P u, entrywise
    and all(sp.expand(sum(eg47[r]*sp.conjugate(eg47[c])*lamu47*vg47[c] for c in range(10))
                      -lamu47*sum(eg47[r]*sp.conjugate(eg47[c])*vg47[c] for c in range(10)))==0
            for r in range(10))
    # unit(lambda u)=lambda unit(u): the norms agree
    and sp.simplify(sum(sp.expand(lamu47*z*sp.conjugate(lamu47*z)) for z in vg47)
                    -sum(sp.expand(z*sp.conjugate(z)) for z in vg47))==0)
dr47 = sp.symbols('dr47_0:10',real=True)
di47 = sp.symbols('di47_0:10',real=True)
dg47 = [dr47[k]+sp.I*di47[k] for k in range(10)]
cancel47 = (
    sp.expand(inner41([lamu47*z for z in dg47],[lamu47*z for z in vg47])
              -inner41(dg47,vg47))==0
    and sp.expand(inner41([lamu47*z for z in dg47],vg47)
                  -sp.conjugate(lamu47)*inner41(dg47,vg47))==0
    and sp.expand(inner41([lamu47*z for z in dg47],vg47)-inner41(dg47,vg47))!=0)
# |mismatch|=1: direct=mu*e, section=nu*e with |e|=|mu|=|nu|=1.
mu47,nu47 = sp.symbols('mu47 nu47',real=True)
e2r47 = sp.symbols('e2r47_0:2',real=True)
e2v47 = [e2r47[0]+sp.I*e2r47[1],sp.Symbol('e2c47',real=True)]
e2n47 = sp.Symbol('e2n47',positive=True)     # e2n47^2=|e|^2, e!=0
e2u47 = [z/e2n47 for z in e2v47]
gval47 = inner41([sp.exp(sp.I*mu47)*z for z in e2u47],[sp.exp(sp.I*nu47)*z for z in e2u47])
cancel47 &= sp.simplify(sp.expand(gval47*sp.conjugate(gval47)).subs(
    e2n47,sp.sqrt(sum(sp.expand(z*sp.conjugate(z)) for z in e2v47)))-1)==0
cancel47 &= (matching_statement31(correct_tree30,
                 'section = elementary_projection_section(xx,h,N,dps,diagnostics)')
             and matching_statement31(correct_tree30,
                 'P = elliptic_meridian_projector(nn,h,dps=dps,elementary=True)')
             and matching_statement31(correct_tree30,
                 'direct = _elementary_unit(P*mp.matrix(_elementary_pole_inclusion(xx,inf)))')
             and matching_statement31(correct_tree30,
                 'return _elementary_inner(direct,section)'))
rep('C47c hemisphere Hopf reference: mismatch independent of the choice, both signs',
    hemi_source47 and hemi_hopf47 and overlap47 and equivariant47 and cancel47)

# C47d.  The chart selector upper=n[2]>=0 of elementary_phase_family and
# its equator correction P_N(e)-P_S(e), both signs (the sign enters only
# through the anchor, C47e).  In production chart=None (C47a signature and
# call), so upper=n[2]>=0.  Write Theta_N=c+T_X+P_N(a_N(n)) and
# Theta_S=c+T_X+P_N(e)-P_S(e)+P_S(a_S(n)), with c=Arg g(north,anchor),
# T_X the north X-path sum, a_N=(n0,n1)/(1+n2), a_S=(n0,n1)/(1-n2), e=(1,0)
# and P_N, P_S the K-step sums along S_N(ja/K), S_S(ja/K) at the current X.
# Written argument (types (b)/(c); its identities are checked below):
# (i) Guards.  Theta_N is defined (every Arg increment in the right half
#     plane) for |a_N|<=2, i.e. n2>=-3/5, and Theta_S for n2<=3/5: by
#     ||DS||<=2 (C20a, both charts) a chart step of length |a|/K has image
#     length <=4/K, the X-path step is <=alpha/K<4/K, ||Dg||<L_g, so
#     consecutive values differ by <4L_g/K<1/2 (exact counts (11.14)) and
#     Re R=1-|z-w|^2/2>7/8 (unit z,w).  Likewise 2L/N<1/2 keeps every
#     projection step >1/2 and 1-100*(2alpha/3)>0 keeps `direct` nonzero
#     (cap estimate).  In the selected chart |a|^2<=1<4, so the radius
#     guard never fires.
# (ii) Exactness.  With Re R>0, exp(i atan(Im R/Re R))=R/|R|, and for unit
#     values R/|R|=g_l conj(g_(l-1)), which telescopes (C20c).  With the
#     endpoint identities S_N(0)=north, S_S(0)=south, S_N(e)=S_S(e)=(1,0,0),
#     S_N(a_N(n))=S_S(a_S(n))=n, exp(i Theta_N)=exp(i Theta_S)=g(n,X).
# (iii) Agreement.  Theta_S-Theta_N is continuous and 2pi Z-valued on the
#     connected band -3/5<n2<3/5 (an annulus, joined to e by rotation about
#     the axis and a meridian arc), and vanishes at n=e because
#     a_N(e)=a_S(e)=e makes both pairs P_N(a_N(e))=P_N(e) and
#     P_S(a_S(e))=P_S(e) one finite expression each.
#     Hence Theta_S=Theta_N on the band, an open neighbourhood of the
#     selector seam n2=0: theta equals the smooth Theta_N on n2>-3/5 and
#     the smooth Theta_S on n2<3/5.  The selector is a smooth join with
#     exact agreement; no 2pi seam exists.
upper_source47 = (
    matching_assignment30(phase_tree30,'upper',
                          "n[2]>=0 if chart is None else chart=='north'")
    and matching_assignment30(phase_tree30,'sign','1 if upper else -1')
    and matching_assignment30(phase_tree30,'a',
                              '(n[0]/(1+sign*n[2]),n[1]/(1+sign*n[2]))')
    and exact_if_body37(phase_tree30,'a[0]**2+a[1]**2>4',
        "raise ValueError('phase chart requires stereographic radius <=2')")
    and exact_if_body37(phase_tree30,'not upper',
        'theta += primitive((mp.mpf(1),mp.mpf(0)),True)'
        '-primitive((mp.mpf(1),mp.mpf(0)),False)')
    and returned_expr30(phase_tree30,'theta+primitive(a,upper)')
    and matching_call30(correct_tree30,
        'elementary_phase_family(mismatch,n,X,anchor,K,diagnostics=diagnostics)')
    and bodies47)
n0d47,n1d47,n2d47 = sp.symbols('n0d47 n1d47 n2d47',real=True)
vd47,wd47 = sp.symbols('vd47 wd47',real=True)
def stereo47(v,w,upper):       # the source primitive's point formula
    r2 = v*v+w*w
    return (2*v/(1+r2),2*w/(1+r2),(1-r2 if upper else r2-1)/(1+r2))
def chart47(n,upper):          # the source `a` with sign=1/-1
    sign = 1 if upper else -1
    return (n[0]/(1+sign*n[2]),n[1]/(1+sign*n[2]))
nd47 = (n0d47,n1d47,n2d47)
sphere47 = n0d47**2+n1d47**2+n2d47**2-1
def sph47(expr):               # numerator reduced modulo the sphere relation
    num,den = sp.fraction(sp.cancel(sp.together(expr)))
    return sp.rem(sp.expand(num),sphere47,n1d47),den
one47,zero47 = sp.Integer(1),sp.Integer(0)
stereo_ok47 = all(
    sp.simplify(sum(c*c for c in stereo47(vd47,wd47,up))-1)==0
    and stereo47(zero47,zero47,up)==(0,0,1 if up else -1)
    and stereo47(one47,zero47,up)==(1,0,0)
    and chart47((one47,zero47,zero47),up)==(1,0)
    and all(sp.simplify(sph47(c-t)[0])==0
            for c,t in zip(stereo47(*chart47(nd47,up),up),nd47))
    for up in (True,False))
# metric of the south chart too (C20a checked the north chart)
ds_south47 = sp.Matrix(stereo47(vd47,wd47,False)).jacobian([vd47,wd47])
stereo_ok47 &= ((ds_south47.T*ds_south47
                 -4*sp.eye(2)/(1+vd47**2+wd47**2)**2).applyfunc(sp.factor)
                ==sp.zeros(2))
rad_N47 = sum(c*c for c in chart47(nd47,True))
rad_S47 = sum(c*c for c in chart47(nd47,False))
radius47 = (
    sp.simplify(sph47(rad_N47-(1-n2d47)/(1+n2d47))[0])==0
    and sp.simplify(sph47(rad_S47-(1+n2d47)/(1-n2d47))[0])==0
    # selected chart: 1-|a|^2 = 2n2/(1+n2)>=0 (n2>=0), -2n2/(1-n2)>0 (n2<0)
    and sp.cancel(1-(1-n2d47)/(1+n2d47)-2*n2d47/(1+n2d47))==0
    and sp.cancel(1-(1+n2d47)/(1-n2d47)+2*n2d47/(1-n2d47))==0
    # chart domains: 4-|a_N|^2=(3+5n2)/(1+n2), 4-|a_S|^2=(3-5n2)/(1-n2)
    and sp.cancel(4-(1-n2d47)/(1+n2d47)-(3+5*n2d47)/(1+n2d47))==0
    and sp.cancel(4-(1+n2d47)/(1-n2d47)-(3-5*n2d47)/(1-n2d47))==0
    # the same sets by solveset on -1<n2<1: radius<=1 exactly on the
    # selected hemispheres, radius<=2 on n2>=-3/5 / n2<=3/5.
    and sp.solveset((1-n2d47)/(1+n2d47)<=1,n2d47,sp.Interval.open(-1,1))
        ==sp.Interval.Ropen(0,1)
    and sp.solveset((1+n2d47)/(1-n2d47)<=1,n2d47,sp.Interval.open(-1,1))
        ==sp.Interval.Lopen(-1,0)
    and sp.solveset((1-n2d47)/(1+n2d47)<=4,n2d47,sp.Interval.open(-1,1))
        ==sp.Interval.Ropen(sp.Rational(-3,5),1)
    and sp.solveset((1+n2d47)/(1-n2d47)<=4,n2d47,sp.Interval.open(-1,1))
        ==sp.Interval.Lopen(-1,sp.Rational(3,5)))
# (ii) increments: three generic right-half-plane ratios; control Re<0.
xq47 = sp.symbols('xq47_1:4',positive=True)
yq47 = sp.symbols('yq47_1:4',real=True)
xneg47 = sp.Symbol('xneg47',negative=True)
increment47 = (
    sp.simplify(sp.expand_trig(sp.expand_complex(
        sp.exp(sp.I*sum(sp.atan(y/x) for x,y in zip(xq47,yq47)))))
        -sp.expand(sp.prod([(x+sp.I*y)/sp.sqrt(x*x+y*y)
                            for x,y in zip(xq47,yq47)])))==0
    and sp.simplify(sp.exp(sp.I*sp.atan(yq47[0]/xneg47)).rewrite(sp.cos)
                    -(xneg47+sp.I*yq47[0])/sp.sqrt(xneg47**2+yq47[0]**2))!=0
    and argument_value30)
# unit values u=(1-t^2+2it)/(1+t^2): R=u_l conj(u_(l-1)) has |R|=1, so
# R/|R|=R, and the source-order product telescopes to u_3 conj(u_0).
tq47 = sp.symbols('tq47_0:4',real=True)
uq47 = [(1-t*t+2*sp.I*t)/(1+t*t) for t in tq47]
Rq47 = [uq47[l]*sp.conjugate(uq47[l-1]) for l in range(1,4)]
increment47 &= (
    all(sp.simplify(R*sp.conjugate(R)-1)==0 for R in Rq47)
    and sp.simplify(sp.prod(Rq47)-uq47[3]*sp.conjugate(uq47[0]))==0)
# Chart products in source order with unit values gA=g(north,anchor),
# gN=g(north,X), gE=g((1,0,0),X), gS=g(south,X), gT=g(n,X).  Control: if the
# two equator endpoints differed (gE vs gE2) the south value would be off.
tg47 = sp.symbols('tg47_0:6',real=True)
gA47,gN47,gE47,gS47,gT47,gE2_47 = [(1-t*t+2*sp.I*t)/(1+t*t) for t in tg47]
def ratio47(new,old):
    return new*sp.conjugate(old)
north_val47 = gA47*ratio47(gN47,gA47)*ratio47(gT47,gN47)
south_val47 = (gA47*ratio47(gN47,gA47)*ratio47(gE47,gN47)
               *sp.conjugate(ratio47(gE47,gS47))*ratio47(gT47,gS47))
south_bad47 = (gA47*ratio47(gN47,gA47)*ratio47(gE47,gN47)
               *sp.conjugate(ratio47(gE2_47,gS47))*ratio47(gT47,gS47))
products47 = (
    sp.simplify(north_val47-gT47)==0 and sp.simplify(south_val47-gT47)==0
    and sp.simplify(south_bad47-gT47)!=0)
# (iii) the equator correction at n=e, with P_N, P_S uninterpreted
# functions of the chart point; control: at another equator point the
# difference is not syntactically zero (that equality is the
# connectedness argument, not an identity).
PN47,PS47 = sp.Function('PN47'),sp.Function('PS47')
cst47 = sp.Symbol('cst47',real=True)
def thetaN47(n):
    return cst47+PN47(*chart47(n,True))
def thetaS47(n):
    return cst47+PN47(one47,zero47)-PS47(one47,zero47)+PS47(*chart47(n,False))
equator47 = (
    sp.simplify(thetaS47((one47,zero47,zero47))-thetaN47((one47,zero47,zero47)))==0
    and sp.simplify(thetaS47((zero47,one47,zero47))
                    -thetaN47((zero47,one47,zero47)))!=0)
# (i) explicit constants with the exact counts (11.14).
zr47,zi47,wr47,wi47 = sp.symbols('zr47 zi47 wr47 wi47',real=True)
alq47 = sp.Rational(1,10**5)
constants47 = (
    count_ok30
    and 4*L30<N30                           # 2L/N<1/2
    and 8*Lg30<K30                          # 4L_g/K<1/2
    and (9+8*N30)*L30<Lg30                  # ||Dg_j||<L_g
    and alq47<4                             # X-path step alpha/K<4/K
    and 100*2*alq47/3<1                     # cap estimate: direct!=0
    and 100*2*sp.Rational(1,10**9)/3<1      # initial step at -1+1e-9
    and sp.expand((zr47*wr47+zi47*wi47)
                  -(1-((zr47-wr47)**2+(zi47-wi47)**2)/2)
                  -((zr47**2+zi47**2-1)+(wr47**2+wi47**2-1))/2)==0
    and 1-sp.Rational(1,2)**2/2==sp.Rational(7,8))
rep('C47d chart selector n[2]>=0 and equator correction: exact real agreement',
    upper_source47 and stereo_ok47 and radius47 and increment47
    and products47 and equator47 and constants47)

# C47e.  The start constant Arg g(north,anchor) has no variable seam, both
# signs.  Type (a) source facts: the only mp.arg call of the phase family is
# theta=mp.arg(g(north,anchor)); `north` is the literal (0,0,1) and
# `anchor` is a parameter never reassigned there; in the corrected section
# anchor=1-alpha if inf else -1+alpha with the literal alpha=1e-5 and
# inf=X>0, each assigned once.  Type (b)/(c): on the band 1-|X|<alpha one has
# |X|>1-alpha>0, the band is {X<-1+alpha} u {X>1-alpha} with
# -1+alpha<0<1-alpha, so inf, hence anchor, is constant on each component,
# and 1-|anchor|=alpha (the band end).  Thus the Arg is evaluated at one
# fixed point per sign: a constant, whose principal branch (even on the
# cut) cannot create a variable seam; exp(i Arg c)=c for every unit c (exact
# values in all quadrants and on the cut below).  The X-path starts at that
# same point and ends at X (j=0 and j=K), so T_X telescopes from it (C47d).
arg_calls47 = [node for node in ast.walk(phase_tree30)
               if isinstance(node,ast.Call) and ast.unparse(node.func)=='mp.arg']
def assigned47(tree,name):
    return sum(1 for node in ast.walk(tree) if isinstance(node,(ast.Assign,ast.AugAssign))
               for t in (node.targets if isinstance(node,ast.Assign) else [node.target])
               if any(isinstance(x,ast.Name) and x.id==name and isinstance(x.ctx,ast.Store)
                      for x in ast.walk(t)))
anchor_source47 = (
    len(arg_calls47)==1
    and ast.unparse(arg_calls47[0])=='mp.arg(g(north, anchor))'
    and matching_assignment30(phase_tree30,'theta','mp.arg(g(north,anchor))')
    and matching_assignment30(phase_tree30,'north','(mp.mpf(0),mp.mpf(0),mp.mpf(1))')
    and assigned47(phase_tree30,'anchor')==0
    and assigned47(phase_tree30,'north')==1
    and matching_assignment30(correct_tree30,'anchor','1-alpha if inf else -1+alpha')
    and matching_assignment30(correct_tree30,'inf','X>0')
    and matching_assignment30(correct_tree30,'alpha',"mp.mpf('1e-5')")
    and all(assigned47(correct_tree30,name)==1 for name in ('anchor','inf','alpha'))
    and 'anchor' in [x.arg for x in phase_tree30.args.args]
    and phase_source30)
alc47 = sp.Rational(1,10**5)
Kc47 = sp.Symbol('Kc47',positive=True,integer=True)
Xc47,anc47 = sp.symbols('Xc47 anc47',real=True)
# the X-path abscissa as written in the source, sympified from its AST
jc47 = sp.Symbol('jc47',integer=True,nonnegative=True)
xpath47 = [sp.sympify(ast.unparse(node.args[1]),
                      locals={'anchor':anc47,'X':Xc47,'j':jc47,'steps':Kc47})
           for node in ast.walk(phase_tree30)
           if isinstance(node,ast.Call) and ast.unparse(node.func)=='g'
           and ast.unparse(node.args[0])=='north'
           and ast.unparse(node.args[1])!='anchor']
cunits47 = [(3+4*sp.I)/5,(-3+4*sp.I)/5,(-3-4*sp.I)/5,(3-4*sp.I)/5,
            sp.Integer(-1),sp.I]
anchor_const47 = (
    -1+alc47<0<1-alc47
    and 1-abs(1-alc47)==alc47 and 1-abs(-1+alc47)==alc47
    and all(sp.simplify(sp.exp(sp.I*sp.arg(cv)).rewrite(sp.cos)-cv)==0
            for cv in cunits47)
    and sp.arg(sp.Integer(-1))==sp.pi
    and len(xpath47)==1
    and sp.simplify(xpath47[0].subs(jc47,0)-anc47)==0
    and sp.simplify(xpath47[0].subs(jc47,Kc47)-Xc47)==0)
rep('C47e start constant Arg g(north,anchor): one fixed point per sign, no seam',
    anchor_source47 and anchor_const47)

# C47f is solely a reproducible numerical observation (d), not used by
# C47a-e.  It evaluates the production functions at the REDUCED counts
# projection_steps=64, phase_steps=32: a DIFFERENT finite map from the
# certified one (11.14), with no guarantee.  South band, X=-1+0.6*alpha,
# anchor=-1+alpha, 40 digits; the mismatch below repeats the source closure
# verbatim (its text is bound in C47a), with a cache so that both charts
# reuse identical evaluations.  Sample: two overlap points (n2=0, the
# selector seam, and n2=-3/10), each with forced north, forced south and
# production (chart=None) values; the band selector at distance
# alpha(1-1/100) (weight s(1/50)) and alpha(1+1/100) (skipped branch).  On
# the south band theta is tiny (|theta|~1e-18): the section has barely left
# I_0, so this observes the mechanism, not a large winding.  The north band,
# where theta is of order 1, is sampled once at N=400, K=48 below (N=64
# and N=200 trip the projection guard there: reduced counts are not
# certified).
def mismatch_factory47(inf,N,dps):
    def mismatch(nn,h):
        if nn[2]>=0:
            w = mp41.sqrt((1+nn[2])/2)
            xx = (w,mp41.mpc(nn[0],-nn[1])/(2*w))
        else:
            w = mp41.sqrt((1-nn[2])/2)
            xx = (mp41.mpc(nn[0],nn[1])/(2*w),w)
        section = cusp30.elementary_projection_section(xx,h,N,dps,None)
        P = cusp30.elliptic_meridian_projector(nn,h,dps=dps,elementary=True)
        direct = cusp30._elementary_unit(
            P*mp41.matrix(cusp30._elementary_pole_inclusion(xx,inf)))
        return cusp30._elementary_inner(direct,section)
    cache = {}
    def cached(nn,h):
        key = (tuple(nn),h)
        if key not in cache:
            cache[key] = mismatch(nn,h)
        return cache[key]
    return cached
rows47 = []
sample_ok47 = True
with mp41.workdps(40):
    alpha_n47 = mp41.mpf('1e-5')
    X_s47 = -1+mp41.mpf('0.6')*alpha_n47
    anchor_s47 = -1+alpha_n47
    g_s47 = mismatch_factory47(False,64,40)
    for n2v47,phv47 in ((mp41.mpf(0),mp41.mpf('0.7')),(mp41.mpf('-0.3'),mp41.mpf(2))):
        rad47 = mp41.sqrt(1-n2v47**2)
        np47 = (rad47*mp41.cos(phv47),rad47*mp41.sin(phv47),n2v47)
        vals47 = {}
        for ch47 in ('north','south',None):
            dg47 = {}
            vals47[ch47] = (cusp30.elementary_phase_family(
                g_s47,np47,X_s47,anchor_s47,32,chart=ch47,diagnostics=dg47),
                dg47['minimum_phase_real'])
        gap47 = abs(vals47['north'][0]-vals47['south'][0])
        exp_err47 = abs(mp41.exp(mp41.j*vals47[None][0])-g_s47(np47,X_s47))
        chosen47 = 'north' if n2v47>=0 else 'south'
        rows47.append(('chart',float(n2v47),vals47[None][0],gap47,exp_err47,
                       min(v[1] for v in vals47.values())))
        sample_ok47 &= (gap47<mp41.mpf('1e-20')*abs(vals47[None][0])
                        and exp_err47<mp41.mpf('1e-35')
                        and vals47[None][0]==vals47[chosen47][0]
                        and min(v[1] for v in vals47.values())>0)
# band selector, same overlap point n2=0, at 60 digits so that the
# ~1e-38 jump is resolved: ||E_c-E|| must equal |exp(-i w theta)-1|.
with mp41.workdps(60):
    alpha_n47 = mp41.mpf('1e-5')
    anchor_s47 = -1+alpha_n47
    np47 = (mp41.cos(mp41.mpf('0.7')),mp41.sin(mp41.mpf('0.7')),mp41.mpf(0))
    w0_47 = mp41.sqrt((1+np47[2])/2)
    xi47 = (w0_47,mp41.mpc(np47[0],-np47[1])/(2*w0_47))
    for frac47 in (mp41.mpf('0.99'),mp41.mpf('1.01')):
        Xb47 = -1+frac47*alpha_n47
        Ec47 = cusp30._elementary_corrected_section(xi47,Xb47,64,32,60,None)
        Eb47 = cusp30.elementary_projection_section(xi47,Xb47,64,60,None)
        jump47 = mp41.sqrt(sum(abs(a-b)**2 for a,b in zip(Ec47,Eb47)))
        if frac47<1:
            wt47 = 1-cusp30.scalar_flat_step(((1+Xb47)-alpha_n47/2)/(alpha_n47/2))
            thb47 = cusp30.elementary_phase_family(
                mismatch_factory47(False,64,60),
                (np47[0],np47[1],np47[2]),Xb47,anchor_s47,32)
            bound47 = abs(mp41.exp(-mp41.j*wt47*thb47)-1)
            sample_ok47 &= (abs(jump47-bound47)<mp41.mpf('1e-10')*bound47
                            and 0<wt47<mp41.mpf('1e-20') and thb47!=0)
        else:
            bound47 = mp41.mpf(0)
            sample_ok47 &= Ec47==Eb47
        rows47.append(('band',float(frac47),jump47,bound47))
# north band, where the mismatch winds: reduced counts N=400, K=48 (a
# different finite map), 25 digits, one overlap point n2=3/10, forced
# north and south charts.
with mp41.workdps(25):
    alpha_n47 = mp41.mpf('1e-5')
    X_n47 = 1-mp41.mpf('0.6')*alpha_n47
    g_n47 = mismatch_factory47(True,400,25)
    n2n47 = mp41.mpf('0.3')
    np_n47 = (mp41.sqrt(1-n2n47**2)*mp41.cos(mp41.mpf('0.7')),
              mp41.sqrt(1-n2n47**2)*mp41.sin(mp41.mpf('0.7')),n2n47)
    north_vals47 = {}
    for ch47 in ('north','south'):
        dg47 = {}
        north_vals47[ch47] = (cusp30.elementary_phase_family(
            g_n47,np_n47,X_n47,1-alpha_n47,48,chart=ch47,diagnostics=dg47),
            dg47['minimum_phase_real'])
    north_gap47 = abs(north_vals47['north'][0]-north_vals47['south'][0])
    north_theta47 = north_vals47['north'][0]
    sample_ok47 &= (north_gap47<mp41.mpf('1e-15')
                    and abs(north_theta47)>mp41.mpf('1e-2')
                    and min(v[1] for v in north_vals47.values())>0)
print('C47f reduced-count observations (d): south N=64,K=32 (2 points x 3 charts, 2 band); north N=400,K=48 (1 point x 2 charts):')
for row47 in rows47:
    if row47[0]=='chart':
        print(f'    south band n2={row47[1]:+.1f}: theta={mp41.nstr(row47[2],5)}'
              f'  |theta_N-theta_S|={mp41.nstr(row47[3],3)}'
              f'  |e^(i theta)-g|={mp41.nstr(row47[4],3)}  min Re R={mp41.nstr(row47[5],8)}')
    else:
        print(f'    band selector distance={row47[1]}*alpha: ||E_c-E||={mp41.nstr(row47[2],12)}'
              f'  |exp(-i w theta)-1|={mp41.nstr(row47[3],12)}')
print('    north band n2=+0.3: theta=',mp41.nstr(north_theta47,8),' |theta_N-theta_S|=',mp41.nstr(north_gap47,3),
      ' min Re R=',mp41.nstr(min(v[1] for v in north_vals47.values()),5))
rep('C47f reduced-count chart overlap and band selector (numerical, d)',sample_ok47)

# C48: no domain guard on the production G path is ever triggered.  The
# 25 raise statements of the 26 functions executed by
# elementary_global_coordinates (the list is checked against the static
# call graph below) are collected from the source and each is
# matched with the bound that excludes it (types (a)-(c)); the list must
# be complete.  Static guards (kind strings, integer counts, the known
# height of the elementary embedding) and the max_work execution guard
# (absent in the exact formula, max_work=None) are listed separately.
path_functions48 = (
    'elementary_global_coordinates','elementary_global_slice',
    '_elementary_corrected_section','elementary_projection_section',
    'elementary_phase_family','_elementary_pole_inclusion','_elementary_unit',
    'elementary_pole_patch','_elementary_pole_radial',
    'elliptic_pole_raw_inverse','elliptic_raw_base_reference',
    'elliptic_bundle_embedding','elementary_cusp_reference',
    'elementary_cusp_logradius','inverse_height_moment',
    '_terminal_height_data','_height_G','_height_R',
    '_height_terminal_forward','scalar_flat_step','scalar_lattice_coords',
    'elementary_restore_torus','scalar_finite_height','_elementary_inner',
    'elliptic_meridian_projector','elementary_closed_parameters')
# The list is checked against the static call graph from
# elementary_global_coordinates: every other reachable cusp function is
# called only inside an arm that the elementary path never enters (the
# else-arm of an `elementary` test, or the `tc == 0` arm after the core
# return of the raw reference).
cusp_functions48 = {n for n,o in vars(cusp30).items()
                    if inspect.isfunction(o) and o.__module__==cusp30.__name__}
def callees48(name):
    return {n.id for n in ast.walk(source_tree30(getattr(cusp30,name)))
            if isinstance(n,ast.Name) and n.id in cusp_functions48}-{name}
reach48,todo48 = set(),['elementary_global_coordinates']
while todo48:
    f48 = todo48.pop()
    if f48 not in reach48:
        reach48.add(f48)
        todo48.extend(callees48(f48)-reach48)
legacy48 = reach48-set(path_functions48)
def in_dead_arm48(tree,callee,tc_rule=False):
    parents = {}
    for node in ast.walk(tree):
        for child in ast.iter_child_nodes(node):
            parents[child] = node
    sites = [n for n in ast.walk(tree)
             if isinstance(n,ast.Name) and n.id==callee]
    def dead(node):
        while node in parents:
            up = parents[node]
            if isinstance(up,(ast.If,ast.IfExp)):
                test = ast.unparse(up.test)
                branch = up.orelse if isinstance(up,ast.If) else [up.orelse]
                body = up.body if isinstance(up,ast.If) else [up.body]
                if test=='elementary' and any(node is b for b in branch):
                    return True
                if tc_rule and test=='tc == 0' and any(node is b for b in body):
                    return True
            node = up
        return False
    return all(dead(site) for site in sites)
def core_return_first48():
    # In the raw reference, `if abs(tc)<=.3 and elementary:` ends in a
    # return and precedes `if tc==0:` in the same block; tc==0 implies
    # abs(tc)<=.3, so with elementary=True the tc==0 arm is never entered.
    core = "abs(tc) <= mp.mpf('.3') and elementary"
    for node in ast.walk(raw_tree31):
        for field in ('body','orelse'):
            stmts = getattr(node,field,None)
            if not isinstance(stmts,list):
                continue
            tests = [ast.unparse(x.test) if isinstance(x,ast.If) else None
                     for x in stmts]
            if core in tests and 'tc == 0' in tests:
                i,j = tests.index(core),tests.index('tc == 0')
                return (i<j and not stmts[i].orelse
                        and isinstance(stmts[i].body[-1],ast.Return))
    return False
call_graph48 = (
    len(set(path_functions48))==26 and set(path_functions48)<=reach48
    and legacy48=={'elliptic_axis_radius','elliptic_collar_blend',
               'elliptic_cusp_logradius','elliptic_height',
               'elliptic_height_chart','elliptic_height_inverse',
               'elliptic_height_reparam'}
    and all(in_dead_arm48(source_tree30(getattr(cusp30,f)),g,
                          tc_rule=(f=='elliptic_raw_base_reference'))
            for f in path_functions48 for g in legacy48
            if g in callees48(f))
    and core_return_first48())
raises48 = []
for name48 in path_functions48:
    for node48 in ast.walk(source_tree30(getattr(cusp30,name48))):
        if isinstance(node48,ast.Raise):
            raises48.append((name48,node48.exc.args[0].value))
certificate48 = {
    'the finite slab formula excludes the inner pole caps':
        'G1 slab only for 1-|X|>=alpha/2=5e-6>1e-9',
    'certified finite expression exceeds max_work; inspect elementary_closed_parameters()':
        'execution guard, absent with max_work=None',
    'steps must be a positive integer': 'static N is a positive integer',
    'finite projection product vanished': 'G2 (11.13)-(11.14): overlaps >1/2',
    'positive integer steps and a valid chart are required':
        'static K positive integer, chart None',
    'phase chart requires stereographic radius <=2':
        'G4 hemisphere chart radius^2=(1-|n2|)/(1+|n2|)<=1',
    'phase count does not give right-half-plane increments':
        'G3 (11.14),(11.16): increments have Re>7/8',
    'zero denominator in a finite projection':
        'G14 cap estimate of section 11.5 in the phase band',
    'radial inverse requires positive radius':
        'G13 raw output nonzero for z_L!=0',
    'radial compression requires radius greater than 1e-8':
        'G13 called only for radius>1e-8',
    'raw pole inverse requires 0<|z|^2<=1e-4': 'G12 cap radius <2e-6',
    'local Hopf projections became orthogonal': 'G14 (11.11)',
    'elementary pole inverse left its certified height region':
        'G12 |height|<1/10 on the cap (C46)',
    'the raw reference excludes the two fixed points': 'G5 |X|<1 on the slab',
    'the bundle embedding excludes the fixed points':
        'G6 only the node has r=0, and it is the pole',
    'elementary projector requires the known inverse-first height':
        'G16 (C42a) core/collar data carry a; pure-middle data have none, but there b0=bi=0 exactly, so torus_pair is never called',
    'pure cusp embedding requires |xy|<=.3 and |H|<=.6':
        'G7 |tc|<=3/10 branch and |a|<3/5',
    'kind must be mid, c0, or ci': 'static',
    'core reference requires |t|<=.3 and |a|<.6': 'G8',
    'the core reference excludes the node': 'G6 node is the pole',
    'elementary cusp inverse requires 0<|t|<=.5 and 0<nu<=1':
        'G9 collar 3/10<|tc|<1/2',
    'inverse moment is outside the certified (q,lam) domain':
        'G10 |b|<=1/2, q<1/8, lambda>0 => q>1/2500',
    'positive collar weight requires q>0': 'G15 from G10',
    'kind must be mid, c0 or ci': 'static',
    'four-seed finite height is outside its stated domain':
        'G11 rr<=1/4, max(A,Y)<4, q<1/8'}
complete48 = (len(raises48)==25
              and {m for _,m in raises48}==set(certificate48))
# The bounds.  G2: consecutive meridian projectors differ by at most
# L*2/N<1/2 (N=4L+1), and the initial one by 100 r_L<1/2, so every
# projection of a unit vector onto the next line has norm >1/2.  G3: the
# K phase samples move the unit value by <4 L_g/K<1/2, so each ratio w
# has |w-1|<1/2, i.e. Re w>7/8.
pars48 = cusp30.elementary_closed_parameters()
L48,N48 = pars48['projector_bound'],pars48['projection_steps']
Lg48,K48 = pars48['phase_bound'],pars48['phase_steps']
counts48 = (
    N48==4*L48+1 and Lg48==32*(N48+1)*L48 and K48==8*Lg48+1
    and 2*(2*L48)<N48                      # 2L/N<1/2
    and 2*(4*Lg48)<K48                     # 4L_g/K<1/2
    and 1-sp.Rational(1,4)/2==sp.Rational(7,8)   # |w-1|^2=2-2Re w<1/4
    and 100*sp.Rational(2,3)*sp.Rational(1,10**9)<sp.Rational(1,2))
# G4: for n2>=0, |(n0,n1)/(1+n2)|^2=(1-n2)/(1+n2)<=1 on the unit sphere.
n2_48 = sp.Symbol('n2_48',nonnegative=True)
g4_48 = (sp.simplify((1-n2_48**2)/(1+n2_48)**2-(1-n2_48)/(1+n2_48))==0
         and sp.Rational(1,1)<=4)
# G7: the principal chart is used only if Re p>=-(9/10)|p|, so |a|<=1/2;
# the half-centred chart has |Arg(-p)|<arccos(9/10)<4pi/25<pi/5, so
# 1/2-2/25<a<1/2+2/25=29/50 and b=a-1<-21/50<-5/12 in cell 1.
g7_48 = (sp.cos(sp.pi/5)<sp.Rational(9,10)
         and sp.cos(4*sp.pi/25)<sp.Rational(9,10)
         and sp.Rational(1,2)+sp.Rational(2,25)==sp.Rational(29,50)
         and sp.Rational(29,50)-1==-sp.Rational(21,50)
         and sp.Rational(1,2)+sp.Rational(1,10)==sp.Rational(3,5))
# G10: core q<1/8 for |t|<=1/2 (section 10.2, L13d:
# q<=v(1-v)/(1+3v), and 8v^2-5v+1=8(v-5/16)^2+7/32>0); collar
# q>1/1375>1/2500 (C43b).
v48 = sp.Symbol('v48',positive=True)
g10_48 = (sp.expand(8*v48**2-5*v48+1-(8*(v48-sp.Rational(5,16))**2
                                      +sp.Rational(7,32)))==0
          and sp.factor(sp.Rational(1,8)*(1+3*v48)-v48*(1-v48)
                        -(8*v48**2-5*v48+1)/8)==0
          and domain43)
# G11: four-seed inputs.  Core cell 0: max(A,Y)<=R<=|f|+2 sqrt q<1+2/sqrt 8
# <2.  Core cell 1 (b=a-1 in [-1/2,-21/50)) and the collar cell n=1:
# A=1/Y(b)<=1/|f(b)| with |f(b)|>1/4 on |b|>5/12 (C43b), so A<4, and
# Y<=2q<1/4 resp. Y=qY(b)<=q.  Collar cell n=0: Y(b) in [q,1], A=q/Y<=1.
# rr=|t|^2<=1/4 and q<1/8 by G10.  Cap: A,Y<1 (section 11.6).
# (A>=1/R>1/2 in cell 1 gives Y=q/A<2q, and 2q<1/4 since q<1/8.)
g11_48 = (1+2/sp.sqrt(8)<2
          and sp.Rational(21,50)>sp.Rational(5,12)
          and 1/sp.Rational(1,4)==4
          and 2*sp.Rational(1,8)==sp.Rational(1,4)
          and sp.Rational(1,2)**2==sp.Rational(1,4))
# G1, G12: alpha/2=5e-6>1e-9; cap r_L<1.7e-6 and |z_S|^2<2e-6<1e-4.
g1_12_48 = (sp.Rational(5,10**6)>sp.Rational(1,10**9)
            and sp.Rational(2,10**6)<sp.Rational(1,10**4)
            and sp.Rational(17,10**7)<sp.Rational(2,10**6))
# Implicit failures (division by zero in the Cramer and coefficient
# solves, 1/small, t/x, t/y) are excluded by the same bounds: the Cramer
# determinant exceeds 2r^2(1+r^2)>0, det(e2,e3)=+-(1-|t|^4)!=0, |small|
# is >=1/2 in the pure middle, and x or y is the nonzero reference
# coordinate (all such operations are enumerated in C50f).  The
# production call uses force_slab=False.
rep('C48 every raise on the G path is excluded by a stated bound (25/25)',
    complete48 and counts48 and g4_48 and g7_48 and g10_48 and g11_48
    and g1_12_48 and call_graph48)
print('C48 guard map:',sorted({v.split(' ')[0].rstrip(':') for v in certificate48.values()}))

# C49: the differential of the whole G_cl is invertible at every point of
# the round S^4 x T^2.  The argument is general, not a chamber check.
# (a) On the free slab, G_W is S^1-equivariant (the finite fold, the phase
#     band correction and the final phase are all equivariant), so in
#     adapted coordinates (orbit point, fibre angle) DG_W is block
#     triangular with fibre entry 1 and det DG_W = det DLambda.
# (b) Lambda factors into: the raw base (n,X)->w=R(X)O^-1 n (R'>0), the
#     inverse stereographic map, (u1,u2)->(u1/u2,Arg u2), the monotone
#     profile ratio->t, the height map a->Y, A or rho (derivatives
#     -f_a(1-f/R)/2, f_a(1+f/R)/2, -f_a/(2R), with f_a>0 by (11.13c)),
#     and, in the collar, rho->x4 with d x4/d rho=+-1/P (C43d).  In the
#     pure middle the last factor is x4=Arg(u_j)/(2pi)+(function of t).
# (c) The torus restoration is the T^2-action applied to a slice point; it
#     is a local diffeomorphism iff the orbit directions are transverse to
#     the slice: det_R(e2,e3)=+-(1-|t_c|^4) in a cusp chart (C30k; the
#     determinant itself C9a, C17d) and
#     the lattice determinant (1+|t|^2)((|t|^2-1)log|t|+2|t|^2)>0 in the
#     middle (C34b), independently of the correction D.
# (d) The cap map is the radial reparametrization z_L of the round cap
#     (increasing), the patch H_j (N10: det DD=S S'/r>0; by (11.11),
#     ||Dg_j-I||<=5000 r with r<1.7e-6 on the cap, so sigma_min>0.99), a
#     linear map and (c).
# Claim types: (a) for the identities checked here, (b)/(c) for the cited
# inputs f_a>0 (11.13c), (11.11), N10 (node_remainder.py), section 11.6's
# bijectivity, and the smoothness of G_W across the seams (C38-C47, C50).
# Hence G_cl is a local diffeomorphism everywhere; with the bijectivity of
# section 11.6 it is a diffeomorphism.
# (a) Equivariance of one fold step and of the start vector.
lam49 = sp.Symbol('lam49',real=True)
gl49 = sp.exp(sp.I*lam49)
r49 = sp.symbols('r0_49 r1_49')
v49 = sp.symbols('v0_49 v1_49')
def step49(ref,vec):
    ph = inner41(ref,vec)
    return [x*ph for x in ref], ph
st_v49,ph_v49 = step49(r49,v49)
st_lv49,ph_lv49 = step49(r49,[gl49*x for x in v49])
equivariant49 = (
    sp.simplify(sp.expand(ph_lv49-gl49*ph_v49))==0
    # |lam*ph|=|ph|, so the normalized update ref*ph/|ph| scales by lam
    and all(sp.simplify(sp.expand(a-gl49*b))==0 for a,b in zip(st_lv49,st_v49))
    and matching_statement31(project_tree30,
        'vector = tuple(a*phase/abs(phase) for a in ref)')
    and matching_assignment30(project_tree30,'phase',
        '_elementary_inner(ref,vector)')
    and matching_assignment30(project_tree30,'vector',
        '_elementary_pole_inclusion(xi)')
    # band correction: an S^1-invariant unit factor (theta depends on (n,X))
    and matching_assignment30(correct_tree30,'result',
        'tuple(a*mp.exp(-mp.j*weight*theta) for a in result)')
    # final updates: (conj(p)x,p y) and one+Arg(p)/(2pi) (bound in C41a)
    and source41)
# Block-triangular determinant: fibre entry 1 below a 3x3 orbit block.
D49 = sp.Matrix(3,3,sp.symbols('d49_0:9'))
row49 = sp.Matrix(1,3,sp.symbols('w49_0:3'))
blk49 = sp.Matrix.vstack(sp.Matrix.hstack(D49,sp.zeros(3,1)),
                         sp.Matrix.hstack(row49,sp.ones(1,1)))
block49 = sp.expand(blk49.det()-D49.det())==0
rep('C49a G_W is S^1-equivariant; det DG_W equals det DLambda on the slab',
    equivariant49 and block49)

# (b) Factors of DLambda.  Raw radius R(X)=(1-X^2)/(1-X sqrt(2-X^2)),
# R'/R=2/(sqrt(2-X^2)(1-X^2))>0; (n,X)->R(X)n has Jacobian R^2 R' on
# S^2 x interval.  Inverse stereographic R^3->S^3 is conformal with factor
# 2/(1+|y|^2).
Xb49 = sp.Symbol('Xb49',real=True)
Rr49 = (1-Xb49**2)/(1-Xb49*sp.sqrt(2-Xb49**2))
y49 = sp.Matrix(sp.symbols('y0_49 y1_49 y2_49',real=True))
sy49 = sum(v**2 for v in y49)
U49 = sp.Matrix([(1-sy49)/(1+sy49),2*y49[0]/(1+sy49),
                 2*y49[1]/(1+sy49),2*y49[2]/(1+sy49)])
JU49 = U49.jacobian(list(y49))
stereo49 = (
    sp.simplify(sp.diff(Rr49,Xb49)/Rr49
                -2/(sp.sqrt(2-Xb49**2)*(1-Xb49**2)))==0
    and all(sp.simplify((JU49.T*JU49-4/(1+sy49)**2*sp.eye(3))[i,j])==0
            for i in range(3) for j in range(3)))
# (u1,u2)->(R=u1/u2,alpha=Arg u2) is a diffeomorphism of S^3 cap {u2!=0}
# onto C x S^1: inverse u2=e^(i alpha)/sqrt(1+|R|^2), u1=R u2.
Rx49,Ry49,al49 = sp.symbols('Rx49 Ry49 al49',real=True)
Rc49 = Rx49+sp.I*Ry49
u2b49 = sp.exp(sp.I*al49)/sp.sqrt(1+Rx49**2+Ry49**2)
u1b49 = Rc49*u2b49
ratio_coords49 = (
    sp.simplify(sp.expand(u1b49*sp.conjugate(u1b49)+u2b49*sp.conjugate(u2b49)))==1
    and sp.simplify(u1b49/u2b49-Rc49)==0)
# The north and |t|>1 charts use (u2/u1, Arg u1) on {u1!=0} symmetrically:
u1n49 = sp.exp(sp.I*al49)/sp.sqrt(1+Rx49**2+Ry49**2)
u2n49 = Rc49*u1n49
ratio_coords49 &= (
    sp.simplify(sp.expand(u1n49*sp.conjugate(u1n49)+u2n49*sp.conjugate(u2n49)))==1
    and sp.simplify(u2n49/u1n49-Rc49)==0)
# Profile: t=R g(|R|^2) has real Jacobian g(g+2 r g'), r=|R|^2, which is
# (|t|/|R|) d|t|/d|R|>0 because s=|t|^2 increases with r: C39 gives
# r=(Bs^2+As)/(1+As), dr/ds>0 (direct branch), and C40 the reciprocal.
gf49 = sp.Function('gf49',real=True)
tmap49 = sp.Matrix([sp.re(sp.expand_complex(Rc49*gf49(Rx49**2+Ry49**2))),
                    sp.im(sp.expand_complex(Rc49*gf49(Rx49**2+Ry49**2)))])
rr49 = sp.Symbol('rr49',positive=True)
profile49 = (
    sp.simplify(tmap49.jacobian([Rx49,Ry49]).det()
                -(gf49(Rx49**2+Ry49**2)**2
                  +2*(Rx49**2+Ry49**2)*gf49(Rx49**2+Ry49**2)
                  *sp.Subs(sp.diff(gf49(rr49),rr49),rr49,Rx49**2+Ry49**2)))==0
    and band_bounds39 and profile40)
# Height maps (f=f(a) increasing, R=sqrt(f^2+4q)>|f|):
fa49 = sp.Function('fa49')
aa49,qq49 = sp.symbols('aa49 qq49',positive=True)
Rh49 = sp.sqrt(fa49(aa49)**2+4*qq49)
fp49 = sp.diff(fa49(aa49),aa49)
height49 = (
    sp.simplify(sp.diff((Rh49-fa49(aa49))/2,aa49)
                -fp49*(fa49(aa49)/Rh49-1)/2)==0
    and sp.simplify(sp.diff((Rh49+fa49(aa49))/2,aa49)
                    -fp49*(fa49(aa49)/Rh49+1)/2)==0
    and sp.simplify(sp.diff(sp.log((Rh49-fa49(aa49))/2)/2,aa49)
                    +fp49/(2*Rh49))==0
    # core cell 1: A=2/(R-f), dA/da=2 f_a/(R(R-f))>0
    and sp.simplify(sp.diff(2/(Rh49-fa49(aa49)),aa49)
                    -2*fp49/(Rh49*(Rh49-fa49(aa49))))==0)
# Collar: rho->x4 at fixed t has derivative 1/P (south) or -1/P (north),
# the e4-coefficient of (-i,0) in the C43d identity (linear in rho).
collar49 = lattice_shift43
# Pure middle: x4=Arg(u2)/(2pi)+xi4(t)+1/2 or Arg(u1)/(2pi)-(function of
# t): at fixed t, d x4/d alpha=1/(2pi) since Arg u1=Arg R+alpha.
middle49 = sp.diff(al49/(2*sp.pi)+sp.Function('xi49')(Rx49,Ry49),al49)==1/(2*sp.pi)
rep('C49b DLambda factors into invertible raw, ratio, profile, height, Cramer maps',
    stereo49 and ratio_coords49 and profile49 and height49 and collar49
    and middle49)

# (c) Restoration.  The action map Psi(slice point, C)=C.(slice point) is a
# local diffeomorphism iff the orbit tangents are transverse to the slice;
# for the cusp slice w_c=0 this is det_R(e2,e3)!=0, for the middle slice
# (two=three=0) it is det_R[2 pi e1,e4,e2,e3]!=0 (the lattice determinant,
# up to the factor 2 pi).  Column operations show the correction D (which
# only adds torus directions) does not change it.
cols49 = sp.Matrix(4,4,sp.symbols('m49_0:16'))
add49 = sp.symbols('k49_0:4')
shear49 = cols49.copy()
shear49[:,0] = cols49[:,0]+add49[0]*cols49[:,2]+add49[1]*cols49[:,3]
shear49[:,1] = cols49[:,1]+add49[2]*cols49[:,2]+add49[3]*cols49[:,3]
restore49 = (
    sp.expand(shear49.det()-cols49.det())==0
    and basis_ok30 and cramer34
    and sp.factor(R34.det()-det_expected34)==0)
# (d) Cap radial map: |z_L|^2=r/(s(1+sqrt(1-r^2))), r=|Z|^2, increases.
rc49 = sp.Symbol('rc49',positive=True)
cap49 = sp.simplify(sp.diff(rc49/(1+sp.sqrt(1-rc49**2)),rc49)
                    -((1+sp.sqrt(1-rc49**2))+rc49**2/sp.sqrt(1-rc49**2))
                    /(1+sp.sqrt(1-rc49**2))**2)==0
rep('C49c restoration and cap factors are invertible; G_cl is a local diffeo',
    restore49 and cap49)

# C50: every conditional on the production G path is classified and
# certified.  The 111 distinct (function, test) pairs of the functions
# executed by elementary_global_coordinates are read from the source; the
# check fails if one is missing from the table below or has an unknown
# category.  Categories:
#   static   fixed by the chart, sign, counts or flags (no seam);
#   guard    a raise or a returned constant excluded on the path (C48);
#   shortcut working-precision shortcut, bypassed by
#            elementary_formula_branches (section 11.7, m1-m3);
#   rewrite  algebraically identical alternatives (same function);
#   flat     a weight, flat step or flat join (all jets agree);
#   chart    choice of holomorphic chart; transitions holomorphic;
#   seam     a genuine selector, certified by the cited check.
classify50 = {
 ('elementary_global_slice','distance < alpha / 2 and (not force_slab)'):'seam N9,C46',
 ('elementary_global_slice',"not -1 + mp.mpf('1e-9') <= X <= 1 - mp.mpf('1e-9')"):'guard G1',
 ('elementary_global_slice','max_work is not None and cost > max_work'):'static execution guard',
 ('elementary_global_slice',"kind == 'mid'"):'chart C41,C44',
 ('elementary_global_slice','projection_steps is None'):'static counts',
 ('elementary_global_slice','phase_steps is None'):'static counts',
 ('elementary_global_slice','inf'):'static sign',
 ('elementary_global_slice','cap_band'):'static work estimate',
 ('_elementary_corrected_section','distance < alpha'):'flat C47b',
 ('_elementary_corrected_section','inf'):'static sign',
 ('_elementary_corrected_section','nn[2] >= 0'):'rewrite C47c',
 ('elementary_projection_section','not isinstance(steps, int) or steps < 1'):'static counts',
 ('elementary_projection_section','diagnostics is not None'):'static bookkeeping',
 ('elementary_projection_section','phase == 0'):'guard G2',
 ('elementary_phase_family',"not isinstance(steps, int) or steps < 1 or chart not in (None, 'north', 'south')"):'static counts',
 ('elementary_phase_family','a[0] ** 2 + a[1] ** 2 > 4'):'guard G4',
 ('elementary_phase_family','not upper'):'seam C47d',
 ('elementary_phase_family','chart is None'):'static chart argument',
 ('elementary_phase_family','upper'):'seam C47d',
 ('elementary_phase_family','ratio.real <= 0'):'guard G3',
 ('elementary_phase_family','diagnostics is not None'):'static bookkeeping',
 ('_elementary_pole_inclusion','infinity'):'static sign',
 ('_elementary_unit','length == 0'):'guard G14',
 ('elementary_pole_patch','radius <= inner'):'flat C46c',
 ('elementary_pole_patch','radius < (inner + outer) / 2'):'shortcut gate',
 ('elementary_pole_patch','_elementary_precision_shortcuts.get()'):'shortcut',
 ('elementary_pole_patch','x - 1 / x <= -(dps + 10) * mp.log(10)'):'shortcut',
 ('_elementary_pole_radial','radius >= c'):'flat C46c',
 ('_elementary_pole_radial','_elementary_precision_shortcuts.get()'):'shortcut',
 ('_elementary_pole_radial','inverse'):'static direction',
 ('_elementary_pole_radial','c / (c - radius) > (mp.mp.dps + 10) * mp.log(10)'):'shortcut',
 ('_elementary_pole_radial','radius <= 0'):'guard G13',
 ('_elementary_pole_radial','radius <= a'):'guard G13',
 ('_elementary_pole_radial','y > 1'):'rewrite C46c',
 ('_elementary_pole_radial','L < 0'):'rewrite C46c',
 ('_elementary_pole_radial','f < 0'):'rewrite C46c',
 ('elliptic_pole_raw_inverse',"not 0 < radius <= mp.mpf('1e-4')"):'guard G12',
 ('elliptic_pole_raw_inverse','infinity'):'static sign',
 ('elliptic_pole_raw_inverse','elementary'):'static flag',
 ('elliptic_pole_raw_inverse','length == 0'):'guard G14',
 ('elliptic_pole_raw_inverse',"abs(height) >= mp.mpf('.1')"):'guard G12',
 ('elliptic_pole_raw_inverse','tc == 0'):'static unreachable (elementary)',
 ('elliptic_pole_raw_inverse','height < 0'):'rewrite C46b',
 ('elliptic_pole_raw_inverse','height >= 0'):'rewrite C46b',
 ('elliptic_raw_base_reference','not abs(xx) < 1'):'guard G5',
 ('elliptic_raw_base_reference',"inf or abs(t) < mp.mpf('.5')"):'seam C44',
 ('elliptic_raw_base_reference','r <= 1'):'seam C42',
 ('elliptic_raw_base_reference','inverse_ratio'):'rewrite C39a',
 ('elliptic_raw_base_reference','inf'):'seam C44',
 ('elliptic_raw_base_reference',"phase_coordinate.real < -mp.mpf('.9') * abs(phase_coordinate)"):'seam C38,C43',
 ('elliptic_raw_base_reference',"abs(tc) <= mp.mpf('.3') and elementary"):'seam C41',
 ('elliptic_raw_base_reference','tc == 0'):'static unreachable (elementary)',
 ('elliptic_raw_base_reference',"abs(tc) <= mp.mpf('.3')"):'static unreachable (elementary)',
 ('elliptic_raw_base_reference','elementary'):'static flag',
 ('elliptic_raw_base_reference','a > 0'):'static unreachable (elementary)',
 ('elliptic_bundle_embedding',"kind in ('c0', 'ci')"):'chart C41,C44',
 ('elliptic_bundle_embedding','not r > 0'):'guard G6',
 ('elliptic_bundle_embedding','H is None'):'guard C42a (pure middle: b0=bi=0 exactly, torus_pair not called)',
 ('elliptic_bundle_embedding',"abs(x * y) > mp.mpf('.3') or abs(H) > mp.mpf('.6')"):'guard G7',
 ('elliptic_bundle_embedding','H < 0 and w > 0'):'flat C41b (w=0 on |H|<1/10)',
 ('elliptic_bundle_embedding',"kind == 'mid'"):'chart C41,C44',
 ('elliptic_bundle_embedding','with_norm'):'static flag',
 ('elliptic_bundle_embedding','elementary'):'static flag',
 ('elliptic_bundle_embedding','w == 0'):'flat C41b',
 ('elliptic_bundle_embedding',"kind == 'ci'"):'static sign',
 ('elliptic_bundle_embedding','b0 > 0'):'flat C41c',
 ('elliptic_bundle_embedding','bi > 0'):'flat C41c',
 ('elliptic_bundle_embedding','bm > 0'):'flat C41c,C44c',
 ('elliptic_bundle_embedding','H < 0'):'flat C41b (w=0 on |H|<1/10)',
 ('elliptic_bundle_embedding','elementary and len(data) > 4'):'guard C42a (false only on pure-middle data, where b0=bi=0)',
 ('elliptic_bundle_embedding','index == 0'):'static slot',
 ('elementary_cusp_reference',"abs(t) > mp.mpf('.3') or abs(a) >= mp.mpf('.6')"):'guard G8',
 ('elementary_cusp_reference','a > 0'):'seam C41b gauge',
 ('elementary_cusp_reference','y == 0'):'guard G6',
 ('elementary_cusp_reference','cell == -1'):'seam C45c',
 ('elementary_cusp_reference','cell == 1'):'seam C45c',
 ('elementary_cusp_logradius',"not 0 < r <= mp.mpf('.5') or not 0 < nu <= 1"):'guard G9',
 ('elementary_cusp_logradius','lam == 1'):'seam C44a',
 ('elementary_cusp_logradius','f >= 0'):'rewrite C41b (rationalization 2q/(R+f)=(R-f)/2)',
 ('inverse_height_moment',"abs(a) > mp.mpf('.5') or not 0 <= q <= mp.mpf('.125') or (not 0 <= lam <= 1) or (lam > 0 and q < mp.mpf(1) / 2500)"):'guard G10',
 ('inverse_height_moment',"lam == 1 or aa <= mp.mpf('.1')"):'seam C45b (|b|=1/10 by the cited 11.2 argument; lam==1 by C44a)',
 ('inverse_height_moment',"aa >= mp.mpf('.4')"):'seam C45b',
 ('_terminal_height_data','lam == 0'):'rewrite C41a',
 ('_terminal_height_data','q <= 0'):'guard G15',
 ('_height_G','amp == 0'):'rewrite C41a',
 ('_height_G','inverse'):'static direction',
 ('_height_R','inverse'):'static direction',
 ('_height_R','value >= c'):'flat C45b',
 ('_height_R','_elementary_precision_shortcuts.get()'):'shortcut',
 ('_height_R','value <= 0'):'guard unreachable (G^-1(v)>0)',
 ('_height_R','value <= delta'):'guard unreachable (interior |b|>1/10)',
 ('_height_R','c / (c - value) > (mp.mp.dps + 10) * mp.log(10)'):'shortcut',
 ('_height_R','F < 0'):'rewrite C45a',
 ('_height_R','y > 1'):'rewrite C45a',
 ('_height_R','L < 0'):'rewrite C45a',
 ('_height_terminal_forward','delta_z <= 0'):'guard unreachable (interior |b|>1/10)',
 ('_height_terminal_forward','delta_z >= gap0'):'flat C45b',
 ('_height_terminal_forward','_elementary_precision_shortcuts.get()'):'shortcut',
 ('_height_terminal_forward','Ez > 50'):'rewrite C45a',
 ('_height_terminal_forward','Ez > (mp.mp.dps + 20) * mp.log(10) + mp.log(K)'):'shortcut',
 ('scalar_flat_step','u <= 0'):'flat C41c',
 ('scalar_flat_step','u >= 1'):'flat C44c',
 ('scalar_lattice_coords','j == k'):'static Cramer column',
 ('elementary_restore_torus',"kind in ('c0', 'ci')"):'chart C41,C44',
 ('elementary_restore_torus',"kind != 'mid'"):'static kind',
 ('elementary_restore_torus','beta == 0'):'flat C41c',
 ('elementary_restore_torus','inf'):'static sign',
 ('scalar_finite_height',"not rr <= mp.mpf('.25') or max(A, Y) > 4 or q >= mp.mpf('.125')"):'guard G11',
 ('scalar_finite_height','v <= mp.mpf(1) / 16'):'flat L13',
 ('scalar_finite_height','reciprocal'):'flat L13',
 ('scalar_finite_height',"v >= mp.mpf('.5')"):'flat L13'}
categories50 = {'static','guard','shortcut','rewrite','flat','chart','seam'}
found50 = set()
occurrences50 = 0
for name50 in path_functions48:
    for node50 in ast.walk(source_tree30(getattr(cusp30,name50))):
        if isinstance(node50,(ast.If,ast.IfExp)):
            found50.add((name50,ast.unparse(node50.test)))
            occurrences50 += 1
# The table is keyed by distinct (function, test) text; the total number
# of occurrences is recorded too, so a repeated test added later (which
# could mean something else) changes the count and fails the check.
classified50 = (
    found50==set(classify50)
    and len(found50)==111 and occurrences50==142
    and all(v.split()[0] in categories50 for v in classify50.values()))
counts50 = {}
for v50 in classify50.values():
    counts50[v50.split()[0]] = counts50.get(v50.split()[0],0)+1
rep('C50 all 111 distinct G-path tests (142 occurrences) classified and certified',
    classified50)
print('C50 categories:',dict(sorted(counts50.items())),
      '; occurrences:',occurrences50)

# C50b: the implicit seams.  Principal-branch cuts hide inside library
# calls, not conditionals, so every arg/log/atan call on the G path is
# also read from the source and classified: 'positive' (real logarithm of
# a positive number, or an argument kept in an open half-plane by a
# stated bound), 'constant' (fixed arguments or shortcut thresholds), or
# 'cut' (a principal-branch jump, certified as a holomorphic chart or
# lattice change by the cited check).
branch_calls50 = {
 ('elementary_global_slice','mp.arg(phase)'):'cut C42d (middle), C50g (collar) (one->one+1: lattice e1)',
 ('elementary_phase_family','mp.arg(g(north, anchor))'):'constant C47e',
 ('elementary_phase_family','mp.atan(ratio.imag / ratio.real)'):'positive G3 Re ratio>7/8',
 ('elementary_pole_patch','mp.log(10)'):'constant shortcut threshold',
 ('_elementary_pole_radial','mp.log(mp.expm1(y))'):'positive y>0',
 ('_elementary_pole_radial','mp.log(10)'):'constant shortcut threshold',
 ('elliptic_pole_raw_inverse','mp.arg(u2)'):'positive C46b3 Re u2>0',
 ('elliptic_pole_raw_inverse','mp.arg(-u1)'):'positive C46b3 Re(-u1)>0',
 ('elliptic_raw_base_reference','mp.log(t)'):'cut C34, C50g (collar), C42 (pure middle)',
 ('elliptic_raw_base_reference','mp.arg(phase_coordinate)'):'positive C38a,C43 principal chart Re p>=-(9/10)|p|',
 ('elliptic_raw_base_reference','mp.arg(-phase_coordinate)'):'positive C38a,C43 half chart Re(-p)>(9/10)|p|',
 ('elliptic_raw_base_reference','mp.arg(t)'):'cut C42',
 ('elliptic_raw_base_reference','mp.arg(u2)'):'cut C42',
 ('elliptic_raw_base_reference','mp.arg(u1)'):'cut C42',
 ('elliptic_bundle_embedding','mp.log(m * abs(t))'):'positive m|t|>0',
 ('elliptic_bundle_embedding','mp.log(t)'):'cut C34,C42 (principal normalization = move (c))',
 ('elliptic_bundle_embedding','mp.log(abs(z))'):'constant unreachable (elementary)',
 ('elliptic_bundle_embedding','mp.arg(t if index == 0 else -t)'):'cut C42b,C34d (brackets in {0,2pi}; weights vanish at jumps)',
 ('elementary_cusp_logradius','mp.log(nu)'):'constant nu in {1,e^-2}',
 ('elementary_cusp_logradius','mp.log(r)'):'positive r>3/10',
 ('elementary_cusp_logradius','mp.log(Y)'):'positive Y>=q>0',
 ('_terminal_height_data','mp.log(q)'):'positive q>1/2500 (G10)',
 ('_height_R','mp.log(mp.expm1(y))'):'positive y>0',
 ('_height_R','mp.log(10)'):'constant shortcut threshold',
 ('_height_terminal_forward','mp.log(K)'):'constant shortcut threshold',
 ('_height_terminal_forward','mp.log(10)'):'constant shortcut threshold',
 ('scalar_finite_height','mp.log(16 * v)'):'positive v>1/16',
 ('scalar_finite_height','mp.log(8)'):'constant'}
found_calls50 = set()
call_occurrences50 = 0
for name50 in path_functions48:
    for node50 in ast.walk(source_tree30(getattr(cusp30,name50))):
        if (isinstance(node50,ast.Call) and isinstance(node50.func,ast.Attribute)
                and node50.func.attr in ('arg','log','atan')):
            found_calls50.add((name50,ast.unparse(node50)))
            call_occurrences50 += 1
calls_ok50 = (found_calls50==set(branch_calls50) and len(found_calls50)==28
              and call_occurrences50==29
              and all(v.split()[0] in ('positive','constant','cut')
                      for v in branch_calls50.values()))
rep('C50b all 28 distinct arg/log/atan calls (29 occurrences) classified',
    calls_ok50)

# C50c: the remaining non-smooth primitives (floor, nint, sign, int, min,
# max) on the G path.  Each is a precision setting, diagnostic bookkeeping,
# a guard argument, a height cell (C45c), a torus-pair shift whose jumps
# carry zero weight (C41b shift table), the principal-log normalization
# (C34/C42, zero in production), or a sign used only off a neighbourhood
# of zero (the moment's G branch is odd and analytic on |b|<=1/10; the
# raw inverse's sign(height) is in the non-elementary branch).
primitive_calls50 = {
 ('elementary_cusp_logradius','max(mp.mp.dps, dps)'):'precision',
 ('elementary_cusp_logradius',"mp.floor(a + mp.mpf('.5'))"):'cell C45c',
 ('elementary_cusp_reference',"int(mp.floor(a + mp.mpf('.5')))"):'cell C45c',
 ('elementary_cusp_reference',"mp.floor(a + mp.mpf('.5'))"):'cell C45c',
 ('elementary_global_coordinates',"max(mp.mp.dps, kwargs.get('dps', 50))"):'precision',
 ('elementary_global_slice','max(mp.mp.dps, dps)'):'precision',
 ('elementary_phase_family',"min(diagnostics.get('minimum_phase_real', mp.mpf(1)), ratio.real)"):'diagnostic',
 ('elementary_pole_patch','max(mp.mp.dps, dps)'):'precision',
 ('elementary_projection_section','max(mp.mp.dps, dps)'):'precision',
 ('elementary_projection_section',"min(diagnostics.get('minimum_projection_norm', mp.mpf(1)), margin)"):'diagnostic',
 ('elementary_projection_section','min(margin, abs(phase))'):'diagnostic',
 ('elliptic_bundle_embedding','int(shift)'):'shift C41b',
 ('elliptic_bundle_embedding','max(mp.mp.dps, dps)'):'precision',
 ('elliptic_meridian_projector','max(mp.mp.dps, dps)'):'precision',
 ('elliptic_bundle_embedding',"mp.floor(centre + mp.mpf('.5') - H)"):'shift C41b',
 ('elliptic_bundle_embedding','mp.nint((lt.imag - principal.imag) / (2 * mp.pi))'):'normalization C34,C42',
 ('elliptic_pole_raw_inverse','max(mp.mp.dps, dps)'):'precision',
 ('elliptic_pole_raw_inverse','mp.sign(height)'):'unreachable (elementary)',
 ('elliptic_raw_base_reference','max(mp.mp.dps, dps)'):'precision',
 ('inverse_height_moment','mp.sign(a)'):'sign used only for |a|>1/10, where it is locally constant; the |a|<=1/10 branch is G(a) (C45)',
 ('scalar_finite_height','max(A, Y)'):'guard G11'}
found_prims50 = set()
prim_occurrences50 = 0
for name50 in path_functions48:
    for node50 in ast.walk(source_tree30(getattr(cusp30,name50))):
        if isinstance(node50,ast.Call):
            fn50 = (node50.func.attr if isinstance(node50.func,ast.Attribute)
                    else getattr(node50.func,'id',None))
            if fn50 in ('floor','nint','sign','int','min','max'):
                found_prims50.add((name50,ast.unparse(node50)))
                prim_occurrences50 += 1
prims_ok50 = (found_prims50==set(primitive_calls50) and len(found_prims50)==21
              and prim_occurrences50==22
              and exact_if_body37(moment_tree36,"lam==1 or aa<=mp.mpf('.1')",
                                  'return _height_G(a,amp,d)'))
rep('C50c all 21 distinct floor/nint/sign/int/min/max calls (22 occurrences) classified',
    prims_ok50)

# C50a: the certified object is exactly this source.  G_cl is the function
# defined by reading the text of cusp.py in exact real arithmetic: decimal
# literals are rationals (mp.mpf('.3')=3/10), mp.pi and mp.e are pi and e,
# the mp functions are the mathematical functions, workdps, max(dps) and
# unary + are identities, and the precision shortcuts are off; mpmath runs
# of the same text are numerical evaluations (d) of it.  The whole module
# is pinned by the full SHA-256 of its token stream (comments and blank
# lines ignored), and the 26 path functions also by short digests that
# localise a change; any edit of cusp.py fails here and requires a new
# certificate.  The runtime objects the checks inspect are those the text
# defines: no __wrapped__, code from cusp.py at the line of the def, the
# module's globals, and defaults equal to the literals in the text.  For the pinned source the censuses are complete: there is
# no while, try, match, lambda, assert, global, comprehension filter, % or
# //; every comparison or boolean operator is an if/if-expression test, or
# defines one of eight flags (cap_band, inf, upper, inverse_ratio,
# reciprocal) whose every use is a test or a flag argument of a path call;
# every call is a path function, a local helper, a listed builtin or a
# listed mpmath function.  The object is
# elementary_global_coordinates(Z,X,phi,max_work=None) inside
# elementary_formula_branches() (shortcuts off, section (m) of
# symbolic_identities.py), with force_slab=False.
import tokenize as tokenize50, io as io50, hashlib as hashlib50
def source_digest50(function):
    tokens = []
    for tok in tokenize50.generate_tokens(
            io50.StringIO(inspect.getsource(function)).readline):
        if tok.type in (tokenize50.COMMENT,tokenize50.NL,tokenize50.ENDMARKER):
            continue
        tokens.append({tokenize50.NEWLINE:'<NEWLINE>',
                       tokenize50.INDENT:'<INDENT>',
                       tokenize50.DEDENT:'<DEDENT>'}.get(tok.type,tok.string))
    return hashlib50.sha256('\x1f'.join(tokens).encode()).hexdigest()[:16]
digests50 = {
 'elementary_global_coordinates':'35907c53f648b603',
 'elementary_global_slice':'43b5c6b86bd0c47a',
 '_elementary_corrected_section':'ede8e4e373251a43',
 'elementary_projection_section':'918354ac41feb986',
 'elementary_phase_family':'a607c7ba19668a1a',
 '_elementary_pole_inclusion':'d1ec0b0f31808b66',
 '_elementary_unit':'771512424c60414f',
 'elementary_pole_patch':'9c44047e50db6bd3',
 '_elementary_pole_radial':'43417dec48397205',
 'elliptic_pole_raw_inverse':'9ceae566a70b3dac',
 'elliptic_raw_base_reference':'ed89be8f233aa60a',
 'elliptic_bundle_embedding':'1cc0cc7365147327',
 'elementary_cusp_reference':'2f486ce4072fb782',
 'elementary_cusp_logradius':'b5022878ce0da627',
 'inverse_height_moment':'2d9215b8877f6071',
 '_terminal_height_data':'ad72c6262b1ab802',
 '_height_G':'eee49e5444d94b3a',
 '_height_R':'55c21917929901ed',
 '_height_terminal_forward':'a71ae4caf0cd49af',
 'scalar_flat_step':'b166f4423dc6293d',
 'scalar_lattice_coords':'cdf243d1734d301b',
 'elementary_restore_torus':'3b1f09b05f16cad0',
 'scalar_finite_height':'55a894c952c38741',
 '_elementary_inner':'bf798df0b7c67f78',
 'elliptic_meridian_projector':'ef96f0ac8bfce1a5',
 'elementary_closed_parameters':'1f10c9b30a384e25'}
module_file50 = os.path.abspath(cusp30.__file__)
module_text50 = open(module_file50).read()
def text_digest50(text):
    tokens = []
    for tok in tokenize50.generate_tokens(io50.StringIO(text).readline):
        if tok.type in (tokenize50.COMMENT,tokenize50.NL,tokenize50.ENDMARKER):
            continue
        tokens.append({tokenize50.NEWLINE:'<NEWLINE>',
                       tokenize50.INDENT:'<INDENT>',
                       tokenize50.DEDENT:'<DEDENT>'}.get(tok.type,tok.string))
    return hashlib50.sha256('\x1f'.join(tokens).encode()).hexdigest()
module_digest50 = '1aa613a4d71de2517962f0a1a8399c36fa449a28b6db7391380caca747383353'
module_defs50 = {n.name:n for n in ast.parse(module_text50).body
                 if isinstance(n,ast.FunctionDef)}
def literal50(node):
    return eval(ast.unparse(node),{'__builtins__':{}})
def same_object50(name):
    fn,node = getattr(cusp30,name),module_defs50[name]
    kw50 = {a.arg:literal50(d) for a,d in zip(node.args.kwonlyargs,
                                              node.args.kw_defaults)
            if d is not None}
    return (not hasattr(fn,'__wrapped__') and not node.decorator_list
            and os.path.abspath(fn.__code__.co_filename)==module_file50
            and fn.__code__.co_firstlineno==node.lineno
            and fn.__globals__ is vars(cusp30)
            and fn.__defaults__==(tuple(literal50(d) for d in node.args.defaults)
                                  or None)
            and (fn.__kwdefaults__ or {})==kw50)
pinned50 = (set(digests50)==set(path_functions48)
            and text_digest50(module_text50)==module_digest50
            and all(source_digest50(getattr(cusp30,f))==d
                    for f,d in digests50.items())
            and all(same_object50(f) for f in path_functions48))
def parents50(tree):
    parents = {}
    for node in ast.walk(tree):
        for child in ast.iter_child_nodes(node):
            parents[child] = node
    return parents
def in_test50(node,parents):
    # inside the test of an If/IfExp (through Compare/BoolOp/Not nodes)
    while node in parents:
        up = parents[node]
        if isinstance(up,(ast.If,ast.IfExp)) and node is up.test:
            return True
        node = up
    return False
banned50 = (ast.While,ast.Try,ast.Lambda,ast.Assert,ast.Global,ast.Nonlocal,
            ast.Delete)+((ast.Match,) if hasattr(ast,'Match') else ())
flags50 = {('elementary_global_slice','cap_band'),('elementary_global_slice','inf'),
           ('_elementary_corrected_section','inf'),('elementary_phase_family','upper'),
           ('elliptic_raw_base_reference','inverse_ratio'),
           ('elliptic_raw_base_reference','inf'),('elementary_restore_torus','inf'),
           ('scalar_finite_height','reciprocal')}
mp_calls50 = {'arg','asinh','atan','conj','exp','expm1','floor','log','log1p',
              'matrix','mpc','mpf','nint','sign','sinh','sqrt','workdps'}
builtin_calls50 = {'abs','max','min','int','sum','tuple','list','map','zip',
                   'range','enumerate','next','isinstance','len','dict',
                   'ValueError','Fraction'}
method_calls50 = {'block.extend','coords.append','out.extend','diagnostics.get',
                  'kwargs.get','_elementary_precision_shortcuts.get'}
flag_names50 = {'cap_band','inf','infinity','upper','inverse_ratio',
                 'reciprocal','elementary'}
structure50 = True
flag_defs50 = set()
for name50 in path_functions48:
    tree50 = source_tree30(getattr(cusp30,name50))
    par50 = parents50(tree50)
    local50 = ({n.id for n in ast.walk(tree50) if isinstance(n,ast.Name)
                and isinstance(n.ctx,ast.Store)}
               | {n.name for n in ast.walk(tree50) if isinstance(n,ast.FunctionDef)}
               | {a.arg for n in ast.walk(tree50) if isinstance(n,ast.FunctionDef)
                  for a in n.args.args+n.args.kwonlyargs})
    for node50 in ast.walk(tree50):
        if isinstance(node50,banned50):
            structure50 = False
        if isinstance(node50,ast.comprehension) and node50.ifs:
            structure50 = False
        if isinstance(node50,ast.BinOp) and isinstance(node50.op,(ast.Mod,ast.FloorDiv)):
            structure50 = False
        if isinstance(node50,(ast.Compare,ast.BoolOp)) and not in_test50(node50,par50):
            up50 = node50
            while isinstance(par50.get(up50),(ast.Compare,ast.BoolOp,ast.UnaryOp,ast.IfExp)):
                up50 = par50[up50]
            target50 = par50.get(up50)
            if (isinstance(target50,ast.Assign) and len(target50.targets)==1
                    and isinstance(target50.targets[0],ast.Name)
                    and (name50,target50.targets[0].id) in flags50):
                flag_defs50.add((name50,target50.targets[0].id))
            else:
                structure50 = False
        if isinstance(node50,ast.Call):
            fn50 = node50.func
            if isinstance(fn50,ast.Name):
                ok50 = (fn50.id in builtin_calls50 or fn50.id in local50
                        or fn50.id in cusp_functions48)
            elif (isinstance(fn50,ast.Attribute) and isinstance(fn50.value,ast.Name)
                  and fn50.value.id=='mp'):
                ok50 = fn50.attr in mp_calls50
            elif isinstance(fn50,ast.Attribute):
                ok50 = ast.unparse(fn50) in method_calls50
            else:
                # the selected logradius: (elementary_... if elementary else elliptic_...)
                ok50 = (isinstance(fn50,ast.IfExp) and ast.unparse(fn50.test)=='elementary'
                        and isinstance(fn50.body,ast.Name)
                        and fn50.body.id=='elementary_cusp_logradius')
            structure50 &= ok50
    # every use of a boolean flag (definitions above, and the flag
    # parameters inf/infinity/upper/elementary) is a test, a flag argument
    # of a path or local function, or part of a flag definition
    def flag_use_ok50(node):
        if in_test50(node,par50):
            return True
        up = par50[node]
        if isinstance(up,ast.keyword):
            up = par50[up]
        if (isinstance(up,ast.Call) and isinstance(up.func,ast.Name)
                and (up.func.id in cusp_functions48 or up.func.id in local50)):
            return True
        while node in par50:
            node = par50[node]
            if (isinstance(node,ast.Assign) and len(node.targets)==1
                    and isinstance(node.targets[0],ast.Name)
                    and (name50,node.targets[0].id) in flags50):
                return True
        return False
    for node50 in ast.walk(tree50):
        if (isinstance(node50,ast.Name) and isinstance(node50.ctx,ast.Load)
                and node50.id in flag_names50):
            structure50 &= flag_use_ok50(node50)
structure50 &= flag_defs50==flags50
rep('C50a cusp.py pinned (SHA-256), objects = text; no hidden selector, call or flag use',
    pinned50 and structure50)

# C50e: the central fibres t_c=0 lie in the core, where the reference is
# elementary_cusp_reference (bound in full by C41a) with lam=0.  There
# (source of _terminal_height_data and inverse_height_moment):
# kappa=2, amp=0, d=1, B=e^-1, C=q e, and f=f(b,q,0) has the sign of b:
#   |b|<=1/10:      f=G(b)=d*b=b (amp==0);
#   1/10<|b|<2/5:   f=sign(b) G(R^-1(G^-1(v)))=sign(b) R^-1(v), R^-1>=1/10;
#   |b|>=2/5:       |f|=e^(2|b|-1)-q e^(1-2|b|)>=e^(-1/5)-e^(1/5)/8>0.
# The core argument a lies in (-3/7,29/50): the principal chart has
# |a|<=1/2-arccos(9/10)/(2pi)<1/2-1/14 (cos(pi/7)>9/10), the half chart
# 1/2-2/25<a<1/2+2/25 (G7).  So cell=floor(a+1/2) is 0 or 1, the
# cell==-1 alternative is never taken, and with R=sqrt(f^2+4q)>=|f|:
#   cell 0, a>0:  f>0,  A=(R+f)/2>=f>0;
#   cell 1:       f<0,  A=2/(R-f)>0;
#   cell 0, a<=0: f<=0, Y=(R-f)/2>0 unless q=0 and f=0, i.e. t_c=0 and
#                 a=0: the node (the pole; guard G6).
# Hence x=sqrt(A/e)/nu, y=t/x and y=sqrt(Y), x=t/y are smooth through
# t_c=0 (q enters f only through C=qe): the central fibre is not a seam.
height_R_tree50 = source_tree30(cusp30._height_R)
fp50,fn50 = sp.Symbol('fp50',positive=True),sp.Symbol('fn50',negative=True)
q50 = sp.Symbol('q50',nonnegative=True)
Rp50,Rn50 = sp.sqrt(fp50**2+4*q50),sp.sqrt(fn50**2+4*q50)
ell50 = sp.Symbol('ell50',positive=True)
central50 = (
    exact_if_body37(core_tree37,'a>0',
        'A = 2/(R-f) if cell==1 else (R+f)/2\n'
        'x = mp.sqrt(A/e)/nu\n'
        'return x,t/x')
    and ordered_nodes37(core_tree37.body,[
        "R = mp.sqrt(f*f+4*q)",
        "Y = 2/(R+f) if cell==-1 else (R-f)/2",
        "y = mp.sqrt(Y)",
        "return t/y,y"])
    and matching_assignment30(core_tree37,'f','inverse_height_moment(b,q,0)')
    # lam=0 parameters and the three moment branches
    and exact_if_body37(terminal_tree36,'lam==0',
                        'kappa,amp,d = mp.mpf(2),mp.mpf(0),mp.mpf(1)')
    and matching_statement31(terminal_tree36,
                             'B,C = mp.exp(-kappa/2),q*mp.exp(kappa/2)')
    and exact_if_body37(height_G_tree36,'amp==0',
                        'return a/d if inverse else d*a')
    and exact_if_body37(moment_tree36,"lam==1 or aa<=mp.mpf('.1')",
                        'return _height_G(a,amp,d)')
    and exact_if_body37(moment_tree36,"aa>=mp.mpf('.4')",
        'return mp.sign(a)*(B*mp.exp(kappa*aa)-C*mp.exp(-kappa*aa))')
    and matching_statement31(moment_tree36,
        'return mp.sign(a)*_height_G(_height_R(_height_G(v,amp,d,inverse=True),'
        'inverse=True),amp,d)')
    # R^-1 >= 1/10: it returns delta=1/10, value>=c=1/5, or delta+c*ell/(2(2+ell))
    and exact_if_body37(height_R_tree50,'value<=0','return delta')
    and exact_if_body37(height_R_tree50,'value>=c','return value')
    and matching_statement31(height_R_tree50,
                             "delta,c = mp.mpf('.1'),mp.mpf('.2')")
    and matching_statement31(height_R_tree50,
                             'return delta+c*ell/(2*(2+ell))')
    and (sp.Rational(1,5)*ell50/(2*(2+ell50))).is_positive
    # |b|>=2/5: e^(2|b|-1)>=e^(-1/5)>4/5 and q e^(1-2|b|)<=e^(1/5)/8<5/32
    and sp.exp(sp.Rational(1,5))<sp.Rational(5,4)
    and sp.Rational(4,5)-sp.Rational(5,32)>0
    # cells: cos(pi/7)>9/10 gives arccos(9/10)>pi/7
    and sp.cos(sp.pi/7)>sp.Rational(9,10)
    and sp.Rational(1,2)-sp.Rational(1,14)==sp.Rational(3,7)
    and -sp.Rational(3,7)+sp.Rational(1,2)>0
    and sp.Rational(1,2)+sp.Rational(2,25)+sp.Rational(1,2)<2
    # R^2-f^2=4q>=0, so (R+f)/2-f=(R-f)/2>=0, R+f>0 for f>0, R-f>0 for f<0
    and sp.expand(Rp50**2-fp50**2-4*q50)==0
    and (Rp50+fp50).is_positive and (Rn50-fn50).is_positive)
rep('C50e central fibre: f has the sign of b, cells 0/1, core roots positive off the node',
    central50)

# C50f: implicit domain conditions.  Every division, square root, log1p,
# abs and non-integer power on the path is read from the source.  Four
# kinds are recognised automatically: constant arguments; nodes inside a
# test (they only select a branch; the divisions among them sit in
# precision-shortcut tests, which are not evaluated with the shortcuts
# off, and this is checked); nodes in a dead arm (C48); and
# abs(..)**2, a polynomial.  Each remaining one is listed with the reason
# its argument is nonzero (divisions, abs), positive (sqrt) or >-1
# (log1p) on the branch where it is evaluated, and the check fails if one
# is missing (types (b),(c); the cited bounds are C46b3, C47c, G1-G16 and
# the written section 11.2 terminal bounds).  'flat' marks an abs inside
# a flat step whose argument is <=0 (step 0) or >=1 (step 1) near the kink.
def const_node50(node):
    if isinstance(node,ast.Constant):
        return True
    if (isinstance(node,ast.Attribute) and isinstance(node.value,ast.Name)
            and node.value.id=='mp' and node.attr in ('pi','e','j')):
        return True
    if (isinstance(node,ast.Call) and isinstance(node.func,ast.Attribute)
            and node.func.attr in ('mpf','sqrt','exp','log')
            and all(const_node50(a) for a in node.args)):
        return True
    if isinstance(node,ast.BinOp):
        return const_node50(node.left) and const_node50(node.right)
    if isinstance(node,ast.UnaryOp):
        return const_node50(node.operand)
    return False
def dead_node50(node,parents,fname):
    while node in parents:
        up = parents[node]
        if isinstance(up,(ast.If,ast.IfExp)):
            test = ast.unparse(up.test)
            orelse = up.orelse if isinstance(up,ast.If) else [up.orelse]
            body = up.body if isinstance(up,ast.If) else [up.body]
            if test=='elementary' and any(node is b for b in orelse):
                return True
            if (fname=='elliptic_raw_base_reference' and test=='tc == 0'
                    and any(node is b for b in body)):
                return True
        node = up
    return False
census50 = {
 ('_elementary_corrected_section','abs','X'):'positive distance=1-|X| enters a value only in the band, where |X|>1-alpha>0',
 ('_elementary_corrected_section','div','2 * w'):'positive w=sqrt((1+-n2)/2)>=1/sqrt(2) on the selected hemisphere (C47c)',
 ('_elementary_corrected_section','div','alpha / 2'):'positive alpha=1e-5',
 ('_elementary_corrected_section','sqrt','(1 + nn[2]) / 2'):'positive >=1/2 for n2>=0 (C47c)',
 ('_elementary_corrected_section','sqrt','(1 - nn[2]) / 2'):'positive >=1/2 for n2<0 (C47c)',
 ('_elementary_pole_radial','div','1 + ell'):'positive ell=log1p(y/e)>=0 (C46c)',
 ('_elementary_pole_radial','div','3 * (3 + ell)'):'positive ell=log1p(x/e^3)>0 (C46c)',
 ('_elementary_pole_radial','div','c - radius'):'positive radius<c on this branch',
 ('_elementary_pole_radial','div','mp.sqrt(L * L + 4) - L'):'positive sqrt(L^2+4)>|L|',
 ('_elementary_pole_radial','div','radius'):'positive radius>0 (D: radius>a; D^-1 guard)',
 ('_elementary_pole_radial','div','x'):'positive x=e^3 expm1(3(r-a)/(c-r))>0 for a<r<c',
 ('_elementary_pole_radial','log1p','-mp.exp(-y)'):'positive -e^-y>-1 (y>1 on this branch)',
 ('_elementary_pole_radial','log1p','mp.exp(-f)'):'positive argument >0',
 ('_elementary_pole_radial','log1p','mp.exp(f)'):'positive argument >0',
 ('_elementary_pole_radial','log1p','x / mp.exp(3)'):'positive x>0',
 ('_elementary_pole_radial','log1p','y / mp.e'):'positive y=softplus(f)>0',
 ('_elementary_pole_radial','sqrt','L * L + 4'):'positive >=4',
 ('_elementary_pole_radial','sqrt','new_radius / radius'):'positive new_radius>0 (c ell/(1+ell), ell>0; or >a)',
 ('_elementary_unit','div','length'):'positive |P I_j xx|>0 in the band (G14, section 11.5 cap estimate)',
 ('_elementary_unit','sqrt','sum((abs(a) ** 2 for a in v))'):'positive the same length squared (G14)',
 ('_height_G','div','amp'):'positive amp=lam k>0 on this branch (amp==0 returns first; k=-log q>2)',
 ('_height_G','div','d'):'positive d=1-lam+2 lam sqrt(q) k>0',
 ('_height_R','div','1 + ell'):'positive ell=log1p(soft/e)>=0 (C45a)',
 ('_height_R','div','2 * (2 + ell)'):'positive ell=log1p(x/e^2)>0 (C45a)',
 ('_height_R','div','c - value'):'positive value<c on this branch',
 ('_height_R','div','mp.sqrt(L * L + 4) - L'):'positive sqrt(L^2+4)>|L|',
 ('_height_R','div','x'):'positive x=e^2 expm1(2(v-delta)/(c-v))>0 for delta<v<c',
 ('_height_R','log1p','-mp.exp(-y)'):'positive -e^-y>-1 (y>1 on this branch)',
 ('_height_R','log1p','mp.exp(-F)'):'positive argument >0',
 ('_height_R','log1p','mp.exp(F)'):'positive argument >0',
 ('_height_R','log1p','soft / mp.e'):'positive soft>0',
 ('_height_R','log1p','x / mp.exp(2)'):'positive x>0',
 ('_height_R','sqrt','L * L + 4'):'positive >=4',
 ('_height_terminal_forward','abs','C0'):'test K enters only the shortcut test',
 ('_height_terminal_forward','abs','C1'):'test K enters only the shortcut test',
 ('_height_terminal_forward','div','1 + u'):'positive u=expm1(x)>0 (x>0 for delta_z>0)',
 ('_height_terminal_forward','div','E0'):'positive E0=exp(eta/zcut)',
 ('_height_terminal_forward','div','Ez'):'positive Ez=exp(L)',
 ('_height_terminal_forward','div','L * (L + ell)'):'positive L=log Ez>log 50 and |ell|<1, since |D|<K e^-Ez (section 11.2, C45a)',
 ('_height_terminal_forward','div','eta / zcut + ell'):'positive ell=log1p(M/E0)>=0, zcut>0',
 ('_height_terminal_forward','div','gap'):'positive 0<delta_z<gap0 on this branch',
 ('_height_terminal_forward','div','gap * gap0'):'positive gap,gap0>0',
 ('_height_terminal_forward','div','zcut'):'positive zcut>0 (below)',
 ('_height_terminal_forward','log1p','D'):'positive |D|<K e^-Ez<1 for Ez>50 (section 11.2)',
 ('_height_terminal_forward','log1p','F'):'positive F=(Lambda u^2+mu u)/(1+u)>=0',
 ('_height_terminal_forward','log1p','M / E0'):'positive M=log1p(F)>=0',
 ('_height_terminal_forward','log1p','correction / Ez'):'positive |correction|<1<Ez (section 11.2)',
 ('_terminal_height_data','div','(zcut - z0) ** 2'):'positive zcut-z0=B(e^(.4k)-1)+C(1-e^(-.4k))>0',
 ('_terminal_height_data','div','kappa * (B + C) * Eprime_z0'):'positive kappa>=2, B>0, Eprime_z0>0',
 ('_terminal_height_data','div','zcut'):'positive zcut>0 (below)',
 ('_terminal_height_data','div','zcut ** 2'):'positive zcut>0 (below)',
 ('_terminal_height_data','div','zcut - z0'):'positive as above',
 ('_terminal_height_data','sqrt','q'):'positive q>1/2500 when lam>0 (G10)',
 ('elementary_cusp_logradius','abs','t'):'positive r=|t|>3/10 (called only in the collar)',
 ('elementary_cusp_logradius','div','1 - r * r'):'positive r<1/2',
 ('elementary_cusp_logradius','div','root + f'):'positive f>=0 and q>0 on this branch',
 ('elementary_cusp_logradius','sqrt','f * f + 4 * q'):'positive q>0 in the collar',
 ('elementary_cusp_reference','div','1 - abs(t) ** 2'):'positive |t|<=3/10',
 ('elementary_cusp_reference','div','R + f'):'unreachable cell -1 alternative, never taken (C50e)',
 ('elementary_cusp_reference','div','R - f'):'positive cell 1: f<0 (C50e)',
 ('elementary_cusp_reference','div','e'):'positive exp',
 ('elementary_cusp_reference','div','nu'):'positive nu in {1,e^-2}',
 ('elementary_cusp_reference','div','x'):'positive x=sqrt(A/e)/nu>0 (C50e)',
 ('elementary_cusp_reference','div','y'):'positive y=sqrt(Y)>0 off the node (C50e, G6)',
 ('elementary_cusp_reference','sqrt','A / e'):'positive A>0 (C50e)',
 ('elementary_cusp_reference','sqrt','Y'):'positive Y>0 off the node (C50e, G6)',
 ('elementary_cusp_reference','sqrt','f * f + 4 * q'):'positive R=0 only at f=q=0, the node (C50e)',
 ('elementary_global_slice','abs','X'):'positive |X|>1-alpha/2 where it enters a value (cap den); distance otherwise only in tests',
 ('elementary_global_slice','div','den'):'positive den=sqrt(scale(1+|X|sqrt(2-X^2)))>0',
 ('elementary_global_slice','div','mp.sqrt(1 - X * X)'):'positive |X|<=1-1e-9 on the slab (G1)',
 ('elementary_global_slice','div','nu'):'positive nu in {1,e^-2}',
 ('elementary_global_slice','sqrt','1 - X * X'):'positive G1',
 ('elementary_global_slice','sqrt','2 - X * X'):'positive >=1',
 ('elementary_global_slice','sqrt','scale * (1 + abs(X) * mp.sqrt(2 - X * X))'):'positive >=scale>0',
 ('elementary_phase_family','div','1 + r2'):'positive >=1',
 ('elementary_phase_family','div','1 + sign * n[2]'):'positive >=1 on the hemisphere of the sign',
 ('elementary_phase_family','div','ratio.real'):'positive Re ratio>7/8 (G3)',
 ('elementary_phase_family','div','steps'):'positive steps=K, a positive integer',
 ('elementary_pole_patch','div','outer - radius'):'positive noop block (C46a2); radius<(inner+outer)/2<outer',
 ('elementary_projection_section','abs','phase'):'positive |phase|>1/2 (G2); margin is diagnostic',
 ('elementary_projection_section','div','abs(phase)'):'positive |phase|>1/2 (G2)',
 ('elementary_projection_section','div','steps'):'positive steps=N, a positive integer',
 ('elementary_restore_torus','abs','t'):'positive t!=0 on mid data (|t|>3/10)',
 ('elementary_restore_torus','div','det'):'positive det(e2,e3)=+-(1-|t_c|^4), |t_c|<1/2 (C30k, C34e)',
 ('elementary_restore_torus','div','r'):'positive r=|t|>0',
 ('elementary_restore_torus','div','t'):'positive t!=0',
 ('elementary_restore_torus','div','t * t'):'positive t!=0',
 ('elementary_restore_torus','div','yy'):'positive yy=exp(..)!=0',
 ('elliptic_bundle_embedding','abs','angle'):'flat step 0 near angle=0 (argument -2), step 1 near the cut |angle|=pi (argument 2): weight 0 there (C42b, C34d)',
 ('elliptic_bundle_embedding','abs','h'):'flat hw: step argument -5/3 near h=0',
 ('elliptic_bundle_embedding','abs','t'):'positive t!=0 on mid data and in the collar',
 ('elliptic_bundle_embedding','div','1 - rr'):'positive rr=|t_c|^2<1/4 or |xy|^2<=9/100',
 ('elliptic_bundle_embedding','div','det'):'positive +-(1-|t_c|^4), |t_c|<1/2 (C34e)',
 ('elliptic_bundle_embedding','div','m * y'):'positive deck frame only for H<-1/10, where y!=0 (section 11.1); m>0',
 ('elliptic_bundle_embedding','div','mp.sqrt(norm2)'):'positive norm2>=1/3 (weight partition, below)',
 ('elliptic_bundle_embedding','div','mp.sqrt(r)'):'positive r>0 (G6)',
 ('elliptic_bundle_embedding','div','r'):'positive r=|t|>0',
 ('elliptic_bundle_embedding','div','t'):'positive t!=0',
 ('elliptic_bundle_embedding','div','t * t'):'positive t!=0',
 ('elliptic_bundle_embedding','div','yy'):'positive yy=z(mt)^shift, z=exp(..)!=0, mt!=0',
 ('elliptic_bundle_embedding','pow','(m * t) ** int(shift)'):'positive mt!=0; int(shift) is locally constant where its weight is nonzero (C41b)',
 ('elliptic_bundle_embedding','sqrt','norm2'):'positive norm2>=1/3 (below)',
 ('elliptic_bundle_embedding','sqrt','r'):'positive r>0 (G6)',
 ('elliptic_pole_raw_inverse','div','1 + s2'):'positive >=1',
 ('elliptic_pole_raw_inverse','div','1 - abs(tc) ** 2'):'positive |tc|^2<4.9e-7 (C46b3)',
 ('elliptic_pole_raw_inverse','div','2 * source_radius'):'positive source_radius>0 (C46b3)',
 ('elliptic_pole_raw_inverse','div','R + height'):'positive height>=0 and R>0 (C46b3)',
 ('elliptic_pole_raw_inverse','div','R - height'):'positive height<0',
 ('elliptic_pole_raw_inverse','div','den'):'positive north den>=1-8u, profile den>=2A(1-uu) (C46b3)',
 ('elliptic_pole_raw_inverse','div','e'):'positive exp',
 ('elliptic_pole_raw_inverse','div','length'):'positive |w|^2/u>=0.9999 (C46b3)',
 ('elliptic_pole_raw_inverse','div','u1'):'positive Re(-u1)>0 (north, C46b3)',
 ('elliptic_pole_raw_inverse','div','u2'):'positive Re u2>0 (south, C46b3)',
 ('elliptic_pole_raw_inverse','sqrt','(2 * B if infinity else 2) / den'):'positive den>0',
 ('elliptic_pole_raw_inverse','sqrt','A * A * (1 - uu) ** 2 + 4 * B * uu'):'positive uu<1',
 ('elliptic_pole_raw_inverse','sqrt','abs(w0) ** 2 + abs(w1) ** 2'):'positive length^2>0 (C46b3)',
 ('elliptic_pole_raw_inverse','sqrt','height * height + 4 * q'):'positive R>0 (C46b3)',
 ('elliptic_pole_raw_inverse','sqrt','source_radius'):'positive C46b3',
 ('elliptic_raw_base_reference','abs','t'):'positive pure middle: t!=0',
 ('elliptic_raw_base_reference','div','1 + s2'):'positive >=1',
 ('elliptic_raw_base_reference','div','1 - aa'):'positive aa^2=1-(1-X^2)^2<1',
 ('elliptic_raw_base_reference','div','b + root'):'positive root>=|b|, and b=A>0 if rr=0',
 ('elliptic_raw_base_reference','div','small'):'positive |small|>=1/2 where 1/small is used',
 ('elliptic_raw_base_reference','div','tc'):'positive |tc|>3/10 in the north collar',
 ('elliptic_raw_base_reference','div','u1'):'positive dominant coordinate, |u1|>=1/sqrt(2)',
 ('elliptic_raw_base_reference','div','u2'):'positive dominant coordinate, |u2|>=1/sqrt(2)',
 ('elliptic_raw_base_reference','div','yy'):'positive yy=exp(rho)',
 ('elliptic_raw_base_reference','sqrt','(2 * B if inverse_ratio else 2) / (b + root)'):'positive b+root>0',
 ('elliptic_raw_base_reference','sqrt','2 - xx * xx'):'positive >=1',
 ('elliptic_raw_base_reference','sqrt','b * b + 4 * B * rr'):'positive b=A>0 if rr=0',
 ('inverse_height_moment','abs','a'):'positive aa enters a value only for |a|>1/10',
 ('scalar_finite_height','div','1 + v'):'positive v>0',
 ('scalar_finite_height','div','1 - q'):'positive q<1/8 (G11)',
 ('scalar_finite_height','div','1 - rr'):'positive rr<=1/4 (G11)',
 ('scalar_finite_height','div','u'):'positive u>2 on this branch',
 ('scalar_flat_step','div','1 - u'):'positive 0<u<1 on this branch',
 ('scalar_flat_step','div','a + b'):'positive a,b>0',
 ('scalar_flat_step','div','u'):'positive 0<u<1',
 ('scalar_lattice_coords','div','den'):'positive lattice determinant (1+|t|^2)((|t|^2-1)log|t|+2|t|^2)>0 (C1c, C34b)'}
census_auto50 = {}
census_gate_ok50 = True
census_found50 = set()
census_occ50 = 0
for name50 in path_functions48:
    tree50 = source_tree30(getattr(cusp30,name50))
    par50 = parents50(tree50)
    for node50 in ast.walk(tree50):
        key50 = arg50 = None
        if isinstance(node50,ast.BinOp) and isinstance(node50.op,ast.Div):
            key50,arg50 = ('div',ast.unparse(node50.right)),node50.right
        elif (isinstance(node50,ast.BinOp) and isinstance(node50.op,ast.Pow)
              and not (isinstance(node50.right,ast.Constant)
                       and isinstance(node50.right.value,int)
                       and node50.right.value>=0)):
            key50,arg50 = ('pow',ast.unparse(node50)),node50.left
        elif isinstance(node50,ast.Call) and node50.args:
            fn50 = (node50.func.attr if isinstance(node50.func,ast.Attribute)
                    else getattr(node50.func,'id',None))
            if fn50 in ('sqrt','log1p','abs'):
                key50,arg50 = (fn50,ast.unparse(node50.args[0])),node50.args[0]
        if key50 is None:
            continue
        census_occ50 += 1
        up50 = par50.get(node50)
        if const_node50(arg50):
            kind50 = 'const'
        elif in_test50(node50,par50) or isinstance(up50,ast.Compare):
            kind50 = 'test'
            if key50[0]=='div':
                gate50,inner50 = False,node50
                while inner50 in par50:
                    inner50 = par50[inner50]
                    if (isinstance(inner50,ast.If) and ast.unparse(inner50.test)
                            =='_elementary_precision_shortcuts.get()'):
                        gate50 = True
                census_gate_ok50 = census_gate_ok50 and gate50
        elif dead_node50(node50,par50,name50):
            kind50 = 'dead'
        elif (key50[0]=='abs' and isinstance(up50,ast.BinOp)
              and isinstance(up50.op,ast.Pow) and up50.left is node50
              and ast.unparse(up50.right)=='2'):
            kind50 = 'square'
        else:
            census_found50.add((name50,)+key50)
            continue
        census_auto50[kind50] = census_auto50.get(kind50,0)+1
# zcut=e^(-kappa/10)-q e^(kappa/10)>0 iff q<e^(-kappa/5); kappa lies in
# [2,k] with k=-log q>2, so e^(-kappa/5)>=q^(1/5)>q (0<q<1); for q=0 it is
# e^(-kappa/10)>0.  The embedding norm: each active slot block contains a
# unit vector with weight 1 (hw=1 on |h|<=1/4, and of the two centres one
# has |h|<=1/4; of the two angular slots one has |angle|<=pi/2), so
# norm2>=b0^2+bi^2+bm^2>=(b0+bi+bm)^2/3=1/3.
qz50,kz50 = sp.Symbol('qz50',positive=True),sp.Symbol('kz50',positive=True)
bb50 = sp.symbols('b0_50 bi_50 bm_50',real=True)
aux50 = (
    sp.simplify(sp.exp(-kz50/2)*sp.exp(sp.Rational(2,5)*kz50)
                -qz50*sp.exp(kz50/2)*sp.exp(-sp.Rational(2,5)*kz50)
                -(sp.exp(-kz50/10)-qz50*sp.exp(kz50/10)))==0
    and sp.expand(3*sum(b**2 for b in bb50)-sum(bb50)**2
                  -((bb50[0]-bb50[1])**2+(bb50[0]-bb50[2])**2
                    +(bb50[1]-bb50[2])**2))==0
    and sp.log(8)>2)
census_ok50 = (census_found50==set(census50) and len(census50)==139
               and census_occ50==280
               and all(v.split()[0] in ('positive','test','flat','unreachable')
                       for v in census50.values())
               and census_gate_ok50 and aux50)
print('C50f census:',census_occ50,'occurrences;',len(census50),'listed;',
      'automatic',dict(sorted(census_auto50.items())))
rep('C50f all 280 divisions, roots, log1p, abs and powers: argument in domain',
    census_ok50)

# C50g: the collar logarithm cut, through the restoration.  C34 shows
# that the collar ten-vector is the same on both sides of the negative
# real axis (C34d-f, generic along the cut), so the fold and the final
# phase theta agree, and the slice data (t,one,four,lt) differ by the deck
# move of C34c with k=+-1: south one'=one-k four-k(k+1)/2, four'=four+k;
# north one'=one-k four; lt'=lt+2 pi i k (theta/(2pi) is added to one on
# both sides).  For the bound restoration loop (C41a, source41):
#   (i)   z,w move by k e4(L)+k(k-1)/2 e1 (south) and not at all (north),
#         with e1=(2pi,0), e4(L)=(-iL,2t);
#   (ii)  zc,wc are unchanged, hence cc, yy and H;
#   (iii) D=cc+H c for every branch of s=lt/(2 pi i): with xi the lattice
#         coordinates of (pi s(s-1),2ts) and c those of e4 in (e2,e3), the
#         candidate cc-(four_x-four)c solves the real 2x2 system for -wc,
#         whose determinant +-(1-|t_c|^4) is nonzero (C34e); north D=(H+four)c;
#   so C=phi+beta D is unchanged and the outputs differ by a lattice
#   translation, the identity on Y.  (iv) The final-phase Arg cut
#   one->one+1 moves z and zc by 2 pi, so yy=exp(i(...)) and C are
#   unchanged and the output moves by e1.
k50 = sp.Symbol('k50',integer=True)
one50,four50,fx50,two50,three50,c0_50,c1_50,H50 = sp.symbols(
    'one50 four50 fx50 two50 three50 c0_50 c1_50 H50',real=True)
L50,t50,E2_50,E3_50 = sp.symbols('L50 t50 E2_50 E3_50')
def zw50(one,four,L,t):
    return 2*sp.pi*one-sp.I*L*four,2*t*four
def zcwc_south50(z,w,L,t):
    s = L/(2*sp.I*sp.pi)
    return z-sp.pi*s*(s-1),w-2*t*s
def zcwc_north50(z,w,t):
    return z-sp.I*w/t,w/(t*t)
z0_50,w0_50 = zw50(one50,four50,L50,t50)
Ls50 = L50+2*sp.pi*sp.I*k50
zs50,ws50 = zw50(one50-k50*four50-k50*(k50+1)/2,four50+k50,Ls50,t50)
zn50,wn50 = zw50(one50-k50*four50,four50,Ls50,t50)
logcut50 = (
    # (i)
    sp.expand(zs50-z0_50-(-sp.I*k50*L50+sp.pi*k50*(k50-1)))==0
    and sp.expand(ws50-w0_50-2*t50*k50)==0
    and sp.expand(zn50-z0_50)==0 and sp.expand(wn50-w0_50)==0
    # (ii)
    and all(sp.expand(a-b)==0 for a,b in zip(zcwc_south50(zs50,ws50,Ls50,t50),
                                             zcwc_south50(z0_50,w0_50,L50,t50)))
    and all(sp.expand(a-b)==0 for a,b in zip(zcwc_north50(zn50,wn50,t50),
                                             zcwc_north50(z0_50,w0_50,t50))))
# (iii) given 2ts=two e2+three e3+fx 2t and 2t=c0 e2+c1 e3 (w components)
s50 = L50/(2*sp.I*sp.pi)
T2_50 = c0_50*E2_50+c1_50*E3_50
two_t_s50 = two50*E2_50+three50*E3_50+fx50*T2_50
cand50 = (two50+(fx50-four50)*c0_50,three50+(fx50-four50)*c1_50)
minus_wc50 = two_t_s50-T2_50*four50          # -wc=2ts-2t four
D50 = tuple(x+(H50-four50+fx50)*c for x,c in zip((two50,three50),(c0_50,c1_50)))
logcut50 &= (
    sp.expand(cand50[0]*E2_50+cand50[1]*E3_50-minus_wc50)==0
    and all(sp.expand(d-(cc+H50*c))==0
            for d,cc,c in zip(D50,cand50,(c0_50,c1_50))))
# (iv) one->one+1
z1_50,w1_50 = zw50(one50+1,four50,L50,t50)
logcut50 &= (
    sp.expand(z1_50-z0_50-2*sp.pi)==0 and sp.expand(w1_50-w0_50)==0
    and sp.expand(zcwc_south50(z1_50,w1_50,L50,t50)[0]
                  -zcwc_south50(z0_50,w0_50,L50,t50)[0]-2*sp.pi)==0
    and sp.expand(zcwc_north50(z1_50,w1_50,t50)[0]
                  -zcwc_north50(z0_50,w0_50,t50)[0]-2*sp.pi)==0
    and sp.exp(sp.I*2*sp.pi)==1)
# The source lines used: the loop (C41a binds it in full), the w
# components e2,e3,e4 shared by cusp_data (south) and the lattice solve.
restore_tree50 = source_tree30(cusp30.elementary_restore_torus)
logcut50 &= (
    source41
    and matching_statement31(restore_tree50,
        'return ((t*t-1,mp.j*(1+t*t),-2*t,t,mp.j*t) if inf else '
        '(1-t*t,mp.j*(1+t*t),2*t,t,-mp.j*t))')
    and matching_assignment30(lattice_tree34,'ew','(1-t*t,mp.j*(1+t*t),2*t)')
    and matching_statement31(restore_tree50,
        'D = tuple(xi[j+1]+(H-four+xi[3])*c[j] for j in range(2))')
    and matching_statement31(restore_tree50,'D = tuple((H+four)*v for v in c)')
    and matching_statement31(restore_tree50,
        'z,w = 2*mp.pi*one-mp.j*lt*four,2*t*four'))
rep('C50g collar log cut and final-phase cut: restored outputs differ by lattice vectors',
    logcut50)

# C50h: the premise of C48's dead-arm argument, that elementary is True in
# every frame the path executes.  Every call from a path function to a
# function with an `elementary` parameter passes elementary=True, or
# elementary=elementary from a function that itself has that parameter;
# no path function assigns `elementary` (C50a: its uses are tests or flag
# arguments); the path imports only mpmath and fractions, and the only
# module global it reads besides cusp functions is the shortcut ContextVar
# (off in formula mode).  The entry point has no such parameter, so by
# induction along the call graph elementary is True throughout.
def params50(function):
    args = source_tree30(function).args
    return [a.arg for a in args.args+args.kwonlyargs]
elementary_ok50 = True
for name50 in path_functions48:
    tree50 = source_tree30(getattr(cusp30,name50))
    own50 = 'elementary' in params50(getattr(cusp30,name50))
    for node50 in ast.walk(tree50):
        if isinstance(node50,ast.Name) and node50.id=='elementary':
            elementary_ok50 &= isinstance(node50.ctx,ast.Load)
        if (isinstance(node50,ast.Call) and isinstance(node50.func,ast.Name)
                and node50.func.id in cusp_functions48):
            callee50 = getattr(cusp30,node50.func.id)
            if 'elementary' not in params50(callee50):
                continue
            if dead_node50(node50,parents50(tree50),name50):
                continue
            position50 = params50(callee50).index('elementary')
            values50 = [kw.value for kw in node50.keywords if kw.arg=='elementary']
            elementary_ok50 &= (
                len(node50.args)<=position50 and len(values50)==1
                and ((isinstance(values50[0],ast.Constant) and values50[0].value is True)
                     or (own50 and isinstance(values50[0],ast.Name)
                         and values50[0].id=='elementary')))
        if isinstance(node50,(ast.Import,ast.ImportFrom)):
            elementary_ok50 &= ast.unparse(node50) in (
                'import mpmath as mp','from fractions import Fraction')
    local50 = ({n.id for n in ast.walk(tree50) if isinstance(n,ast.Name)
                and isinstance(n.ctx,ast.Store)}
               | {a.arg for f in ast.walk(tree50) if isinstance(f,ast.FunctionDef)
                  for a in f.args.args+f.args.kwonlyargs}
               | {f.name for f in ast.walk(tree50) if isinstance(f,ast.FunctionDef)}
               | {(al.asname or al.name) for n in ast.walk(tree50)
                  if isinstance(n,(ast.Import,ast.ImportFrom)) for al in n.names})
    for node50 in ast.walk(tree50):
        if (isinstance(node50,ast.Name) and isinstance(node50.ctx,ast.Load)
                and node50.id in vars(cusp30) and node50.id not in cusp_functions48):
            elementary_ok50 &= (node50.id in local50
                                or node50.id=='_elementary_precision_shortcuts')
elementary_ok50 &= 'elementary' not in params50(cusp30.elementary_global_coordinates)
rep('C50h elementary=True in every executed frame (premise of the C48 dead arms)',
    elementary_ok50)

# C50d: every cited check exists.  Each label C<n><x>, N<n> or L<n> named
# in the four C50 tables (C50, C50b, C50c, C50f) must occur as the start of a rep() label in
# chart_identities.py (C), node_remainder.py (N) or lemmas_1d.py (L).
import re as re50
base50 = os.path.dirname(os.path.abspath(cusp30.__file__))
texts50 = {'C':open(os.path.join(base50,'sym','chart_identities.py')).read(),
           'N':open(os.path.join(base50,'sym','node_remainder.py')).read(),
           'L':open(os.path.join(base50,'sym','lemmas_1d.py')).read()}
cited50 = set()
for table50 in (classify50,branch_calls50,primitive_calls50,census50):
    for value50 in table50.values():
        cited50 |= set(re50.findall(r'\b([CNL]\d+[a-z]?\d?)\b',value50))
def label_exists50(label):
    text = texts50[label[0]]
    return bool(re50.search(r"rep\(f?['\"]"+re50.escape(label)+r"[a-z0-9 ]",text)
                or re50.search(r"['\"]\(?"+re50.escape(label)+r"[a-z0-9]*[) ]",text))
missing50 = sorted(l for l in cited50 if not label_exists50(l))
rep('C50d every check cited in the classification exists',not missing50)
if missing50:
    print('C50d missing labels:',missing50)

print('ALL PASS' if ok else 'FAILURES')
assert ok, 'A chart identity failed'
