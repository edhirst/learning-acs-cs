"""The two cusps of Y (paper: Prop. 'Mumford-Tate filling'), on the T^2-slice w = 0 of the filling.

Cusp data (paper, Sections 'The cusp at t = 0' and Lemma 'atinfty'), in the base coordinate t (resp. s = 1/t):
    CUSP0:   alpha2 = t,  alpha3 = -i t,  e2bar = 1 - t^2,  e3bar = i(1 + t^2),  e4bar =  2t,  nu = 1
    CUSPINF: alpha2 = s,  alpha3 =  i s,  e2bar = s^2 - 1,  e3bar = i(1 + s^2),  e4bar = -2s,  nu = e^-2
The T^2 element c multiplies the X-factor by exp(c2 alpha2 + c3 alpha3). On the slice, G4 followed by the
T^2-correction restoring w = 0 is
    psi_hat(zeta, t) = (nu t exp(-(c2 alpha2 + c3 alpha3)) zeta, t),   c = R(t)^{-1} e4bar(t),
and log|zeta| o psi_hat = log|zeta| + P(t), with the paper's height P = log|nu| + log|t| - phi_t(e4bar).

Height (docs/s4t2_explicit_diffeo.md, 5a):  mu = sum_n n w_n,  w_n ~ exp(-2 n log|zeta| + n (n + 1) P),
so mu o psi_hat = mu + 1. In the toric chart U_k (zeta = x^k y^(k+1), t = x y), with j = n - k,
    w_n ~ |x|^(j (j + 1)) |y|^(j (j - 1)) exp((2 k j + j^2 + j) (log|nu| - phi(t))),
smooth across the central fibre; mu = k at the node P_k.

Weight-1 function (notes, section 7):  E = sum_n chi(mu + n) E0(psi_hat^n p),  E0 = y0 + conj(x0) in U_0.
"""
from __future__ import annotations

import numpy as np


def bump_step(u):
    """C-infinity step: 0 for u <= 0, 1 for u >= 1."""
    u = np.clip(u, 0.0, 1.0)
    a = np.where(u > 0, np.exp(-1.0 / np.maximum(u, 1e-300)), 0.0)
    b = np.where(u < 1, np.exp(-1.0 / np.maximum(1 - u, 1e-300)), 0.0)
    return a / (a + b)


def chi(m, inner=0.4, outer=0.9):
    """Cutoff in the height: 1 on |m| <= inner, 0 on |m| >= outer (so sum_n chi(m + n) > 0)."""
    return bump_step((outer - np.abs(m)) / (outer - inner))


def real_coeffs(u, a, b):
    """Real (c2, c3) with c2 a + c3 b = u (a, b R-independent)."""
    det = (a.conj() * b).imag
    return (u.conj() * b).imag / det, (a.conj() * u).imag / det


class Cusp:
    def __init__(self, name, al2, al3, e2, e3, e4, nu, e0coef=1.0):
        self.name, self.al2, self.al3, self.e2, self.e3, self.e4, self.nu = name, al2, al3, e2, e3, e4, nu
        self.e0coef = e0coef        # local model E0 = y0 + e0coef * conj(x0); zero set {mu_node = 0, arg t = pi}

    def _c(self, t):
        t = np.asarray(t, complex)
        return real_coeffs(self.e4(t), self.e2(t), self.e3(t))

    def phi(self, t):
        """phi_t(e4bar), phi_t the R-linear functional with phi_t(e_j bar) = Re alpha_j (j = 2, 3)."""
        t = np.asarray(t, complex)
        c2, c3 = self._c(t)
        return c2 * self.al2(t).real + c3 * self.al3(t).real

    def P(self, t):
        t = np.asarray(t, complex)
        return np.log(abs(self.nu)) + np.log(np.abs(t)) - self.phi(t)

    def multiplier(self, t):
        """psi_hat(zeta, t) = (multiplier(t) zeta, t)."""
        t = np.asarray(t, complex)
        c2, c3 = self._c(t)
        return self.nu * t * np.exp(-(c2 * self.al2(t) + c3 * self.al3(t)))

    def psi_hat(self, zeta, t):
        return self.multiplier(t) * zeta, t

    def mu_torus(self, zeta, t, J=12):
        rho = np.log(np.abs(zeta))
        p = self.P(t)
        nstar = np.round(rho / p - 0.5)
        n = nstar[..., None] + np.arange(-J, J + 1)
        ex = -2 * n * rho[..., None] + n * (n + 1) * p[..., None]
        w = np.exp(ex - ex.max(-1, keepdims=True))
        return (n * w).sum(-1) / w.sum(-1)

    def mu_chart(self, x, y, k, J=8):
        t = x * y
        g = np.log(abs(self.nu)) - self.phi(t)
        X, Y = np.abs(x) ** 2, np.abs(y) ** 2
        num, den = np.zeros_like(X), np.zeros_like(X)
        for j in range(-J, J + 1):
            term = X ** (j * (j + 1) // 2) * Y ** (j * (j - 1) // 2) * np.exp((2 * k * j + j * j + j) * g)
            num, den = num + j * term, den + term
        return k + num / den

    def E0_torus(self, zeta, t):
        """y0 + e0coef conj(x0) in U_0 (y0 = zeta, x0 = t / zeta). With e0coef = |nu| the zero set
        {|y0| = |nu| |x0|, arg t = pi} is exactly where the node height a0|x0|^2 - |y0|^2 vanishes."""
        return zeta * (1 + self.e0coef * np.conj(t) / np.abs(zeta) ** 2)

    def E_torus(self, zeta, t):
        mu = self.mu_torus(zeta, t)
        mt = self.multiplier(t)
        out = np.zeros_like(zeta)
        n0 = -np.round(mu).astype(int)
        for dn in (-1, 0, 1):
            n = n0 + dn
            out = out + chi(mu + n) * self.E0_torus(mt ** n * zeta, t)
        return out


CUSP0 = Cusp('0', lambda t: t, lambda t: -1j * t, lambda t: 1 - t ** 2, lambda t: 1j * (1 + t ** 2),
             lambda t: 2 * t, 1.0)
CUSPINF = Cusp('inf', lambda s: s, lambda s: 1j * s, lambda s: s ** 2 - 1, lambda s: 1j * (1 + s ** 2),
               lambda s: -2 * s, np.exp(-2.0), e0coef=np.exp(-2.0))


def elliptic_height_reparam(h,inverse=False,scale=None):
    """Degree-one real lift of atan(c*tan(pi*h))/pi, c=pi/2.

    mpmath scalar arithmetic; half-cell reduction makes the branch explicit.
    This normalization makes the nodal linear term a0|x|^2-|y|^2.
    """
    import mpmath as mp
    h = mp.mpf(h)
    n = mp.floor(h+mp.mpf('.5'))
    delta = h-n
    c = mp.pi/2 if scale is None else mp.mpf(scale)
    if inverse:
        c = 1/c
    return n+mp.atan2(c*mp.sin(mp.pi*delta),mp.cos(mp.pi*delta))/mp.pi


def elliptic_height(rho,P,normalized=True,dps=50,blend=1):
    """Closed special-function cusp height, programme 10.5 (scalar reference).

    Uses Jacobi amplitude via sn/cn and an explicitly reduced real branch.
    P<0 and t!=0; central axes use elliptic_height_chart. These independent
    helpers do NOT change Cusp.mu_* or the legacy numerical construction.
    """
    import mpmath as mp
    P_check,blend_check = mp.mpf(P),mp.mpf(blend)
    if P_check >= 0:
        raise ValueError('elliptic_height requires P<0')
    if not 0 <= blend_check <= 1:
        raise ValueError('elliptic height blend must lie in [0,1]')
    # Resolve 1-k^2 when the complementary nome is extremely small.
    precision = max(mp.mp.dps,dps)+max(0,int(-2*P_check/mp.log(10)))+10
    with mp.workdps(precision):
        # Convert inside the guarded context: mp.mpf rounds an existing mpf.
        # Premature rounding loses inverse accuracy on nearly flat node fibres.
        rho,P,blend = mp.mpf(rho),mp.mpf(P),mp.mpf(blend)
        q = mp.exp(2*P)
        B = mp.jtheta(3,0,q)**2
        kp = mp.jtheta(2,0,q)**2/B
        parameter = blend*blend*(1-kp*kp)
        if blend != 1:
            B = -mp.ellipk(parameter)/P
        n = mp.floor(rho/P)
        u = -2*B*(rho-n*P)
        amplitude = mp.atan2(mp.ellipfun('sn',u,parameter),mp.ellipfun('cn',u,parameter))
        h = n+amplitude/mp.pi-mp.mpf('.5')
        scale = 1+(mp.pi/2-1)*blend
        return +elliptic_height_reparam(h,scale=scale) if normalized else +h


def elliptic_height_inverse(a,P,normalized=True,dps=50,blend=1):
    """Closed inverse rho(a,P), using standard incomplete elliptic F.

    No root finder or user-defined inverse is involved. F uses m=k^2 in
    mpmath, whereas the document/DLMF use the modulus k.
    """
    import mpmath as mp
    P_check,blend_check = mp.mpf(P),mp.mpf(blend)
    if P_check >= 0:
        raise ValueError('elliptic_height_inverse requires P<0')
    if not 0 <= blend_check <= 1:
        raise ValueError('elliptic height blend must lie in [0,1]')
    precision = max(mp.mp.dps,dps)+max(0,int(-2*P_check/mp.log(10)))+10
    with mp.workdps(precision):
        a,P,blend = mp.mpf(a),mp.mpf(P),mp.mpf(blend)
        q = mp.exp(2*P)
        B = mp.jtheta(3,0,q)**2
        kp = mp.jtheta(2,0,q)**2/B
        parameter = blend*blend*(1-kp*kp)
        if blend != 1:
            B = -mp.ellipk(parameter)/P
        scale = 1+(mp.pi/2-1)*blend
        h = elliptic_height_reparam(a,inverse=True,scale=scale) if normalized else a
        n = mp.floor(h+mp.mpf('.5'))
        phi = mp.pi*(h-n+mp.mpf('.5'))
        return +(n*P-mp.ellipf(phi,parameter)/(2*B))


def elliptic_height_chart(x,y,nu=1,normalized=True,dps=50):
    """The same special-function height in U0, including both exact axes."""
    import mpmath as mp
    with mp.workdps(max(mp.mp.dps,dps)):
        x,y,nu = mp.mpc(x),mp.mpc(y),mp.mpf(nu)
        if x == 0 or y == 0:
            h = 2/mp.pi*(mp.atan(nu*nu*abs(x)**2)-mp.atan(abs(y)**2))
            return +elliptic_height_reparam(h) if normalized else +h
        r2 = abs(x*y)**2
        if r2 >= 1:
            raise ValueError('elliptic_height_chart requires |xy|<1')
        P = mp.log(nu)+mp.log(abs(x*y))-2*r2/(1-r2)
        return elliptic_height(mp.log(abs(y)),P,normalized,dps)


def elliptic_collar_blend(r):
    """Exact scalar flat cutoff: one on r<=.3, zero on r>=.5."""
    import mpmath as mp
    u = (mp.mpf(r)-mp.mpf('.3'))/mp.mpf('.2')
    if u <= 0:
        return mp.mpf(1)
    if u >= 1:
        return mp.mpf(0)
    a,b = mp.exp(-1/u),mp.exp(-1/(1-u))
    return b/(a+b)


def elliptic_axis_radius(a):
    """Closed central-axis inverse a=+/-R((2/pi)atan(rad)), |a|<1."""
    import mpmath as mp
    a = abs(mp.mpf(a))
    if a >= 1:
        raise ValueError('central U0 axis height must satisfy |a|<1')
    return mp.tan(mp.pi*elliptic_height_reparam(a,inverse=True)/2)


def elliptic_cusp_angle(rho,t,nu=1,dps=50):
    """Complete closed height/collar block, without the separate node correction."""
    import mpmath as mp
    with mp.workdps(max(mp.mp.dps,dps)):
        r = abs(mp.mpc(t))
        if not 0 < r <= mp.mpf('.5'):
            raise ValueError('elliptic_cusp_angle requires 0<|t|<=.5')
        P = mp.log(mp.mpf(nu))+mp.log(r)-2*r*r/(1-r*r)
        return elliptic_height(rho,P,dps=dps,blend=elliptic_collar_blend(r))


def elliptic_cusp_logradius(a,t,nu=1,dps=50):
    """Closed inverse of elliptic_cusp_angle, including the entire collar."""
    import mpmath as mp
    with mp.workdps(max(mp.mp.dps,dps)):
        r = abs(mp.mpc(t))
        if not 0 < r <= mp.mpf('.5'):
            raise ValueError('elliptic_cusp_logradius requires 0<|t|<=.5')
        P = mp.log(mp.mpf(nu))+mp.log(r)-2*r*r/(1-r*r)
        return elliptic_height_inverse(a,P,dps=dps,blend=elliptic_collar_blend(r))


def elliptic_pole_raw_inverse(z,infinity=False,dps=50,elementary=False):
    """Closed local inverse lift before smoothing at the pole (programme 11.3).

    z is the two-complex-coordinate linear pole model, not a sphere point.
    Return the slice pair (sqrt(a0)*conj(x),y). Defined on 0<|z|^2<=1e-4;
    no implicit root or iteration is used. Smoothness at z=0 is NOT asserted.
    """
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
        # Inverse O, and a stable inverse stereographic expression at each pole.
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
            u1 = mp.mpc(1-s2,2*aa)/(1+s2)
            u2 = mp.mpc(2*(bb+1),2*cc)/(1+s2)
            ratio = u1/u2
            height = mp.arg(u2)/(2*mp.pi)
        A,B = (4*mp.pi)**2,mp.exp(4)
        uu = abs(ratio)**2
        den = A*(1-uu)+mp.sqrt(A*A*(1-uu)**2+4*B*uu)
        tc = ratio*mp.sqrt((2*B if infinity else 2)/den)
        nu = mp.exp(-2) if infinity else mp.mpf(1)
        if elementary:
            # This cap lies in |a|<.1, where f(a,q)=a. Recover u and Y
            # directly; this is smooth on both nonnode central axes and does
            # not form the polar 0/0 quotient nu^2*|tc|^2/Y there.
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


def elliptic_pole_patch(z,infinity=False,dps=50,elementary=False):
    """Inverse-first closed pole patch: identity on |z|^2<=1e-8.

    The flat blend reaches elliptic_pole_raw_inverse at |z|^2>=2e-8.
    elementary=True is the cap used by the global map of programme 11.
    The default retains the earlier, separate elliptic design.
    """
    import mpmath as mp
    with mp.workdps(max(mp.mp.dps,dps)+20):
        zz = tuple(map(mp.mpc,z))
        radius = sum(abs(w)**2 for w in zz)
        u = radius/mp.mpf('1e-8')-1
        if u<=0:
            return zz
        raw = elliptic_pole_raw_inverse(zz,infinity,dps=dps,elementary=elementary)
        if u>=1:
            return raw
        a,b = mp.exp(-1/u),mp.exp(-1/(1-u))
        weight = a/(a+b)
        return tuple(+(w+weight*(v-w)) for w,v in zip(zz,raw))


def _elementary_pole_radial(pair,inverse=False):
    """The radial conjugation D (or D inverse) in programme (11.12g--h)."""
    import mpmath as mp
    z = tuple(map(mp.mpc,pair))
    radius = sum(abs(w)**2 for w in z)
    a,c = mp.mpf('1e-8'),mp.mpf('1.5e-8')
    if radius>=c:
        return z
    # In (11.12i), q=exp(-c/(c-radius)) and |S(radius)-radius|/radius
    # (likewise for S^{-1}) is <176q here. Avoid forming 1/q as an exp once
    # this is below working precision.
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


def elementary_pole_raw_forward(w,infinity=False,dps=50):
    """Finite elementary inverse of the raw local lift (11.10), away from zero."""
    import mpmath as mp
    with mp.workdps(max(mp.mp.dps,dps)+20):
        pair = tuple(map(mp.mpc,w))
        radius = sum(abs(a)**2 for a in pair)
        if radius<=0:
            raise ValueError('raw pole forward map requires nonzero spinor')
        nu = mp.exp(-2) if infinity else mp.mpf(1)
        tc = mp.conj(pair[0])*pair[1]/nu
        u,v = abs(pair[0])**2,abs(pair[1])**2
        e = mp.exp(-4*abs(tc)**2/(1-abs(tc)**2))
        height = e*u-v
        A,B = (4*mp.pi)**2,mp.exp(4)
        if infinity:
            ratio = tc*mp.sqrt((abs(tc)**2+A)/(A*abs(tc)**2+B))
            u1 = mp.exp(2*mp.pi*mp.j*(-height+mp.mpf('.5')))/mp.sqrt(1+abs(ratio)**2)
            u2 = ratio*u1
        else:
            ratio = tc*mp.sqrt((A+B*abs(tc)**2)/(1+A*abs(tc)**2))
            u2 = mp.exp(2*mp.pi*mp.j*height)/mp.sqrt(1+abs(ratio)**2)
            u1 = ratio*u2
        if infinity:
            # Rational inverted stereographic coordinates stay defined at
            # the exact north axis u1=-1, u2=0.
            den = 2*(1-u2.real)
            stereo = (u1.imag/den,(u2.real-1-u1.real)/den,u2.imag/den)
        else:
            den = 1+u1.real
            stereo = (u1.imag/den,u2.real/den-1,u2.imag/den)
        aa,bb,cc = stereo[0],-stereo[2],-stereo[1]
        scale = mp.pi if infinity else 2*mp.pi
        p = tuple(a/scale for a in (aa,bb,cc))
        # The local linear pole model is v=s_j O_j Hopf(z).
        pL = ((p[2],-p[1],p[0]) if infinity else (p[2],-p[0],-p[1]))
        rL = mp.sqrt(sum(a*a for a in pL))
        n = tuple(a/rL for a in pL)
        projected = ((1+n[2])*pair[0]+(n[0]+mp.j*n[1])*pair[1],
                     (n[0]-mp.j*n[1])*pair[0]+(1-n[2])*pair[1])
        length = mp.sqrt(sum(abs(a)**2 for a in projected))
        if length==0:
            raise ValueError('raw pole forward projector is orthogonal')
        factor = mp.sqrt(rL)/length
        return tuple(+(factor*a) for a in projected)


def elementary_pole_patch(z,infinity=False,dps=50):
    """Two-sided elementary pole patch H_j of programme (11.12h)."""
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
            if x-1/x<=-(dps+10)*mp.log(10):
                # u=S(radius)<(outer/e)*10^(-dps-10). For u<c/10,
                # J(t)=log(expm1(psi(t))) has J'(t)<15/(2t).
                # The raw bounds (11.9),(11.11) then give
                # |H(z)-z|/|z|<1000u, below the requested precision.
                return pair
        compressed = _elementary_pole_radial(pair)
        raw = elliptic_pole_raw_inverse(compressed,infinity,dps=dps,elementary=True)
        return _elementary_pole_radial(raw,inverse=True)


def elementary_pole_forward(w,infinity=False,dps=50):
    """Finite elementary inverse H_j^{-1} on the certified cap image."""
    import mpmath as mp
    with mp.workdps(max(mp.mp.dps,dps)+50):
        pair = tuple(map(mp.mpc,w))
        inner = mp.mpf('1e-8')
        outer = mp.mpf('1.5e-8')
        radius = sum(abs(a)**2 for a in pair)
        if radius<=inner:
            return pair
        if radius<(inner+outer)/2:
            x = mp.exp(3)*mp.expm1(3*(radius-inner)/(outer-radius))
            if x-1/x<=-(dps+10)*mp.log(10):
                # The same estimate gives |H^{-1}(w)-w|/|w|<3000S(radius).
                return pair
        compressed = _elementary_pole_radial(pair)
        raw = elementary_pole_raw_forward(compressed,infinity,dps=dps)
        return _elementary_pole_radial(raw,inverse=True)


def scalar_flat_step(u):
    """mpmath reference for the existing flat step, with exact endpoint pieces."""
    import mpmath as mp
    u = mp.mpf(u)
    if u<=0:
        return mp.mpf(0)
    if u>=1:
        return mp.mpf(1)
    a,b = mp.exp(-1/u),mp.exp(-1/(1-u))
    return a/(a+b)


def inverse_height_moment(a,q):
    """Elementary inverse-first moment f(a,q), |a|<=1/2, 0<=q<=1/8.

    f_a>=8/9. It equals a near zero and the prescribed exponential near each
    cell endpoint. This defines the inverse height, not a new named inverse.
    """
    import mpmath as mp
    a,q = mp.mpf(a),mp.mpf(q)
    if abs(a)>mp.mpf('.5') or not 0<=q<=mp.mpf('.125'):
        raise ValueError('inverse moment requires |a|<=1/2 and 0<=q<=1/8')
    aa = abs(a)
    if aa<=mp.mpf('.1'):
        return a
    weight = scalar_flat_step((aa-mp.mpf('.1'))/mp.mpf('.2'))
    E = mp.exp(2*(aa-mp.mpf('.5')))
    return mp.sign(a)*((1-weight)*aa+weight*(E-q/E))


def elementary_cusp_logradius(a,t,nu=1,dps=50):
    """Elementary inverse cusp height including its full collar.

    The output has rho(a+1)=rho(a)+P, with all jets matching at cell seams.
    No special function, implicit solve, series or quadrature is evaluated.
    """
    import mpmath as mp
    with mp.workdps(max(mp.mp.dps,dps)+20):
        a,t,nu = mp.mpf(a),mp.mpc(t),mp.mpf(nu)
        r = abs(t)
        if not 0<r<=mp.mpf('.5') or not 0<nu<=1:
            raise ValueError('elementary cusp inverse requires 0<|t|<=.5 and 0<nu<=1')
        P = mp.log(nu)+mp.log(r)-2*r*r/(1-r*r)
        q = mp.exp(2*P)
        n = mp.floor(a+mp.mpf('.5'))
        b = a-n
        f = inverse_height_moment(b,q)
        root = mp.sqrt(f*f+4*q)
        Y = 2*q/(root+f) if f>=0 else (root-f)/2
        rho0 = mp.log(Y)/2
        lam = scalar_flat_step((r-mp.mpf('.3'))/mp.mpf('.2'))
        return +(n*P+(1-lam)*rho0+lam*P*(b+mp.mpf('.5')))


def elementary_cusp_reference(t,a,nu=1):
    """Smooth core reference in the nonzero-coordinate gauge, including axes.

    |t|<=.3, |a|<.6, excluding the node. The neighbouring height cells give
    explicit open chart extensions across the central-circle representative
    seam; none of these formulas divides |t|^2 by a vanishing coordinate.
    """
    import mpmath as mp
    t,a,nu = mp.mpc(t),mp.mpf(a),mp.mpf(nu)
    if abs(t)>mp.mpf('.3') or abs(a)>=mp.mpf('.6'):
        raise ValueError('core reference requires |t|<=.3 and |a|<.6')
    e = mp.exp(-4*abs(t)**2/(1-abs(t)**2))
    q = e*nu*nu*abs(t)**2
    cell = int(mp.floor(a+mp.mpf('.5')))
    b = a-cell
    f = inverse_height_moment(b,q)
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


def scalar_lattice_coords(z,w,t,logt):
    """Fixed-size real Cramer's rule, matching tcore.lattice_coords."""
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


def elliptic_raw_base_reference(v,X,dps=50,elementary=False):
    """Closed reference over an UNMODIFIED algebraic-profile orbit map.

    Input is (v,X) in the S4 circle quotient, |v|+X^2=1, |X|<1.
    elementary=True uses programme 11's elementary inverse-first height;
    False uses the earlier elliptic height. Returns (kind,data) for the embedding.
    Reference phases are
    local gauges, not a global section. Pole smoothing is a separate block.
    """
    import mpmath as mp
    with mp.workdps(max(mp.mp.dps,dps)+30):
        vv = tuple(map(mp.mpf,v))
        xx = mp.mpf(X)
        if not abs(xx)<1:
            raise ValueError('the raw reference excludes the two fixed points')
        aa = xx*mp.sqrt(2-xx*xx)
        w = (vv[0]/(1-aa),-vv[2]/(1-aa),-vv[1]/(1-aa))
        y = (w[0],w[1]+1,w[2])
        s2 = sum(z*z for z in y)
        u1 = mp.mpc(1-s2,2*y[0])/(1+s2)
        u2 = mp.mpc(2*y[1],2*y[2])/(1+s2)
        A,B = (4*mp.pi)**2,mp.exp(4)
        # Always divide by the dominant S3 coordinate. The positive square-root
        # coefficient is smooth even when the numerator is exactly zero.
        inverse_ratio = abs(u1)>abs(u2)
        ratio = u2/u1 if inverse_ratio else u1/u2
        rr = abs(ratio)**2
        b = A*(1-rr)
        root = mp.sqrt(b*b+4*B*rr)
        small = ratio*mp.sqrt((2*B if inverse_ratio else 2)/(b+root))
        inf = inverse_ratio and abs(small)<mp.mpf('.5')
        t = None if inf else (1/small if inverse_ratio else small)
        if inf or abs(t)<mp.mpf('.5'):
            tc = small if inf else t
            phase_coordinate = -mp.conj(u1) if inf else u2
            # Two open argument charts; the half-centred chart avoids an
            # axis-only branch at the height +/-1/2 seam. On the selector
            # overlap they describe the same W point by the deck action.
            if phase_coordinate.real<-mp.mpf('.9')*abs(phase_coordinate):
                a = mp.mpf('.5')+mp.arg(-phase_coordinate)/(2*mp.pi)
            else:
                a = mp.arg(phase_coordinate)/(2*mp.pi)
            nu = mp.exp(-2) if inf else mp.mpf(1)
            if abs(tc)<=mp.mpf('.3') and elementary:
                xs,yy = elementary_cusp_reference(tc,a,nu)
                return ('ci' if inf else 'c0',(xs,yy,a))
            if tc==0:
                rad = elliptic_axis_radius(a)
                x,y = (mp.sqrt(rad)/nu,mp.mpf(0)) if a>0 else (mp.mpf(0),mp.sqrt(rad))
                return ('ci' if inf else 'c0',(mp.mpc(x),mp.mpc(y)))
            rho = (elementary_cusp_logradius if elementary else elliptic_cusp_logradius)(a,tc,nu,dps=dps)
            yy = mp.exp(rho)
            xs = tc/yy
            if abs(tc)<=mp.mpf('.3'):
                data = (xs,mp.mpc(yy),a) if elementary else (xs,mp.mpc(yy))
                return ('ci' if inf else 'c0',data)
            if inf:
                t = 1/tc
            lt = mp.log(t)
            if inf:
                zm,wm = -mp.j*rho,mp.mpc(0)
            else:
                ss = lt/(2*mp.j*mp.pi)
                zm,wm = -mp.j*rho+mp.pi*ss*(ss-1),2*t*ss
            lattice = scalar_lattice_coords(zm,wm,t,lt)
            data = (t,lattice[0],lattice[3],lt,a) if elementary else (t,lattice[0],lattice[3],lt)
            return ('mid',data)
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


def elliptic_meridian_projector(n,X,dps=50,elementary=False):
    """Closed projector on a raw-base meridian; no transport/integration.

    For elementary=True the fixed global bound is programme (11.13).
    The default elliptic alternative has no such global bound asserted here.
    """
    import mpmath as mp
    with mp.workdps(max(mp.mp.dps,dps)+20):
        h = mp.mpf(X)
        v = tuple((1-h*h)*mp.mpf(a) for a in n)
        kind,data = elliptic_raw_base_reference(v,h,dps=dps,elementary=elementary)
        e = elliptic_bundle_embedding(kind,data,dps=dps,elementary=elementary)
        return mp.matrix([[a*mp.conj(b) for b in e] for a in e])


def elliptic_bundle_embedding(kind,data,dps=50,with_norm=False,elementary=False):
    """Ten-component scalar embedding for the closed elliptic height.

    mid data=(t,x1,x4,logt); c0/ci data=(x,y), |xy|<=.3, |H|<=.6.
    elementary=True adds a known height to cusp/collar data. In middle data
    this height refers to the principal-log cusp representative; normalize
    logt before converting to that representative, preserving the Y point.
    Values are in S19, with the source circle acting by scalar multiplication.
    This has the same finite chart prescription as tlift.bundle_embedding_*.
    """
    import mpmath as mp
    with mp.workdps(max(mp.mp.dps,dps)+20):
        def hw(h):
            return 1-scalar_flat_step((abs(h)-mp.mpf('.25'))/mp.mpf('.15'))
        def frame(x,y,nu):
            r = nu*nu*abs(x)**2+abs(y)**2
            if not r>0:
                raise ValueError('the bundle embedding excludes the fixed points')
            return (nu*mp.conj(x)/mp.sqrt(r),y/mp.sqrt(r))
        def torus_pair(z,t,nu,known_height=None):
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
        vec = [mp.mpc(0) for _ in range(10)]
        if kind in ('c0','ci'):
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
        elif kind=='mid':
            t,one,four,lt = data[:4]
            known_height = mp.mpf(data[4]) if elementary and len(data)>4 else None
            t,lt,one,four = mp.mpc(t),mp.mpc(lt),mp.mpf(one),mp.mpf(four)
            principal = mp.log(t)
            branch = mp.nint((lt.imag-principal.imag)/(2*mp.pi))
            one,lt = one+branch*four,principal
            r = abs(t)
            b0 = 1-scalar_flat_step((r-mp.mpf('.3'))/mp.mpf('.2'))
            bi = 1-scalar_flat_step((1/r-mp.mpf('.3'))/mp.mpf('.2'))
            bm = 1-b0-bi
            z,w = 2*mp.pi*one-mp.j*lt*four,2*t*four
            def coeffs(rhs,e2,e3):
                det = e2.real*e3.imag-e3.real*e2.imag
                return ((rhs.real*e3.imag-e3.real*rhs.imag)/det,
                        (e2.real*rhs.imag-rhs.real*e2.imag)/det)
            if b0>0:
                s = lt/(2*mp.j*mp.pi)
                zc,wc = z-mp.pi*s*(s-1),w-2*t*s
                c2,c3 = coeffs(-wc,1-t*t,mp.j*(1+t*t))
                zz = mp.exp(mp.j*(zc-mp.j*c2*t-c3*t))
                vec[:4] = [b0*q for q in torus_pair(zz,t,mp.mpf(1),known_height)]
            if bi>0:
                s = 1/t
                zc,wc = z-mp.j*w/t,w/(t*t)
                c2,c3 = coeffs(-wc,s*s-1,mp.j*(1+s*s))
                zz = mp.exp(mp.j*(zc-mp.j*c2*s+c3*s))
                vec[4:8] = [bi*q for q in torus_pair(zz,s,mp.exp(-2),known_height)]
            if bm>0:
                for index,centre in enumerate((mp.mpf(0),mp.pi)):
                    angle = mp.arg(t if index==0 else -t)
                    weight = 1-scalar_flat_step((abs(angle)-mp.pi/2)/(mp.pi/4))
                    vec[8+index] = bm*weight*mp.exp(mp.j*(2*mp.pi*one-(angle+centre-lt.imag)*four))
        else:
            raise ValueError('kind must be mid, c0, or ci')
        norm2 = sum(abs(q)**2 for q in vec)
        result = tuple(+(q/mp.sqrt(norm2)) for q in vec)
        return (result,+norm2) if with_norm else result


def elementary_closed_parameters():
    """Exact integers/rationals in programme 11; no enormous powers evaluated.

    These counts specify a finite expression, not a practical numerical method.
    The executable evaluator has an explicit work guard for that reason.
    """
    from fractions import Fraction
    L = 10**114
    N = 4*L+1
    M = 32*(N+1)*L
    return dict(delta=Fraction(1,10**9),alpha=Fraction(1,10**5),
                epsilon=Fraction(1,10**8),projector_bound=L,
                projection_steps=N,phase_bound=M,phase_steps=8*M+1)


def _elementary_unit(v):
    import mpmath as mp
    length = mp.sqrt(sum(abs(a)**2 for a in v))
    if length==0:
        raise ValueError('zero denominator in a finite projection')
    return tuple(a/length for a in v)


def _elementary_inner(u,v):
    import mpmath as mp
    return sum(mp.conj(a)*b for a,b in zip(u,v))


def _elementary_pole_inclusion(xi,infinity=False):
    """I_j xi: U_j^* xi in the indicated C2 block of C10."""
    import mpmath as mp
    a,b = xi
    pair = ((a+b)/mp.sqrt(2),(a-b)/mp.sqrt(2)) if infinity else (
        (a-mp.j*b)/mp.sqrt(2),(a+mp.j*b)/mp.sqrt(2))
    result = [mp.mpc(0)]*10
    start = 4 if infinity else 0
    result[start:start+2] = pair
    return tuple(result)


def elementary_projection_section(xi,X,steps,dps=50,diagnostics=None):
    """The exact finite product section for a specified positive integer count.

    At the certified count this is the programme's global slab lift. Smaller
    counts are independently defined numerical instances, not certified ones.
    diagnostics records the smallest evaluated projection norm (not a proof).
    """
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


def elementary_phase_family(g,n,X,anchor,steps,chart=None,diagnostics=None):
    """Finite real logarithm of a specified S2 x interval -> S1 family.

    The constant Arg is evaluated only at (north,anchor). All variable increments
    must lie in the right half-plane. This is programme 10.6 with an extra finite
    path in X; no unwrapping algorithm or quadrature is used.
    """
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


def elementary_global_slice(Z,X,*,projection_steps=None,phase_steps=None,dps=50,
                            max_work=10**6,diagnostics=None,force_slab=False):
    """Finite elementary G_W of programme 11, as a W-chart representative.

    Defaults are the proved, extremely large fixed counts. max_work prevents an
    accidental astronomical evaluation; None removes that execution guard.
    Supplying smaller counts evaluates a different finite map and provides NO
    global guarantee. Such counts are useful for independent local checks.
    force_slab is only for checking the overlap with the pole formula.
    """
    import mpmath as mp
    with mp.workdps(max(mp.mp.dps,dps)+20):
        Z = tuple(map(mp.mpc,Z))
        X = mp.mpf(X)
        alpha = mp.mpf('1e-5')
        distance = 1-abs(X)
        if distance<alpha/2 and not force_slab:
            inf = X>0
            scale = mp.pi if inf else 2*mp.pi
            den = mp.sqrt(scale*(1+abs(X)*mp.sqrt(2-X*X)))
            pair = _elementary_pole_inclusion(Z,inf)[4:6] if inf else _elementary_pole_inclusion(Z)[:2]
            zz = elementary_pole_patch(tuple(a/den for a in pair),inf,dps=dps)
            nu = mp.exp(-2) if inf else mp.mpf(1)
            # A third entry is unnecessary: the torus restoration uses H_f.
            return ('ci' if inf else 'c0',(mp.conj(zz[0])/nu,zz[1]))
        if not -1+mp.mpf('1e-9')<=X<=1-mp.mpf('1e-9'):
            raise ValueError('the finite slab formula excludes the inner pole caps')
        params = elementary_closed_parameters()
        N = params['projection_steps'] if projection_steps is None else projection_steps
        K = params['phase_steps'] if phase_steps is None else phase_steps
        cap_band = distance<alpha
        cost = (N+1)*(12*(K+1) if cap_band else 1)
        if max_work is not None and cost>max_work:
            raise ValueError('certified finite expression exceeds max_work; inspect elementary_closed_parameters()')
        xi = tuple(a/mp.sqrt(1-X*X) for a in Z)
        pc = 2*xi[0]*mp.conj(xi[1])
        n = (pc.real,pc.imag,abs(xi[0])**2-abs(xi[1])**2)
        result = elementary_projection_section(xi,X,N,dps,diagnostics)
        if cap_band:
            inf = X>0
            anchor = 1-alpha if inf else -1+alpha
            def mismatch(nn,h):
                # Smooth reference on either hemisphere; mismatch is phase invariant.
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
        kind,data = elliptic_raw_base_reference(tuple((1-X*X)*a for a in n),X,dps=dps,elementary=True)
        ref = elliptic_bundle_embedding(kind,data,dps=dps,elementary=True)
        phase = _elementary_inner(ref,result)
        if kind=='mid':
            t,one,four,lt = data[:4]
            return kind,(t,one+mp.arg(phase)/(2*mp.pi),four,lt)
        x,y = data[:2]
        return kind,(mp.conj(phase)*x,phase*y)


def scalar_finite_height(x,y,nu=1):
    """Four-seed H_f for the bounded representatives used by the closed G.

    The more general torch implementation is separate. Here A,Y<=4 suffices
    for the exact q^2 tail; an explicit assertion guards this formula's domain.
    """
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


def elementary_restore_torus(kind,data,phi):
    """Closed tau^{-1}, using finite H_f, in a holomorphic Y chart.

    Output (kind,(z1,z2,z3)); the real Jacobian gives the coframe in (8.14).
    The inverse-first base height is not needed here.
    """
    import mpmath as mp
    phi = tuple(map(mp.mpf,phi))
    def coeffs(rhs,e2,e3):
        det = e2.real*e3.imag-e3.real*e2.imag
        return ((rhs.real*e3.imag-e3.real*rhs.imag)/det,
                (e2.real*rhs.imag-rhs.real*e2.imag)/det)
    def cusp_data(t,inf):
        return ((t*t-1,mp.j*(1+t*t),-2*t,t,mp.j*t) if inf else
                (1-t*t,mp.j*(1+t*t),2*t,t,-mp.j*t))
    if kind in ('c0','ci'):
        x,y = data[:2]
        inf = kind=='ci'
        t = x*y
        e2,e3,e4,a2,a3 = cusp_data(t,inf)
        c = coeffs(e4,e2,e3)
        H = scalar_finite_height(x,y,mp.exp(-2) if inf else mp.mpf(1))
        C = tuple(v+H*w for v,w in zip(phi,c))
        A = C[0]*a2+C[1]*a3
        return kind,(x*mp.exp(-A),y*mp.exp(A),C[0]*e2+C[1]*e3)
    if kind!='mid':
        raise ValueError('kind must be mid, c0 or ci')
    t,one,four,lt = data[:4]
    r = abs(t)
    z,w = 2*mp.pi*one-mp.j*lt*four,2*t*four
    C = list(phi)
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
    return 'mid',(z-mp.j*t*C[0]-t*C[1],w+(1-t*t)*C[0]+mp.j*(1+t*t)*C[1],t)


def elementary_global_coordinates(Z,X,phi,**kwargs):
    """Finite elementary G, returning a holomorphic Y-chart representative."""
    import mpmath as mp
    with mp.workdps(max(mp.mp.dps,kwargs.get('dps',50))+20):
        kind,data = elementary_global_slice(Z,X,**kwargs)
        return elementary_restore_torus(kind,data,phi)


def elementary_tensor_numeric(u,missing=4,sign=1,**kwargs):
    """Independent high-precision local check of (11.19), using mp.diff.

    u consists of four retained R5 coordinates and two torus angles. This is a
    numerical derivative evaluator, not the definition of the symbolic tensor.
    The latter uses the finite chain rules in the programme. As for G, the
    certified default counts are protected by max_work; reduced counts test
    different finite maps on the regions where their denominators stay nonzero.
    """
    import mpmath as mp
    if missing not in range(5) or sign not in (-1,1):
        raise ValueError('invalid round-sphere chart')
    with mp.workdps(max(mp.mp.dps,kwargs.get('dps',50))+20):
        u = mp.matrix(u)
        def coordinates(v):
            q = list(v[:4])
            q.insert(missing,sign*mp.sqrt(1-sum(a*a for a in q)))
            kind,h = elementary_global_coordinates((mp.mpc(q[0],q[1]),mp.mpc(q[2],q[3])),q[4],v[4:6],**kwargs)
            return kind,mp.matrix([b for a in h for b in (a.real,a.imag)])
        kind,_ = coordinates(u)
        matrix = mp.matrix(6)
        for j in range(6):
            direction = mp.eye(6)[:,j]
            def evaluate(t):
                current,result = coordinates(u+t*direction)
                if current!=kind:
                    raise ValueError('choose an interior source chart for the numerical derivative')
                return result
            matrix[:,j] = mp.diff(evaluate,0)
        standard = mp.matrix(6)
        for j in range(3):
            standard[2*j,2*j+1] = -1
            standard[2*j+1,2*j] = 1
        return matrix**-1*standard*matrix,matrix
