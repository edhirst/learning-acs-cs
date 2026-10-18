"""Explicit direct complex structure on S4 x (R/Z)^2.

The definition is programme section 12, with exact formula branches.
No historical test script is executed at import. Numerical evaluation
approximates elementary real functions; see proof_contracts() for the
written global premises. The public matrix uses the global frame described
in complex_structure(), not a coordinate frame.
"""
from __future__ import annotations
from contextlib import nullcontext
from dataclasses import dataclass
import mpmath as mp
import numpy as np

PHASE_STEPS = 65536

def real_array(a):
    return np.vectorize(mp.re, otypes=[object])(a)

def imag_array(a):
    return np.vectorize(mp.im, otypes=[object])(a)

class Jet:
    """Six real input directions; complex value, gradient and Hessian."""

    def __init__(self, value, gradient=None, hessian=None):
        self.value = mp.mpc(value)
        self.gradient = np.zeros(6, object) if gradient is None else np.asarray(gradient, object)
        self.hessian = np.zeros((6, 6), object) if hessian is None else np.asarray(hessian, object)

    @staticmethod
    def coerce(x):
        return x if isinstance(x, Jet) else Jet(x)

    def unary(self, value, first, second):
        return Jet(value, first * self.gradient, first * self.hessian + second * np.outer(self.gradient, self.gradient))

    def __add__(self, other):
        if isinstance(other, Matrix):
            return NotImplemented
        other = self.coerce(other)
        return Jet(self.value + other.value, self.gradient + other.gradient, self.hessian + other.hessian)
    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.value, -self.gradient, -self.hessian)

    def __pos__(self):
        return self

    def __sub__(self, other):
        return self + -self.coerce(other)

    def __rsub__(self, other):
        return self.coerce(other) + -self

    def __mul__(self, other):
        if isinstance(other, Matrix):
            return NotImplemented
        other = self.coerce(other)
        return Jet(self.value * other.value, self.gradient * other.value + self.value * other.gradient, self.hessian * other.value + self.value * other.hessian + np.outer(self.gradient, other.gradient) + np.outer(other.gradient, self.gradient))
    __rmul__ = __mul__

    def __pow__(self, power):
        if isinstance(power, Jet):
            return arithmetic.exp(power * arithmetic.log(self))
        if power == 0:
            return Jet(1)
        if power == 1:
            return self
        if power == 2:
            return self * self
        if self.value == 0:
            if power > 2:
                return Jet(0)
            if power < 0:
                raise ZeroDivisionError('negative power of a zero jet')
            return Jet(0, np.full(6, np.nan), np.full((6, 6), np.nan))
        return self.unary(self.value ** power, power * self.value ** (power - 1), power * (power - 1) * self.value ** (power - 2))

    def __truediv__(self, other):
        return self * self.coerce(other) ** (-1)

    def __rtruediv__(self, other):
        return self.coerce(other) * self ** (-1)

    @property
    def real(self):
        return Jet(self.value.real, real_array(self.gradient), real_array(self.hessian))

    @property
    def imag(self):
        return Jet(self.value.imag, imag_array(self.gradient), imag_array(self.hessian))

    def conjugate(self):
        return Jet(self.value.conjugate(), self.gradient.conj(), self.hessian.conj())

    def __abs__(self):
        return norm_squared(self) ** 0.5

    def __complex__(self):
        return complex(self.value)

    def __float__(self):
        if self.value.imag != 0:
            raise TypeError('nonreal jet comparison/conversion')
        return float(self.value.real)

    def __int__(self):
        return int(float(self))

    def __bool__(self):
        return self.value != 0

    def __lt__(self, other):
        return self.value.real < self.coerce(other).value.real

    def __le__(self, other):
        return self.value.real <= self.coerce(other).value.real

    def __gt__(self, other):
        return self.value.real > self.coerce(other).value.real

    def __ge__(self, other):
        return self.value.real >= self.coerce(other).value.real

    def __eq__(self, other):
        try:
            return self.value == self.coerce(other).value
        except (TypeError, ValueError):
            return False

    def __ne__(self, other):
        return not self == other

def norm_squared(x):
    return x.real * x.real + x.imag * x.imag

class Matrix:

    def __init__(self, *args):
        if len(args) > 1 or isinstance(args[0], int):
            rows = args[0]
            cols = args[1] if len(args) > 1 else rows
            self.data = np.zeros((rows, cols), object)
        else:
            self.data = np.array(args[0], dtype=object)
            if self.data.ndim == 1:
                self.data = self.data[:, None]

    def __iter__(self):
        return iter(self.data.flat)

    def __getitem__(self, index):
        if not isinstance(index, tuple):
            return self.data.flat[index]
        result = self.data[index]
        if not isinstance(result, np.ndarray):
            return result
        if result.ndim == 1:
            result = result[None, :] if isinstance(index[0], int) else result[:, None]
        return Matrix(result)

    def __setitem__(self, index, value):
        self.data[index] = value.data if isinstance(value, Matrix) else value

    @property
    def H(self):
        return Matrix([[arithmetic.conj(x) for x in row] for row in self.data.T])

    @property
    def T(self):
        return Matrix(self.data.T)

    def __mul__(self, other):
        if isinstance(other, Matrix):
            return Matrix(self.data @ other.data)
        return Matrix([[x * other for x in row] for row in self.data])

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        return Matrix([[x / other for x in row] for row in self.data])

    def __add__(self, other):
        return Matrix(self.data + other.data)

    def __sub__(self, other):
        return Matrix(self.data - other.data)

class Arithmetic:
    pi = mp.pi
    e = mp.e
    j = mp.j
    inf = mp.inf
    mp = mp.mp

    @staticmethod
    def mpf(x):
        if isinstance(x, Jet):
            if x.value.imag != 0:
                raise TypeError('mpf received a nonreal jet')
            return x.real
        return mp.mpf(x)

    @staticmethod
    def mpc(x=0, y=None):
        if y is not None:
            return x + 1j * y
        return x if isinstance(x, Jet) else mp.mpc(x)

    @staticmethod
    def workdps(dps):
        return nullcontext()

    @staticmethod
    def conj(x):
        return x.conjugate() if hasattr(x, 'conjugate') else x

    @staticmethod
    def re(x):
        return x.real

    @staticmethod
    def im(x):
        return x.imag

    @staticmethod
    def floor(x):
        return mp.floor(x.value.real if isinstance(x, Jet) else x)

    @staticmethod
    def nint(x):
        return mp.nint(x.value.real if isinstance(x, Jet) else x)

    @staticmethod
    def sign(x):
        value = x.value.real if isinstance(x, Jet) else x
        if value == 0 and isinstance(x, Jet):
            return Jet(0, np.full(6, np.nan), np.full((6, 6), np.nan))
        return mp.sign(value)

    @staticmethod
    def sqrt(x):
        return x ** 0.5 if isinstance(x, Jet) else mp.sqrt(x)

    @staticmethod
    def exp(x):
        if not isinstance(x, Jet):
            return mp.exp(x)
        value = mp.exp(x.value)
        return x.unary(value, value, value)

    @staticmethod
    def log(x):
        if not isinstance(x, Jet):
            return mp.log(x)
        return x.unary(mp.log(x.value), 1 / x.value, -1 / x.value ** 2)

    @staticmethod
    def log1p(x):
        if not isinstance(x, Jet):
            return mp.log1p(x)
        return x.unary(mp.log1p(x.value), 1 / (1 + x.value), -1 / (1 + x.value) ** 2)

    @staticmethod
    def expm1(x):
        if not isinstance(x, Jet):
            return mp.expm1(x)
        value = mp.exp(x.value)
        return x.unary(mp.expm1(x.value), value, value)

    @staticmethod
    def sin(x):
        if not isinstance(x, Jet):
            return mp.sin(x)
        return x.unary(mp.sin(x.value), mp.cos(x.value), -mp.sin(x.value))

    @staticmethod
    def cos(x):
        if not isinstance(x, Jet):
            return mp.cos(x)
        return x.unary(mp.cos(x.value), -mp.sin(x.value), -mp.cos(x.value))

    @staticmethod
    def sinh(x):
        if not isinstance(x, Jet):
            return mp.sinh(x)
        return x.unary(mp.sinh(x.value), mp.cosh(x.value), mp.sinh(x.value))

    @staticmethod
    def cosh(x):
        if not isinstance(x, Jet):
            return mp.cosh(x)
        return x.unary(mp.cosh(x.value), mp.sinh(x.value), mp.cosh(x.value))

    @staticmethod
    def asinh(x):
        if not isinstance(x, Jet):
            return mp.asinh(x)
        return x.unary(mp.asinh(x.value), 1 / mp.sqrt(1 + x.value * x.value), -x.value / (1 + x.value * x.value) ** 1.5)

    @staticmethod
    def arg(x):
        if not isinstance(x, Jet):
            return mp.arg(x)
        return Jet(mp.arg(x.value), imag_array(x.gradient / x.value), imag_array(x.hessian / x.value - np.outer(x.gradient, x.gradient) / x.value ** 2))
    matrix = Matrix

    @staticmethod
    def zeros(rows, cols=None):
        return Matrix(rows, rows if cols is None else cols)

    @staticmethod
    def diag(entries):
        return Matrix(np.diag(np.array(entries, dtype=object)))

    @staticmethod
    def norm(entries):
        return arithmetic.sqrt(sum((norm_squared(x) for x in entries)))

arithmetic = Arithmetic()

def boundary_transition(theta, phi, infinity=False):
    """Vectorized two-row formula for delta at the EXACT cusp radius 1/2.

    Float64 evaluation is numerical evidence. The finite elementary formula
    and all normalization domains are certified in programme 12.7--12.8.
    Inactive local rows are not divided by their possibly zero norms.
    """
    (th, ph) = np.broadcast_arrays(np.asarray(theta, dtype=float), np.asarray(phi, dtype=float))

    def step(x):
        xx = np.where((x > 0) & (x < 1), x, 0.5)
        aa = np.exp(-1 / xx)
        bb = np.exp(-1 / (1 - xx))
        return np.where(x <= 0, 0.0, np.where(x >= 1, 1.0, aa / (aa + bb)))

    def unit(row, active):
        norm = np.sqrt(np.sum(np.abs(row) ** 2, axis=-1))
        if np.any(active & (norm == 0)):
            raise ValueError('an active boundary row vanished')
        return row / np.where(active, norm, 1.0)[..., None]
    r = 2.0 if infinity else 0.5
    ratio2 = r * r * ((4 * np.pi) ** 2 + np.exp(4) * r * r) / (1 + (4 * np.pi) ** 2 * r * r)
    (a, b) = (np.sqrt(ratio2 / (1 + ratio2)), 1 / np.sqrt(1 + ratio2))
    u1 = a * np.exp(1j * (th + ph))
    u2 = b * np.exp(1j * ph)
    p = np.stack((u1.imag, -u2.imag, 1 + u1.real - u2.real), axis=-1)
    n = p / np.linalg.norm(p, axis=-1)[..., None]
    Pi = np.empty(th.shape + (2, 2), complex)
    Pi[..., 0, 0] = (1 + n[..., 2]) / 2
    Pi[..., 1, 1] = (1 - n[..., 2]) / 2
    Pi[..., 0, 1] = (n[..., 0] + 1j * n[..., 1]) / 2
    Pi[..., 1, 0] = Pi[..., 0, 1].conj()

    def project(v):
        return np.einsum('...i,...ij->...j', v.conj(), Pi)
    H = np.angle(-u1.conj() if infinity else u2) / (2 * np.pi)
    tc = -th if infinity else th
    k = np.log(2) + 2 / 3 + (2 if infinity else 0)
    den = np.sqrt(2 * np.cosh(2 * k * H))
    q0 = np.exp(k * H - 1j * tc) / den
    q1 = np.exp(-k * H) / den
    Uphase = 1 if infinity else 1j
    F = np.stack(((q0 + q1) / np.sqrt(2), Uphase * (q0 - q1) / np.sqrt(2)), axis=-1)
    ch = step((np.abs(H) - 0.3) / 0.05)
    c0 = unit(project(F), ch < 1)
    H1 = np.where(H >= 0, H, H + 1)
    tau = np.pi / 2 * step((H1 - 0.35) / 0.3)
    v = np.stack(((np.cos(tau) + np.sin(tau)) / np.sqrt(2), Uphase * (np.cos(tau) - np.sin(tau)) / np.sqrt(2)), axis=-1)
    c1 = np.exp(1j * np.where(H >= 0, tc, 0))[..., None] * unit(project(v), ch > 0)
    cr = unit((1 - ch)[..., None] * c0 + ch[..., None] * c1, np.ones(th.shape, bool))
    th1 = np.where(th >= 0, th, th + 2 * np.pi)
    cs = 1 / (3 * np.log(2) + 2)
    C = th1 / (2 * np.pi) + (0 if infinity else cs)
    w = step((th1 - 5 * np.pi / 8) / (3 * np.pi / 4))
    four = -H - 0.5 if infinity else H + th / (2 * np.pi) + cs
    z = np.exp(1j * ph)
    vm = np.stack((1j * (a + 1j * b) * (z - 1) / 2, (z + 1) / 2), axis=-1)
    cm = step((np.abs(th) - 3 * np.pi / 8) / (np.pi / 4))
    m0 = unit(Pi[..., 0, :], cm < 1)
    phase = (b - 1j * a) * np.exp(1j * (w * (np.pi - 2 * np.pi * C) + (th1 - th) * four))
    m1 = phase[..., None] * unit(project(vm), cm > 0)
    mr = unit((1 - cm)[..., None] * m0 + cm[..., None] * m1, np.ones(th.shape, bool))
    ell = -np.log(2)
    P = ell - 2 / 3
    xi1 = -(th * th + ell * ell) / (8 * np.pi ** 2) - th / (4 * np.pi) - ell / (10 * np.pi ** 2) + th * ell / (4 * np.pi * P)
    one = th * (H + 0.5) / (2 * np.pi) if infinity else xi1 - th * (H + 0.5) / (2 * np.pi)
    return np.exp(2j * np.pi * one) * np.sum(mr.conj() * cr, axis=-1)

def phase_value(theta, phi, steps=65536):
    """Certified-count finite phase sum on [-pi,pi]^2, in float64.

    Any other count denotes a different, uncertified numerical experiment.
    The proof uses the exact arithmetic expression, not float64 residuals.
    No count-dependent source projection chain occurs here.
    """
    if not isinstance(steps, int) or steps <= 0:
        raise ValueError('steps must be a positive integer')
    if abs(theta) > np.pi or abs(phi) > np.pi:
        raise ValueError('use principal input angles in [-pi,pi]')
    s = np.arange(steps + 1, dtype=float) / steps
    g = boundary_transition(s * theta, s * phi, True) / boundary_transition(s * theta, s * phi)
    adjacent = g[1:] / g[:-1]
    if np.any(adjacent.real <= 0):
        raise ValueError('phase increment left its certified right half-plane')
    h = np.log(2)
    alpha = h * (1 / 5 - h / 4) / np.pi
    phase = -alpha + np.sum(np.arctan(adjacent.imag / adjacent.real), dtype=np.longdouble)
    return (float(phase), float(adjacent.real.min()))

class Formulas:
    def __init__(self, ops=arithmetic):
        self.ops = ops

    def pole_raw_inverse(self, z, infinity=False, dps=50):
        """Closed local inverse lift before smoothing at the pole (programme 11.3).

        z is the two-complex-coordinate linear pole model, not a sphere point.
        Return the slice pair (sqrt(a0)*conj(x),y). Defined on 0<|z|^2<=1e-4;
        no implicit root or iteration is used. Smoothness at z=0 is NOT asserted.
        """
        mp = self.ops
        with mp.workdps(max(mp.mp.dps, dps) + 20):
            (z0, z1) = map(mp.mpc, z)
            radius = norm_squared(z0) + norm_squared(z1)
            if not 0 < radius <= mp.mpf('1e-4'):
                raise ValueError('raw pole inverse requires 0<|z|^2<=1e-4')
            pc = 2 * z0 * mp.conj(z1)
            p = (pc.real, pc.imag, norm_squared(z0) - norm_squared(z1))
            scale = mp.pi if infinity else 2 * mp.pi
            v = tuple((scale * w for w in ((p[2], -p[1], p[0]) if infinity else (-p[1], -p[2], p[0]))))
            (aa, bb, cc) = (v[0], -v[2], -v[1])
            if infinity:
                s2 = aa * aa + bb * bb + cc * cc
                den = 1 + 2 * bb + 2 * s2
                u1 = mp.mpc(-1 - 2 * bb, 2 * aa) / den
                u2 = mp.mpc(2 * (bb + s2), 2 * cc) / den
                ratio = u2 / u1
                height = -mp.arg(-u1) / (2 * mp.pi)
            else:
                s2 = aa * aa + (bb + 1) ** 2 + cc * cc
                u1 = mp.mpc(-(aa * aa + bb * bb + cc * cc + 2 * bb), 2 * aa) / (1 + s2)
                u2 = mp.mpc(2 * (bb + 1), 2 * cc) / (1 + s2)
                ratio = u1 / u2
                height = mp.arg(u2) / (2 * mp.pi)
            (A, B) = ((4 * mp.pi) ** 2, mp.exp(4))
            uu = norm_squared(ratio)
            den = A * (1 - uu) + mp.sqrt(A * A * (1 - uu) ** 2 + 4 * B * uu)
            tc = ratio * mp.sqrt((2 * B if infinity else 2) / den)
            nu = mp.exp(-2) if infinity else mp.mpf(1)
            if abs(height) >= mp.mpf('.1'):
                raise ValueError('elementary pole inverse left its certified height region')
            e = mp.exp(-4 * norm_squared(tc) / (1 - norm_squared(tc)))
            q = e * nu * nu * norm_squared(tc)
            R = mp.sqrt(height * height + 4 * q)
            AA = 2 * q / (R - height) if height < 0 else (R + height) / 2
            yy = 2 * q / (R + height) if height >= 0 else (R - height) / 2
            xx = AA / e
            (moment, source_radius) = (xx - yy, xx + yy)
            source_pc = 2 * nu * mp.conj(tc)
            w0 = ((source_radius + moment) * z0 + source_pc * z1) / (2 * source_radius)
            w1 = (mp.conj(source_pc) * z0 + (source_radius - moment) * z1) / (2 * source_radius)
            length = mp.sqrt(norm_squared(w0) + norm_squared(w1))
            if length == 0:
                raise ValueError('local Hopf projections became orthogonal')
            factor = mp.sqrt(source_radius) / length
            return (+(factor * w0), +(factor * w1))

    def _elementary_pole_radial(self, pair, inverse=False):
        """The radial conjugation D (or D inverse) in programme (11.12g--h)."""
        mp = self.ops
        z = tuple(map(mp.mpc, pair))
        radius = sum((norm_squared(w) for w in z))
        (a, c) = (mp.mpf('1e-8'), mp.mpf('1.5e-8'))
        if radius >= c:
            return z
        if inverse:
            if radius <= 0:
                raise ValueError('radial inverse requires positive radius')
            y = mp.e * mp.expm1(radius / (c - radius))
            L = y + mp.log1p(-mp.exp(-y)) if y > 1 else mp.log(mp.expm1(y))
            x = 2 / (mp.sqrt(L * L + 4) - L) if L < 0 else (L + mp.sqrt(L * L + 4)) / 2
            ell = mp.log1p(x / mp.exp(3))
            new_radius = a + c * ell / (3 * (3 + ell))
        else:
            if radius <= a:
                raise ValueError('radial compression requires radius greater than 1e-8')
            x = mp.exp(3) * mp.expm1(3 * (radius - a) / (c - radius))
            f = x - 1 / x
            y = mp.log1p(mp.exp(f)) if f < 0 else f + mp.log1p(mp.exp(-f))
            ell = mp.log1p(y / mp.e)
            new_radius = c * ell / (1 + ell)
        factor = mp.sqrt(new_radius / radius)
        return tuple((+(factor * w) for w in z))

    def elementary_pole_raw_forward(self, w, infinity=False, dps=50):
        """Finite elementary inverse of the raw local lift (11.10), away from zero."""
        mp = self.ops
        with mp.workdps(max(mp.mp.dps, dps) + 20):
            pair = tuple(map(mp.mpc, w))
            radius = sum((norm_squared(a) for a in pair))
            if radius <= 0:
                raise ValueError('raw pole forward map requires nonzero spinor')
            nu = mp.exp(-2) if infinity else mp.mpf(1)
            tc = mp.conj(pair[0]) * pair[1] / nu
            (u, v) = (norm_squared(pair[0]), norm_squared(pair[1]))
            e = mp.exp(-4 * norm_squared(tc) / (1 - norm_squared(tc)))
            height = e * u - v
            (A, B) = ((4 * mp.pi) ** 2, mp.exp(4))
            if infinity:
                ratio = tc * mp.sqrt((norm_squared(tc) + A) / (A * norm_squared(tc) + B))
                u1 = mp.exp(2 * mp.pi * mp.j * (-height + mp.mpf('.5'))) / mp.sqrt(1 + norm_squared(ratio))
                u2 = ratio * u1
            else:
                ratio = tc * mp.sqrt((A + B * norm_squared(tc)) / (1 + A * norm_squared(tc)))
                u2 = mp.exp(2 * mp.pi * mp.j * height) / mp.sqrt(1 + norm_squared(ratio))
                u1 = ratio * u2
            if infinity:
                den = 2 * (1 - u2.real)
                stereo = (u1.imag / den, (u2.real - 1 - u1.real) / den, u2.imag / den)
            else:
                den = 1 + u1.real
                stereo = (u1.imag / den, u2.real / den - 1, u2.imag / den)
            (aa, bb, cc) = (stereo[0], -stereo[2], -stereo[1])
            scale = mp.pi if infinity else 2 * mp.pi
            p = tuple((a / scale for a in (aa, bb, cc)))
            pL = (p[2], -p[1], p[0]) if infinity else (p[2], -p[0], -p[1])
            rL = mp.sqrt(sum((a * a for a in pL)))
            n = tuple((a / rL for a in pL))
            projected = ((1 + n[2]) * pair[0] + (n[0] + mp.j * n[1]) * pair[1], (n[0] - mp.j * n[1]) * pair[0] + (1 - n[2]) * pair[1])
            length = mp.sqrt(sum((norm_squared(a) for a in projected)))
            if length == 0:
                raise ValueError('raw pole forward projector is orthogonal')
            factor = mp.sqrt(rL) / length
            return tuple((+(factor * a) for a in projected))

    def elementary_pole_patch(self, z, infinity=False, dps=50):
        """Two-sided elementary pole patch H_j of programme (11.12h)."""
        mp = self.ops
        with mp.workdps(max(mp.mp.dps, dps) + 50):
            pair = tuple(map(mp.mpc, z))
            inner = mp.mpf('1e-8')
            outer = mp.mpf('1.5e-8')
            radius = sum((norm_squared(a) for a in pair))
            if radius <= inner:
                return pair
            if radius < (inner + outer) / 2:
                x = mp.exp(3) * mp.expm1(3 * (radius - inner) / (outer - radius))
            compressed = self._elementary_pole_radial(pair)
            raw = self.pole_raw_inverse(compressed, infinity, dps=dps)
            return self._elementary_pole_radial(raw, inverse=True)

    def elementary_pole_forward(self, w, infinity=False, dps=50):
        """Finite elementary inverse H_j^{-1} on the certified cap image."""
        mp = self.ops
        with mp.workdps(max(mp.mp.dps, dps) + 50):
            pair = tuple(map(mp.mpc, w))
            inner = mp.mpf('1e-8')
            outer = mp.mpf('1.5e-8')
            radius = sum((norm_squared(a) for a in pair))
            if radius <= inner:
                return pair
            if radius < (inner + outer) / 2:
                x = mp.exp(3) * mp.expm1(3 * (radius - inner) / (outer - radius))
            compressed = self._elementary_pole_radial(pair)
            raw = self.elementary_pole_raw_forward(compressed, infinity, dps=dps)
            return self._elementary_pole_radial(raw, inverse=True)

    def scalar_flat_step(self, u):
        """mpmath reference for the existing flat step, with exact endpoint pieces."""
        mp = self.ops
        u = mp.mpf(u)
        if u <= 0:
            return mp.mpf(0)
        if u >= 1:
            return mp.mpf(1)
        (a, b) = (mp.exp(-1 / u), mp.exp(-1 / (1 - u)))
        return a / (a + b)

    def _terminal_height_data(self, q, lam):
        """Parameters of the terminal-cut height conjugation."""
        mp = self.ops
        if lam == 0:
            (kappa, amp, d) = (mp.mpf(2), mp.mpf(0), mp.mpf(1))
        else:
            if q <= 0:
                raise ValueError('positive collar weight requires q>0')
            k = -mp.log(q)
            kappa = 2 + lam * (k - 2)
            amp = lam * k
            d = 1 - lam + 2 * lam * mp.sqrt(q) * k
        (B, C) = (mp.exp(-kappa / 2), q * mp.exp(kappa / 2))
        z0 = B - C
        zcut = B * mp.exp(mp.mpf('.4') * kappa) - C * mp.exp(-mp.mpf('.4') * kappa)
        eta = mp.mpf('.01')
        E0 = mp.exp(eta / zcut)
        Acut = mp.exp(eta / (zcut - z0))
        Lambda = mp.exp(Acut - E0)
        Eprime0 = E0 * eta / zcut ** 2
        Eprime_z0 = Acut * eta / (zcut - z0) ** 2
        mu = d * Eprime0 / (kappa * (B + C) * Eprime_z0)
        return (kappa, amp, d, B, C, z0, zcut, eta, E0, Acut, Lambda, mu)

    def _height_G(self, a, amp, d, inverse=False):
        mp = self.ops
        if amp == 0:
            return a / d if inverse else d * a
        return mp.asinh(amp * a / d) / amp if inverse else d * mp.sinh(amp * a) / amp

    def _height_R(self, value, inverse=False):
        """Shifted flat map R:(.1,infinity)->(0,infinity), or its inverse."""
        mp = self.ops
        (delta, c) = (mp.mpf('.1'), mp.mpf('.2'))
        if inverse:
            if value <= 0:
                return delta
        elif value <= delta:
            return mp.mpf(0)
        if value >= c:
            return value
        if inverse:
            y = mp.e * mp.expm1(value / (c - value))
            L = y + mp.log1p(-mp.exp(-y)) if y > 1 else mp.log(mp.expm1(y))
            x = 2 / (mp.sqrt(L * L + 4) - L) if L < 0 else (L + mp.sqrt(L * L + 4)) / 2
            ell = mp.log1p(x / mp.exp(2))
            return delta + c * ell / (2 * (2 + ell))
        x = mp.exp(2) * mp.expm1(2 * (value - delta) / (c - value))
        F = x - 1 / x
        soft = mp.log1p(mp.exp(F)) if F < 0 else F + mp.log1p(mp.exp(-F))
        ell = mp.log1p(soft / mp.e)
        return c * ell / (1 + ell)

    def _height_terminal_forward(self, delta_z, data):
        """T(z0+delta_z), using shifted differences at z0."""
        mp = self.ops
        (_, _, _, _, _, z0, zcut, eta, E0, Acut, Lambda, mu) = data
        gap0 = zcut - z0
        if delta_z <= 0:
            return mp.mpf(0)
        if delta_z >= gap0:
            return z0 + delta_z
        gap = gap0 - delta_z
        L = eta / gap
        Ez = mp.exp(L)
        (C0, C1) = (1 + mu - 2 * Lambda, Lambda - mu)
        K = 1 + abs(C0) * mp.exp(E0) + abs(C1) * mp.exp(E0 + Acut)
        if Ez > 50:
            D = C0 * mp.exp(E0 - Ez) + C1 * mp.exp(E0 + Acut - 2 * Ez)
            correction = mp.log1p(D)
            ell = mp.log1p(correction / Ez)
            z = z0 + delta_z
            return z + eta * ell / (L * (L + ell))
        x = Acut * mp.expm1(eta * delta_z / (gap * gap0))
        u = mp.expm1(x)
        F = (Lambda * u * u + mu * u) / (1 + u)
        M = mp.log1p(F)
        ell = mp.log1p(M / E0)
        return zcut * ell / (eta / zcut + ell)

    def _height_terminal_inverse_delta(self, y, data):
        """T^{-1}(y)-z0, including the flat terminal-cut seam."""
        mp = self.ops
        (_, _, _, _, _, z0, zcut, eta, E0, Acut, Lambda, mu) = data
        if y <= 0:
            return mp.mpf(0)
        if y >= zcut:
            return y - z0
        L = eta / (zcut - y)
        Ey = mp.exp(L)
        (C0, C1) = (1 + mu - 2 * Lambda, Lambda - mu)
        K = 1 + abs(C0) + 2 * Lambda * abs(C1)
        if Ey > 50:
            tau = mp.exp(E0 - Ey)
            h = (1 - C0 * tau + mp.sqrt((1 - C0 * tau) ** 2 - 4 * Lambda * C1 * tau * tau)) / 2
            ell = mp.log1p(mp.log(h) / Ey)
            z = y + eta * ell / (L * (L + ell))
            return z - z0
        V = E0 * mp.expm1(eta * y / ((zcut - y) * zcut))
        w = mp.expm1(V)
        coeff = mu - w
        disc = mp.sqrt(coeff * coeff + 4 * Lambda * w)
        u = 2 * w / (disc + coeff) if coeff >= 0 else (w - mu + disc) / (2 * Lambda)
        x = mp.log1p(u)
        ell = mp.log1p(x / Acut)
        L0 = eta / (zcut - z0)
        return eta * ell / (L0 * (L0 + ell))

    def _height_g_inverse_delta(self, delta_z, data):
        """g^{-1}(z0+delta_z), rationalized at the inner endpoint."""
        mp = self.ops
        (kappa, _, _, B, C, _, _, _, _, _, _, _) = data
        if delta_z <= 0:
            return mp.mpf(0)
        coeff = B + C - delta_z
        disc = mp.sqrt(coeff * coeff + 4 * B * delta_z)
        xi = 2 * delta_z / (disc + coeff) if coeff >= 0 else (-coeff + disc) / (2 * B)
        return mp.log1p(xi) / kappa

    def inverse_height_moment(self, a, q, lam):
        """Terminal-cut elementary moment f(a,q,lam), |a|<=1/2."""
        mp = self.ops
        (a, q, lam) = (mp.mpf(a), mp.mpf(q), mp.mpf(lam))
        if abs(a) > mp.mpf('.5') or not 0 <= q <= mp.mpf('.125') or (not 0 <= lam <= 1) or (lam > 0 and q < mp.mpf(1) / 2500):
            raise ValueError('inverse moment is outside the certified (q,lam) domain')
        data = self._terminal_height_data(q, lam)
        (kappa, amp, d, B, C, z0, _, _, _, _, _, _) = data
        aa = abs(a)
        if lam == 1 or aa <= mp.mpf('.1'):
            return self._height_G(a, amp, d)
        if aa >= mp.mpf('.4'):
            return mp.sign(a) * (B * mp.exp(kappa * aa) - C * mp.exp(-kappa * aa))
        u = self._height_R(aa)
        delta_z = B * mp.expm1(kappa * u) - C * mp.expm1(-kappa * u)
        v = self._height_terminal_forward(delta_z, data)
        return mp.sign(a) * self._height_G(self._height_R(self._height_G(v, amp, d, inverse=True), inverse=True), amp, d)

    def elementary_forward_moment(self, value, q, lam):
        """Finite elementary inverse of inverse_height_moment on its range."""
        mp = self.ops
        (value, q, lam) = (mp.mpf(value), mp.mpf(q), mp.mpf(lam))
        if abs(value) > 1 - q or not 0 <= q <= mp.mpf('.125') or (not 0 <= lam <= 1) or (lam > 0 and q < mp.mpf(1) / 2500):
            raise ValueError('forward moment is outside the certified (q,lam) domain')
        data = self._terminal_height_data(q, lam)
        (_, amp, d, _, _, _, _, _, _, _, _, _) = data
        target = abs(value)
        if lam == 1 or target <= self._height_G(mp.mpf('.1'), amp, d):
            return mp.sign(value) * self._height_G(target, amp, d, inverse=True)
        v = self._height_G(target, amp, d, inverse=True)
        u = self._height_R(v)
        y = self._height_G(u, amp, d)
        delta_z = self._height_terminal_inverse_delta(y, data)
        raw = self._height_g_inverse_delta(delta_z, data)
        return mp.sign(value) * self._height_R(raw, inverse=True)

    def elementary_cusp_logradius(self, a, t, nu=1, dps=50):
        """Elementary inverse cusp height with terminal-cut two-way moment.

        The output has rho(a+1)=rho(a)+P, with all jets matching at cell seams.
        No special function, implicit solve, series or quadrature is evaluated.
        """
        mp = self.ops
        with mp.workdps(max(mp.mp.dps, dps) + 50):
            (a, t, nu) = (mp.mpf(a), mp.mpc(t), mp.mpf(nu))
            r = abs(t)
            if not 0 < r <= mp.mpf('.5') or not 0 < nu <= 1:
                raise ValueError('elementary cusp inverse requires 0<|t|<=.5 and 0<nu<=1')
            P = mp.log(nu) + mp.log(r) - 2 * r * r / (1 - r * r)
            q = mp.exp(2 * P)
            n = mp.floor(a + mp.mpf('.5'))
            b = a - n
            lam = self.scalar_flat_step((r - mp.mpf('.3')) / mp.mpf('.2'))
            if lam == 1:
                return +(P * (a + mp.mpf('.5')))
            f = self.inverse_height_moment(b, q, lam)
            root = mp.sqrt(f * f + 4 * q)
            Y = 2 * q / (root + f) if f >= 0 else (root - f) / 2
            return +(n * P + mp.log(Y) / 2)

    def elementary_cusp_height(self, rho, t, nu=1, dps=50):
        """Finite elementary forward height for t!=0, inverse to logradius."""
        mp = self.ops
        with mp.workdps(max(mp.mp.dps, dps) + 50):
            (rho, t, nu) = (mp.mpf(rho), mp.mpc(t), mp.mpf(nu))
            radius = abs(t)
            if not 0 < radius <= mp.mpf('.5') or not 0 < nu <= 1:
                raise ValueError('elementary cusp height requires 0<|t|<=.5 and 0<nu<=1')
            P = mp.log(nu) + mp.log(radius) - 2 * radius * radius / (1 - radius * radius)
            lam = self.scalar_flat_step((radius - mp.mpf('.3')) / mp.mpf('.2'))
            if lam == 1:
                return +(rho / P - mp.mpf('.5'))
            q = mp.exp(2 * P)
            n = mp.floor(rho / P)
            Y = mp.exp(2 * (rho - n * P))
            value = q / Y - Y
            return +(n + self.elementary_forward_moment(value, q, lam))

    def elementary_cusp_height_chart(self, x, y, nu=1, dps=50):
        """Source-coordinate forward height on a chart, including the two axes.

        On an axis this uses the U_0 core domain |a|<.6, or squared nonzero
        coordinate <exp(.2) in the normalized radial coordinate. At the node,
        a=0 is the smooth local value of (11.7). The global map uses its
        pole-cap formula there before any phase or Hopf normalization.
        """
        mp = self.ops
        with mp.workdps(max(mp.mp.dps, dps) + 50):
            (x, y, nu) = (mp.mpc(x), mp.mpc(y), mp.mpf(nu))
            if not 0 < nu <= 1:
                raise ValueError('elementary cusp height requires 0<nu<=1')
            if x != 0 and y != 0:
                return +self.elementary_cusp_height(mp.log(abs(y)), x * y, nu, dps)
            if x == 0 and y == 0:
                return mp.mpf(0)
            radius = nu * nu * norm_squared(x) if y == 0 else norm_squared(y)
            if radius >= mp.exp(mp.mpf('.2')):
                raise ValueError('axis height chart requires |a|<.6')
            if radius > 1:
                aa = mp.mpf('.5') + mp.log(radius) / 2
            else:
                aa = self.elementary_forward_moment(radius, mp.mpf(0), mp.mpf(0))
            return +(aa if y == 0 else -aa)

    def elementary_cusp_reference(self, t, a, nu=1):
        """Smooth core reference in the nonzero-coordinate gauge, including axes.

        |t|<=.3, |a|<.6, excluding the node. The neighbouring height cells give
        explicit open chart extensions across the central-circle representative
        seam; none of these formulas divides |t|^2 by a vanishing coordinate.
        """
        mp = self.ops
        (t, a, nu) = (mp.mpc(t), mp.mpf(a), mp.mpf(nu))
        if abs(t) > mp.mpf('.3') or abs(a) >= mp.mpf('.6'):
            raise ValueError('core reference requires |t|<=.3 and |a|<.6')
        e = mp.exp(-4 * norm_squared(t) / (1 - norm_squared(t)))
        q = e * nu * nu * norm_squared(t)
        cell = int(mp.floor(a + mp.mpf('.5')))
        b = a - cell
        f = self.inverse_height_moment(b, q, 0)
        R = mp.sqrt(f * f + 4 * q)
        if a > 0:
            A = 2 / (R - f) if cell == 1 else (R + f) / 2
            x = mp.sqrt(A / e) / nu
            return (x, t / x)
        Y = 2 / (R + f) if cell == -1 else (R - f) / 2
        y = mp.sqrt(Y)
        if y == 0:
            raise ValueError('the core reference excludes the node')
        return (t / y, y)

    def scalar_lattice_coords(self, z, w, t, logt):
        """Fixed-size real Cramer's rule, matching tcore.lattice_coords."""
        mp = self.ops
        ez = (-mp.j * t, -t, -mp.j * logt)
        ew = (1 - t * t, mp.j * (1 + t * t), 2 * t)
        rows = [[v.imag for v in ez], [v.real for v in ew], [v.imag for v in ew]]
        rhs = (z.imag, w.real, w.imag)

        def det(a):
            return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
        den = det(rows)
        coords = []
        for k in range(3):
            replaced = [[rhs[i] if j == k else rows[i][j] for j in range(3)] for i in range(3)]
            coords.append(det(replaced) / den)
        one = (z.real - sum((a * b.real for (a, b) in zip(coords, ez)))) / (2 * mp.pi)
        return (one, *coords)

    def base_reference(self, v, X, dps=50):
        """Closed reference over an UNMODIFIED algebraic-profile orbit map.

        Input is (v,X) in the S4 circle quotient, |v|+X^2=1, |X|<1.
        Uses programme 11's elementary inverse-first height.
        Returns (kind,data) for the embedding.
        Reference phases are
        local gauges, not a global section. Pole smoothing is a separate block.
        """
        mp = self.ops
        with mp.workdps(max(mp.mp.dps, dps) + 30):
            vv = tuple(map(mp.mpf, v))
            xx = mp.mpf(X)
            if not abs(xx) < 1:
                raise ValueError('the raw reference excludes the two fixed points')
            aa = xx * mp.sqrt(2 - xx * xx)
            w = (vv[0] / (1 - aa), -vv[2] / (1 - aa), -vv[1] / (1 - aa))
            y = (w[0], w[1] + 1, w[2])
            s2 = sum((z * z for z in y))
            u1 = mp.mpc(1 - s2, 2 * y[0]) / (1 + s2)
            u2 = mp.mpc(2 * y[1], 2 * y[2]) / (1 + s2)
            (A, B) = ((4 * mp.pi) ** 2, mp.exp(4))
            inverse_ratio = abs(u1) > abs(u2)
            ratio = u2 / u1 if inverse_ratio else u1 / u2
            rr = norm_squared(ratio)
            b = A * (1 - rr)
            root = mp.sqrt(b * b + 4 * B * rr)
            small = ratio * mp.sqrt((2 * B if inverse_ratio else 2) / (b + root))
            inf = inverse_ratio and abs(small) < mp.mpf('.5')
            t = None if inf else 1 / small if inverse_ratio else small
            if inf or abs(t) < mp.mpf('.5'):
                tc = small if inf else t
                phase_coordinate = -mp.conj(u1) if inf else u2
                if phase_coordinate.real < -mp.mpf('.9') * abs(phase_coordinate):
                    a = mp.mpf('.5') + mp.arg(-phase_coordinate) / (2 * mp.pi)
                else:
                    a = mp.arg(phase_coordinate) / (2 * mp.pi)
                nu = mp.exp(-2) if inf else mp.mpf(1)
                if abs(tc) <= mp.mpf('.3'):
                    (xs, yy) = self.elementary_cusp_reference(tc, a, nu)
                    return ('ci' if inf else 'c0', (xs, yy, a))
                rho = self.elementary_cusp_logradius(a, tc, nu, dps=dps)
                yy = mp.exp(rho)
                xs = tc / yy
                if abs(tc) <= mp.mpf('.3'):
                    data = (xs, mp.mpc(yy), a)
                    return ('ci' if inf else 'c0', data)
                if inf:
                    t = 1 / tc
                lt = mp.log(t)
                if inf:
                    (zm, wm) = (-mp.j * rho, mp.mpc(0))
                else:
                    ss = lt / (2 * mp.j * mp.pi)
                    (zm, wm) = (-mp.j * rho + mp.pi * ss * (ss - 1), 2 * t * ss)
                lattice = self.scalar_lattice_coords(zm, wm, t, lt)
                data = (t, lattice[0], lattice[3], lt, a)
                return ('mid', data)
            r = abs(t)
            lt = mp.log(t)
            ss = lt / (2 * mp.j * mp.pi)
            xi4 = self.scalar_lattice_coords(mp.pi * ss * (ss - 1), 2 * t * ss, t, lt)[3]
            if r <= 1:
                four = mp.arg(u2) / (2 * mp.pi) + xi4 + mp.mpf('.5')
            else:
                lam = 1 - self.scalar_flat_step((r - mp.mpf('1.2')) / mp.mpf('.4'))
                delta = mp.arg(t) / (2 * mp.pi) - xi4
                a2 = -mp.arg(u1) / (2 * mp.pi) + lam * delta + (1 - lam) / 2
                four = -a2 - mp.mpf('.5')
            return ('mid', (t, mp.mpf(0), four, lt))

    def bundle_embedding(self, kind, data, dps=50, with_norm=False):
        """Ten-component scalar embedding with the known elementary height.

        mid data=(t,x1,x4,logt); c0/ci data=(x,y), |xy|<=.3, |H|<=.6.
        Cusp/collar data include their known height. In middle data
        this height refers to the principal-log cusp representative; normalize
        logt before converting to that representative, preserving the Y point.
        Values are in S19, with the source circle acting by scalar multiplication.
        This has the same finite chart prescription as tlift.bundle_embedding_*.
        """
        mp = self.ops
        with mp.workdps(max(mp.mp.dps, dps) + 20):

            def hw(h):
                return 1 - self.scalar_flat_step((abs(h) - mp.mpf('.25')) / mp.mpf('.15'))

            def frame(x, y, nu):
                r = nu * nu * norm_squared(x) + norm_squared(y)
                if not r > 0:
                    raise ValueError('the bundle embedding excludes the fixed points')
                return (nu * mp.conj(x) / mp.sqrt(r), y / mp.sqrt(r))

            def torus_pair(z, t, nu, known_height=None):
                rr = norm_squared(t)
                m = nu * mp.exp(-2 * rr / (1 - rr))
                P = mp.log(m * abs(t))
                H = known_height
                if H is None:
                    raise ValueError('elementary projector requires the known inverse-first height')
                out = []
                for centre in (mp.mpf(0), mp.mpf('.5')):
                    shift = mp.floor(centre + mp.mpf('.5') - H)
                    w = hw(H + shift - centre)
                    if w == 0:
                        out.extend([mp.mpc(0), mp.mpc(0)])
                    else:
                        yy = z * (m * t) ** int(shift)
                        out.extend((w * q for q in frame(t / yy, yy, nu)))
                return out
            vec = [mp.mpc(0) for _ in range(10)]
            if kind in ('c0', 'ci'):
                (x, y) = map(mp.mpc, data[:2])
                nu = mp.exp(-2) if kind == 'ci' else mp.mpf(1)
                H = mp.mpf(data[2])
                if abs(x * y) > mp.mpf('.3') or abs(H) > mp.mpf('.6'):
                    raise ValueError('pure cusp embedding requires |xy|<=.3 and |H|<=.6')
                f = frame(x, y, nu)
                block = [hw(H) * q for q in f]
                w = hw(H + mp.mpf('.5') if H < 0 else H - mp.mpf('.5'))
                if H < 0 and w > 0:
                    rr = norm_squared(x * y)
                    m = nu * mp.exp(-2 * rr / (1 - rr))
                    f = frame(1 / (m * y), m * x * y * y, nu)
                block.extend((w * q for q in f))
                start = 4 if kind == 'ci' else 0
                vec[start:start + 4] = block
            elif kind == 'mid':
                (t, one, four, lt) = data[:4]
                known_height = mp.mpf(data[4]) if len(data) > 4 else None
                (t, lt, one, four) = (mp.mpc(t), mp.mpc(lt), mp.mpf(one), mp.mpf(four))
                principal = mp.log(t)
                branch = mp.nint((lt.imag - principal.imag) / (2 * mp.pi))
                (one, lt) = (one + branch * four, principal)
                r = abs(t)
                b0 = 1 - self.scalar_flat_step((r - mp.mpf('.3')) / mp.mpf('.2'))
                bi = 1 - self.scalar_flat_step((1 / r - mp.mpf('.3')) / mp.mpf('.2'))
                bm = 1 - b0 - bi
                (z, w) = (2 * mp.pi * one - mp.j * lt * four, 2 * t * four)

                def coeffs(rhs, e2, e3):
                    det = e2.real * e3.imag - e3.real * e2.imag
                    return ((rhs.real * e3.imag - e3.real * rhs.imag) / det, (e2.real * rhs.imag - rhs.real * e2.imag) / det)
                if b0 > 0:
                    s = lt / (2 * mp.j * mp.pi)
                    (zc, wc) = (z - mp.pi * s * (s - 1), w - 2 * t * s)
                    (c2, c3) = coeffs(-wc, 1 - t * t, mp.j * (1 + t * t))
                    zz = mp.exp(mp.j * (zc - mp.j * c2 * t - c3 * t))
                    vec[:4] = [b0 * q for q in torus_pair(zz, t, mp.mpf(1), known_height)]
                if bi > 0:
                    s = 1 / t
                    (zc, wc) = (z - mp.j * w / t, w / (t * t))
                    (c2, c3) = coeffs(-wc, s * s - 1, mp.j * (1 + s * s))
                    zz = mp.exp(mp.j * (zc - mp.j * c2 * s + c3 * s))
                    vec[4:8] = [bi * q for q in torus_pair(zz, s, mp.exp(-2), known_height)]
                if bm > 0:
                    for (index, centre) in enumerate((mp.mpf(0), mp.pi)):
                        angle = mp.arg(t if index == 0 else -t)
                        weight = 1 - self.scalar_flat_step((abs(angle) - mp.pi / 2) / (mp.pi / 4))
                        vec[8 + index] = bm * weight * mp.exp(mp.j * (2 * mp.pi * one - (angle + centre - lt.imag) * four))
            else:
                raise ValueError('kind must be mid, c0, or ci')
            norm2 = sum((norm_squared(q) for q in vec))
            result = tuple((+(q / mp.sqrt(norm2)) for q in vec))
            return (result, +norm2) if with_norm else result

    def _elementary_pole_inclusion(self, xi, infinity=False):
        """I_j xi: U_j^* xi in the indicated C2 block of C10."""
        mp = self.ops
        (a, b) = xi
        pair = ((a + b) / mp.sqrt(2), (a - b) / mp.sqrt(2)) if infinity else ((a - mp.j * b) / mp.sqrt(2), (a + mp.j * b) / mp.sqrt(2))
        result = [mp.mpc(0)] * 10
        start = 4 if infinity else 0
        result[start:start + 2] = pair
        return tuple(result)

    def pole_slice(self, Z, X, *, dps=50):
        mp = self.ops
        with mp.workdps(max(mp.mp.dps, dps) + 20):
            Z = tuple(map(mp.mpc, Z))
            X = mp.mpf(X)
            inf = X > 0
            scale = mp.pi if inf else 2 * mp.pi
            den = mp.sqrt(scale * (1 + abs(X) * mp.sqrt(2 - X * X)))
            pair = self._elementary_pole_inclusion(Z, inf)[4:6] if inf else self._elementary_pole_inclusion(Z)[:2]
            zz = self.elementary_pole_patch(tuple((a / den for a in pair)), inf, dps=dps)
            nu = mp.exp(-2) if inf else mp.mpf(1)
            return ('ci' if inf else 'c0', (mp.conj(zz[0]) / nu, zz[1]))

    def scalar_finite_height(self, x, y, nu=1):
        """Four-seed H_f for the bounded representatives used by the closed G.

        The more general torch implementation is separate. Here A,Y<=4 suffices
        for the exact q^2 tail; an explicit assertion guards this formula's domain.
        """
        mp = self.ops
        rr = norm_squared(x * y)
        e = mp.exp(-4 * rr / (1 - rr))
        (A, Y) = (nu * nu * e * norm_squared(x), norm_squared(y))
        q = A * Y
        if not rr <= mp.mpf('.25') or max(A, Y) > 4 or q >= mp.mpf('.125'):
            raise ValueError('four-seed finite height is outside its stated domain')

        def seed(u):
            reciprocal = u > 2
            v = 1 / u if reciprocal else u
            if v <= mp.mpf(1) / 16:
                value = v
            elif v >= mp.mpf('.5'):
                value = v / (1 + v)
            else:
                weight = self.scalar_flat_step(mp.log(16 * v) / mp.log(8))
                value = v + weight * (v / (1 + v) - v)
            return 1 - value if reciprocal else value
        return seed(A) + seed(q * A) - seed(Y) - seed(q * Y) + q * q * (A - Y) / (1 - q)

    def elementary_restore_torus(self, kind, data, phi):
        """Closed tau^{-1}, using finite H_f, in a holomorphic Y chart.

        Output (kind,(z1,z2,z3)); the real Jacobian gives the coframe in (8.14).
        The inverse-first base height is not needed here.
        """
        mp = self.ops
        phi = tuple(map(mp.mpf, phi))

        def coeffs(rhs, e2, e3):
            det = e2.real * e3.imag - e3.real * e2.imag
            return ((rhs.real * e3.imag - e3.real * rhs.imag) / det, (e2.real * rhs.imag - rhs.real * e2.imag) / det)

        def cusp_data(t, inf):
            return (t * t - 1, mp.j * (1 + t * t), -2 * t, t, mp.j * t) if inf else (1 - t * t, mp.j * (1 + t * t), 2 * t, t, -mp.j * t)
        if kind in ('c0', 'ci'):
            (x, y) = data[:2]
            inf = kind == 'ci'
            t = x * y
            (e2, e3, e4, a2, a3) = cusp_data(t, inf)
            c = coeffs(e4, e2, e3)
            H = self.scalar_finite_height(x, y, mp.exp(-2) if inf else mp.mpf(1))
            C = tuple((v + H * w for (v, w) in zip(phi, c)))
            A = C[0] * a2 + C[1] * a3
            return (kind, (x * mp.exp(-A), y * mp.exp(A), C[0] * e2 + C[1] * e3))
        if kind != 'mid':
            raise ValueError('kind must be mid, c0 or ci')
        (t, one, four, lt) = data[:4]
        r = abs(t)
        (z, w) = (2 * mp.pi * one - mp.j * lt * four, 2 * t * four)
        C = list(phi)
        for (inf, tc, beta) in ((False, t, 1 - self.scalar_flat_step((r - mp.mpf('.3')) / mp.mpf('.2'))), (True, 1 / t, 1 - self.scalar_flat_step((1 / r - mp.mpf('.3')) / mp.mpf('.2')))):
            if beta == 0:
                continue
            (e2, e3, e4, a2, a3) = cusp_data(tc, inf)
            c = coeffs(e4, e2, e3)
            if inf:
                (zc, wc) = (z - mp.j * w / t, w / (t * t))
            else:
                s = lt / (2 * mp.j * mp.pi)
                (zc, wc) = (z - mp.pi * s * (s - 1), w - 2 * t * s)
            cc = coeffs(-wc, e2, e3)
            yy = mp.exp(mp.j * (zc - mp.j * (cc[0] * a2 + cc[1] * a3)))
            H = self.scalar_finite_height(tc / yy, yy, mp.exp(-2) if inf else mp.mpf(1))
            if inf:
                D = tuple(((H + four) * v for v in c))
            else:
                s = lt / (2 * mp.j * mp.pi)
                xi = self.scalar_lattice_coords(mp.pi * s * (s - 1), 2 * t * s, t, lt)
                D = tuple((xi[j + 1] + (H - four + xi[3]) * c[j] for j in range(2)))
            C = [v + beta * d for (v, d) in zip(C, D)]
        return ('mid', (z - mp.j * t * C[0] - t * C[1], w + (1 - t * t) * C[0] + mp.j * (1 + t * t) * C[1], t))

    def _elementary_real_coeffs(self, rhs, e2, e3):
        """The two real coefficients of a complex vector in a real basis."""
        mp = self.ops
        det = e2.real * e3.imag - e3.real * e2.imag
        if det == 0:
            raise ValueError('degenerate cusp torus basis')
        return ((rhs.real * e3.imag - e3.real * rhs.imag) / det, (e2.real * rhs.imag - rhs.real * e2.imag) / det)

    def _elementary_cusp_basis(self, t, inf):
        mp = self.ops
        return (t * t - 1, mp.j * (1 + t * t), -2 * t, t, mp.j * t) if inf else (1 - t * t, mp.j * (1 + t * t), 2 * t, t, -mp.j * t)

    def _elementary_reduce_cusp(self, x, y, nu):
        """Node-adapted deck representative and original deck exponent."""
        mp = self.ops
        t = x * y
        if t != 0:
            radius2 = norm_squared(t)
            h = nu * mp.exp(-2 * radius2 / (1 - radius2)) * t
            n = int(mp.floor(mp.log(abs(y)) / mp.log(abs(h))))
            factor = h ** n
            return (x * factor, y / factor, n)
        if x == 0 and y != 0 and (abs(y) > 1):
            return (1 / (nu * y), mp.mpc(0), -1)
        if y == 0 and x != 0 and (nu * abs(x) > 1):
            return (mp.mpc(0), 1 / (nu * x), 1)
        return (x, y, 0)

    def _elementary_source_coordinates(self, kind, data, dps):
        """Undo tau_f^-1: return a W representative and its real torus angles."""
        mp = self.ops
        if kind in ('c0', 'ci'):
            (x, y, w) = map(mp.mpc, data)
            inf = kind == 'ci'
            nu = mp.exp(-2) if inf else mp.mpf(1)
            t = x * y
            if abs(t) > mp.mpf('.3'):
                raise ValueError('cusp source chart requires |xy|<=.3')
            (e2, e3, e4, a2, a3) = self._elementary_cusp_basis(t, inf)
            C = self._elementary_real_coeffs(w, e2, e3)
            A = C[0] * a2 + C[1] * a3
            (xs, ys) = (x * mp.exp(A), y * mp.exp(-A))
            c = self._elementary_real_coeffs(e4, e2, e3)
            (x0, y0, deck) = self._elementary_reduce_cusp(xs, ys, nu)
            H = self.scalar_finite_height(x0, y0, nu) + deck
            phi = tuple((C[j] - H * c[j] for j in range(2)))
            a = self.elementary_cusp_height_chart(x0, y0, nu, dps=dps)
            return (kind, (x0, y0, a), phi)
        if kind != 'mid':
            raise ValueError('source chart must be c0, ci, or mid')
        (z, w, t) = map(mp.mpc, data)
        if t == 0:
            raise ValueError('middle source chart requires t!=0')
        lt = mp.log(t)
        (one, two, three, four) = self.scalar_lattice_coords(z, w, t, lt)
        phi = [two, three]
        radius = abs(t)
        b0 = 1 - self.scalar_flat_step((radius - mp.mpf('.3')) / mp.mpf('.2'))
        bi = 1 - self.scalar_flat_step((1 / radius - mp.mpf('.3')) / mp.mpf('.2'))
        (zw, ww) = (2 * mp.pi * one - mp.j * lt * four, 2 * t * four)
        height = None
        cusp_candidate = None
        for (inf, tc, beta) in ((False, t, b0), (True, 1 / t, bi)):
            if beta == 0:
                continue
            (e2, e3, e4, a2, a3) = self._elementary_cusp_basis(tc, inf)
            c = self._elementary_real_coeffs(e4, e2, e3)
            if inf:
                (zc, wc) = (zw - mp.j * ww / t, ww / (t * t))
            else:
                s = lt / (2 * mp.j * mp.pi)
                (zc, wc) = (zw - mp.pi * s * (s - 1), ww - 2 * t * s)
            cc = self._elementary_real_coeffs(-wc, e2, e3)
            yy = mp.exp(mp.j * (zc - mp.j * (cc[0] * a2 + cc[1] * a3)))
            nu = mp.exp(-2) if inf else mp.mpf(1)
            (x0, y0, deck) = self._elementary_reduce_cusp(tc / yy, yy, nu)
            if abs(tc) <= mp.mpf('.3'):
                cusp_candidate = (inf, x0, y0)
            H = self.scalar_finite_height(x0, y0, nu) + deck
            height = self.elementary_cusp_height(mp.log(abs(yy)), tc, nu, dps=dps)
            if inf:
                D = tuple(((H + four) * v for v in c))
            else:
                s = lt / (2 * mp.j * mp.pi)
                xi = self.scalar_lattice_coords(mp.pi * s * (s - 1), 2 * t * s, t, lt)
                D = tuple((xi[j + 1] + (H - four + xi[3]) * c[j] for j in range(2)))
            phi = [phi[j] - beta * D[j] for j in range(2)]
        if cusp_candidate is not None:
            (inf, x0, y0) = cusp_candidate
            nu = mp.exp(-2) if inf else mp.mpf(1)
            a = self.elementary_cusp_height_chart(x0, y0, nu, dps=dps)
            return ('ci' if inf else 'c0', (x0, y0, a), tuple(phi))
        wdata = (t, one, four, lt) if height is None else (t, one, four, lt, height)
        return (kind, wdata, tuple(phi))

    def _elementary_raw_orbit(self, kind, data):
        """Algebraic-profile raw W/S1 -> (v,X), before pole conjugation."""
        mp = self.ops
        inf = kind == 'ci'
        if kind == 'mid':
            (t, one, four, lt) = data[:4]
            inf = abs(t) > 1
            tc = 1 / t if inf else t
            if len(data) > 4:
                a = -data[4] + mp.mpf('.5') if inf else data[4]
            else:
                s = lt / (2 * mp.j * mp.pi)
                xi4 = self.scalar_lattice_coords(mp.pi * s * (s - 1), 2 * t * s, t, lt)[3]
                if inf:
                    lam = 1 - self.scalar_flat_step((abs(t) - mp.mpf('1.2')) / mp.mpf('.4'))
                    delta = mp.arg(t) / (2 * mp.pi) - xi4
                    a = four + mp.mpf('.5') + lam * delta + (1 - lam) / 2
                else:
                    a = four - xi4 - mp.mpf('.5')
        else:
            (x, y, a) = data
            tc = x * y
        (A, B) = ((4 * mp.pi) ** 2, mp.exp(4))
        if inf:
            ratio = tc * mp.sqrt((norm_squared(tc) + A) / (A * norm_squared(tc) + B))
            u1 = mp.exp(2 * mp.pi * mp.j * (a if kind == 'mid' else -a + mp.mpf('.5'))) / mp.sqrt(1 + norm_squared(ratio))
            u2 = ratio * u1
        else:
            ratio = tc * mp.sqrt((A + B * norm_squared(tc)) / (1 + A * norm_squared(tc)))
            u2 = mp.exp(2 * mp.pi * mp.j * a) / mp.sqrt(1 + norm_squared(ratio))
            u1 = ratio * u2
        den = 1 + u1.real
        if den <= 0:
            raise ValueError('raw orbit reached the north stereographic pole')
        w = (u1.imag / den, u2.real / den - 1, u2.imag / den)
        r2 = sum((v * v for v in w))
        v = (2 * w[0] / (1 + r2), -2 * w[2] / (1 + r2), -2 * w[1] / (1 + r2))
        X = (mp.sqrt(r2) - 1) / mp.sqrt(1 + r2)
        return (v, X)

    def _elementary_hopf_reference(self, n):
        mp = self.ops
        if n[2] >= 0:
            w = mp.sqrt((1 + n[2]) / 2)
            return (w, mp.mpc(n[0], -n[1]) / (2 * w))
        w = mp.sqrt((1 - n[2]) / 2)
        return (mp.mpc(n[0], n[1]) / (2 * w), w)

    def middle_lift(self, u1, u2, t, Ps, phase_offset=None):
        """The candidate of (12.11), with the C61 collar offset when supplied.

        Inputs use the existing raw base and its exact source projector.
        Returns the normalized 10x2 map, its unnormalized norm square,
        and the corrected overlap real part when both local maps are used.
        Numerical evaluation uses mpmath at the caller's working precision.
        """
        mp = self.ops
        (a, b) = (abs(u1), abs(u2))
        p = mp.matrix([u1.imag, -u2.imag, 1 + u1.real - u2.real])
        n = p / mp.norm(p)
        Pt = mp.matrix([[1 + n[2], n[0] + mp.j * n[1]], [n[0] - mp.j * n[1], 1 - n[2]]]) / 2
        th = mp.arg(t)
        chi = self.scalar_flat_step((abs(th) - 3 * mp.pi / 8) / (mp.pi / 4))

        def local(slot, frame):
            B = Ps[:, slot] * frame.H * Pt
            return B / mp.sqrt(sum((norm_squared(z) for z in B)))
        out = mp.zeros(10, 2)
        gamma_real = None
        if chi < 1:
            L0 = local(8, mp.matrix([1, 0]))
            out += (1 - chi) * L0
        if chi > 0:
            z = u2 / b
            v = mp.matrix([mp.j * (a + mp.j * b) * (z - 1) / 2, (z + 1) / 2])
            th1 = mp.arg(-t) + mp.pi
            log1 = mp.log(abs(t)) + mp.j * th1
            s = log1 / (2 * mp.j * mp.pi)
            xi4 = self.scalar_lattice_coords(mp.pi * s * (s - 1), 2 * t * s, t, log1)[3]
            lam = 1 - self.scalar_flat_step((abs(t) - mp.mpf('1.2')) / mp.mpf('.4'))
            C = (1 - lam) * th1 / (2 * mp.pi) + lam * xi4 + lam / 2
            if phase_offset is not None:
                C = phase_offset(th1)
            w = self.scalar_flat_step((th1 - 5 * mp.pi / 8) / (3 * mp.pi / 4))
            k = mp.j * (a - mp.j * b)
            d = mp.conj(k) * mp.exp(mp.j * w * (mp.pi - 2 * mp.pi * C))
            L1 = d * local(9, v)
            out += chi * L1
            if chi < 1:
                gamma_real = mp.re(sum((mp.conj(v) * w for (v, w) in zip(L1, L0))))
        normsq = mp.re(sum((norm_squared(z) for z in out)))
        return (out / mp.sqrt(normsq), normsq, gamma_real)

    def cusp_lift(self, u1, u2, tc, Ps, infinity=False, pole_coordinate=None):
        """Finite regional lift on |tc|<1/2, excluding the node.

        All source phases come from the actual projector Ps. No new source
        section is assumed. Optional pole_coordinate is the raw round X;
        C60 then matches the minimal raw pole lift. The actual pole replacement
        remains the existing elementary_pole_patch, not this matrix function.
        """
        mp = self.ops
        p = mp.matrix([u1.imag, -u2.imag, 1 + u1.real - u2.real])
        n = p / mp.norm(p)
        Pt = mp.matrix([[1 + n[2], n[0] + mp.j * n[1]], [n[0] - mp.j * n[1], 1 - n[2]]]) / 2
        phase = -mp.conj(u1) if infinity else u2
        H = mp.arg(phase) / (2 * mp.pi)
        H1 = mp.mpf('.5') + mp.arg(-phase) / (2 * mp.pi)
        chi = self.scalar_flat_step((abs(H) - mp.mpf('.3')) / mp.mpf('.05'))
        U = (mp.matrix([[1, 1], [1, -1]]) if infinity else mp.matrix([[1, 1], [mp.j, -mp.j]])) / mp.sqrt(2)
        start = 4 if infinity else 0
        out = mp.zeros(10, 2)
        gamma = None
        norms = []
        if chi < 1:
            d = mp.exp(-2 * norm_squared(tc) / (1 - norm_squared(tc)))
            if pole_coordinate is not None:
                distance = 1 - (pole_coordinate if infinity else -pole_coordinate)
                lam = self.scalar_flat_step((distance - mp.mpf('1e-5')) / mp.mpf('1e-5'))
                d = 1 - lam * (1 - d)
            D = mp.diag([d, 1])
            B = Ps[:, start:start + 2] * D * U.H * Pt
            q0 = sum((norm_squared(z) for z in B))
            L0 = B / mp.sqrt(q0)
            out += (1 - chi) * L0
            norms.append(q0)
        if chi > 0:
            tau = mp.pi / 2 * self.scalar_flat_step((H1 - mp.mpf('.35')) / mp.mpf('.3'))
            v = U * mp.matrix([mp.cos(tau), mp.sin(tau)])
            B = Ps[:, start + 2] * v.H * Pt
            q1 = sum((norm_squared(z) for z in B))
            L1 = B / mp.sqrt(q1)
            out += chi * L1
            norms.append(q1)
            if chi < 1:
                gamma = mp.re(sum((mp.conj(v) * w for (v, w) in zip(L1, L0))))
        normsq = sum((norm_squared(z) for z in out))
        return (out / mp.sqrt(normsq), normsq, gamma, min(norms))

    def collar_offset(self, t, u2, th1):
        mp = self.ops
        r = abs(t)
        ph = mp.arg(u2)
        if r < mp.mpf('.5'):
            H = ph / (2 * mp.pi)
            P = mp.log(r) - 2 * r * r / (1 - r * r)
            rho = self.elementary_cusp_logradius(H, t, dps=60)
            log1 = mp.log(r) + mp.j * th1
            s = log1 / (2 * mp.j * mp.pi)
            xi4 = self.scalar_lattice_coords(mp.pi * s * (s - 1), 2 * t * s, t, log1)[3]
            return xi4 + rho / P - H
        if r > 2:
            rc = 1 / r
            H = (mp.pi - th1 - ph) / (2 * mp.pi)
            P = -2 + mp.log(rc) - 2 * rc * rc / (1 - rc * rc)
            rho = self.elementary_cusp_logradius(H, 1 / t, mp.exp(-2), dps=60)
            return 1 - rho / P - ph / (2 * mp.pi)
        raise ValueError('the collar offset is only for r<1/2 or r>2')

    def extended_cusp_lift(self, th, ph, rc, infinity=False):
        """Extended cusp lift with source weights cancelled; 0.49<=rc<=0.51.

        Returns the extended matrix and the independent middle matrix. Above
        rc=1/2 use rho=P(H+1/2). Below it retain the actual elementary height.
        The global circle section uses this formula on its boundary collars.
        """
        mp = self.ops
        r = 1 / rc if infinity else rc
        t = r * mp.exp(mp.j * th)
        tc = 1 / t if infinity else t
        nu = mp.exp(-2) if infinity else mp.mpf(1)
        ratio2 = r * r * ((4 * mp.pi) ** 2 + mp.exp(4) * r * r) / (1 + (4 * mp.pi) ** 2 * r * r)
        (a, b) = (mp.sqrt(ratio2 / (1 + ratio2)), 1 / mp.sqrt(1 + ratio2))
        u1 = a * mp.exp(mp.j * (th + ph))
        u2 = b * mp.exp(mp.j * ph)
        p = mp.matrix([u1.imag, -u2.imag, 1 + u1.real - u2.real])
        n = p / mp.norm(p)
        Pi = mp.matrix([[1 + n[2], n[0] + mp.j * n[1]], [n[0] - mp.j * n[1], 1 - n[2]]]) / 2
        H = mp.arg(-mp.conj(u1) if infinity else u2) / (2 * mp.pi)
        m = nu * mp.exp(-2 * rc * rc / (1 - rc * rc))
        P = mp.log(m * rc)
        rho = P * (H + mp.mpf('.5')) if rc >= mp.mpf('.5') else self.elementary_cusp_logradius(H, tc, nu, dps=60)
        y = mp.exp(rho)
        x = tc / y
        lt = mp.log(r) + mp.j * th
        s = lt / (2 * mp.j * mp.pi)
        z = -mp.j * rho + (0 if infinity else mp.pi * s * (s - 1))
        w = 0 if infinity else 2 * t * s
        (one, _, _, four) = self.scalar_lattice_coords(z, w, t, lt)
        e = mp.matrix(self.bundle_embedding('mid', (t, one, four, lt, H), dps=60))
        Ps = e * e.H
        U = (mp.matrix([[1, 1], [1, -1]]) if infinity else mp.matrix([[1, 1], [mp.j, -mp.j]])) / mp.sqrt(2)
        H1 = H if H >= 0 else H + 1
        chi = self.scalar_flat_step((abs(H) - mp.mpf('.3')) / mp.mpf('.05'))
        row = mp.zeros(1, 2)
        if chi < 1:
            F = U * mp.matrix([m * mp.conj(x), y])
            b0 = F.H * Pi
            row += (1 - chi) * b0 / mp.norm(b0)
        if chi > 0:
            tau = mp.pi / 2 * self.scalar_flat_step((H1 - mp.mpf('.35')) / mp.mpf('.3'))
            v = U * mp.matrix([mp.cos(tau), mp.sin(tau)])
            phase = x / abs(x) if H >= 0 else mp.mpf(1)
            b1 = phase * v.H * Pi
            row += chi * b1 / mp.norm(b1)
        Lc = e * row / mp.norm(row)
        offset = (lambda th1: self.collar_offset(t, u2, th1)) if rc < mp.mpf('.5') else None
        Lm = self.middle_lift(u1, u2, t, Ps, phase_offset=offset)[0]
        return (Lc, Lm, Ps, Pi, (u1, u2, tc))

    def circle_section(self, xi, X, dps=60, phase_local=None):
        mp = self.ops
        pc = 2 * xi[0] * mp.conj(xi[1])
        n = (pc.real, pc.imag, norm_squared(xi[0]) - norm_squared(xi[1]))
        vv = tuple(((1 - X * X) * v for v in n))
        (kind, data) = self.base_reference(vv, X, dps=dps)
        e = mp.matrix(self.bundle_embedding(kind, data, dps=dps))
        Ps = e * e.H
        aa = X * mp.sqrt(2 - X * X)
        yy = (vv[0] / (1 - aa), 1 - vv[2] / (1 - aa), -vv[1] / (1 - aa))
        ss = sum((v * v for v in yy))
        u1 = mp.mpc(1 - ss, 2 * yy[0]) / (1 + ss)
        u2 = mp.mpc(2 * yy[1], 2 * yy[2]) / (1 + ss)
        if kind in ('c0', 'ci'):
            tc = data[0] * data[1]
            inf = kind == 'ci'
            L = self.cusp_lift(u1, u2, tc, Ps, inf, pole_coordinate=X)[0]
        else:
            t = data[0]
            r = abs(t)
            if r < mp.mpf('.5') or r > 2:
                inf = r > 2
                tc = 1 / t if inf else t
                if abs(tc) >= mp.mpf('.49'):
                    L = self.extended_cusp_lift(mp.arg(t), mp.arg(u2), abs(tc), inf)[0]
                else:
                    L = self.cusp_lift(u1, u2, tc, Ps, inf, pole_coordinate=X)[0]
            else:
                (th, ph) = (mp.arg(t), mp.arg(u2))
                L = self.middle_lift(u1, u2, t, Ps)[0]

                def transition(rc, inf):
                    (LC, LM, *_) = self.extended_cusp_lift(th, ph, rc, inf)
                    return sum((mp.conj(v) * w for (v, w) in zip(LM, LC)))

                def retreat(rc):
                    return mp.mpf('.5') + (rc - mp.mpf('.5')) * (1 - self.scalar_flat_step((rc - mp.mpf('.5')) / mp.mpf('.01')))
                if r <= mp.mpf('.75'):
                    correction = transition(retreat(r), False)
                elif r >= mp.mpf('1.25'):
                    correction = transition(retreat(1 / r), True)
                else:
                    phase = mp.mpf(phase_value(float(th), float(ph))[0]) if phase_local is None else phase_local(th, ph)
                    weight = self.scalar_flat_step((r - mp.mpf('.75')) / mp.mpf('.5'))
                    correction = transition(mp.mpf('.5'), False) * mp.exp(mp.j * weight * phase)
                L = L * correction
        return (L * mp.matrix(xi), kind, data)

    def coordinates(self, Z, X, phi, dps=60, phase_local=None):
        """G in a selected holomorphic chart, with the fixed finite phase."""
        mp = self.ops
        with mp.workdps(max(mp.mp.dps, dps) + 20):
            X = mp.mpf(X)
            Z = tuple(map(mp.mpc, Z))
            if 1 - abs(X) < mp.mpf('5e-6'):
                (kind, data) = self.pole_slice(Z, X, dps=dps)
            else:
                xi = tuple((z / mp.sqrt(1 - X * X) for z in Z))
                (result, kind, data) = self.circle_section(xi, X, dps, phase_local)
                ref = self.bundle_embedding(kind, data, dps=dps)
                phase = sum((mp.conj(v) * w for (v, w) in zip(ref, result)))
                if kind == 'mid':
                    (t, one, four, lt) = data[:4]
                    data = (t, one + mp.arg(phase) / (2 * mp.pi), four, lt)
                else:
                    (x, y) = data[:2]
                    data = (mp.conj(phase) * x, phase * y)
            return self.elementary_restore_torus(kind, data, phi)

    def forward(self, kind, data, dps=60):
        """Explicit inverse using the same circle map and base/cap inverses."""
        mp = self.ops
        with mp.workdps(max(mp.mp.dps, dps) + 30):
            (wk, wd, phi) = self._elementary_source_coordinates(kind, data, dps)
            if wk in ('c0', 'ci'):
                inf = wk == 'ci'
                nu = mp.exp(-2) if inf else mp.mpf(1)
                (x, y) = wd[:2]
                if nu * nu * norm_squared(x) + norm_squared(y) < mp.mpf('2e-6'):
                    zz = self.elementary_pole_forward((nu * mp.conj(x), y), inf, dps=dps)
                    scale = mp.pi if inf else 2 * mp.pi
                    h = scale * sum((norm_squared(v) for v in zz))
                    X = (1 if inf else -1) * (1 - h) / mp.sqrt(1 + h * h)
                    if 1 - abs(X) < mp.mpf('5e-6'):
                        pair = ((zz[0] + zz[1]) / mp.sqrt(2), (zz[0] - zz[1]) / mp.sqrt(2)) if inf else ((zz[0] + zz[1]) / mp.sqrt(2), mp.j * (zz[0] - zz[1]) / mp.sqrt(2))
                        return (tuple((mp.sqrt(2 * scale / (1 + h * h)) * v for v in pair)), X, phi)
            (v, X) = self._elementary_raw_orbit(wk, wd)
            radius = 1 - X * X
            n = tuple((a / radius for a in v))
            xi = self._elementary_hopf_reference(n)
            (section, _, _) = self.circle_section(xi, X, dps)
            embedding = self.bundle_embedding(wk, wd, dps=dps)
            phase = sum((mp.conj(a) * b for (a, b) in zip(section, embedding)))
            if abs(abs(phase) - 1) > mp.mpf('1e-25'):
                raise ValueError('new forward map recovered the wrong source line')
            phase /= abs(phase)
            return (tuple((mp.sqrt(radius) * phase * a for a in xi)), X, phi)

    def boundary_mismatch(self, theta, phi):
        mp = self.ops
        deltas = []
        for inf in (False, True):
            (Lc, Lm, *_) = self.extended_cusp_lift(theta, phi, mp.mpf('.5'), inf)
            deltas.append(sum((mp.conj(a) * b for (a, b) in zip(Lm, Lc))))
        return deltas[1] / deltas[0]


def phase_jet(theta, phi, formulas):
    """The fixed finite phase value and its exact endpoint chain rule.

    Section 12.9 proves exp(i Theta)=g and periodic smoothness. Hence
    d Theta=Im(dg/g), d2 Theta=Im(d2g/g-dg dg/g2), checked independently
    by the focused verifier. This is a derivative rule for that finite sum.
    """
    g = formulas.boundary_mismatch(theta, phi)
    value = phase_value(float(theta), float(phi), PHASE_STEPS)[0]
    if not isinstance(g, Jet):
        return mp.mpf(value)
    return Jet(value, imag_array(g.gradient / g.value),
               imag_array(g.hessian / g.value
                          - np.outer(g.gradient, g.gradient) / g.value**2))


@dataclass(frozen=True)
class CoordinateTensor:
    """Numerical coordinate data; derivative[k,i,j] means partial_k J^i_j."""
    matrix: np.ndarray
    derivative: np.ndarray
    jacobian: np.ndarray
    hessian: np.ndarray
    nijenhuis: np.ndarray
    chart: str
    missing: int
    sign: int


def nijenhuis_components(J, dJ):
    """Coordinate formula, with the repository's index convention."""
    return np.array([[[sum(J[a, i]*dJ[a, k, j] - J[a, j]*dJ[a, k, i]
                          - J[k, a]*(dJ[i, a, j] - dJ[j, a, i])
                          for a in range(6))
                       for j in range(6)] for i in range(6)] for k in range(6)])


def coordinate_tensor(coordinates, *, missing=4, sign=1, dps=60):
    """Evaluate J and dJ in a specified round hemisphere coordinate chart.

    Inputs are the four retained sphere coordinates, then the two torus
    coordinates of period one. This function uses finite chain rules.
    Numerical differentiation and the historical certificate are not called.
    """
    if len(coordinates) != 6 or missing not in range(5) or sign not in (-1, 1):
        raise ValueError('six coordinates, missing in 0..4, sign in {-1,1} required')
    with mp.workdps(dps):
        u = [mp.mpf(v) for v in coordinates]
        if sum(v*v for v in u[:4]) >= 1:
            raise ValueError('coordinates must lie in the open hemisphere chart')
        q = [Jet(v, np.eye(6)[i]) for i, v in enumerate(u)]
        sphere = q[:4]
        sphere.insert(missing, sign*arithmetic.sqrt(1-sum(v*v for v in sphere)))
        formulas = NumericalFormulas()
        kind, values = formulas.coordinates(
            (sphere[0]+1j*sphere[1], sphere[2]+1j*sphere[3]), sphere[4], q[4:],
            dps=dps, phase_local=lambda t, p: phase_jet(t, p, formulas))
        real_values = [part for z in values for part in (z.real, z.imag)]
        A = np.array([real_array(z.gradient) for z in real_values], float)
        H = np.array([real_array(z.hessian) for z in real_values], float)
        if not (np.isfinite(A).all() and np.isfinite(H).all()):
            raise ArithmeticError('a singular intermediate reached the output derivatives')
        J0 = np.kron(np.eye(3), [[0., -1.], [1., 0.]])
        J = np.linalg.solve(A, J0 @ A)
        dJ = np.array([np.linalg.solve(A, J0 @ H[:, :, k] - H[:, :, k] @ J)
                       for k in range(6)])
        return CoordinateTensor(J, dJ, A, H, nijenhuis_components(J, dJ),
                                kind, missing, sign)


def frame_matrices(sphere, missing, *, dtype=float):
    """Coordinate-to-global frame and its inverse; sphere[missing] != 0.

    Global identification: (v,a,b) -> (v+a*x,b), for x.v=0.
    Inverse: (w,b) -> (w-(x.w)x, x.w, b).
    Torus tangent components are measured in coordinates of period one.
    """
    x = np.asarray(sphere, dtype=dtype)
    keep = [i for i in range(5) if i != missing]
    T = np.zeros((6, 6), dtype=dtype)
    for j, i in enumerate(keep):
        T[i, j] = 1
        T[missing, j] = -x[i]/x[missing]
    T[:5, 4] = x
    T[5, 5] = 1
    inverse = np.zeros((6, 6), dtype=dtype)
    inverse[:4, :5] = (np.eye(5, dtype=dtype)-np.outer(x, x))[keep, :]
    inverse[4, :5] = x
    inverse[5, 5] = 1
    return T, inverse


def frame_brackets(sphere, *, dtype=float):
    """[E_i,E_j]=x_i E_j-x_j E_i for i,j<5; [E_i,E_5]=0."""
    x = np.asarray(sphere,dtype=dtype)
    c = np.zeros((6,6,6),dtype=dtype)
    for i in range(5):
        for j in range(5):
            c[j,i,j] += x[i]
            c[i,i,j] -= x[j]
    return c


def frame_derivatives(sphere, missing, *, dtype=float):
    """The six coordinate derivatives of the coordinate-to-global matrix."""
    x = np.asarray(sphere,dtype=dtype)
    keep = [i for i in range(5) if i!=missing]
    T,_ = frame_matrices(x,missing,dtype=dtype)
    derivatives = np.zeros((6,6,6),dtype=dtype)
    for k in range(4):
        derivatives[k,missing,:4] = (-np.eye(4,dtype=dtype)[k]/x[missing]
                                    -x[keep]*x[keep[k]]/x[missing]**3)
        derivatives[k,:5,4] = T[:5,k]
    return derivatives


def global_frame_data(sphere, coordinate):
    """Global matrix, frame directional derivatives, and frame Nijenhuis.

    This exposes the non-coordinate bracket terms when checking the public
    matrix directly. It uses exact differentiation formulas for the frame
    change; numerical arrays approximate the corresponding exact expressions.
    """
    x = np.asarray(sphere,float)
    missing = coordinate.missing
    T,inverse = frame_matrices(x,missing)
    J = T @ coordinate.matrix @ inverse
    derivatives = []
    for k,dT in enumerate(frame_derivatives(x,missing)):
        derivatives.append(dT @ coordinate.matrix @ inverse
                           + T @ coordinate.derivative[k] @ inverse
                           - J @ dT @ inverse)
    directions = np.einsum('ia,ikj->akj',inverse,np.asarray(derivatives))
    c = frame_brackets(x)
    N = nijenhuis_components(J,directions)
    N += np.einsum('ai,bj,kab->kij',J,J,c)
    N -= np.einsum('kb,ai,baj->kij',J,J,c)
    N -= np.einsum('kb,aj,bia->kij',J,J,c)
    N += np.einsum('kb,bij->kij',J @ J,c)
    return J,directions,N


def complex_structure(point, *, dps=60, symbolic=False):
    """Return the 6x6 matrix of J at (x0,...,x4,t1,t2).

    sum(x_i**2)=1; t1,t2 are real representatives modulo integers.
    The output uses the smooth GLOBAL tangent identification
        (v,a,b) -> (v+a*x,b) in R5 x R,  x.v=0.
    Consequently its square is -I6 and matrices at different points use
    the same specified frame. This frame is not a coordinate frame:
    use coordinate_tensor() for the coordinate Nijenhuis calculation.

    With symbolic=True, return a SymbolicTensor describing the exact
    formula branch at this witness; its .matrix is in the same global
    frame and its graph records the named elementary assignments.
    Global joins are supplied by the proof contracts.

    Float/multiprecision inputs and the returned float64 matrix are
    approximate evaluations of the exact definition. An input is checked
    to floating-point tolerance; the selected hemisphere reconstructs
    its omitted coordinate on the sphere. No claim of uniform rounding
    accuracy is made. Global mathematical premises are in proof_contracts().
    """
    if len(point) != 7:
        raise ValueError('expected five sphere and two period-one torus coordinates')
    with mp.workdps(dps):
        p = [mp.mpf(v) for v in point]
        if not all(mp.isfinite(v) for v in p):
            raise ValueError('all coordinates must be finite and real')
        if abs(sum(v*v for v in p[:5])-1) > mp.mpf('1e-12'):
            raise ValueError('sphere coordinates must have squared norm one')
        missing = max(range(5), key=lambda i: abs(p[i]))
        sign = 1 if p[missing] > 0 else -1
        u = [p[i] for i in range(5) if i != missing]+p[5:]
        if symbolic:
            return symbolic_definition(u, missing=missing, sign=sign, dps=dps)
        result = coordinate_tensor(u, missing=missing, sign=sign, dps=dps)
        # Reconstruct identically to the coordinate evaluator, rather than
        # mixing the input rounding with a different tangent plane.
        x = p[:5]
        x[missing] = sign*mp.sqrt(1-sum(v*v for v in u[:4]))
        T, inverse = frame_matrices(x, missing)
        return T @ result.matrix @ inverse


@dataclass(frozen=True)
class ScalarAssignment:
    symbol: object
    expression: object
    differential_expression: object
    contract: str


class ExpressionGraph:
    """Real elementary assignments, with numerical witnesses only for branching.

    Expressions and constants are exact SymPy objects. A witness selects an
    open formula branch; it is never used to prove an identity or inequality.
    Global coverage and joins are the separate contracts below. The phase
    node has the finite-sum value and the certified local argument derivatives.
    """
    def __init__(self):
        import sympy as sp
        self.sp = sp
        self.assignments = []
        self.choices = []
        self.inputs = sp.symbols('q0:6', real=True)

    def real_node(self, expression, differential=None, contract='elementary'):
        expression = self.sp.sympify(expression)
        if not expression.free_symbols:
            return expression
        if isinstance(expression, self.sp.Symbol):
            return expression
        symbol = self.sp.Symbol('v'+str(len(self.assignments)), real=True)
        self.assignments.append(ScalarAssignment(
            symbol, expression, expression if differential is None else differential, contract))
        return symbol

    def number(self, value, re=None, im=None):
        if isinstance(value, SymbolicNumber):
            return value
        value = mp.mpc(value)
        if re is None:
            re = self.sp.Rational(str(value.real))
        if im is None:
            im = self.sp.Rational(str(value.imag))
        return SymbolicNumber(self, value, self.sp.sympify(re), self.sp.sympify(im))

    def emit(self, value, re, im=0):
        return self.number(value, self.real_node(re), self.real_node(im))

    def active_assignments(self, outputs):
        needed = set().union(*(e.free_symbols for e in outputs))
        active = []
        for node in reversed(self.assignments):
            if node.symbol in needed:
                active.append(node)
                needed.update(node.differential_expression.free_symbols)
                needed.update(node.expression.free_symbols)
        return tuple(reversed(active))

    def finite_phase_definition(self):
        """Definition of the named finite-sum primitive, without expansion.

        boundary_mismatch is the concrete Formulas method, not an input
        function to be chosen by the caller. All 65536 summands are specified.
        """
        sp = self.sp
        theta, phi = sp.symbols('theta phi', real=True)
        j = sp.Symbol('j', integer=True)
        g = sp.Function('boundary_mismatch')
        ratio = g(j*theta/PHASE_STEPS,j*phi/PHASE_STEPS)/g(
            (j-1)*theta/PHASE_STEPS,(j-1)*phi/PHASE_STEPS)
        h = sp.log(2)
        return -h*(sp.Rational(1,5)-h/4)/sp.pi + sp.Sum(
            sp.atan(sp.im(ratio)/sp.re(ratio)),(j,1,PHASE_STEPS))

    def differentiate(self, outputs):
        """Exact first/second derivatives as a straight-line program.

        Local partials are obtained by SymPy differentiation of each actual
        assignment. The product/chain rule keeps earlier derivatives named.
        Mixed partial symmetry follows from the separately checked chain rule;
        every local scalar Hessian is also checked here, before composition.
        """
        sp = self.sp
        gradient = {q: tuple(sp.Integer(i == k) for k in range(6))
                    for i, q in enumerate(self.inputs)}
        hessian = {q: sp.zeros(6) for q in self.inputs}
        derivative_assignments = []
        for node in self.active_assignments(outputs):
            parents = sorted(node.differential_expression.free_symbols, key=str)
            if any(p not in gradient for p in parents):
                raise ValueError('unbound symbolic dependency')
            f = node.differential_expression
            first = [sp.diff(f, p) for p in parents]
            second = [[sp.diff(first[i], p) for p in parents] for i in range(len(parents))]
            if any(sp.simplify(second[i][j]-second[j][i]) != 0
                   for i in range(len(parents)) for j in range(i)):
                raise AssertionError('a scalar primitive has inconsistent mixed derivatives')
            grad = []
            H = sp.zeros(6)
            for i in range(6):
                expression = sum(c*gradient[p][i] for c, p in zip(first, parents))
                symbol = sp.Symbol(str(node.symbol)+'d'+str(i), real=True)
                if expression == 0:
                    grad.append(sp.S.Zero)
                else:
                    derivative_assignments.append((symbol, expression))
                    grad.append(symbol)
                for j in range(i, 6):
                    expression = sum(c*hessian[p][i, j] for c, p in zip(first, parents))
                    expression += sum(second[a][b]*gradient[p][i]*gradient[q][j]
                                      for a, p in enumerate(parents)
                                      for b, q in enumerate(parents))
                    if expression == 0:
                        H[i, j] = H[j, i] = 0
                    else:
                        symbol = sp.Symbol(str(node.symbol)+'h'+str(i)+str(j), real=True)
                        derivative_assignments.append((symbol, expression))
                        H[i, j] = H[j, i] = symbol
            gradient[node.symbol], hessian[node.symbol] = tuple(grad), H
        A = sp.Matrix([gradient.get(v, (0,)*6) for v in outputs])
        H = tuple(hessian.get(v, sp.zeros(6)) for v in outputs)
        return A, H, tuple(derivative_assignments)


class SymbolicNumber:
    """A complex pair of real expression references and a branch witness."""
    def __init__(self, graph, value, re, im):
        self.graph, self.value, self.re, self.im = graph, value, re, im

    def coerce(self, other):
        return self.graph.number(other)

    def __add__(self, other):
        if isinstance(other, Matrix):
            return NotImplemented
        other = self.coerce(other)
        return self.graph.emit(self.value+other.value, self.re+other.re, self.im+other.im)
    __radd__ = __add__

    def __neg__(self):
        return self.graph.emit(-self.value, -self.re, -self.im)

    def __pos__(self):
        return self

    def __sub__(self, other):
        return self + -self.coerce(other)

    def __rsub__(self, other):
        return self.coerce(other) + -self

    def __mul__(self, other):
        if isinstance(other, Matrix):
            return NotImplemented
        other = self.coerce(other)
        return self.graph.emit(self.value*other.value,
                               self.re*other.re-self.im*other.im,
                               self.re*other.im+self.im*other.re)
    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.coerce(other)
        denominator = other.re**2+other.im**2
        return self.graph.emit(self.value/other.value,
                               (self.re*other.re+self.im*other.im)/denominator,
                               (self.im*other.re-self.re*other.im)/denominator)

    def __rtruediv__(self, other):
        return self.coerce(other)/self

    def __pow__(self, power):
        if isinstance(power, SymbolicNumber):
            return symbolic_unary('exp', power*symbolic_unary('log', self))
        if power == 0:
            return self.coerce(1)
        if power == 1:
            return self
        if power == 2:
            return self*self
        if self.im == 0 and self.value.real >= 0:
            p = self.graph.sp.Rational(str(power))
            return self.graph.emit(self.value**power, self.re**p)
        if int(power) == power:
            if power < 0:
                return 1/(self**(-power))
            result = self.coerce(1)
            for _ in range(int(power)):
                result *= self
            return result
        raise ValueError('nonreal noninteger power has no declared branch')

    @property
    def real(self):
        return self.graph.number(self.value.real, self.re, 0)

    @property
    def imag(self):
        return self.graph.number(self.value.imag, self.im, 0)

    def conjugate(self):
        return self.graph.emit(self.value.conjugate(), self.re, -self.im)

    def __abs__(self):
        return norm_squared(self)**0.5

    def __float__(self):
        if self.value.imag:
            raise TypeError('a real scalar is required')
        return float(self.value.real)

    def __int__(self):
        value = int(self.value.real)
        self.graph.choices.append(('int', self.re, value))
        return value

    def _compare(self, other, operation):
        other = self.coerce(other)
        if self.value.imag or other.value.imag:
            raise TypeError('ordered comparisons require real scalars')
        value = operation(self.value.real, other.value.real)
        self.graph.choices.append((operation.__name__, self.re-other.re, bool(value)))
        return bool(value)

    def __lt__(self, other):
        import operator
        return self._compare(other, operator.lt)

    def __le__(self, other):
        import operator
        return self._compare(other, operator.le)

    def __gt__(self, other):
        import operator
        return self._compare(other, operator.gt)

    def __ge__(self, other):
        import operator
        return self._compare(other, operator.ge)

    def __eq__(self, other):
        other = self.coerce(other)
        result = self.value == other.value
        self.graph.choices.append(('eq', (self.re-other.re, self.im-other.im), result))
        return result

    def __ne__(self, other):
        return not self == other

    def __bool__(self):
        return self != 0


def symbolic_unary(name, x):
    """Exact elementary formulas, separated into real components."""
    g, sp = x.graph, x.graph.sp
    a, b = x.re, x.im
    value = getattr(mp, name)(x.value)
    if name == 'exp':
        return g.emit(value, sp.exp(a)*sp.cos(b), sp.exp(a)*sp.sin(b))
    if name == 'arg':
        return g.emit(value, sp.atan2(b, a))
    if name == 'log':
        if b == 0 and x.value.real > 0:
            return g.emit(value, sp.log(a))
        return g.emit(value, sp.log(a*a+b*b)/2, sp.atan2(b, a))
    if name == 'expm1':
        return g.emit(value, sp.exp(a)*sp.cos(b)-1, sp.exp(a)*sp.sin(b))
    if name == 'log1p':
        if b == 0 and x.value.real > -1:
            return g.emit(value, sp.log(1+a))
        return g.emit(value, sp.log((1+a)**2+b*b)/2, sp.atan2(b, 1+a))
    if name == 'sin':
        return g.emit(value, sp.sin(a)*sp.cosh(b), sp.cos(a)*sp.sinh(b))
    if name == 'cos':
        return g.emit(value, sp.cos(a)*sp.cosh(b), -sp.sin(a)*sp.sinh(b))
    if b != 0:
        raise ValueError('this primitive is used only on real inputs')
    return g.emit(value, getattr(sp, name)(a))


class SymbolicArithmetic(Arithmetic):
    """Same Formulas methods; exact expressions and symbolic constants."""
    def __init__(self, graph):
        self.graph = graph
        self.pi = graph.number(mp.pi, graph.sp.pi, 0)
        self.e = graph.number(mp.e, graph.sp.E, 0)
        self.j = graph.number(1j, 0, 1)

    def mpf(self, x):
        if isinstance(x, SymbolicNumber):
            if x.im != 0:
                raise TypeError('mpf received a nonreal symbolic scalar')
            return x
        return self.graph.number(x, self.graph.sp.Rational(str(x)), 0)

    def mpc(self, x=0, y=None):
        return self.graph.number(x) if y is None else self.graph.number(x)+self.j*y

    def sqrt(self, x):
        return self.graph.number(x)**0.5

    def norm(self, values):
        return self.sqrt(sum(norm_squared(x) for x in values))

    def floor(self, x):
        value = mp.floor(x.value.real)
        x.graph.choices.append(('floor', x.re, int(value)))
        return int(value)

    def nint(self, x):
        value = mp.nint(x.value.real)
        x.graph.choices.append(('nint', x.re, int(value)))
        return int(value)

    def sign(self, x):
        value = mp.sign(x.value.real)
        x.graph.choices.append(('sign', x.re, int(value)))
        return int(value)

    def __getattribute__(self, name):
        if name in ('exp', 'log', 'expm1', 'log1p', 'sin', 'cos', 'sinh', 'cosh', 'asinh', 'arg'):
            return lambda x: symbolic_unary(name, self.graph.number(x))
        return super().__getattribute__(name)


@dataclass(frozen=True)
class SymbolicTensor:
    graph: ExpressionGraph
    outputs: tuple
    chart: str
    missing: int
    sign: int

    def differentiate(self):
        """Return exact A, six Hessians and their named assignments."""
        return self.graph.differentiate(self.outputs)

    def tensor_expressions(self):
        """J and dJ with named exact derivatives; inverse means adj(A)/det(A)."""
        sp = self.graph.sp
        A, H, assignments = self.differentiate()
        A = sp.ImmutableMatrix(A)
        J0 = sp.diag(*([sp.Matrix([[0, -1], [1, 0]])]*3))
        inverse = sp.Inverse(A)
        J = sp.MatMul(inverse, J0, A, evaluate=False)
        dJ = []
        for k in range(6):
            Ak = sp.ImmutableMatrix(6, 6, lambda a, i: H[a][i, k])
            dJ.append(sp.MatMul(inverse, sp.MatAdd(
                sp.MatMul(J0, Ak, evaluate=False),
                sp.MatMul(-1, Ak, J, evaluate=False), evaluate=False), evaluate=False))
        return J, tuple(dJ), assignments

    @property
    def matrix(self):
        """Exact global-frame matrix expression, with named scalar definitions."""
        sp = self.graph.sp
        x = list(self.graph.inputs[:4])
        x.insert(self.missing,self.sign*sp.sqrt(1-sum(v*v for v in x)))
        T,inverse = frame_matrices(x,self.missing,dtype=object)
        J,_,_ = self.tensor_expressions()
        return sp.MatMul(sp.ImmutableMatrix(T),J,sp.ImmutableMatrix(inverse),evaluate=False)


def symbolic_definition(coordinates, *, missing=4, sign=1, dps=60):
    """Trace the actual coordinate map with six free symbolic inputs.

    The supplied numerical witness chooses a formula branch. The result is
    a local symbolic definition; the contracts supply global coverage/joins.
    No floating-point residual is accepted as a symbolic proof.
    """
    if len(coordinates)!=6 or missing not in range(5) or sign not in (-1,1):
        raise ValueError('six coordinates and a valid hemisphere chart are required')
    with mp.workdps(dps):
        if sum(mp.mpf(v)**2 for v in coordinates[:4])>=1:
            raise ValueError('the witness must lie in the open hemisphere chart')
        graph = ExpressionGraph()
        ops = SymbolicArithmetic(graph)
        formulas = Formulas(ops)
        q = [graph.number(v, s, 0) for v, s in zip(coordinates, graph.inputs)]
        sphere = q[:4]
        sphere.insert(missing, sign*ops.sqrt(1-sum(v*v for v in sphere)))

        def phase(theta, phi):
            g = formulas.boundary_mismatch(theta, phi)
            value = phase_value(float(theta), float(phi), PHASE_STEPS)[0]
            exact_sum = graph.sp.Function('FinitePhase65536')(theta.re, phi.re)
            symbol = graph.real_node(exact_sum, graph.sp.atan2(g.im, g.re), 'finite_phase')
            return graph.number(value, symbol, 0)

        kind, values = formulas.coordinates(
            (sphere[0]+ops.j*sphere[1], sphere[2]+ops.j*sphere[3]), sphere[4], q[4:],
            dps=dps, phase_local=phase)
        outputs = tuple(part for z in values for part in (z.re, z.im))
        return SymbolicTensor(graph, outputs, kind, missing, sign)


@dataclass(frozen=True)
class ProofContract:
    name: str
    methods: tuple
    guarantee: str
    reference: str


def proof_contracts():
    """Explicit written proof obligations, not assertions proved by this list.

    References are to docs/s4t2_symbolic_programme.md and its existing C/L/N
    checks. Those checks verify identities and arithmetic; the written domain,
    inequality, all-order smoothness and winding arguments require review.
    The focused verifier additionally checks the actual formula correspondence.
    """
    return (
        ProofContract('flat step', ('scalar_flat_step',),
            'Smooth, constant outside [0,1], flat at both ends; values in [0,1].',
            'L1; C40/C45/C47; sections 11.2 and 11.9'),
        ProofContract('two-way height', (
            '_terminal_height_data','_height_G','_height_R','_height_terminal_forward',
            '_height_terminal_inverse_delta','_height_g_inverse_delta','inverse_height_moment',
            'elementary_forward_moment','elementary_cusp_logradius','elementary_cusp_height',
            'elementary_cusp_height_chart','elementary_cusp_reference'),
            'On the stated q/weight domains, positive denominators, increasing moment, '
            'inverse identities, smooth central axes and all-order joins/deck compatibility.',
            'section 11.2; C29, C43-C45, C49-C50; written terminal bounds and inner-seam argument'),
        ProofContract('pole inverse pair', (
            'pole_raw_inverse','_elementary_pole_radial','elementary_pole_raw_forward',
            'elementary_pole_patch','elementary_pole_forward','_elementary_pole_inclusion','pole_slice'),
            'The raw inverse is used off zero; radial conjugation extends smoothly as the identity '
            'on the inner ball, has an explicit inverse and agrees with the raw lift on the outer overlap.',
            'section 11.3 and cap domains in 11.6; C25-C28, C46; N10; all-order blow-up argument'),
        ProofContract('orbit coordinates', (
            'base_reference','_elementary_raw_orbit','_elementary_hopf_reference','scalar_lattice_coords'),
            'The dominant-coordinate formulas and their inverses cover the free orbit space; '
            'all argument choices represent the same quotient point and the orbit differential is invertible.',
            'sections 11.2/11.6/11.9; C38-C45, C49-C50; N11'),
        ProofContract('source line', ('bundle_embedding',),
            'A nonzero equivariant reference defines the same line through every chart/deck change; '
            'normalization has squared norm at least 1/2 before its positive common factors.',
            'section 11.1; C21-C22, C38-C45, C50'),
        ProofContract('torus restoration', (
            'scalar_finite_height','elementary_restore_torus','_elementary_real_coeffs',
            '_elementary_cusp_basis','_elementary_reduce_cusp','_elementary_source_coordinates'),
            'The finite height has the exact deck shift; restoration and its inverse respect '
            'integer torus representatives and holomorphic changes of the target charts.',
            'sections 10.2/11.6/11.9; C30k, C38, C41-C44, C50g'),
        ProofContract('regional lifts', ('middle_lift','cusp_lift'),
            'Every active projected row is nonzero; corrected overlaps have positive real part; '
            'the normalized blends are unit circle isomorphisms with adjoint inverses.',
            'sections 12.1-12.5; C56-C60; explicit positive-polynomial and hemisphere bounds'),
        ProofContract('boundary extension', ('collar_offset','extended_cusp_lift'),
            'Source slot weights cancel exactly; the two-row formula extends through the boundary '
            'and respects angular/height deck changes, with the same all-order endpoint jets.',
            'sections 12.6-12.7; C61-C62; C45 flat height join'),
        ProofContract('periodic finite phase', ('boundary_mismatch',),
            'The mismatch has unit modulus, both winding numbers zero, and derivative bounds '
            '1900/2700. The 65536-term real lift is smooth and periodic with positive increment real parts.',
            'sections 12.8-12.9; C63-C64; exact endpoint derivative identity C66'),
        ProofContract('global assembly', ('circle_section','coordinates','forward'),
            'Flat radial retreats match the capped cusp and middle isomorphisms to all orders. '
            'The explicit inverse reverses the orbit, circle and torus maps on a covering of the manifold.',
            'sections 12.5/12.9/12.11; equations (12.20), (12.35); inherited cap/orbit inverse domains'),
    )


class NumericalFormulas(Formulas):
    """Numerical tail evaluation, separate from the exact Formulas definition.

    These are the existing cusp.py precision shortcuts. They avoid forming
    double exponentials in flat tails. The value cutoffs are described in
    programme 11.7 and the cusp.py source. They are NOT identities, are never
    used by the symbolic definition, and give no uniform bound here for
    the numerical error in J or dJ.
    All remaining calculations execute the exact formula bodies above.
    """
    def _height_R(self, value, inverse=False):
        c = mp.mpf('.2')
        if 0 < value < c and c/(c-value) > (mp.mp.dps+10)*mp.log(10):
            return value
        return super()._height_R(value,inverse)

    def _elementary_pole_radial(self, pair, inverse=False):
        radius = sum(norm_squared(z) for z in pair)
        c = mp.mpf('1.5e-8')
        if 0 < radius < c and c/(c-radius) > (mp.mp.dps+10)*mp.log(10):
            return tuple(pair)
        return super()._elementary_pole_radial(pair,inverse)

    def _height_terminal_forward(self, delta_z, data):
        _,_,_,_,_,z0,zcut,eta,E0,Acut,Lambda,mu = data
        if 0 < delta_z < zcut-z0:
            Ez = self.ops.exp(eta/(zcut-z0-delta_z))
            C0,C1 = 1+mu-2*Lambda,Lambda-mu
            K = 1+abs(C0)*self.ops.exp(E0)+abs(C1)*self.ops.exp(E0+Acut)
            if Ez > (mp.mp.dps+20)*mp.log(10)+self.ops.log(K):
                return z0+delta_z
        return super()._height_terminal_forward(delta_z,data)

    def _height_terminal_inverse_delta(self, y, data):
        _,_,_,_,_,z0,zcut,eta,E0,Acut,Lambda,mu = data
        if 0 < y < zcut:
            Ey = self.ops.exp(eta/(zcut-y))
            C0,C1 = 1+mu-2*Lambda,Lambda-mu
            K = 1+abs(C0)+2*Lambda*abs(C1)
            if Ey > E0+(mp.mp.dps+20)*mp.log(10)+self.ops.log(K):
                return y-z0
        return super()._height_terminal_inverse_delta(y,data)

    def _inner_tail(self, pair, dps):
        radius = sum(norm_squared(z) for z in pair)
        inner,outer = mp.mpf('1e-8'),mp.mpf('1.5e-8')
        if inner < radius < (inner+outer)/2:
            x = self.ops.exp(3)*self.ops.expm1(3*(radius-inner)/(outer-radius))
            return x-1/x <= -(dps+10)*mp.log(10)
        return False

    def elementary_pole_patch(self, z, infinity=False, dps=50):
        if self._inner_tail(z,dps):
            return tuple(z)
        return super().elementary_pole_patch(z,infinity,dps)

    def elementary_pole_forward(self, w, infinity=False, dps=50):
        if self._inner_tail(w,dps):
            return tuple(w)
        return super().elementary_pole_forward(w,infinity,dps)

