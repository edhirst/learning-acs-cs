"""Numerical reconstruction, on its fibration model, of the complex structure
given by the supplied S^4 x T^2 period model (the explicit period model of the
paper), in the style of extra/s6_reconstruction.py.

Model over C*_t. The fibres are the 2-tori X_t = C^2/Lambda_t (additive
coordinates (z, w), zeta = e^{iz}) with
    e1 = (2 pi, 0),
    e2 = (-i alpha_2(t), eb_2(t)),   e3 = (-i alpha_3(t), eb_3(t)),
    e4 = (beta - i log t, eb_4(t)).
In the real coordinates (Re t, Im t, x_1..x_4), x the coordinates in the basis
e_j, the complex structure J has (1,0)-forms dt and Pi(t) dx, Pi = [e1 e2 e3 e4]
(2 x 4). J is an ACS iff the e_j are R-independent (det_R > 0) and integrable
iff Pi is holomorphic in t (it is, by construction: N = 0 up to round-off).

Period data. Holomorphy of the chart change at t = inf,
Theta_lam(zeta, w, t) = (zeta e^{lam w/t}, w/t^2, 1/t), together with the
Mumford-Tate hypothesis (H1) at both cusps, forces the shape
    eb_j = a_j + c_j t^2,   alpha_j = -lam c_j t   (j = 2, 3),   eb_4 = b_4 t,
with lam > 0 and a_j, c_j, b_4, beta in C; the supplied model is lam = 1,
(a_2, a_3) = (1, i), (c_2, c_3) = (-1, i), b_4 = 2, beta = 0 (and its
deformations in (lam, beta)). The remaining requirements are
    (H2) omega(eb_2, eb_3) != 0 near both cusps, with OPPOSITE chirality:
         omega(a_2, a_3) > 0 > omega(c_2, c_3);
    (ND) det_R(e1, e2, e3, e4)(t) > 0 for every t in C*.
`learn_period_data` finds such data from a random start by Levenberg-Marquardt
on a barrier for (ND) and (H2) at resampled points of C* (a feasibility
problem: no target values, no formula from the supplied model); `certify` checks the
resulting structure on the model.
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import math

import numpy as np
import torch

from extra.s6_reconstruction import J_of_Pi, nijenhuis

DTYPE = torch.float64
C128 = torch.complex128
PAPER = {'lam': 1.0, 'a2': 1.0 + 0j, 'a3': 1j, 'c2': -1.0 + 0j, 'c3': 1j, 'b4': 2.0 + 0j, 'beta': 0j}
KEYS = ['a2', 'a3', 'c2', 'c3', 'b4', 'beta']


def omega(u, v):
    """omega(u, v) = Im(conj(u) v)."""
    return (u.conj() * v).imag


class PeriodData:
    """Holds (lam, a2, a3, c2, c3, b4, beta); `theta` is the real parameter vector."""

    def __init__(self, lam=1.0, **kw):
        self.lam = float(lam)
        for k in KEYS:
            setattr(self, k, complex(kw.get(k, 0j)))

    @classmethod
    def paper(cls):
        return cls(**PAPER)

    @classmethod
    def from_theta(cls, th, lam=1.0):
        th = [float(v) for v in th]
        kw = {k: complex(th[2 * i], th[2 * i + 1]) for i, k in enumerate(KEYS)}
        return cls(lam=lam, **kw)

    def theta(self):
        return np.array([p for k in KEYS for p in (getattr(self, k).real, getattr(self, k).imag)])

    def as_dict(self):
        d = {'lam': self.lam}
        d.update({k: [getattr(self, k).real, getattr(self, k).imag] for k in KEYS})
        return d

    # ---- periods -----------------------------------------------------------------------
    def Pi(self, t, k=0):
        """(B, 2, 4) period matrix at complex t (torch complex128); k selects the
        branch log t + 2 pi i k."""
        t = t.to(C128)
        lam = self.lam
        a2, a3, c2, c3, b4, beta = (torch.tensor(getattr(self, n), dtype=C128) for n in KEYS)
        lg = torch.log(t) + 2j * math.pi * k
        z = torch.zeros_like(t)
        r1 = torch.stack([z + 2 * math.pi, 1j * lam * c2 * t, 1j * lam * c3 * t, beta - 1j * lg], -1)
        r2 = torch.stack([z, a2 + c2 * t ** 2, a3 + c3 * t ** 2, b4 * t], -1)
        return torch.stack([r1, r2], -2)

    def Pi_inf(self, s, k=0):
        """Period matrix of the chart at infinity (coordinates (z'', w', s), basis
        (e1, e2, e3, -e4)), from the formulas: alpha'_j = lam a_j s, eb'_j = c_j + a_j s^2,
        -e4 -> (-i log(nu s), -b4 s), nu = exp(-i beta - lam b4)."""
        s = s.to(C128)
        lam = self.lam
        a2, a3, c2, c3, b4, beta = (torch.tensor(getattr(self, n), dtype=C128) for n in KEYS)
        lg = torch.log(s) - 2j * math.pi * k
        log_nu = -1j * beta - lam * b4
        z = torch.zeros_like(s)
        r1 = torch.stack([z + 2 * math.pi, -1j * lam * a2 * s, -1j * lam * a3 * s, -1j * (log_nu + lg)], -1)
        r2 = torch.stack([z, c2 + a2 * s ** 2, c3 + a3 * s ** 2, -b4 * s], -1)
        return torch.stack([r1, r2], -2)

    def det(self, t):
        """det_R(e1..e4) in the complex orientation (Re z, Im z, Re w, Im w)."""
        Pi = self.Pi(t)
        R = torch.stack([Pi[:, 0].real, Pi[:, 0].imag, Pi[:, 1].real, Pi[:, 1].imag], 1)   # (B, 4, 4) rows
        return torch.linalg.det(R)

    def det_normalised(self, t):
        """det / (2 pi (1 + r^4)(1 + |log r|)): O(1) uniformly on C* (the asymptotic
        forms are 2 pi (-log r) omega(a2, a3) at 0 and 2 pi (-log r) r^4 omega(c2, c3) at inf)."""
        r = t.abs()
        return self.det(t) / (2 * math.pi * (1 + r ** 4) * (1 + torch.log(r).abs()))

    def chirality(self):
        a2, a3, c2, c3 = (torch.tensor(getattr(self, n), dtype=C128) for n in ('a2', 'a3', 'c2', 'c3'))
        return omega(a2, a3).item(), omega(c2, c3).item()


# ---------------------------------------------------------------------------
# J on the model and its certificate
# ---------------------------------------------------------------------------

def J_field(pd, ab, k=0):
    t = torch.complex(ab[:, 0], ab[:, 1])
    return J_of_Pi(pd.Pi(t, k))


def sample_annulus(n, rmin, rmax, rng, avoid_cut=0.02):
    u = rng.uniform(np.log(rmin), np.log(rmax), n)
    th = rng.uniform(-np.pi + avoid_cut, np.pi - avoid_cut, n)
    return np.exp(u) * np.exp(1j * th)


def min_det_adversarial(pd, n=4000, steps=60, L=7.0, seed=0):
    """Minimum of the normalised determinant over C* (|t| in [e^-L, e^L]):
    dense sample, then descent from the 64 smallest values in (log r, arg t)."""
    rng = np.random.default_rng(seed)
    t = torch.tensor(sample_annulus(n, np.exp(-L), np.exp(L), rng, 0.0), dtype=C128)
    q = pd.det_normalised(t)
    idx = torch.topk(-q, 64).indices
    p = torch.stack([torch.log(t[idx].abs()), torch.angle(t[idx])], 1)
    step = torch.full((64, 1), 0.05, dtype=DTYPE)
    f = lambda p: pd.det_normalised(torch.exp(p[:, 0]) * torch.exp(1j * p[:, 1].to(C128)))
    cur = f(p).detach()
    for _ in range(steps):
        pg = p.clone().requires_grad_(True)
        g, = torch.autograd.grad(f(pg).sum(), pg)
        pn = p - step * g / (g.norm(dim=1, keepdim=True) + 1e-30)
        pn[:, 0] = pn[:, 0].clamp(-L, L)
        vn = f(pn).detach()
        down = vn < cur
        p = torch.where(down.unsqueeze(1), pn.detach(), p)
        cur = torch.where(down, vn, cur)
        step = torch.where(down.unsqueeze(1), step * 1.2, step * 0.5)
    return min(q.min().item(), cur.min().item())


def height_ratio(pd, t, at_inf=False):
    """P(t)/log|t| at the cusp (the filling condition): P = -det_R/(2 pi omega(eb2, eb3))
    for the lattice of the filling at that cusp."""
    Pi = pd.Pi_inf(t) if at_inf else pd.Pi(t)
    R = torch.stack([Pi[:, 0].real, Pi[:, 0].imag, Pi[:, 1].real, Pi[:, 1].imag], 1)
    det = torch.linalg.det(R)
    om = omega(Pi[:, 1, 1], Pi[:, 1, 2])
    P = -det / (2 * math.pi * om)
    return P / torch.log(t.abs())


def filling_radius(pd, at_inf=False, n_ang=64):
    """Largest eps (on a log grid in (1e-6, 1/16]) such that omega(eb2, eb3) has one sign
    on 0 < |t| <= eps and P/log|t| lies in [1/2, 3/2] there (Step 0 of the filling)."""
    radii = np.logspace(-6, np.log10(1 / 16), 120)
    th = np.linspace(-np.pi + 1e-3, np.pi - 1e-3, n_ang)
    ok_eps = 0.0
    for r in radii:
        t = torch.tensor((r * np.exp(1j * th)), dtype=C128)
        Pi = pd.Pi_inf(t) if at_inf else pd.Pi(t)
        om = omega(Pi[:, 1, 1], Pi[:, 1, 2])
        hr = height_ratio(pd, t, at_inf)
        one_sign = bool((om > 0).all()) or bool((om < 0).all())
        if one_sign and bool(((hr >= 0.5) & (hr <= 1.5)).all()):
            ok_eps = float(r)
        else:
            break
    return ok_eps


def log_transform_xi(pd, ab, k=0):
    """Real fibre translation xi(t) with Pi(t) xi = sigma(t), sigma = (pi s(s-1) + beta s,
    b4 t s), s = log t/(2 pi i): the multiplicity-one logarithmic transform."""
    t = torch.complex(ab[:, 0], ab[:, 1])
    s = (torch.log(t) + 2j * math.pi * k) / (2j * math.pi)
    beta = torch.tensor(pd.beta, dtype=C128); b4 = torch.tensor(pd.b4, dtype=C128)
    sig = torch.stack([math.pi * s * (s - 1) + beta * s, b4 * t * s], -1)
    Pi = pd.Pi(t, k)
    R = torch.stack([Pi[:, 0].real, Pi[:, 0].imag, Pi[:, 1].real, Pi[:, 1].imag], 1)
    rhs = torch.stack([sig[:, 0].real, sig[:, 0].imag, sig[:, 1].real, sig[:, 1].imag], 1)
    return (torch.linalg.inv(R) @ rhs.unsqueeze(-1)).squeeze(-1)


def certify(pd, eps_frac=0.5, seed=20261004, eps=None):
    """Certificate of the period data `pd`. The regular-part tests run on the annulus
    eps*eps_frac <= |t| <= 1/(eps*eps_frac) and the disc tests inside it, with eps the
    data's own filling radius, or `eps` if given: pass a common eps to compare two sets
    of period data on the same annulus (their own radii differ, and the transition
    discrepancies grow with the tested range)."""
    rng = np.random.default_rng(seed)
    rep = {'data': pd.as_dict()}
    rep['chirality_0'], rep['chirality_inf'] = pd.chirality()
    rep['det_normalised_min'] = min_det_adversarial(pd)
    eps0, epsi = filling_radius(pd, False), filling_radius(pd, True)
    rep['filling_eps_0'], rep['filling_eps_inf'] = eps0, epsi
    eps = max(1e-6, min(eps0, epsi)) if eps is None else float(eps)
    rep['annulus_eps'] = eps
    # complex structure on the regular part (annulus between the two filling discs)
    tt = sample_annulus(3000, eps * eps_frac, 1 / (eps * eps_frac), rng)
    ab = torch.tensor(np.stack([tt.real, tt.imag], 1), dtype=DTYPE)
    Jf = lambda p: J_field(pd, p)
    J = Jf(ab)
    rep['acs_residual_max'] = (J @ J + torch.eye(6, dtype=DTYPE)).abs().max().item()
    N = torch.cat([nijenhuis(Jf, ab[i:i + 500]).detach().flatten(1).norm(dim=1) for i in range(0, len(ab), 500)])
    rep['N_mean'], rep['N_max'] = N.mean().item(), N.max().item()
    # monodromy across the branch cut of log t: e4 -> e4 + e1, x -> T x
    T = torch.eye(4, dtype=DTYPE); T[0, 3] = 1.0
    D = torch.block_diag(torch.eye(2, dtype=DTYPE), torch.linalg.inv(T))
    J0, J1 = J_field(pd, ab, 0), J_field(pd, ab, 1)
    rep['monodromy_consistency_max'] = (J1 - D @ J0 @ torch.linalg.inv(D)).flatten(1).norm(dim=1).max().item()
    Dw = torch.block_diag(torch.eye(2, dtype=DTYPE), T)                      # wrong transport must fail
    rep['monodromy_wrong_transport'] = (J1 - Dw @ J0 @ torch.linalg.inv(Dw)).flatten(1).norm(dim=1).max().item()
    # chart at infinity: Pi_inf(s) must equal L(t) Pi(t) diag(1,1,1,-1), L = [[1, -i lam/t], [0, 1/t^2]]
    t = torch.complex(ab[:, 0], ab[:, 1])
    Lt = torch.zeros(len(t), 2, 2, dtype=C128)
    Lt[:, 0, 0] = 1; Lt[:, 0, 1] = -1j * pd.lam / t; Lt[:, 1, 1] = 1 / t ** 2
    flip = torch.diag(torch.tensor([1, 1, 1, -1], dtype=C128))
    lhs = pd.Pi_inf(1 / t)
    rhs = Lt @ pd.Pi(t) @ flip
    # equal up to adding a multiple of e1 to the last column (branch of log s vs log t)
    diff = lhs - rhs
    k = torch.round((diff[:, 0, 3].real) / (2 * math.pi))
    diff[:, 0, 3] = diff[:, 0, 3] - 2 * math.pi * k
    rep['chart_infinity_consistency_max'] = diff.abs().max().item()
    # logarithmic transform: x -> x + xi(t) preserves J iff sigma is holomorphic
    def phi_J(p):
        xi_fn = lambda q: log_transform_xi(pd, q.unsqueeze(0)).squeeze(0)
        Dxi = torch.func.vmap(torch.func.jacfwd(xi_fn))(p)                    # (B, 4, 2)
        DPhi = torch.eye(6, dtype=DTYPE).repeat(len(p), 1, 1)
        DPhi[:, 2:, :2] = Dxi
        return DPhi @ J_field(pd, p) @ torch.linalg.inv(DPhi)
    near0 = sample_annulus(800, eps * eps_frac * 0.5, eps * eps_frac, rng)
    ab0 = torch.tensor(np.stack([near0.real, near0.imag], 1), dtype=DTYPE)
    rep['log_transform_consistency_max'] = (phi_J(ab0) - J_field(pd, ab0)).flatten(1).norm(dim=1).max().item()
    # N on the two filling discs (regular points; the central fibres are the toric models)
    # in the holomorphic coordinate u = log t (resp. log s), so that the size of the
    # coordinate tensor does not blow up like 1/|t| as the central fibre is approached
    for nm, inf in (('0', False), ('inf', True)):
        uu = rng.uniform(np.log(1e-8), np.log(eps * eps_frac), 800) + 1j * rng.uniform(-np.pi + 0.02, np.pi - 0.02, 800)
        abd = torch.tensor(np.stack([uu.real, uu.imag], 1), dtype=DTYPE)
        def Jd(p, inf=inf):
            tt_ = torch.exp(torch.complex(p[:, 0], p[:, 1]))
            return J_of_Pi(pd.Pi_inf(tt_) if inf else pd.Pi(tt_))
        Nd = nijenhuis(Jd, abd).detach().flatten(1).norm(dim=1)
        rep[f'N_disc_{nm}_max'] = Nd.max().item()
    keys = ['N_max', 'monodromy_consistency_max', 'chart_infinity_consistency_max',
            'log_transform_consistency_max', 'N_disc_0_max', 'N_disc_inf_max']
    rep['certified_sup'] = max(rep[k] for k in keys)
    rep['admissible'] = bool(rep['det_normalised_min'] > 0 and rep['chirality_0'] > 0 > rep['chirality_inf']
                             and eps0 > 0 and epsi > 0)
    return rep


# ---------------------------------------------------------------------------
# Physics-informed recovery of admissible period data
# ---------------------------------------------------------------------------

def learn_period_data(seed=0, delta=0.05, rounds=40, n_pts=512, L=7.0, lam=1.0, log=print,
                      maximise=True, grow=1.15, max_steps=60):
    """Feasibility by Levenberg-Marquardt: residuals relu(delta - q(t)) on resampled
    t in C* (q = normalised det_R), relu(delta -/+ chirality) and a scale gauge
    |a|^2 + |c|^2 = 4. Random start; no values from the supplied model. With `maximise`, the
    margin delta is then raised geometrically while the problem stays feasible
    (the robust, well-conditioned admissible data), and the last feasible data
    are returned."""
    from src.optim.lm import lm_minimise
    g = np.random.default_rng(seed)
    th = torch.tensor(g.normal(0, 1, 2 * len(KEYS)), dtype=DTYPE)

    def resid(theta, t, dl):
        pd = _TorchPD(theta, lam)
        q = pd.det_normalised(t)
        w0, wi = pd.chirality_t()
        gauge = (theta[:8] ** 2).sum() - 4.0
        return torch.cat([torch.relu(dl - q) / math.sqrt(len(t)),
                          torch.stack([torch.relu(dl - w0), torch.relu(dl + wi), 0.1 * gauge])])

    def solve(th, dl, n_rounds, lam_lm):
        for rnd in range(n_rounds):
            t = torch.tensor(sample_annulus(n_pts, np.exp(-L), np.exp(L), g, 0.0), dtype=C128)
            th, lg = lm_minimise(lambda p, c: resid(p, c, dl), th, [t], n_iter=10, lam0=lam_lm)
            lam_lm = min(max(lg[-1]['lam'], 1e-12), 1e6)
            if lg[-1]['loss'] < 1e-26:
                pd = PeriodData.from_theta(th.numpy(), lam)
                md = min_det_adversarial(pd, n=2000, steps=30, L=L, seed=rnd)
                c0, ci = pd.chirality()
                if md > 0.9 * dl and c0 >= 0.999 * dl and ci <= -0.999 * dl:
                    return th, True, rnd + 1, md, lam_lm
        return th, False, n_rounds, None, lam_lm

    hist = []
    th, ok, nr, md, lam_lm = solve(th, delta, rounds, 1e-2)
    hist.append({'delta': delta, 'feasible': ok, 'rounds': nr, 'min_det': md})
    if log:
        log(f"  seed {seed}: delta {delta:.3f} feasible={ok} after {nr} rounds (min det {md})")
    if not ok:
        return PeriodData.from_theta(th.numpy(), lam), hist
    best = (delta, th.clone())
    dl = delta
    for _ in range(max_steps if maximise else 0):
        dl = dl * grow
        th2, ok, nr, md, lam_lm = solve(best[1].clone(), dl, 12, 1e-2)
        hist.append({'delta': dl, 'feasible': ok, 'rounds': nr, 'min_det': md})
        if not ok:
            break
        best = (dl, th2.clone())
        if log:
            log(f"  seed {seed}: margin {dl:.3f} feasible (min det {md:.3f})")
    return PeriodData.from_theta(best[1].numpy(), lam), hist


class _TorchPD(PeriodData):
    """Differentiable view of PeriodData on a torch parameter vector."""

    def __init__(self, theta, lam):
        self.lam = lam
        self._th = theta

    def _c(self, i):
        return torch.complex(self._th[2 * i], self._th[2 * i + 1])

    def Pi(self, t, k=0):
        a2, a3, c2, c3, b4, beta = (self._c(i) for i in range(6))
        lg = torch.log(t) + 2j * math.pi * k
        z = torch.zeros_like(t)
        r1 = torch.stack([z + 2 * math.pi, 1j * self.lam * c2 * t, 1j * self.lam * c3 * t, beta - 1j * lg], -1)
        r2 = torch.stack([z, a2 + c2 * t ** 2, a3 + c3 * t ** 2, b4 * t], -1)
        return torch.stack([r1, r2], -2)

    def chirality_t(self):
        return omega(self._c(0), self._c(1)), omega(self._c(2), self._c(3))
