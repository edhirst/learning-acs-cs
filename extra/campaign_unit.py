"""Uniform chart-free PINN campaign — one (case, task, seed) unit.

Formulation is chosen by topology (src/geometry/explicit_acs.py):
  parallelizable  (s1xs1, s3xs1, s3xs3, s5xs1)  → global-frame Nijenhuis, no charts
  non-parallelizable (s2, s4, s6, s2xs2, s4xs2) → ambient-projector Nijenhuis, no charts

Tasks (unsupervised PINN, identical across formulations):
  acs  free endomorphism field K (frame: K = net(x); ambient: K = P·net(x)·P),
       loss = mean ‖K² + Id_T‖²_F. From scratch, no base point. Certified a
       posteriori by the polar correction J* = K(−K²_T)^{−1/2} when
       sup ‖K²+Id_T‖ < 1.
  cs   exact-ACS parametrization J = e^{B_s}⋯e^{B_1} J_base e^{−B_1}⋯e^{−B_s}
       (B_k tangential net outputs, initialised at O(1) ⇒ a generic
       NON-integrable start), loss = mean ‖N(J)‖² only. J_base is the
       canonical J₀ (frame cases) or a classical explicit ACS (ambient cases);
       s4 has no ACS, so its cs unit runs the free-field two-term loss
       (w·‖K²+P‖² + ‖N‖²) as the "nothing exists" control.

Every unit is restartable (state.pt), time-guarded (--wall-minutes), logs
progress.jsonl, and ends with a dense evaluation + the ORIGINAL chart-protocol
scores (frame fields are bridged to ambient fields, then to ball charts), in
result.json.

Usage:
  python extra/campaign_unit.py --case s3xs3 --task cs --seed 7 --out DIR \
      [--steps 40000] [--hidden 256,256,256,256] [--stages 2] [--rar] \
      [--wall-minutes 4200] [--eval-every 2000] [--ckpt-every 1000]
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import os
import sys
import json
import time
import argparse
import torch
import torch.nn as nn

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from src.geometry.frames import FRAMES
from src.geometry.explicit_acs import (ALL_CASES, factor_dims, is_frame_case,
                                       ambient_manifold, base_acs)
from src.loss.frame_loss import frame_nij_norm
from src.loss.ambient_loss import (ambient_nij_norm, compute_J_and_DJ,
                                   compute_P_and_DP)
from src.geometry.sphere import EvenSphere, RADIAL_MIDPOINT
from src.geometry.sphere_product import SphereProduct
from src.loss.sphere_loss import _nij_tensor, OverlapLoss
from src.loss.product_loss import ProductOverlapLoss

DTYPE = torch.float64
BATCH = 192
OVERLAP_SIZE = 0.2
R_CUT = RADIAL_MIDPOINT + OVERLAP_SIZE

# Results table of the paper (before this work), for side-by-side reporting.
TABLE_BEFORE = {
    's2':    {'acs': 1.63e-4, 'cs': 5.27e-4},
    's4':    {'acs': 2.83,    'cs': 8.50},
    's6':    {'acs': 3.57e-2, 'cs': 9.87},
    's1xs1': {'acs': 5.81e-5, 'cs': 1.90e-4},
    's2xs2': {'acs': 9.27e-3, 'cs': 4.04e-2},
    's3xs1': {'acs': 1.11e-2, 'cs': 2.00e-1},
    's3xs3': {'acs': 1.61e-1, 'cs': 7.63},
    's4xs2': {'acs': 1.91,    'cs': 6.28},
    's5xs1': {'acs': 7.98e-2, 'cs': 6.42},
}


# ---------------------------------------------------------------------------
# Geometry wrapper
# ---------------------------------------------------------------------------

class Geometry:
    def __init__(self, case):
        self.case = case
        self.dims = factor_dims(case)
        self.d = sum(self.dims)
        if is_frame_case(case):
            self.kind = 'frame'
            self.frame = FRAMES[case]()
            self.N = self.frame.ambient_dim
            probe = self.frame.sample(200)
            c = self.frame.structure_functions(probe)
            self.c_const = c[:1] if (c - c.mean(0, keepdim=True)).abs().max() < 1e-10 else None
        else:
            self.kind = 'ambient'
            self.manifold = ambient_manifold(case)
            self.N = self.manifold.ambient_dim
            self.c_const = None

    def sample(self, n):
        return self.frame.sample(n) if self.kind == 'frame' else self.manifold.sample(n)

    def structure(self, x):
        if self.c_const is not None:
            return self.c_const.expand(len(x), self.d, self.d, self.d)
        return self.frame.structure_functions(x)


def canonical_J0(d):
    J0 = torch.zeros(d, d, dtype=DTYPE)
    for k in range(d // 2):
        J0[2 * k, 2 * k + 1] = -1.0
        J0[2 * k + 1, 2 * k] = 1.0
    return J0


# ---------------------------------------------------------------------------
# Network + field construction
# ---------------------------------------------------------------------------

class Net(nn.Module):
    def __init__(self, n_in, n_mat, stages, hidden, out_scale):
        super().__init__()
        layers, prev = [], n_in
        for h in hidden:
            layers += [nn.Linear(prev, h), nn.GELU()]
            prev = h
        layers += [nn.Linear(prev, stages * n_mat * n_mat)]
        self.net = nn.Sequential(*layers)
        with torch.no_grad():
            self.net[-1].weight.mul_(out_scale)
            self.net[-1].bias.normal_(0.0, out_scale)
        self.n_mat, self.stages = n_mat, stages

    def forward(self, x):
        return self.net(x).reshape(-1, self.stages, self.n_mat, self.n_mat)


class PolyNet(nn.Module):
    """Degree-≤2 polynomial field in the ambient coordinates: features
    [1, x_i, x_i x_j (i≤j)] → n_mat×n_mat. The lowest 'harmonics' of the
    spectral curriculum: a low-degree field cannot form many orientation
    domain walls, so the ACS flow settles into one global orientation."""

    def __init__(self, n_in, n_mat, out_scale, feat='all2', blocks=None):
        """feat: harmonic content of the class —
             'all2'  : 1, x, x⊗x (degree ≤ 2)        'lin'   : x only (ℓ = 1)
             'lin+c' : 1, x                          'quad_h': x⊗x only
             'bilin' : x_A ⊗ x_B across the two product factors (bidegree (1,1);
                       `blocks` = ambient block sizes of the factors)
             'lin+bilin': x, x_A ⊗ x_B"""
        super().__init__()
        self.feat, self.n_in = feat, n_in
        iu = torch.triu_indices(n_in, n_in)
        self.register_buffer('iu', iu)
        if blocks is not None and len(blocks) == 2:
            na = blocks[0]
            ia, ib = torch.meshgrid(torch.arange(na), torch.arange(na, n_in), indexing='ij')
            self.register_buffer('bil', torch.stack([ia.reshape(-1), ib.reshape(-1)]))
        else:
            self.bil = None
        n_feat = self.features(torch.zeros(1, n_in)).shape[1]
        self.lin = nn.Linear(n_feat, n_mat * n_mat, bias=False)
        with torch.no_grad():
            self.lin.weight.normal_(0.0, out_scale / n_feat ** 0.5)
        self.n_mat, self.stages = n_mat, 1

    def features(self, x):
        one = torch.ones_like(x[:, :1])
        quad = x[:, self.iu[0]] * x[:, self.iu[1]]
        if self.feat == 'all2':
            return torch.cat([one, x, quad], dim=1)
        if self.feat == 'c':
            return one
        if self.feat == 'lin':
            return x
        if self.feat == 'lin+c':
            return torch.cat([one, x], dim=1)
        if self.feat == 'quad_h':
            return quad
        bil = x[:, self.bil[0]] * x[:, self.bil[1]]
        if self.feat == 'bilin':
            return bil
        if self.feat == 'lin+bilin':
            return torch.cat([x, bil], dim=1)
        raise ValueError(self.feat)

    def forward(self, x):
        return self.lin(self.features(x)).reshape(-1, 1, self.n_mat, self.n_mat)


class CurriculumNet(nn.Module):
    """Spectral curriculum field: K̃ = poly(x) + mlp(x). Phase 1 trains the
    polynomial alone (mlp frozen, initialised tiny) — climbing the HARMONIC
    LADDER of feature classes (see LADDER) until the field certifies; phase 2
    trains both."""

    def __init__(self, n_in, n_mat, hidden, out_scale_poly=0.3, out_scale_mlp=0.02,
                 feat='lin', blocks=None):
        super().__init__()
        self.poly = PolyNet(n_in, n_mat, out_scale_poly, feat=feat, blocks=blocks)
        self.mlp = Net(n_in, n_mat, 1, hidden, out_scale_mlp)
        self.n_mat, self.stages = n_mat, 1

    def forward(self, x):
        return self.poly(x) + self.mlp(x)


# Harmonic ladder: lowest spherical-harmonic content first. A low class has no
# room for defect cores, so the residual flow settles on a global structure
# (or, on S⁴, stays loudly non-zero). Product-only classes are skipped on
# single spheres. Acceptance = adversarial sup ‖K²+Id_T‖ < 1 (certifiable).
LADDER = ['c', 'lin', 'bilin', 'lin+c', 'lin+bilin', 'all2']


def make_field(geo, task, model, base_fn=None):
    """Returns (J_fn, meta). J_fn maps ambient points → matrix field:
    frame cases: (B, d, d) in frame components; ambient cases: (B, N, N).
    `base_fn` (CS task): the base ACS field, e.g. the PINN's own certified
    ACS from the acs unit; default = canonical J₀ / classical explicit ACS."""
    case = geo.case
    if geo.kind == 'frame':
        d = geo.d
        if task == 'acs':
            return (lambda x: model(x)[:, 0]), {'exact_acs': False}
        J0 = canonical_J0(d)

        def J_fn(x):
            raw = model(x)
            J = base_fn(x) if base_fn is not None else J0.unsqueeze(0).expand(len(x), d, d)
            for s in range(model.stages):
                Bs = raw[:, s]
                J = torch.matrix_exp(Bs) @ J @ torch.matrix_exp(-Bs)
            return J
        return J_fn, {'exact_acs': True}

    proj = geo.manifold.projector
    if task == 'acs' or (case == 's4' and base_fn is None):
        def K_fn(x):
            P = proj(x)
            return P @ model(x)[:, 0] @ P
        return K_fn, {'exact_acs': False}

    def J_fn(x):
        P = proj(x)
        raw = model(x)
        J = base_fn(x) if base_fn is not None else base_acs(case, x, proj)
        for s in range(model.stages):
            Bs = P @ raw[:, s] @ P
            J = torch.matrix_exp(Bs) @ J @ torch.matrix_exp(-Bs)
        return J
    return J_fn, {'exact_acs': True}


def load_acs_field(geo, run_dir):
    """Rebuild the frozen free ACS field K(x) of a finished acs unit."""
    ck = torch.load(os.path.join(run_dir, 'final.pt'), map_location='cpu', weights_only=False)
    n_mat = geo.d if geo.kind == 'frame' else geo.N
    if ck.get('model_type') == 'curriculum':
        model = CurriculumNet(geo.N, n_mat, ck['hidden'], feat=ck.get('feat', 'all2'),
                              blocks=[m + 1 for m in geo.dims]).to(DTYPE)
    else:
        model = Net(geo.N, n_mat, ck['stages'], ck['hidden'], 0.3).to(DTYPE)
    model.load_state_dict(ck['model'])
    model.eval()
    for p in model.parameters():
        p.requires_grad_(False)
    K_fn, _ = make_field(geo, 'acs', model)
    return K_fn


# ---------------------------------------------------------------------------
# Residuals
# ---------------------------------------------------------------------------

def identity_T(geo, x):
    if geo.kind == 'frame':
        return torch.eye(geo.d, dtype=DTYPE).expand(len(x), geo.d, geo.d)
    return geo.manifold.projector(x)


def acs_residual(geo, J_fn, x):
    K = J_fn(x)
    return torch.sqrt(((K @ K + identity_T(geo, x)) ** 2).sum(dim=(1, 2)) + 1e-20)


def nij_norm(geo, J_fn, x):
    """Per-sample intrinsic Nijenhuis norm of the field (frame or ambient)."""
    if geo.kind == 'frame':
        def single(xs):
            return J_fn(xs.unsqueeze(0)).squeeze(0)
        J = J_fn(x)
        dJ_amb = torch.func.vmap(torch.func.jacfwd(single))(x)            # (B,d,d,N)
        E = geo.frame.frame_fields(x)
        dJ_frame = torch.einsum('beaN,bcN->beac', dJ_amb, E)
        return frame_nij_norm(J, dJ_frame, geo.structure(x))
    J, DJ = compute_J_and_DJ(J_fn, x)
    P, DP = compute_P_and_DP(geo.manifold, x)
    return ambient_nij_norm(J, DJ, P, DP)


def dirichlet_energy(geo, K_fn, x):
    """INTRINSIC Dirichlet energy of a tangential endomorphism field:
    ambient: Σ_a ‖P (D_{Pe_a} K) P‖²_F (Levi-Civita covariant derivative on a
    submanifold = projected ambient derivative); frame: Σ_a ‖E_a(K)‖²_F.
    Used as an ANNEALED homotopy term in the ACS flow: the free residual
    flow has orientation domain walls as local minima at every resolution;
    walls carry covariant-derivative energy while the round structures are
    parallel (S²) or the smoothest in their class, so λ(t)·E with λ → 0
    removes the walls without any reference solution."""
    def single(xs):
        return K_fn(xs.unsqueeze(0))[0]
    DK = torch.func.vmap(torch.func.jacfwd(single))(x)          # (B, n, n, N)
    if geo.kind == 'frame':
        E = geo.frame.frame_fields(x)
        DKa = torch.einsum('bijN,baN->baij', DK, E)
        return (DKa ** 2).sum(dim=(1, 2, 3)).mean()
    P = geo.manifold.projector(x)
    DKa = torch.einsum('bijN,bNa->baij', DK, P)
    DKa = torch.einsum('bik,bakl,blj->baij', P, DKa, P)
    return (DKa ** 2).sum(dim=(1, 2, 3)).mean()


def unit_loss(geo, task, J_fn, x, w_c=10.0):
    if task == 'acs':
        r = acs_residual(geo, J_fn, x)
        return (r ** 2).mean(), {'acs': (r ** 2).mean().item()}
    nij = (nij_norm(geo, J_fn, x) ** 2).mean()
    if geo.case == 's4':                    # no ACS exists: two-term control
        r = acs_residual(geo, J_fn, x)
        acs = (r ** 2).mean()
        return w_c * acs + nij, {'acs': acs.item(), 'nij': nij.item()}
    return nij, {'nij': nij.item()}


# ---------------------------------------------------------------------------
# Polar correction (ACS certification)
# ---------------------------------------------------------------------------

def inv_sqrt_newton_schulz(S, iters=30):
    """Principal S^{-1/2} by the coupled Newton–Schulz iteration on S/2
    (converges quadratically when the spectrum of S lies in (0, 4), i.e.
    whenever sup ‖K²+Id_T‖ < 1). Pure matrix products: differentiable,
    forward-AD and vmap safe, unlike torch.linalg.eig."""
    n = S.shape[-1]
    I = torch.eye(n, dtype=S.dtype, device=S.device).expand_as(S)
    Y, Z = S / 2.0, I
    for _ in range(iters):
        T = 0.5 * (3.0 * I - Z @ Y)
        Y, Z = Y @ T, T @ Z
    return Z / 2.0 ** 0.5


def polar_corrected(geo, K_fn, method='ns'):
    """J* = K · (−K²_T)^{−1/2}: exact J*² = −Id_T whenever K has no real
    eigenvalues on T. S = −K² + (I − Id_T) is NOT symmetric for a generic
    (non-normal) K, so a symmetric eigensolver is wrong; use Newton–Schulz
    (default, AD/vmap safe → usable as a CS base point) or the general
    eigendecomposition ('eig', certification cross-check)."""
    def J_star(x):
        K = K_fn(x)
        IdT = identity_T(geo, x)
        eye = torch.eye(K.shape[-1], dtype=K.dtype).expand_as(K)
        S = -K @ K + (eye - IdT)
        if method == 'eig':
            w, V = torch.linalg.eig(S)
            return K @ (V @ torch.diag_embed(w ** -0.5) @ torch.linalg.inv(V)).real
        return K @ inv_sqrt_newton_schulz(S)
    return J_star


# ---------------------------------------------------------------------------
# Adversarial sup (defect finder): existence is a sup statement
# ---------------------------------------------------------------------------

def _blocks(geo):
    return [m + 1 for m in geo.dims]          # every case is a product of unit spheres


def retract(geo, x):
    parts, off = [], 0
    for n in _blocks(geo):
        blk = x[:, off:off + n]
        parts.append(blk / blk.norm(dim=1, keepdim=True))
        off += n
    return torch.cat(parts, dim=1)


def tangential(geo, x, g):
    parts, off = [], 0
    for n in _blocks(geo):
        xb, gb = x[:, off:off + n], g[:, off:off + n]
        parts.append(gb - (gb * xb).sum(dim=1, keepdim=True) * xb)
        off += n
    return torch.cat(parts, dim=1)


def adversarial_sup(geo, r_fn, n_starts=256, steps=150, lr=0.02):
    """Riemannian gradient ASCENT of a per-sample residual r(x) from many
    starts: returns (max r found, argmax point). A sampled sup can miss a
    small defect core (S⁴!); the ascent finds it: at a degenerate point K
    has a real eigenvalue and ‖K²+Id_T‖ ≥ 1."""
    # starts = the worst points of a large sample (sampling + ascent combined)
    xs = geo.sample(32 * n_starts)
    with torch.no_grad():
        rs = _eval_r(r_fn, xs)
    x = xs[torch.topk(rs, n_starts).indices].detach()
    r_cur = _eval_r(r_fn, x).detach()
    step = torch.full((n_starts, 1), lr, dtype=DTYPE)
    for _ in range(steps):                     # per-sample line-search ascent
        xg = x.detach().requires_grad_(True)
        r = r_fn(xg)
        g, = torch.autograd.grad(r.sum(), xg)
        g = tangential(geo, xg, g)
        x_new = retract(geo, x + step * g / (g.norm(dim=1, keepdim=True) + 1e-12))
        r_new = _eval_r(r_fn, x_new).detach()
        improved = r_new > r_cur
        x = torch.where(improved.unsqueeze(1), x_new.detach(), x)
        r_cur = torch.where(improved, r_new, r_cur)
        step = torch.where(improved.unsqueeze(1), step * 1.2, step * 0.5).clamp(1e-6, 0.3)
    k = int(r_cur.argmax())
    return r_cur[k].item(), x[k].detach()


def _eval_r(r_fn, x):
    """Evaluate a residual that may need autograd internally (Nijenhuis via
    jacfwd) without building a graph on x."""
    return r_fn(x.detach())


# ---------------------------------------------------------------------------
# Chart-protocol evaluation (the paper's original metric)
# ---------------------------------------------------------------------------

def ambient_field(geo, J_fn):
    """Frame fields are bridged to ambient tangential fields:
    J_amb = Eᵀ J_frame (E Eᵀ)⁻¹ E."""
    if geo.kind == 'ambient':
        return J_fn

    def J_amb(p):
        E = geo.frame.frame_fields(p)                                       # (B,d,N)
        G = E @ E.transpose(1, 2)
        Jf = J_fn(p)
        return E.transpose(1, 2) @ Jf @ torch.linalg.inv(G) @ E
    return J_amb


def factor_embedding(x, chart):
    if chart == 0:
        return EvenSphere.ball_to_sphere(x, from_north=True)
    p = EvenSphere.ball_to_sphere(x, from_north=False)
    return torch.cat([-p[:, :-1], p[:, -1:]], dim=1)


def chart_field(J_amb, dims, charts):
    def emb_single(xc):
        parts, off = [], 0
        for m, ch in zip(dims, charts):
            parts.append(factor_embedding(xc[off:off + m].unsqueeze(0), ch).squeeze(0))
            off += m
        return torch.cat(parts)

    def J_of(x_chart):
        p = torch.func.vmap(emb_single)(x_chart)
        D = torch.func.vmap(torch.func.jacfwd(emb_single))(x_chart)
        Ja = J_amb(p)
        D_pinv = torch.linalg.inv(D.transpose(1, 2) @ D) @ D.transpose(1, 2)
        return D_pinv @ Ja @ D
    return J_of


def chart_nij(J_of, x_chart, dims, chunk=400):
    off, mask = 0, torch.ones(len(x_chart), dtype=torch.bool)
    for m in dims:
        mask &= x_chart[:, off:off + m].norm(dim=1) <= R_CUT
        off += m
    xa = x_chart[mask]
    d = sum(dims)
    norms = []
    for s in range(0, len(xa), chunk):
        xr = xa[s:s + chunk].detach().clone().requires_grad_(True)
        J = J_of(xr)
        dJ = torch.zeros(len(xr), d, d, d, dtype=DTYPE)
        for i in range(d):
            for j in range(d):
                g, = torch.autograd.grad(J[:, i, j].sum(), xr, retain_graph=True)
                dJ[:, i, j, :] = g
        N = _nij_tensor(J.detach(), dJ)
        norms.append(torch.sqrt((N ** 2).sum(dim=(1, 2, 3)) + 1e-20))
    v = torch.cat(norms) if norms else torch.zeros(1, dtype=DTYPE)
    return v.mean().item(), v.max().item()


def chart_protocol(geo, J_fn, n_test):
    """Returns (chart_cs = mean hard-cutoff Nij over patches, chart_acs = overlap)."""
    J_amb = ambient_field(geo, J_fn)
    dims = geo.dims
    if len(dims) == 1:
        m = dims[0]
        sphere = EvenSphere(dimension=m)
        x1 = sphere.sample_points_patch1(n_test).to(DTYPE)
        x2 = EvenSphere.patch_change_coords(x1)
        J1f, J2f = chart_field(J_amb, dims, (0,)), chart_field(J_amb, dims, (1,))
        n1, _ = chart_nij(J1f, x1, dims)
        n2, _ = chart_nij(J2f, x2, dims)
        with torch.no_grad():
            ov, _ = OverlapLoss(overlap_size=OVERLAP_SIZE).to(DTYPE)(
                J1f(x1), J2f(x2), x1, x2, hard_cutoff=True)
        return 0.5 * (n1 + n2), ov.item()
    ma, mb = dims
    sp = SphereProduct(ma, mb)
    xa0, xb0 = sp.sample_points(n_test)
    xa0, xb0 = xa0.to(DTYPE), xb0.to(DTYPE)
    pa = {0: xa0, 1: EvenSphere.patch_change_coords(xa0)}
    pb = {0: xb0, 1: EvenSphere.patch_change_coords(xb0)}
    nijs, J_dict = [], {}
    for ia, ib in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        J_of = chart_field(J_amb, dims, (ia, ib))
        xc = torch.cat([pa[ia], pb[ib]], dim=1)
        n_mean, _ = chart_nij(J_of, xc, dims)
        nijs.append(n_mean)
        with torch.no_grad():
            J_dict[(ia, ib)] = J_of(xc)
    with torch.no_grad():
        ov, _ = ProductOverlapLoss(ma, mb, OVERLAP_SIZE, OVERLAP_SIZE).to(DTYPE)(
            J_dict, xa0, xb0, pa[1], pb[1], hard_cutoff=True)
    return sum(nijs) / 4.0, ov.item()


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def dense_eval(geo, task, J_fn, n, chunk=500):
    accs, nijs = [], []
    for s in range(0, n, chunk):
        x = geo.sample(min(chunk, n - s))
        with torch.no_grad():
            accs.append(acs_residual(geo, J_fn, x))
        nijs.append(nij_norm(geo, J_fn, x).detach())
    acs, nij = torch.cat(accs), torch.cat(nijs)
    return {'acs_mean': acs.mean().item(), 'acs_sup': acs.max().item(),
            'nij_mean': nij.mean().item(), 'nij_sup': nij.max().item()}


def full_report(geo, task, J_fn, meta, n_dense, n_chart, adv_starts=256, adv_steps=150):
    """Dense + ADVERSARIAL evaluation. Certification of an ACS requires the
    adversarial sup of ‖K²+Id_T‖ to be < 1 (then the polar correction is an
    exact ACS); the adversarial sup of ‖N‖ is reported for the final field."""
    rep = {'dense': dense_eval(geo, task, J_fn, n_dense)}
    field_for_chart = J_fn
    if task == 'acs' or not meta['exact_acs']:
        adv, _ = adversarial_sup(geo, lambda x: acs_residual(geo, J_fn, x), adv_starts, adv_steps)
        rep['adv_acs_sup'] = adv
        if adv < 1.0:
            J_star = polar_corrected(geo, J_fn)
            corr = dense_eval(geo, task, J_star, n_dense)
            corr['adv_acs_sup'], _ = adversarial_sup(
                geo, lambda x: acs_residual(geo, J_star, x), adv_starts // 2, adv_steps // 2)
            rep['corrected'] = corr
            rep['certified_exact_acs'] = corr['adv_acs_sup'] < 1e-8
            field_for_chart = J_star
        else:
            rep['certified_exact_acs'] = False
    adv_n, _ = adversarial_sup(geo, lambda x: nij_norm(geo, field_for_chart, x),
                               adv_starts // 2, adv_steps // 2)
    rep['adv_nij_sup'] = adv_n
    cs, acs = chart_protocol(geo, field_for_chart, n_chart)
    rep['chart_cs'], rep['chart_acs'] = cs, acs
    return rep


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', required=True, choices=ALL_CASES)
    ap.add_argument('--task', required=True, choices=['acs', 'cs'])
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--out', required=True)
    ap.add_argument('--steps', type=int, default=None)
    ap.add_argument('--hidden', default=None)
    ap.add_argument('--stages', type=int, default=None)
    ap.add_argument('--rar', action='store_true')
    ap.add_argument('--lbfgs', type=int, default=200)
    ap.add_argument('--wall-minutes', type=float, default=4200.0)
    ap.add_argument('--eval-every', type=int, default=2000)
    ap.add_argument('--ckpt-every', type=int, default=1000)
    ap.add_argument('--n-dense', type=int, default=6000)
    ap.add_argument('--n-chart', type=int, default=3000)
    ap.add_argument('--poly-steps', type=int, default=3000,
                    help='acs: spectral-curriculum phase 1 (degree-2 polynomial field) steps; 0 = plain MLP')
    ap.add_argument('--base-from', default=None,
                    help='cs: directory of a finished acs unit; its certified polar-corrected field is the base ACS')
    ap.add_argument('--feat', default='auto',
                    help="acs phase 1 feature class: 'auto' climbs the harmonic ladder " + str(LADDER))
    ap.add_argument('--lam0', type=float, default=0.0, help='acs phase 1: initial Dirichlet weight (0 = off)')
    ap.add_argument('--lam1', type=float, default=1e-3, help='acs phase 1: final Dirichlet weight')
    ap.add_argument('--adv-starts', type=int, default=256)
    ap.add_argument('--adv-steps', type=int, default=150)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    torch.manual_seed(args.seed)
    geo = Geometry(args.case)
    d, N = geo.d, geo.N
    steps = args.steps or (10000 if args.task == 'acs' else 40000)
    hidden = [int(h) for h in args.hidden.split(',')] if args.hidden else \
        ([256, 256, 256, 256] if d >= 6 else [128, 128, 128])
    stages = args.stages or (2 if args.task == 'cs' else 1)
    n_mat = d if geo.kind == 'frame' else N
    out_scale = 0.3 if args.task == 'acs' else 0.5
    tag = f"{args.case}/{args.task}/s{args.seed} [{geo.kind}]"

    curriculum = args.task == 'acs' and args.poly_steps > 0
    blocks = [m + 1 for m in geo.dims]
    feat_used = None
    if curriculum:
        feat_used = args.feat if args.feat != 'auto' else LADDER[0]
        model = CurriculumNet(N, n_mat, hidden, out_scale_poly=out_scale,
                              feat=feat_used, blocks=blocks).to(DTYPE)
    else:
        model = Net(N, n_mat, stages, hidden, out_scale).to(DTYPE)

    base_fn, base_info = None, None
    if args.task == 'cs' and args.base_from:
        # Unsupervised end to end: the base ACS is the acs unit's own field.
        K_base = load_acs_field(geo, args.base_from)
        adv, _ = adversarial_sup(geo, lambda x: acs_residual(geo, K_base, x),
                                 args.adv_starts, args.adv_steps)
        base_info = {'dir': args.base_from, 'adv_acs_sup_raw': adv}
        if adv >= 1.0:
            base_info['status'] = 'no_certified_base'
            result = {'case': args.case, 'task': args.task, 'seed': args.seed,
                      'formulation': geo.kind, 'status': 'no_certified_base',
                      'base': base_info, 'table_before': TABLE_BEFORE[args.case][args.task]}
            with open(os.path.join(args.out, 'result.json'), 'w') as f:
                json.dump(result, f, indent=2)
            print(f"[{tag}] NO CERTIFIED BASE ACS (adversarial sup {adv:.3f} >= 1): CS stage not run",
                  flush=True)
            return
        base_fn = polar_corrected(geo, K_base)
        xb = geo.sample(400)
        base_info['status'] = 'certified'
        base_info['base_acs_residual_sup'] = acs_residual(geo, base_fn, xb).max().item()
        base_info['base_nij_mean'] = nij_norm(geo, base_fn, xb).mean().item()
        print(f"[{tag}] base ACS from {args.base_from}: adversarial raw sup {adv:.3f}, "
              f"polar-corrected residual {base_info['base_acs_residual_sup']:.1e}, "
              f"Nij mean {base_info['base_nij_mean']:.3e}", flush=True)
    J_fn, meta = make_field(geo, args.task, model, base_fn=base_fn)

    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=steps, eta_min=1e-5)
    state_path = os.path.join(args.out, 'state.pt')
    start_step, elapsed_prev, init_info, poly_done = 0, 0.0, None, False
    if os.path.exists(state_path):
        st = torch.load(state_path, map_location='cpu', weights_only=False)
        if curriculum and st.get('feat'):
            feat_used = st['feat']
            model.poly = PolyNet(N, n_mat, out_scale, feat=feat_used, blocks=blocks).to(DTYPE)
            J_fn, meta = make_field(geo, args.task, model, base_fn=base_fn)
            opt = torch.optim.Adam(model.parameters(), lr=1e-3)
            sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=steps, eta_min=1e-5)
        model.load_state_dict(st['model']); opt.load_state_dict(st['opt'])
        sched.load_state_dict(st['sched'])
        start_step, elapsed_prev, init_info = st['step'], st['elapsed'], st.get('init')
        poly_done = st.get('poly_done', False)
        print(f"[{tag}] resumed at step {start_step} ({elapsed_prev/60:.1f} min used)", flush=True)
    if init_info is None:
        x0 = geo.sample(300)
        _, init_info = unit_loss(geo, args.task, J_fn, x0)
        print(f"[{tag}] init: {init_info}", flush=True)

    progress = open(os.path.join(args.out, 'progress.jsonl'), 'a')
    t0 = time.time()
    elapsed = lambda: elapsed_prev + (time.time() - t0)
    wall_s = args.wall_minutes * 60.0
    nonfinite, step = 0, start_step
    stopped_by_wall = False

    def save_state(step_):
        torch.save({'model': model.state_dict(), 'opt': opt.state_dict(),
                    'sched': sched.state_dict(), 'step': step_, 'poly_done': poly_done,
                    'feat': feat_used,
                    'elapsed': elapsed(), 'init': init_info}, state_path)

    # Phase 1 (ACS): polynomial field only (MLP frozen), climbing the harmonic
    # ladder; optional annealed intrinsic Dirichlet term λ(t): λ0 → λ1.
    ladder_log = []
    if curriculum and not poly_done and start_step == 0:
        classes = [args.feat] if args.feat != 'auto' else \
            [f for f in LADDER if ('bilin' not in f or len(geo.dims) == 2)]
        best = None
        for feat in classes:
            model.poly = PolyNet(N, n_mat, out_scale, feat=feat, blocks=blocks).to(DTYPE)
            opt_p = torch.optim.Adam(model.poly.parameters(), lr=5e-3)
            sched_p = torch.optim.lr_scheduler.CosineAnnealingLR(opt_p, T_max=args.poly_steps, eta_min=2e-4)
            for ps in range(args.poly_steps):
                x = geo.sample(BATCH)
                loss, info = unit_loss(geo, args.task, J_fn, x)
                if args.lam0 > 0:
                    lam = args.lam0 * (args.lam1 / args.lam0) ** (ps / max(1, args.poly_steps - 1))
                    dir_e = dirichlet_energy(geo, J_fn, x)
                    info['dirichlet'], info['lam'] = dir_e.item(), lam
                    loss = loss + lam * dir_e
                if not bool(torch.isfinite(loss)):
                    continue
                opt_p.zero_grad(); loss.backward(); opt_p.step(); sched_p.step()
            ev = dense_eval(geo, args.task, J_fn, 1000)
            adv, _ = adversarial_sup(geo, lambda x: acs_residual(geo, J_fn, x),
                                     max(32, args.adv_starts // 4), max(30, args.adv_steps // 3))
            entry = {'phase': 'poly', 'feat': feat, 'n_params': sum(p.numel() for p in model.poly.parameters()),
                     'elapsed_min': elapsed() / 60.0, 'adv_acs_sup': adv, **ev}
            ladder_log.append(entry)
            progress.write(json.dumps(entry) + '\n'); progress.flush()
            print(f"[{tag}] ladder {feat:10s} ({entry['n_params']} params) | held-out {ev} | adversarial sup {adv:.3f}"
                  f" -> {'ACCEPT' if adv < 1 else 'reject'}", flush=True)
            if best is None or adv < best[0]:
                best = (adv, feat, {k: v.clone() for k, v in model.poly.state_dict().items()})
            if adv < 1.0:
                break
        if best is not None and best[0] >= 1.0:          # nothing certified: keep the least bad
            model.poly = PolyNet(N, n_mat, out_scale, feat=best[1], blocks=blocks).to(DTYPE)
            model.poly.load_state_dict(best[2])
        feat_used = model.poly.feat
        J_fn, meta = make_field(geo, args.task, model, base_fn=base_fn)
        opt = torch.optim.Adam(model.parameters(), lr=1e-3)
        sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=steps, eta_min=1e-5)
        poly_done = True
        save_state(0)

    for step in range(start_step, steps):
        if elapsed() > 0.85 * wall_s:
            stopped_by_wall = True
            break
        if args.rar and args.task == 'cs':
            xc = geo.sample(2 * BATCH)
            r = nij_norm(geo, J_fn, xc).detach()
            idx = torch.topk(r, BATCH // 2).indices
            x = torch.cat([xc[idx], geo.sample(BATCH - BATCH // 2)], dim=0)
        else:
            x = geo.sample(BATCH)
        loss, info = unit_loss(geo, args.task, J_fn, x)
        if not bool(torch.isfinite(loss)):
            nonfinite += 1
            if nonfinite >= 200:
                print(f"[{tag}] aborting: 200 consecutive non-finite losses", flush=True)
                break
            continue
        nonfinite = 0
        opt.zero_grad(); loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step(); sched.step()
        if (step + 1) % args.ckpt_every == 0:
            save_state(step + 1)
        if (step + 1) % args.eval_every == 0 or step == steps - 1:
            ev = dense_eval(geo, args.task, J_fn, 1500)
            rec = {'step': step + 1, 'elapsed_min': elapsed() / 60.0, 'batch': info, **ev}
            progress.write(json.dumps(rec) + '\n'); progress.flush()
            print(f"[{tag}] step {step+1:6d} | batch {info} | held-out {ev}", flush=True)
    save_state(step + 1 if not stopped_by_wall else step)

    # L-BFGS polish (time-guarded)
    if elapsed() < 0.92 * wall_s and args.lbfgs > 0:
        xf = geo.sample(768)
        lb = torch.optim.LBFGS(model.parameters(), lr=1.0, max_iter=args.lbfgs,
                               history_size=100, line_search_fn='strong_wolfe')

        def closure():
            lb.zero_grad()
            l, _ = unit_loss(geo, args.task, J_fn, xf)
            l.backward()
            return l
        try:
            lb.step(closure)
        except Exception as e:                     # keep the Adam result on failure
            print(f"[{tag}] LBFGS failed ({type(e).__name__}); keeping Adam weights", flush=True)
        save_state(steps)

    torch.save({'model': model.state_dict(), 'case': args.case, 'task': args.task,
                'seed': args.seed, 'hidden': hidden, 'stages': stages, 'kind': geo.kind,
                'model_type': 'curriculum' if curriculum else 'mlp', 'feat': feat_used,
                'base_from': args.base_from},
               os.path.join(args.out, 'final.pt'))

    rep = full_report(geo, args.task, J_fn, meta, args.n_dense, args.n_chart,
                      args.adv_starts, args.adv_steps)
    result = {'case': args.case, 'task': args.task, 'seed': args.seed,
              'formulation': geo.kind, 'hidden': hidden, 'stages': stages,
              'model_type': 'curriculum' if curriculum else 'mlp', 'feat': feat_used,
              'ladder': ladder_log, 'base': base_info,
              'steps_done': min(step + 1, steps), 'steps_planned': steps,
              'stopped_by_wall': stopped_by_wall, 'elapsed_min': elapsed() / 60.0,
              'init': init_info, 'exact_acs_by_construction': meta['exact_acs'],
              'table_before': TABLE_BEFORE[args.case][args.task], **rep}
    with open(os.path.join(args.out, 'result.json'), 'w') as f:
        json.dump(result, f, indent=2)
    verdict = ('CERTIFIED exact ACS' if rep.get('certified_exact_acs') else
               ('exact by construction' if meta['exact_acs'] else 'NOT certified (defects)'))
    print(f"[{tag}] RESULT: {verdict} | adv acs sup {rep.get('adv_acs_sup', float('nan')):.3e} "
          f"| adv Nij sup {rep['adv_nij_sup']:.3e} | dense {rep['dense']} "
          f"| chart {args.task} = {rep['chart_cs'] if args.task == 'cs' else rep['chart_acs']:.3e} "
          f"(table before {TABLE_BEFORE[args.case][args.task]})", flush=True)


if __name__ == '__main__':
    main()
