"""Uniform chart-free PINN campaign — one (case, task, seed) unit.

Formulation (--formulation, default 'ambient'): ONE common description for
all registered cases — every case is a product of unit spheres in ℝ^N, fields are
J = P·M(p)·P with P the tangential projector, and the Nijenhuis tensor is
assembled from ambient brackets (src/loss/ambient_loss.py). No charts. The
global-frame Nijenhuis remains available for the parallelizable cases
(--formulation frame / auto).

Tasks (unsupervised PINN, identical across formulations):
  acs  free endomorphism field K (frame: K = net(x); ambient: K = P·net(x)·P),
       loss = mean ‖K² + Id_T‖²_F. From scratch, no base point. Harmonic ladder
       of polynomial classes, each tried plain → annealed Dirichlet → orientation
       penalty (--attempts). Certified a posteriori by the polar correction
       J* = K(−K²_T)^{−1/2} when sup ‖K²+Id_T‖ < 1.
  cs   --cs-mode poly (recommended): free polynomial field K = P·M_W(p)·P from
       the unit's own certified ACS, loss E‖N_K‖² + w E‖K²+P‖², degree ladder,
       Levenberg–Marquardt finisher (extra/poly_cs.py); certified and
       polar-corrected like an ACS unit.
  cs   --cs-mode conj (original): exact-ACS parametrization
       J = g_s⋯g_1 J_base g_1⁻¹⋯g_s⁻¹, g_k = e^{B_k} (or --conj cayley / affine)
       (B_k tangential net outputs, initialised at O(1) ⇒ a generic
       NON-integrable start), loss = mean ‖N(J)‖² only. J_base is normally the
       certified ACS learned by an earlier unit (--base-from), with canonical
       or classical explicit structures retained as optional fallbacks;
       s4 has no ACS, so its cs unit runs the free-field two-term loss
       (w·‖K²+P‖² + ‖N‖²) as the "nothing exists" control.

Every unit is restartable (state.pt), time-guarded (--wall-minutes), logs
progress.jsonl, and ends with a dense evaluation + the ORIGINAL chart-protocol
scores (frame fields are bridged to ambient fields, then to ball charts), in
result.json.

Usage:
  python extra/campaign_unit.py --case s3xs3 --task acs --seed 7 --out ACS_DIR
  python extra/campaign_unit.py --case s3xs3 --task cs --seed 7 --cs-mode poly \
      --base-from ACS_DIR --out DIR [--poly-degrees 2,3]
  python extra/campaign_unit.py --case s3xs3 --task cs --seed 7 --out DIR \
      [--steps 40000] [--hidden 256,256,256,256] [--stages 2] [--rar] [--tail] \
      [--conj exp] [--fast-cs] [--wall-minutes 4200] [--eval-every 2000] [--ckpt-every 1000]
  bash extra/local_campaign.sh [SEED]      # every case, ACS then poly-mode CS, then the table
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import os
import sys
import json
import time
import argparse
import itertools
import numpy as np
import torch
import torch.nn as nn

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from src.geometry.frames import FRAMES
from src.geometry.explicit_acs import (ALL_CASES, factor_dims, is_frame_case,
                                       ambient_manifold, base_acs)
from src.loss.frame_loss import frame_nij_norm
from src.loss.ambient_loss import (ambient_nij_norm, compute_J_and_DJ,
                                   compute_P_and_DP)
from src.geometry.sphere import EvenSphere, RADIAL_MIDPOINT, sphere_uniform_sample
from src.loss.sphere_loss import _nij_tensor

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
    """formulation: 'ambient' (default; every case as a product of unit spheres
    with the tangential projector — one common description for all cases),
    'frame' (global-frame Nijenhuis, parallelizable cases only) or 'auto'
    (frame where available)."""

    def __init__(self, case, formulation='ambient'):
        self.case = case
        self.dims = factor_dims(case)
        self.d = sum(self.dims)
        use_frame = is_frame_case(case) and formulation in ('frame', 'auto')
        if formulation == 'frame' and not is_frame_case(case):
            raise ValueError(f"'{case}' is not parallelizable: no global frame; use --formulation ambient")
        if use_frame:
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
             'bilin' : x_A ⊗ x_B across every pair of distinct product factors
                       (`blocks` = ambient block sizes of all factors)
             'lin+bilin': x, x_A ⊗ x_B"""
        super().__init__()
        self.feat, self.n_in = feat, n_in
        iu = torch.triu_indices(n_in, n_in)
        self.register_buffer('iu', iu)
        if blocks is not None and sum(blocks) != n_in:
            raise ValueError(f"ambient block sizes must sum to n_in={n_in}, got {blocks!r}")
        if blocks is not None and len(blocks) >= 2:
            offsets = [0]
            for width in blocks:
                offsets.append(offsets[-1] + width)
            pairs = []
            for a, b in itertools.combinations(range(len(blocks)), 2):
                ia, ib = torch.meshgrid(
                    torch.arange(offsets[a], offsets[a + 1]),
                    torch.arange(offsets[b], offsets[b + 1]),
                    indexing='ij',
                )
                pairs.append(torch.stack([ia.reshape(-1), ib.reshape(-1)]))
            self.register_buffer('bil', torch.cat(pairs, dim=1))
        else:
            self.bil = None
        if feat in ('bilin', 'lin+bilin') and self.bil is None:
            raise ValueError(f"feature class {feat!r} requires at least two ambient factor blocks")
        # Three-factor products x_A,a x_B,b x_C,c (one coordinate from each of three
        # distinct factors): the multidegree-(1,1,1) harmonics of a product with
        # >= 3 factors. Same convention as the s4t2 branch ('trilin', 'all2+trilin').
        triples = []
        if blocks is not None and len(blocks) >= 3:
            for factors in itertools.combinations(range(len(blocks)), 3):
                grids = torch.meshgrid(*[torch.arange(offsets[a], offsets[a + 1]) for a in factors],
                                       indexing='ij')
                triples.append(torch.stack([g.reshape(-1) for g in grids]))
        self.register_buffer('tri', torch.cat(triples, dim=1) if triples else
                             torch.empty(3, 0, dtype=torch.long), persistent=False)
        if feat in ('trilin', 'all2+trilin') and not triples:
            raise ValueError(f"feature class {feat!r} requires at least three ambient factor blocks")
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
        if self.feat in ('trilin', 'all2+trilin'):
            tri = x[:, self.tri[0]] * x[:, self.tri[1]] * x[:, self.tri[2]]
            return tri if self.feat == 'trilin' else torch.cat([one, x, quad, tri], dim=1)
        if self.feat == 'c':
            return one
        if self.feat == 'lin':
            return x
        if self.feat == 'lin+c':
            return torch.cat([one, x], dim=1)
        if self.feat == 'quad_h':
            return quad
        if self.feat == 'bilin':
            return x[:, self.bil[0]] * x[:, self.bil[1]]
        if self.feat == 'lin+bilin':
            bil = x[:, self.bil[0]] * x[:, self.bil[1]]
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


def conjugator(B, conj='exp'):
    """(g, g⁻¹) from a generator field B (B, n, n). 'exp': matrix exponential
    (the original parametrization). 'cayley': g = (I − B/2)⁻¹(I + B/2),
    g⁻¹ = (I + B/2)⁻¹(I − B/2). Both commute with anything B commutes with
    (in particular with P), so the conjugated field stays an exact ACS. The
    Cayley map is a rational function: forward-mode AD through it is ~25x
    cheaper than through matrix_exp (a batched solve instead of the
    Padé/scaling-squaring graph), which is what makes the CS step affordable."""
    if conj == 'exp':
        return torch.matrix_exp(B), torch.matrix_exp(-B)
    # NOTE: torch.linalg.inv, not torch.linalg.solve — under vmap + forward AD
    # (torch 2.14) solve returns wrong derivatives for batch elements > 0.
    I = torch.eye(B.shape[-1], dtype=B.dtype, device=B.device).expand_as(B)
    if conj == 'affine':
        # g = I + B: the WHOLE group GL⁺ (g − I is unrestricted), no Lie-algebra
        # chart at all; det g = 0 is a barrier of the loss (‖N‖² ~ ‖g⁻¹‖²), so
        # the flow stays in the identity component. Needs a small init.
        g = I + B
        return g, torch.linalg.inv(g)
    A, C = I - 0.5 * B, I + 0.5 * B
    Ai = torch.linalg.inv(A)
    return Ai @ C, torch.linalg.inv(C) @ A


# Harmonic ladder: lowest spherical-harmonic content first. A low class has no
# room for defect cores, so the residual flow settles on a global structure
# (or, on S⁴, stays loudly non-zero). Product-only classes are skipped on
# single spheres. Acceptance = adversarial sup ‖K²+Id_T‖ < 1 (certifiable).
LADDER = ['c', 'lin', 'bilin', 'quad_h', 'lin+c', 'lin+bilin', 'all2', 'all2+trilin']


def ladder_classes(n_factors):
    """Ladder rungs that exist for a product of `n_factors` spheres."""
    return [f for f in LADDER
            if ('bilin' not in f or n_factors >= 2) and ('trilin' not in f or n_factors >= 3)]


def make_field(geo, task, model, base_fn=None, conj='exp'):
    """Returns (J_fn, meta). J_fn maps ambient points → matrix field:
    frame cases: (B, d, d) in frame components; ambient cases: (B, N, N).
    `base_fn` (CS task): the base ACS field, e.g. the PINN's own certified
    ACS from the acs unit; default = canonical J₀ / classical explicit ACS.
    `conj`: 'exp' or 'cayley' conjugation (see `conjugator`)."""
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
                g, gi = conjugator(raw[:, s], conj)
                J = g @ J @ gi
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
            g, gi = conjugator(P @ raw[:, s] @ P, conj)
            J = g @ J @ gi
        return J
    return J_fn, {'exact_acs': True}


def make_conjugation(geo, model, conj='exp'):
    """Total conjugator field g(x) = g_S ⋯ g_1 (ambient CS task), so that
    J = g J_base g⁻¹. Used by the fast CS loss."""
    proj = geo.manifold.projector

    def g_fn(x):
        P = proj(x)
        raw = model(x)
        g = None
        for s in range(model.stages):
            gs, _ = conjugator(P @ raw[:, s] @ P, conj)
            g = gs if g is None else gs @ g
        return g
    return g_fn


def cs_nij_norm_fast(geo, g_fn, base_fn, x):
    """Per-sample ‖N_J‖ for J = g J_b g⁻¹ with the derivative decomposed:
    DJ = Dg J_b g⁻¹ + g DJ_b g⁻¹ − g J_b g⁻¹ Dg g⁻¹. The frozen base (which
    contains the Newton–Schulz polar correction) is differentiated once
    WITHOUT an autograd graph; forward-mode AD only runs through the
    (cheap) conjugator net. Identical to `nij_norm` up to round-off, but
    only valid for parameter gradients (no x-gradient through J_b), so it is
    used for the training loss only."""
    with torch.no_grad():
        Jb, DJb = compute_J_and_DJ(base_fn, x)
        P, DP = compute_P_and_DP(geo.manifold, x)
    g, Dg = compute_J_and_DJ(g_fn, x)
    gi = torch.linalg.inv(g)
    J = g @ Jb @ gi
    DJ = torch.einsum('bikc,bkl,blj->bijc', Dg, Jb, gi) \
        + torch.einsum('bik,bklc,blj->bijc', g, DJb, gi) \
        - torch.einsum('bik,bklc,blj->bijc', J, Dg, gi)
    return ambient_nij_norm(J, DJ, P, DP)


def load_acs_field(geo, run_dir):
    """Rebuild the frozen free ACS field K(x) of a finished acs unit."""
    ck = torch.load(os.path.join(run_dir, 'final.pt'), map_location='cpu', weights_only=False)
    if ck.get('case') not in (None, geo.case):
        raise ValueError(f"base ACS in {run_dir} is for case {ck['case']!r}, expected {geo.case!r}")
    if ck.get('task') not in (None, 'acs'):
        raise ValueError(f"base checkpoint in {run_dir} has task {ck['task']!r}, expected 'acs'")
    if ck.get('kind', geo.kind) != geo.kind:
        raise ValueError(f"base ACS in {run_dir} was trained in the '{ck['kind']}' formulation, "
                         f"this unit runs '{geo.kind}': use the same --formulation for both stages")
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


def oriented_tangent_frame(geo, x):
    """(B, N, d) orthonormal frame of T_pM, positively oriented for the fixed
    global orientation in which (n_1, ..., n_k, E) is a positive basis of R^N
    (n_A = unit outward normal of factor A). Ambient formulation only."""
    P = geo.manifold.projector(x)
    _, evecs = torch.linalg.eigh(P)
    E = evecs[:, :, -geo.d:]
    normals, off = [], 0
    for n in _blocks(geo):
        v = torch.zeros_like(x)
        v[:, off:off + n] = x[:, off:off + n] / x[:, off:off + n].norm(dim=1, keepdim=True)
        normals.append(v.unsqueeze(2))
        off += n
    sgn = torch.sign(torch.linalg.det(torch.cat(normals + [E], dim=2)))
    return torch.cat([E[:, :, :1] * sgn.view(-1, 1, 1), E[:, :, 1:]], dim=2)


def pfaffian(A):
    """Batched Pfaffian of antisymmetric (..., 2n, 2n) matrices (expansion along
    the first row; 15 terms for 2n = 6)."""
    d = A.shape[-1]
    if d == 0:
        return torch.ones(A.shape[:-2], dtype=A.dtype)
    if d == 2:
        return A[..., 0, 1]
    out = torch.zeros(A.shape[:-2], dtype=A.dtype)
    for j in range(1, d):
        keep = [k for k in range(1, d) if k != j]
        sub = A[..., keep, :][..., :, keep]
        out = out + (-1) ** (j + 1) * A[..., 0, j] * pfaffian(sub)
    return out


def orientation_index(geo, K_fn, x):
    """o(p) = Pf(E^T (K^T - K) E) / 2^{d/2}. For an exact ACS A, omega = A^T - A
    = h A with h = I + A^T A positive and A-invariant, so omega is nondegenerate,
    |o| = sqrt(det(I + A^T A)) / 2^{d/2} >= 1 (singular values of A pair as s, 1/s)
    and sign(o) is the complex orientation of A. At an orientation wall K ~ 0 and
    o ~ 0."""
    E = oriented_tangent_frame(geo, x)
    K = K_fn(x)
    om = torch.einsum('bpa,bpq,bqc->bac', E, K.transpose(1, 2) - K, E)
    return pfaffian(om) / 2.0 ** (geo.d // 2)


def unit_loss(geo, task, J_fn, x, w_c=10.0, fast=None):
    """`fast`: optional (g_fn, base_fn) pair — CS training loss through the
    decomposed derivative (cs_nij_norm_fast) instead of AD through the base."""
    if task == 'acs':
        r = acs_residual(geo, J_fn, x)
        return (r ** 2).mean(), {'acs': (r ** 2).mean().item()}
    if fast is not None:
        nij = (cs_nij_norm_fast(geo, fast[0], fast[1], x) ** 2).mean()
        return nij, {'nij': nij.item()}
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
    """Evaluate a global field in all ``2**k`` product charts.

    Returns ``(chart_cs, chart_acs)``: the mean hard-cutoff coordinate
    Nijenhuis norm over every chart and the mean consistency residual over
    every pair of overlapping charts.  The implementation is deliberately
    factor-count agnostic; in particular, individual factors may be odd
    dimensional even though the total dimension of an ACS candidate is even.
    """
    J_amb = ambient_field(geo, J_fn)
    dims = geo.dims
    factor_coords = []
    for m in dims:
        x0 = torch.as_tensor(sphere_uniform_sample(n_test, m), dtype=DTYPE)
        factor_coords.append((x0, EvenSphere.patch_change_coords(x0)))

    patches = list(itertools.product((0, 1), repeat=len(dims)))
    nijs, J_dict = [], {}
    for patch in patches:
        J_of = chart_field(J_amb, dims, patch)
        xc = torch.cat([factor_coords[a][patch[a]] for a in range(len(dims))], dim=1)
        n_mean, _ = chart_nij(J_of, xc, dims)
        nijs.append(n_mean)
        with torch.no_grad():
            J_dict[patch] = J_of(xc)
    ov = product_chart_overlap(J_dict, factor_coords, dims, patches)
    return sum(nijs) / len(nijs), ov.item()


def _batch_block_diag(blocks):
    """Batched block diagonal assembly without in-place tensor mutation."""
    batch = blocks[0].shape[0]
    total = sum(block.shape[1] for block in blocks)
    rows, off = [], 0
    for block in blocks:
        width = block.shape[1]
        left = block.new_zeros(batch, width, off)
        right = block.new_zeros(batch, width, total - off - width)
        rows.append(torch.cat([left, block, right], dim=2))
        off += width
    return torch.cat(rows, dim=1)


def product_chart_overlap(J_dict, factor_coords, dims, patches=None,
                          overlap_size=OVERLAP_SIZE):
    """Mean tensor-transition residual over every product-chart pair.

    A pair may flip one or several sphere factors.  Its transition Jacobian
    is assembled block diagonally from the per-factor patch changes, so this
    is the direct ``k``-factor generalisation of ``ProductOverlapLoss``.
    """
    if patches is None:
        patches = list(itertools.product((0, 1), repeat=len(dims)))
    batch = next(iter(J_dict.values())).shape[0]
    dtype = next(iter(J_dict.values())).dtype
    device = next(iter(J_dict.values())).device
    r_outer = RADIAL_MIDPOINT + overlap_size
    r_inner = (1.0 - r_outer) / (1.0 + r_outer)
    r_floor = max(0.05, r_inner - 2.0 * overlap_size)
    pair_losses = []

    for source, target in itertools.combinations(patches, 2):
        flipped = [a for a in range(len(dims)) if source[a] != target[a]]
        src_blocks, tgt_blocks = [], []
        weight = torch.ones(batch, dtype=dtype, device=device)
        valid = torch.ones(batch, dtype=torch.bool, device=device)
        for a, m in enumerate(dims):
            if a in flipped:
                xs = factor_coords[a][source[a]].to(device=device, dtype=dtype)
                xt = factor_coords[a][target[a]].to(device=device, dtype=dtype)
                src_blocks.append(EvenSphere.patch_change_jacobian(xs))
                tgt_blocks.append(EvenSphere.patch_change_jacobian(xt))
                rho_s, rho_t = xs.norm(dim=1), xt.norm(dim=1)
                valid &= (rho_s >= r_floor) & (rho_t >= r_floor)
                weight *= ((rho_s >= r_inner) & (rho_s <= r_outer)).to(dtype)
            else:
                eye = torch.eye(m, dtype=dtype, device=device).unsqueeze(0).expand(batch, m, m)
                src_blocks.append(eye)
                tgt_blocks.append(eye)

        F_src = _batch_block_diag(src_blocks)
        F_tgt = _batch_block_diag(tgt_blocks)
        J_src, J_tgt = J_dict[source], J_dict[target]
        J_src_to_tgt = F_src @ J_src @ torch.linalg.inv(F_src)
        J_tgt_to_src = F_tgt @ J_tgt @ torch.linalg.inv(F_tgt)
        norm_fwd = torch.sqrt(((J_src_to_tgt - J_tgt) ** 2).sum(dim=(1, 2)) + 1e-20)
        norm_back = torch.sqrt(((J_tgt_to_src - J_src) ** 2).sum(dim=(1, 2)) + 1e-20)
        weight = weight * valid.to(dtype)
        w_sum = weight.sum()
        if w_sum >= 1.0:
            pair_losses.append((weight * (norm_fwd + norm_back)).sum() / (2.0 * w_sum))

    if not pair_losses:
        return next(iter(J_dict.values())).sum() * 0.0
    return torch.stack(pair_losses).mean()


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

def run_poly_unit(args, geo, base_fn, base_info, tag, table_before):
    """CS unit in polynomial free-field mode (see extra/poly_cs.py)."""
    from extra.poly_cs import run_poly_cs, run_polar_ladder, run_conjpoly_ladder
    progress = open(os.path.join(args.out, 'progress.jsonl'), 'a')
    t0 = time.time()
    elapsed = lambda: time.time() - t0

    def log_fn(rec):
        progress.write(json.dumps(rec) + '\n'); progress.flush()
    runner = {'polar': run_polar_ladder, 'conjpoly': run_conjpoly_ladder}.get(args.cs_mode, run_poly_cs)
    field, K_fn, ladder = runner(args, geo, base_fn, base_info, tag, log_fn, elapsed)
    torch.save({'model': field.state_dict(), 'degree': field.degree, 'trilin': field.trilin, 'blocks': [m + 1 for m in geo.dims], 'case': args.case, 'task': 'cs',
                'seed': args.seed, 'kind': geo.kind, 'model_type': 'poly_field', 'base_from': args.base_from},
               os.path.join(args.out, 'final.pt'))
    meta = {'exact_acs': args.cs_mode == 'conjpoly'}      # conjpoly returns J itself
    rep = full_report(geo, 'cs', K_fn, meta, args.n_dense, args.n_chart, args.adv_starts, args.adv_steps)
    result = {'case': args.case, 'task': 'cs', 'seed': args.seed, 'formulation': geo.kind,
              'cs_mode': args.cs_mode, 'degree': field.label, 'n_params': field.W.numel(),
              'model_type': 'poly_field', 'feat': f'poly{field.label}', 'ladder': ladder,
              'base': base_info, 'steps_done': len(ladder),
              'steps_planned': len(str(args.poly_degrees).split(',')), 'stopped_by_wall': False,
              'elapsed_min': elapsed() / 60.0, 'exact_acs_by_construction': meta['exact_acs'],
              'table_before': table_before, **rep}
    with open(os.path.join(args.out, 'result.json'), 'w') as f:
        json.dump(result, f, indent=2)
    corr = rep.get('corrected', {})
    if meta['exact_acs']:
        corr = rep['dense']
        rep['certified_exact_acs'] = True                     # exact by construction
    print(f"[{tag}] RESULT ({args.cs_mode} {field.label}): "
          f"{'exact ACS by construction' if meta['exact_acs'] else ('ACS CERTIFIED' if rep.get('certified_exact_acs') else 'ACS NOT certified')} | "
          f"CS {'YES' if rep['adv_nij_sup'] < args.poly_adv_tol else 'no'} (adv sup |N| vs {args.poly_adv_tol:g}) | adv |K^2+P| "
          f"{rep.get('adv_acs_sup', float('nan')):.2e} | corrected N mean {corr.get('nij_mean', float('nan')):.3e} "
          f"max {corr.get('nij_sup', float('nan')):.3e} | adv N sup {rep['adv_nij_sup']:.3e} | "
          f"chart cs {rep['chart_cs']:.3e} (table before {table_before})", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', required=True, choices=ALL_CASES)
    ap.add_argument('--formulation', default='ambient', choices=['ambient', 'frame', 'auto'],
                    help="'ambient' (default): one common ambient-projector description for ALL cases; "
                         "'frame': global-frame Nijenhuis (parallelizable cases only); 'auto': frame where available")
    ap.add_argument('--task', required=True, choices=['acs', 'cs'])
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--out', required=True)
    ap.add_argument('--steps', type=int, default=None)
    ap.add_argument('--hidden', default=None)
    ap.add_argument('--stages', type=int, default=None)
    ap.add_argument('--rar', action='store_true',
                    help='cs: residual-adaptive half batches (worst of a 2x candidate pool)')
    ap.add_argument('--tail', action='store_true',
                    help='cs: tail-aware batches = uniform half + worst of a pool + a hard set of '
                         'adversarially ascended points (refreshed every --tail-refresh steps)')
    ap.add_argument('--tail-refresh', type=int, default=500)
    ap.add_argument('--tail-hard', type=int, default=24, help='number of adversarial points per batch')
    ap.add_argument('--lr', type=float, default=1e-3, help='initial Adam lr (use ~3e-4 for warm starts)')
    ap.add_argument('--init-from', default=None,
                    help='cs: warm start from the final.pt of a previous unit of the same case/architecture')
    ap.add_argument('--lbfgs', type=int, default=200, help='total L-BFGS finisher iterations')
    ap.add_argument('--lbfgs-batch', type=int, default=768)
    ap.add_argument('--lbfgs-chunk', type=int, default=200,
                    help='iterations per L-BFGS call; a fresh batch is drawn for each chunk')
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
    ap.add_argument('--lam0', type=float, default=0.0,
                    help='acs phase 1: force the Dirichlet homotopy with this initial weight on every attempt (0 = ladder default)')
    ap.add_argument('--lam-retry', type=float, default=3.0,
                    help='acs phase 1: initial Dirichlet weight of the retry attempt after a rejected plain attempt (0 = no retry)')
    ap.add_argument('--lam1', type=float, default=1e-3, help='acs phase 1: final Dirichlet weight of the anneal')
    ap.add_argument('--attempts', default='plain,dirichlet,orient',
                    help="acs phase 1: attempt kinds per ladder class, in order (plain, dirichlet, orient)")
    ap.add_argument('--mu-orient', type=float, default=1.0,
                    help="acs phase 1: weight of the orientation penalty in 'orient' attempts")
    ap.add_argument('--adv-starts', type=int, default=256)
    ap.add_argument('--adv-steps', type=int, default=150)
    ap.add_argument('--conj', default='exp', choices=['exp', 'cayley', 'affine'],
                    help="cs: conjugation map g(B): matrix exponential (original), Cayley transform or "
                         "affine g = I + B (rational maps: ~25x cheaper forward-mode derivative)")
    ap.add_argument('--cs-scale', type=float, default=None,
                    help='cs: output scale of the generator nets at init; default 0.5 for exp (O(1) generic '
                         'start), 0.05 for cayley/affine (their g is singular for O(1) generators)')
    ap.add_argument('--fast-cs', action='store_true',
                    help='cs (ambient, with a base): training loss via the decomposed derivative '
                         '(the frozen base is differentiated without an autograd graph)')
    ap.add_argument('--cs-mode', default='conj', choices=['conj', 'poly', 'polar', 'conjpoly', 'atlas'],
                    help="cs: 'conj' = J = g J_base g^-1 with MLP generators (original); 'poly' = free "
                         "polynomial field K = P M(p) P from the base ACS, loss E|N|^2 + w E|K^2+P|^2, degree "
                         "ladder, Levenberg-Marquardt finisher (extra/poly_cs.py)")
    ap.add_argument('--dense-limit', type=int, default=16000,
                    help='cs polar: above this many coefficients use matrix-free CG-LM instead of dense LM')
    ap.add_argument('--cg-iters', type=int, default=200, help='cs polar: CG iterations per LM step')
    ap.add_argument('--polar-lm-iters', type=int, default=3, help='cs polar: LM iterations per round')
    ap.add_argument('--patience', type=int, default=5,
                    help='cs polar/conjpoly: stop a degree after this many rounds without a 10% held-out gain')
    ap.add_argument('--polar-barrier', type=float, default=1e4,
                    help='cs polar: weight of the admissibility barrier relu(|K^2+P| - r0)^2')
    ap.add_argument('--polar-r0', type=float, default=0.5, help='cs polar: barrier threshold (< 1)')
    ap.add_argument('--atlas-charts', type=int, default=40, help='cs atlas: number of charts (farthest points)')
    ap.add_argument('--atlas-overlap', type=float, default=1.3,
                    help='cs atlas: chart radius = this x covering radius (> 1 for overlaps)')
    ap.add_argument('--atlas-pool', type=int, default=40000, help='cs atlas: sample pool per round')
    ap.add_argument('--atlas-pts', type=int, default=40, help='cs atlas: points per overlap per round')
    ap.add_argument('--atlas-barrier', type=float, default=10.0, help='cs atlas: invertibility barrier weight')
    ap.add_argument('--atlas-reg', type=float, default=1e-4,
                    help='cs atlas: Tikhonov weight on nonlinear chart coefficients (raised on rollback)')
    ap.add_argument('--lp-max', type=float, default=2.0,
                    help='cs conjpoly: L^p objective by IRLS, p ramped from 2 to this (2 = mean square)')
    ap.add_argument('--lp-ramp', type=int, default=10, help='cs conjpoly: rounds over which p ramps up')
    ap.add_argument('--init-poly-from', default=None,
                    help='cs conjpoly: warm start the first degree from a finished conjpoly unit')
    ap.add_argument('--conj-reg', type=float, default=0.0,
                    help='cs conjpoly: Tikhonov weight on the generator B (keeps e^B well conditioned)')
    ap.add_argument('--polar-lm-points', type=int, default=384, help='cs polar: minimum LM collocation points')
    ap.add_argument('--polar-gauge', type=float, default=0.3,
                    help="cs polar: weight of the gauge term E|K^2+P|^2 (J* = polar(K) is exact regardless)")
    ap.add_argument('--poly-adv-tol', type=float, default=1e-6,
                    help='cs poly/polar: a class is a CS only if the adversarial sup of |N(J*)| is below this')
    ap.add_argument('--poly-degrees', default='2,3',
                    help="cs poly: class ladder, e.g. '2,3' or '2t,3' (t = + three-factor products)")
    ap.add_argument('--poly-adam', type=int, default=2000, help='cs poly: Adam warm-up steps per degree')
    ap.add_argument('--poly-lr', type=float, default=3e-3)
    ap.add_argument('--poly-w', type=float, default=10.0, help='cs poly: weight of E|K^2+P|^2')
    ap.add_argument('--poly-lm-points', type=int, default=128, help='cs poly: minimum LM collocation points')
    ap.add_argument('--poly-lm-oversample', type=float, default=2.0,
                    help='cs poly: LM points >= oversample * n_params / residual rows per point')
    ap.add_argument('--poly-lm-rounds', type=int, default=30,
                    help='cs poly: maximum LM rounds per class (fresh points each)')
    ap.add_argument('--poly-stall', type=float, default=1.5,
                    help='cs poly: a round must cut the held-out loss by this factor, else the class stalls')
    ap.add_argument('--poly-lm-iters', type=int, default=8, help='cs poly: LM iterations per round')
    ap.add_argument('--poly-tol', type=float, default=1e-20,
                    help='cs poly: stop LM and accept a class once the held-out loss is below this')
    ap.add_argument('--poly-tol-stall', type=float, default=1e-16,
                    help='cs poly: a class whose LM stalls below this held-out loss is also accepted '
                         '(|N| rms ~ 1e-8; LM converges only linearly near degenerate zero sets)')
    args = ap.parse_args()
    table_before = TABLE_BEFORE.get(args.case, {}).get(args.task)

    os.makedirs(args.out, exist_ok=True)
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    geo = Geometry(args.case, args.formulation)
    d, N = geo.d, geo.N
    steps = args.steps or (10000 if args.task == 'acs' else 40000)
    hidden = [int(h) for h in args.hidden.split(',')] if args.hidden else \
        ([256, 256, 256, 256] if d >= 6 else [128, 128, 128])
    stages = args.stages or (2 if args.task == 'cs' else 1)
    n_mat = d if geo.kind == 'frame' else N
    if args.cs_scale is None:
        args.cs_scale = 0.5 if args.conj == 'exp' else 0.05
    out_scale = 0.3 if args.task == 'acs' else args.cs_scale
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
                      'base': base_info, 'table_before': table_before}
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
    if args.task == 'cs' and args.cs_mode == 'atlas':
        if base_fn is None or geo.kind != 'ambient':
            raise ValueError('--cs-mode atlas needs the ambient formulation and --base-from (own certified ACS)')
        from extra.atlas_cs import run_atlas
        prog = open(os.path.join(args.out, 'progress.jsonl'), 'a')
        t0 = time.time()
        log = lambda rec: (prog.write(json.dumps(rec) + '\n'), prog.flush())
        atlas, ladder = run_atlas(args, geo, base_fn, tag, log, lambda: time.time() - t0)
        torch.save({'W': atlas.W.detach(), 'centres': atlas.c, 'rho': atlas.rho, 'degree': atlas.degree,
                    'case': args.case, 'base_from': args.base_from}, os.path.join(args.out, 'final.pt'))
        best = min(ladder, key=lambda e: e.get('adv_N_sup', float('inf')))
        result = {'case': args.case, 'task': 'cs', 'seed': args.seed, 'formulation': geo.kind,
                  'cs_mode': 'atlas', 'ladder': ladder, 'base': base_info, 'table_before': table_before,
                  'adv_nij_sup': best.get('adv_N_sup'), 'cs_found': (best.get('adv_N_sup') or 1e9) < args.poly_adv_tol,
                  'elapsed_min': (time.time() - t0) / 60}
        with open(os.path.join(args.out, 'result.json'), 'w') as f:
            json.dump(result, f, indent=2)
        print(f"[{tag}] RESULT (atlas): best adversarial sup |N(J*)| {best.get('adv_N_sup')} | "
              f"CS {'YES' if result['cs_found'] else 'no'}", flush=True)
        return
    if args.task == 'cs' and args.cs_mode in ('poly', 'polar', 'conjpoly'):
        if base_fn is None or geo.kind != 'ambient':
            raise ValueError('--cs-mode poly needs the ambient formulation and --base-from (own certified ACS)')
        run_poly_unit(args, geo, base_fn, base_info, tag, table_before)
        return
    J_fn, meta = make_field(geo, args.task, model, base_fn=base_fn, conj=args.conj)
    fast = None
    if args.fast_cs:
        if args.task != 'cs' or geo.kind != 'ambient' or base_fn is None:
            raise ValueError('--fast-cs needs an ambient cs unit with --base-from')
        fast = (make_conjugation(geo, model, args.conj), base_fn)
        xq = geo.sample(64)
        with torch.no_grad():
            n_ref = nij_norm(geo, J_fn, xq)
        n_fast = cs_nij_norm_fast(geo, fast[0], fast[1], xq).detach()
        dev = ((n_fast - n_ref).abs() / (n_ref.abs() + 1e-12)).max().item()
        if dev > 1e-6:
            raise RuntimeError(f'fast CS loss disagrees with the reference ({dev:.2e})')
        print(f"[{tag}] fast CS loss verified (max rel. dev. {dev:.1e})", flush=True)

    init_from_info = None
    state_path = os.path.join(args.out, 'state.pt')
    if args.init_from and not os.path.exists(state_path):
        ck = torch.load(os.path.join(args.init_from, 'final.pt'), map_location='cpu', weights_only=False)
        for k in ('case', 'task', 'hidden', 'stages', 'kind'):
            if ck.get(k) != {'case': args.case, 'task': args.task, 'hidden': hidden,
                             'stages': stages, 'kind': geo.kind}[k]:
                raise ValueError(f"--init-from mismatch on '{k}': checkpoint {ck.get(k)!r} vs unit "
                                 f"{ {'case': args.case, 'task': args.task, 'hidden': hidden, 'stages': stages, 'kind': geo.kind}[k]!r}")
        model.load_state_dict(ck['model'])
        prev = os.path.join(args.init_from, 'result.json')
        prev_res = json.load(open(prev)) if os.path.exists(prev) else {}
        init_from_info = {'dir': args.init_from, 'prev_steps': prev_res.get('steps_done'),
                          'prev_nij_mean': (prev_res.get('dense') or {}).get('nij_mean'),
                          'prev_base_from': ck.get('base_from')}
        print(f"[{tag}] warm start from {args.init_from} (prev steps {init_from_info['prev_steps']}, "
              f"prev Nij mean {init_from_info['prev_nij_mean']})", flush=True)

    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=steps, eta_min=1e-5)
    start_step, elapsed_prev, init_info, poly_done = 0, 0.0, None, False
    ladder_log = []
    if os.path.exists(state_path):
        st = torch.load(state_path, map_location='cpu', weights_only=False)
        expected_state = {'case': args.case, 'task': args.task, 'seed': args.seed,
                          'hidden': hidden, 'stages': stages, 'kind': geo.kind}
        for key, expected in expected_state.items():
            if key in st and st[key] != expected:
                raise ValueError(f"resume state mismatch on {key!r}: {st[key]!r} vs {expected!r}")
        if curriculum and st.get('feat'):
            feat_used = st['feat']
            model.poly = PolyNet(N, n_mat, out_scale, feat=feat_used, blocks=blocks).to(DTYPE)
            J_fn, meta = make_field(geo, args.task, model, base_fn=base_fn, conj=args.conj)
            opt = torch.optim.Adam(model.parameters(), lr=args.lr)
            sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=steps, eta_min=1e-5)
        model.load_state_dict(st['model']); opt.load_state_dict(st['opt'])
        sched.load_state_dict(st['sched'])
        start_step, elapsed_prev, init_info = st['step'], st['elapsed'], st.get('init')
        poly_done = st.get('poly_done', False)
        ladder_log = st.get('ladder_log', [])
        if st.get('torch_rng_state') is not None:
            torch.set_rng_state(st['torch_rng_state'])
        if st.get('numpy_rng_state') is not None:
            np.random.set_state(st['numpy_rng_state'])
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
                    'feat': feat_used, 'case': args.case, 'task': args.task,
                    'seed': args.seed, 'hidden': hidden, 'stages': stages, 'kind': geo.kind,
                    'elapsed': elapsed(), 'init': init_info, 'ladder_log': ladder_log,
                    'torch_rng_state': torch.get_rng_state(),
                    'numpy_rng_state': np.random.get_state()}, state_path)

    # Phase 1 (ACS): polynomial field only (MLP frozen), climbing the harmonic
    # ladder; optional annealed intrinsic Dirichlet term λ(t): λ0 → λ1.
    if curriculum and not poly_done and start_step == 0:
        classes = [args.feat] if args.feat != 'auto' else ladder_classes(len(geo.dims))
        best = None
        # Per class: attempt 1 = plain residual flow; attempt 2 (only if
        # rejected) = the same flow with the annealed intrinsic-Dirichlet
        # homotopy λ: lam_retry → lam1 (removes orientation walls where the
        # structure is parallel/Kähler-type, e.g. S¹×S¹, S², S²×S²).
        # Attempt kinds per class: 'plain' residual flow; 'dirichlet' = annealed
        # intrinsic-Dirichlet homotopy (lam_retry -> lam1); 'orient' = the residual
        # flow plus the consistent orientation penalty mu*relu(1 - s*o)^2 (s = the
        # majority orientation of the initial field), which vanishes on every ACS
        # of that orientation and penalises the K ~ 0 cores of orientation walls.
        kinds = [k.strip() for k in args.attempts.split(',') if k.strip()]
        attempts = [(k, {'plain': 0.0, 'dirichlet': args.lam_retry, 'orient': args.mu_orient}[k])
                    for k in kinds]
        attempts = [(k, v) for k, v in attempts if k == 'plain' or v > 0]
        if args.lam0 > 0:
            attempts = [('dirichlet', args.lam0)]
        accepted = False
        for feat in classes:
            for kind, lam0 in attempts:
                model.poly = PolyNet(N, n_mat, out_scale, feat=feat, blocks=blocks).to(DTYPE)
                opt_p = torch.optim.Adam(model.poly.parameters(), lr=5e-3)
                sched_p = torch.optim.lr_scheduler.CosineAnnealingLR(opt_p, T_max=args.poly_steps, eta_min=2e-4)
                orient_sign = None
                if kind == 'orient':
                    with torch.no_grad():
                        o0 = orientation_index(geo, J_fn, geo.sample(2000))
                    orient_sign = 1.0 if o0.mean().item() >= 0 else -1.0
                for ps in range(args.poly_steps):
                    x = geo.sample(BATCH)
                    loss, info = unit_loss(geo, args.task, J_fn, x)
                    if kind == 'dirichlet':
                        lam = lam0 * (args.lam1 / lam0) ** (ps / max(1, args.poly_steps - 1))
                        dir_e = dirichlet_energy(geo, J_fn, x)
                        info['dirichlet'], info['lam'] = dir_e.item(), lam
                        loss = loss + lam * dir_e
                    elif kind == 'orient':
                        o = orientation_index(geo, J_fn, x)
                        pen = torch.relu(1.0 - orient_sign * o).pow(2).mean()
                        info['orient_pen'] = pen.item()
                        loss = loss + lam0 * pen
                    if not bool(torch.isfinite(loss)):
                        continue
                    opt_p.zero_grad(); loss.backward(); opt_p.step(); sched_p.step()
                ev = dense_eval(geo, args.task, J_fn, 1000)
                adv, _ = adversarial_sup(geo, lambda x: acs_residual(geo, J_fn, x),
                                         max(32, args.adv_starts // 4), max(30, args.adv_steps // 3))
                entry = {'phase': 'poly', 'feat': feat, 'attempt': kind, 'lam0': lam0,
                         'orient_sign': orient_sign,
                         'n_params': sum(p.numel() for p in model.poly.parameters()),
                         'elapsed_min': elapsed() / 60.0, 'adv_acs_sup': adv, **ev}
                ladder_log.append(entry)
                progress.write(json.dumps(entry) + '\n'); progress.flush()
                print(f"[{tag}] ladder {feat:10s} {kind:9s} λ0={lam0:<4g} ({entry['n_params']} params) | held-out {ev} "
                      f"| adversarial sup {adv:.3f} -> {'ACCEPT' if adv < 1 else 'reject'}", flush=True)
                if best is None or adv < best[0]:
                    best = (adv, feat, {k: v.clone() for k, v in model.poly.state_dict().items()})
                if adv < 1.0:
                    accepted = True
                    break
            if accepted:
                break
        if best is not None and best[0] >= 1.0:          # nothing certified: keep the least bad
            model.poly = PolyNet(N, n_mat, out_scale, feat=best[1], blocks=blocks).to(DTYPE)
            model.poly.load_state_dict(best[2])
        feat_used = model.poly.feat
        J_fn, meta = make_field(geo, args.task, model, base_fn=base_fn, conj=args.conj)
        opt = torch.optim.Adam(model.parameters(), lr=args.lr)
        sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=steps, eta_min=1e-5)
        poly_done = True
        save_state(0)

    # Tail-aware sampling state (cs): a hard set of adversarially ascended
    # points, refreshed periodically; the sup is 100–1000x the mean on the
    # unfinished cases, so the mean loss alone does not press on the tail.
    hard_set = None
    cs_residual = (lambda x: nij_norm(geo, J_fn, x))
    # candidate-pool scoring needs no x-gradient: use the decomposed residual
    # when available (~4x cheaper than AD through the Newton-Schulz base)
    pool_score = (lambda x: cs_nij_norm_fast(geo, fast[0], fast[1], x).detach()) if fast is not None \
        else (lambda x: nij_norm(geo, J_fn, x).detach())

    def refresh_hard_set():
        pts = []
        for _ in range(max(1, args.tail_hard // 8)):
            _, xa = adversarial_sup(geo, cs_residual, n_starts=8, steps=20)
            pts.append(xa.unsqueeze(0))
        return torch.cat(pts, dim=0)

    for step in range(start_step, steps):
        if elapsed() > 0.85 * wall_s:
            stopped_by_wall = True
            break
        if args.task == 'cs' and args.tail:
            if hard_set is None or step % args.tail_refresh == 0:
                hard_set = refresh_hard_set()
            n_hard = len(hard_set)
            n_pool = BATCH // 2 - n_hard
            xc = geo.sample(4 * n_pool)
            r = pool_score(xc)
            idx = torch.topk(r, n_pool).indices
            x = torch.cat([geo.sample(BATCH - BATCH // 2), xc[idx], hard_set], dim=0)
        elif args.rar and args.task == 'cs':
            xc = geo.sample(2 * BATCH)
            r = pool_score(xc)
            idx = torch.topk(r, BATCH // 2).indices
            x = torch.cat([xc[idx], geo.sample(BATCH - BATCH // 2)], dim=0)
        else:
            x = geo.sample(BATCH)
        loss, info = unit_loss(geo, args.task, J_fn, x, fast=fast)
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

    # Second-order finisher: L-BFGS in chunks, a fresh large batch per chunk
    # (time-guarded; a chunk that makes the held-out loss worse is reverted).
    lbfgs_log = []
    if elapsed() < 0.92 * wall_s and args.lbfgs > 0:
        done = 0
        x_val = geo.sample(1500)
        def val_loss():
            return (nij_norm(geo, J_fn, x_val).detach() ** 2).mean().item() if args.task == 'cs' \
                else (acs_residual(geo, J_fn, x_val) ** 2).mean().item()
        best_val = val_loss()
        while done < args.lbfgs and elapsed() < 0.95 * wall_s:
            n_it = min(args.lbfgs_chunk, args.lbfgs - done)
            xf = geo.sample(args.lbfgs_batch)
            if args.task == 'cs' and args.tail and hard_set is not None:
                xf = torch.cat([xf, hard_set], dim=0)
            snapshot = {k: v.clone() for k, v in model.state_dict().items()}
            lb = torch.optim.LBFGS(model.parameters(), lr=1.0, max_iter=n_it,
                                   history_size=100, line_search_fn='strong_wolfe')

            def closure():
                lb.zero_grad()
                l, _ = unit_loss(geo, args.task, J_fn, xf, fast=fast)
                l.backward()
                return l
            try:
                lb.step(closure)
            except Exception as e:                 # keep the previous weights on failure
                print(f"[{tag}] LBFGS chunk failed ({type(e).__name__}); reverting", flush=True)
                model.load_state_dict(snapshot)
                break
            val = val_loss()
            reverted = val > best_val
            if reverted:
                model.load_state_dict(snapshot)
            else:
                best_val = val
            done += n_it
            rec = {'phase': 'lbfgs', 'iters_done': done, 'val_loss': val, 'best_val': best_val,
                   'reverted': reverted, 'elapsed_min': elapsed() / 60.0}
            lbfgs_log.append(rec)
            progress.write(json.dumps(rec) + '\n'); progress.flush()
            print(f"[{tag}] lbfgs {done:5d}/{args.lbfgs} | held-out loss {val:.3e} (best {best_val:.3e})"
                  f"{' REVERTED' if reverted else ''}", flush=True)
        save_state(steps)

    model_type = 'curriculum' if curriculum else 'mlp'
    torch.save({'model': model.state_dict(), 'case': args.case, 'task': args.task,
                'seed': args.seed, 'hidden': hidden, 'stages': stages, 'kind': geo.kind,
                'model_type': model_type, 'feat': feat_used, 'conj': args.conj,
                'base_from': args.base_from},
               os.path.join(args.out, 'final.pt'))

    rep = full_report(geo, args.task, J_fn, meta, args.n_dense, args.n_chart,
                      args.adv_starts, args.adv_steps)
    result = {'case': args.case, 'task': args.task, 'seed': args.seed,
              'formulation': geo.kind, 'hidden': hidden, 'stages': stages,
              'model_type': model_type, 'feat': feat_used, 'conj': args.conj,
              'cs_scale': args.cs_scale, 'fast_cs': bool(args.fast_cs),
              'ladder': ladder_log, 'base': base_info, 'init_from': init_from_info,
              'tail': bool(args.tail), 'lr': args.lr, 'lbfgs': lbfgs_log,
              'steps_done': min(step + 1, steps), 'steps_planned': steps,
              'stopped_by_wall': stopped_by_wall, 'elapsed_min': elapsed() / 60.0,
              'init': init_info, 'exact_acs_by_construction': meta['exact_acs'],
              'table_before': table_before, **rep}
    with open(os.path.join(args.out, 'result.json'), 'w') as f:
        json.dump(result, f, indent=2)
    verdict = ('CERTIFIED exact ACS' if rep.get('certified_exact_acs') else
               ('exact by construction' if meta['exact_acs'] else 'NOT certified (defects)'))
    print(f"[{tag}] RESULT: {verdict} | adv acs sup {rep.get('adv_acs_sup', float('nan')):.3e} "
          f"| adv Nij sup {rep['adv_nij_sup']:.3e} | dense {rep['dense']} "
          f"| chart {args.task} = {rep['chart_cs'] if args.task == 'cs' else rep['chart_acs']:.3e} "
          f"(table before {table_before if table_before is not None else 'n/a'})", flush=True)


if __name__ == '__main__':
    main()
