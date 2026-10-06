"""
Free tangential polynomial endomorphism fields on a product of spheres, with
explicit derivatives, for the CS stage of the ambient pipeline.

    K(p) = P(p) M(p) P(p),      M(p) = sum_f phi_f(p) W_f,

phi = all monomials of total degree <= `degree` in the ambient coordinates
(degree 0, 1, 2, ... = the lowest spherical harmonics first; the ACS ladder's
classes are sub-classes of these). K is LINEAR in the weights W, so the
Nijenhuis tensor N_K is a QUADRATIC polynomial map of W and the CS search in a
fixed class is a small nonlinear least-squares problem: residual vector r(W),
objective |r|^2 = E|N_K|^2 + w E|K^2 + P|^2. It is zero exactly at an
integrable almost complex structure of the class.

Everything that does not depend on W (P, DP, phi, Dphi, an orthonormal tangent
frame E) is precomputed once per batch (`PolyCache`), after which K, DK and N
are plain einsums: no autograd through x, and exact parameter Jacobians by
forward-mode AD over W (used by `src/optim/lm.py`).

Residual layout (intrinsic, orthonormal): with E(p) an orthonormal basis of
T_pM, N is tangent-valued and antisymmetric, so

    sum_{i,j} |N(Pe_i, Pe_j)|^2 = sum_{a,b} |N(E_a, E_b)|^2 = 2 sum_{a<b} |E^T N(E_a, E_b)|^2,

and K^2 + P = E (E^T K^2 E + I_d) E^T. The residual vector is therefore
[sqrt(2) E^T N(E_a, E_b) (a<b) ; sqrt(w) (E^T K^2 E + I_d)] / sqrt(B):
d * d(d-1)/2 + d^2 rows per point instead of the N^3 + N^2 overcomplete ones,
with the same squared norm.
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import itertools
import math

import torch
import torch.nn as nn

from src.loss.ambient_loss import ambient_nijenhuis, compute_P_and_DP


def monomial_keys(n_in: int, degree: int, blocks=None, trilin: bool = False):
    """Monomials as sorted index tuples: all of total degree <= `degree`
    (degree 0, 1, 2, ... in combinations-with-replacement order), plus, with
    `trilin`, the three-factor products x_A,a x_B,b x_C,c of degree 3 over
    every triple of distinct factor blocks (the multidegree-(1,1,1) harmonics
    of a product with >= 3 factors) that are not already present."""
    keys = []
    for k in range(degree + 1):
        keys += list(itertools.combinations_with_replacement(range(n_in), k))
    if trilin:
        if blocks is None or len(blocks) < 3:
            raise ValueError("trilin features need >= 3 factor blocks")
        offs = [0]
        for b in blocks:
            offs.append(offs[-1] + b)
        have = set(keys)
        for A, B_, C in itertools.combinations(range(len(blocks)), 3):
            for t in itertools.product(range(offs[A], offs[A + 1]), range(offs[B_], offs[B_ + 1]),
                                       range(offs[C], offs[C + 1])):
                t = tuple(sorted(t))
                if t not in have:
                    keys.append(t); have.add(t)
    return keys


def monomial_index(keys, n_in: int) -> torch.Tensor:
    """(n_feat, width) index matrix; entry n_in is a sentinel for the factor 1."""
    width = max(1, max(len(k) for k in keys))
    return torch.tensor([list(k) + [n_in] * (width - len(k)) for k in keys], dtype=torch.long)


def parse_class(spec):
    """'2' -> (2, False); '2t' -> (2, True) (degree <= 2 plus three-factor products)."""
    spec = str(spec).strip()
    return (int(spec[:-1]), True) if spec.endswith('t') else (int(spec), False)


class PolyField(nn.Module):
    """K = P M P with M a degree-<=`degree` matrix polynomial on R^N.
    forward(x) returns M(x) with shape (B, 1, N, N), compatible with the
    `model(x)[:, 0]` convention of extra/campaign_unit.make_field."""

    def __init__(self, n_in: int, degree: int, init_scale: float = 0.3, blocks=None,
                 trilin: bool = False):
        super().__init__()
        self.n_in, self.degree, self.trilin = n_in, degree, trilin
        self.keys = monomial_keys(n_in, degree, blocks, trilin)
        self.register_buffer('idx', monomial_index(self.keys, n_in))
        self.n_feat = self.idx.shape[0]
        self.W = nn.Parameter(torch.randn(n_in * n_in, self.n_feat, dtype=torch.float64)
                              * init_scale / math.sqrt(self.n_feat))
        self.n_mat, self.stages = n_in, 1

    @property
    def label(self) -> str:
        return f"{self.degree}{'t' if self.trilin else ''}"

    def raw_features(self, x: torch.Tensor) -> torch.Tensor:
        xe = torch.cat([x, torch.ones_like(x[:, :1])], dim=1)
        return xe[:, self.idx].prod(dim=2)

    def features(self, x: torch.Tensor) -> torch.Tensor:
        F = self.raw_features(x)
        return F if getattr(self, 'T', None) is None else F @ self.T

    def orthonormalize(self, manifold, n: int = 20000, rtol: float = 1e-9) -> int:
        """Replace the monomials by an orthonormal basis of their span in
        L^2(M) (SVD over n samples). Removes the redundancy from |x_A| = 1 and
        conditions J^T J; the function class is unchanged. Resets W (call
        before initialising). Returns the number of basis functions."""
        with torch.no_grad():
            F = self.raw_features(manifold.sample(n)) / math.sqrt(n)
            _, sv, Vh = torch.linalg.svd(F, full_matrices=False)
            keep = sv > rtol * sv[0]
            self.register_buffer('T', (Vh[keep].T / sv[keep]).contiguous())
        self.n_feat = self.T.shape[1]
        self.W = nn.Parameter(torch.zeros(self.n_in * self.n_in, self.n_feat, dtype=torch.float64))
        return self.n_feat

    def load_from(self, other: 'PolyField') -> None:
        """Warm start from a smaller class: copy coefficients by monomial identity,
        zero elsewhere."""
        pos = {k: i for i, k in enumerate(self.keys)}
        with torch.no_grad():
            self.W.zero_()
            for j, k in enumerate(other.keys):
                self.W[:, pos[k]] = other.W[:, j]

    def features_and_grad(self, x: torch.Tensor):
        """phi (B, F) and its ambient gradient Dphi (B, F, N)."""
        def single(xs):
            return self.features(xs.unsqueeze(0)).squeeze(0)
        with torch.no_grad():
            phi = self.features(x)
            dphi = torch.func.vmap(torch.func.jacfwd(single))(x)
        return phi, dphi

    def matrix(self, phi: torch.Tensor, W: torch.Tensor = None) -> torch.Tensor:
        W = self.W if W is None else W
        return (phi @ W.T).reshape(-1, self.n_in, self.n_in)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.matrix(self.features(x)).unsqueeze(1)


class PolyCache:
    """Everything W-independent for a batch of points on the manifold."""

    def __init__(self, manifold, field: PolyField, x: torch.Tensor):
        self.x = x
        self.P, self.DP = compute_P_and_DP(manifold, x)
        self.P, self.DP = self.P.detach(), self.DP.detach()
        self.phi, self.dphi = field.features_and_grad(x)
        d = manifold.dim
        evals, evecs = torch.linalg.eigh(self.P)              # ascending: normal dirs first
        self.E = evecs[:, :, -d:]                              # (B, N, d) orthonormal tangent frame
        self.iu = torch.triu_indices(d, d, 1)
        self.B, self.d = x.shape[0], d


def poly_K_DK(field: PolyField, cache: PolyCache, W: torch.Tensor = None):
    """K (B,N,N) and DK (B,N,N,N) [b, a, r, c] = d_c K_ar, explicitly."""
    W = field.W if W is None else W
    n = field.n_in
    P, DP = cache.P, cache.DP
    M = (cache.phi @ W.T).reshape(-1, n, n)
    DM = torch.einsum('bfc,kf->bkc', cache.dphi, W).reshape(-1, n, n, n)
    K = P @ M @ P
    DK = torch.einsum('bikc,bkl,blj->bijc', DP, M, P) \
        + torch.einsum('bik,bklc,blj->bijc', P, DM, P) \
        + torch.einsum('bik,bkl,bljc->bijc', P, M, DP)
    return K, DK


def poly_residual(field: PolyField, cache: PolyCache, W: torch.Tensor = None,
                  w_acs: float = 10.0, per_point: bool = False):
    """Intrinsic residual vector r with |r|^2 = mean_b (|N_K|^2 + w |K^2+P|^2).
    per_point=True returns (B, rows) instead of the flat normalised vector."""
    K, DK = poly_K_DK(field, cache, W)
    Nt = ambient_nijenhuis(K, DK, cache.P, cache.DP)              # (B, N, N, N) [b, a, i, j]
    E = cache.E
    N_int = torch.einsum('bpk,bpij,bia,bjc->bkac', E, Nt, E, E)   # (B, d, d, d)
    rN = N_int[:, :, cache.iu[0], cache.iu[1]].reshape(cache.B, -1) * math.sqrt(2.0)
    C = torch.einsum('bpa,bpq,bqc->bac', E, K @ K, E) + torch.eye(cache.d, dtype=K.dtype)
    rC = C.reshape(cache.B, -1) * math.sqrt(w_acs)
    r = torch.cat([rN, rC], dim=1)
    if per_point:
        return r
    return r.reshape(-1) / math.sqrt(cache.B)


def poly_split_losses(field: PolyField, cache: PolyCache, W: torch.Tensor = None):
    """(mean |N_K|^2, mean |K^2+P|^2) on the cached batch."""
    r = poly_residual(field, cache, W, w_acs=1.0, per_point=True)
    nN = cache.d * (cache.d * (cache.d - 1) // 2)
    return (r[:, :nN] ** 2).sum(1).mean(), (r[:, nN:] ** 2).sum(1).mean()


def fit_field(field: PolyField, target_fn, manifold, n: int = 4000) -> float:
    """Least-squares projection of a (tangential) matrix field into the class:
    W = argmin sum |M(p) - target(p)|^2 over n samples. Returns the sup of the
    tangential misfit |P M P - target| on the samples."""
    x = manifold.sample(n)
    with torch.no_grad():
        T = target_fn(x)
        F = field.features(x)
        W = torch.linalg.lstsq(F, T.reshape(n, -1)).solution        # (F, N*N)
        field.W.copy_(W.T)
        P = manifold.projector(x)
        mis = (P @ field.matrix(F) @ P - T).flatten(1).norm(dim=1).max().item()
    return mis


def _residual_rows_KDK(K, DK, P, DP, E, iu, w_acs):
    """Per-point residual rows from (K, DK) — the map whose derivative G is
    factorised against the (linear) feature map in `poly_residual_and_jacobian`."""
    d = E.shape[-1]
    Nt = ambient_nijenhuis(K, DK, P, DP)
    N_int = torch.einsum('bpk,bpij,bia,bjc->bkac', E, Nt, E, E)
    rN = N_int[:, :, iu[0], iu[1]].reshape(K.shape[0], -1) * math.sqrt(2.0)
    C = torch.einsum('bpa,bpq,bqc->bac', E, K @ K, E) + torch.eye(d, dtype=K.dtype)
    return torch.cat([rN, C.reshape(K.shape[0], -1) * math.sqrt(w_acs)], dim=1)


def _factorised_jacobian(field, cache, W, rows_fn):
    """Flat residual r and exact Jacobian dr/dW for ANY per-point residual
    rows_fn(K, DK, P, DP, E) (batched), using that K, DK are linear in W:
        dK_ij/dW_(ar)f   = phi_f P_ia P_rj
        dDK_ijc/dW_(ar)f = phi_f (dP_ia/dx_c P_rj + P_ia dP_rj/dx_c) + d_c phi_f P_ia P_rj.
    G = d rows / d(K, DK) by forward AD over the N^2 + N^3 entries per point,
    independently of the number of features."""
    K, DK = poly_K_DK(field, cache, W)
    P, DP, E = cache.P, cache.DP, cache.E

    def rows_pt(Kp, DKp, Pp, DPp, Ep):
        return rows_fn(Kp.unsqueeze(0), DKp.unsqueeze(0), Pp.unsqueeze(0),
                       DPp.unsqueeze(0), Ep.unsqueeze(0)).squeeze(0)
    with torch.no_grad():
        r = rows_fn(K, DK, P, DP, E)                                             # (B, R)
        GK, GD = torch.func.vmap(torch.func.jacfwd(rows_pt, argnums=(0, 1)))(
            K.detach(), DK.detach(), P, DP, E)                                  # (B,R,N,N), (B,R,N,N,N)
        alpha = torch.einsum('bxij,bia,brj->bxar', GK, P, P) \
            + torch.einsum('bxijc,biac,brj->bxar', GD, DP, P) \
            + torch.einsum('bxijc,bia,brjc->bxar', GD, P, DP)
        beta = torch.einsum('bxijc,bia,brj->bxarc', GD, P, P)
        Jac = torch.einsum('bxar,bf->bxarf', alpha, cache.phi) \
            + torch.einsum('bxarc,bfc->bxarf', beta, cache.dphi)
    s = 1.0 / math.sqrt(cache.B)
    return r.reshape(-1) * s, Jac.reshape(r.numel(), -1) * s


def poly_residual_and_jacobian(field: PolyField, cache: PolyCache, W: torch.Tensor = None,
                               w_acs: float = 10.0):
    """Two-term residual (N_K, sqrt(w) (K^2+P)) and its exact Jacobian."""
    iu = cache.iu
    return _factorised_jacobian(field, cache, W,
                                lambda K, DK, P, DP, E: _residual_rows_KDK(K, DK, P, DP, E, iu, w_acs))


# ---------------------------------------------------------------------------
# Exact-ACS residual: N of the polar-corrected field J* = K S^{-1/2}
# ---------------------------------------------------------------------------

def _inv_sqrt_ns(S, iters):
    n = S.shape[-1]
    I = torch.eye(n, dtype=S.dtype).expand_as(S)
    Y, Z = S / 2.0, I
    for _ in range(iters):
        T = 0.5 * (3.0 * I - Z @ Y)
        Y, Z = Y @ T, T @ Z
    return Z / math.sqrt(2.0)


def polar_JDJ(K, DK, P, DP, iters: int = 30):
    """J* = K S^{-1/2}, S = -K^2 + (I - P), and its ambient derivative
    DJ*_c = d polar[K, P] (DK_c, DP_c) by forward AD through Newton-Schulz.
    J*^2 = -P exactly on T_pM whenever K has no real eigenvalue there."""
    n = K.shape[-1]
    eye = torch.eye(n, dtype=K.dtype)

    def polar(Kb, Pb):
        return Kb @ _inv_sqrt_ns(-Kb @ Kb + (eye - Pb), iters)
    J = polar(K, P)
    DJ = torch.func.vmap(lambda dk, dp: torch.func.jvp(polar, (K, P), (dk, dp))[1],
                         in_dims=(-1, -1), out_dims=-1)(DK, DP)
    return J, DJ


def _polar_rows_KDK(K, DK, P, DP, E, iu, w_gauge, iters, mu_bar=0.0, r0=0.5):
    """Rows [sqrt(2) E^T N_{J*}(E_a,E_b) (a<b); sqrt(w_gauge) (E^T K^2 E + I)].
    The second block is only a weak gauge term (polar(K) = polar(K f(K)) for
    positive f): it keeps K near an ACS and away from the degenerate set."""
    J, DJ = polar_JDJ(K, DK, P, DP, iters)
    d = E.shape[-1]
    Nt = ambient_nijenhuis(J, DJ, P, DP)
    N_int = torch.einsum('bpk,bpij,bia,bjc->bkac', E, Nt, E, E)
    rN = N_int[:, :, iu[0], iu[1]].reshape(K.shape[0], -1) * math.sqrt(2.0)
    C = torch.einsum('bpa,bpq,bqc->bac', E, K @ K, E) + torch.eye(d, dtype=K.dtype)
    out = [rN, C.reshape(K.shape[0], -1) * math.sqrt(w_gauge)]
    if mu_bar > 0:          # admissibility barrier: penalise |K^2+P| above r0 (< 1)
        rr = torch.sqrt((C ** 2).sum(dim=(1, 2)) + 1e-30)
        out.append((math.sqrt(mu_bar) * torch.relu(rr - r0)).unsqueeze(1))
        # exactness rows: J*^2 + P = 0 wherever Newton-Schulz returned the true
        # polar factor; outside its convergence region the iterate is not an ACS
        # and its "Nijenhuis tensor" must not be rewarded
        X = torch.einsum('bpa,bpq,bqc->bac', E, J @ J, E) + torch.eye(d, dtype=K.dtype)
        out.append(X.reshape(K.shape[0], -1) * math.sqrt(mu_bar))
    return torch.cat(out, dim=1)


def polar_residual(field: PolyField, cache: PolyCache, W: torch.Tensor = None,
                   w_gauge: float = 1e-2, iters: int = 30, mu_bar: float = 0.0, r0: float = 0.5):
    K, DK = poly_K_DK(field, cache, W)
    r = _polar_rows_KDK(K, DK, cache.P, cache.DP, cache.E, cache.iu, w_gauge, iters, mu_bar, r0)
    return r.reshape(-1) / math.sqrt(cache.B)


def polar_residual_and_jacobian(field: PolyField, cache: PolyCache, W: torch.Tensor = None,
                                w_gauge: float = 1e-2, iters: int = 30, mu_bar: float = 0.0,
                                r0: float = 0.5):
    iu = cache.iu
    return _factorised_jacobian(field, cache, W,
                                lambda K, DK, P, DP, E: _polar_rows_KDK(K, DK, P, DP, E, iu, w_gauge,
                                                                         iters, mu_bar, r0))


def polar_split(field: PolyField, cache: PolyCache, W: torch.Tensor = None, iters: int = 30):
    """(mean |N_{J*}|^2, mean |K^2+P|^2) on the cached batch."""
    K, DK = poly_K_DK(field, cache, W)
    r = _polar_rows_KDK(K, DK, cache.P, cache.DP, cache.E, cache.iu, 1.0, iters)
    nN = cache.d * (cache.d * (cache.d - 1) // 2)
    return (r[:, :nN] ** 2).sum(1).mean(), (r[:, nN:] ** 2).sum(1).mean()


# ---------------------------------------------------------------------------
# Matrix-free factorised Jacobian (for classes too large for dense J^T J)
# ---------------------------------------------------------------------------

class FactorisedJacobian:
    """J = dr/dW stored as per-point factors (alpha, beta) with
        J[(b,x), (a,r,f)] = alpha[b,x,a,r] phi[b,f] + sum_c beta[b,x,a,r,c] dphi[b,f,c],
    so J v and J^T u cost O(B R N^3 + B N^3 F) and memory O(B R N^3), independent
    of the number of features F. Used by the CG-based LM in src/optim/lm.py."""

    def __init__(self, field, cache, W, rows_fn):
        K, DK = poly_K_DK(field, cache, W)
        P, DP, E = cache.P, cache.DP, cache.E

        def rows_pt(Kp, DKp, Pp, DPp, Ep):
            return rows_fn(Kp.unsqueeze(0), DKp.unsqueeze(0), Pp.unsqueeze(0),
                           DPp.unsqueeze(0), Ep.unsqueeze(0)).squeeze(0)
        with torch.no_grad():
            self.r = rows_fn(K, DK, P, DP, E)
            GK, GD = torch.func.vmap(torch.func.jacfwd(rows_pt, argnums=(0, 1)))(
                K.detach(), DK.detach(), P, DP, E)
            self.alpha = torch.einsum('bxij,bia,brj->bxar', GK, P, P) \
                + torch.einsum('bxijc,biac,brj->bxar', GD, DP, P) \
                + torch.einsum('bxijc,bia,brjc->bxar', GD, P, DP)
            self.beta = torch.einsum('bxijc,bia,brj->bxarc', GD, P, P)
        self.phi, self.dphi = cache.phi, cache.dphi
        self.n = field.n_in
        self.s = 1.0 / math.sqrt(cache.B)

    def matvec(self, v):                      # v: (N*N*F,) -> (B*R,)
        V = v.view(self.n, self.n, -1)
        a = torch.einsum('bf,arf->bar', self.phi, V)
        c = torch.einsum('bfc,arf->barc', self.dphi, V)
        out = torch.einsum('bxar,bar->bx', self.alpha, a) + torch.einsum('bxarc,barc->bx', self.beta, c)
        return out.reshape(-1) * self.s

    def rmatvec(self, u):                     # u: (B*R,) -> (N*N*F,)
        U = u.view(self.alpha.shape[0], self.alpha.shape[1])
        A = torch.einsum('bxar,bx->bar', self.alpha, U)
        C = torch.einsum('bxarc,bx->barc', self.beta, U)
        out = torch.einsum('bar,bf->arf', A, self.phi) + torch.einsum('barc,bfc->arf', C, self.dphi)
        return out.reshape(-1) * self.s

    def residual(self):
        return self.r.reshape(-1) * self.s


def polar_factorised(field, cache, W=None, w_gauge: float = 1e-2, iters: int = 30,
                     mu_bar: float = 0.0, r0: float = 0.5):
    iu = cache.iu
    return FactorisedJacobian(field, cache, W,
                              lambda K, DK, P, DP, E: _polar_rows_KDK(K, DK, P, DP, E, iu, w_gauge,
                                                                       iters, mu_bar, r0))


# ---------------------------------------------------------------------------
# Conjugation ansatz with a polynomial generator: J = e^B J_b e^{-B}, B = P M_W P
# ---------------------------------------------------------------------------

def _expm_ss(B, s: int = 4, order: int = 10):
    """exp(B) by scaling and squaring of a Taylor polynomial: matmuls only, so
    forward AD and vmap are safe (no solve)."""
    n = B.shape[-1]
    I = torch.eye(n, dtype=B.dtype).expand_as(B)
    A = B / (2.0 ** s)
    E, T = I, I
    for k in range(1, order + 1):
        T = T @ A / k
        E = E + T
    for _ in range(s):
        E = E @ E
    return E


class ConjCache(PolyCache):
    """PolyCache plus the frozen base ACS J_b and its derivative DJ_b."""

    def __init__(self, manifold, field, x, base_fn):
        super().__init__(manifold, field, x)
        from src.loss.ambient_loss import compute_J_and_DJ
        with torch.no_grad():
            self.Jb, self.DJb = compute_J_and_DJ(base_fn, x)
            self.Jb, self.DJb = self.Jb.detach(), self.DJb.detach()


def conj_JDJ(B, DB, Jb, DJb, s: int = 4, order: int = 10):
    """J = g Jb g^-1 with g = exp(B), g^-1 = exp(-B), and DJ by forward AD."""
    def J_of(Bm, Jm):
        return _expm_ss(Bm, s, order) @ Jm @ _expm_ss(-Bm, s, order)
    J = J_of(B, Jb)
    DJ = torch.func.vmap(lambda db, dj: torch.func.jvp(J_of, (B, Jb), (db, dj))[1],
                         in_dims=(-1, -1), out_dims=-1)(DB, DJb)
    return J, DJ


def _conj_rows(B, DB, P, DP, E, Jb, DJb, iu, w_reg):
    J, DJ = conj_JDJ(B, DB, Jb, DJb)
    Nt = ambient_nijenhuis(J, DJ, P, DP)
    N_int = torch.einsum('bpk,bpij,bia,bjc->bkac', E, Nt, E, E)
    rN = N_int[:, :, iu[0], iu[1]].reshape(B.shape[0], -1) * math.sqrt(2.0)
    if w_reg > 0:      # mild Tikhonov on the generator keeps g well conditioned
        return torch.cat([rN, B.reshape(B.shape[0], -1) * math.sqrt(w_reg)], dim=1)
    return rN


def _point_weights(cache, r):
    w = getattr(cache, 'w', None)
    return r if w is None else r * w.sqrt().unsqueeze(1)


def conj_residual(field, cache, W=None, w_reg: float = 0.0):
    """Flat residual; per-point weights `cache.w` (IRLS for an L^p objective) if set."""
    B, DB = poly_K_DK(field, cache, W)
    r = _conj_rows(B, DB, cache.P, cache.DP, cache.E, cache.Jb, cache.DJb, cache.iu, w_reg)
    return _point_weights(cache, r).reshape(-1) / math.sqrt(cache.B)


def conj_residual_and_jacobian(field, cache, W=None, w_reg: float = 0.0):
    B, DB = poly_K_DK(field, cache, W)
    P, DP, E, iu = cache.P, cache.DP, cache.E, cache.iu

    def rows_pt(Bp, DBp, Pp, DPp, Ep, Jp, DJp):
        u = lambda t: t.unsqueeze(0)
        return _conj_rows(u(Bp), u(DBp), u(Pp), u(DPp), u(Ep), u(Jp), u(DJp), iu, w_reg).squeeze(0)
    with torch.no_grad():
        r = _point_weights(cache, _conj_rows(B, DB, P, DP, E, cache.Jb, cache.DJb, iu, w_reg))
        GK, GD = torch.func.vmap(torch.func.jacfwd(rows_pt, argnums=(0, 1)))(
            B.detach(), DB.detach(), P, DP, E, cache.Jb, cache.DJb)
        if getattr(cache, 'w', None) is not None:
            sw = cache.w.sqrt()
            GK, GD = GK * sw.view(-1, 1, 1, 1), GD * sw.view(-1, 1, 1, 1, 1)
        alpha = torch.einsum('bxij,bia,brj->bxar', GK, P, P) \
            + torch.einsum('bxijc,biac,brj->bxar', GD, DP, P) \
            + torch.einsum('bxijc,bia,brjc->bxar', GD, P, DP)
        beta = torch.einsum('bxijc,bia,brj->bxarc', GD, P, P)
        Jac = torch.einsum('bxar,bf->bxarf', alpha, cache.phi) \
            + torch.einsum('bxarc,bfc->bxarf', beta, cache.dphi)
    s = 1.0 / math.sqrt(cache.B)
    return r.reshape(-1) * s, Jac.reshape(r.numel(), -1) * s
