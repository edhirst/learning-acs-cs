"""Integrable-by-construction atlas search for a complex structure (`--cs-mode atlas`).

Unsupervised and unbiased: the manifold is covered by geodesic caps around
centres chosen by farthest-point sampling; on each cap a polynomial map
phi_i : R^N -> R^d = C^{d/2} of the local ambient offset u = (p - c_i)/rho is a
candidate holomorphic coordinate system. The pulled-back structure

    J_i(p) = P Dphi_i^T (Dphi_i P Dphi_i^T)^{-1} J_std Dphi_i P

is an exact ACS on T_pM (J_i^2 = -P) and INTEGRABLE wherever Dphi_i|_T is
invertible (pullback of the standard structure by a local diffeomorphism), so
N(J_i) = 0 identically. The only equation left is that the charts define the
same structure on overlaps:  r_ij(p) = w_i(p) w_j(p) (J_i(p) - J_j(p)),  plus a
barrier keeping Dphi_i|_T uniformly invertible. A zero of this system is a
complex structure (the transition maps are then holomorphic).

Why this localises the problem: in a small cap the holomorphic coordinates of
any complex structure are real-analytic, so degree-k polynomial charts
approximate them with an error decaying geometrically in k. Global complexity,
which defeats global polynomial fields in round ambient coordinates on S^6,
is carried only by the gluing.

Initialisation from the unit's own certified ACS J_b: phi_i linear, mapping a
J_b(c_i)-complex frame at the centre to the standard one, so J_i(c_i) = J_b(c_i).
Training uses first derivatives only. Certification: the partition-of-unity
blend K = sum_i w_i J_i / sum_i w_i, polar-corrected to an exact ACS, with the
true Nijenhuis tensor measured densely and by adversarial ascent.
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import itertools
import math
import time

import torch

DTYPE = torch.float64


def j_std(d):
    J = torch.zeros(d, d, dtype=DTYPE)
    for k in range(d // 2):
        J[2 * k, 2 * k + 1], J[2 * k + 1, 2 * k] = -1.0, 1.0
    return J


def farthest_points(manifold, m, n_pool=20000, seed=0):
    g = torch.Generator().manual_seed(seed)
    torch.manual_seed(seed)
    pool = manifold.sample(n_pool)
    idx = [int(torch.randint(n_pool, (1,), generator=g))]
    dmin = (pool - pool[idx[0]]).norm(dim=1)
    for _ in range(m - 1):
        i = int(dmin.argmax())
        idx.append(i)
        dmin = torch.minimum(dmin, (pool - pool[i]).norm(dim=1))
    return pool[idx], dmin.max().item()          # centres, covering radius (chordal)


class Atlas(torch.nn.Module):
    """Polynomial chart maps phi_i(p) = W_i feat((p - c_i)/rho), one per centre."""

    def __init__(self, manifold, centres, rho, degree):
        super().__init__()
        self.man, self.rho, self.degree = manifold, rho, degree
        self.N, self.d = manifold.ambient_dim, manifold.dim
        self.register_buffer('c', centres.clone())
        keys = []
        for k in range(1, degree + 1):                       # constants are irrelevant (Dphi only)
            keys += list(itertools.combinations_with_replacement(range(self.N), k))
        width = max(len(k) for k in keys)
        self.register_buffer('idx', torch.tensor([list(k) + [self.N] * (width - len(k)) for k in keys]))
        self.n_feat = len(keys)
        self.n_lin = self.N                                   # first N features are the linear ones
        self.W = torch.nn.Parameter(torch.zeros(len(centres), self.d, self.n_feat, dtype=DTYPE))
        self.Jstd = j_std(self.d)

    def feat(self, u):                                        # (B, N) -> (B, F)
        ue = torch.cat([u, torch.ones_like(u[:, :1])], dim=1)
        return ue[:, self.idx].prod(dim=2)

    def dfeat(self, u):                                       # (B, F, N)
        def single(us):
            return self.feat(us.unsqueeze(0)).squeeze(0)
        return torch.func.vmap(torch.func.jacfwd(single))(u)

    def weight(self, p, i):
        """Smooth bump of chart i: (1 - |p - c_i|^2 / rho^2)_+^3."""
        s = ((p - self.c[i]) ** 2).sum(-1) / self.rho ** 2
        return torch.relu(1.0 - s) ** 3

    def Dphi(self, p, i, W=None):
        """Ambient Jacobian of chart i at points p: (B, d, N)."""
        W = self.W if W is None else W
        u = (p - self.c[i]) / self.rho
        return torch.einsum('df,bfn->bdn', W[i], self.dfeat(u)) / self.rho

    def J_chart(self, p, i, W=None, P=None):
        P = self.man.projector(p) if P is None else P
        D = self.Dphi(p, i, W) @ P                            # (B, d, N), rows span T*
        G = D @ D.transpose(1, 2)
        Ginv = torch.linalg.inv(G)
        return D.transpose(1, 2) @ Ginv @ self.Jstd @ D, D

    def init_from(self, base_fn):
        """Linear charts with J_i(c_i) = J_b(c_i): rows of A_i = a J_b-adapted
        complex coframe at the centre, extended by the tangent projector."""
        P = self.man.projector(self.c)
        Jb = base_fn(self.c)
        evals, evecs = torch.linalg.eigh(P)
        E = evecs[:, :, -self.d:]                             # (M, N, d)
        with torch.no_grad():
            self.W.zero_()
            for i in range(len(self.c)):
                Jf = E[i].T @ Jb[i] @ E[i]                   # d x d ACS in the frame
                cols, basis = [], torch.eye(self.d, dtype=DTYPE)
                for v in basis:                              # greedy complex basis v, Jv
                    M = torch.stack(cols, 1) if cols else torch.zeros(self.d, 0, dtype=DTYPE)
                    cand = torch.cat([M, v.unsqueeze(1), (Jf @ v).unsqueeze(1)], 1)
                    if torch.linalg.matrix_rank(cand) == cand.shape[1]:
                        cols += [v, Jf @ v]
                    if len(cols) == self.d:
                        break
                V = torch.stack(cols, 1)                     # V^{-1} Jf V = J_std
                A = torch.linalg.inv(V) @ E[i].T             # d x N
                lin = [k for k in range(self.n_feat) if self.idx[k, 1:].eq(self.N).all()]
                for k in lin:
                    self.W[i, :, k] = A[:, int(self.idx[k, 0])] * self.rho


class OverlapBatch:
    """All overlap collocation points of a round in one batch, with the
    W-independent pieces precomputed: chart indices (ci, cj), P, feature
    Jacobians at both charts, and the bump weight sqrt(w_i w_j)."""

    def __init__(self, atlas, n_pool, max_pts):
        pool = atlas.man.sample(n_pool)
        inside = torch.cdist(pool, atlas.c) < atlas.rho
        cc = torch.cdist(atlas.c, atlas.c)
        xs, ci, cj = [], [], []
        M = len(atlas.c)
        for i in range(M):
            for j in range(i + 1, M):
                if cc[i, j] >= 2 * atlas.rho:
                    continue
                m = inside[:, i] & inside[:, j]
                if m.sum() >= 4:
                    x = pool[m][:max_pts]
                    xs.append(x); ci += [i] * len(x); cj += [j] * len(x)
        self.n_pairs = len(xs)
        with torch.no_grad():
            self.x = torch.cat(xs)
            self.ci, self.cj = torch.tensor(ci), torch.tensor(cj)
            self.P = atlas.man.projector(self.x)
            self.dfi = atlas.dfeat((self.x - atlas.c[self.ci]) / atlas.rho) / atlas.rho
            self.dfj = atlas.dfeat((self.x - atlas.c[self.cj]) / atlas.rho) / atlas.rho
            si = ((self.x - atlas.c[self.ci]) ** 2).sum(-1) / atlas.rho ** 2
            sj = ((self.x - atlas.c[self.cj]) ** 2).sum(-1) / atlas.rho ** 2
            self.w = (torch.relu(1 - si) ** 3 * torch.relu(1 - sj) ** 3).sqrt()

    def __len__(self):
        return len(self.x)


def pair_samples(atlas, n_pool, max_pts):
    return OverlapBatch(atlas, n_pool, max_pts)


def _J_batch(W, idx, df, P, Jstd):
    D = torch.einsum('bdf,bfn->bdn', W[idx], df) @ P
    G = D @ D.transpose(1, 2)
    return D.transpose(1, 2) @ torch.linalg.inv(G) @ Jstd @ D, D


def gluing_residual(atlas, batch, W, tau, mu, reg=0.0):
    """Flat residual: overlap mismatch rows + invertibility barrier rows +
    Tikhonov rows on the nonlinear chart coefficients (keeps each chart a
    near-linear diffeomorphism of its cap, so Dphi cannot go singular between
    collocation points)."""
    n = max(1, len(batch))
    Ji, Di = _J_batch(W, batch.ci, batch.dfi, batch.P, atlas.Jstd)
    Jj, Dj = _J_batch(W, batch.cj, batch.dfj, batch.P, atlas.Jstd)
    out = [((Ji - Jj) * batch.w.view(-1, 1, 1)).reshape(-1)]
    if mu > 0:              # invertibility barrier on the Gram determinant of Dphi|_T
        d = atlas.d
        for D in (Di, Dj):
            g = torch.linalg.det(D @ D.transpose(1, 2)).clamp_min(0) ** (0.5 / d)
            out.append(math.sqrt(mu) * torch.relu(tau - g))
    if reg > 0:
        out.append(math.sqrt(reg * n) * W[:, :, atlas.n_lin:].reshape(-1))
    return torch.cat(out) / math.sqrt(n)


def lm_cg(fn, theta0, n_iter, lam0=1e-3, cg_iters=150, log=None):
    """Matrix-free LM on |fn(theta)|^2 with autograd jvp / vjp and CG."""
    theta = theta0.detach().clone()
    r = fn(theta).detach()
    L, lam = (r ** 2).sum().item(), lam0
    for it in range(n_iter):
        t0 = time.time()
        r, vjp_fn = torch.func.vjp(fn, theta)
        r = r.detach()
        g = vjp_fn(r)[0]
        Jv = lambda v: torch.func.jvp(fn, (theta,), (v,))[1]
        JTJ = lambda v: vjp_fn(Jv(v))[0]
        z = torch.randn(3, *theta.shape, dtype=DTYPE)
        dscale = sum((Jv(zz) ** 2).sum().item() for zz in z) / (3 * theta.numel())
        accepted = False
        while lam < 1e10:
            mu = lam * dscale
            A = lambda v: JTJ(v) + mu * v
            x = torch.zeros_like(theta); res = -g; p = res.clone(); rs = (res * res).sum()
            g0 = rs.sqrt()
            for _ in range(cg_iters):
                Ap = A(p)
                al = rs / (p * Ap).sum()
                x = x + al * p; res = res - al * Ap
                rs_new = (res * res).sum()
                if rs_new.sqrt() < 1e-8 * g0:
                    break
                p = res + (rs_new / rs) * p; rs = rs_new
            L_new = (fn(theta + x).detach() ** 2).sum().item()
            if math.isfinite(L_new) and L_new < L:
                theta, L, lam, accepted = theta + x, L_new, max(lam / 3, 1e-12), True
                break
            lam *= 4
        if log is not None:
            log(it, L, lam, accepted, time.time() - t0)
        if not accepted:
            break
    return theta, L, lam


def blended_field(atlas):
    """Global candidate: K = sum_i w_i J_i / sum_i w_i (exact where one chart
    dominates); certification polar-corrects it."""
    def K(p):
        P = atlas.man.projector(p)
        num = torch.zeros(p.shape[0], atlas.N, atlas.N, dtype=DTYPE)
        den = torch.zeros(p.shape[0], dtype=DTYPE)
        for i in range(len(atlas.c)):
            w = atlas.weight(p, i)
            if bool((w > 0).any()):
                Ji, _ = atlas.J_chart(p, i, None, P)
                num = num + w.view(-1, 1, 1) * Ji
                den = den + w
        return num / den.clamp_min(1e-30).view(-1, 1, 1)
    return K


def run_atlas(args, geo, base_fn, tag, log_fn, elapsed):
    from extra.campaign_unit import adversarial_sup, nij_norm, polar_corrected, acs_residual
    man = geo.manifold
    centres, cover = farthest_points(man, args.atlas_charts, seed=args.seed)
    rho = args.atlas_overlap * cover
    ladder, prev = [], None
    for deg in [int(k) for k in str(args.poly_degrees).split(',')]:
        atlas = Atlas(man, centres, rho, deg)
        if prev is None:
            atlas.init_from(base_fn)
        else:
            with torch.no_grad():
                atlas.W[:, :, :prev.n_feat] = prev.W            # monomial order nests by degree
        D0 = atlas.Dphi(centres[:1], 0) @ man.projector(centres[:1])
        tau = 0.3 * torch.linalg.det(D0 @ D0.transpose(1, 2)).item() ** (0.5 / atlas.d)
        entry = {'phase': 'atlas', 'degree': deg, 'charts': len(centres), 'rho': rho,
                 'cover_radius': cover, 'n_params': atlas.W.numel(), 'rounds': []}
        print(f"[{tag}] atlas deg {deg}: {len(centres)} charts, rho {rho:.3f} (cover {cover:.3f}), "
              f"{atlas.W.numel()} params", flush=True)
        lam, best, since, reg = 1e-2, float('inf'), 0, args.atlas_reg
        for rnd in range(args.poly_lm_rounds):
            pairs = pair_samples(atlas, args.atlas_pool, args.atlas_pts)
            fn = lambda th: gluing_residual(atlas, pairs, th, tau, args.atlas_barrier, reg)
            W_prev = atlas.W.detach().clone()
            th, L, lam = lm_cg(fn, atlas.W.detach(), args.polar_lm_iters, lam0=lam, cg_iters=args.cg_iters)
            with torch.no_grad():
                atlas.W.copy_(th)
            hold = pair_samples(atlas, args.atlas_pool // 2, args.atlas_pts // 2)
            with torch.no_grad():
                Lh = (gluing_residual(atlas, hold, atlas.W, tau, 0.0) ** 2).sum().item()
                Lf = (gluing_residual(atlas, pairs, atlas.W, tau, 0.0) ** 2).sum().item()
            rolled = (not math.isfinite(Lh)) or Lh > 3.0 * Lf + 1e-30
            if rolled:                    # a singular hypersurface between collocation points
                with torch.no_grad():
                    atlas.W.copy_(W_prev)
                lam, reg = lam * 30.0, reg * 3.0 + 1e-6
            rec = {'round': rnd + 1, 'fit': L, 'held_out': Lh, 'fit_mismatch': Lf, 'rolled_back': rolled,
                   'reg': reg, 'pairs': pairs.n_pairs, 'elapsed_min': elapsed() / 60}
            entry['rounds'].append(rec); log_fn(dict(phase='atlas_lm', degree=deg, **rec))
            print(f"[{tag}] atlas deg {deg} round {rnd+1} ({pairs.n_pairs} overlaps) | fit {Lf:.3e} | "
                  f"held-out mismatch {Lh:.3e}{' ROLLED BACK' if rolled else ''} | reg {reg:.1e} | "
                  f"{elapsed()/60:.1f} min", flush=True)
            if rolled:
                continue
            if Lh < best / 1.1:
                best, since = Lh, 0
            else:
                since += 1
            if since >= args.patience or Lh < 1e-24:
                break
        K = blended_field(atlas)
        x = man.sample(1000)
        with torch.no_grad():
            r = acs_residual(geo, K, x)
        entry['blend_acs_sup'] = r.max().item()
        if r.max().item() < 1.0:
            Js = polar_corrected(geo, K)
            nn_ = torch.cat([nij_norm(geo, Js, x[i:i + 100]).detach() for i in range(0, 1000, 100)])
            advN, _ = adversarial_sup(geo, lambda y: nij_norm(geo, Js, y), 32, 40)
            entry.update({'N_mean': nn_.mean().item(), 'N_sup_dense': nn_.max().item(), 'adv_N_sup': advN})
            print(f"[{tag}] atlas deg {deg} RESULT: blend |K^2+P| sup {r.max():.2e} | N(J*) mean "
                  f"{nn_.mean():.3e} dense sup {nn_.max():.3e} adversarial sup {advN:.3e}", flush=True)
        else:
            print(f"[{tag}] atlas deg {deg}: blended field not admissible (|K^2+P| sup {r.max():.2f})", flush=True)
        ladder.append(entry); log_fn(entry)
        prev = atlas
        if entry.get('adv_N_sup', 1e9) < args.poly_adv_tol or elapsed() > 0.85 * args.wall_minutes * 60:
            break
    return prev, ladder
