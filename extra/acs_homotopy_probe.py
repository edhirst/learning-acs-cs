"""Probe: does an annealed smoothness homotopy remove orientation domain walls
in the free ACS flow?  loss(t) = mean ‖K²+Id_T‖² + λ(t)·mean ‖D_x K̃‖²,
λ geometric from λ0 → λ1. Polynomial (degree ≤ 2) field. Reports final
residual mean / sampled sup / adversarial sup per (case, seed).
Usage: python extra/acs_homotopy_probe.py s2 s4 s6 [--steps 3000] [--lam0 3] [--lam1 1e-3] [--seeds 3]
"""
from __future__ import annotations
import os, sys, argparse, time, torch
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from extra.campaign_unit import (Geometry, PolyNet, make_field, acs_residual, adversarial_sup,
                                 polar_corrected, identity_T, DTYPE, BATCH)


def dirichlet(geo, K_fn, x):
    """INTRINSIC Dirichlet energy of the tangential field K along the manifold:
    ambient: Σ_a ‖P (D_{Pe_a} K) P‖²_F  (Levi-Civita covariant derivative of an
    endomorphism field on a submanifold = projected ambient derivative);
    frame:   Σ_a ‖E_a(K)‖²_F (derivatives along the global frame).
    The round S² structure [x]_× is parallel (energy 0); a domain-wall field
    ⟨x,n⟩·J_std is not — so annealing this term removes the walls."""
    def single(xs):
        return K_fn(xs.unsqueeze(0))[0]
    DK = torch.func.vmap(torch.func.jacfwd(single))(x)          # (B, n, n, N)
    if geo.kind == 'frame':
        E = geo.frame.frame_fields(x)                            # (B, d, N)
        DKa = torch.einsum('bijN,baN->baij', DK, E)
        return (DKa ** 2).sum(dim=(1, 2, 3)).mean()
    P = geo.manifold.projector(x)
    DKa = torch.einsum('bijN,bNa->baij', DK, P)                  # derivative along P e_a
    DKa = torch.einsum('bik,bakl,blj->baij', P, DKa, P)          # project to T⊗T*
    return (DKa ** 2).sum(dim=(1, 2, 3)).mean()


def dirichlet_U(geo, K_fn, x, eps=1e-2):
    """Intrinsic Dirichlet energy of the SCALE-INVARIANT orthogonal polar factor
    U = K (KᵀK + εP + (I−P))^{-1/2}. Rescaling K costs nothing; an orientation
    wall (K → −K across a surface) is a sign flip of U with energy ~ 1/ε.
    Needed because on S⁶ / S⁴×S² no ACS is parallel: the energy of K itself
    prefers ⟨x,n⟩·J_oct (60/7) over J_oct (24)."""
    def U_fn(xs):
        K = K_fn(xs)
        IdT = identity_T(geo, xs)
        eye = torch.eye(K.shape[-1], dtype=K.dtype).expand_as(K)
        S = K.transpose(1, 2) @ K + eps * IdT + (eye - IdT)
        # SPD inverse sqrt by scaled Newton–Schulz (eigh's forward-AD breaks on
        # repeated eigenvalues); spectrum of S/c lies in (0, 1].
        c = S.norm(dim=(1, 2), keepdim=True)
        Y, Z = S / c, eye
        for _ in range(20):
            T = 0.5 * (3.0 * eye - Z @ Y)
            Y, Z = Y @ T, T @ Z
        return K @ (Z / c.sqrt())
    return dirichlet(geo, U_fn, x)


def run(case, seed, steps, lam0, lam1, lr, norm=False, eps=1e-2, feat='all2'):
    torch.manual_seed(seed)
    geo = Geometry(case)
    n_mat = geo.d if geo.kind == 'frame' else geo.N
    model = PolyNet(geo.N, n_mat, 0.3, feat=feat, blocks=[m + 1 for m in geo.dims]).to(DTYPE)
    K_fn, _ = make_field(geo, 'acs', model)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=steps, eta_min=lr / 30)
    hist = []
    for t in range(steps):
        lam = lam0 * (lam1 / lam0) ** (t / max(1, steps - 1)) if lam0 > 0 else 0.0
        x = geo.sample(BATCH)
        res = (acs_residual(geo, K_fn, x) ** 2).mean()
        dir_ = (dirichlet_U(geo, K_fn, x, eps) if norm else dirichlet(geo, K_fn, x)) if lam > 0 else torch.zeros(())
        loss = res + lam * dir_
        opt.zero_grad(); loss.backward(); opt.step(); sched.step()
        if (t + 1) % (steps // 5) == 0:
            hist.append((t + 1, lam, res.item(), dir_.item()))
    xs = geo.sample(4000)
    r = acs_residual(geo, K_fn, xs)
    adv, _ = adversarial_sup(geo, lambda x: acs_residual(geo, K_fn, x), 128, 100)
    return r.mean().item(), r.max().item(), adv, hist


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('cases', nargs='+')
    ap.add_argument('--steps', type=int, default=3000)
    ap.add_argument('--lam0', type=float, default=3.0)
    ap.add_argument('--lam1', type=float, default=1e-3)
    ap.add_argument('--lr', type=float, default=5e-3)
    ap.add_argument('--seeds', type=int, default=3)
    ap.add_argument('--norm', action='store_true', help='Dirichlet energy of the polar factor U instead of K')
    ap.add_argument('--eps', type=float, default=1e-2)
    ap.add_argument('--feat', default='all2')
    a = ap.parse_args()
    for case in a.cases:
        for seed in range(a.seeds):
            t0 = time.time()
            m, s, adv, hist = run(case, seed, a.steps, a.lam0, a.lam1, a.lr, a.norm, a.eps, a.feat)
            traj = '  '.join(f"{st}:{res:.2e}/{d:.1f}" for st, lam, res, d in hist)
            print(f"{case:6s} [{a.feat}] seed {seed} {'U' if a.norm else 'K'}-λ {a.lam0}→{a.lam1}: residual mean {m:.2e} sup {s:.3f} ADV {adv:.3f} "
                  f"{'CERTIFIABLE' if adv < 1 else 'defects'} [{time.time()-t0:.0f}s]  traj(res/dir) {traj}", flush=True)
