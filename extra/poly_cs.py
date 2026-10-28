"""CS stage, polynomial free-field mode (`campaign_unit.py --task cs --cs-mode poly`).

Unsupervised end to end, like the conjugation mode: the only inputs are the
case and the unit's OWN certified ACS (--base-from); the only training signals
are the residuals of J^2 = -Id and N_J = 0 at random points.

  field   K = P M_W(p) P,  M_W = degree-<=k matrix polynomial on R^N
          (src/models/poly_field.py). K is linear in W, N_K quadratic.
  init    least-squares projection of the certified base ACS J* = polar(K_acs)
          into the class (it must itself be admissible: adversarial sup < 1).
  loss    E|N_K|^2 + w E|K^2 + P|^2  (zero iff K is an integrable ACS).
  ladder  degree k = 2, 3, ...: Adam warm-up on the fast explicit residual, then
          dense Levenberg-Marquardt rounds on freshly resampled point sets with
          the factorised exact Jacobian; accept the first degree whose
          held-out loss is below --poly-tol, otherwise climb (warm start by
          zero-padding the coefficients).
  output  the free field K; full_report certifies it (adversarial sup of
          |K^2+P| < 1), polar-corrects it to an exact ACS J* and reports N of J*
          (dense, adversarial) and the chart protocol.

Why: the conjugation mode J = g J_base g^-1 must represent the integrable
structure through a complicated g (J_base is an MLP-refined learned field),
while the lowest-degree integrable structures live in small polynomial classes
(an offline check shows the Hopf and Calabi-Eckmann structures are exactly
degree-2 ambient fields). In a fixed class the problem is a small zero-residual
least-squares problem, where LM converges quadratically.
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import math
import os

import torch

from src.models.poly_field import (PolyField, PolyCache, poly_residual, parse_class,
                                   poly_residual_and_jacobian, poly_split_losses, fit_field)
from src.optim.lm import lm_minimise


def _held_out(field, manifold, n=512, W=None, seed_points=None):
    x = manifold.sample(n) if seed_points is None else seed_points
    c = PolyCache(manifold, field, x)
    with torch.no_grad():
        n2, c2 = poly_split_losses(field, c, W)
    return n2.item(), c2.item()


def run_poly_cs(args, geo, base_fn, base_info, tag, log_fn, elapsed):
    """Returns (field, K_fn, ladder_log). `log_fn(dict)` appends to progress.jsonl."""
    from extra.campaign_unit import adversarial_sup, acs_residual, BATCH   # lazy: avoid cycle
    man = geo.manifold
    N = geo.N
    classes = [parse_class(k) for k in str(args.poly_degrees).split(',')]
    blocks = [m + 1 for m in geo.dims]
    ladder = []
    prev = None
    wall_s = args.wall_minutes * 60.0

    def K_of(fld):
        return lambda x: man.projector(x) @ fld(x)[:, 0] @ man.projector(x)

    best = None
    for deg, trilin in classes:
        fld = PolyField(N, deg, blocks=blocks, trilin=trilin)
        if prev is None:
            fit_sup = fit_field(fld, base_fn, man, n=max(4000, 3 * fld.n_feat))
            init = 'lsq(base)'
        else:
            fld.load_from(prev)
            fit_sup, init = None, f'pad({prev.label})'
        K_fn = K_of(fld)
        adv0, _ = adversarial_sup(geo, lambda x: acs_residual(geo, K_fn, x), 64, 60)
        n2, c2 = _held_out(fld, man)
        entry = {'phase': 'poly_cs', 'degree': fld.label, 'n_params': fld.W.numel(), 'init': init,
                 'init_fit_sup': fit_sup, 'init_adv_acs_sup': adv0, 'init_N2': n2, 'init_C2': c2}
        print(f"[{tag}] poly class {fld.label}: {fld.W.numel()} params, init {init} "
              f"(fit sup {fit_sup if fit_sup is None else f'{fit_sup:.2e}'}), adversarial |K^2+P| {adv0:.3f}, "
              f"held-out N2 {n2:.3e} C2 {c2:.3e}", flush=True)

        # --- Adam warm-up on the fast explicit residual (first rung only: a
        # warm-started rung already sits at the previous optimum, and Adam's
        # noise floor is far above what LM reaches) ------------------------
        adam_steps = args.poly_adam if prev is None else 0
        opt = torch.optim.Adam(fld.parameters(), lr=args.poly_lr)
        sch = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=max(1, adam_steps),
                                                         eta_min=args.poly_lr / 30)
        for s in range(adam_steps):
            c = PolyCache(man, fld, man.sample(BATCH))
            loss = (poly_residual(fld, c, w_acs=args.poly_w) ** 2).sum()
            if not bool(torch.isfinite(loss)):
                continue
            opt.zero_grad(); loss.backward()
            torch.nn.utils.clip_grad_norm_(fld.parameters(), 1.0)
            opt.step(); sch.step()
            if (s + 1) % max(1, adam_steps // 4) == 0:
                n2, c2 = _held_out(fld, man)
                rec = {'phase': 'poly_adam', 'degree': fld.label, 'step': s + 1, 'N2': n2, 'C2': c2,
                       'elapsed_min': elapsed() / 60}
                log_fn(rec)
                print(f"[{tag}] poly {fld.label} adam {s+1:5d} | held-out N2 {n2:.3e} C2 {c2:.3e}", flush=True)

        # --- LM rounds on freshly resampled point sets --------------------
        rows = geo.d * (geo.d * (geo.d - 1) // 2) + geo.d ** 2
        n_pts = max(args.poly_lm_points, int(args.poly_lm_oversample * fld.W.numel() / rows))
        n_pts = 32 * ((n_pts + 31) // 32)
        lam = 1e-2
        lm_log = []
        best_round = float('inf')
        for rnd in range(args.poly_lm_rounds):
            if elapsed() > 0.8 * wall_s:
                break
            caches = [PolyCache(man, fld, man.sample(32)) for _ in range(n_pts // 32)]
            res = lambda W, c: poly_residual(fld, c, W, w_acs=args.poly_w)
            jac = lambda W, c: poly_residual_and_jacobian(fld, c, W, w_acs=args.poly_w)
            W, lg = lm_minimise(res, fld.W.detach(), caches, n_iter=args.poly_lm_iters,
                                lam0=lam, jac_fn=jac)
            with torch.no_grad():
                fld.W.copy_(W)
            # carry the damping to the next round, capped: a round that ended saturated
            # (lam >= 1e12 takes no step in lm_minimise) must not freeze every later round
            lam = min(max(lg[-1]['lam'], 1e-8), 1e6)
            n2, c2 = _held_out(fld, man)
            rec = {'phase': 'poly_lm', 'degree': fld.label, 'round': rnd + 1, 'n_points': n_pts,
                   'train_loss': lg[-1]['loss'], 'lm_iters': len(lg) - 1, 'N2': n2, 'C2': c2,
                   'elapsed_min': elapsed() / 60}
            lm_log.append(rec); log_fn(rec)
            print(f"[{tag}] poly {fld.label} LM round {rnd+1} ({n_pts} pts, {len(lg)-1} it) | train "
                  f"{lg[-1]['loss']:.3e} | held-out N2 {n2:.3e} C2 {c2:.3e} | {elapsed()/60:.1f} min", flush=True)
            score_r = n2 + args.poly_w * c2
            if score_r < args.poly_tol:
                break
            # keep going while the held-out loss still falls by >= poly_stall per
            # round (LM converges linearly on these degenerate zero sets);
            # climb once it stalls
            if rnd >= 2 and score_r > best_round / args.poly_stall:
                print(f"[{tag}] poly {fld.label}: LM stalled (held-out {score_r:.2e}, best {best_round:.2e})",
                      flush=True)
                break
            best_round = min(best_round, score_r)
        n2, c2 = _held_out(fld, man, n=2000)
        entry.update({'final_N2': n2, 'final_C2': c2, 'lm': lm_log, 'elapsed_min': elapsed() / 60})
        ladder.append(entry); log_fn(entry)
        score = n2 + args.poly_w * c2
        if best is None or score < best[0]:
            best = (score, fld)
        prev = fld
        if score < args.poly_tol or score < args.poly_tol_stall:
            # acceptance rule: the polar-corrected exact ACS must pass the
            # ADVERSARIAL sup test on |N|, not just the sampled mean
            from extra.campaign_unit import polar_corrected, nij_norm
            Js = polar_corrected(geo, K_of(fld))
            advN, _ = adversarial_sup(geo, lambda y: nij_norm(geo, Js, y), 64, 60)
            entry['adv_N_sup'] = advN
            if advN < args.poly_adv_tol:
                print(f"[{tag}] poly {fld.label}: ACCEPT (held-out loss {score:.2e}, adversarial sup "
                      f"|N(J*)| {advN:.2e} < {args.poly_adv_tol:g})", flush=True)
                entry['accepted'] = True
                best = (score, fld)      # return the ACCEPTED rung, even if an earlier
                break                    # rung had a lower score but failed the sup test
            print(f"[{tag}] poly {fld.label}: held-out loss {score:.2e} but adversarial sup |N(J*)| "
                  f"{advN:.2e} >= {args.poly_adv_tol:g} -> climb", flush=True)
            prev = fld
            continue
        print(f"[{tag}] poly {fld.label}: held-out loss {score:.2e} >= {args.poly_tol_stall:g} -> climb", flush=True)
        if elapsed() > 0.8 * wall_s:
            break
    field = best[1]
    return field, K_of(field), ladder


# ---------------------------------------------------------------------------
# Degree-ladder MEASUREMENT on an exact ACS (`--cs-mode polar`)
# ---------------------------------------------------------------------------

def run_polar_ladder(args, geo, base_fn, base_info, tag, log_fn, elapsed):
    """Per degree k, the floor of E|N(J*)|^2 over the class, J* = polar(K),
    K = P M_W P with M_W in an orthonormalised degree-<=k basis, subject to
    admissibility (adversarial sup |K^2+P| < 1, so J* is an exact ACS).

    Unlike the two-term loss, the field is an exact ACS at every iterate, so the
    recorded floor measures how far the class is from containing a CS, with no
    trade-off against ACS error. Init: LSQ projection of the certified base J*
    (first degree) or of the previous degree's M (exact: the spans nest).
    Solver: dense LM rounds on fresh collocation sets with the factorised exact
    Jacobian; a round that makes K inadmissible is rolled back and the gauge
    weight doubled. A degree ends when LM stalls; then the ladder climbs."""
    from extra.campaign_unit import adversarial_sup, acs_residual, nij_norm, polar_corrected
    from src.models.poly_field import (PolyField, PolyCache, polar_residual, polar_factorised,
                                       polar_residual_and_jacobian, polar_split)
    from src.optim.lm import lm_cg_minimise
    man, N = geo.manifold, geo.N
    degrees = [int(k) for k in str(args.poly_degrees).split(',')]
    rows = geo.d * (geo.d * (geo.d - 1) // 2) + geo.d ** 2
    ladder, prev = [], None
    wall_s = args.wall_minutes * 60.0

    def K_of(fld):
        return lambda x: man.projector(x) @ fld(x)[:, 0] @ man.projector(x)

    def held_out(fld, n=1000):
        c = PolyCache(man, fld, man.sample(n))
        with torch.no_grad():
            n2, c2 = polar_split(fld, c)
        return n2.item(), c2.item()

    def adv_acs(fld):
        Kf = K_of(fld)
        return adversarial_sup(geo, lambda x: acs_residual(geo, Kf, x), 64, 60)[0]

    for deg in degrees:
        fld = PolyField(N, deg)
        nf = fld.orthonormalize(man)
        src = base_fn if prev is None else K_of(prev)
        fit_sup = fit_field(fld, src, man, n=max(4000, 4 * nf))
        a0 = adv_acs(fld)
        n2, c2 = held_out(fld)
        entry = {'phase': 'polar_ladder', 'degree': deg, 'n_feat': nf, 'n_params': fld.W.numel(),
                 'init': 'lsq(base)' if prev is None else f'lsq(deg{prev.degree})',
                 'init_fit_sup': fit_sup, 'init_adv_acs_sup': a0, 'init_N2': n2, 'lm': []}
        print(f"[{tag}] polar deg {deg}: {nf} basis fns, {fld.W.numel()} params, init fit {fit_sup:.1e}, "
              f"adv |K^2+P| {a0:.3f}, held-out |N(J*)|^2 {n2:.4e}", flush=True)
        if a0 >= 1.0:
            entry['status'] = 'inadmissible init'
            ladder.append(entry); log_fn(entry)
            break
        n_pts = max(args.polar_lm_points, int(args.poly_lm_oversample * fld.W.numel() / rows))
        n_pts = 32 * ((n_pts + 31) // 32)
        lam, wg, best, since = 1.0, args.polar_gauge, n2, 0
        for rnd in range(args.poly_lm_rounds):
            if elapsed() > 0.85 * wall_s:
                entry['stopped_by_wall'] = True
                break
            W_prev = fld.W.detach().clone()
            # collocation = uniform points + the worst |K^2+P| points of a large
            # sample (tail-aware), so the barrier sees where K approaches degeneracy
            with torch.no_grad():
                xp = man.sample(20000)
                rp = acs_residual(geo, K_of(fld), xp)
                hard = xp[torch.topk(rp, 64).indices]
            xs = torch.cat([man.sample(n_pts - 64), hard], dim=0)
            xs = xs[torch.randperm(len(xs))]
            caches = [PolyCache(man, fld, xs[i:i + 32]) for i in range(0, len(xs), 32)]
            mb, r0 = args.polar_barrier, args.polar_r0
            res = lambda W, c: polar_residual(fld, c, W, w_gauge=wg, mu_bar=mb, r0=r0)
            if fld.W.numel() <= args.dense_limit:
                jac = lambda W, c: polar_residual_and_jacobian(fld, c, W, w_gauge=wg, mu_bar=mb, r0=r0)
                W, lg = lm_minimise(res, fld.W.detach(), caches, n_iter=args.polar_lm_iters, lam0=lam,
                                    jac_fn=jac)
            else:                      # matrix-free: factorised J, CG normal equations
                fac = lambda W, c: polar_factorised(fld, c, W, w_gauge=wg, mu_bar=mb, r0=r0)
                W, lg = lm_cg_minimise(res, fac, fld.W.detach(), caches, n_iter=args.polar_lm_iters,
                                       lam0=lam, cg_iters=args.cg_iters)
            with torch.no_grad():
                fld.W.copy_(W)
            lam = min(max(lg[-1]['lam'], 1e-10), 1e6)   # capped, see run_poly_cs
            n2, c2 = held_out(fld)
            with torch.no_grad():                       # cheap dense admissibility check
                xd = man.sample(4000)
                dense_r = acs_residual(geo, K_of(fld), xd).max().item()
            a = adv_acs(fld) if (rnd % 3 == 2 and dense_r < 0.9) else float('nan')
            rolled = (not math.isfinite(n2)) or dense_r >= 0.9 or (math.isfinite(a) and a >= 0.95)
            if rolled:                                  # smaller steps, firmer gauge
                with torch.no_grad():
                    fld.W.copy_(W_prev)
                lam = min(max(lam, 1e-3) * 100.0, 1e6)
                wg = min(wg * 2.0, 10.0)
            rec = {'degree': deg, 'round': rnd + 1, 'n_points': n_pts, 'train_loss': lg[-1]['loss'],
                   'N2': n2, 'C2': c2, 'adv_acs_sup': a, 'gauge': wg, 'rolled_back': rolled,
                   'elapsed_min': elapsed() / 60}
            entry['lm'].append(rec); log_fn(dict(phase='polar_lm', **rec))
            print(f"[{tag}] polar deg {deg} round {rnd+1} ({n_pts} pts) | train {lg[-1]['loss']:.3e} | "
                  f"held-out |N(J*)|^2 {n2:.4e} |K^2+P|^2 {c2:.2e} | adv |K^2+P| {a:.3f}"
                  f"{' ROLLED BACK' if rolled else ''} | {elapsed()/60:.1f} min", flush=True)
            if rolled:
                continue
            if n2 < best / 1.1:
                best, since = n2, 0
            else:
                since += 1
            if since >= args.patience:
                break
        # floor measurement on the exact ACS J*
        Js = polar_corrected(geo, K_of(fld))
        x = man.sample(2000)
        nn_ = torch.cat([nij_norm(geo, Js, x[i:i + 250]).detach() for i in range(0, 2000, 250)])
        advN, _ = adversarial_sup(geo, lambda y: nij_norm(geo, Js, y), 64, 60)
        entry.update({'floor_N2': (nn_ ** 2).mean().item(), 'N_mean': nn_.mean().item(),
                      'N_sup_dense': nn_.max().item(), 'adv_N_sup': advN,
                      'final_adv_acs_sup': adv_acs(fld), 'elapsed_min': elapsed() / 60})
        ladder.append(entry); log_fn(entry)
        print(f"[{tag}] polar deg {deg} FLOOR: E|N(J*)|^2 {entry['floor_N2']:.4e}, mean |N| {entry['N_mean']:.3e}, "
              f"adv sup |N| {advN:.3e}, adv |K^2+P| {entry['final_adv_acs_sup']:.3f}", flush=True)
        prev = fld
        if advN < args.poly_adv_tol:
            print(f"[{tag}] polar deg {deg}: adversarial sup |N(J*)| < {args.poly_adv_tol:g} -> CS found", flush=True)
            break
        if elapsed() > 0.85 * wall_s:
            break
    return prev, K_of(prev), ladder


# ---------------------------------------------------------------------------
# Conjugation with a POLYNOMIAL generator + LM (`--cs-mode conjpoly`)
# ---------------------------------------------------------------------------

def run_conjpoly_ladder(args, geo, base_fn, base_info, tag, log_fn, elapsed):
    """J = e^B J_b e^{-B}, B = P M_W P with M_W in an orthonormalised degree-<=k
    basis, J_b the unit's own certified ACS. Exact ACS by construction and
    g = e^B always invertible, so there is no degenerate set to exploit. W = 0
    at the first degree (J = J_b); later degrees copy the previous generator.
    Dense LM (or CG-LM above --dense-limit) on resampled collocation sets plus
    the 64 worst-|N| points of a 2000-point sample (tail-aware). Returns the
    field, the J function, and the per-degree log with dense and adversarial
    sup of |N|."""
    from extra.campaign_unit import adversarial_sup, nij_norm
    from src.models.poly_field import (PolyField, ConjCache, conj_residual,
                                       conj_residual_and_jacobian, _expm_ss)
    man, N = geo.manifold, geo.N
    degrees = [int(k) for k in str(args.poly_degrees).split(',')]
    rows = geo.d * (geo.d * (geo.d - 1) // 2)
    ladder, prev = [], None
    wall_s = args.wall_minutes * 60.0

    def J_of(fld):
        def J(x):
            P = man.projector(x)
            B = P @ fld(x)[:, 0] @ P
            return _expm_ss(B) @ base_fn(x) @ _expm_ss(-B)
        return J

    def held_out(fld, n=600):
        c = ConjCache(man, fld, man.sample(n), base_fn)
        with torch.no_grad():
            r = conj_residual(fld, c)
        return (r ** 2).sum().item()

    for deg in degrees:
        fld = PolyField(N, deg)
        nf = fld.orthonormalize(man)
        if prev is not None:
            fit_field(fld, lambda x: prev(x)[:, 0], man, n=max(4000, 4 * nf))
        elif args.init_poly_from:
            ck = torch.load(os.path.join(args.init_poly_from, 'final.pt'), map_location='cpu',
                            weights_only=False)
            if ck['degree'] != deg:
                raise ValueError(f"--init-poly-from has degree {ck['degree']}, ladder starts at {deg}")
            fld.load_state_dict(ck['model'])
            nf = fld.n_feat
        n2 = held_out(fld)
        entry = {'phase': 'conjpoly_ladder', 'degree': deg, 'n_feat': nf, 'n_params': fld.W.numel(),
                 'init_N2': n2, 'lm': []}
        print(f"[{tag}] conjpoly deg {deg}: {nf} basis fns, {fld.W.numel()} params, held-out |N|^2 {n2:.4e}",
              flush=True)
        n_pts = max(args.polar_lm_points, int(args.poly_lm_oversample * fld.W.numel() / rows))
        n_pts = 32 * ((n_pts + 31) // 32)
        lam, best, since = 1.0, float('inf'), 0
        for rnd in range(args.poly_lm_rounds):
            if elapsed() > 0.85 * wall_s:
                entry['stopped_by_wall'] = True
                break
            with torch.no_grad():
                xp = man.sample(2000)
                cp = ConjCache(man, fld, xp, base_fn)
                rp = conj_residual(fld, cp).view(2000, -1).pow(2).sum(1)
                hard = xp[torch.topk(rp, 64).indices]
            xs = torch.cat([man.sample(n_pts - 64), hard], dim=0)[torch.randperm(n_pts)]
            caches = [ConjCache(man, fld, xs[i:i + 32], base_fn) for i in range(0, n_pts, 32)]
            # L^p objective by IRLS: weights (|N(x)| / rms)^(p-2), p rising 2 -> lp_max
            p_now = 2.0 + (args.lp_max - 2.0) * min(1.0, rnd / max(1, args.lp_ramp))
            if p_now > 2.0:
                with torch.no_grad():
                    nrm = [conj_residual(fld, c).view(c.B, -1).pow(2).sum(1) * c.B for c in caches]
                    rms = torch.cat(nrm).mean().sqrt()
                    for c, n_ in zip(caches, nrm):
                        c.w = (n_.sqrt() / rms).clamp_min(1e-3) ** (p_now - 2.0)
                    wmean = torch.cat([c.w for c in caches]).mean()
                    for c in caches:
                        c.w = c.w / wmean
            res = lambda W, c: conj_residual(fld, c, W, w_reg=args.conj_reg)
            jac = lambda W, c: conj_residual_and_jacobian(fld, c, W, w_reg=args.conj_reg)
            W, lg = lm_minimise(res, fld.W.detach(), caches, n_iter=args.polar_lm_iters, lam0=lam, jac_fn=jac)
            with torch.no_grad():
                fld.W.copy_(W)
            lam = min(max(lg[-1]['lam'], 1e-10), 1e6)   # capped, see run_poly_cs
            n2 = held_out(fld)
            with torch.no_grad():                       # sampled-sup proxy on 2000 fresh points
                xq = man.sample(2000)
                cq = ConjCache(man, fld, xq, base_fn)
                supq = conj_residual(fld, cq).view(2000, -1).pow(2).sum(1).mul(2000).sqrt().max().item()
            rec = {'degree': deg, 'round': rnd + 1, 'n_points': n_pts, 'train_loss': lg[-1]['loss'],
                   'N2': n2, 'N_sup_sampled': supq, 'p': p_now, 'W_norm': fld.W.norm().item(),
                   'elapsed_min': elapsed() / 60}
            entry['lm'].append(rec); log_fn(dict(phase='conjpoly_lm', **rec))
            print(f"[{tag}] conjpoly deg {deg} round {rnd+1} ({n_pts} pts, p={p_now:.1f}) | train {lg[-1]['loss']:.3e} | "
                  f"held-out |N|^2 {n2:.4e} | sampled sup |N| {supq:.3e} | |W| {rec['W_norm']:.2f} | "
                  f"{elapsed()/60:.1f} min", flush=True)
            score = supq if args.lp_max > 2 else n2      # an L^p run is judged by the sup
            if score < best / 1.1:
                best, since = score, 0
            else:
                since += 1
            if since >= args.patience:
                break
        J = J_of(fld)
        x = man.sample(2000)
        nn_ = torch.cat([nij_norm(geo, J, x[i:i + 250]).detach() for i in range(0, 2000, 250)])
        advN, _ = adversarial_sup(geo, lambda y: nij_norm(geo, J, y), 64, 60)
        entry.update({'floor_N2': (nn_ ** 2).mean().item(), 'N_mean': nn_.mean().item(),
                      'N_sup_dense': nn_.max().item(), 'adv_N_sup': advN, 'elapsed_min': elapsed() / 60})
        ladder.append(entry); log_fn(entry)
        print(f"[{tag}] conjpoly deg {deg} FLOOR: E|N|^2 {entry['floor_N2']:.4e}, mean |N| {entry['N_mean']:.3e}, "
              f"dense sup {entry['N_sup_dense']:.3e}, adv sup |N| {advN:.3e}", flush=True)
        prev = fld
        if advN < args.poly_adv_tol or elapsed() > 0.85 * wall_s:
            break
    return prev, J_of(prev), ladder
