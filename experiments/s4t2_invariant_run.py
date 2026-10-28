"""Unbiased search for a T^2-invariant complex structure on the round S^4 x T^2
(reduced model of extra/s4t2_invariant.py). ACS stage: polynomial classes of
degree 1, 2 in x from a random start, Levenberg-Marquardt on |K^2 + P|^2 (plain;
retry with the orientation penalty), certified by the adversarial sup < 1 and
polar-corrected. CS stage: degree 2, 3, 4 free fields from the own ACS, LM on
E|N_K|^2 + 10 E|K^2+P|^2 with resampled collocation; certificate of the
polar-corrected field (dense and adversarial sup of |N|).

    .venv/bin/python experiments/s4t2_invariant_run.py --out outputs/paper_rerun/runs/s4t2inv_s7 --seed 7
"""
from __future__ import annotations

import os
import sys
import json
import time
import math
import argparse

import torch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from extra.s4t2_invariant import (MAN, XPolyField, K_of, acs_residual, nij_norm, polar, adversarial_sup,
                                  orientation_rows, oriented_frame)
from src.models.poly_field import (PolyCache, fit_field, poly_residual_and_jacobian, poly_split_losses,
                                   _factorised_jacobian)
from src.optim.lm import lm_minimise

DTYPE = torch.float64


def caches_for(field, n_pts):
    out = []
    x = MAN.sample(n_pts)
    for i in range(0, n_pts, 32):
        c = PolyCache(MAN, field, x[i:i + 32])
        c.E = oriented_frame(c.E, c.x)
        out.append(c)
    return out


def held_out(field, n=600):
    c = PolyCache(MAN, field, MAN.sample(n))
    with torch.no_grad():
        n2, c2 = poly_split_losses(field, c)
    return n2.item(), c2.item()


def acs_stage(seed, log, rounds=12, w=10.0, mu=1.0):
    ladder = []
    for deg in (1, 2):
        for attempt in ('plain', 'orient'):
            torch.manual_seed(seed * 100 + deg * 10 + (attempt == 'orient'))
            fld = XPolyField(deg)
            s = None
            if attempt == 'orient':
                c0 = caches_for(fld, 512)
                with torch.no_grad():
                    o = torch.cat([orientation_rows(K_of(fld)(c.x), c.E, 1.0, 1.0)[1] for c in c0])
                s = 1.0 if o.mean() >= 0 else -1.0

            def rows(K, DK, P, DP, E, s=s):
                d = E.shape[-1]
                C = torch.einsum('bpa,bpq,bqc->bac', E, K @ K, E) + torch.eye(d, dtype=K.dtype)
                r = [C.reshape(K.shape[0], -1) * math.sqrt(w)]
                if s is not None:
                    r.append(orientation_rows(K, E, s, mu)[0])
                return torch.cat(r, 1)
            lam = 1e-2
            for rnd in range(rounds):
                cs = caches_for(fld, 256)
                W, lg = lm_minimise(lambda W, c: _factorised_jacobian(fld, c, W, rows)[0], fld.W.detach(), cs,
                                    n_iter=6, lam0=lam,
                                    jac_fn=lambda W, c: _factorised_jacobian(fld, c, W, rows))
                with torch.no_grad():
                    fld.W.copy_(W)
                lam = min(max(lg[-1]['lam'], 1e-12), 1e6)
            Kf = K_of(fld)
            adv = adversarial_sup(lambda p: acs_residual(Kf, p), 64, 50)
            with torch.no_grad():
                xs = MAN.sample(2000)
                r = acs_residual(Kf, xs)
            e = {'degree': deg, 'attempt': attempt, 'n_params': fld.W.numel(), 'acs_mean': r.mean().item(),
                 'acs_max': r.max().item(), 'adv_acs_sup': adv}
            ladder.append(e)
            log(f"ACS deg {deg} {attempt:6s}: mean |K^2+P| {e['acs_mean']:.2e} max {e['acs_max']:.2e} adv {adv:.3e} "
                f"-> {'ACCEPT' if adv < 1 else 'reject'}")
            if adv < 1:
                return fld, ladder
    return None, ladder


def certify_field(J_fn, n=3000):
    xs = MAN.sample(n)
    N = torch.cat([nij_norm(J_fn, xs[i:i + 300]).detach() for i in range(0, n, 300)])
    adv = adversarial_sup(lambda p: nij_norm(J_fn, p), 64, 60)
    with torch.no_grad():
        Jx = J_fn(xs[:500]); P = MAN.projector(xs[:500])
        acs = (Jx @ Jx + P).flatten(1).norm(dim=1).max().item()
    return {'N_mean': N.mean().item(), 'N_max': N.max().item(), 'N_adv': adv, 'acs_residual_max': acs}


def cs_stage(Jstar, log, degrees=(2, 3, 4), rounds=30, w=10.0, stall=1.5, patience=4, elapsed=None):
    prev, out = None, []
    for deg in degrees:
        fld = XPolyField(deg)
        if prev is None:
            fit = fit_field(fld, Jstar, MAN, n=max(4000, 3 * fld.n_feat))
        else:
            fld.load_from(prev); fit = None
        n2, c2 = held_out(fld)
        rows_pp = 6 * 15 + 36
        n_pts = max(128, int(2 * fld.W.numel() / rows_pp)); n_pts = 32 * ((n_pts + 31) // 32)
        e = {'degree': deg, 'n_params': fld.W.numel(), 'init_fit_sup': fit, 'init_N2': n2, 'init_C2': c2, 'lm': []}
        log(f"CS deg {deg}: {fld.W.numel()} params, init held-out N2 {n2:.3e} C2 {c2:.3e} ({n_pts} pts/round)")
        lam, best, since = 1e-2, float('inf'), 0
        for rnd in range(rounds):
            cs = caches_for(fld, n_pts)
            W, lg = lm_minimise(lambda W, c: poly_residual_and_jacobian(fld, c, W, w)[0], fld.W.detach(), cs,
                                n_iter=8, lam0=lam, jac_fn=lambda W, c: poly_residual_and_jacobian(fld, c, W, w))
            with torch.no_grad():
                fld.W.copy_(W)
            lam = min(max(lg[-1]['lam'], 1e-12), 1e6)
            n2, c2 = held_out(fld)
            score = n2 + w * c2
            rec = {'round': rnd + 1, 'train': lg[-1]['loss'], 'N2': n2, 'C2': c2, 'elapsed_min': elapsed()}
            e['lm'].append(rec)
            log(f"  deg {deg} round {rnd+1}: train {lg[-1]['loss']:.3e} held-out N2 {n2:.3e} C2 {c2:.3e} "
                f"({elapsed():.1f} min)")
            if score < 1e-20:
                break
            if score < best / stall:
                best, since = score, 0
            else:
                since += 1
            if since >= patience:
                break
        Kf = K_of(fld)
        e['adv_acs_sup'] = adversarial_sup(lambda p: acs_residual(Kf, p), 64, 50)
        if e['adv_acs_sup'] < 1:
            e.update(certify_field(polar(Kf)))
            log(f"CS deg {deg} RESULT: adv |K^2+P| {e['adv_acs_sup']:.2e}; N(J*) mean {e['N_mean']:.3e} "
                f"max {e['N_max']:.3e} adv {e['N_adv']:.3e}")
        else:
            log(f"CS deg {deg} RESULT: free field not admissible (adv |K^2+P| {e['adv_acs_sup']:.2f})")
        out.append(e)
        prev = fld
        if e.get('N_adv', 1e9) < 1e-6:
            break
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--degrees', default='2,3,4')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    torch.set_num_threads(int(os.environ.get('OMP_NUM_THREADS', 2)))
    torch.manual_seed(a.seed)
    t0 = time.time()
    elapsed = lambda: (time.time() - t0) / 60
    logf = open(os.path.join(a.out, 'log.txt'), 'w')

    def log(msg):
        print(msg, flush=True); logf.write(msg + '\n'); logf.flush()
    res = {'seed': a.seed}
    fld, ladder = acs_stage(a.seed, log)
    res['acs_ladder'] = ladder
    if fld is None:
        res['status'] = 'no certified ACS'
    else:
        Jstar = polar(K_of(fld))
        res['acs'] = certify_field(Jstar)
        res['acs']['adv_acs_raw'] = ladder[-1]['adv_acs_sup']
        log(f"own ACS J*: |J*^2+P| {res['acs']['acs_residual_max']:.1e}, N mean {res['acs']['N_mean']:.3f} "
            f"adv {res['acs']['N_adv']:.3f}")
        res['cs'] = cs_stage(Jstar, log, degrees=[int(k) for k in a.degrees.split(',')], elapsed=elapsed)
    res['elapsed_min'] = elapsed()
    json.dump(res, open(os.path.join(a.out, 'done.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
