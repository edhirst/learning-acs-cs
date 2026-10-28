"""T^2-invariant conjugation with a POLYNOMIAL generator, solved by
Levenberg-Marquardt, on the reduced model of extra/s4t2_invariant.py:
J = e^B J_b e^-B, B = P M_W(x) P with M_W of degree k in x in S^4 (k = 1, 2, 3),
J_b the unit's own exact T^2-invariant ACS. Exact ACS by construction; the same
protocol as the polynomial-generator runs on S^6 (25 residual rows per unknown,
the 64 worst points of a 2000-point sample added every round, patience 5).

    .venv/bin/python experiments/s4t2_invariant_conjpoly.py --out outputs/paper_rerun/runs/s4t2invcp_s7 --seed 7
"""
from __future__ import annotations

import os
import sys
import json
import time
import argparse

import torch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extra.s4t2_invariant import MAN, XPolyField, K_of, polar, nij_norm, adversarial_sup
from src.models.poly_field import ConjCache, conj_residual, conj_residual_and_jacobian, _expm_ss, fit_field
from src.optim.lm import lm_minimise
from s4t2_invariant_run import acs_stage

DTYPE = torch.float64


def J_of(fld, Jb):
    def J(p):
        P = MAN.projector(p)
        B = P @ fld(p)[:, 0] @ P
        return _expm_ss(B) @ Jb(p) @ _expm_ss(-B)
    return J


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--degrees', default='1,2,3')
    ap.add_argument('--rounds', type=int, default=40)
    ap.add_argument('--patience', type=int, default=5)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    torch.set_num_threads(int(os.environ.get('OMP_NUM_THREADS', 2)))
    t0 = time.time()
    el = lambda: (time.time() - t0) / 60
    logf = open(os.path.join(a.out, 'log.txt'), 'w')

    def log(msg):
        print(msg, flush=True); logf.write(msg + '\n'); logf.flush()
    fld0, ladder = acs_stage(a.seed, log)
    for p in fld0.parameters():
        p.requires_grad_(False)
    Jb = polar(K_of(fld0))
    out, prev = [], None
    for deg in [int(k) for k in a.degrees.split(',')]:
        fld = XPolyField(deg)
        with torch.no_grad():
            fld.W.zero_()
        if prev is not None:
            fld.load_from(prev)
        n_pts = max(384, int(25 * fld.W.numel() / 90)); n_pts = 32 * ((n_pts + 31) // 32)
        e = {'degree': deg, 'n_params': fld.W.numel(), 'lm': []}
        log(f"conjpoly deg {deg}: {fld.W.numel()} params, {n_pts} pts/round")
        lam, best, since = 1.0, float('inf'), 0
        for rnd in range(a.rounds):
            with torch.no_grad():
                xp = MAN.sample(2000)
                cp = ConjCache(MAN, fld, xp, Jb)
                rp = conj_residual(fld, cp).view(2000, -1).pow(2).sum(1)
                hard = xp[torch.topk(rp, 64).indices]
            xs = torch.cat([MAN.sample(n_pts - 64), hard], 0)[torch.randperm(n_pts)]
            caches = [ConjCache(MAN, fld, xs[i:i + 32], Jb) for i in range(0, n_pts, 32)]
            W, lg = lm_minimise(lambda W, c: conj_residual(fld, c, W), fld.W.detach(), caches, n_iter=3, lam0=lam,
                                jac_fn=lambda W, c: conj_residual_and_jacobian(fld, c, W))
            with torch.no_grad():
                fld.W.copy_(W)
            lam = min(max(lg[-1]['lam'], 1e-10), 1e6)
            with torch.no_grad():
                xq = MAN.sample(1000)
                cq = ConjCache(MAN, fld, xq, Jb)
                rq = conj_residual(fld, cq).view(1000, -1).pow(2).sum(1).mul(1000)
            n2, supq = rq.mean().item(), rq.sqrt().max().item()
            e['lm'].append({'round': rnd + 1, 'train': lg[-1]['loss'], 'N2': n2, 'N_sup_sampled': supq, 'elapsed_min': el()})
            log(f"  deg {deg} round {rnd+1}: train {lg[-1]['loss']:.3e} held-out E|N|^2 {n2:.3e} sampled sup {supq:.3e} "
                f"({el():.1f} min)")
            if n2 < best / 1.1:
                best, since = n2, 0
            else:
                since += 1
            if since >= a.patience:
                break
        J = J_of(fld, Jb)
        xs = MAN.sample(3000)
        N = torch.cat([nij_norm(J, xs[i:i + 300]).detach() for i in range(0, 3000, 300)])
        e.update({'N_mean': N.mean().item(), 'N_max': N.max().item(),
                  'N_adv': adversarial_sup(lambda p: nij_norm(J, p), 64, 60), 'elapsed_min': el()})
        log(f"conjpoly deg {deg} RESULT: N mean {e['N_mean']:.3e} max {e['N_max']:.3e} adv {e['N_adv']:.3e}")
        out.append(e)
        prev = fld
    json.dump({'seed': a.seed, 'acs_ladder': ladder, 'conjpoly': out, 'elapsed_min': el()},
              open(os.path.join(a.out, 'done.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
