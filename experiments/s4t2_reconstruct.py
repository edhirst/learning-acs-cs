"""Reconstruction of the complex structure on S^4 x T^2 of [HSSS26] on its
fibration model, and physics-informed recovery of admissible period data
(extra/s4t2_reconstruction.py).

    .venv/bin/python experiments/s4t2_reconstruct.py --out outputs/paper_rerun/runs/s4t2_reconstruction

Writes certificate.json (certificate of the data of [HSSS26] and of the best
learned data), learned.json (all random starts) and grids.npz (normalised
determinant over C*, chirality and height-ratio curves) for the figures.
"""
from __future__ import annotations

import os
import sys
import json
import time
import argparse

import numpy as np
import torch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from extra.s4t2_reconstruction import PeriodData, certify, learn_period_data, height_ratio, omega

C128 = torch.complex128


def det_grid(pd, n_r=241, n_th=181, L=5.0):
    u = np.linspace(-L, L, n_r); th = np.linspace(-np.pi, np.pi, n_th)
    U, TH = np.meshgrid(u, th, indexing='ij')
    t = torch.tensor((np.exp(U) * np.exp(1j * TH)).ravel(), dtype=C128)
    return U, TH, pd.det_normalised(t).numpy().reshape(U.shape)


def ratio_curves(pd, at_inf):
    r = np.logspace(-6, np.log10(0.2), 200)
    out = []
    for rr in r:
        t = torch.tensor(rr * np.exp(1j * np.linspace(-np.pi + 1e-3, np.pi - 1e-3, 64)), dtype=C128)
        hr = height_ratio(pd, t, at_inf)
        out.append((hr.min().item(), hr.max().item()))
    return r, np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--seeds', type=int, default=20)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    torch.set_num_threads(int(os.environ.get('OMP_NUM_THREADS', 2)))
    t0 = time.time()
    paper = PeriodData.paper()
    report = {'paper': certify(paper)}
    print(f"[HSSS26] data: certified sup {report['paper']['certified_sup']:.2e}, "
          f"min normalised det {report['paper']['det_normalised_min']:.3f}, "
          f"filling eps {report['paper']['filling_eps_0']:.3g} / {report['paper']['filling_eps_inf']:.3g}", flush=True)
    learned = []
    for seed in range(a.seeds):
        ts = time.time()
        pd, hist = learn_period_data(seed, log=None)
        feas = [h for h in hist if h['feasible']]
        rec = {'seed': seed, 'feasible': bool(feas), 'margin': feas[-1]['delta'] if feas else None,
               'steps': len(hist), 'seconds': time.time() - ts, 'data': pd.as_dict(), 'hist': hist}
        if feas:
            c = certify(pd)
            rec['certificate'] = {k: c[k] for k in ('det_normalised_min', 'chirality_0', 'chirality_inf',
                                                     'filling_eps_0', 'filling_eps_inf', 'certified_sup', 'admissible')}
            a2, a3, c2, c3, b4 = pd.a2, pd.a3, pd.c2, pd.c3, pd.b4
            tt = np.exp(np.linspace(-2, 2, 7)) * np.exp(1j * np.linspace(0, 6, 7))
            e = np.stack([a2 + c2 * tt ** 2, a3 + c3 * tt ** 2, b4 * tt])
            rec['null_conic_defect'] = float(np.abs((e ** 2).sum(0)).max() / (np.abs(e) ** 2).sum(0).max())
            rec['relation_c2_plus_a2'] = abs(c2 + a2) / (abs(a2) + abs(c2))
            rec['relation_c3_minus_a3'] = abs(c3 - a3) / (abs(a3) + abs(c3))
        learned.append(rec)
        print(f"seed {seed:2d}: feasible={rec['feasible']} margin {rec['margin']} "
              f"certified {rec.get('certificate', {}).get('certified_sup')}", flush=True)
    ok = [r for r in learned if r['feasible'] and r['certificate']['admissible']]
    best = max(ok, key=lambda r: (r['margin'], r['certificate']['det_normalised_min'])) if ok else None
    report['learned_summary'] = {'starts': len(learned), 'admissible': len(ok),
                                 'margins': [r['margin'] for r in ok], 'best_seed': best['seed'] if best else None}
    if best:
        pd_best = PeriodData(lam=best['data']['lam'], **{k: complex(*best['data'][k]) for k in
                                                         ('a2', 'a3', 'c2', 'c3', 'b4', 'beta')})
        report['learned_best'] = certify(pd_best)
        print(f"best learned (seed {best['seed']}): certified sup {report['learned_best']['certified_sup']:.2e}, "
              f"min normalised det {report['learned_best']['det_normalised_min']:.3f}", flush=True)
        # like-for-like comparison: both data sets on the same (the wider) annulus
        eps_c = min(report['paper']['annulus_eps'], report['learned_best']['annulus_eps'])
        keys = ('annulus_eps', 'acs_residual_max', 'N_max', 'monodromy_consistency_max',
                'chart_infinity_consistency_max', 'log_transform_consistency_max', 'N_disc_0_max',
                'N_disc_inf_max', 'certified_sup')
        report['common_annulus'] = {nm: {k: c[k] for k in keys} for nm, c in
                                    (('paper', certify(paper, eps=eps_c)), ('learned_best', certify(pd_best, eps=eps_c)))}
        print(f"common annulus (eps {eps_c:.3g}): certified sup [HSSS26] "
              f"{report['common_annulus']['paper']['certified_sup']:.2e}, learned "
              f"{report['common_annulus']['learned_best']['certified_sup']:.2e}", flush=True)
    report['elapsed_min'] = (time.time() - t0) / 60
    json.dump(report, open(os.path.join(a.out, 'certificate.json'), 'w'), indent=1)
    json.dump(learned, open(os.path.join(a.out, 'learned.json'), 'w'), indent=1)
    U, TH, Dp = det_grid(paper)
    grids = {'U': U, 'TH': TH, 'det_paper': Dp}
    if best:
        grids['det_learned'] = det_grid(pd_best)[2]
    for nm, pdx in (('paper', paper),) + ((('learned', pd_best),) if best else ()):
        for end, inf in (('0', False), ('inf', True)):
            r, rc = ratio_curves(pdx, inf)
            grids[f'ratio_r'] = r
            grids[f'ratio_{nm}_{end}'] = rc
    np.savez(os.path.join(a.out, 'grids.npz'), **grids)
    json.dump({'ok': True}, open(os.path.join(a.out, 'done.json'), 'w'))


if __name__ == '__main__':
    main()
