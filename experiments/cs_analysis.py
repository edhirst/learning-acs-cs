"""Post-hoc analysis of the CS runs of the report (never used in training).

  (i)  defect concentration: |N| of the exact ACS J of each CS run on fresh
       points, sorted; the share of E|N|^2 carried by the worst q-fraction of
       the points (the L^2 loophole of the report, Section 7.1);
  (ii) the complex structures found by the polynomial mode: orthogonality
       defect |J^T J - P|, singular values of J on T, homogeneity of the defect.

    .venv/bin/python experiments/cs_analysis.py --runs outputs/paper_rerun/runs --out outputs/paper_rerun/runs/cs_analysis
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
from extra.campaign_unit import (Geometry, Net, load_acs_field, polar_corrected, nij_norm, make_field,
                                 acs_residual, DTYPE)
from src.models.poly_field import PolyField, _expm_ss

FRACS = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]


def field_of(geo, run):
    """Exact ACS J of a finished CS run (None if the run has no exact field)."""
    res = json.load(open(os.path.join(run, 'result.json')))
    ck = torch.load(os.path.join(run, 'final.pt'), map_location='cpu', weights_only=False)
    man = geo.manifold
    mode = res.get('cs_mode', 'conj')
    if mode == 'conj':
        model = Net(geo.N, geo.N, ck['stages'], ck['hidden'], 0.5).to(DTYPE)
        model.load_state_dict(ck['model']); model.eval()
        for p in model.parameters():
            p.requires_grad_(False)
        base = polar_corrected(geo, load_acs_field(geo, ck['base_from']))
        J, _ = make_field(geo, 'cs', model, base_fn=base, conj=ck.get('conj', 'exp'))
        return J, mode, res
    if mode == 'conjpoly':
        fld = PolyField(geo.N, ck['degree'])
        fld.orthonormalize(man)
        fld.load_state_dict(ck['model'])
        for p in fld.parameters():
            p.requires_grad_(False)
        base = polar_corrected(geo, load_acs_field(geo, ck['base_from']))

        def J(x):
            P = man.projector(x)
            B = P @ fld(x)[:, 0] @ P
            return _expm_ss(B) @ base(x) @ _expm_ss(-B)
        return J, mode, res
    if mode in ('poly', 'polar'):
        if not res.get('certified_exact_acs'):
            return None, mode, res
        fld = PolyField(geo.N, ck['degree'], blocks=ck.get('blocks'), trilin=ck.get('trilin', False))
        if mode == 'polar':
            fld.orthonormalize(man)
        fld.load_state_dict(ck['model'])
        for p in fld.parameters():
            p.requires_grad_(False)
        K = lambda x: man.projector(x) @ fld(x)[:, 0] @ man.projector(x)
        return polar_corrected(geo, K), mode, res
    return None, mode, res


def nij_values(geo, J, n, chunk=250, seed=0):
    torch.manual_seed(seed)
    x = geo.sample(n)
    vals = [nij_norm(geo, J, x[i:i + chunk]).detach() for i in range(0, n, chunk)]
    return x, torch.cat(vals).numpy()


def concentration(v):
    s = np.sort(v ** 2)[::-1]
    c = np.cumsum(s) / s.sum()
    n = len(s)
    return {f'{q:g}': float(c[max(0, int(round(q * n)) - 1)]) for q in FRACS}


def structure_stats(geo, J, x):
    with torch.no_grad():
        Jx = J(x)
        P = geo.manifold.projector(x)
        orth = (Jx.transpose(1, 2) @ Jx - P).flatten(1).norm(dim=1)
        ev, evec = torch.linalg.eigh(P)
        E = evec[:, :, -geo.d:]
        sv = torch.linalg.svdvals(E.transpose(1, 2) @ Jx @ E)
        acs = acs_residual(geo, J, x)
    return {'orth_defect_mean': orth.mean().item(), 'orth_defect_max': orth.max().item(),
            'orth_defect_std': orth.std().item(), 'sv_min': sv.min().item(), 'sv_max': sv.max().item(),
            'acs_residual_max': acs.max().item()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--runs', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--n', type=int, default=12000)
    ap.add_argument('--n6', type=int, default=6000, help='points for the six-dimensional MLP runs')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    torch.set_num_threads(int(os.environ.get('OMP_NUM_THREADS', 2)))
    names = sorted(n for n in os.listdir(a.runs)
                   if n.startswith(('cspoly_', 'csconj', 'conjpoly_', 'lp_', 'polar_')) and
                   os.path.exists(os.path.join(a.runs, n, 'result.json')) and
                   os.path.exists(os.path.join(a.runs, n, 'final.pt')))
    # reference: the exact ACS J* of the ACS units (the starting point of every CS run)
    refs = [f'acs_{c}_s7' for c in ('s6', 's4xs2', 's4xs1xs1', 's2xs2', 's3xs1', 's3xs3', 's5xs1')]
    out = {}
    for name in refs + names:
        run = os.path.join(a.runs, name)
        if not os.path.exists(os.path.join(run, 'result.json')):
            continue
        res = json.load(open(os.path.join(run, 'result.json')))
        geo = Geometry(res['case'], 'ambient')
        t0 = time.time()
        try:
            if name.startswith('acs_'):
                if not res.get('certified_exact_acs'):
                    continue
                J, mode = polar_corrected(geo, load_acs_field(geo, run)), 'acs'
            else:
                J, mode, res = field_of(geo, run)
            if J is None:
                out[name] = {'case': res['case'], 'mode': mode, 'skipped': 'no exact field'}
                continue
            n = a.n6 if (geo.d == 6 and mode in ('conj', 'acs', 'conjpoly')) else a.n
            x, v = nij_values(geo, J, n)
        except Exception as e:                       # keep going: one bad run must not stop the table
            out[name] = {'case': res.get('case'), 'error': f'{type(e).__name__}: {e}'}
            print(name, 'ERROR', e, flush=True)
            continue
        rec = {'case': res['case'], 'mode': mode, 'n': int(n), 'mean': float(v.mean()),
               'rms': float(np.sqrt((v ** 2).mean())), 'max': float(v.max()),
               'quantiles': {f'{q:g}': float(np.quantile(v, q)) for q in (0.5, 0.9, 0.99, 0.999)},
               'share_top': concentration(v), 'adv_nij_sup': res.get('adv_nij_sup')}
        if mode in ('poly', 'polar') or name.startswith('cspoly_'):
            rec.update(structure_stats(geo, J, x[:2000]))
        rec['seconds'] = time.time() - t0
        out[name] = rec
        print(f"{name:28s} mean {rec['mean']:.3e} max {rec['max']:.3e} top5% share {rec['share_top']['0.05']:.2f}"
              f" ({rec['seconds']:.0f}s)", flush=True)
        json.dump(out, open(os.path.join(a.out, 'partial.json'), 'w'), indent=1)
    json.dump(out, open(os.path.join(a.out, 'done.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
