"""Re-run the a-posteriori certification of campaign ACS units from their
saved final.pt. The polar inverse square root uses the differentiable
Newton–Schulz implementation from campaign_unit, valid under the certified
residual bound even when −K² is not symmetric.

Usage:  python extra/campaign_recertify.py RESULTS_ROOT [--n 6000]
Updates result.json in place (keys: corrected, certified_exact_acs,
chart_cs/chart_acs recomputed on the corrected field when admissible) and
prints a summary. Also reports the corrected field's Nijenhuis, i.e. what
the CS column would read for the certified ACS.
"""
from __future__ import annotations

import os
import sys
import json
import glob
import argparse
import numpy as np
import torch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from extra.campaign_unit import (Geometry, load_acs_field, dense_eval,
                                 polar_corrected, chart_protocol)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('root')
    ap.add_argument('--n', type=int, default=6000)
    ap.add_argument('--n-chart', type=int, default=2000)
    args = ap.parse_args()
    torch.manual_seed(123)
    np.random.seed(123)

    for p in sorted(glob.glob(os.path.join(args.root, '*_acs_*', 'result.json'))):
        d = os.path.dirname(p)
        r = json.load(open(p))
        if r.get('schema_version', 1) >= 2:
            print(f"{d}: skipped new ACS experiment; use extra/validate_ambient_acs.py")
            continue
        geo = Geometry(r['case'], r.get('formulation', 'ambient'))
        K_fn = load_acs_field(geo, d)
        raw = dense_eval(geo, 'acs', K_fn, args.n)
        line = f"{os.path.basename(d):16s} raw acs {raw['acs_mean']:.2e}/{raw['acs_sup']:.2e}"
        r['dense'] = raw
        if raw['acs_sup'] < 1.0:
            J_star = polar_corrected(geo, K_fn)
            corr = dense_eval(geo, 'acs', J_star, args.n)
            r['corrected'] = corr
            r['certified_exact_acs'] = corr['acs_sup'] < 1e-8
            cs, acs = chart_protocol(geo, J_star, args.n_chart)
            r['chart_cs'], r['chart_acs'] = cs, acs
            line += (f" | corrected acs {corr['acs_mean']:.1e}/{corr['acs_sup']:.1e}"
                     f"  N(J*) {corr['nij_mean']:.2e}/{corr['nij_sup']:.2e}"
                     f"  certified={r['certified_exact_acs']}  chart_cs(J*) {cs:.3e}")
        else:
            r['corrected'] = None
            r['certified_exact_acs'] = False
            line += "  | not admissible (sup >= 1): defects remain"
        with open(p, 'w') as f:
            json.dump(r, f, indent=2)
        print(line, flush=True)


if __name__ == '__main__':
    main()
