"""Aggregate campaign_unit results into the before/after table.

Usage:  python extra/campaign_summary.py RESULTS_ROOT [--all]
Scans RESULTS_ROOT/*/result.json; per (case, task) keeps the best seed by the
chart-protocol score (--all lists every run).
"""
from __future__ import annotations

import os
import sys
import json
import glob

ORDER = ['s2', 's4', 's6', 's1xs1', 's2xs2', 's3xs1', 's3xs3', 's4xs2', 's5xs1']
NAMES = {'s2': 'S2', 's4': 'S4', 's6': 'S6', 's1xs1': 'S1xS1', 's2xs2': 'S2xS2',
         's3xs1': 'S3xS1', 's3xs3': 'S3xS3', 's4xs2': 'S4xS2', 's5xs1': 'S5xS1'}


def main():
    root = sys.argv[1]
    show_all = '--all' in sys.argv
    runs = []
    for p in glob.glob(os.path.join(root, '*', 'result.json')):
        with open(p) as f:
            r = json.load(f)
        r['_dir'] = os.path.basename(os.path.dirname(p))
        if 'dense' not in r:                      # e.g. status 'no_certified_base' (S⁴ CS)
            r['_score'] = float('inf')
            r['dense'] = {'acs_mean': float('nan'), 'acs_sup': float('nan'),
                          'nij_mean': float('nan'), 'nij_sup': float('nan')}
        else:
            r['_score'] = r['dense']['nij_mean'] if r['task'] == 'cs' else \
                (r.get('adv_acs_sup') if r.get('adv_acs_sup') is not None else r['dense']['acs_sup'])
        runs.append(r)
    if not runs:
        print("no result.json found under", root)
        return

    best = {}
    for r in runs:
        key = (r['case'], r['task'])
        if key not in best or r['_score'] < best[key]['_score']:
            best[key] = r

    print("| case | task | formulation | before (chart) | chart protocol | intrinsic mean / sup | **adversarial sup** | verdict | class | steps | seed |")
    print("|---|---|---|---:|---:|---:|---:|:-:|:-:|---:|---:|")
    for case in ORDER:
        for task in ('acs', 'cs'):
            r = best.get((case, task))
            if r is None:
                print(f"| {NAMES[case]} | {task.upper()} | — | — | (no result) | | | | | | |")
                continue
            if r.get('status') == 'no_certified_base':
                print(f"| {NAMES[case]} | {task.upper()} | {r['formulation']} | {r['table_before']:.3g} | — | — | "
                      f"{r['base']['adv_acs_sup_raw']:.2f} (base) | **refused: no certified ACS** | | | {r['seed']} |")
                continue
            d = r['dense']
            if task == 'acs':
                verdict = '**certified exact ACS**' if r.get('certified_exact_acs') else 'NOT certified (defects)'
                adv = r.get('adv_acs_sup', float('nan'))
                chart = r.get('chart_acs', float('nan'))
                intr, sup = d['acs_mean'], d['acs_sup']
            else:
                verdict = 'exact ACS by constr.' if r.get('exact_acs_by_construction') else 'free field'
                adv = r.get('adv_nij_sup', float('nan'))
                chart = r.get('chart_cs', float('nan'))
                intr, sup = d['nij_mean'], d['nij_sup']
            print(f"| {NAMES[case]} | {task.upper()} | {r['formulation']} | {r['table_before']:.3g} | "
                  f"{chart:.2e} | {intr:.2e} / {sup:.2e} | **{adv:.2e}** | {verdict} | "
                  f"{r.get('feat') or '—'} | {r['steps_done']} | {r['seed']} |")

    if show_all:
        print("\nAll runs:")
        for r in sorted(runs, key=lambda r: (ORDER.index(r['case']), r['task'], r['seed'])):
            print(f"  {r['_dir']}: score {r['_score']:.3e}  dense {r['dense']}  "
                  f"wall-stopped={r['stopped_by_wall']}  {r['elapsed_min']:.0f} min")


if __name__ == '__main__':
    main()
