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
        r['_score'] = r['chart_cs'] if r['task'] == 'cs' else r['chart_acs']
        runs.append(r)
    if not runs:
        print("no result.json found under", root)
        return

    best = {}
    for r in runs:
        key = (r['case'], r['task'])
        if key not in best or r['_score'] < best[key]['_score']:
            best[key] = r

    print("| case | task | formulation | before | **after (chart protocol)** | intrinsic mean / sup | certified exact ACS | steps | seed |")
    print("|---|---|---|---:|---:|---:|:-:|---:|---:|")
    for case in ORDER:
        for task in ('acs', 'cs'):
            r = best.get((case, task))
            if r is None:
                print(f"| {NAMES[case]} | {task.upper()} | — | — | (no result) | | | | |")
                continue
            d = r['dense']
            cert = r.get('certified_exact_acs')
            cert_s = '✓' if cert else ('exact by constr.' if r.get('exact_acs_by_construction') else '—')
            intr = d['nij_mean'] if task == 'cs' else d['acs_mean']
            sup = d['nij_sup'] if task == 'cs' else d['acs_sup']
            print(f"| {NAMES[case]} | {task.upper()} | {r['formulation']} | {r['table_before']:.3g} | "
                  f"**{r['_score']:.3e}** | {intr:.3e} / {sup:.3e} | {cert_s} | "
                  f"{r['steps_done']} | {r['seed']} |")

    if show_all:
        print("\nAll runs:")
        for r in sorted(runs, key=lambda r: (ORDER.index(r['case']), r['task'], r['seed'])):
            print(f"  {r['_dir']}: score {r['_score']:.3e}  dense {r['dense']}  "
                  f"wall-stopped={r['stopped_by_wall']}  {r['elapsed_min']:.0f} min")


if __name__ == '__main__':
    main()
