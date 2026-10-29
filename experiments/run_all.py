"""Re-run every experiment of the paper from scratch (dependency-aware scheduler).

    .venv/bin/python experiments/run_all.py [--budget 9] [--dry] [--only PREFIX,...]

Each job is a node with dependencies, a thread allowance and a rough duration;
ready jobs are launched longest-first (critical path) while the total thread
budget allows. A job counts as done when its output marker exists, so the
script can be re-run to resume. Outputs go to $PAPER_OUT (default
outputs/paper_rerun; its records are tracked, its checkpoints are not):
runs/<job>/, logs/<job>.txt, scheduler_status.json; experiments/make_figures.py
then writes figures/ and tables/ there.
"""
from __future__ import annotations

import os
import sys
import json
import time
import argparse
import subprocess

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
PY = os.path.join(ROOT, '.venv', 'bin', 'python')
OUT = os.path.abspath(os.environ.get('PAPER_OUT', os.path.join(ROOT, 'outputs', 'paper_rerun')))
RUNS = os.path.join(OUT, 'runs')
LOGS = os.path.join(OUT, 'logs')
EXP = os.path.join(ROOT, 'experiments')
UNIT = os.path.join(ROOT, 'extra', 'campaign_unit.py')

CASES = ['s2', 's4', 's6', 's1xs1', 's2xs2', 's3xs1', 's3xs3', 's4xs2', 's5xs1', 's4xs1xs1']
SEEDS = [7, 8]
SIX_D = {'s6', 's3xs3', 's4xs2', 's5xs1', 's4xs1xs1'}


def run_dir(name):
    return os.path.join(RUNS, name)


def unit(name, *args, deps=(), threads=1, minutes=10):
    return {'name': name, 'cmd': [PY, UNIT, *map(str, args), '--out', run_dir(name)],
            'deps': list(deps), 'threads': threads, 'minutes': minutes,
            'marker': os.path.join(run_dir(name), 'result.json')}


def script(name, path, *args, deps=(), threads=1, minutes=10, marker=None, soft=False):
    return {'name': name, 'cmd': [PY, path, *map(str, args)], 'deps': list(deps), 'threads': threads,
            'minutes': minutes, 'marker': marker or os.path.join(run_dir(name), 'done.json'), 'soft': soft}


def supervised_jobs():
    """Paper Table 4: supervised fits of the Hopf structure of (C^m - 0)/<z -> e z> (c = 1/2pi,
    extra/ambient_supervised_hopf.py) with the conjugation network of the unsupervised units, from
    the campaign's own certified ACS units (the base orientation is matched automatically). Recipe
    of docs/ambient_supervised_hopf_log.md: 40k Adam steps on S^3 x S^1, 80k on S^5 x S^1, then
    4 x 100 L-BFGS iterations on 2048 points; values + first derivatives ('sobolev') on two bases
    x network seeds, values only on base 7. The 13-h S^5 x S^1 runs have low priority, so that the
    rest of the campaign is not starved of threads.
    On S^3 x S^1 the ACS homotopy classes relative to the target carry a Hopf invariant
    (experiments/hopf_class.py): bases 7 and 10 are in the target's class, base 8 is not (invariant
    +-1), so no conjugation reaches the target from it; its run is kept as the example of this
    obstruction, and Table 4 uses bases 7 and 10. On S^5 x S^1 there is one class per orientation."""
    S = os.path.join(ROOT, 'extra', 'ambient_supervised_hopf.py')
    spec = [  # name, case, base seed, network seed, loss, threads, minutes, priority
        ('hopf_s3xs1_sob_b7', 's3xs1', 7, 7, 'sobolev', 2, 130, 0),
        ('hopf_s3xs1_sob_b8', 's3xs1', 8, 8, 'sobolev', 2, 130, 0),
        ('hopf_s3xs1_sob_b10', 's3xs1', 10, 10, 'sobolev', 2, 130, 0),
        ('hopf_s3xs1_val_b7', 's3xs1', 7, 7, 'value', 1, 15, 0),
        ('hopf_s5xs1_val_b7', 's5xs1', 7, 7, 'value', 2, 70, 0),
        ('hopf_s5xs1_sob_b7', 's5xs1', 7, 7, 'sobolev', 2, 840, -1),
        ('hopf_s5xs1_sob_b8', 's5xs1', 8, 8, 'sobolev', 2, 840, -1),
        ('hopf_s5xs1_sob_b7_s17', 's5xs1', 7, 17, 'sobolev', 2, 840, -1),
        ('hopf_s5xs1_sob_b8_s18', 's5xs1', 8, 18, 'sobolev', 2, 840, -1),
    ]
    J = []
    for name, c, b, seed, loss, th, mins, prio in spec:
        j = script(name, S, '--case', c, '--base-from', run_dir(f'acs_{c}_s{b}'), '--seed', seed,
                   '--loss', loss, '--out', run_dir(name), deps=[f'acs_{c}_s{b}'], threads=th, minutes=mins,
                   marker=os.path.join(run_dir(name), 'result.json'))
        j['prio'] = prio
        J.append(j)
    return J


def jobs():
    J = []
    # ---- A. ACS, harmonic ladder, every case and seed (default settings, finer logging)
    acs_minutes = {'s2': 4, 's1xs1': 6, 's4': 20, 's2xs2': 5, 's3xs1': 12, 's6': 25, 's3xs3': 30,
                   's4xs2': 35, 's5xs1': 30, 's4xs1xs1': 80}
    for c in CASES:
        for s in SEEDS:
            J.append(unit(f'acs_{c}_s{s}', '--case', c, '--task', 'acs', '--seed', s, '--eval-every', 500,
                          minutes=acs_minutes[c]))
    # extra S^3 x S^1 bases for Table 4 (added 2026-10-25): a supervised fit by conjugation needs a base
    # in the homotopy class of the Hopf target, and seed 8 is not in it (experiments/hopf_class.py)
    for s in (9, 10):
        J.append(unit(f'acs_s3xs1_s{s}', '--case', 's3xs1', '--task', 'acs', '--seed', s, '--eval-every', 500,
                      minutes=acs_minutes['s3xs1']))
    # ---- B. CS, unbiased polynomial mode with LM, from the unit's own ACS
    for c in CASES:
        for s in (SEEDS if c not in ('s6', 's4xs2', 's4xs1xs1') else [7]):   # hard cases: failure mode, seed 7
            degrees = '2t' if c == 's4xs1xs1' else '2,3'
            minutes = {'s6': 90, 's4xs2': 150, 's4xs1xs1': 60, 's3xs3': 15, 's5xs1': 20}.get(c, 8)
            J.append(unit(f'cspoly_{c}_s{s}', '--case', c, '--task', 'cs', '--seed', s, '--cs-mode', 'poly',
                          '--poly-degrees', degrees, '--base-from', run_dir(f'acs_{c}_s{s}'),
                          deps=[f'acs_{c}_s{s}'], threads=2 if c in SIX_D else 1, minutes=minutes))
    # ---- C. CS, original conjugation pipeline (baselines), seed 7
    for c, steps, mins in (('s2xs2', 6000, 60), ('s3xs1', 6000, 60), ('s3xs3', 4000, 75),
                           ('s5xs1', 4000, 75), ('s6', 4000, 75), ('s4xs2', 4000, 75)):
        J.append(unit(f'csconjexp_{c}', '--case', c, '--task', 'cs', '--seed', 7, '--conj', 'exp',
                      '--fast-cs', '--steps', steps, '--eval-every', 250,
                      '--lbfgs', 200, '--lbfgs-chunk', 100, '--lbfgs-batch', 1024,
                      '--base-from', run_dir(f'acs_{c}_s7'),
                      deps=[f'acs_{c}_s7'], threads=2 if c in SIX_D else 1, minutes=mins))
    # affine conjugator on every case that has an exponential run (S^3 x S^3 and S^6 added 2026-10-25)
    for c, steps, mins in (('s2xs2', 10000, 30), ('s3xs1', 10000, 30), ('s3xs3', 4000, 60), ('s5xs1', 4000, 45),
                           ('s6', 4000, 45), ('s4xs2', 4000, 45)):
        J.append(unit(f'csconjaff_{c}', '--case', c, '--task', 'cs', '--seed', 7, '--conj', 'affine',
                      '--fast-cs', '--steps', steps, '--eval-every', 250,
                      '--lbfgs', 200, '--lbfgs-chunk', 100, '--lbfgs-batch', 1024,
                      '--base-from', run_dir(f'acs_{c}_s7'),
                      deps=[f'acs_{c}_s7'], threads=2 if c in SIX_D else 1, minutes=mins))
    for c, steps in (('s6', 2500), ('s4xs2', 2500), ('s4xs1xs1', 1500)):
        J.append(unit(f'csconjtail_{c}', '--case', c, '--task', 'cs', '--seed', 7, '--conj', 'exp',
                      '--cs-scale', 0.1, '--fast-cs', '--tail', '--steps', steps, '--eval-every', 250,
                      '--lbfgs', 200, '--lbfgs-chunk', 100, '--lbfgs-batch', 1024,
                      '--base-from', run_dir(f'acs_{c}_s7'), deps=[f'acs_{c}_s7'], threads=2, minutes=150))
    # continuation of the S4xS1xS1 tail run (the only hard case whose supremum kept falling)
    J.append(unit('csconjtail_s4xs1xs1_cont', '--case', 's4xs1xs1', '--task', 'cs', '--seed', 7, '--conj', 'exp',
                  '--cs-scale', 0.1, '--fast-cs', '--tail', '--steps', 4000, '--eval-every', 250, '--lr', 3e-4,
                  '--init-from', run_dir('csconjtail_s4xs1xs1'),
                  '--lbfgs', 200, '--lbfgs-chunk', 100, '--lbfgs-batch', 1024,
                  '--base-from', run_dir('acs_s4xs1xs1_s7'), deps=['csconjtail_s4xs1xs1'], threads=2, minutes=200))
    # ---- D. polynomial-generator conjugation on the hard cases (Table 5)
    J.append(unit('conjpoly_s6', '--case', 's6', '--task', 'cs', '--seed', 7, '--cs-mode', 'conjpoly',
                  '--poly-degrees', '1,2,3', '--poly-lm-rounds', 60, '--poly-lm-oversample', 25,
                  '--base-from', run_dir('acs_s6_s7'), deps=['acs_s6_s7'], threads=2, minutes=200))
    J.append(unit('conjpoly_s4xs2', '--case', 's4xs2', '--task', 'cs', '--seed', 7, '--cs-mode', 'conjpoly',
                  '--poly-degrees', '1,2', '--poly-lm-rounds', 25, '--poly-lm-oversample', 25, '--patience', 3,
                  '--base-from', run_dir('acs_s4xs2_s7'), deps=['acs_s4xs2_s7'], threads=2, minutes=240))
    # ---- E. reconstruction and its convergence study, offline checks
    J.append(script('reconstruction', os.path.join(ROOT, 'extra', 'reconstruct_s6.py'),
                    '--out', run_dir('reconstruction'), threads=2, minutes=5,
                    marker=os.path.join(run_dir('reconstruction'), 'certificate.json')))
    J.append(script('reconstruction_K', os.path.join(EXP, 'period_map_study.py'),
                    '--out', run_dir('reconstruction_K'), threads=1, minutes=10))
    J.append(script('offline_checks', os.path.join(EXP, 'offline_checks.py'),
                    '--out', run_dir('offline_checks'), '--acs-root', RUNS,
                    deps=[f'acs_{c}_s7' for c in CASES if c != 's4'], minutes=5))
    # ---- F. S^4 x T^2: reconstruction of the supplied period model (Figure 5, Table 6) and
    # the T^2-invariant reduction to S^4 (Appendix B): ACS + free-field CS, MLP conjugation, polynomial
    # generator. Commands as run by hand on 2026-10-04.
    J.append(script('s4t2_reconstruction', os.path.join(EXP, 's4t2_reconstruct.py'),
                    '--out', run_dir('s4t2_reconstruction'), threads=2, minutes=5))
    J.append(script('s4t2inv_s7', os.path.join(EXP, 's4t2_invariant_run.py'),
                    '--out', run_dir('s4t2inv_s7'), '--seed', 7, threads=1, minutes=30))
    J.append(script('s4t2invconj_s7', os.path.join(EXP, 's4t2_invariant_conj.py'),
                    '--out', run_dir('s4t2invconj_s7'), '--seed', 7, '--steps', 4000, threads=1, minutes=70))
    J.append(script('s4t2invconj_s7_cont', os.path.join(EXP, 's4t2_invariant_conj.py'),
                    '--out', run_dir('s4t2invconj_s7_cont'), '--seed', 7, '--steps', 8000, '--lr', 3e-4,
                    '--init-from', run_dir('s4t2invconj_s7'), deps=['s4t2invconj_s7'], threads=1, minutes=60))
    J.append(script('s4t2invcp_s7', os.path.join(EXP, 's4t2_invariant_conjpoly.py'),
                    '--out', run_dir('s4t2invcp_s7'), '--seed', 7, threads=1, minutes=100))
    # ---- G. supervised Hopf benchmark (Table 4): see SUPERVISED below
    J.extend(supervised_jobs())
    # ---- H. post-hoc analysis and figures: soft dependencies (run once every
    # dependency has finished, successfully or not)
    J.append(script('hopf_class', os.path.join(EXP, 'hopf_class.py'), '--runs', RUNS, '--out', run_dir('hopf_class'),
                    deps=[j['name'] for j in J if j['name'].startswith(('acs_s3xs1_', 'cspoly_s3xs1_'))]
                    + ['csconjexp_s3xs1', 'csconjaff_s3xs1'], minutes=25, soft=True))
    cs_jobs = [j['name'] for j in J if j['name'].startswith(('cspoly_', 'csconj', 'conjpoly_'))]
    J.append(script('cs_analysis', os.path.join(EXP, 'cs_analysis.py'),
                    '--runs', RUNS, '--out', run_dir('cs_analysis'), deps=cs_jobs, threads=2, minutes=30, soft=True))
    J.append(script('figures', os.path.join(EXP, 'make_figures.py'),
                    deps=[j['name'] for j in J], minutes=2, soft=True,
                    marker=os.path.join(OUT, 'figures_done.txt')))
    return J


def done(job):
    return os.path.exists(job['marker'])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--budget', type=int, default=9, help='total threads in use at any time')
    ap.add_argument('--dry', action='store_true')
    ap.add_argument('--only', default=None, help='comma-separated job-name prefixes')
    args = ap.parse_args()
    os.makedirs(RUNS, exist_ok=True); os.makedirs(LOGS, exist_ok=True)
    J = jobs()
    if args.only:
        pre = args.only.split(',')
        J = [j for j in J if any(j['name'].startswith(p) for p in pre)]
    byname = {j['name']: j for j in J}
    # priority = critical-path length: own minutes + the longest chain of dependants
    children = {j['name']: [] for j in J}
    for j in J:
        for d in j['deps']:
            if d in children:
                children[d].append(j['name'])
    crit = {}

    def cp(n):
        if n not in crit:
            crit[n] = byname[n]['minutes'] + max((cp(c) for c in children[n] if not byname[c].get('soft')), default=0)
        return crit[n]
    for j in J:
        cp(j['name'])
    if args.dry:
        for j in J:
            print(f"{j['name']:28s} deps={j['deps']} threads={j['threads']} ~{j['minutes']}min done={done(j)}")
        print(f"{len(J)} jobs, ~{sum(j['minutes'] * j['threads'] for j in J) / 60:.1f} thread-hours")
        return
    running, failed = {}, set()
    t0 = time.time()
    status_path = os.path.join(OUT, 'scheduler_status.json')
    while True:
        for name, (p, j, ts) in list(running.items()):
            rc = p.poll()
            if rc is not None:
                del running[name]
                ok = rc == 0 and done(j)
                if not ok:
                    failed.add(name)
                print(f"[{(time.time()-t0)/60:6.1f} min] {'DONE' if ok else 'FAILED'} {name} "
                      f"({(time.time()-ts)/60:.1f} min, rc={rc})", flush=True)
        used = sum(j['threads'] for _, j, _ in running.values())
        pending = [j for j in J if not done(j) and j['name'] not in running and j['name'] not in failed]
        blocked = [j for j in pending if not j.get('soft') and any(d in failed for d in j['deps'])]
        for j in blocked:
            failed.add(j['name'])
            print(f"SKIP {j['name']} (dependency failed)", flush=True)
        alljobs = {j['name']: j for j in jobs()}           # deps outside --only: judged by their markers
        ready = [j for j in pending if j not in blocked and
                 all(done(alljobs[d]) or (j.get('soft') and d in failed) for d in j['deps'] if d in alljobs)]
        ready.sort(key=lambda j: (-j.get('prio', 0), -crit[j['name']]))
        for j in ready:
            if used + j['threads'] > args.budget and running:
                continue
            env = dict(os.environ, OMP_NUM_THREADS=str(j['threads']), MKL_NUM_THREADS=str(j['threads']),
                       VECLIB_MAXIMUM_THREADS=str(j['threads']), PYTHONUNBUFFERED='1')
            log = open(os.path.join(LOGS, j['name'] + '.txt'), 'w')
            p = subprocess.Popen(j['cmd'], stdout=log, stderr=subprocess.STDOUT, env=env, cwd=ROOT)
            running[j['name']] = (p, j, time.time())
            used += j['threads']
            print(f"[{(time.time()-t0)/60:6.1f} min] START {j['name']} (threads {j['threads']}, ~{j['minutes']} min)",
                  flush=True)
        with open(status_path, 'w') as f:
            json.dump({'elapsed_min': (time.time() - t0) / 60, 'running': sorted(running),
                       'failed': sorted(failed), 'done': sorted(j['name'] for j in J if done(j)),
                       'pending': sorted(j['name'] for j in J if not done(j) and j['name'] not in running)}, f, indent=1)
        if not running and not ready:
            break
        time.sleep(10)
    print(f"finished in {(time.time()-t0)/60:.1f} min; failed: {sorted(failed)}", flush=True)


if __name__ == '__main__':
    main()
