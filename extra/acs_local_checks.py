"""Bounded CPU learning checks, not a substitute for the HPC pilot campaign.

Usage: python -m extra.acs_local_checks --out scratch/acs_local_checks/controls
At most 5,616 total Adam updates including target smokes; no analytic labels.
"""
import argparse
from pathlib import Path

from dataclasses import replace

from extra.acs_training import Config, Runner, atomic_json, read_config


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    results = []
    for case in ('s2', 's1xs1', 's3xs1', 's6', 's4xs2', 's4'):
        feature = {'s2': 'lin', 's1xs1': 'bilin'}.get(case, 'all2')
        cfg = Config(case=case, feat=feature, hidden=(16,),
                     poly_steps=2000 if case in ('s2', 's1xs1') else 250,
                     refine_steps=0, batch=64, pool_size=128, validation_size=256,
                     eval_every=500, audit_every=1000, audit_starts=16, audit_steps=10,
                     boundary_starts=32, boundary_steps=20, lbfgs_steps=0,
                     wall_minutes=5)
        result = Runner(cfg, Path(args.out) / case).run()
        results.append(result)
    atomic_json(results, Path(args.out) / 'controls.json')
    for result in results:
        print(result['case'], result['best_raw'], flush=True)
    for result in results[:2]:
        if result['best_raw']['max_found'] >= .05:
            raise SystemExit(f"Positive learning control {result['case']} did not reach residual < 0.05")
    root = Path(__file__).resolve().parents[1]
    configs = root / 'configs' / 'ambient' / 's4xs1xs1_acs'
    for seed in (7, 8, 9):
        cfg = read_config(configs / 'smoke.yaml')
        cfg.seed = seed
        Runner(cfg, Path(args.out) / f'target_smoke_s{seed}').run()
    # Exercise the actual 4x256 MLP widths and all four objectives, without
    # running the 10,000-step pilots locally.
    for path in sorted(configs.glob('pilot_*.yaml')):
        cfg = replace(read_config(path), poly_steps=2, refine_steps=2, batch=12,
                      pool_size=24, validation_size=20, replay_size=8, replay_starts=4,
                      replay_ascent_steps=2, boundary_starts=4, boundary_steps=2,
                      lbfgs_steps=0, wall_minutes=5)
        Runner(cfg, Path(args.out) / (path.stem + '_wiring')).run()


if __name__ == '__main__':
    main()
