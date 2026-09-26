"""Run the successful ACS recipe with explicit best-checkpoint restarts.

This replays an experimental protocol, not a guarantee of convergence or bitwise
equivalence across PyTorch/platform versions. Every stage keeps its own logs.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import signal
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import yaml

from extra.acs_training import Config, Runner, atomic_json, sha256


def read_recipe(path, seed=None):
    with open(path) as stream:
        raw = yaml.safe_load(stream)
    keys = {'seed', 'common', 'initial', 'polish', 'polish_seed_offsets'}
    if not isinstance(raw, dict) or set(raw) != keys:
        raise ValueError(f'Recipe must contain exactly {sorted(keys)}')
    seed = raw['seed'] if seed is None else seed
    offsets = raw['polish_seed_offsets']
    if not isinstance(offsets, list) or any(type(v) is not int or v <= 0 for v in offsets):
        raise ValueError('polish_seed_offsets must be a list of positive integers')
    if len(set(offsets)) != len(offsets):
        raise ValueError('polish seeds must be distinct')
    for name in ('common', 'initial', 'polish'):
        if not isinstance(raw[name], dict) or 'seed' in raw[name]:
            raise ValueError(f'{name} must be a mapping without a seed override')
    configs = [Config(**{**raw['common'], **raw['initial'], 'seed': seed}).validate()]
    configs += [Config(**{**raw['common'], **raw['polish'], 'seed': seed + offset}).validate()
                for offset in offsets]
    first = configs[0]
    if first.poly_steps <= 0:
        raise ValueError('Initial stage must learn the polynomial from scratch')
    for cfg in configs[1:]:
        if cfg.poly_steps or (cfg.case, cfg.feat, cfg.hidden) != (first.case, first.feat, first.hidden):
            raise ValueError('Polishing must preserve case/architecture and skip polynomial training')
    return configs


def run_recipe(configs, out, resume=False, wall_minutes=180.0):
    if not 0 < wall_minutes < float('inf'):
        raise ValueError('wall_minutes must be positive and finite')
    out = Path(out)
    # Normalize tuples for exact JSON comparison. Duration is per allocation.
    specification = json.loads(json.dumps([asdict(cfg) for cfg in configs]))
    manifest = out / 'recipe_manifest.json'
    if resume:
        if not manifest.is_file() or json.loads(manifest.read_text())['stages'] != specification:
            raise ValueError('Resume requires the same recipe and seed')
    else:
        if out.exists() and any(out.iterdir()):
            raise ValueError('Output directory must be empty; use --resume for this recipe')
        out.mkdir(parents=True, exist_ok=True)
        atomic_json({'stages': specification, 'note': 'Numerical recipe; independent audit required.'}, manifest)
    ledger_path = out / 'recipe_result.json'
    ledger = json.loads(ledger_path.read_text()) if resume and ledger_path.exists() else {'stages': []}
    started = time.monotonic()
    source = None
    completed = []
    for index, original in enumerate(configs):
        directory = out / f'{index:02d}_s{original.seed}'
        final = directory / 'final.pt'
        if final.exists():
            result = json.loads((directory / 'result.json').read_text())
            if result['status'] != 'completed_pending_independent_audit':
                raise ValueError('Stage final.pt exists without completed result')
            digest = sha256(final)
            if index < len(ledger['stages']) and ledger['stages'][index]['sha256'] != digest:
                raise ValueError('Completed stage checkpoint changed since the recipe recorded it')
            # A crash after stage completion but before updating the ledger is recoverable.
        else:
            remaining = wall_minutes - (time.monotonic() - started) / 60
            if remaining <= 0:
                break
            cfg = Config(**{**asdict(original), 'wall_minutes': remaining})
            continuing = (directory / 'state.pt').exists()
            runner = Runner(cfg, directory, resume=continuing,
                            init_from=source if not continuing else None,
                            expected_source_sha256=sha256(source) if source and not continuing else None)
            old = {sig: signal.signal(sig, lambda *_: setattr(runner, 'stop_requested', True))
                   for sig in (signal.SIGTERM, signal.SIGINT)}
            try:
                result = runner.run()
            finally:
                for sig, handler in old.items():
                    signal.signal(sig, handler)
            if result['status'] == 'incomplete':
                break
            digest = sha256(final)
        source = final
        completed.append({'directory': str(directory.resolve()), 'sha256': digest,
                          'seed': original.seed, 'counts': result['counts']})
        atomic_json({'status': 'incomplete', 'stages': completed}, ledger_path)
    done = len(completed) == len(configs)
    report = {'status': 'completed_pending_independent_audit' if done else 'incomplete',
              'stages': completed, 'final_checkpoint': str(source.resolve()) if done else None,
              'numerically_validated': False}
    atomic_json(report, ledger_path)
    print(json.dumps(report, indent=2), flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--seed', type=int)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--wall-minutes', type=float, default=180.0)
    args = parser.parse_args()
    report = run_recipe(read_recipe(args.config, args.seed), args.out, args.resume, args.wall_minutes)
    if report['status'] == 'incomplete':
        raise SystemExit(75)


if __name__ == '__main__':
    main()
