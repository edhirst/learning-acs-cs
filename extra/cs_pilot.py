"""Bounded CS pilot starting at a full-audited, checksum-matched learned ACS.

Wraps the existing ambient conjugation/Nijenhuis trainer. Deliberately no
automatic resume, long-run promotion, or claim of integrability. A CS-specific
independent audit and robust long-run restart are subsequent work.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import torch
import yaml

from extra.acs_training import atomic_json, sha256
from extra import campaign_unit


def checked_base(base_dir, audit_path):
    checkpoint = Path(base_dir) / 'final.pt'
    report = json.loads(Path(audit_path).read_text())
    digest = sha256(checkpoint)
    if (report.get('case') != 's4xs1xs1' or report.get('audit_tier') != 'full'
            or report.get('numerically_validated') is not True
            or report.get('screening_targets_passed') is not True):
        raise ValueError('CS pilot requires a passing full S4xS1xS1 ACS audit')
    if report.get('sha256') != digest:
        raise ValueError('ACS checkpoint SHA256 differs from the full audit')
    return {'checkpoint': str(checkpoint.resolve()), 'sha256': digest,
            'audit': str(Path(audit_path).resolve()), 'audit_sha256': sha256(audit_path)}


def read_config(path):
    with open(path) as stream:
        cfg = yaml.safe_load(stream)
    return validate_config(cfg)


def validate_config(cfg):
    keys = {'seed', 'steps', 'hidden', 'stages', 'batch', 'lr', 'wall_minutes',
            'eval_every', 'ckpt_every', 'n_dense', 'n_chart', 'adv_starts', 'adv_steps',
            'n_initial', 'n_validation', 'eval_chunk'}
    if not isinstance(cfg, dict) or set(cfg) != keys:
        raise ValueError(f'CS config must contain exactly {sorted(keys)}')
    for key in keys - {'hidden', 'lr', 'wall_minutes'}:
        if type(cfg[key]) is not int or cfg[key] < (0 if key == 'seed' else 1):
            raise ValueError(f'Invalid {key}')
    if cfg['adv_starts'] < 2 or cfg['adv_steps'] < 2:
        raise ValueError('Adversarial budgets must be >=2 (legacy CS report halves them)')
    if not isinstance(cfg['hidden'], list) or not cfg['hidden'] or any(type(h) is not int or h < 1 for h in cfg['hidden']):
        raise ValueError('hidden must contain positive integer widths')
    if not 1e-5 <= cfg['lr'] < float('inf') or not 0 < cfg['wall_minutes'] < float('inf'):
        raise ValueError('Invalid learning rate or wall duration')
    return cfg


def run(cfg, base_dir, audit_path, out):
    evidence = checked_base(base_dir, audit_path)
    out = Path(out)
    if out.exists() and any(out.iterdir()):
        raise ValueError('CS pilots require a new output directory; automatic resume is disabled')
    out.mkdir(parents=True, exist_ok=True)
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    atomic_json({'config': cfg, 'base': evidence, 'revision': revision,
                 'dirty_worktree': bool(subprocess.check_output(
                     ['git', 'status', '--porcelain'], text=True).strip()),
                 'torch': torch.__version__, 'identity_init': True, 'sampling': 'uniform',
                 'lbfgs': 0, 'note': 'Exploratory pilot, not a CS certificate.'},
                out / 'cs_pilot_manifest.json')
    argv = ['--case', 's4xs1xs1', '--formulation', 'ambient', '--task', 'cs',
            '--base-from', str(Path(base_dir).resolve()), '--out', str(out),
            '--identity-init', '--lbfgs', '0', '--poly-steps', '0']
    for key, value in cfg.items():
        if key == 'hidden':
            value = ','.join(map(str, value))
        argv.extend(['--' + key.replace('_', '-'), str(value)])
    torch.set_num_threads(1)
    campaign_unit.main(argv)
    if sha256(evidence['checkpoint']) != evidence['sha256']:
        raise RuntimeError('Base checkpoint changed during CS training; results are invalid')
    if sha256(evidence['audit']) != evidence['audit_sha256']:
        raise RuntimeError('Base audit changed during CS training; results are invalid')
    result_path = out / 'result.json'
    result = json.loads(result_path.read_text())
    result['audited_base'] = evidence
    result['numerically_validated_cs'] = False
    completed = (not result.get('stopped_by_wall', True)
                 and result.get('steps_done') == cfg['steps']
                 and result.get('nonfinite_batches', 0) == 0)
    result['status'] = 'pilot_completed_pending_review' if completed else 'incomplete'
    if not completed and (out / 'final.pt').exists():
        (out / 'final.pt').rename(out / 'partial.pt')
    atomic_json(result, result_path)
    if completed:
        ck = torch.load(out / 'final.pt', map_location='cpu', weights_only=False)
        ck['audited_base'] = evidence
        torch.save(ck, out / 'final.pt')
    print(f"CS {result['status']}: {result_path}", flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    parser.add_argument('--base-from', required=True)
    parser.add_argument('--base-audit', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--seed', type=int)
    parser.add_argument('--stages', type=int, choices=[1, 2])
    args = parser.parse_args()
    cfg = read_config(args.config)
    for key in ('seed', 'stages'):
        value = getattr(args, key)
        if value is not None:
            if value < 0:
                parser.error('seed must be nonnegative')
            cfg[key] = value
    if run(cfg, args.base_from, args.base_audit, args.out)['status'] == 'incomplete':
        raise SystemExit(75)


if __name__ == '__main__':
    main()
