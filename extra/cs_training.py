"""Checked, resumable CS continuation with fixed validation and best selection.

Reuses the ambient conjugation and Nijenhuis definitions, but not the legacy
trainer's automatic resume. Finite experiments never certify integrability.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import signal
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import torch

from extra.acs_training import atomic_json, atomic_torch_save, sha256
from extra.cs_pilot import checked_base, read_config
from extra.campaign_unit import (Geometry, Net, load_acs_field, make_field,
                                 polar_corrected, nij_norm, adversarial_sup)

DTYPE = torch.float64
ACS_LIMIT = 1e-10
CONDITION_LIMIT = 1e6


class Runner:
    def __init__(self, cfg, base_dir, audit_path, out, init_from=None, resume=False,
                 expected_source_sha256=None):
        self.cfg = cfg
        self.out = Path(out)
        self.base = checked_base(base_dir, audit_path)
        if resume == bool(init_from):
            raise ValueError('Choose exactly one of init_from and resume')
        if resume and expected_source_sha256:
            raise ValueError('Source hash applies only to a new continuation')
        if not resume and self.out.exists() and any(self.out.iterdir()):
            raise ValueError('New continuation requires an empty output directory')
        if resume and (self.out / 'final.pt').exists():
            raise ValueError('Already completed; continue into a new directory')
        torch.set_num_threads(1)
        self.geo = Geometry('s4xs1xs1', 'ambient')
        with torch.random.fork_rng():
            torch.manual_seed(cfg['seed'])
            self.model = Net(9, 9, cfg['stages'], cfg['hidden'], .5).to(DTYPE)
            self.base_field = polar_corrected(self.geo, load_acs_field(self.geo, base_dir))
        self.field = make_field(self.geo, 'cs', self.model, base_fn=self.base_field)[0]
        self.train_rng = torch.Generator().manual_seed(cfg['seed'] + 10000)
        validation_rng = torch.Generator().manual_seed(cfg['seed'] + 20000)
        self.validation = self.geo.manifold.sample(cfg['n_validation'], generator=validation_rng)
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=cfg['lr'])
        self.scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
            self.optimizer, T_max=cfg['steps'], eta_min=1e-5)
        self.step = 0
        self.best = self.initial = None
        self.previous_elapsed = 0.0
        self.stop_requested = False
        self.started = time.monotonic()
        if resume:
            state = torch.load(self.out / 'state.pt', map_location='cpu', weights_only=False)
            previous = state['config'].copy()
            expected = cfg.copy()
            previous.pop('wall_minutes'); expected.pop('wall_minutes')
            if state.get('schema_version') != 1 or previous != expected:
                raise ValueError('Resume configuration differs from saved state')
            if state['base']['sha256'] != self.base['sha256'] or state['base']['audit_sha256'] != self.base['audit_sha256']:
                raise ValueError('Resume base/audit differs from saved state')
            self.model.load_state_dict(state['model'])
            self.optimizer.load_state_dict(state['optimizer'])
            self.scheduler.load_state_dict(state['scheduler'])
            self.train_rng.set_state(state['train_rng'])
            self.validation = state['validation']
            self.step, self.best, self.initial = state['step'], state['best'], state['initial']
            self.source = state['source']
            self.previous_elapsed = state['elapsed_seconds']
        else:
            path = Path(init_from)
            path = path / 'final.pt' if path.is_dir() else path
            digest = sha256(path)
            if expected_source_sha256 and digest != expected_source_sha256:
                raise ValueError('CS source checksum differs from requested SHA256')
            ck = torch.load(path, map_location='cpu', weights_only=False)
            for key, value in {'case': 's4xs1xs1', 'task': 'cs', 'kind': 'ambient',
                               'model_type': 'mlp', 'hidden': cfg['hidden'], 'stages': cfg['stages']}.items():
                if ck.get(key) != value:
                    raise ValueError(f'CS source architecture mismatch: {key}')
            if ck.get('audited_base', {}).get('sha256') != self.base['sha256']:
                raise ValueError('CS source uses a different or unrecorded ACS base')
            self.model.load_state_dict(ck['model'])
            if sha256(path) != digest:
                raise ValueError('CS source changed during loading')
            self.source = {'path': str(path.resolve()), 'sha256': digest,
                           'optimizer_restarted': True}
        self.out.mkdir(parents=True, exist_ok=True)
        if not resume:
            atomic_json({'config': cfg, 'base': self.base, 'source': self.source,
                         'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
                         'dirty_worktree': bool(subprocess.check_output(['git', 'status', '--porcelain'], text=True).strip()),
                         'torch': torch.__version__, 'acs_limit': ACS_LIMIT,
                         'condition_limit': CONDITION_LIMIT, 'lr_min': 1e-5,
                         'note': 'Fixed validation for selection, not independent certification.'},
                        self.out / 'manifest.json')

    def payload(self):
        return {'case': 's4xs1xs1', 'task': 'cs', 'kind': 'ambient', 'model_type': 'mlp',
                'hidden': self.cfg['hidden'], 'stages': self.cfg['stages'],
                'seed': self.cfg['seed'], 'model': copy.deepcopy(self.model.state_dict()),
                'audited_base': self.base, 'source': self.source, 'step': self.step}

    def evaluate(self, points):
        norms, acs, conditions, jnorms = [], [], [], []
        for x in points.split(self.cfg['eval_chunk']):
            with torch.no_grad():
                P = self.geo.manifold.projector(x)
                J = self.field(x)
                acs.append((J @ J + P).norm(dim=(1, 2)))
                jnorms.append(J.norm(dim=(1, 2)))
                raw = self.model(x)
                G = torch.eye(9, dtype=DTYPE).expand(len(x), 9, 9)
                for stage in range(self.cfg['stages']):
                    G = torch.matrix_exp(P @ raw[:, stage] @ P) @ G
                singular = torch.linalg.svdvals(G)
                conditions.append(singular[:, 0] / singular[:, -1])
                norms.append(nij_norm(self.geo, self.field, x))
        n, a, c, j = map(torch.cat, (norms, acs, conditions, jnorms))
        if not all(torch.isfinite(t).all() for t in (n, a, c, j)):
            raise FloatingPointError('Nonfinite CS validation diagnostics')
        return {'nij_mean': n.mean().item(), 'nij_rms': n.square().mean().sqrt().item(),
                'nij_sup': n.max().item(), 'acs_floorfree_sup': a.max().item(),
                'acs_mean': a.mean().item(), 'acs_sup': a.max().item(),
                'conjugation_condition_max': c.max().item(), 'J_frobenius_max': j.max().item()}

    def observe(self):
        metrics = self.evaluate(self.validation)
        if self.initial is None:
            self.initial = metrics
        admissible = metrics['acs_floorfree_sup'] <= ACS_LIMIT and metrics['conjugation_condition_max'] <= CONDITION_LIMIT
        if admissible and (self.best is None or metrics['nij_rms'] < self.best['metrics']['nij_rms']):
            self.best = {**self.payload(), 'metrics': metrics}
        with (self.out / 'progress.jsonl').open('a') as stream:
            stream.write(json.dumps({'step': self.step, 'validation': metrics,
                                     'elapsed_min': (self.previous_elapsed + time.monotonic() - self.started) / 60,
                                     'lr': self.optimizer.param_groups[0]['lr'],
                                     'eligible': admissible, 'best_step': self.best['step'] if self.best else None}) + '\n')
        print(f"step {self.step}: {metrics}", flush=True)

    def save(self):
        atomic_torch_save({**self.payload(), 'schema_version': 1, 'config': self.cfg,
                          'base': self.base, 'optimizer': self.optimizer.state_dict(),
                          'scheduler': self.scheduler.state_dict(), 'train_rng': self.train_rng.get_state(),
                          'validation': self.validation, 'initial': self.initial, 'best': self.best,
                          'elapsed_seconds': self.previous_elapsed + time.monotonic() - self.started},
                         self.out / 'state.pt')
        if self.best:
            atomic_torch_save(self.best, self.out / 'best.pt')

    def run(self, max_updates=None):
        updates = 0
        if self.initial is None:
            self.observe()
            self.save()
        while self.step < self.cfg['steps']:
            if self.stop_requested or time.monotonic() - self.started >= 60*self.cfg['wall_minutes']:
                break
            x = self.geo.manifold.sample(self.cfg['batch'], generator=self.train_rng)
            self.optimizer.zero_grad()
            loss = nij_norm(self.geo, self.field, x).square().mean()
            if not torch.isfinite(loss):
                raise FloatingPointError('Nonfinite CS loss; last atomic state remains available')
            loss.backward()
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), 1.0, error_if_nonfinite=True)
            self.optimizer.step(); self.scheduler.step()
            self.step += 1
            updates += 1
            if self.step % self.cfg['eval_every'] == 0 or self.step == self.cfg['steps']:
                self.observe()
            if self.step % self.cfg['ckpt_every'] == 0:
                self.save()
            if max_updates is not None and updates >= max_updates:
                break
        self.save()  # last state, not the selected best; resume follows the same trajectory
        report = {'case': 's4xs1xs1', 'task': 'cs', 'formulation': 'ambient', 'seed': self.cfg['seed'],
                  'hidden': self.cfg['hidden'], 'stages': self.cfg['stages'], 'model_type': 'mlp',
                  'counts': {'adam': self.step},
                  'status': 'incomplete', 'steps_done': self.step, 'steps_planned': self.cfg['steps'],
                  'source': self.source, 'audited_base': self.base, 'initial_validation': self.initial,
                  'best_validation': self.best['metrics'] if self.best else None,
                  'selected_step': self.best['step'] if self.best else None,
                  'stopped_by_wall': self.step < self.cfg['steps'], 'numerically_validated_cs': False,
                  'exact_acs_by_construction': True,
                  'elapsed_min': (self.previous_elapsed + time.monotonic() - self.started) / 60}
        if self.step == self.cfg['steps']:
            if self.best is None:
                raise RuntimeError('No candidate passed algebraic/conditioning selection gates')
            self.model.load_state_dict(self.best['model'])
            rng = torch.Generator().manual_seed(self.cfg['seed'] + 30000)
            report['dense'] = self.evaluate(self.geo.manifold.sample(self.cfg['n_dense'], generator=rng))
            with torch.random.fork_rng():
                torch.manual_seed(self.cfg['seed'] + 40000)
                report['adv_nij_sup'], _ = adversarial_sup(
                    self.geo, lambda x: nij_norm(self.geo, self.field, x),
                    self.cfg['adv_starts'], self.cfg['adv_steps'])
            if (sha256(self.base['checkpoint']) != self.base['sha256']
                    or sha256(self.base['audit']) != self.base['audit_sha256']):
                raise RuntimeError('ACS base/audit changed during training')
            atomic_torch_save(self.best, self.out / 'final.pt')
            report['status'] = 'completed_pending_independent_audit'
        report['elapsed_min'] = (self.previous_elapsed + time.monotonic() - self.started) / 60
        atomic_json(report, self.out / 'result.json')
        return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    parser.add_argument('--base-from', required=True)
    parser.add_argument('--base-audit', required=True)
    parser.add_argument('--out', required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--init-from')
    mode.add_argument('--resume', action='store_true')
    parser.add_argument('--expected-source-sha256')
    args = parser.parse_args()
    runner = Runner(read_config(args.config), args.base_from, args.base_audit, args.out,
                    args.init_from, args.resume, args.expected_source_sha256)
    old = {sig: signal.signal(sig, lambda *_: setattr(runner, 'stop_requested', True))
           for sig in (signal.SIGTERM, signal.SIGINT)}
    try:
        report = runner.run()
    finally:
        for sig, handler in old.items():
            signal.signal(sig, handler)
    print(f"{report['status']}: {runner.out}", flush=True)
    if report['status'] == 'incomplete':
        raise SystemExit(75)


if __name__ == '__main__':
    main()
