"""Allocation-resumable long CS search with separate mean/tail candidates.

Selection witnesses never enter training replay. Finite diagnostics are not
integrability certificates. Stop thresholds are operational, not obstructions.
"""
from __future__ import annotations

import argparse
import copy
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import signal
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import torch
import yaml

from extra.acs_training import Stopped, atomic_json, atomic_torch_save, sha256
from extra.cs_matched import Runner as ReplayRunner, Sampling
from extra.cs_pilot import validate_config
from extra.cs_training import ACS_LIMIT, CONDITION_LIMIT
from extra.campaign_unit import nij_norm
from extra.diagnose_cs import ascend, values


@dataclass
class Controls:
    monitor_every: int = 2500
    monitor_pool: int = 2048
    monitor_starts: int = 64
    monitor_steps: int = 60
    plateau_window: int = 10000
    minimum_relative_gain: float = .005
    norm_growth_limit: float = 2.
    condition_growth_limit: float = 3.
    candidate_regression_limit: float = 1.05

    def validate(self, eval_every):
        for key, value in asdict(self).items():
            if key in ('monitor_every', 'monitor_pool', 'monitor_starts', 'monitor_steps', 'plateau_window'):
                if type(value) is not int or value < 1:
                    raise ValueError(f'Invalid {key}')
            elif type(value) not in (int, float) or not 0 < value < float('inf'):
                raise ValueError(f'Invalid {key}')
        if (self.monitor_pool < self.monitor_starts or self.monitor_every % eval_every
                or self.plateau_window % self.monitor_every or self.minimum_relative_gain >= 1
                or min(self.norm_growth_limit, self.condition_growth_limit, self.candidate_regression_limit) < 1):
            raise ValueError('Inconsistent long-run monitoring controls')
        return self


def read_config(path):
    with open(path) as stream:
        raw = yaml.safe_load(stream)
    if not isinstance(raw, dict) or set(raw) != {'training', 'sampling', 'long_run'}:
        raise ValueError('Expected training, sampling, long_run sections')
    cfg = validate_config(raw['training'])
    if cfg['lr'] != 1e-5:
        raise ValueError('Long recipe requires constant 1e-5 LR (equal to scheduler floor)')
    sampling = Sampling(**raw['sampling']).validate(cfg['batch'])
    control = Controls(**raw['long_run']).validate(cfg['eval_every'])
    return {**cfg, 'matched': {'mode': 'replay', 'sampling': asdict(sampling)}, 'long_run': asdict(control)}


class Runner(ReplayRunner):
    def __init__(self, cfg, *args, **kwargs):
        self.controls = Controls(**cfg['long_run']).validate(cfg['eval_every'])
        self.tail_best = None
        self.monitor_points = torch.empty((0, 9), dtype=torch.float64)
        self.monitor_rng = torch.Generator().manual_seed(cfg['seed'] + 70000)
        self.last_monitor = 0
        self.pending_observation = False
        self.history = []
        self.terminal_reason = None
        super().__init__(cfg, *args, **kwargs)
        if kwargs.get('resume', False):
            extra = torch.load(self.out / 'state.pt', map_location='cpu', weights_only=False)['training_extension']['long']
            for key in ('tail_best', 'monitor_points', 'last_monitor', 'pending_observation', 'history', 'terminal_reason'):
                setattr(self, key, extra[key])
            self.monitor_rng.set_state(extra['monitor_rng'])
        else:
            manifest = json.loads((self.out / 'manifest.json').read_text())
            manifest.update(selection='Separate minimum uniform RMS and minimum monitored tail; source-relative regression gates',
                            note='Growing selection-only witness set; final diagnostics remain separate. No integrability certificate.')
            atomic_json(manifest, self.out / 'manifest.json')

    def training_state(self):
        return {**super().training_state(), 'long': {
            'tail_best': self.tail_best, 'monitor_points': self.monitor_points,
            'monitor_rng': self.monitor_rng.get_state(), 'last_monitor': self.last_monitor,
            'pending_observation': self.pending_observation, 'history': self.history,
            'terminal_reason': self.terminal_reason}}

    def measure(self):
        metrics = self.evaluate(self.validation)
        metrics['challenge'] = self.evaluate(self.challenge)
        groups = [metrics, metrics['challenge']]
        if len(self.monitor_points):
            metrics['monitor'] = self.evaluate(self.monitor_points)
            groups.append(metrics['monitor'])
        metrics['tail_max'] = max(m['nij_sup'] for m in groups)
        metrics['numerically_eligible'] = all(
            m['acs_floorfree_sup'] <= ACS_LIMIT and m['tangency_max'] <= 1e-10
            and m['conjugation_condition_max'] <= CONDITION_LIMIT for m in groups)
        return metrics

    def refresh_monitor(self):
        rng = torch.Generator().set_state(self.monitor_rng.get_state())
        pool = self.geo.manifold.sample(self.controls.monitor_pool, generator=rng)
        fn = lambda x: nij_norm(self.geo, self.field, x)
        scores = values(fn, pool, self.cfg['eval_chunk'], self.check_stop)
        points, _ = ascend(self.geo, fn, pool[scores.topk(self.controls.monitor_starts).indices],
                           self.controls.monitor_steps, self.cfg['eval_chunk'], self.check_stop)
        self.check_stop()
        self.monitor_points = torch.cat((self.monitor_points, points))
        self.monitor_rng.set_state(rng.get_state())
        self.last_monitor = self.step
        # Re-evaluate incumbents on the expanded set so candidate tail scores use
        # exactly the same points. Uniform sets and model weights never change here.
        current = copy.deepcopy(self.model.state_dict())
        try:
            for candidate in (self.best, self.tail_best):
                if candidate is not None:
                    self.model.load_state_dict(candidate['model'])
                    candidate['metrics'] = self.measure()
        finally:
            self.model.load_state_dict(current)

    def observe(self):
        if self.challenge is None:
            self.prepare_challenge()
        if self.step > self.last_monitor and self.step % self.controls.monitor_every == 0:
            self.refresh_monitor()
        metrics = self.measure()
        if self.initial is None:
            self.initial = copy.deepcopy(metrics)
        # Compare growth on the same fixed uniform/challenge points, never on a
        # changing sample maximum. These are safety stops, not geometry constraints.
        c = self.controls
        fixed_pairs = [(metrics, self.initial), (metrics['challenge'], self.initial['challenge'])]
        if not metrics['numerically_eligible']:
            self.terminal_reason = 'numerical_guard'
        elif any(m['J_frobenius_max'] > c.norm_growth_limit*b['J_frobenius_max']
                 or m['conjugation_condition_max'] > c.condition_growth_limit*b['conjugation_condition_max']
                 for m, b in fixed_pairs):
            self.terminal_reason = 'growth_guard'
        eligible = metrics['numerically_eligible'] and self.terminal_reason is None
        # Gate against the original source's measured errors. No claim that these
        # source values bound the field globally or that gates are mathematically necessary.
        candidates = [ck for ck in (self.best, self.tail_best) if ck is not None]
        if eligible:
            candidates.append({**self.payload(), 'metrics': copy.deepcopy(metrics)})
        valid = [ck for ck in candidates if ck['metrics']['numerically_eligible']]
        average_choices = [ck for ck in valid if ck['metrics']['tail_max'] <= c.candidate_regression_limit*self.initial['tail_max']]
        tail_choices = [ck for ck in valid if ck['metrics']['nij_rms'] <= c.candidate_regression_limit*self.initial['nij_rms']]
        # An incumbent can change rank or become ineligible when new monitoring
        # points are added. Compare all retained candidates on that same expanded set.
        if average_choices:
            self.best = min(average_choices, key=lambda ck: ck['metrics']['nij_rms'])
        if tail_choices:
            self.tail_best = min(tail_choices, key=lambda ck: ck['metrics']['tail_max'])
        if not average_choices or not tail_choices:
            self.terminal_reason = self.terminal_reason or 'candidate_gate_review'
        if self.step % c.monitor_every == 0 and self.best and self.tail_best:
            row = {'step': self.step, 'rms': self.best['metrics']['nij_rms'], 'tail': self.tail_best['metrics']['tail_max']}
            previous = next((h for h in self.history if h['step'] == self.step-c.plateau_window), None)
            if previous is not None and not self.terminal_reason:
                gains = [1-row[k]/max(previous[k], 1e-14) for k in ('rms', 'tail')]
                if max(gains) < c.minimum_relative_gain:
                    self.terminal_reason = 'plateau_review'
            self.history.append(row)
        self.pending_observation = False
        row = {'step': self.step, 'validation': metrics, 'lr': self.optimizer.param_groups[0]['lr'],
               'average_step': self.best['step'] if self.best else None,
               'tail_step': self.tail_best['step'] if self.tail_best else None,
               'monitor_count': len(self.monitor_points), 'terminal_reason': self.terminal_reason}
        with (self.out / 'progress.jsonl').open('a') as stream:
            stream.write(json.dumps(row) + '\n')
        print(f"step {self.step}: RMS={metrics['nij_rms']:.6g}, monitored tail={metrics['tail_max']:.6g}, "
              f"Jmax={metrics['J_frobenius_max']:.6g}, stop={self.terminal_reason}", flush=True)

    def run(self, max_updates=10000):
        if type(max_updates) is not int or max_updates < 1:
            raise ValueError('max_updates must be positive')
        updates = 0
        try:
            if self.initial is None or self.pending_observation:
                self.observe()
                self.save()
            while self.step < self.cfg['steps'] and self.terminal_reason is None and updates < max_updates:
                self.check_stop()
                x = self.training_batch()
                self.optimizer.zero_grad()
                loss = nij_norm(self.geo, self.field, x).square().mean()
                if not torch.isfinite(loss):
                    raise FloatingPointError('Nonfinite loss; resume only after review of last atomic state')
                loss.backward()
                torch.nn.utils.clip_grad_norm_(self.model.parameters(), 1., error_if_nonfinite=True)
                self.optimizer.step()
                self.scheduler.step()
                self.step += 1
                updates += 1
                if self.step % self.cfg['eval_every'] == 0 or self.step == self.cfg['steps']:
                    # Persist if a monitor search is interrupted: resume completes
                    # this observation before another update (same candidate selection).
                    self.pending_observation = True
                    self.observe()
                if self.step % self.cfg['ckpt_every'] == 0:
                    self.save()
        except Stopped:
            pass
        self.save()
        self.check_inputs()
        finished = self.step == self.cfg['steps'] and not self.pending_observation
        terminal = (finished or self.terminal_reason is not None) and self.best is not None and self.tail_best is not None
        if terminal:
            atomic_torch_save(self.best, self.out / 'final_average.pt')
            atomic_torch_save(self.tail_best, self.out / 'final_tail.pt')
        status = ('stopped_for_review' if self.terminal_reason else
                  'completed_pending_independent_audit' if finished else 'incomplete')
        report = {'case': 's4xs1xs1', 'task': 'cs', 'status': status, 'terminal_reason': self.terminal_reason,
                  'steps_done': self.step, 'steps_planned': self.cfg['steps'], 'updates_this_allocation': updates,
                  'pending_observation': self.pending_observation, 'config': self.cfg,
                  'source': self.source, 'audited_base': self.base, 'witness_input': self.witness_input,
                  'initial_validation': self.initial, 'monitor_count': len(self.monitor_points),
                  'elapsed_min': (self.previous_elapsed + time.monotonic()-self.started)/60,
                  'numerically_validated_cs': False, 'ready_for_diagnostics': terminal,
                  'candidates': {name: {'step': ck['step'], 'metrics': ck['metrics']} if ck else None
                                 for name, ck in [('average', self.best), ('tail', self.tail_best)]}}
        atomic_json(report, self.out / 'result.json')
        return report


def summary(root):
    root = Path(root)
    path = root / 'result.json'
    training = json.loads(path.read_text()) if path.exists() else {}
    report = {'training': training, 'diagnostics': {}, 'numerically_validated_cs': False}
    lines = ['Long CS continuation (not an integrability certificate)',
             f"Training: {training.get('status', 'MISSING')}; steps={training.get('steps_done')}/{training.get('steps_planned')}; stop={training.get('terminal_reason')}"]
    for name in ('average', 'tail'):
        path = root / f'diagnostic_{name}/report.json'
        diagnostic = json.loads(path.read_text()) if path.exists() else {}
        checkpoint = root / f'final_{name}.pt'
        matches = checkpoint.is_file() and diagnostic.get('sha256') == sha256(checkpoint)
        report['diagnostics'][name] = {'report': diagnostic, 'matches_checkpoint': bool(matches)}
        lines += [f"{name}: selected={training.get('candidates', {}).get(name)}",
                  f"  diagnostic={diagnostic.get('status', 'MISSING')}; checkpoint matches={bool(matches)}",
                  f"  fresh uniform={diagnostic.get('uniform')}",
                  f"  worst found={diagnostic.get('combined_completed_groups', {}).get('nij_max_found')}"]
    report['diagnostics_complete'] = bool(training.get('ready_for_diagnostics')) and all(
        d['matches_checkpoint'] and d['report'].get('status') == 'completed_diagnostic' for d in report['diagnostics'].values())
    lines.append(f"Both diagnostics complete: {report['diagnostics_complete']}. No CS certification.")
    atomic_json(report, root / 'comparison.json')
    (root / 'comparison.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('summary').add_argument('--root', required=True)
    train = sub.add_parser('train')
    for name in ('config', 'base-from', 'base-audit', 'out', 'witness-dir', 'init-from', 'expected-source-sha256'):
        train.add_argument('--' + name, required=True)
    train.add_argument('--allocation', type=int, required=True)
    train.add_argument('--max-updates', type=int, default=10000)
    args = parser.parse_args()
    if args.command == 'summary':
        summary(args.root)
        return
    if args.allocation < 1:
        parser.error('allocation must be positive')
    out = Path(args.out)
    resume = args.allocation > 1
    if sha256(args.init_from) != args.expected_source_sha256:
        parser.error('Source checkpoint hash mismatch')
    if resume:
        receipt_path = out / f'allocation_{args.allocation-1}.json'
        if not receipt_path.exists():
            parser.error('Previous allocation did not finish safely; inspect its log before resubmitting')
        previous_receipt = json.loads(receipt_path.read_text())
        if (previous_receipt.get('source_sha256') != args.expected_source_sha256
                or previous_receipt.get('state_sha256') != sha256(out / 'state.pt')):
            parser.error('Allocation chain source/state mismatch')
    def receipt():
        atomic_json({'allocation': args.allocation, 'source_sha256': args.expected_source_sha256,
                     'state_sha256': sha256(out / 'state.pt')}, out / f'allocation_{args.allocation}.json')
    # Later chained slots become quick no-ops after a deliberate early stop or completion.
    if resume and (out / 'result.json').exists():
        result = json.loads((out / 'result.json').read_text())
        if result.get('ready_for_diagnostics'):
            print('Training already terminal; retaining selected candidates.', flush=True)
            receipt()
            return
    runner = Runner(read_config(args.config), args.base_from, args.base_audit, out, args.witness_dir,
                    init_from=None if resume else args.init_from, resume=resume,
                    expected_source_sha256=None if resume else args.expected_source_sha256)
    old = {sig: signal.signal(sig, lambda *_: setattr(runner, 'stop_requested', True))
           for sig in (signal.SIGTERM, signal.SIGINT)}
    try:
        report = runner.run(args.max_updates)
    finally:
        for sig, handler in old.items():
            signal.signal(sig, handler)
    # A productive, atomically saved allocation may continue in the next chained job.
    # Signals and zero-progress stops fail rather than quietly consuming the whole chain.
    if runner.stop_requested or (not report['ready_for_diagnostics'] and not report['updates_this_allocation']):
        raise SystemExit(75)
    receipt()


if __name__ == '__main__':
    main()
