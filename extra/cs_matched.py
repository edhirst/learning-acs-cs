"""Matched uniform/replay CS experiments; numerical search, never certification."""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import signal
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import torch
import yaml

from extra.acs_training import Stopped, atomic_json, sha256
from extra.cs_pilot import validate_config
from extra.cs_training import Runner as Continuation, ACS_LIMIT, CONDITION_LIMIT
from extra.campaign_unit import nij_norm, retract, tangential
from extra.diagnose_cs import ascend, values, point_metrics, summarize


@dataclass
class Sampling:
    target_count: int = 8
    capacity: int = 256
    refresh_every: int = 250
    pool_size: int = 512
    starts: int = 32
    ascent_steps: int = 20
    jitter: float = .05
    separation: float = .05
    challenge_pool: int = 2048
    challenge_starts: int = 64
    challenge_steps: int = 40

    def validate(self, batch):
        for key, value in asdict(self).items():
            if key in ('jitter', 'separation'):
                if type(value) not in (float, int) or not 0 < value < 1:
                    raise ValueError(f'Invalid {key}')
            elif type(value) is not int or value < 1:
                raise ValueError(f'Invalid {key}')
        if (self.target_count >= batch or self.pool_size < self.starts
                or self.challenge_pool < self.challenge_starts):
            raise ValueError('Invalid targeted batch or search sizes')
        return self


def read_config(path, mode):
    if mode not in ('uniform', 'replay'):
        raise ValueError('mode must be uniform or replay')
    with open(path) as stream:
        raw = yaml.safe_load(stream)
    if not isinstance(raw, dict) or set(raw) != {'training', 'sampling'}:
        raise ValueError('Expected training and sampling config sections')
    cfg = validate_config(raw['training'])
    sampling = Sampling(**raw['sampling']).validate(cfg['batch'])
    return {**cfg, 'matched': {'mode': mode, 'sampling': asdict(sampling)}}


def diverse(points, scores, capacity, separation):
    """Greedy residual-ranked diversity in ambient chord distance; no duplicate fill."""
    chosen = []
    for index in scores.argsort(descending=True).tolist():
        if not chosen or (points[chosen] - points[index]).norm(dim=1).min() >= separation:
            chosen.append(index)
            if len(chosen) == capacity:
                break
    if not chosen:
        raise ValueError('No replay witnesses')
    return points[chosen].detach().clone()


def load_witnesses(directory, source_hash):
    directory = Path(directory)
    report_path, points_path = directory / 'report.json', directory / 'diagnostic_points.pt'
    before = [sha256(report_path), sha256(points_path)]
    report = json.loads(report_path.read_text())
    data = torch.load(points_path, map_location='cpu', weights_only=False)
    if (report.get('status') != 'completed_diagnostic' or data.get('status') != 'completed_diagnostic'
            or report.get('sha256') != source_hash or data.get('checkpoint_sha256') != source_hash
            or report.get('points_sha256') != before[1]):
        raise ValueError('Witness diagnostic status/source/checksum mismatch')
    points, groups, scores = data['points'], data['groups'], data['metrics']['nij']
    if (not isinstance(points, torch.Tensor) or points.dtype != torch.float64
            or points.ndim != 2 or points.shape[1] != 9 or not len(points)
            or not isinstance(scores, torch.Tensor) or scores.shape != (len(points),)
            or len(groups) != len(points) or not all(isinstance(g, str) for g in groups)
            or not torch.isfinite(points).all() or not torch.isfinite(scores).all() or (scores < 0).any()):
        raise ValueError('Malformed witness arrays')
    for block in points.split([5, 2, 2], dim=1):
        if not torch.allclose(block.norm(dim=1), torch.ones(len(points), dtype=points.dtype), atol=1e-10, rtol=0):
            raise ValueError('Witness points are off the product of spheres')
    mask = torch.tensor([g.startswith('search_') for g in groups])
    if not mask.any():
        raise ValueError('No adversarial endpoints in witness file')
    if before != [sha256(report_path), sha256(points_path)]:
        raise ValueError('Witness inputs changed during loading')
    evidence = {'report': str(report_path.resolve()), 'report_sha256': before[0],
                'points': str(points_path.resolve()), 'points_sha256': before[1],
                'checkpoint_sha256': source_hash}
    return points[mask], scores[mask], evidence


class Runner(Continuation):
    def __init__(self, cfg, base_dir, audit_path, out, witness_dir, init_from=None,
                 resume=False, expected_source_sha256=None):
        self.mode = cfg['matched']['mode']
        if self.mode not in ('uniform', 'replay'):
            raise ValueError('Invalid matched mode')
        self.sampling = Sampling(**cfg['matched']['sampling']).validate(cfg['batch'])
        state = torch.load(Path(out) / 'state.pt', map_location='cpu', weights_only=False) if resume else None
        source = Path(state['source']['path'] if resume else init_from)
        source = source / 'final.pt' if source.is_dir() else source
        digest = sha256(source)
        if expected_source_sha256 and digest != expected_source_sha256:
            raise ValueError('CS source checksum mismatch')
        points, scores, self.witness_input = load_witnesses(witness_dir, digest)
        self.replay = diverse(points, scores, self.sampling.capacity, self.sampling.separation)
        self.replay_rng = torch.Generator().manual_seed(cfg['seed'] + 50000)
        self.challenge = None
        self.last_refresh = 0
        super().__init__(cfg, base_dir, audit_path, out, init_from, resume, expected_source_sha256)
        if state is not None:
            extra = state.get('training_extension')
            if not extra or extra['witness_input'] != self.witness_input or self.source['sha256'] != digest:
                raise ValueError('Resume witness/source provenance mismatch')
            self.replay = extra['replay']
            self.replay_rng.set_state(extra['replay_rng'])
            self.challenge = extra['challenge']
            self.last_refresh = extra['last_refresh']
        else:
            manifest = json.loads((self.out / 'manifest.json').read_text())
            manifest.update(witness_input=self.witness_input,
                            selection='minimax of relative uniform RMS and fixed-challenge max; each <= 1.05 of source',
                            note='Witness points are training data, not independent validation. No CS certificate.')
            atomic_json(manifest, self.out / 'manifest.json')

    def check_stop(self):
        if self.stop_requested or time.monotonic() - self.started >= 60*self.cfg['wall_minutes']:
            raise Stopped()

    def training_state(self):
        return {'witness_input': self.witness_input, 'replay': self.replay,
                'replay_rng': self.replay_rng.get_state(), 'challenge': self.challenge,
                'last_refresh': self.last_refresh}

    def prepare_challenge(self):
        # Identical fresh source-based points in both arms. Never used for gradients
        # of the training loss; repeated selection still makes these validation data.
        rng = torch.Generator().manual_seed(self.cfg['seed'] + 60000)
        pool = self.geo.manifold.sample(self.sampling.challenge_pool, generator=rng)
        fn = lambda x: nij_norm(self.geo, self.field, x)
        scores = values(fn, pool, self.cfg['eval_chunk'], self.check_stop)
        points, _ = ascend(self.geo, fn, pool[scores.topk(self.sampling.challenge_starts).indices],
                           self.sampling.challenge_steps, self.cfg['eval_chunk'], self.check_stop)
        self.challenge = points

    def refresh(self):
        # All changes are committed together. An interrupted search must not advance
        # either sampler or alter the buffer before the corresponding Adam update.
        rng = torch.Generator().set_state(self.replay_rng.get_state())
        pool = self.geo.manifold.sample(self.sampling.pool_size, generator=rng)
        candidates = torch.cat((self.replay, pool))
        fn = lambda x: nij_norm(self.geo, self.field, x)
        scores = values(fn, candidates, self.cfg['eval_chunk'], self.check_stop)
        starts = diverse(candidates, scores, self.sampling.starts, self.sampling.separation)
        endpoints, _ = ascend(self.geo, fn, starts, self.sampling.ascent_steps,
                              self.cfg['eval_chunk'], self.check_stop)
        # Keep old and new regions, ranked using the *current* field, with diversity.
        candidates = torch.cat((self.replay, endpoints, pool))
        scores = values(fn, candidates, self.cfg['eval_chunk'], self.check_stop)
        replay = diverse(candidates, scores, self.sampling.capacity, self.sampling.separation)
        self.check_stop()
        self.replay = replay
        self.replay_rng.set_state(rng.get_state())
        self.last_refresh = self.step
        print(f'refreshed replay at step {self.step}: {len(replay)} points', flush=True)

    def training_batch(self):
        if self.mode == 'replay' and self.step > self.last_refresh and self.step % self.sampling.refresh_every == 0:
            self.refresh()
        # Both arms draw the same full uniform batch. Replay substitutes the last
        # quarter, preserving common random numbers for the uniform component.
        uniform = super().training_batch()
        if self.mode == 'uniform':
            return uniform
        count = self.sampling.target_count
        indices = torch.randint(len(self.replay), (count,), generator=self.replay_rng)
        target = self.replay[indices].clone()
        # Half exact witnesses, half nearby points. Scale is a tangent displacement
        # in the fixed product metric; retraction normalizes all three factors.
        jitter_count = count // 2
        if jitter_count:
            x = target[:jitter_count]
            noise = torch.randn(x.shape, dtype=x.dtype, generator=self.replay_rng)
            noise = tangential(self.geo, x, noise)
            noise = noise / noise.norm(dim=1, keepdim=True).clamp_min(1e-12)
            target[:jitter_count] = retract(self.geo, x + self.sampling.jitter * noise)
        return torch.cat((uniform[:-count], target))

    def evaluate(self, points):
        m = summarize(point_metrics(self.geo, self.model, self.field, points, self.cfg['eval_chunk'], lambda: None))
        return {**m, 'nij_sup': m['nij_max_found'], 'acs_floorfree_sup': m['acs_max'],
                'acs_sup': m['acs_max'], 'conjugation_condition_max': m['G_condition_max']}

    def observe(self):
        if self.challenge is None:
            self.prepare_challenge()
        metrics = self.evaluate(self.validation)
        challenge = self.evaluate(self.challenge)
        metrics['challenge'] = challenge
        if self.initial is None:
            self.initial = metrics
        rms_ratio = metrics['nij_rms'] / max(self.initial['nij_rms'], 1e-14)
        tail_ratio = challenge['nij_sup'] / max(self.initial['challenge']['nij_sup'], 1e-14)
        metrics['selection_score'] = max(rms_ratio, tail_ratio)
        eligible = all(m['acs_floorfree_sup'] <= ACS_LIMIT and m['tangency_max'] <= 1e-10
                       and m['conjugation_condition_max'] <= CONDITION_LIMIT for m in (metrics, challenge))
        eligible = eligible and rms_ratio <= 1.05 and tail_ratio <= 1.05
        if eligible and (self.best is None or metrics['selection_score'] < self.best['metrics']['selection_score']):
            self.best = {**self.payload(), 'metrics': metrics}
        with (self.out / 'progress.jsonl').open('a') as stream:
            stream.write(json.dumps({'step': self.step, 'validation': metrics, 'eligible': eligible,
                                     'lr': self.optimizer.param_groups[0]['lr'],
                                     'best_step': self.best['step'] if self.best else None,
                                     'replay_size': len(self.replay), 'last_refresh': self.last_refresh}) + '\n')
        print(f"step {self.step} ({self.mode}): RMS={metrics['nij_rms']:.6g}, "
              f"challenge max={challenge['nij_sup']:.6g}, score={metrics['selection_score']:.6g}", flush=True)

    def check_inputs(self):
        super().check_inputs()
        if (sha256(self.source['path']) != self.source['sha256']
                or sha256(self.witness_input['report']) != self.witness_input['report_sha256']
                or sha256(self.witness_input['points']) != self.witness_input['points_sha256']):
            raise RuntimeError('CS source or witness inputs changed during training')

    def run(self, max_updates=None):
        try:
            report = super().run(max_updates)
        except Stopped:
            # Only initial challenge creation raises through the parent. It has not
            # changed model/optimizer/RNG; a resume rebuilds the same challenge.
            self.save()
            report = {'case': 's4xs1xs1', 'task': 'cs', 'status': 'incomplete',
                      'steps_done': self.step, 'steps_planned': self.cfg['steps'],
                      'source': self.source, 'numerically_validated_cs': False}
        report.update(mode=self.mode, witness_input=self.witness_input,
                      selection='minimax_relative_uniform_rms_and_fixed_challenge_max',
                      config=self.cfg)
        atomic_json(report, self.out / 'result.json')
        return report


def comparison(root, experiment='uniform-replay'):
    """Read both arms, including failed/missing runs; never promote to CS success."""
    root = Path(root)
    if experiment not in ('uniform-replay', 'replay-mixtures'):
        raise ValueError('Unknown comparison experiment')
    modes = ('uniform', 'replay') if experiment == 'uniform-replay' else ('replay25', 'replay12p5')
    report = {'experiment': experiment, 'numerically_validated_cs': False, 'arms': {}, 'comparable': False}
    lines = [f'Matched CS comparison: {experiment} (not an integrability certificate)', '']
    for mode in modes:
        result_path, diagnostic_path = root / mode / 'result.json', root / mode / 'diagnostic' / 'report.json'
        result = json.loads(result_path.read_text()) if result_path.exists() else {}
        diagnostic = json.loads(diagnostic_path.read_text()) if diagnostic_path.exists() else {}
        final = root / mode / 'final.pt'
        diagnostic_matches = bool(final.is_file() and diagnostic.get('sha256') == sha256(final))
        report['arms'][mode] = {'training': result, 'diagnostic': diagnostic, 'diagnostic_matches_final': diagnostic_matches}
        lines.append(f"{mode}: training={result.get('status', 'MISSING')}, diagnostic={diagnostic.get('status', 'MISSING')}")
        lines.append(f"  selected step={result.get('selected_step')}; fixed validation={result.get('best_validation')}")
        if diagnostic:
            lines.append(f"  fresh uniform={diagnostic.get('uniform')}")
            lines.append(f"  worst found={diagnostic.get('combined_completed_groups', {}).get('nij_max_found')}")
    a, b = (report['arms'][m] for m in modes)
    if all(arm['training'].get('status') == 'completed_pending_independent_audit'
           and arm['diagnostic'].get('status') == 'completed_diagnostic'
           and arm['diagnostic_matches_final'] for arm in (a, b)):
        ca, cb = (dict(arm['training']['config']) for arm in (a, b))
        if experiment == 'uniform-replay':
            variants_valid = ca['matched']['mode'] == 'uniform' and cb['matched']['mode'] == 'replay'
            ca['matched'] = {**ca['matched'], 'mode': 'comparison'}
            cb['matched'] = {**cb['matched'], 'mode': 'comparison'}
        else:
            variants_valid = (ca['matched']['mode'] == cb['matched']['mode'] == 'replay'
                              and ca['batch'] == cb['batch'] == 32
                              and ca['matched']['sampling']['target_count'] == 8
                              and cb['matched']['sampling']['target_count'] == 4)
            for cfg in (ca, cb):
                cfg['matched'] = {**cfg['matched'], 'sampling': {**cfg['matched']['sampling'], 'target_count': 'comparison'}}
        report['comparable'] = (variants_valid and ca == cb and a['training']['source'] == b['training']['source']
                                and a['training']['audited_base'] == b['training']['audited_base']
                                and a['training']['initial_validation'] == b['training']['initial_validation']
                                and a['training']['witness_input'] == b['training']['witness_input']
                                and a['diagnostic']['config'] == b['diagnostic']['config'])
    lines.append(f"Comparable completed pair: {report['comparable']}")
    lines.append('Review mean/tail tradeoffs and conditioning. No automatic winner or CS certification.')
    atomic_json(report, root / 'comparison.json')
    text = '\n'.join(lines) + '\n'
    (root / 'comparison.txt').write_text(text)
    print(text)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    summary = sub.add_parser('summary')
    summary.add_argument('--root', required=True)
    summary.add_argument('--experiment', choices=['uniform-replay', 'replay-mixtures'], default='uniform-replay')
    train = sub.add_parser('train')
    for key in ('config', 'base-from', 'base-audit', 'out', 'witness-dir'):
        train.add_argument('--' + key, required=True)
    train.add_argument('--mode', choices=['uniform', 'replay'], required=True)
    source = train.add_mutually_exclusive_group(required=True)
    source.add_argument('--init-from')
    source.add_argument('--resume', action='store_true')
    train.add_argument('--expected-source-sha256')
    args = parser.parse_args()
    if args.command == 'summary':
        comparison(args.root, args.experiment)
        return
    runner = Runner(read_config(args.config, args.mode), args.base_from, args.base_audit, args.out,
                    args.witness_dir, args.init_from, args.resume, args.expected_source_sha256)
    old = {sig: signal.signal(sig, lambda *_: setattr(runner, 'stop_requested', True))
           for sig in (signal.SIGTERM, signal.SIGINT)}
    try:
        report = runner.run()
    finally:
        for sig, handler in old.items():
            signal.signal(sig, handler)
    if report['status'] == 'incomplete':
        raise SystemExit(75)


if __name__ == '__main__':
    main()
