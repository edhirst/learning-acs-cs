"""Reproducible ambient ACS experiments; no reference field or CS objective.

Run through ``python extra/campaign_unit.py --config CONFIG --out NEW_DIR``.
The legacy campaign CLI is unchanged. All acceptance here is *numerical*, not
a proof of a uniform bound on the manifold. Independent audit is a later step.
"""
from __future__ import annotations

import argparse
import copy
from dataclasses import asdict, dataclass, fields
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import time

import torch
import yaml

from extra.campaign_unit import CurriculumNet, Geometry, make_field, retract, tangential

DTYPE = torch.float64
SCHEMA = 2


@dataclass
class Config:
    case: str = 's4xs1xs1'
    task: str = 'acs'
    formulation: str = 'ambient'
    seed: int = 7
    feat: str = 'all2+trilin'
    hidden: tuple = (128, 128)
    poly_steps: int = 5000
    refine_steps: int = 5000
    batch: int = 192
    sampling: str = 'uniform'
    pool_size: int = 768
    replay_size: int = 256
    replay_refresh: int = 250
    replay_starts: int = 48
    replay_ascent_steps: int = 30
    replay_min_distance: float = 0.15
    replay_jitter: float = 0.03
    tail_weight: float = 0.0
    tail_fraction: float = 0.1
    tail_ramp: float = 0.2
    poly_lr: float = 0.005
    refine_lr: float = 0.001
    lr_min: float = 0.00001
    grad_clip: float = 1.0
    eval_every: int = 500
    validation_size: int = 4096
    audit_every: int = 2000
    audit_starts: int = 128
    audit_steps: int = 50
    boundary_starts: int = 512
    boundary_steps: int = 150
    lbfgs_steps: int = 200
    lbfgs_chunk: int = 50
    lbfgs_batch: int = 1024
    admissible_max: float = 0.5
    target_max: float = 0.01
    target_rms: float = 0.0001
    checkpoint_every: int = 250
    wall_minutes: float = 45.0
    threads: int = 1

    def validate(self):
        if self.task != 'acs' or self.formulation != 'ambient':
            raise ValueError('This runner supports task: acs and formulation: ambient only')
        geo = Geometry(self.case, 'ambient')
        if geo.d % 2:
            raise ValueError('ACS experiments require even total tangent dimension')
        self.hidden = tuple(self.hidden)
        if not self.hidden or any(type(h) is not int or h < 1 for h in self.hidden):
            raise ValueError('hidden must be a nonempty list of positive integers')
        for name in ('poly_steps', 'refine_steps', 'lbfgs_steps', 'seed'):
            if type(getattr(self, name)) is not int or getattr(self, name) < 0:
                raise ValueError(f'{name} must be a nonnegative integer')
        for name in ('batch', 'pool_size', 'replay_size', 'replay_refresh',
                     'replay_starts', 'replay_ascent_steps', 'eval_every',
                     'validation_size', 'audit_every', 'audit_starts', 'audit_steps',
                     'boundary_starts', 'boundary_steps', 'lbfgs_chunk',
                     'lbfgs_batch', 'checkpoint_every', 'threads'):
            if type(getattr(self, name)) is not int or getattr(self, name) < 1:
                raise ValueError(f'{name} must be a positive integer')
        for f in fields(self):
            if isinstance(f.default, float):
                value = getattr(self, f.name)
                if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                    raise ValueError(f'{f.name} must be a finite nonnegative number')
        if self.sampling not in ('uniform', 'mixed') or self.batch < 4:
            raise ValueError('sampling must be uniform or mixed; batch must be >= 4')
        if self.pool_size < max(self.batch, self.replay_starts):
            raise ValueError('pool_size must cover batch and replay_starts')
        if not 0 < self.tail_fraction <= 1 or not 0 < self.tail_ramp <= 1:
            raise ValueError('tail_fraction and tail_ramp must be in (0, 1]')
        if not 0 < self.target_max < self.admissible_max < 1 or self.target_rms <= 0:
            raise ValueError('require 0 < target_max < admissible_max < 1 and target_rms > 0')
        if min(self.poly_lr, self.refine_lr, self.wall_minutes, self.grad_clip) <= 0:
            raise ValueError('learning rates, wall_minutes and grad_clip must be positive')
        if self.lr_min > min(self.poly_lr, self.refine_lr):
            raise ValueError('lr_min cannot exceed an initial learning rate')
        # Validate feature availability without consuming the experiment RNG.
        with torch.random.fork_rng():
            build_model(self)
        return self


def read_config(path):
    with open(path) as stream:
        raw = yaml.safe_load(stream)
    if not isinstance(raw, dict):
        raise ValueError('Configuration must be a YAML mapping')
    unknown = set(raw) - {f.name for f in fields(Config)}
    if unknown:
        raise ValueError(f'Unknown configuration keys: {sorted(unknown)}')
    return Config(**raw).validate()


def build_model(cfg):
    geo = Geometry(cfg.case, 'ambient')
    return CurriculumNet(geo.N, geo.N, cfg.hidden, feat=cfg.feat,
                         blocks=[d + 1 for d in geo.dims]).to(DTYPE)


def load_checkpoint(path):
    # Checkpoints are pickle-based; load only artifacts you trust.
    path = Path(path)
    return torch.load(path / 'final.pt' if path.is_dir() else path,
                      map_location='cpu', weights_only=False)


def restore_model(checkpoint, case):
    if checkpoint.get('case') != case or checkpoint.get('task') != 'acs':
        raise ValueError('Checkpoint case/task does not match this ACS experiment')
    if checkpoint.get('kind', checkpoint.get('formulation')) != 'ambient':
        raise ValueError('Continuation requires an ambient checkpoint')
    if checkpoint.get('model_type') != 'curriculum':
        raise ValueError('Continuation currently requires a curriculum checkpoint')
    cfg = Config(case=case, hidden=tuple(checkpoint['hidden']), feat=checkpoint['feat'])
    model = build_model(cfg)
    model.load_state_dict(checkpoint['model'])
    if not checkpoint.get('mlp_active', True):
        model.polynomial_phase()
    return model


def migrate_features(model, feat, blocks):
    """Extend a polynomial basis by zero coefficients, preserving the field."""
    from extra.campaign_unit import PolyNet
    old = model.poly
    new = PolyNet(old.n_in, old.n_mat, 0.0, feat=feat, blocks=blocks).to(DTYPE)
    destination = {key: i for i, key in enumerate(new.feature_keys())}
    if any(key not in destination for key in old.feature_keys()):
        raise ValueError('Feature migration would discard existing monomials')
    with torch.no_grad():
        for i, key in enumerate(old.feature_keys()):
            new.lin.weight[:, destination[key]] = old.lin.weight[:, i]
    model.poly = new


def residual_sq(geo, field, x):
    """Floor-free squared Frobenius residual, differentiable even at zero."""
    K = field(x)
    return (K @ K + geo.manifold.projector(x)).square().sum(dim=(-1, -2))


def values(geo, field, points, chunk=512):
    with torch.no_grad():
        result = torch.cat([residual_sq(geo, field, x) for x in points.split(chunk)])
    if not torch.isfinite(result).all():
        raise FloatingPointError('Nonfinite ACS residual')
    return result


def metrics(squared):
    return {'rms': squared.mean().sqrt().item(), 'mean': squared.sqrt().mean().item(),
            'max_found': squared.max().sqrt().item()}


def diagnostics(squared):
    r = squared.sqrt()
    return {**metrics(squared), 'q50': r.quantile(.5).item(),
            'q90': r.quantile(.9).item(), 'q99': r.quantile(.99).item(),
            'fraction_ge_1': (r >= 1).to(DTYPE).mean().item()}


def rank(metrics_, cfg):
    return (metrics_['max_found'] > cfg.admissible_max,
            metrics_['max_found'], metrics_['rms'])


def tail_loss(uniform_sq, hard_sq, weight, fraction):
    n = max(1, math.ceil(fraction * hard_sq.numel()))
    return uniform_sq.mean() + weight * hard_sq.topk(n).values.mean()


class Stopped(Exception):
    pass


def ascent(geo, field, starts, steps, check_stop=lambda: None, lr=0.08):
    """Detached, monotone, tangent ascent with per-point backtracking.

    Returns every endpoint, not just a single maximum. No parameter .grad is
    accumulated, and rejected/nonfinite proposals cannot replace a point.
    """
    x = starts.detach().clone()
    for _ in range(steps):
        check_stop()
        x.requires_grad_(True)
        r = residual_sq(geo, field, x)
        if not torch.isfinite(r).all():
            raise FloatingPointError('Nonfinite adversarial starting residual')
        g, = torch.autograd.grad(r.sum(), x)
        g = tangential(geo, x, g).detach()
        if not torch.isfinite(g).all():
            raise FloatingPointError('Nonfinite adversarial gradient')
        g = g / g.norm(dim=1, keepdim=True).clamp_min(1e-30)
        x, r = x.detach(), r.detach()
        step = x.new_full((len(x), 1), lr)
        with torch.no_grad():
            for _ in range(8):
                proposal = retract(geo, x + step * g)
                rp = residual_sq(geo, field, proposal)
                take = torch.isfinite(rp) & (rp > r)
                x = torch.where(take[:, None], proposal, x)
                r = torch.where(take, rp, r)
                step = torch.where(take[:, None], torch.zeros_like(step), step / 2)
    return x.detach()


def diverse_worst(points, scores, size, min_distance):
    """Greedy score-ordered chord-distance deduplication across sphere blocks."""
    chosen = []
    for index in scores.argsort(descending=True).tolist():
        point = points[index]
        if not chosen or (torch.stack(chosen) - point).norm(dim=1).min() >= min_distance:
            chosen.append(point)
            if len(chosen) == size:
                break
    return torch.stack(chosen).detach()


def atomic_torch_save(value, path):
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp')
    with open(temporary, 'wb') as stream:
        torch.save(value, stream)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def atomic_json(value, path):
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp')
    with open(temporary, 'w') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


class Runner:
    """CPU/float64 staged optimizer with explicit resumable phase state."""

    def __init__(self, cfg, out, resume=False, init_from=None, migrate=False,
                 expected_source_sha256=None):
        self.cfg = cfg.validate()
        self.out = Path(out)
        if resume and init_from:
            raise ValueError('--resume and --init-from are mutually exclusive')
        if expected_source_sha256 and not init_from:
            raise ValueError('--expected-source-sha256 requires --init-from')
        if resume and not (self.out / 'state.pt').is_file():
            raise ValueError('--resume requires state.pt in --out')
        if resume and (self.out / 'final.pt').exists():
            raise ValueError('This run is already complete; use --init-from in a new output directory')
        if not resume and self.out.exists() and any(self.out.iterdir()):
            raise ValueError('Output directory is not empty; use --resume or a new --out')
        torch.set_num_threads(cfg.threads)
        self.geo = Geometry(cfg.case, 'ambient')
        with torch.random.fork_rng():
            torch.manual_seed(cfg.seed)
            self.model = build_model(cfg)
        self.model.polynomial_phase()
        self.generators = {name: torch.Generator().manual_seed(cfg.seed + offset)
                           for name, offset in [('train', 101), ('replay', 202),
                                                ('audit', 303), ('validation', 404)]}
        self.validation = self.sample(cfg.validation_size, 'validation')
        self.replay = torch.empty(0, self.geo.N, dtype=DTYPE)
        self.challenge = torch.empty_like(self.replay)
        self.phase = 'poly'
        self.counts = {'poly': 0, 'refine': 0, 'lbfgs': 0}
        self.best = None
        self.best_poly = None
        self.lbfgs_log = []
        self.pending_observation = None
        self.source = None
        self.previous_elapsed = 0.0
        self.stop_requested = False
        self.start = time.monotonic()
        self.last_save = self.start
        if init_from:
            source_path = Path(init_from)
            source_path = source_path / 'final.pt' if source_path.is_dir() else source_path
            source_digest = sha256(source_path)
            if expected_source_sha256 is not None and source_digest != expected_source_sha256.lower():
                raise ValueError('Source checkpoint SHA256 does not match the audited candidate')
            checkpoint = load_checkpoint(source_path)
            if sha256(source_path) != source_digest:
                raise ValueError('Source checkpoint changed during loading')
            with torch.random.fork_rng():
                self.model = restore_model(checkpoint, cfg.case)
                if tuple(checkpoint['hidden']) != cfg.hidden:
                    raise ValueError('Continuation hidden widths must match the source')
                if checkpoint['feat'] != cfg.feat:
                    if not migrate:
                        raise ValueError('Continuation feat must match source (or use --migrate-features)')
                    migrate_features(self.model, cfg.feat, [d + 1 for d in self.geo.dims])
            if cfg.poly_steps:
                raise ValueError('Continuation requires poly_steps: 0; it must not rerun the ladder')
            if not self.model.mlp_active:
                self.model.refinement_phase(zero_output=True)
            self.phase = 'refine'
            self.source = {'path': str(source_path.resolve()), 'sha256': source_digest,
                           'feat': checkpoint['feat'], 'migrated': bool(migrate),
                           'counts': checkpoint.get('counts'),
                           'optimizer_restarted': True}
        self.set_optimizer()
        if resume:
            self.restore()
        self.out.mkdir(parents=True, exist_ok=True)
        if not resume:
            try:
                revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
                dirty = bool(subprocess.check_output(['git', 'status', '--porcelain'], text=True).strip())
            except (OSError, subprocess.CalledProcessError):
                revision, dirty = None, None
            atomic_json({'schema_version': SCHEMA, 'config': asdict(cfg), 'revision': revision,
                         'dirty_worktree': dirty, 'source': self.source,
                         'pbs_jobid': os.environ.get('PBS_JOBID'), 'torch': torch.__version__},
                        self.out / 'manifest.json')

    @property
    def field(self):
        return make_field(self.geo, 'acs', self.model)[0]

    @property
    def adam_steps(self):
        return self.counts['poly'] + self.counts['refine']

    def sample(self, n, stream):
        return self.geo.manifold.sample(n, generator=self.generators[stream])

    def check_stop(self):
        if self.stop_requested or (time.monotonic() - self.start) >= self.cfg.wall_minutes * 60:
            raise Stopped()

    def set_optimizer(self):
        if self.phase not in ('poly', 'refine'):
            self.optimizer = self.scheduler = None
            return
        lr = self.cfg.poly_lr if self.phase == 'poly' else self.cfg.refine_lr
        params = self.model.poly.parameters() if self.phase == 'poly' else self.model.parameters()
        self.optimizer = torch.optim.Adam(params, lr=lr)
        self.scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
            self.optimizer, T_max=max(1, getattr(self.cfg, self.phase + '_steps')), eta_min=self.cfg.lr_min)

    def payload(self):
        return {'schema_version': SCHEMA, 'case': self.cfg.case, 'task': 'acs', 'kind': 'ambient',
                'model_type': 'curriculum', 'hidden': list(self.cfg.hidden), 'stages': 1,
                'feat': self.model.poly.feat, 'mlp_active': self.model.mlp_active,
                'feature_keys': self.model.poly.feature_keys(),
                'blocks': [d + 1 for d in self.geo.dims], 'feature_scaling': 'none',
                'source': self.source,
                'model': copy.deepcopy(self.model.state_dict()), 'counts': dict(self.counts),
                'phase': self.phase}

    def save(self):
        state = {**self.payload(), 'config': asdict(self.cfg),
                 'optimizer': self.optimizer.state_dict() if self.optimizer else None,
                 'scheduler': self.scheduler.state_dict() if self.scheduler else None,
                 'rng': {k: g.get_state() for k, g in self.generators.items()},
                 'replay': self.replay, 'challenge': self.challenge,
                 'best': self.best, 'best_poly': self.best_poly, 'source': self.source,
                 'elapsed_seconds': self.previous_elapsed + time.monotonic() - self.start,
                 'lbfgs_log': self.lbfgs_log}
        state['pending_observation'] = self.pending_observation
        atomic_torch_save(state, self.out / 'state.pt')
        if self.best is not None:
            atomic_torch_save(self.best, self.out / 'best.pt')
        if self.best_poly is not None:
            atomic_torch_save(self.best_poly, self.out / 'best_poly.pt')
        self.last_save = time.monotonic()

    def restore(self):
        state = load_checkpoint(self.out / 'state.pt')
        if state.get('schema_version') != SCHEMA:
            raise ValueError('Cannot exactly resume this checkpoint schema; use --init-from')
        expected, actual = asdict(self.cfg), state['config'].copy()
        # Allocation duration may change; numerical trajectory settings may not.
        for key in ('wall_minutes',):
            expected.pop(key)
            actual.pop(key)
        actual['hidden'] = tuple(actual['hidden'])
        if expected != actual:
            raise ValueError('Resume configuration differs; use a new output with --init-from for continuation')
        self.model.load_state_dict(state['model'])
        if state['mlp_active']:
            self.model.refinement_phase(zero_output=False)
        else:
            self.model.polynomial_phase()
        self.phase, self.counts = state['phase'], state['counts']
        self.best, self.best_poly = state['best'], state['best_poly']
        self.replay, self.challenge = state['replay'], state['challenge']
        self.source, self.lbfgs_log = state['source'], state['lbfgs_log']
        self.pending_observation = state.get('pending_observation')
        self.previous_elapsed = state['elapsed_seconds']
        for key, value in state['rng'].items():
            self.generators[key].set_state(value)
        self.set_optimizer()
        if self.optimizer:
            self.optimizer.load_state_dict(state['optimizer'])
            self.scheduler.load_state_dict(state['scheduler'])

    def append_log(self, record):
        record = {'phase': self.phase, 'counts': dict(self.counts), **record}
        with open(self.out / 'progress.jsonl', 'a') as stream:
            stream.write(json.dumps(record, allow_nan=False) + '\n')
        print(json.dumps(record, allow_nan=False), flush=True)

    def refresh_replay(self):
        cfg = self.cfg
        rng_before = self.generators['replay'].get_state()
        pool = self.sample(cfg.pool_size, 'replay')
        candidates = torch.cat([pool, self.replay])
        score = values(self.geo, self.field, candidates)
        starts = diverse_worst(candidates, score, cfg.replay_starts, cfg.replay_min_distance)
        try:
            endpoints = ascent(self.geo, self.field, starts, cfg.replay_ascent_steps, self.check_stop)
        except Stopped:
            self.generators['replay'].set_state(rng_before)
            raise
        candidates = torch.cat([endpoints, self.replay, pool])
        score = values(self.geo, self.field, candidates)
        self.replay = diverse_worst(candidates, score, cfg.replay_size, cfg.replay_min_distance)

    def training_batch(self):
        cfg = self.cfg
        if cfg.sampling == 'uniform':
            x = self.sample(cfg.batch, 'train')
            return x, len(x)
        if not len(self.replay) or self.adam_steps % cfg.replay_refresh == 0:
            self.refresh_replay()
        n_uniform, n_pool = cfg.batch // 2, cfg.batch // 4
        n_replay = cfg.batch - n_uniform - n_pool
        uniform = self.sample(n_uniform, 'train')
        pool = self.sample(cfg.pool_size, 'train')
        worst = pool[values(self.geo, self.field, pool).topk(n_pool).indices]
        index = torch.randint(len(self.replay), (n_replay,), generator=self.generators['train'])
        replay = self.replay[index]
        noise = torch.randn(replay.shape, dtype=DTYPE, generator=self.generators['train'])
        # Keep half the selected hard examples exact; jitter the others on M.
        mask = (torch.arange(n_replay) % 2 == 0)[:, None]
        replay = retract(self.geo, replay + mask * cfg.replay_jitter * tangential(self.geo, replay, noise))
        return torch.cat([uniform, worst, replay]).detach(), n_uniform

    def objective(self, x, n_uniform):
        squared = residual_sq(self.geo, self.field, x)
        if self.cfg.tail_weight == 0:
            return squared.mean()
        budget = max(1, getattr(self.cfg, self.phase + '_steps', 1))
        ramp = min(1.0, (self.counts.get(self.phase, 0) + 1) / (self.cfg.tail_ramp * budget))
        hard = squared[n_uniform:] if n_uniform < len(x) else squared
        return tail_loss(squared[:n_uniform], hard, self.cfg.tail_weight * ramp, self.cfg.tail_fraction)

    def observe(self, audit=False, boundary=False):
        cfg = self.cfg
        field = self.field
        if audit:
            self.check_stop()
            rng_before = self.generators['audit'].get_state()
            count = cfg.boundary_starts if boundary else cfg.audit_starts
            steps = cfg.boundary_steps if boundary else cfg.audit_steps
            pool = self.sample(max(cfg.pool_size, 8 * count), 'audit')
            candidates = torch.cat([pool, self.challenge, self.replay])
            starts = diverse_worst(candidates, values(self.geo, field, candidates), count,
                                   cfg.replay_min_distance)
            try:
                endpoints = ascent(self.geo, field, starts, steps, self.check_stop)
            except Stopped:
                self.generators['audit'].set_state(rng_before)
                raise
            # Historical maxima remain challengers for future candidates.
            combined = torch.cat([self.challenge, endpoints])
            self.challenge = diverse_worst(combined, values(self.geo, field, combined),
                                           max(cfg.replay_size, cfg.boundary_starts),
                                           cfg.replay_min_distance)
        challenge = torch.cat([self.validation, self.challenge, self.replay])
        raw = diagnostics(values(self.geo, field, self.validation))
        raw['uniform_max'] = raw['max_found']
        raw['max_found'] = values(self.geo, field, challenge).max().sqrt().item()
        current = {**self.payload(), 'metrics': raw}
        # Compare both on IDENTICAL current points; never compare stale maxima.
        for attr in ('best', 'best_poly') if self.phase == 'poly' else ('best',):
            previous = getattr(self, attr)
            if previous is not None:
                with torch.random.fork_rng():
                    old_model = restore_model(previous, cfg.case)
                old_field = make_field(self.geo, 'acs', old_model)[0]
                old = diagnostics(values(self.geo, old_field, self.validation))
                old['uniform_max'] = old['max_found']
                old['max_found'] = values(self.geo, old_field, challenge).max().sqrt().item()
                previous['metrics'] = old
                if rank(old, cfg) <= rank(raw, cfg):
                    continue
            setattr(self, attr, copy.deepcopy(current))
        self.append_log({'raw': raw, 'best': self.best['metrics'], 'audit': audit,
                         'replay_points': len(self.replay), 'challenge_points': len(self.challenge),
                         'elapsed_seconds': self.previous_elapsed + time.monotonic() - self.start,
                         'parameter_norm': sum(p.detach().square().sum() for p in self.model.parameters()).sqrt().item()})
        return raw

    def transition(self):
        # Re-run a boundary audit before selecting a phase endpoint.
        self.observe(audit=True, boundary=True)
        if self.phase == 'poly':
            self.model.load_state_dict(self.best_poly['model'])
            m = self.best_poly['metrics']
            if m['max_found'] <= self.cfg.target_max and m['rms'] <= self.cfg.target_rms:
                self.phase = 'lbfgs'
                self.append_log({'skipped_refinement': 'polynomial candidate meets screening targets'})
            else:
                self.model.refinement_phase(zero_output=True)
                self.phase = 'refine'
        else:
            self.model.load_state_dict(self.best['model'])
            if self.best['mlp_active']:
                self.model.refinement_phase(zero_output=False)
            else:
                self.model.polynomial_phase()
            self.phase = 'lbfgs'
        self.set_optimizer()
        self.save()

    def adam_step(self):
        x, n_uniform = self.training_batch()
        self.optimizer.zero_grad(set_to_none=True)
        loss = self.objective(x, n_uniform)
        if not torch.isfinite(loss):
            raise FloatingPointError('Nonfinite ACS loss')
        loss.backward()
        torch.nn.utils.clip_grad_norm_(self.model.parameters(), self.cfg.grad_clip, error_if_nonfinite=True)
        self.optimizer.step()
        self.scheduler.step()
        self.counts[self.phase] += 1

    def lbfgs_step(self):
        cfg = self.cfg
        remaining = cfg.lbfgs_steps - self.counts['lbfgs']
        if remaining <= 0 or self.best['metrics']['max_found'] > cfg.admissible_max:
            self.phase = 'done'
            return
        previous = copy.deepcopy(self.model.state_dict())
        challenge = torch.cat([self.validation, self.challenge, self.replay])
        before_sq = values(self.geo, self.field, challenge)
        before = metrics(before_sq)
        rng_before = self.generators['train'].get_state()
        x = torch.cat([self.sample(cfg.lbfgs_batch, 'train'), self.replay])
        opt = torch.optim.LBFGS([p for p in self.model.parameters() if p.requires_grad],
                              max_iter=min(cfg.lbfgs_chunk, remaining), line_search_fn='strong_wolfe')

        def closure():
            self.check_stop()
            opt.zero_grad(set_to_none=True)
            loss = self.objective(x, cfg.lbfgs_batch)
            if not torch.isfinite(loss):
                raise FloatingPointError('Nonfinite L-BFGS loss')
            loss.backward()
            if any(p.grad is not None and not torch.isfinite(p.grad).all() for p in self.model.parameters()):
                raise FloatingPointError('Nonfinite L-BFGS gradient')
            return loss

        try:
            opt.step(closure)
            after = metrics(values(self.geo, self.field, challenge))
            accept = after['max_found'] <= before['max_found'] and after['rms'] <= before['rms']
        except Stopped:
            self.model.load_state_dict(previous)
            self.generators['train'].set_state(rng_before)
            raise
        except FloatingPointError:
            after, accept = None, False
        n_iter = next(iter(opt.state.values()), {}).get('n_iter', 0)
        self.counts['lbfgs'] += n_iter
        if not accept:
            self.model.load_state_dict(previous)
        self.lbfgs_log.append({'iterations': n_iter, 'accepted': accept, 'before': before, 'after': after})
        # Record the terminal decision before an interruptible audit.
        if not accept or n_iter == 0:
            self.phase = 'done'
        self.pending_observation = {'audit': True}
        self.observe(audit=True)
        self.pending_observation = None
        self.save()

    def run(self, max_updates=None):
        """max_updates is a deterministic test interruption at an Adam boundary."""
        updates = 0
        status = 'incomplete'
        try:
            self.check_stop()
            if self.best is None:
                self.observe()
                self.save()
            if self.pending_observation is not None:
                self.observe(**self.pending_observation)
                self.pending_observation = None
            while self.phase != 'done':
                self.check_stop()
                if self.pending_observation is not None:
                    self.observe(**self.pending_observation)
                    self.pending_observation = None
                if self.phase == 'lbfgs':
                    self.lbfgs_step()
                    continue
                if self.counts[self.phase] >= getattr(self.cfg, self.phase + '_steps'):
                    self.transition()
                    continue
                self.adam_step()
                updates += 1
                audit = self.adam_steps % self.cfg.audit_every == 0
                if audit or self.adam_steps % self.cfg.eval_every == 0:
                    self.pending_observation = {'audit': audit}
                    self.observe(audit=audit)
                    self.pending_observation = None
                if self.adam_steps % self.cfg.checkpoint_every == 0 or time.monotonic() - self.last_save >= 300:
                    self.save()
                if max_updates is not None and updates >= max_updates:
                    raise Stopped()
            self.observe(audit=True, boundary=True)
            self.save()
            # final.pt is the selected candidate, state.pt is the exact last state.
            atomic_torch_save(self.best, self.out / 'final.pt')
            status = 'completed_pending_independent_audit'
        except Stopped:
            self.save()
        except Exception:
            # Retain the most recent atomically saved state, not a possibly corrupt update.
            raise
        result = {'schema_version': SCHEMA, 'case': self.cfg.case, 'task': 'acs',
                  'formulation': 'ambient', 'seed': self.cfg.seed, 'feat': self.cfg.feat,
                  'status': status, 'phase': self.phase, 'counts': self.counts,
                  'source': self.source,
                  'steps_done': self.adam_steps, 'lbfgs': self.lbfgs_log,
                  'stopped_by_wall': status == 'incomplete',
                  'elapsed_min': (self.previous_elapsed + time.monotonic() - self.start) / 60,
                  'best_raw': self.best['metrics'] if self.best else None,
                  'certified_exact_acs': False,
                  'note': 'Sampled/adversarial maxima are lower bounds on the true supremum; not a certificate.'}
        if self.best:
            m = self.best['metrics']
            result.update(dense={'acs_mean': m['mean'], 'acs_sup': m['uniform_max'], 'acs_rms': m['rms']},
                          adv_acs_sup=m['max_found'])
        atomic_json(result, self.out / 'result.json')
        print(f"{status}: {self.out}", flush=True)
        return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--seed', type=int)
    parser.add_argument('--wall-minutes', type=float)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--resume', action='store_true')
    mode.add_argument('--init-from')
    parser.add_argument('--migrate-features', action='store_true')
    parser.add_argument('--expected-source-sha256', help='Refuse continuation unless this audited source checksum matches')
    args = parser.parse_args(argv)
    if args.migrate_features and not args.init_from:
        parser.error('--migrate-features requires --init-from')
    cfg = read_config(args.config)
    if args.seed is not None:
        cfg.seed = args.seed
    if args.wall_minutes is not None:
        cfg.wall_minutes = args.wall_minutes
    runner = Runner(cfg, args.out, resume=args.resume, init_from=args.init_from,
                    migrate=args.migrate_features, expected_source_sha256=args.expected_source_sha256)
    old_handlers = {}
    for sig in (signal.SIGTERM, signal.SIGINT):
        old_handlers[sig] = signal.signal(sig, lambda *_: setattr(runner, 'stop_requested', True))
    try:
        result = runner.run()
    finally:
        for sig, handler in old_handlers.items():
            signal.signal(sig, handler)
    # A wall stop is not success; afterok must not imply training completed.
    if result['status'] == 'incomplete':
        raise SystemExit(75)


if __name__ == '__main__':
    main()
