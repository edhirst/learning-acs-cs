"""Hash-bound, finite-budget CS precision campaign. Never certifies integrability."""
from __future__ import annotations

import argparse
from dataclasses import dataclass, asdict
import hashlib
import json
import math
from pathlib import Path
import signal
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import torch
import yaml

from extra.acs_training import Stopped, atomic_json, atomic_torch_save, sha256
from extra.campaign_unit import Geometry, load_acs_field, polar_corrected, make_field, nij_norm
from extra.cs_pilot import checked_base
from extra.cs_matched import diverse
from extra.cs_precision_inputs import verify
from extra.cs_precision_fields import (load_model, payload, tensor, Correction,
                                       crosscheck, six_torus_control)
from extra.cs_precision_optim import Parameters, weights, gn_step, backward_objective, loss_value
from extra.diagnose_cs import point_metrics, summarize, values, ascend
from extra.validate_ambient_acs import structured_points

REPO = Path(__file__).resolve().parents[1]
METHODS = ('adam', 'lbfgs', 'gn')
DTYPE = torch.float64
TARGET = .01
ALGEBRA_LIMIT = 1e-10


@dataclass
class Config:
    seed: int = 507
    rounds: int = 6
    work_minutes: float = 150.
    wall_minutes: float = 210.
    chunk: int = 8
    batch: int = 128
    replay_fraction: float = .125
    collocation: int = 2048
    hard_points: int = 256
    validation: int = 4096
    circle_points: int = 8
    search_pool: int = 2048
    search_starts: int = 64
    search_steps: int = 60
    adam_steps: int = 5000
    lbfgs_steps: int = 500
    gn_steps: int = 100
    eval_every: int = 25
    refresh_every: int = 100
    lr: float = 1e-5
    lr_min: float = 1e-7
    cg_steps: int = 8
    damping: float = .001
    minimum_gain: float = .01
    norm_growth_limit: float = 1.5
    condition_growth_limit: float = 2.
    tail_weight: float = .1
    tail_weight_max: float = 1.
    tail_fraction: float = .1
    correction_hidden: tuple = (64, 64)
    audit_uniform: int = 100000
    audit_pool: int = 16384
    audit_starts: int = 512
    audit_steps: int = 200
    audit_circle_points: int = 32

    def validate(self):
        self.correction_hidden = list(self.correction_hidden)
        for key, value in asdict(self).items():
            if key == 'correction_hidden':
                if not value or any(type(v) is not int or v < 1 for v in value):
                    raise ValueError('Invalid correction_hidden')
            elif isinstance(getattr(Config(), key), int):
                if type(value) is not int or value < (0 if key == 'seed' else 1):
                    raise ValueError(f'Invalid {key}')
            elif type(value) not in (float, int) or not math.isfinite(value) or value <= 0:
                raise ValueError(f'Invalid {key}')
        if (self.work_minutes >= self.wall_minutes or self.lr_min > self.lr
                or not 0 < self.replay_fraction < 1 or not 0 < self.tail_fraction <= 1
                or self.minimum_gain >= 1 or self.search_starts > self.search_pool
                or self.audit_starts > self.audit_pool or self.tail_weight > self.tail_weight_max
                or min(self.norm_growth_limit, self.condition_growth_limit) < 1):
            raise ValueError('Inconsistent precision configuration')
        return self


def read_config(path):
    raw = yaml.safe_load(Path(path).read_text())
    if not isinstance(raw, dict) or set(raw) != set(asdict(Config())):
        raise ValueError('Precision config must explicitly specify all fields')
    return Config(**raw).validate()


def code_hash():
    h = hashlib.sha256()
    for folder in ('extra', 'src'):
        for path in sorted((REPO/folder).rglob('*.py')):
            h.update(str(path.relative_to(REPO)).encode())
            h.update(path.read_bytes())
    return h.hexdigest()


class Budget:
    def __init__(self, cfg):
        self.started = time.monotonic()
        self.cfg, self.signalled = cfg, False

    def check(self, work=False):
        minutes = self.cfg.work_minutes if work else self.cfg.wall_minutes
        if self.signalled or time.monotonic()-self.started >= minutes*60:
            raise Stopped()


def receipt(root, name, data, artifacts=()):
    report = {**data, 'artifacts': {str(Path(p).resolve()): sha256(p) for p in artifacts}}
    atomic_json(report, Path(root)/f'{name}.json')
    return report


def required(root, name):
    path = Path(root)/f'{name}.json'
    report = json.loads(path.read_text())
    if report.get('status') not in ('complete', 'review', 'target_on_selection', 'budget_exhausted'):
        raise ValueError(f'{name}: predecessor did not finish safely')
    if not report.get('artifacts'):
        raise ValueError(f'{name}: missing artifact receipt')
    for file, digest in report['artifacts'].items():
        if sha256(file) != digest:
            raise ValueError(f'{name}: artifact changed: {file}')
    return report


def context(root, budget):
    manifest = required(root, 'prepare')
    if manifest['code_hash'] != code_hash():
        raise ValueError('Code changed since preparation; start a new campaign after review')
    for file, digest in manifest['inputs']['files'].items():
        if sha256(file) != digest:
            raise ValueError(f'Original input changed: {file}')
    cfg = Config(**manifest['config']).validate()
    if asdict(cfg) != asdict(budget.cfg):
        raise ValueError('Configuration changed since preparation')
    geo = Geometry('s4xs1xs1', 'ambient')
    inputs = manifest['inputs']
    evidence = checked_base(inputs['base_dir'], inputs['base_audit'])
    base = polar_corrected(geo, load_acs_field(geo, inputs['base_dir']))
    probes = torch.load(Path(root)/'probes.pt', weights_only=False)
    return manifest, geo, base, evidence, probes


def restore(path, geo, base, evidence):
    ck = torch.load(path, map_location='cpu', weights_only=False)
    for key in ('sha256', 'audit_sha256'):
        if ck.get('audited_base', {}).get(key) != evidence[key]:
            raise ValueError('Candidate ACS provenance mismatch')
    model, arch = load_model(ck, geo, base)
    return model, arch, make_field(geo, 'cs', model, base_fn=base)[0]


def metric(geo, model, field, points, cfg, budget):
    return summarize(point_metrics(geo, model, field, points, cfg.chunk, budget.check))


def conditioning(geo, model, field, points, cfg, budget):
    """Cheap same-point baseline comparison; no redundant Nijenhuis derivatives."""
    jmax, gmax = 0., 0.
    for x in points.split(cfg.chunk):
        budget.check()
        with torch.no_grad():
            P, raw, J = geo.manifold.projector(x), model(x), field(x)
            G = torch.eye(geo.N, dtype=DTYPE).expand(len(x), geo.N, geo.N)
            for stage in range(model.stages):
                G = torch.matrix_exp(P @ raw[:, stage] @ P) @ G
            sv = torch.linalg.svdvals(G)
            if not torch.isfinite(J).all() or not torch.isfinite(sv).all():
                raise FloatingPointError('Nonfinite baseline conditioning')
            jmax = max(jmax, J.norm(dim=(1, 2)).max().item())
            gmax = max(gmax, (sv[:, 0]/sv[:, -1]).max().item())
    return {'J_max': jmax, 'G_condition': gmax}


def measure(geo, model, field, probes, cfg, budget):
    uniform = metric(geo, model, field, probes['uniform'], cfg, budget)
    challenge = metric(geo, model, field, probes['challenge'], cfg, budget)
    return {'uniform': uniform, 'challenge': challenge,
            'rms': uniform['nij_rms'],
            'tail': max(uniform['nij_max_found'], challenge['nij_max_found']),
            'J_max': max(uniform['J_frobenius_max'], challenge['J_frobenius_max']),
            'G_condition': max(uniform['G_condition_max'], challenge['G_condition_max']),
            'valid': all(m['acs_max'] <= ALGEBRA_LIMIT and m['tangency_max'] <= ALGEBRA_LIMIT
                         and m['G_condition_max'] <= 1e6 for m in (uniform, challenge))}


def admissible(m, anchor, cfg):
    return (m['valid'] and m['J_max'] <= cfg.norm_growth_limit*anchor['J_max']
            and m['G_condition'] <= cfg.condition_growth_limit*anchor['G_condition'])


def score(m, anchor):
    return max(m['rms']/max(anchor['rms'], 1e-15), m['tail']/max(anchor['tail'], 1e-15))


def search(geo, field, rng, cfg, budget, audit=False, work=False):
    count, starts, steps = ((cfg.audit_pool, cfg.audit_starts, cfg.audit_steps) if audit else
                            (cfg.search_pool, cfg.search_starts, cfg.search_steps))
    check = lambda: budget.check(work)
    pool = geo.manifold.sample(count, generator=rng)
    fn = lambda x: nij_norm(geo, field, x)
    norms = values(fn, pool, cfg.chunk, check)
    initial = diverse(pool, norms, starts, .05)
    endpoints, _ = ascend(geo, fn, initial, steps, cfg.chunk, check)
    # Preserve the entire pool: its maximum is also evidence, even if diversity
    # excludes a high point from the starts.
    return torch.cat((pool, endpoints))


def prepare(root, source_run, base_dir, cfg, expected_tail, expected_revision, budget):
    root = Path(root)
    initial_code_hash = code_hash()
    if (root/'prepare.json').exists():
        raise ValueError('Campaign already prepared; never overwrite provenance')
    allowed = {'config.yaml', 'login_preflight.json', 'submission.tsv', 'logs'}
    if any(p.name not in allowed for p in root.iterdir()):
        raise ValueError('Preparation requires a new campaign output directory')
    inputs = verify(source_run, base_dir)
    tail = Path(source_run)/'final_tail.pt'
    if inputs['files'][str(tail.resolve())] != expected_tail:
        raise ValueError('Source differs from login-node preflight')
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip()
    if revision != expected_revision:
        raise ValueError('Git revision differs from the submitted revision')
    geo = Geometry('s4xs1xs1', 'ambient')
    evidence = checked_base(base_dir, inputs['base_audit'])
    base = polar_corrected(geo, load_acs_field(geo, base_dir))
    rng = torch.Generator().manual_seed(cfg.seed + 10000)
    uniform = geo.manifold.sample(cfg.validation, generator=rng)
    challenge = structured_points(geo, rng, cfg.circle_points)
    models = {}
    old_witnesses = []
    for name in ('tail', 'average'):
        path = Path(source_run)/f'final_{name}.pt'
        models[name] = restore(path, geo, base, evidence)
        cloud = torch.load(Path(source_run)/f'diagnostic_{name}/diagnostic_points.pt', weights_only=False)
        points = cloud['points']
        if cloud.get('checkpoint_sha256') != sha256(path) or cloud.get('status') != 'completed_diagnostic':
            raise ValueError('Previous witness checkpoint/status mismatch')
        if (points.dtype != DTYPE or points.ndim != 2 or points.shape[1] != 9
                or not torch.isfinite(points).all()):
            raise ValueError('Invalid previous witness points')
        for block in points.split((5, 2, 2), dim=1):
            if not torch.allclose(block.norm(dim=1), torch.ones(len(points), dtype=DTYPE), atol=1e-10, rtol=0):
                raise ValueError('Previous witnesses are off the manifold')
        indices = cloud['metrics']['nij'].topk(min(cfg.hard_points, len(points))).indices
        old_witnesses.append(points[indices])
        # New searches are selection-only; old diagnostic witnesses may train.
        challenge = torch.cat((challenge, search(geo, models[name][2], rng, cfg, budget)))
    probes = {'uniform': uniform, 'challenge': challenge}
    observations, checks = {}, {}
    for name in ('previous', 'average', 'tail'):
        if name == 'previous':
            models[name] = restore(inputs['previous'], geo, base, evidence)
        model, _, field = models[name]
        observations[name] = measure(geo, model, field, probes, cfg, budget)
        # All fields see the same random points and historical tail witnesses.
        test_points = torch.cat((uniform[:4], old_witnesses[0][:4], old_witnesses[1][:4]))
        checks[name] = crosscheck(geo, field, test_points, budget.check)
    control_geo, control = six_torus_control()
    control_points = control_geo.manifold.sample(3, generator=rng)
    controls = crosscheck(control_geo, control, control_points, budget.check)
    controls['integrable_max'] = tensor(control_geo, control, control_points).flatten(1).norm(dim=1).max().item()
    def negative(x):
        P = control_geo.manifold.projector(x)
        A = torch.zeros((12, 12), dtype=DTYPE)
        A[0, 0] = .3
        B = P @ (x[:, 4, None, None]*A) @ P
        return torch.matrix_exp(B) @ control(x) @ torch.matrix_exp(-B)
    controls['nonintegrable'] = crosscheck(control_geo, negative, control_points, budget.check)
    controls['nonintegrable_max'] = tensor(control_geo, negative, control_points).flatten(1).norm(dim=1).max().item()
    # Exercise the actual higher-order GN path on the HPC Torch build before
    # releasing three expensive pilots. This probe never mutates the source.
    probe_model = models['tail'][0]
    parameters = Parameters(probe_model)
    probe_residual = functional_residual(probe_model, parameters, geo, base)
    _, controls['gn_kernel_probe'] = gn_step(
        probe_residual, parameters.vector(), uniform[:2],
        torch.ones(min(2, len(uniform)), dtype=DTYPE)/min(2, len(uniform)),
        damping=cfg.damping, chunk=1, cg_steps=2, check=budget.check)
    if not (all(c['passed'] for c in checks.values()) and controls['passed']
            and controls['integrable_max'] < 1e-10 and controls['nonintegrable']['passed']
            and controls['nonintegrable_max'] > 1e-3 and observations['tail']['valid']):
        atomic_json({'checks': checks, 'controls': controls, 'metrics': observations}, root/'diagnostic_failure.json')
        raise ValueError('Independent checks failed; inspect diagnostic_failure.json')
    atomic_torch_save(probes, root/'probes.pt')
    atomic_torch_save(torch.cat(old_witnesses), root/'training_witnesses.pt')
    if code_hash() != initial_code_hash:
        raise ValueError('Code changed during preparation')
    for path, fingerprint in inputs['files'].items():
        if sha256(path) != fingerprint:
            raise ValueError(f'Input changed during preparation: {path}')
    report = {'status': 'complete', 'inputs': inputs, 'config': asdict(cfg),
              'revision': revision, 'code_hash': initial_code_hash, 'metrics': observations,
              'crosschecks': checks, 'controls': controls,
              'targets': {'rms': TARGET, 'max_found': TARGET, 'acs_tangency': ALGEBRA_LIMIT},
              'integrability_certified': False}
    receipt(root, 'prepare', report, [root/'probes.pt', root/'training_witnesses.pt'])


def functional_residual(model, params, geo, base):
    def residual(vector, x):
        def field(z):
            raw = torch.func.functional_call(model, params.dictionary(vector), (z,))
            P, J = geo.manifold.projector(z), base(z)
            for stage in range(model.stages):
                B = P @ raw[:, stage] @ P
                J = torch.matrix_exp(B) @ J @ torch.matrix_exp(-B)
            return J
        return tensor(geo, field, x).flatten(1)
    return residual


def train_block(root, name, source, method, cfg, budget, *, round_index=0, correction=False):
    manifest, geo, base, evidence, probes = context(root, budget)
    directory = Path(root)/name
    if directory.exists():
        raise ValueError(f'{name} already exists; inspect failure before retrying into a new campaign')
    directory.mkdir()
    model, arch, field = restore(source, geo, base, evidence)
    torch.manual_seed(cfg.seed + 30000 + round_index)
    if correction:
        model = Correction(model, base, geo, cfg.correction_hidden)
        arch = {'kind': 'correction', 'parent': arch, 'hidden': list(cfg.correction_hidden)}
        field = make_field(geo, 'cs', model, base_fn=base)[0]
    params = Parameters(model)
    residual = functional_residual(model, params, geo, base)
    source_hash = sha256(source)
    source_info = {'path': str(Path(source).resolve()), 'sha256': source_hash,
                   'method': method, 'optimizer_restarted': True, 'correction': correction}
    rng = torch.Generator().manual_seed(cfg.seed + 40000 + round_index)
    uniform = geo.manifold.sample(cfg.collocation, generator=rng)
    historical = torch.load(Path(root)/'training_witnesses.pt', weights_only=False)
    # Separate from the immutable selection set. Only this cloud enters training.
    discovered = search(geo, field, rng, cfg, budget, work=True)
    pool = torch.cat((historical, discovered))
    norms = values(lambda x: nij_norm(geo, field, x), pool, cfg.chunk, lambda: budget.check(True))
    hard = diverse(pool, norms, cfg.hard_points, .05)
    collocation = torch.cat((uniform, hard))
    initial = measure(geo, model, field, probes, cfg, budget)
    anchor = manifest['metrics']['tail']
    best, best_metric = payload(model, arch, evidence, source_info), initial
    atomic_torch_save(best, directory/'best.pt')
    parameters = [p for _, p in params.items]
    optimizer = (torch.optim.Adam(parameters, lr=cfg.lr) if method == 'adam' else
                 torch.optim.LBFGS(parameters, lr=1., max_iter=1, max_eval=8,
                                   history_size=20, line_search_fn='strong_wolfe') if method == 'lbfgs' else None)
    limit = getattr(cfg, method+'_steps')
    # A GN step can cost many Adam steps. Validate each one, and each small
    # L-BFGS block, so wall-limited jobs still retain useful intermediate fields.
    observe_every = (cfg.eval_every if method == 'adam' else
                     min(5, cfg.eval_every) if method == 'lbfgs' else 1)
    refresh_every = cfg.refresh_every if method == 'adam' else min(10, cfg.refresh_every)
    reweight_every = min(5, cfg.eval_every)
    damping, done, accepted, last_weights = cfg.damping, 0, 0, None
    # Pilots share a mean-square mixture; explicit tail weighting starts in rounds.
    tail_weight = 0. if round_index == 0 else min(cfg.tail_weight_max, cfg.tail_weight*2**(round_index-1))
    stop_reason = 'iteration_budget'
    for step in range(limit):
        before = params.vector()
        try:
            budget.check(True)
            if step and step % refresh_every == 0:
                discovered = search(geo, field, rng, cfg, budget, work=True)
                pool = torch.cat((hard, discovered))
                norms = values(lambda x: nij_norm(geo, field, x), pool, cfg.chunk, lambda: budget.check(True))
                hard = diverse(pool, norms, cfg.hard_points, .05)
                uniform = geo.manifold.sample(cfg.collocation, generator=rng)
                collocation = torch.cat((uniform, hard))
                last_weights = None
            if method == 'adam':
                n_hard = max(1, round(cfg.batch*cfg.replay_fraction))
                xu = geo.manifold.sample(cfg.batch-n_hard, generator=rng)
                xh = hard[torch.randint(len(hard), (n_hard,), generator=rng)]
                points = torch.cat((xu, xh))
                norms = values(lambda x: nij_norm(geo, field, x), points, cfg.chunk, lambda: budget.check(True))
                w = (torch.full_like(norms, 1/len(points)) if not tail_weight else
                     weights(norms, len(xu), tail_weight, cfg.tail_fraction))
                optimizer.param_groups[0]['lr'] = cfg.lr*(cfg.lr_min/cfg.lr)**(step/max(1, limit-1))
            else:
                points = collocation
                if last_weights is None or step % reweight_every == 0:
                    norms = values(lambda x: nij_norm(geo, field, x), points, cfg.chunk, lambda: budget.check(True))
                    if tail_weight:
                        next_weights = weights(norms, len(uniform), tail_weight, cfg.tail_fraction)
                    else:
                        next_weights = torch.cat((torch.full((len(uniform),), (1-cfg.replay_fraction)/len(uniform), dtype=DTYPE),
                                                  torch.full((len(hard),), cfg.replay_fraction/len(hard), dtype=DTYPE)))
                    if method == 'lbfgs' and (last_weights is None or not torch.equal(last_weights, next_weights)):
                        optimizer.state.clear()  # objective/weights changed
                    last_weights = next_weights
                w = last_weights
            check = lambda: budget.check(True)
            if method == 'gn':
                vector, info = gn_step(residual, before, points, w, damping=damping,
                                       chunk=cfg.chunk, cg_steps=cfg.cg_steps, check=check)
                params.assign(vector)
                damping = info['next_damping']
            else:
                def closure():
                    optimizer.zero_grad()
                    value = backward_objective(residual, params, points, w, cfg.chunk, check)
                    return torch.tensor(value, dtype=DTYPE)
                if method == 'adam':
                    value = closure()
                    gradient = torch.nn.utils.clip_grad_norm_(parameters, 1., error_if_nonfinite=True)
                    optimizer.step()
                    info = {'before': value.item(), 'gradient_norm': gradient.item(), 'accepted': True,
                            'lr': optimizer.param_groups[0]['lr']}
                else:
                    before_loss = loss_value(residual, before, points, w, cfg.chunk, check)
                    optimizer.step(closure)
                    after_loss = loss_value(residual, params.vector(), points, w, cfg.chunk, check)
                    good = after_loss < before_loss
                    if not good:
                        params.assign(before)
                        optimizer.state.clear()
                    info = {'before': before_loss, 'after': after_loss, 'accepted': good}
            if not torch.isfinite(params.vector()).all():
                raise FloatingPointError('Nonfinite parameters')
            done += 1
            accepted += int(info['accepted'])
        except Stopped:
            params.assign(before)
            if budget.signalled:
                raise  # no success receipt after external interruption
            stop_reason = 'work_time_budget'
            break
        except Exception:
            params.assign(before)
            raise
        if done % observe_every == 0 or done == limit:
            m = measure(geo, model, field, probes, cfg, budget)
            if admissible(m, anchor, cfg) and score(m, initial) < score(best_metric, initial):
                best_metric, best = m, payload(model, arch, evidence, source_info)
                atomic_torch_save(best, directory/'best.pt')
            with (directory/'progress.jsonl').open('a') as stream:
                stream.write(json.dumps({'step': done, 'metrics': m, 'optimizer': info,
                                         'best_score': score(best_metric, initial)})+'\n')
            print(f'{name} step={done} RMS={m["rms"]:.6g} tail={m["tail"]:.6g}', flush=True)
            if not admissible(m, anchor, cfg):
                stop_reason = 'numerical_or_growth_review'
                break
            if m['rms'] <= TARGET and m['tail'] <= TARGET:
                stop_reason = 'target_on_selection'
                break
    if done == 0:
        raise RuntimeError('No optimizer steps completed: reduce collocation/CG budget after review')
    m = measure(geo, model, field, probes, cfg, budget)
    if admissible(m, anchor, cfg) and score(m, initial) < score(best_metric, initial):
        best_metric, best = m, payload(model, arch, evidence, source_info)
    atomic_torch_save(best, directory/'final.pt')
    if sha256(source) != source_hash:
        raise ValueError('Source changed during optimization')
    # Recheck original files/code before issuing a successful receipt.
    context(root, budget)
    return receipt(root, name, {'status': 'complete', 'method': method, 'source': source_info,
        'checkpoint': str((directory/'final.pt').resolve()), 'steps': done,
        'accepted_steps': accepted, 'stop_reason': stop_reason, 'initial': initial,
        'metrics': best_metric, 'tail_weight': tail_weight,
        'elapsed_minutes': (time.monotonic()-budget.started)/60,
        'integrability_certified': False}, [directory/'final.pt'])


def decide(root, name, previous, candidates, cfg, budget, *, correction_used=False,
           confirmation=False):
    """Compare all fields on the SAME expanded selection set, never audit data."""
    if (Path(root)/f'{name}.json').exists():
        raise ValueError('Decision already exists; refusing to overwrite')
    manifest, geo, base, evidence, probes = context(root, budget)
    old = required(root, previous) if previous else None
    retained = (old['checkpoint'] if old else
                str(Path(manifest['inputs']['source_run'])/'final_tail.pt'))
    old_monitor = torch.load(old['monitor_file'], weights_only=False) if old else probes['challenge']
    rng = torch.Generator().manual_seed(cfg.seed + 60000 + (old['round']+1 if old else 0))
    entries = [('retained', retained, old['method'] if old else 'adam')]
    for candidate_name in candidates:
        r = required(root, candidate_name)
        entries.append((candidate_name, r['checkpoint'], r['method']))
    models, monitor = {}, old_monitor
    for label, path, method in entries:
        models[label] = restore(path, geo, base, evidence)
        monitor = torch.cat((monitor, search(geo, models[label][2], rng, cfg, budget)))
    probes = {**probes, 'challenge': monitor}
    measurements = {label: measure(geo, m, f, probes, cfg, budget)
                    for label, (m, _, f) in models.items()}
    # Re-evaluate the original source on these points; norm gates never reset
    # after rebasing or after the selection set grows.
    m, _, f = restore(Path(manifest['inputs']['source_run'])/'final_tail.pt', geo, base, evidence)
    anchor = measure(geo, m, f, probes, cfg, budget)
    reference = measurements['retained']
    if not admissible(reference, anchor, cfg):
        raise ValueError('Retained candidate fails expanded-set numerical/growth gates; review required')
    eligible = [(label, path, method) for label, path, method in entries
                if admissible(measurements[label], anchor, cfg)]
    if not eligible:
        raise ValueError('No numerically admissible retained candidate; manual review required')
    winner = min(eligible, key=lambda entry: score(measurements[entry[0]], reference))
    gain = 1-score(measurements[winner[0]], reference)
    promoted = winner[0] != 'retained' and gain >= cfg.minimum_gain
    if not promoted:
        winner = entries[0]
    label, path, method = winner
    selected = measurements[label]
    used = correction_used or bool(old and old['correction_used'])
    # An average-parent confirmation is informative, not a plateau test on the
    # retained tail lineage. Other unproductive blocks increment the counter.
    stalls = (0 if promoted else (old['stalls'] if confirmation and old else
                                  (old['stalls']+1 if old else 1)))
    status = 'target_on_selection' if selected['rms'] <= TARGET and selected['tail'] <= TARGET else 'complete'
    if stalls >= 2 and used and status == 'complete':
        status = 'review'
    pending = status == 'complete' and not used and ((old is None and not promoted) or stalls >= 2)
    monitor_file = Path(root)/f'{name}.monitor.pt'
    atomic_torch_save(monitor, monitor_file)
    artifacts = [monitor_file, Path(path), Path(root)/'prepare.json']
    artifacts += [Path(root)/f'{n}.json' for n in candidates]
    if old:
        artifacts.append(Path(root)/f'{previous}.json')
    report = {'status': status, 'checkpoint': path, 'method': method,
              'round': 0 if old is None else old['round']+1,
              'metrics': selected, 'comparisons': measurements, 'gain': gain,
              'promoted': promoted, 'stalls': stalls, 'correction_used': used,
              'pending_correction': pending, 'monitor_file': str(monitor_file.resolve()),
              'integrability_certified': False,
              'reason': ('target requires independent audit' if status == 'target_on_selection' else
                         'plateau after correction; review required' if status == 'review' else
                         'joint gain accepted' if promoted else 'retained prior checkpoint')}
    return receipt(root, name, report, artifacts)


def continuation(root, index, cfg, budget):
    if not 1 <= index <= cfg.rounds:
        raise ValueError('Invalid continuation round')
    if (Path(root)/f'round_{index}.json').exists():
        raise ValueError('Round receipt already exists')
    previous = 'select' if index == 1 else f'round_{index-1}'
    old = required(root, previous)
    context(root, budget)
    if old['status'] in ('review', 'target_on_selection', 'budget_exhausted'):
        return receipt(root, f'round_{index}', {**old, 'round': index,
                       'reason': 'safe terminal predecessor; training skipped'},
                       [Path(root)/f'{previous}.json', old['checkpoint'], old['monitor_file']])
    correction = old['pending_correction']
    source, method = old['checkpoint'], old['method']
    # The first normal round confirms the chosen optimizer on the average
    # parent. Later rounds use the retained winner, with independent training RNG.
    confirmation = index == 1 and not correction
    if confirmation:
        manifest = required(root, 'prepare')
        source = str(Path(manifest['inputs']['source_run'])/'final_average.pt')
    if correction:
        # Fit only the new head with a damped, identity-start correction solve.
        method = 'gn'
    block = f'continuation_{index}'
    train_block(root, block, source, method, cfg, budget, round_index=index, correction=correction)
    result = decide(root, f'round_{index}', previous, [block], cfg, budget,
                    correction_used=correction, confirmation=confirmation)
    if index == cfg.rounds and result['status'] == 'complete':
        result.update(status='budget_exhausted', reason='Finite automatic budget exhausted; audit retained candidate')
        atomic_json(result, Path(root)/f'round_{index}.json')
    return result


def audit(root, part, cfg, budget):
    if (Path(root)/f'audit_{part}.json').exists():
        raise ValueError('Audit receipt already exists')
    manifest, geo, base, evidence, _ = context(root, budget)
    final = required(root, f'round_{cfg.rounds}')
    path = final['checkpoint']
    fingerprint = sha256(path)
    model, _, field = restore(path, geo, base, evidence)
    model.requires_grad_(False)
    rng = torch.Generator().manual_seed(cfg.seed + 900000 + part*1000)
    if part == 0:
        uniform = geo.manifold.sample(cfg.audit_uniform, generator=rng)
        structured = structured_points(geo, rng, cfg.audit_circle_points)
        groups = [('uniform', uniform), ('structured', structured)]
    elif 1 <= part <= 3:
        groups = [('adversarial_pool_and_endpoints', search(geo, field, rng, cfg, budget, audit=True))]
    else:
        raise ValueError('Audit part must be 0 through 3')
    reports, clouds, metrics = {}, [], []
    for label, points in groups:
        m = point_metrics(geo, model, field, points, cfg.chunk, budget.check)
        reports[label] = summarize(m)
        clouds.append(points)
        metrics.append(m)
    points = torch.cat(clouds)
    combined = {key: torch.cat([m[key] for m in metrics]) for key in metrics[0]}
    worst = points[combined['nij'].topk(min(8, len(points))).indices]
    check = crosscheck(geo, field, worst, budget.check)
    original = Path(manifest['inputs']['source_run'])/'final_tail.pt'
    ref_model, _, ref_field = restore(original, geo, base, evidence)
    ref_model.requires_grad_(False)
    reference = conditioning(geo, ref_model, ref_field, points, cfg, budget)
    overall = summarize(combined)
    growth_passed = (overall['J_frobenius_max'] <= cfg.norm_growth_limit*reference['J_max']
                     and overall['G_condition_max'] <= cfg.condition_growth_limit*reference['G_condition']
                     and overall['G_condition_max'] <= 1e6)
    points_file = Path(root)/f'audit_{part}.points.pt'
    atomic_torch_save({'points': points, 'metrics': combined, 'checkpoint_sha256': fingerprint}, points_file)
    context(root, budget)
    if sha256(path) != fingerprint:
        raise ValueError('Audited checkpoint changed')
    return receipt(root, f'audit_{part}', {'status': 'complete', 'part': part,
        'checkpoint': path, 'sha256': fingerprint, 'groups': reports,
        'combined': overall, 'crosscheck': check, 'reference_conditioning': reference,
        'growth_passed': growth_passed,
        'seed': cfg.seed+900000+part*1000, 'integrability_certified': False},
        [path, points_file, Path(root)/f'round_{cfg.rounds}.json'])


def summary(root, cfg):
    root = Path(root)
    lines = ['CS precision campaign — finite numerical evidence, not an integrability certificate',
             'Targets: uniform RMS <= 0.01; maximum found <= 0.01; ACS/tangency <= 1e-10.']
    names = ['prepare'] + [f'pilot_{m}' for m in METHODS] + ['select']
    names += [f'round_{i}' for i in range(1, cfg.rounds+1)] + [f'audit_{i}' for i in range(4)]
    records, errors = {}, {}
    for name in names:
        try:
            records[name] = required(root, name)
            r = records[name]
            m = r.get('metrics', {})
            lines.append(f'{name}: {r["status"]}; {r.get("reason", r.get("stop_reason", ""))}'
                         + (f'; RMS={m["rms"]:.6g}; tail={m["tail"]:.6g}' if 'rms' in m else ''))
        except (OSError, ValueError, KeyError) as exc:
            detail = str(exc)
            failure = root/f'failure_{name}.json'
            if failure.exists():
                try:
                    failed = json.loads(failure.read_text())
                    detail = failed.get('error', failed.get('status', detail))
                except (OSError, ValueError):
                    pass
            errors[name] = detail
            lines.append(f'{name}: MISSING/FAILED/INVALID — {detail}')
    if 'prepare' in records:
        prepared = records['prepare']
        try:
            if prepared['code_hash'] != code_hash() or asdict(Config(**prepared['config']).validate()) != asdict(cfg):
                raise ValueError('Code/config changed since preparation')
            for file, fingerprint in prepared['inputs']['files'].items():
                if sha256(file) != fingerprint:
                    raise ValueError(f'Original input changed: {file}')
        except (OSError, ValueError, KeyError) as exc:
            errors['provenance'] = str(exc)
            lines.append(f'Provenance: INVALID — {exc}')
    audits = [records.get(f'audit_{i}') for i in range(4)]
    all_complete = not errors and all(audits) and len({a['sha256'] for a in audits}) == 1
    passed = False
    full_protocol = (all_complete and cfg.audit_uniform >= 100000 and cfg.audit_pool >= 16384
                     and cfg.audit_starts >= 512 and cfg.audit_steps >= 200
                     and cfg.audit_circle_points >= 32)
    if all_complete:
        rms = audits[0]['groups']['uniform']['nij_rms']
        worst = max(a['combined']['nij_max_found'] for a in audits)
        full_protocol = (full_protocol and audits[0]['groups']['uniform']['count'] >= 100000
                         and all(a['combined']['count'] >= cfg.audit_pool+cfg.audit_starts for a in audits[1:]))
        passed = (full_protocol and rms <= TARGET and worst <= TARGET
                  and all(a['crosscheck']['passed'] and a['growth_passed'] for a in audits)
                  and all(a['combined']['acs_max'] <= ALGEBRA_LIMIT
                          and a['combined']['tangency_max'] <= ALGEBRA_LIMIT for a in audits))
        lines.append(f'Independent audit: RMS={rms:.9g}; maximum found={worst:.9g}')
    selected = records.get(f'round_{cfg.rounds}')
    if selected:
        lines.append(f'Retained checkpoint: {selected["checkpoint"]}')
    lines += [f'Full audit budget used: {full_protocol}', f'Agreed sampled numerical targets passed: {passed}',
              'No rigorous global supremum bound or proof of an exact CS is claimed.']
    report = {'all_jobs_valid': all_complete, 'full_audit_protocol': full_protocol,
              'sampled_targets_passed': bool(passed),
              'integrability_certified': False, 'errors': errors,
              'selected': records.get(f'round_{cfg.rounds}'), 'audits': audits}
    atomic_json(report, root/'comparison.json')
    (root/'comparison.txt').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines), flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=['prepare', 'pilot', 'select', 'continue', 'audit', 'summary'])
    parser.add_argument('--root', required=True)
    parser.add_argument('--config', required=True)
    parser.add_argument('--method', choices=METHODS)
    parser.add_argument('--index', type=int)
    parser.add_argument('--source-run')
    parser.add_argument('--base-dir')
    parser.add_argument('--expected-tail')
    parser.add_argument('--revision')
    args = parser.parse_args()
    cfg = read_config(args.config)
    root = Path(args.root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    if args.phase == 'summary':
        summary(root, cfg)
        return
    budget = Budget(cfg)
    def stop(signum, frame):
        budget.signalled = True
    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    torch.set_num_threads(1)
    torch.manual_seed(cfg.seed)
    name = args.phase
    try:
        if args.phase == 'prepare':
            if not all((args.source_run, args.base_dir, args.expected_tail, args.revision)):
                raise ValueError('Preparation requires source, ACS base, expected tail hash and revision')
            prepare(root, args.source_run, args.base_dir, cfg, args.expected_tail, args.revision, budget)
        elif args.phase == 'pilot':
            if not args.method:
                raise ValueError('Pilot requires --method')
            name = f'pilot_{args.method}'
            manifest = required(root, 'prepare')
            source = Path(manifest['inputs']['source_run'])/'final_tail.pt'
            train_block(root, name, source, args.method, cfg, budget)
        elif args.phase == 'select':
            decide(root, 'select', None, [f'pilot_{m}' for m in METHODS], cfg, budget)
        elif args.phase == 'continue':
            name = f'round_{args.index}'
            continuation(root, args.index, cfg, budget)
        elif args.phase == 'audit':
            name = f'audit_{args.index}'
            audit(root, args.index, cfg, budget)
    except Stopped:
        atomic_json({'status': 'interrupted', 'phase': name}, root/f'failure_{name}.json')
        raise SystemExit(75)
    except Exception as exc:
        # Do not replace a valid prior receipt when a duplicate invocation fails.
        atomic_json({'status': 'failed', 'phase': name, 'error': repr(exc)}, root/f'failure_{name}.json')
        raise


if __name__ == '__main__':
    main()
