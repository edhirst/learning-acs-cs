"""Resumable, fixed-objective CS branches with independent final audits.

Training state is distinct from nomination: a temporary RMS/tail tradeoff never
silently resets Adam. Audits remain finite numerical evidence, not a theorem.
"""
from __future__ import annotations

import argparse
import copy
from dataclasses import asdict, dataclass, fields
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

from extra import cs_precision as p
from extra.cs_refine_inputs import verify
from extra.cs_precision_fields import Correction, payload, tensor, crosscheck
from extra.cs_precision_optim import weights
from extra.acs_training import atomic_json, atomic_torch_save, sha256, Stopped

BRANCHES = ('low', 'tail', 'correction')
ROLES = ('rms', 'tail', 'balanced')


@dataclass
class Config:
    seed: int = 607
    blocks: int = 4
    steps_per_block: int = 10000
    schedule_steps: int = 30000
    lr: float = 3e-6
    lr_min: float = 1e-7
    work_minutes: float = 170.
    wall_minutes: float = 210.
    batch: int = 128
    chunk: int = 16
    replay_fraction: float = .125
    low_weight: float = .05
    tail_weight: float = .1
    correction_weight: float = .05
    tail_fraction: float = .1
    eval_every: int = 250
    refresh_every: int = 1000
    monitor_uniform: int = 512
    monitor_challenge: int = 128
    validation: int = 4096
    circle_points: int = 8
    search_pool: int = 1024
    search_starts: int = 32
    search_steps: int = 40
    hard_points: int = 256
    correction_hidden: tuple = (64, 64)
    norm_growth_limit: float = 1.5
    condition_growth_limit: float = 2.
    audit_uniform: int = 100000
    audit_pool: int = 16384
    audit_starts: int = 512
    audit_steps: int = 200
    audit_circle_points: int = 32

    def validate(self):
        self.correction_hidden = list(self.correction_hidden)
        defaults = Config()
        for key, value in asdict(self).items():
            if key == 'correction_hidden':
                if not value or any(type(v) is not int or v < 1 for v in value):
                    raise ValueError('Invalid correction_hidden')
            elif type(getattr(defaults, key)) is int:
                if type(value) is not int or value < (0 if key == 'seed' else 1):
                    raise ValueError(f'Invalid {key}')
            elif type(value) not in (int, float) or not math.isfinite(value) or value <= 0:
                raise ValueError(f'Invalid {key}')
        replay = round(self.batch*self.replay_fraction)
        if (not 0 < self.replay_fraction < 1 or not 0 < replay < self.batch
                or self.lr_min > self.lr or self.work_minutes >= self.wall_minutes
                or self.wall_minutes > 210 or self.monitor_uniform > self.validation
                or self.tail_fraction > 1 or self.search_starts > self.search_pool
                or self.audit_starts > self.audit_pool
                or min(self.norm_growth_limit, self.condition_growth_limit) < 1):
            raise ValueError('Inconsistent refinement configuration')
        return self

    def protocol(self):
        # Reuse the already-tested geometry/provenance/audit implementation.
        return p.Config(seed=self.seed, rounds=self.blocks, work_minutes=self.work_minutes,
                        wall_minutes=self.wall_minutes, chunk=self.chunk, batch=self.batch,
                        validation=self.validation, circle_points=self.circle_points,
                        search_pool=self.search_pool, search_starts=self.search_starts,
                        search_steps=self.search_steps, hard_points=self.hard_points,
                        norm_growth_limit=self.norm_growth_limit,
                        condition_growth_limit=self.condition_growth_limit,
                        audit_uniform=self.audit_uniform, audit_pool=self.audit_pool,
                        audit_starts=self.audit_starts, audit_steps=self.audit_steps,
                        audit_circle_points=self.audit_circle_points).validate()


def read_config(path):
    raw = yaml.safe_load(Path(path).read_text())
    if not isinstance(raw, dict) or set(raw) != {f.name for f in fields(Config)}:
        raise ValueError('Refinement config must explicitly specify all fields')
    return Config(**raw).validate()


def context(root, cfg, budget):
    result = p.context(root, budget)
    if result[0]['refinement_config'] != asdict(cfg):
        raise ValueError('Refinement configuration changed')
    return result


def prepare(root, cfg, source_root, expected_sha, revision, budget):
    root = Path(root)
    if any(x.name not in {'config.yaml', 'login_preflight.json', 'submission.tsv', 'logs'} for x in root.iterdir()):
        raise ValueError('Preparation requires a new directory')
    start_hash = p.code_hash()
    inputs = verify(source_root)
    if inputs['checkpoint_sha256'] != expected_sha:
        raise ValueError('Source checkpoint differs from preflight')
    actual_revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=p.REPO, text=True).strip()
    if actual_revision != revision:
        raise ValueError('Git revision changed')
    geo = p.Geometry('s4xs1xs1', 'ambient')
    evidence = p.checked_base(inputs['base_dir'], inputs['base_audit'])
    base = p.polar_corrected(geo, p.load_acs_field(geo, inputs['base_dir']))
    model, _, field = p.restore(inputs['checkpoint'], geo, base, evidence)
    rng = torch.Generator().manual_seed(cfg.seed+10000)
    historical = []
    for path in inputs['clouds']:
        cloud = torch.load(path, weights_only=False, map_location='cpu')
        x, scores = cloud['points'], cloud['metrics']['nij']
        if (cloud['checkpoint_sha256'] != expected_sha or x.dtype != p.DTYPE
                or x.ndim != 2 or x.shape[1] != 9 or scores.shape != (len(x),)
                or not torch.isfinite(x).all() or not torch.isfinite(scores).all()):
            raise ValueError('Invalid historical audit cloud')
        for block in x.split((5, 2, 2), dim=1):
            if not torch.allclose(block.norm(dim=1), torch.ones(len(x), dtype=p.DTYPE), atol=1e-10, rtol=0):
                raise ValueError('Historical cloud is off manifold')
        historical.append(p.diverse(x, scores, min(cfg.hard_points, len(x)), .05))
    witnesses = torch.cat(historical)
    full = {'uniform': geo.manifold.sample(cfg.validation, generator=rng),
            'challenge': torch.cat((p.structured_points(geo, rng, cfg.circle_points),
                                     p.search(geo, field, rng, cfg.protocol(), budget)))}
    # Construct monitoring separately from the historical training cloud.
    # Deterministic structured grids can overlap earlier grids; random samples
    # and fresh-start searches supply the independent components.
    scores = p.values(lambda x: p.nij_norm(geo, field, x), full['challenge'], cfg.chunk, budget.check)
    small = {'uniform': full['uniform'][:cfg.monitor_uniform],
             'challenge': p.diverse(full['challenge'], scores, cfg.monitor_challenge, .05)}
    source_metrics = p.measure(geo, model, field, full, cfg.protocol(), budget)
    source_small = p.measure(geo, model, field, small, cfg.protocol(), budget)
    original = Path(inputs['source_run'])/'final_tail.pt'
    ref_model, _, ref_field = p.restore(original, geo, base, evidence)
    anchor = p.measure(geo, ref_model, ref_field, full, cfg.protocol(), budget)
    anchor_small = p.measure(geo, ref_model, ref_field, small, cfg.protocol(), budget)
    checks = crosscheck(geo, field, torch.cat((small['uniform'][:4], witnesses[:4])), budget.check)
    if not checks['passed'] or not p.admissible(source_metrics, anchor, cfg.protocol()):
        raise ValueError('Source fails independent checks or growth gates')
    for name, obj in [('probes.pt', full), ('monitor.pt', small), ('training_witnesses.pt', witnesses)]:
        atomic_torch_save(obj, root/name)
    if p.code_hash() != start_hash:
        raise ValueError('Code changed during preparation')
    for path, fingerprint in inputs['files'].items():
        if sha256(path) != fingerprint:
            raise ValueError(f'Input changed: {path}')
    return p.receipt(root, 'prepare', {
        'status': 'complete', 'inputs': inputs, 'config': asdict(cfg.protocol()),
        'refinement_config': asdict(cfg), 'revision': revision, 'code_hash': start_hash,
        'runtime': {'torch': str(torch.__version__), 'python': sys.version,
                    'device': 'cpu', 'torch_threads': torch.get_num_threads()},
        'metrics': {'tail': anchor}, 'anchor_small': anchor_small,
        'source_metrics': source_metrics, 'source_small': source_small, 'crosschecks': checks,
        'integrability_certified': False},
        [root/'probes.pt', root/'monitor.pt', root/'training_witnesses.pt', root/'config.yaml'])


def learning_rate(cfg, completed):
    fraction = min(completed/max(1, cfg.schedule_steps-1), 1.)
    return cfg.lr*(cfg.lr_min/cfg.lr)**fraction


def nomination_score(role, metrics, anchor):
    if role == 'rms':
        return metrics['rms']
    if role == 'tail':
        return metrics['tail']
    return p.score(metrics, anchor)


def snapshot(model, arch, evidence, source, optimizer, completed, rng, hard):
    result = payload(model, arch, evidence, source)
    result['training'] = {'optimizer': copy.deepcopy(optimizer.state_dict()),
                          'completed': completed, 'rng': rng.get_state().clone(),
                          'torch_rng': torch.get_rng_state().clone(), 'hard': hard.clone()}
    return result


def restore_training(ck, optimizer, rng):
    state = ck['training']
    optimizer.load_state_dict(state['optimizer'])
    rng.set_state(state['rng'])
    torch.set_rng_state(state['torch_rng'])
    return state['completed'], state['hard'].clone()


def train(root, cfg, branch, index, budget):
    if branch not in BRANCHES or not 1 <= index <= cfg.blocks:
        raise ValueError('Invalid branch/block')
    root = Path(root)
    name = f'{branch}_{index}'
    folder = root/name
    if folder.exists() or (root/f'{name}.json').exists():
        raise ValueError('Block already exists; never overwrite a training allocation')
    manifest, geo, base, evidence, full = context(root, cfg, budget)
    prior = p.required(root, f'{branch}_{index-1}') if index > 1 else None
    if prior and prior['status'] == 'review':
        return p.receipt(root, name, {**prior, 'block': index, 'steps': 0, 'elapsed_minutes': 0.,
                         'reason': 'Numerical review predecessor; skipped'},
                         [root/f'{branch}_{index-1}.json', prior['resume']])
    folder.mkdir()
    source_path = prior['resume'] if prior else manifest['inputs']['checkpoint']
    source_sha = sha256(source_path)
    model, arch, field = p.restore(source_path, geo, base, evidence)
    torch.manual_seed(cfg.seed+20000)
    if not prior and branch == 'correction':
        model = Correction(model, base, geo, cfg.correction_hidden)
        arch = {'kind': 'correction', 'parent': arch, 'hidden': cfg.correction_hidden}
        field = p.make_field(geo, 'cs', model, base_fn=base)[0]
    source = {'path': source_path, 'sha256': source_sha, 'branch': branch,
              'correction': branch == 'correction', 'optimizer_restarted': prior is None}
    parameters = [v for v in model.parameters() if v.requires_grad]
    optimizer = torch.optim.Adam(parameters, lr=cfg.lr)
    rng = torch.Generator().manual_seed(cfg.seed+30000)
    if prior:
        ck = torch.load(source_path, weights_only=False, map_location='cpu')
        completed, hard = restore_training(ck, optimizer, rng)
    else:
        completed = 0
        hard = torch.load(root/'training_witnesses.pt', weights_only=False)
    initial_step = completed
    monitor = torch.load(root/'monitor.pt', weights_only=False)
    nominations = copy.deepcopy(prior['nominations']) if prior else {}
    initial = p.measure(geo, model, field, monitor, cfg.protocol(), budget)
    if not p.admissible(initial, manifest['anchor_small'], cfg.protocol()):
        raise ValueError('Initial field fails numerical gates')
    timings = {key: 0. for key in ('optimization', 'monitoring', 'search', 'saving')}
    initial_state = snapshot(model, arch, evidence, source, optimizer, completed, rng, hard)
    safe = initial_state
    atomic_torch_save(safe, folder/'last.pt')
    if not nominations:
        for role in ROLES:
            nominations[role] = {'checkpoint': str((folder/'last.pt').resolve()),
                                 'step': completed, 'metrics': initial}
        # Initial nominees must be immutable, separate from mutable last.pt.
        atomic_torch_save(initial_state, folder/'initial.pt')
        for item in nominations.values():
            item['checkpoint'] = str((folder/'initial.pt').resolve())
    full_nominations = copy.deepcopy(prior['full_nominations']) if prior else {
        role: {'checkpoint': str((folder/'initial.pt').resolve()), 'step': completed,
               'metrics': manifest['source_metrics']} for role in ROLES}
    weight = getattr(cfg, f'{branch}_weight')
    status, reason = 'complete', 'iteration_budget'

    def observe():
        nonlocal safe, status, reason
        start = time.monotonic()
        metric = p.measure(geo, model, field, monitor, cfg.protocol(), budget)
        timings['monitoring'] += time.monotonic()-start
        if not p.admissible(metric, manifest['anchor_small'], cfg.protocol()):
            status, reason = 'review', 'Numerical/growth gate; rolled back to validated training state'
            return metric, False
        safe = snapshot(model, arch, evidence, source, optimizer, completed, rng, hard)
        start = time.monotonic()
        for role in ROLES:
            if nomination_score(role, metric, manifest['source_small']) < nomination_score(
                    role, nominations[role]['metrics'], manifest['source_small']):
                path = folder/f'best_{role}.pt'
                atomic_torch_save(safe, path)
                nominations[role] = {'checkpoint': str(path.resolve()), 'step': completed, 'metrics': metric}
        atomic_torch_save(safe, folder/'last.pt')
        timings['saving'] += time.monotonic()-start
        return metric, True

    last_observed = completed
    for _ in range(cfg.steps_per_block):
        iteration_rng = rng.get_state().clone()
        iteration_torch_rng = torch.get_rng_state().clone()
        iteration_hard = hard
        try:
            budget.check(True)
            if completed and completed % cfg.refresh_every == 0:
                start = time.monotonic()
                discovered = p.search(geo, field, rng, cfg.protocol(), budget, work=True)
                pool = torch.cat((hard, discovered))
                norms = p.values(lambda x: p.nij_norm(geo, field, x), pool, cfg.chunk, lambda: budget.check(True))
                hard = p.diverse(pool, norms, cfg.hard_points, .05)
                timings['search'] += time.monotonic()-start
            start = time.monotonic()
            n_hard = round(cfg.batch*cfg.replay_fraction)
            xu = geo.manifold.sample(cfg.batch-n_hard, generator=rng)
            xh = hard[torch.randint(len(hard), (n_hard,), generator=rng)]
            points = torch.cat((xu, xh))
            norms = p.values(lambda x: p.nij_norm(geo, field, x), points, cfg.chunk, lambda: budget.check(True))
            w = weights(norms, len(xu), weight, cfg.tail_fraction)
            optimizer.zero_grad(set_to_none=True)
            objective = 0.
            # Direct parameter backward avoids flattening/rebuilding the entire
            # network for every chunk. Frozen parents KEEP spatial derivatives.
            for x, part_weight in zip(points.split(cfg.chunk), w.split(cfg.chunk)):
                budget.check(True)
                residual = tensor(geo, field, x).flatten(1)
                loss = (residual.square().sum(1)*part_weight).sum()
                if not torch.isfinite(loss):
                    raise FloatingPointError('Nonfinite training objective')
                loss.backward()
                objective += loss.detach().item()
            gradient = torch.nn.utils.clip_grad_norm_(parameters, 1., error_if_nonfinite=True)
            lr = learning_rate(cfg, completed)
            optimizer.param_groups[0]['lr'] = lr
            # No budget checks inside/after the atomic optimizer update: a
            # partial backward can stop, but cannot change Adam state/weights.
            budget.check(True)
            optimizer.step()
            completed += 1
            timings['optimization'] += time.monotonic()-start
        except Stopped:
            # Retry an interrupted batch/search from the same RNG/replay state
            # next allocation, even if completed updates await observation.
            rng.set_state(iteration_rng)
            torch.set_rng_state(iteration_torch_rng)
            hard = iteration_hard
            optimizer.zero_grad(set_to_none=True)
            if budget.signalled:
                raise  # prior last.pt remains recoverable; no success receipt
            reason = 'work_time_budget'
            break
        if completed % cfg.eval_every == 0:
            metric, good = observe()
            last_observed = completed
            with (folder/'progress.jsonl').open('a') as stream:
                stream.write(json.dumps({'step': completed, 'metrics': metric, 'loss': objective,
                    'gradient_norm': gradient.item(), 'lr': lr, 'timings_seconds': timings})+'\n')
            print(f'{name} step={completed} RMS={metric["rms"]:.6g} tail={metric["tail"]:.6g}', flush=True)
            if not good:
                break
    if status != 'review' and last_observed != completed:
        observe()
    if completed == initial_step:
        raise RuntimeError('No updates completed; inspect timing/budget')
    # Always resume the last VALIDATED iterate, not the nominated winner.
    atomic_torch_save(safe, folder/'last.pt')
    finalist_metrics = {}
    # Validate all saved extrema on the larger common set at block boundaries.
    for role, nominee in nominations.items():
        m, _, f = p.restore(nominee['checkpoint'], geo, base, evidence)
        finalist_metrics[role] = p.measure(geo, m, f, full, cfg.protocol(), budget)
    m, _, f = p.restore(folder/'last.pt', geo, base, evidence)
    latest_full = p.measure(geo, m, f, full, cfg.protocol(), budget)
    if not p.admissible(latest_full, manifest['metrics']['tail'], cfg.protocol()):
        # Small monitors cannot replace boundary gates. Resume only from the
        # previously full-validated block input after a boundary violation.
        status, reason = 'review', 'Full-validation numerical/growth gate; rolled back to block input'
        safe = initial_state
        atomic_torch_save(safe, folder/'last.pt')
    # Keep boundary-validated extrema across ALL allocations, not just extrema
    # of the cheap monitor in the last allocation. These common probes never
    # change, so prior full metrics remain comparable and hash-bound.
    boundary_candidates = [(nominations[role], metric) for role, metric in finalist_metrics.items()]
    if status != 'review':
        boundary_candidates.append(({'checkpoint': str((folder/'last.pt').resolve()),
                                     'step': safe['training']['completed']}, latest_full))
    for nominee, metric in boundary_candidates:
        if p.admissible(metric, manifest['metrics']['tail'], cfg.protocol()):
            for role in ROLES:
                if nomination_score(role, metric, manifest['source_metrics']) < nomination_score(
                        role, full_nominations[role]['metrics'], manifest['source_metrics']):
                    full_nominations[role] = {**nominee, 'metrics': metric}
    if sha256(source_path) != source_sha:
        raise ValueError('Parent changed during training')
    context(root, cfg, budget)
    paths = [folder/'last.pt', root/'prepare.json']+[Path(n['checkpoint']) for n in nominations.values()]
    paths += [Path(n['checkpoint']) for n in full_nominations.values()]
    if prior:
        paths.append(root/f'{branch}_{index-1}.json')
    return p.receipt(root, name, {'status': status, 'reason': reason, 'branch': branch,
        'block': index, 'steps': completed-initial_step, 'completed': safe['training']['completed'],
        'optimizer_resumed': prior is not None, 'resume': str((folder/'last.pt').resolve()),
        'nominations': nominations, 'full_metrics': finalist_metrics, 'tail_weight': weight,
        'full_nominations': full_nominations,
        'attempted_latest_full_metrics': latest_full,
        'timings_seconds': timings, 'elapsed_minutes': (time.monotonic()-budget.started)/60,
        'integrability_certified': False}, paths)


def select(root, cfg, budget):
    root = Path(root)
    if (root/f'round_{cfg.blocks}.json').exists():
        raise ValueError('Selection already exists')
    manifest, geo, base, evidence, probes = context(root, cfg, budget)
    candidates = {'source': manifest['inputs']['checkpoint']}
    receipts = []
    for branch in BRANCHES:
        name = f'{branch}_{cfg.blocks}'
        report = p.required(root, name)
        receipts.append(root/f'{name}.json')
        candidates.update({f'{branch}_{role}': entry['checkpoint'] for role, entry in report['full_nominations'].items()})
        candidates[f'{branch}_latest'] = report['resume']
    # Selection points are new and common to ALL candidates. Final audits use
    # separate seeds and are never fed back into this campaign's training.
    rng = torch.Generator().manual_seed(cfg.seed+70000)
    monitor = probes['challenge']
    seen = set()
    for path in candidates.values():
        h = sha256(path)
        if h in seen:
            continue
        seen.add(h)
        _, _, f = p.restore(path, geo, base, evidence)
        monitor = torch.cat((monitor, p.search(geo, f, rng, cfg.protocol(), budget)))
    probes = {**probes, 'challenge': monitor}
    anchor_model, _, anchor_field = p.restore(Path(manifest['inputs']['source_run'])/'final_tail.pt', geo, base, evidence)
    anchor = p.measure(geo, anchor_model, anchor_field, probes, cfg.protocol(), budget)
    metrics = {}
    for name, path in candidates.items():
        model, _, field = p.restore(path, geo, base, evidence)
        metrics[name] = p.measure(geo, model, field, probes, cfg.protocol(), budget)
    eligible = [n for n, m in metrics.items() if p.admissible(m, anchor, cfg.protocol())]
    if not eligible:
        raise ValueError('No admissible candidate on expanded selection set')
    # Equal ABSOLUTE targets: choose closest to satisfying BOTH, while preserving
    # other nominees for review. No simultaneous-improvement veto is used.
    winner = min(eligible, key=lambda n: (max(metrics[n]['rms'], metrics[n]['tail']), metrics[n]['rms']))
    atomic_torch_save(monitor, root/'selection.monitor.pt')
    context(root, cfg, budget)
    return p.receipt(root, f'round_{cfg.blocks}', {'status': 'complete', 'winner': winner,
        'checkpoint': candidates[winner], 'metrics': metrics[winner], 'comparisons': metrics,
        'candidates': candidates, 'integrability_certified': False},
        receipts+[root/'prepare.json', root/'selection.monitor.pt']+[Path(v) for v in candidates.values()])


def summary(root, cfg):
    root = Path(root)
    errors, records = {}, {}
    names = ['prepare']+[f'{b}_{i}' for b in BRANCHES for i in range(1, cfg.blocks+1)]
    names += [f'round_{cfg.blocks}']+[f'audit_{i}' for i in range(4)]
    for name in names:
        try:
            records[name] = p.required(root, name)
        except (OSError, ValueError, KeyError) as exc:
            errors[name] = str(exc)
    try:
        context(root, cfg, p.Budget(cfg.protocol()))
    except (OSError, ValueError, KeyError) as exc:
        errors['provenance'] = str(exc)
    audits = [records.get(f'audit_{i}') for i in range(4)]
    selected = records.get(f'round_{cfg.blocks}')
    valid = (not errors and all(audits) and selected is not None
             and all(a['sha256'] == sha256(selected['checkpoint'])
                     and Path(a['checkpoint']).resolve() == Path(selected['checkpoint']).resolve()
                     for a in audits))
    full, passed = False, False
    lines = ['Resumable CS refinement — finite numerical evidence, NOT an integrability certificate']
    for name in names:
        r = records.get(name)
        if not r:
            lines.append(f'{name}: MISSING/INVALID: {errors[name]}')
        elif 'nominations' in r:
            lines.append(f'{name}: {r["status"]}; steps={r["steps"]}; total={r["completed"]}; {r["reason"]}')
    if valid:
        rms = audits[0]['groups']['uniform']['nij_rms']
        worst = max(a['combined']['nij_max_found'] for a in audits)
        full = (cfg.audit_uniform >= 100000 and cfg.audit_pool >= 16384 and cfg.audit_starts >= 512
                and cfg.audit_steps >= 200 and cfg.audit_circle_points >= 32
                and audits[0]['groups']['uniform']['count'] >= 100000
                and audits[0]['groups']['structured']['count'] >= 20480
                and all(a['combined']['count'] >= cfg.audit_pool+cfg.audit_starts for a in audits[1:]))
        passed = (full and rms <= p.TARGET and worst <= p.TARGET
                  and all(a['crosscheck']['passed'] and a['growth_passed']
                          and a['combined']['acs_max'] <= p.ALGEBRA_LIMIT
                          and a['combined']['tangency_max'] <= p.ALGEBRA_LIMIT for a in audits))
        lines.append(f'Independent audit: RMS={rms:.9g}; maximum found={worst:.9g}')
    if selected:
        lines.append(f'Selected {selected["winner"]}: {selected["checkpoint"]}')
    lines += [f'All receipts valid: {bool(valid)}', f'Full audit budget used: {bool(full)}',
              f'Agreed sampled targets passed: {bool(passed)}',
              'Targets: RMS <= 0.01 AND maximum found <= 0.01; ACS/tangency <= 1e-10.',
              'No rigorous global bound or proof of an exact CS is claimed.']
    report = {'all_jobs_valid': bool(valid), 'full_audit_protocol': bool(full),
              'sampled_targets_passed': bool(passed), 'integrability_certified': False,
              'errors': errors, 'selected': selected, 'audits': audits}
    atomic_json(report, root/'comparison.json')
    (root/'comparison.txt').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines), flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=['prepare', 'train', 'select', 'audit', 'summary'])
    parser.add_argument('--root', required=True)
    parser.add_argument('--config', required=True)
    parser.add_argument('--branch', choices=BRANCHES)
    parser.add_argument('--index', type=int)
    parser.add_argument('--source-root')
    parser.add_argument('--expected-sha')
    parser.add_argument('--revision')
    args = parser.parse_args()
    root, cfg = Path(args.root).resolve(), read_config(args.config)
    if args.phase == 'summary':
        summary(root, cfg)
        return
    budget = p.Budget(cfg.protocol())
    def stop(signum, frame):
        budget.signalled = True
    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    torch.set_num_threads(1)
    name = f'{args.branch}_{args.index}' if args.phase == 'train' else (
        f'audit_{args.index}' if args.phase == 'audit' else args.phase)
    try:
        if args.phase == 'prepare':
            if not all((args.source_root, args.expected_sha, args.revision)):
                raise ValueError('Preparation requires source, hash and revision')
            prepare(root, cfg, args.source_root, args.expected_sha, args.revision, budget)
        elif args.phase == 'train':
            train(root, cfg, args.branch, args.index, budget)
        elif args.phase == 'select':
            select(root, cfg, budget)
        elif args.phase == 'audit':
            context(root, cfg, budget)
            p.audit(root, args.index, cfg.protocol(), budget)
        else:
            raise ValueError('Unknown phase')
    except Stopped:
        atomic_json({'status': 'interrupted', 'phase': name}, root/f'failure_{name}.json')
        raise SystemExit(75)
    except Exception as exc:
        atomic_json({'status': 'failed', 'phase': name, 'error': repr(exc)}, root/f'failure_{name}.json')
        raise


if __name__ == '__main__':
    main()
