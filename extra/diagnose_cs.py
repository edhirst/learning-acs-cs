"""Read-only CS residual/conditioning diagnostic with saved adversarial witnesses."""
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

from extra.acs_training import Stopped, atomic_json, atomic_torch_save, sha256
from extra.cs_pilot import checked_base
from extra.campaign_unit import (Geometry, Net, load_acs_field, make_field,
                                 polar_corrected, nij_norm, tangential, retract)
from extra.validate_ambient_acs import structured_points

DTYPE = torch.float64


@dataclass
class Config:
    seed: int = 920107
    uniform_count: int = 16384
    circle_points: int = 16
    searches: int = 3
    pool_size: int = 4096
    starts: int = 128
    steps: int = 100
    chunk: int = 16
    wall_minutes: float = 180.0

    def validate(self):
        for key, value in asdict(self).items():
            if key != 'wall_minutes' and (type(value) is not int or value < (0 if key == 'seed' else 1)):
                raise ValueError(f'Invalid {key}')
        if self.pool_size < self.starts or not 0 < self.wall_minutes < float('inf'):
            raise ValueError('Require pool_size >= starts and positive finite wall duration')
        return self


def read_config(path):
    with open(path) as stream:
        return Config(**yaml.safe_load(stream)).validate()


def values(fn, points, chunk, check):
    output = []
    for x in points.split(chunk):
        check()
        with torch.no_grad():
            r = fn(x).detach()
        if not torch.isfinite(r).all():
            raise FloatingPointError('Nonfinite Nijenhuis diagnostic')
        output.append(r)
    return torch.cat(output)


def ascend(geo, fn, starts, steps, chunk, check=lambda: None):
    """Chunked input gradients; preserve every monotonically improving endpoint."""
    x = starts.detach().clone()
    current = values(fn, x, chunk, check)
    step = x.new_full((len(x), 1), .02)
    for _ in range(steps):
        gradients = []
        for block in x.split(chunk):
            check()
            probe = block.detach().requires_grad_(True)
            gradient, = torch.autograd.grad(fn(probe).sum(), probe)
            if not torch.isfinite(gradient).all():
                raise FloatingPointError('Nonfinite adversarial input gradient')
            gradients.append(gradient.detach())
        g = tangential(geo, x, torch.cat(gradients))
        proposal = retract(geo, x + step * g / g.norm(dim=1, keepdim=True).clamp_min(1e-12))
        score = values(fn, proposal, chunk, check)
        improved = score > current
        x = torch.where(improved[:, None], proposal, x).detach()
        current = torch.where(improved, score, current)
        step = torch.where(improved[:, None], step * 1.2, step * .5).clamp(1e-6, .3)
    return x, current


def point_metrics(geo, model, field, points, chunk, check):
    records = {key: [] for key in ('nij', 'acs', 'tangency', 'G_condition', 'J_frobenius')}
    for x in points.split(chunk):
        check()
        with torch.no_grad():
            P, J, raw = geo.manifold.projector(x), field(x), model(x)
            G = torch.eye(geo.N, dtype=DTYPE).expand(len(x), geo.N, geo.N)
            for stage in range(model.stages):
                G = torch.matrix_exp(P @ raw[:, stage] @ P) @ G
            sv = torch.linalg.svdvals(G)
            row = {'nij': nij_norm(geo, field, x), 'acs': (J @ J + P).norm(dim=(1, 2)),
                   'tangency': (P @ J @ P - J).norm(dim=(1, 2)),
                   'G_condition': sv[:, 0] / sv[:, -1], 'J_frobenius': J.norm(dim=(1, 2))}
            for key, val in row.items():
                if not torch.isfinite(val).all():
                    raise FloatingPointError(f'Nonfinite {key} diagnostic')
                records[key].append(val.detach())
    return {key: torch.cat(parts) for key, parts in records.items()}


def summarize(metrics):
    r = metrics['nij']
    return {'count': len(r), 'nij_mean': r.mean().item(), 'nij_rms': r.square().mean().sqrt().item(),
            'nij_max_found': r.max().item(), 'nij_q90': r.quantile(.9).item(),
            'nij_q99': r.quantile(.99).item(), 'fraction_nij_gt_1': (r > 1).to(DTYPE).mean().item(),
            'acs_max': metrics['acs'].max().item(), 'tangency_max': metrics['tangency'].max().item(),
            'G_condition_max': metrics['G_condition'].max().item(),
            'J_frobenius_max': metrics['J_frobenius'].max().item()}


def diagnose(cfg, checkpoint, base_dir, audit_path, out, expected_sha256=None, stop=lambda: False):
    cfg.validate()
    out, checkpoint = Path(out), Path(checkpoint)
    if out.exists() and any(out.iterdir()):
        raise ValueError('Diagnostic output must be empty; inputs/reports are never overwritten')
    evidence = checked_base(base_dir, audit_path)
    fingerprint = sha256(checkpoint)
    if expected_sha256 and fingerprint != expected_sha256:
        raise ValueError('CS checkpoint checksum differs from requested SHA256')
    ck = torch.load(checkpoint, map_location='cpu', weights_only=False)
    if any(ck.get(k) != v for k, v in {'case': 's4xs1xs1', 'task': 'cs', 'kind': 'ambient', 'model_type': 'mlp'}.items()):
        raise ValueError('Expected an ambient S4xS1xS1 CS checkpoint')
    if ck.get('audited_base', {}).get('sha256') != evidence['sha256']:
        raise ValueError('CS checkpoint uses a different or unrecorded ACS base')
    torch.set_num_threads(1)
    geo = Geometry('s4xs1xs1', 'ambient')
    with torch.random.fork_rng():
        model = Net(9, 9, ck['stages'], ck['hidden'], .5).to(DTYPE)
        model.load_state_dict(ck['model'])
        base_field = polar_corrected(geo, load_acs_field(geo, base_dir))
    model.eval().requires_grad_(False)
    field = make_field(geo, 'cs', model, base_fn=base_field)[0]
    fn = lambda x: nij_norm(geo, field, x)
    out.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    report = {'config': asdict(cfg), 'checkpoint': str(checkpoint.resolve()), 'sha256': fingerprint,
              'audited_base': evidence, 'status': 'incomplete', 'torch': torch.__version__,
              'numerically_validated_cs': False, 'searches': [],
              'note': 'Exploratory finite diagnostic, not a global bound or CS certificate.'}
    def check():
        if stop() or time.monotonic() - started >= cfg.wall_minutes * 60:
            raise Stopped()
    clouds, metrics, labels = [], [], []
    try:
        rng = torch.Generator().manual_seed(cfg.seed)
        samples = [('uniform', geo.manifold.sample(cfg.uniform_count, generator=rng)),
                   ('structured', structured_points(geo, rng, cfg.circle_points))]
        for name, points in samples:
            m = point_metrics(geo, model, field, points, cfg.chunk, check)
            clouds.append(points); metrics.append(m); labels += [name] * len(points)
            report[name] = summarize(m)
            print(f'{name}: {report[name]}', flush=True)
        for index in range(cfg.searches):
            seed = cfg.seed + 1000 * (index + 1)
            rng = torch.Generator().manual_seed(seed)
            pool = geo.manifold.sample(cfg.pool_size, generator=rng)
            scores = values(fn, pool, cfg.chunk, check)
            starts = pool[scores.topk(cfg.starts).indices]
            points, _ = ascend(geo, fn, starts, cfg.steps, cfg.chunk, check)
            m = point_metrics(geo, model, field, points, cfg.chunk, check)
            clouds.append(points); metrics.append(m); labels += [f'search_{index+1}'] * len(points)
            report['searches'].append({'seed': seed, 'pool_max_found': scores.max().item(), **summarize(m)})
            print(f'search {index+1}/{cfg.searches}: {report["searches"][-1]}', flush=True)
        report['status'] = 'completed_diagnostic'
    except Stopped:
        print('Diagnostic interrupted; reporting only completed sample groups.', flush=True)
    # Verify all inputs again even on a cooperative stop.
    if (sha256(checkpoint) != fingerprint or sha256(evidence['checkpoint']) != evidence['sha256']
            or sha256(evidence['audit']) != evidence['audit_sha256']):
        raise RuntimeError('An input changed during diagnosis; results are invalid')
    if clouds:
        points = torch.cat(clouds)
        combined = {k: torch.cat([m[k] for m in metrics]) for k in metrics[0]}
        report['combined_completed_groups'] = summarize(combined)
        worst = combined['nij'].topk(min(16, len(points))).indices.tolist()
        report['worst_points'] = [{'group': labels[i], 'point': points[i].tolist(),
                                   **{k: v[i].item() for k, v in combined.items()}} for i in worst]
        atomic_torch_save({'checkpoint_sha256': fingerprint, 'config': asdict(cfg),
                           'status': report['status'], 'points': points, 'groups': labels,
                           'metrics': combined}, out / 'diagnostic_points.pt')
        report['points_file'] = str((out / 'diagnostic_points.pt').resolve())
        report['points_sha256'] = sha256(out / 'diagnostic_points.pt')
    report['elapsed_min'] = (time.monotonic() - started) / 60
    atomic_json(report, out / 'report.json')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('config', 'checkpoint', 'base-from', 'base-audit', 'out'):
        parser.add_argument('--' + key, required=True)
    parser.add_argument('--expected-source-sha256')
    args = parser.parse_args()
    stopped = [False]
    def request_stop(*_):
        stopped[0] = True
    old = {sig: signal.signal(sig, request_stop) for sig in (signal.SIGTERM, signal.SIGINT)}
    try:
        report = diagnose(read_config(args.config), args.checkpoint, args.base_from, args.base_audit,
                          args.out, args.expected_source_sha256, lambda: stopped[0])
    finally:
        for sig, handler in old.items():
            signal.signal(sig, handler)
    print(f'{report["status"]}: {args.out}/report.json', flush=True)
    if report['status'] != 'completed_diagnostic':
        raise SystemExit(75)


if __name__ == '__main__':
    main()
