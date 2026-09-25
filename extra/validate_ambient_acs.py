"""Independent, read-only numerical audit of a learned ambient ACS checkpoint.

This is not a rigorous global certificate. Results are written to a NEW report,
never over a training result. --quick is diagnostic and cannot grant final
acceptance. Full defaults reserve independent random seeds and search budgets.
"""
from __future__ import annotations

import argparse
import itertools
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import torch
from scipy.linalg import sqrtm

from extra.acs_training import (
    DTYPE, Config, ascent, atomic_json, diagnostics, diverse_worst,
    load_checkpoint, restore_model, sha256, values,
)
from extra.campaign_unit import (
    Geometry, chart_field, make_field, polar_corrected, product_chart_overlap, tangential,
)
from src.geometry.sphere import EvenSphere


def structured_points(geo, generator, circle_points=16):
    """Poles × circle grids on the target; mixed pole slices on other products."""
    if geo.dims == (4, 1, 1):
        poles = torch.cat([torch.eye(5, dtype=DTYPE), -torch.eye(5, dtype=DTYPE)])
        angle = torch.arange(circle_points, dtype=DTYPE) * (2 * math.pi / circle_points)
        circle = torch.stack([angle.cos(), angle.sin()], dim=1)
        grid = torch.tensor(list(itertools.product(range(10), range(circle_points), range(circle_points))))
        exact = torch.cat([poles[grid[:, 0]], circle[grid[:, 1]], circle[grid[:, 2]]], dim=1)
        nearby = exact + .001 * torch.randn(exact.shape, dtype=DTYPE, generator=generator)
        blocks = nearby.split([5, 2, 2], dim=1)
        nearby = torch.cat([b / b.norm(dim=1, keepdim=True) for b in blocks], dim=1)
        return torch.cat([exact, nearby])
    points, offset = [], 0
    for width in (d + 1 for d in geo.dims):
        x = geo.manifold.sample(2 * width, generator=generator)
        x[:, offset:offset + width] = torch.cat([torch.eye(width), -torch.eye(width)]).to(DTYPE)
        points.append(x)
        offset += width
    return torch.cat(points)


def tangent_spectrum(geo, field, x):
    """Exclude the ambient normal kernel; frame is diagnostic only."""
    with torch.no_grad():
        P = geo.manifold.projector(x)
        _, basis = torch.linalg.eigh(P)
        E = basis[:, :, -geo.d:]
        A = E.transpose(1, 2) @ field(x) @ E
        singular = torch.linalg.svdvals(A)
        eigen = torch.linalg.eigvals(A)
    return {'minimum_tangent_singular_value': singular.min().item(),
            'maximum_tangent_condition': (singular[:, 0] / singular[:, -1].clamp_min(1e-300)).max().item(),
            'minimum_abs_imaginary_eigenvalue': eigen.imag.abs().min().item(),
            'tangent_dimension': geo.d, 'excluded_normal_modes': geo.N - geo.d}


def derivative_check(geo, field, x, generator):
    direction = torch.randn(x.shape, dtype=DTYPE, generator=generator)
    direction = tangential(geo, x, direction)
    direction /= direction.norm(dim=1, keepdim=True).clamp_min(1e-30)
    _, derivative = torch.func.jvp(field, (x,), (direction,))
    eps = 1e-5
    with torch.no_grad():
        center = field(x)
        plus, minus = field(x + eps * direction), field(x - eps * direction)
        finite_difference = (plus - minus) / (2 * eps)
        second_difference = (plus - 2 * center + minus) / eps**2
        error = (derivative - finite_difference).norm() / derivative.norm().clamp_min(1.)
    _, second = torch.func.jvp(
        lambda y: torch.func.jvp(field, (y,), (direction,))[1], (x,), (direction,))
    return {'first_derivative_fd_relative_error': error.item(),
            'first_derivative_finite': bool(torch.isfinite(derivative).all()),
            'second_derivative_finite': bool(torch.isfinite(second).all()),
            'second_difference_finite': bool(torch.isfinite(second_difference).all()),
            'second_derivative_fd_relative_error': (
                (second-second_difference).norm() / second.norm().clamp_min(1.)).item()}


def overlap_check(geo, field, generator, count=16):
    coords = []
    for d in geo.dims:
        p = torch.randn(count, d, dtype=DTYPE, generator=generator)
        p = .4 * p / p.norm(dim=1, keepdim=True)
        coords.append((p, EvenSphere.patch_change_coords(p)))
    charts = list(itertools.product((0, 1), repeat=len(geo.dims)))
    with torch.no_grad():
        tensors = {chart: chart_field(field, geo.dims, chart)(
            torch.cat([coords[i][chart[i]] for i in range(len(chart))], dim=1)) for chart in charts}
        return product_chart_overlap(tensors, coords, geo.dims, charts).item()


def audit(checkpoint_path, quick=False, seed=900001):
    """Full audit budgets intentionally fixed; quick never sets numerical success."""
    checkpoint_path = Path(checkpoint_path)
    if checkpoint_path.is_dir():
        checkpoint_path = checkpoint_path / 'final.pt'
    fingerprint = sha256(checkpoint_path)
    checkpoint = load_checkpoint(checkpoint_path)
    geo = Geometry(checkpoint['case'], 'ambient')
    with torch.random.fork_rng():
        model = restore_model(checkpoint, geo.case)
    model.eval()
    model.requires_grad_(False)
    field = make_field(geo, 'acs', model)[0]
    count, n_starts, n_steps, repeats = (512, 16, 20, 1) if quick else (100000, 1024, 200, 3)
    generator = torch.Generator().manual_seed(seed)
    uniform = geo.manifold.sample(count, generator=generator)
    structured = structured_points(geo, generator, 4 if quick else 16)
    uniform_sq, structured_sq = values(geo, field, uniform), values(geo, field, structured)
    adv_points, searches = [], []
    for i in range(repeats):
        search_seed = seed + 1000 * (i + 1)
        search_rng = torch.Generator().manual_seed(search_seed)
        pool = geo.manifold.sample(32*n_starts, generator=search_rng)
        starts = diverse_worst(pool, values(geo, field, pool), n_starts, .1)
        endpoints = ascent(geo, field, starts, n_steps)
        adv_points.append(endpoints)
        searches.append({'seed': search_seed, 'pool': len(pool), 'starts': len(starts),
                         'steps': n_steps, **diagnostics(values(geo, field, endpoints))})
        print(f"Audit search {i+1}/{repeats}: maximum found {searches[-1]['max_found']:.6g}", flush=True)
    points = torch.cat([uniform, structured, *adv_points])
    squared = values(geo, field, points)
    worst = points[squared.topk(min(16, len(points))).indices]
    report = {'schema_version': 2, 'checkpoint': str(checkpoint_path.resolve()),
              'sha256': fingerprint, 'case': geo.case, 'seed': seed,
              'audit_tier': 'quick' if quick else 'full', 'uniform_count': count,
              'structured_count': len(structured), 'raw_uniform': diagnostics(uniform_sq),
              'raw_structured': diagnostics(structured_sq), 'raw_combined': diagnostics(squared),
              'searches': searches, 'worst_points': worst.tolist(),
              'tangent_spectrum_at_worst': tangent_spectrum(geo, field, worst),
              'numerically_validated': False, 'certified_exact_acs': False,
              'note': 'Finite numerical audit, not an upper bound on the true supremum.'}
    if squared.max().sqrt().item() <= Config.admissible_max:
        corrected = polar_corrected(geo, field)
        corrected_sq = values(geo, corrected, points)
        tangency = 0.0
        with torch.no_grad():
            for block in points.split(512):
                P, J = geo.manifold.projector(block), corrected(block)
                tangency = max(tangency, (P @ J @ P - J).norm(dim=(1, 2)).max().item())
            K, P = field(worst), geo.manifold.projector(worst)
            J = corrected(worst).numpy()
        agreement, imaginary = 0.0, 0.0
        for i in range(len(worst)):
            S = (-K[i] @ K[i] + torch.eye(geo.N, dtype=DTYPE) - P[i]).numpy()
            reference = K[i].numpy() @ np.linalg.inv(sqrtm(S))
            agreement = max(agreement, float(np.linalg.norm(reference.real - J[i])))
            imaginary = max(imaginary, float(np.abs(reference.imag).max()))
        derivative = derivative_check(geo, corrected, worst[:4], generator)
        overlap = overlap_check(geo, corrected, generator)
        report['corrected'] = {**diagnostics(corrected_sq), 'tangency_max': tangency,
                               'scipy_sqrtm_max_difference': agreement, 'scipy_imaginary_max': imaginary,
                               'chart_overlap_mean': overlap, **derivative}
        passed = (uniform_sq.mean().sqrt().item() <= Config.target_rms
                  and squared.max().sqrt().item() <= Config.target_max
                  and corrected_sq.max().sqrt().item() <= 1e-10 and tangency <= 1e-12
                  and agreement <= 1e-10 and imaginary <= 1e-12 and overlap <= 1e-9
                  and derivative['first_derivative_finite'] and derivative['second_derivative_finite']
                  and derivative['second_difference_finite']
                  and derivative['first_derivative_fd_relative_error'] <= 1e-5
                  and derivative['second_derivative_fd_relative_error'] <= 1e-3)
        report['numerically_validated'] = bool(passed and not quick)
        report['screening_targets_passed'] = bool(passed)
    else:
        report['corrected'] = None
        report['screening_targets_passed'] = False
    if sha256(checkpoint_path) != fingerprint:
        raise RuntimeError('Checkpoint changed during audit; report is invalid')
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('checkpoint')
    parser.add_argument('--report', required=True)
    parser.add_argument('--quick', action='store_true')
    parser.add_argument('--seed', type=int, default=900001)
    args = parser.parse_args(argv)
    report_path = Path(args.report)
    if report_path.exists():
        parser.error('--report already exists; choose a new immutable report path')
    torch.set_num_threads(1)
    report = audit(args.checkpoint, quick=args.quick, seed=args.seed)
    atomic_json(report, report_path)
    print(f"numerically_validated={report['numerically_validated']}; report={report_path}")


if __name__ == '__main__':
    main()
