"""Representation, optimization and restart tests for the ambient ACS runner.

Analytic coefficients/oracles are confined to this test module. Production
training sees only sampled coordinates and the algebraic residual K² + P.
"""
import copy
import itertools
import json

import pytest
import torch

from extra.acs_training import (
    Config, Runner, ascent, build_model, diverse_worst, load_checkpoint,
    metrics, migrate_features, rank, read_config, residual_sq, restore_model,
    tail_loss, values,
)
from extra.campaign_unit import (
    Geometry, PolyNet, canonical_J0, load_acs_field, make_field, polar_corrected,
)
from extra.test_ambient_products import _exact_s4xt2_acs

DTYPE = torch.float64


@pytest.fixture(autouse=True)
def one_thread():
    previous = torch.get_num_threads()
    torch.set_num_threads(1)
    yield
    torch.set_num_threads(previous)


def tiny(**overrides):
    args = dict(hidden=(8,), poly_steps=4, refine_steps=4, batch=12, pool_size=24,
                replay_size=8, replay_refresh=2, replay_starts=4, replay_ascent_steps=2,
                validation_size=20, eval_every=2, audit_every=4, audit_starts=4,
                audit_steps=2, boundary_starts=4, boundary_steps=2, lbfgs_steps=0,
                checkpoint_every=2, wall_minutes=5)
    args.update(overrides)
    return Config(**args).validate()


def oracle_coefficients(model):
    """Explicit ambient polynomial M with PMP equal to the test-only oracle."""
    keys = {key: i for i, key in enumerate(model.feature_keys())}
    weights = torch.zeros_like(model.lin.weight).reshape(9, 9, -1)
    Q = canonical_J0(6)[:5, :5]
    R = torch.tensor([[0., -1.], [1., 0.]], dtype=DTYPE)

    def add(i, j, key, value):
        weights[i, j, keys[tuple(sorted(key))]] += value
        weights[j, i, keys[tuple(sorted(key))]] -= value

    for i in range(5):
        for j in range(i + 1, 5):
            add(i, j, (), Q[i, j])
    for i, a, j, b in itertools.product(range(5), range(5), range(2), range(2)):
        add(i, 5 + j, (a, 5 + b), Q[i, a] * R[j, b])
    for j, b in itertools.product(range(2), repeat=2):
        add(4, 7 + j, (7 + b,), -R[j, b])
    for i, j, a, b in itertools.product(range(2), repeat=4):
        add(5 + i, 7 + j, (4, 5 + a, 7 + b), -R[i, a] * R[j, b])
    with torch.no_grad():
        model.lin.weight.copy_(weights.reshape_as(model.lin.weight))


def test_features_counts_order_and_generic_four_factors():
    p = PolyNet(9, 9, .1, feat='all2+trilin', blocks=[5, 2, 2]).to(DTYPE)
    assert p.lin.in_features == 75
    assert len(p.feature_keys()) == len(set(p.feature_keys())) == 75
    assert p.feature_keys()[-20:] == list(itertools.product(range(5), range(5, 7), range(7, 9)))
    x = torch.randn(3, 9, dtype=DTYPE)
    expected = torch.stack([x[:, key].prod(1) if key else torch.ones(3, dtype=DTYPE)
                            for key in p.feature_keys()], dim=1)
    torch.testing.assert_close(p.features(x), expected, atol=0, rtol=0)
    q = PolyNet(8, 8, .1, feat='trilin', blocks=[2] * 4)
    assert q.lin.in_features == 4 * 2**3
    with pytest.raises(ValueError, match='three'):
        PolyNet(6, 6, .1, feat='all2+trilin', blocks=[3, 3])


def test_exact_acs_representable_by_explicit_coefficients():
    geo = Geometry('s4xs1xs1')
    poly = PolyNet(9, 9, .1, feat='all2+trilin', blocks=[5, 2, 2]).to(DTYPE)
    oracle_coefficients(poly)
    field = make_field(geo, 'acs', poly)[0]
    x = geo.manifold.sample(512, generator=torch.Generator().manual_seed(712))
    torch.testing.assert_close(field(x), _exact_s4xt2_acs(x), atol=2e-15, rtol=0)
    assert values(geo, field, x).max().sqrt() < 5e-15


def test_triple_parity_is_missing_from_quadratic_space():
    geo = Geometry('s4xs1xs1')
    x = geo.manifold.sample(5, generator=torch.Generator().manual_seed(33))
    poly = PolyNet(9, 9, .2, feat='all2', blocks=[5, 2, 2]).to(DTYPE)
    field = make_field(geo, 'acs', poly)[0]
    odd_poly, odd_oracle = torch.zeros(5, 2, 2), torch.zeros(5, 2, 2)
    for a, b, c in itertools.product((-1, 1), repeat=3):
        p = x * x.new_tensor([a] * 5 + [b] * 2 + [c] * 2)
        odd_poly = odd_poly + a*b*c * field(p)[:, 5:7, 7:9]
        odd_oracle = odd_oracle + a*b*c * _exact_s4xt2_acs(p)[:, 5:7, 7:9]
    assert odd_poly.abs().max() < 1e-14
    assert odd_oracle.abs().max() > .1


def test_polynomial_phase_and_zero_output_refinement():
    model = build_model(tiny())
    x = torch.randn(4, 9, dtype=DTYPE)
    model.polynomial_phase()
    before = model(x).detach()
    torch.testing.assert_close(before, model.poly(x), atol=0, rtol=0)
    model(x).square().mean().backward()
    assert all(p.grad is None for p in model.mlp.parameters())
    model.refinement_phase()
    torch.testing.assert_close(model(x), before, atol=0, rtol=0)
    model.zero_grad(set_to_none=True)
    model(x).square().mean().backward()
    assert model.mlp.net[-1].weight.grad.abs().max() > 0
    assert model.mlp.net[0].weight.abs().max() > 0


def test_migration_preserves_old_field_and_legacy_state_dict():
    model = build_model(tiny(feat='all2'))
    state = copy.deepcopy(model.state_dict())
    assert not any('tri' in key for key in state)
    build_model(tiny(feat='all2')).load_state_dict(state, strict=True)
    x = torch.randn(6, 9, dtype=DTYPE)
    before = model(x).detach()
    migrate_features(model, 'all2+trilin', [5, 2, 2])
    torch.testing.assert_close(model(x), before, atol=2e-16, rtol=1e-15)
    assert torch.count_nonzero(model.poly.lin.weight[:, -20:]) == 0
    with pytest.raises(ValueError, match='discard'):
        migrate_features(model, 'lin', [5, 2, 2])


def test_tangent_and_loss_gradients_finite_difference():
    runner_geo = Geometry('s4xs1xs1')
    model = build_model(tiny())
    model.polynomial_phase()
    field = make_field(runner_geo, 'acs', model)[0]
    x = runner_geo.sample(3).requires_grad_(True)
    P = runner_geo.manifold.projector(x)
    K = field(x)
    torch.testing.assert_close(P @ K @ P, K, atol=1e-15, rtol=1e-14)
    loss = residual_sq(runner_geo, field, x).mean()
    dx, dw = torch.autograd.grad(loss, (x, model.poly.lin.weight))
    assert torch.isfinite(dx).all() and torch.isfinite(dw).all()
    eps = 1e-6
    direction = torch.randn_like(x)
    fd = (residual_sq(runner_geo, field, x + eps*direction).mean()
          - residual_sq(runner_geo, field, x - eps*direction).mean()) / (2*eps)
    torch.testing.assert_close(fd, (dx*direction).sum(), rtol=1e-5, atol=1e-8)


def test_tail_loss_uses_ceiling_and_uniform_anchor():
    uniform = torch.tensor([1., 3.], dtype=DTYPE, requires_grad=True)
    hard = torch.arange(1., 12., dtype=DTYPE, requires_grad=True)
    loss = tail_loss(uniform, hard, .5, .1)
    assert loss == 2 + .5 * 10.5
    loss.backward()
    assert torch.count_nonzero(hard.grad) == 2
    assert torch.all(uniform.grad > 0)


def test_adversarial_search_finds_manufactured_local_defect():
    geo = Geometry('s4xs1xs1')
    target = geo.manifold.sample(1, generator=torch.Generator().manual_seed(7))
    def field(x):
        bump = torch.exp(-8 * (x-target).square().sum(1))
        return (1 + bump[:, None, None]) * _exact_s4xt2_acs(x)
    x = geo.manifold.sample(256, generator=torch.Generator().manual_seed(19))
    v = values(geo, field, x)
    starts = x[v.topk(12).indices]
    endpoints = ascent(geo, field, starts, 80)
    end = values(geo, field, endpoints)
    assert torch.all(end >= values(geo, field, starts) - 1e-12)
    assert end.max().sqrt() > 7.0
    assert v.max().sqrt() < 1.0  # Uniform samples miss this constructed defect.
    for block in endpoints.split([5, 2, 2], dim=1):
        torch.testing.assert_close(block.norm(dim=1), torch.ones(len(block), dtype=DTYPE))


def test_diverse_replay_and_parameter_grad_isolation():
    cfg = tiny(sampling='mixed')
    geo = Geometry(cfg.case)
    model = build_model(cfg)
    field = make_field(geo, 'acs', model)[0]
    x = geo.sample(5)
    points = torch.cat([x, x])
    selected = diverse_worst(points, torch.arange(10.), 10, .01)
    assert len(selected) == 5
    ascent(geo, field, x, 2)
    assert all(p.grad is None for p in model.parameters())


@pytest.mark.parametrize('interrupt_after', [2, 4, 6])
def test_exact_restart_across_phases_with_replay(tmp_path, interrupt_after):
    cfg = tiny(sampling='mixed', tail_weight=1)
    full = Runner(cfg, tmp_path / 'full')
    full.run()
    part = Runner(cfg, tmp_path / 'part')
    result = part.run(max_updates=interrupt_after)
    assert result['status'] == 'incomplete'
    assert not (tmp_path / 'part' / 'final.pt').exists()
    resumed = Runner(cfg, tmp_path / 'part', resume=True)
    resumed.run()
    assert resumed.counts == full.counts
    for key, value in full.model.state_dict().items():
        torch.testing.assert_close(resumed.model.state_dict()[key], value, atol=0, rtol=0)
    for key in full.generators:
        assert torch.equal(resumed.generators[key].get_state(), full.generators[key].get_state())
    torch.testing.assert_close(full.replay, resumed.replay, atol=0, rtol=0)
    assert full.best['metrics'] == resumed.best['metrics']


def test_evaluation_does_not_consume_training_rng(tmp_path):
    runner = Runner(tiny(), tmp_path)
    state = runner.generators['train'].get_state().clone()
    global_state = torch.get_rng_state().clone()
    runner.observe(audit=True)
    assert torch.equal(state, runner.generators['train'].get_state())
    assert torch.equal(global_state, torch.get_rng_state())


def test_checkpoint_roundtrip_legacy_loader_and_best_poly(tmp_path):
    cfg = tiny()
    runner = Runner(cfg, tmp_path)
    runner.run()
    checkpoint = load_checkpoint(tmp_path)
    model = restore_model(checkpoint, cfg.case)
    x = runner.geo.sample(5)
    torch.testing.assert_close(load_acs_field(runner.geo, str(tmp_path))(x),
                               make_field(runner.geo, 'acs', model)[0](x), atol=0, rtol=0)
    assert load_checkpoint(tmp_path / 'best_poly.pt')['mlp_active'] is False
    assert not json.loads((tmp_path / 'result.json').read_text())['certified_exact_acs']


def test_legacy_continuation_preserves_field_and_restarts_optimizer(tmp_path):
    cfg = tiny(feat='all2')
    model = build_model(cfg)
    checkpoint = dict(model=model.state_dict(), case=cfg.case, task='acs', kind='ambient',
                      model_type='curriculum', hidden=list(cfg.hidden), stages=1, feat='all2')
    source = tmp_path / 'legacy.pt'
    torch.save(checkpoint, source)
    cfg.poly_steps = 0
    runner = Runner(cfg, tmp_path / 'continued', init_from=source)
    x = runner.geo.sample(5)
    torch.testing.assert_close(runner.model(x), model(x), atol=0, rtol=0)
    assert runner.phase == 'refine' and not runner.optimizer.state
    assert runner.scheduler.T_max == cfg.refine_steps
    assert runner.source['optimizer_restarted']
    cfg.feat = 'all2+trilin'
    with pytest.raises(ValueError, match='feat must match'):
        Runner(cfg, tmp_path / 'bad', init_from=source)
    migrated = Runner(cfg, tmp_path / 'migrated', init_from=source, migrate=True)
    torch.testing.assert_close(migrated.model(x), model(x), atol=1e-15, rtol=1e-15)


def test_wall_stop_and_configuration_guards(tmp_path):
    cfg = tiny()
    runner = Runner(cfg, tmp_path / 'stopped')
    runner.stop_requested = True
    assert runner.run()['status'] == 'incomplete'
    assert (runner.out / 'state.pt').exists()
    assert not (runner.out / 'final.pt').exists()
    cfg.wall_minutes = 10
    Runner(cfg, runner.out, resume=True)  # Allocation limit can change.
    changed = copy.deepcopy(cfg)
    changed.refine_steps += 1
    with pytest.raises(ValueError, match='configuration differs'):
        Runner(changed, runner.out, resume=True)
    with pytest.raises(ValueError, match='not empty'):
        Runner(cfg, runner.out)
    with pytest.raises(ValueError, match='only'):
        tiny(task='cs')
    with pytest.raises(ValueError, match='nonnegative integer'):
        tiny(poly_steps=True)


def test_raw_metrics_have_no_reporting_floor_and_rank_worst_first():
    assert metrics(torch.zeros(5)) == dict(rms=0., mean=0., max_found=0.)
    a = dict(rms=.01, max_found=.6)
    b = dict(rms=.02, max_found=.4)
    assert rank(b, tiny()) < rank(a, tiny())


def test_polar_root_against_scipy_on_nonnormal_field():
    import numpy as np
    from scipy.linalg import sqrtm
    geo = Geometry('s4xs1xs1')
    x = geo.manifold.sample(6, generator=torch.Generator().manual_seed(42))
    def field(p):
        P = geo.manifold.projector(p)
        scale = torch.diag(p.new_tensor([1., 1.1, 1.2, 1., .9, 1., .8, 1., 1.]))
        A = P @ scale @ P
        I = torch.eye(9, dtype=DTYPE).expand_as(P)
        A = A + I - P
        return 1.02 * A @ _exact_s4xt2_acs(p) @ torch.linalg.inv(A)
    K = field(x)
    assert (K @ K.transpose(1, 2) - K.transpose(1, 2) @ K).abs().max() > .1
    J = polar_corrected(geo, field)(x)
    assert values(geo, polar_corrected(geo, field), x).max().sqrt() < 1e-12
    P = geo.manifold.projector(x)
    for i in range(len(x)):
        S = (-K[i] @ K[i] + torch.eye(9, dtype=DTYPE) - P[i]).numpy()
        ref = K[i].numpy() @ np.linalg.inv(sqrtm(S))
        np.testing.assert_allclose(J[i].numpy(), ref, atol=1e-12, rtol=1e-12)


def test_all_shipped_configs_validate():
    from pathlib import Path
    root = Path(__file__).resolve().parents[1] / 'configs' / 'ambient' / 's4xs1xs1_acs'
    paths = list(root.glob('*.yaml'))
    assert paths
    for path in paths:
        read_config(path)


def test_short_audit_never_claims_final_success_even_for_exact_fixture(tmp_path):
    from extra.validate_ambient_acs import audit
    runner = Runner(tiny(), tmp_path / 'fixture')
    oracle_coefficients(runner.model.poly)
    source = tmp_path / 'oracle.pt'
    torch.save(runner.payload(), source)
    before = source.read_bytes()
    report = audit(source, quick=True, seed=812345)
    assert report['raw_combined']['max_found'] < 1e-12
    assert report['corrected']['max_found'] < 1e-12
    assert report['corrected']['tangency_max'] < 1e-12
    assert report['corrected']['scipy_sqrtm_max_difference'] < 1e-12
    assert report['tangent_spectrum_at_worst']['excluded_normal_modes'] == 3
    assert report['tangent_spectrum_at_worst']['minimum_tangent_singular_value'] > .999
    assert report['screening_targets_passed']
    assert not report['numerically_validated']
    assert not report['certified_exact_acs']
    assert source.read_bytes() == before


def test_lbfgs_rollback_and_real_iteration_count(tmp_path, monkeypatch):
    runner = Runner(tiny(lbfgs_steps=3, lbfgs_chunk=3), tmp_path)
    oracle_coefficients(runner.model.poly)
    with torch.no_grad():
        runner.model.poly.lin.weight.mul_(1.01)
    runner.observe()
    runner.phase = 'lbfgs'
    runner.set_optimizer()
    before = copy.deepcopy(runner.model.state_dict())
    def bad_step(opt, closure):
        closure()
        with torch.no_grad():
            next(runner.model.parameters()).mul_(5)
        opt.state[next(iter(opt.param_groups[0]['params']))]['n_iter'] = 2
    monkeypatch.setattr(torch.optim.LBFGS, 'step', bad_step)
    runner.lbfgs_step()
    assert not runner.lbfgs_log[-1]['accepted']
    assert runner.counts['lbfgs'] == 2
    for key, value in before.items():
        torch.testing.assert_close(value, runner.model.state_dict()[key], atol=0, rtol=0)


def test_lbfgs_runs_and_preserves_an_improvement(tmp_path):
    runner = Runner(tiny(lbfgs_steps=3, lbfgs_chunk=3), tmp_path)
    oracle_coefficients(runner.model.poly)
    with torch.no_grad():
        runner.model.poly.lin.weight.mul_(1.01)
    runner.observe()
    before = runner.best['metrics']['max_found']
    runner.phase = 'lbfgs'
    runner.set_optimizer()
    runner.lbfgs_step()
    assert 0 < runner.counts['lbfgs'] <= 3
    assert runner.best['metrics']['max_found'] < before
    assert runner.lbfgs_log[-1]['accepted']


def test_interrupted_ascent_rolls_back_rng_and_retries_pending_observation(tmp_path, monkeypatch):
    from extra import acs_training
    cfg = tiny(audit_every=2)
    full = Runner(cfg, tmp_path / 'full')
    full.run()
    interrupted = Runner(cfg, tmp_path / 'interrupted')
    def stop_ascent(*args, **kwargs):
        raise acs_training.Stopped()
    with monkeypatch.context() as patch:
        patch.setattr(acs_training, 'ascent', stop_ascent)
        interrupted.run()
    assert interrupted.pending_observation == {'audit': True}
    resumed = Runner(cfg, interrupted.out, resume=True)
    resumed.run()
    for key, value in full.model.state_dict().items():
        torch.testing.assert_close(value, resumed.model.state_dict()[key], atol=0, rtol=0)
    assert full.best['metrics'] == resumed.best['metrics']
    assert torch.equal(full.generators['audit'].get_state(), resumed.generators['audit'].get_state())


def test_nan_loss_is_rejected_before_parameter_update(tmp_path):
    runner = Runner(tiny(), tmp_path)
    before = copy.deepcopy(runner.model.state_dict())
    runner.objective = lambda *args: torch.tensor(float('nan'), dtype=DTYPE)
    with pytest.raises(FloatingPointError, match='Nonfinite ACS loss'):
        runner.adam_step()
    for key, value in before.items():
        torch.testing.assert_close(value, runner.model.state_dict()[key], atol=0, rtol=0)


def test_pbs_training_wrapper_with_mock_environment(tmp_path):
    """Real bash, stub module/conda/Python: tests paths, quoting and pipefail."""
    import os
    from pathlib import Path
    import subprocess
    repo = Path(__file__).resolve().parents[1]
    script = repo / 'scratch' / 'submit_S4xS1xS1_ACS_v2.sh'
    if not script.exists():
        pytest.skip('scratch scripts are intentionally gitignored; use the documented local copy')
    work = tmp_path / 'parent with spaces'
    work.mkdir()
    (work / 'github').symlink_to(repo, target_is_directory=True)
    binaries = tmp_path / 'bin'
    binaries.mkdir()
    for name in ('module', 'conda'):
        stub = binaries / name
        stub.write_text('#!/bin/bash\nexit 0\n')
        stub.chmod(0o755)
    stub = binaries / 'python3'
    stub.write_text('#!/bin/bash\nif [[ "$1" == "-c" ]]; then exit 0; fi\nprintf "%s\\n" "$@"\nexit "${FAKE_EXIT:-0}"\n')
    stub.chmod(0o755)
    env = {**os.environ, 'PATH': str(binaries) + ':' + os.environ['PATH'],
           'PBS_O_WORKDIR': str(work), 'PBS_JOBID': '123.server',
           'RECIPE': 'pilot_trilin_tail', 'SEED': '8'}
    result = subprocess.run(['bash', str(script)], env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert '--config\nconfigs/ambient/s4xs1xs1_acs/pilot_trilin_tail.yaml' in result.stdout
    assert str(work / 'outputs/s4xs1xs1_acs/pilots/pilot_trilin_tail/s8/123.server') in result.stdout
    result = subprocess.run(['bash', str(script)], env={**env, 'FAKE_EXIT': '75'}, capture_output=True)
    assert result.returncode == 75
    result = subprocess.run(['bash', str(script)], env={**env, 'RECIPE': '../bad'}, capture_output=True)
    assert result.returncode == 2
