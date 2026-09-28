"""Real AD/optimizer/campaign tests, plus PBS mocks (no remote submissions)."""
import copy
from dataclasses import asdict
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import torch

from extra import cs_precision as campaign
from extra import cs_precision_inputs as inputs
from extra.cs_precision import Config, Budget
from extra.cs_precision_fields import (six_torus_control, tensor, crosscheck, Correction,
                                       payload, load_model)
from extra.cs_precision_optim import Parameters, gn_step, weights, cg
from extra.campaign_unit import Geometry, load_acs_field, polar_corrected, make_field
from extra.cs_pilot import checked_base
from extra.acs_training import sha256, atomic_json, atomic_torch_save, Stopped
from extra.test_cs_training import setup
from extra.test_acs_training import one_thread
from extra.test_acs_pbs import shell_env, calls

REPO = Path(__file__).resolve().parents[1]


def tiny():
    return Config(rounds=2, chunk=2, batch=4, collocation=4, hard_points=2, validation=4,
                  circle_points=2, search_pool=4, search_starts=2, search_steps=2,
                  adam_steps=2, lbfgs_steps=2, gn_steps=2, eval_every=1, refresh_every=1,
                  cg_steps=2, correction_hidden=[4], audit_uniform=4, audit_pool=4,
                  audit_starts=2, audit_steps=2, audit_circle_points=2).validate()


def source_fixture(tmp_path, monkeypatch):
    torch.manual_seed(7123)
    _, base, audit, previous = setup(tmp_path)
    target_audit = base/'final.audit.full.seed900001.json'
    shutil.copy(audit, target_audit)
    monkeypatch.setattr(inputs, 'ACS_SHA256', sha256(base/'final.pt'))
    source = tmp_path/'long_run'
    source.mkdir()
    evidence = checked_base(base, target_audit)
    ck = torch.load(previous, weights_only=False)
    ck['audited_base'] = evidence
    atomic_torch_save(ck, previous)
    geo = Geometry('s4xs1xs1', 'ambient')
    for name in ('tail', 'average'):
        shutil.copy(previous, source/f'final_{name}.pt')
        folder = source/f'diagnostic_{name}'
        folder.mkdir()
        points = geo.sample(4)
        atomic_torch_save({'points': points, 'metrics': {'nij': torch.arange(4, dtype=torch.float64)},
                          'status': 'completed_diagnostic', 'checkpoint_sha256': sha256(source/f'final_{name}.pt')},
                          folder/'diagnostic_points.pt')
        atomic_json({'status': 'completed_diagnostic', 'sha256': sha256(source/f'final_{name}.pt'),
                     'points_sha256': sha256(folder/'diagnostic_points.pt'), 'audited_base': evidence}, folder/'report.json')
    atomic_json({'status': 'completed_pending_independent_audit', 'steps_done': 50000,
                 'ready_for_diagnostics': True, 'source': {'path': str(previous), 'sha256': sha256(previous)}},
                source/'result.json')
    return source, base, target_audit


def prepare_fixture(tmp_path, monkeypatch, real=False):
    source, base, audit = source_fixture(tmp_path, monkeypatch)
    cfg, root = tiny(), tmp_path/'campaign'
    root.mkdir()
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip()
    if real:
        campaign.prepare(root, source, base, cfg, sha256(source/'final_tail.pt'), revision, Budget(cfg))
    else:
        geo = Geometry('s4xs1xs1', 'ambient')
        field = polar_corrected(geo, load_acs_field(geo, base))
        model, _, fn = campaign.restore(source/'final_tail.pt', geo, field, checked_base(base, audit))
        probes = {'uniform': geo.sample(4), 'challenge': geo.sample(4)}
        m = campaign.measure(geo, model, fn, probes, cfg, Budget(cfg))
        atomic_torch_save(probes, root/'probes.pt')
        atomic_torch_save(geo.sample(4), root/'training_witnesses.pt')
        campaign.receipt(root, 'prepare', {'status': 'complete', 'config': asdict(cfg), 'code_hash': campaign.code_hash(),
             'inputs': inputs.verify(source, base), 'metrics': {'tail': m}},
             [root/'probes.pt', root/'training_witnesses.pt'])
    return cfg, root, source, base, audit


def test_coordinate_and_finite_difference_controls_detect_wrong_tensor():
    geo, field = six_torus_control()
    points = geo.manifold.sample(2, generator=torch.Generator().manual_seed(22))
    assert tensor(geo, field, points).norm() < 1e-12
    assert crosscheck(geo, field, points)['passed']
    def negative(x):
        P = geo.manifold.projector(x)
        A = torch.zeros((12, 12), dtype=torch.float64)
        A[0, 0] = .3
        B = P @ (x[:, 4, None, None]*A) @ P
        return torch.matrix_exp(B) @ field(x) @ torch.matrix_exp(-B)
    assert tensor(geo, negative, points).norm() > .01
    assert crosscheck(geo, negative, points)['passed']
    from extra import cs_precision_fields
    with pytest.MonkeyPatch.context() as patch:
        original = cs_precision_fields.tensor
        patch.setattr(cs_precision_fields, 'tensor', lambda *a: original(*a)*0)
        assert not crosscheck(geo, negative, points)['passed']


def test_correction_exact_start_spatial_derivatives_and_roundtrip(tmp_path):
    _, base, audit, source = setup(tmp_path)
    geo = Geometry('s4xs1xs1', 'ambient')
    fixed = polar_corrected(geo, load_acs_field(geo, base))
    model, arch, field = campaign.restore(source, geo, fixed, checked_base(base, audit))
    new = Correction(model, fixed, geo, [4])
    corrected = make_field(geo, 'cs', new, base_fn=fixed)[0]
    x = geo.sample(2)
    torch.testing.assert_close(corrected(x), field(x), atol=1e-13, rtol=0)
    torch.testing.assert_close(tensor(geo, corrected, x), tensor(geo, field, x), atol=1e-12, rtol=0)
    tensor(geo, corrected, x).square().sum().backward()
    assert new.head.net[-1].weight.grad.norm() > 0
    assert all(p.grad is None for p in new.parent.parameters())
    optimizer = torch.optim.Adam(new.head.parameters(), lr=1e-4)
    optimizer.step()
    J = corrected(x)
    assert (J@J+geo.manifold.projector(x)).norm() < 1e-11
    spec = {'kind': 'correction', 'parent': arch, 'hidden': [4]}
    ck = payload(new, spec, checked_base(base, audit), {})
    restored, _ = load_model(ck, geo, fixed)
    torch.testing.assert_close(restored(x), new(x), atol=0, rtol=0)
    assert all(not p.requires_grad for p in restored.parent.parameters())


def test_gn_real_tensor_derivatives_and_reduction(tmp_path):
    _, base, audit, source = setup(tmp_path)
    geo = Geometry('s4xs1xs1', 'ambient')
    fixed = polar_corrected(geo, load_acs_field(geo, base))
    model, _, field = campaign.restore(source, geo, fixed, checked_base(base, audit))
    params = Parameters(model)
    residual = campaign.functional_residual(model, params, geo, fixed)
    v, x = params.vector(), geo.sample(2)
    torch.testing.assert_close(residual(v, x), tensor(geo, field, x).flatten(1))
    direction = torch.randn_like(v); direction /= direction.norm()
    _, jvp = torch.func.jvp(lambda z: residual(z, x), (v,), (direction,))
    fd = (residual(v+1e-5*direction, x)-residual(v-1e-5*direction, x))/(2e-5)
    torch.testing.assert_close(jvp, fd, atol=1e-7, rtol=1e-5)
    _, pullback = torch.func.vjp(lambda z: residual(z, x), v)
    witness = torch.randn_like(jvp)
    torch.testing.assert_close((jvp*witness).sum(), pullback(witness)[0].dot(direction))
    after, info = gn_step(residual, v, x, torch.ones(2, dtype=torch.float64)/2,
                          damping=.1, chunk=1, cg_steps=2, check=lambda: None)
    assert info['accepted'] and info['after'] < info['before']
    torch.testing.assert_close(params.vector(), v, atol=0, rtol=0)  # transactional
    assert torch.isfinite(after).all()


def test_weight_scaling_cg_and_interruption():
    scores = torch.tensor([1., 2., 3., 4.], dtype=torch.float64)
    w = weights(scores, 3, .5, .25)
    torch.testing.assert_close(w, torch.tensor([1/3, 1/3, 1/3, .5], dtype=torch.float64))
    rhs = torch.tensor([1., 2.], dtype=torch.float64)
    answer, count = cg(lambda v: 2*v, rhs, 4, .01, lambda: None)
    torch.testing.assert_close(answer, rhs/2)
    assert count == 1
    with pytest.raises(FloatingPointError):
        cg(lambda v: -v, rhs, 2, .01, lambda: None)
    def stop():
        raise Stopped()
    with pytest.raises(Stopped):
        gn_step(lambda v, x: v.expand(len(x), 2), rhs, torch.ones(2, 1), torch.ones(2),
                damping=1., chunk=1, cg_steps=2, check=stop)


def test_real_prepare_all_pilots_selection_continuation_and_audit(tmp_path, monkeypatch):
    cfg, root, source, base, audit = prepare_fixture(tmp_path, monkeypatch, real=True)
    original = {str(p): sha256(p) for p in (source/'final_tail.pt', source/'final_average.pt', base/'final.pt')}
    for method in campaign.METHODS:
        result = campaign.train_block(root, f'pilot_{method}', source/'final_tail.pt', method, cfg, Budget(cfg))
        assert 1 <= result['steps'] <= 2 and result['status'] == 'complete'
        assert result['metrics']['valid']
    decision = campaign.decide(root, 'select', None, [f'pilot_{m}' for m in campaign.METHODS], cfg, Budget(cfg))
    assert Path(decision['checkpoint']).is_file()
    for index in range(1, cfg.rounds+1):
        campaign.continuation(root, index, cfg, Budget(cfg))
    for part in range(4):
        report = campaign.audit(root, part, cfg, Budget(cfg))
        assert report['crosscheck']['passed']
    report = campaign.summary(root, cfg)
    assert report['all_jobs_valid'] and not report['integrability_certified']
    assert {p: sha256(p) for p in original} == original
    with pytest.raises(ValueError, match='already exists'):
        campaign.continuation(root, 1, cfg, Budget(cfg))
    final = campaign.required(root, f'round_{cfg.rounds}')
    Path(final['checkpoint']).write_bytes(b'corrupted')
    assert not campaign.summary(root, cfg)['all_jobs_valid']


def test_failure_gates_changed_inputs_config_and_receipts(tmp_path, monkeypatch):
    cfg, root, source, base, audit = prepare_fixture(tmp_path, monkeypatch)
    with pytest.raises(FileNotFoundError):
        campaign.continuation(root, 1, cfg, Budget(cfg))
    with pytest.raises(ValueError, match='Configuration changed'):
        other = copy.deepcopy(cfg); other.lr_min = 1e-8
        campaign.context(root, Budget(other))
    with pytest.raises(ValueError, match='Code changed'):
        with monkeypatch.context() as patch:
            patch.setattr(campaign, 'code_hash', lambda: 'changed')
            campaign.context(root, Budget(cfg))
    (root/'training_witnesses.pt').write_bytes(b'bad')
    with pytest.raises(ValueError, match='artifact changed'):
        campaign.context(root, Budget(cfg))
    (source/'final_tail.pt').write_bytes(b'changed')
    with pytest.raises(ValueError, match='provenance mismatch'):
        inputs.verify(source, base)


def test_review_stop_propagates_without_training(tmp_path, monkeypatch):
    cfg, root, source, base, audit = prepare_fixture(tmp_path, monkeypatch)
    campaign.receipt(root, 'select', {'status': 'review', 'checkpoint': str(source/'final_tail.pt'),
                    'monitor_file': str(root/'probes.pt'), 'round': 0}, [source/'final_tail.pt', root/'probes.pt'])
    monkeypatch.setattr(campaign, 'train_block', lambda *a, **kw: pytest.fail('must not train after review'))
    for index in range(1, cfg.rounds+1):
        assert campaign.continuation(root, index, cfg, Budget(cfg))['status'] == 'review'


def test_real_correction_gn_block_preserves_parent(tmp_path, monkeypatch):
    cfg, root, source, base, audit = prepare_fixture(tmp_path, monkeypatch)
    cfg.gn_steps = 1
    # Update this synthetic preparation to the exact test configuration.
    prepared = campaign.required(root, 'prepare')
    prepared['config'] = asdict(cfg)
    atomic_json(prepared, root/'prepare.json')
    digest = sha256(source/'final_tail.pt')
    report = campaign.train_block(root, 'correction_test', source/'final_tail.pt', 'gn', cfg,
                                   Budget(cfg), round_index=1, correction=True)
    assert report['steps'] == 1 and report['metrics']['valid']
    assert sha256(source/'final_tail.pt') == digest
    ck = torch.load(report['checkpoint'], weights_only=False)
    assert ck['architecture']['kind'] == 'correction'


def test_configuration_and_invalid_ranges():
    cfg = campaign.read_config(REPO/'configs/ambient/s4xs1xs1_cs/precision.yaml')
    assert cfg.lr_min == 1e-7 and cfg.audit_uniform == 100000 and cfg.rounds == 6
    for key, value in [('lr_min', 1.), ('batch', True), ('work_minutes', 999), ('minimum_gain', 2.), ('rounds', 0)]:
        bad = copy.deepcopy(cfg); setattr(bad, key, value)
        with pytest.raises(ValueError):
            bad.validate()


def test_summary_missing_work_never_claims_audit_or_targets(tmp_path):
    report = campaign.summary(tmp_path, Config().validate())
    assert not report['all_jobs_valid']
    assert not report['full_audit_protocol']
    assert not report['sampled_targets_passed']
    assert not report['integrability_certified']


def test_six_dimensional_integrable_control_recovered_by_gn():
    geo, base = six_torus_control()
    rng = torch.Generator().manual_seed(712)
    x = geo.manifold.sample(4, generator=rng)
    A = torch.zeros((12, 12), dtype=torch.float64); A[0, 0] = 1.
    def residual(v, points):
        def field(z):
            P = geo.manifold.projector(z)
            B = P @ (v[0]*z[:, 4, None, None]*A) @ P
            return torch.matrix_exp(B) @ base(z) @ torch.matrix_exp(-B)
        return tensor(geo, field, points).flatten(1)
    v = torch.tensor([.1], dtype=torch.float64)
    for _ in range(4):
        v, info = gn_step(residual, v, x, torch.ones(4, dtype=torch.float64)/4,
                          damping=1e-5, chunk=2, cg_steps=2, check=lambda: None)
    fresh = geo.manifold.sample(16, generator=rng)
    assert residual(v, fresh).norm(dim=1).max() < 1e-6


@pytest.mark.parametrize('failure', [None, 'FAIL_CALL', 'EMPTY_CALL'])
def test_launcher_dependencies_and_partial_submission(shell_env, tmp_path, failure):
    script = REPO/'scratch/launch_S4xS1xS1_CS_precision.sh'
    if not script.exists():
        pytest.skip('scratch intentionally ignored')
    work = tmp_path/'submit'; work.mkdir()
    for name in (script.name, 'submit_S4xS1xS1_CS_precision.sh'):
        shutil.copy(REPO/'scratch'/name, work/name)
    python = Path(shell_env['PATH'].split(':')[0])/'python3'
    python.write_text(f'#!{sys.executable}\n' + '''import json,os,sys
if sys.argv[1].endswith('cs_precision_inputs.py'):
    source=sys.argv[sys.argv.index('--source-run')+1]
    print(json.dumps({'source_run':source, 'files':{source+'/final_tail.pt':'a'*64}}))
else:
    os.execv(sys.executable, [sys.executable]+sys.argv[1:])
''')
    env = dict(shell_env)
    if failure:
        env[failure] = '3'
    result = subprocess.run(['bash', str(work/script.name)], cwd=work, env=env, capture_output=True, text=True)
    assert result.returncode == (1 if failure else 0), result.stderr
    submitted = calls(env)
    assert len(submitted) == (3 if failure else 16)
    if not failure:
        assert submitted[4][submitted[4].index('-W')+1] == 'depend=afterany:102.server:103.server:104.server'
        assert submitted[-1][submitted[-1].index('-W')+1] == 'depend=afterany:112.server:113.server:114.server:115.server'
        for args in submitted[11:15]:
            assert args[args.index('-W')+1] == 'depend=afterany:111.server'


def test_pbs_arguments_and_failure_exit_code(shell_env, tmp_path):
    script = REPO/'scratch/submit_S4xS1xS1_CS_precision.sh'
    if not script.exists():
        pytest.skip('scratch intentionally ignored')
    python = Path(shell_env['PATH'].split(':')[0])/'python3'
    python.write_text('#!/bin/bash\nprintf "%s\\n" "$@"\nexit "${FAKE_EXIT:-0}"\n')
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip()
    env = {**shell_env, 'ROOT': str(tmp_path/'output'), 'PHASE': 'continue', 'INDEX': '2',
           'CONFIG': '/absolute/config.yaml', 'REVISION': revision}
    result = subprocess.run(['bash', str(script)], env=env, capture_output=True, text=True)
    assert result.returncode == 0 and '--index\n2' in result.stdout
    result = subprocess.run(['bash', str(script)], env={**env, 'FAKE_EXIT': '75'}, capture_output=True, text=True)
    assert result.returncode == 75
