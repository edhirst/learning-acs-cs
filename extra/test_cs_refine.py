"""Real small-geometry refinement tests and mocked HPC submission, no remote work."""
import copy
from dataclasses import asdict
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import torch
import yaml

from extra import cs_refine as r
from extra import cs_refine_inputs as inputs
from extra import cs_precision as p
from extra.acs_training import atomic_json, atomic_torch_save, sha256, Stopped
from extra.test_cs_precision import source_fixture
from extra.test_acs_training import one_thread
from extra.test_acs_pbs import shell_env, calls

REPO = Path(__file__).resolve().parents[1]


def tiny(**kwargs):
    values = dict(blocks=2, steps_per_block=2, schedule_steps=10, chunk=2, batch=4,
                  replay_fraction=.25, eval_every=1, refresh_every=2, monitor_uniform=4,
                  monitor_challenge=2, validation=4, circle_points=1, search_pool=4,
                  search_starts=2, search_steps=2, hard_points=2, correction_hidden=[4],
                  audit_uniform=4, audit_pool=4, audit_starts=2, audit_steps=2,
                  audit_circle_points=1)
    return r.Config(**{**values, **kwargs}).validate()


def fixture(tmp_path, monkeypatch, cfg=None):
    source, base, audit = source_fixture(tmp_path, monkeypatch)
    from extra.cs_precision_inputs import verify
    data = verify(source, base)
    data.update(checkpoint=str(source/'final_tail.pt'), checkpoint_sha256=sha256(source/'final_tail.pt'),
                clouds=[str(source/'diagnostic_tail/diagnostic_points.pt')], source_root=str(source))
    monkeypatch.setattr(r, 'verify', lambda _: data)
    cfg = cfg or tiny()
    root = tmp_path/'refine'
    root.mkdir()
    (root/'config.yaml').write_text(yaml.safe_dump(asdict(cfg)))
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip()
    r.prepare(root, cfg, source, data['checkpoint_sha256'], revision, p.Budget(cfg.protocol()))
    return cfg, root


def test_real_branches_resume_select_and_audit(tmp_path, monkeypatch):
    cfg, root = fixture(tmp_path, monkeypatch)
    for branch in r.BRANCHES:
        first = r.train(root, cfg, branch, 1, p.Budget(cfg.protocol()))
        assert first['completed'] == 2 and not first['optimizer_resumed']
        second = r.train(root, cfg, branch, 2, p.Budget(cfg.protocol()))
        assert second['completed'] == 4 and second['optimizer_resumed']
        assert set(second['nominations']) == set(r.ROLES)
        assert set(second['full_nominations']) == set(r.ROLES)
        anchor = p.required(root, 'prepare')['source_metrics']
        for role in r.ROLES:
            assert r.nomination_score(role, second['full_nominations'][role]['metrics'], anchor) <= (
                r.nomination_score(role, first['full_nominations'][role]['metrics'], anchor) + 1e-12)
        assert second['tail_weight'] == first['tail_weight']
        ck = torch.load(second['resume'], weights_only=False)
        assert ck['training']['optimizer']['state']
        assert all(v['step'] == 4 for v in ck['training']['optimizer']['state'].values())
        if branch == 'correction':
            initial = torch.load(root/'correction_1/initial.pt', weights_only=False)
            assert ck['architecture']['kind'] == 'correction'
            for key in ck['model']:
                if key.startswith('parent.'):
                    torch.testing.assert_close(ck['model'][key], initial['model'][key], rtol=0, atol=0)
    selection = r.select(root, cfg, p.Budget(cfg.protocol()))
    assert Path(selection['checkpoint']).exists()
    assert len(selection['comparisons']) == 13
    for i in range(4):
        assert p.audit(root, i, cfg.protocol(), p.Budget(cfg.protocol()))['crosscheck']['passed']
    summary = r.summary(root, cfg)
    assert summary['all_jobs_valid']
    assert not summary['full_audit_protocol']
    assert not summary['sampled_targets_passed']
    assert not summary['integrability_certified']
    Path(selection['checkpoint']).write_bytes(b'bad')
    assert not r.summary(root, cfg)['all_jobs_valid']


def test_exact_optimizer_sampling_schedule_across_allocations(tmp_path, monkeypatch):
    cfg, root = fixture(tmp_path, monkeypatch, tiny(steps_per_block=2, eval_every=2, refresh_every=2))
    fullroot = tmp_path/'uninterrupted'
    shutil.copytree(root, fullroot)
    fullcfg = copy.deepcopy(cfg)
    fullcfg.steps_per_block = 4
    (fullroot/'config.yaml').write_text(yaml.safe_dump(asdict(fullcfg)))
    # Bind the copied fixture to its own config/artifacts, as a fresh campaign.
    manifest = json.loads((fullroot/'prepare.json').read_text())
    manifest['refinement_config'] = asdict(fullcfg)
    manifest['artifacts'] = {str(fullroot/Path(path).name): sha256(fullroot/Path(path).name)
                             for path in manifest['artifacts']}
    atomic_json(manifest, fullroot/'prepare.json')
    r.train(root, cfg, 'low', 1, p.Budget(cfg.protocol()))
    split = r.train(root, cfg, 'low', 2, p.Budget(cfg.protocol()))
    full = r.train(fullroot, fullcfg, 'low', 1, p.Budget(fullcfg.protocol()))
    a, b = [torch.load(v['resume'], weights_only=False) for v in (split, full)]
    for key in a['model']:
        torch.testing.assert_close(a['model'][key], b['model'][key], atol=0, rtol=0)
    for key in ('rng', 'torch_rng', 'hard'):
        torch.testing.assert_close(a['training'][key], b['training'][key], atol=0, rtol=0)
    assert a['training']['completed'] == b['training']['completed'] == 4
    for key, state in a['training']['optimizer']['state'].items():
        for name, val in state.items():
            torch.testing.assert_close(val, b['training']['optimizer']['state'][key][name], atol=0, rtol=0)
    assert a['training']['optimizer']['param_groups'] == b['training']['optimizer']['param_groups']


def test_tradeoff_nominations_are_independent_and_schedule_persists():
    initial = {'rms': .286719, 'tail': .747038}
    tradeoff = {'rms': .289626, 'tail': .600297}
    assert r.nomination_score('tail', tradeoff, initial) < r.nomination_score('tail', initial, initial)
    assert r.nomination_score('rms', tradeoff, initial) > r.nomination_score('rms', initial, initial)
    cfg = r.Config().validate()
    assert r.learning_rate(cfg, 0) == cfg.lr
    assert r.learning_rate(cfg, 2000) < r.learning_rate(cfg, 1000)
    assert r.learning_rate(cfg, cfg.schedule_steps-1) == pytest.approx(cfg.lr_min)
    assert r.learning_rate(cfg, cfg.schedule_steps+100) == pytest.approx(cfg.lr_min)


def test_direct_backward_matches_existing_functional_tensor(tmp_path, monkeypatch):
    source, base_dir, audit = source_fixture(tmp_path, monkeypatch)
    geo = p.Geometry('s4xs1xs1', 'ambient')
    base = p.polar_corrected(geo, p.load_acs_field(geo, base_dir))
    model, _, field = p.restore(source/'final_tail.pt', geo, base, p.checked_base(base_dir, audit))
    from extra.cs_precision_optim import Parameters, backward_objective
    params = Parameters(model)
    x = geo.sample(2)
    w = torch.tensor([.5, .7], dtype=p.DTYPE)
    residual = p.functional_residual(model, params, geo, base)
    old = backward_objective(residual, params, x, w, 1, lambda: None)
    expected = {n: v.grad.clone() for n, v in params.items}
    model.zero_grad(set_to_none=True)
    direct = (r.tensor(geo, field, x).flatten(1).square().sum(1)*w).sum()
    direct.backward()
    assert direct.item() == pytest.approx(old, abs=1e-12)
    for name, value in params.items:
        torch.testing.assert_close(value.grad, expected[name], atol=1e-10, rtol=1e-10)


def test_time_limit_saves_completed_iterate_but_signal_fails_closed(tmp_path, monkeypatch):
    cfg, root = fixture(tmp_path, monkeypatch, tiny(steps_per_block=5, eval_every=1))
    class WorkBudget(p.Budget):
        def check(self, work=False):
            if work and (root/'low_1/progress.jsonl').exists():
                raise Stopped()
            return super().check(work)
    report = r.train(root, cfg, 'low', 1, WorkBudget(cfg.protocol()))
    assert report['steps'] == 1 and report['reason'] == 'work_time_budget'
    ck = torch.load(report['resume'], weights_only=False)
    assert ck['training']['completed'] == 1
    class SignalBudget(WorkBudget):
        def check(self, work=False):
            if work:
                self.signalled = True
                raise Stopped()
            return p.Budget.check(self, work)
    with pytest.raises(Stopped):
        r.train(root, cfg, 'tail', 1, SignalBudget(cfg.protocol()))
    assert not (root/'tail_1.json').exists()
    assert (root/'tail_1/last.pt').exists()


def test_partial_batch_budget_restores_sampling_before_unobserved_update(tmp_path, monkeypatch):
    cfg, root = fixture(tmp_path, monkeypatch, tiny(steps_per_block=2, eval_every=2, refresh_every=100))
    control = tmp_path/'control'
    shutil.copytree(root, control)
    one = copy.deepcopy(cfg); one.steps_per_block = 1
    (control/'config.yaml').write_text(yaml.safe_dump(asdict(one)))
    manifest = json.loads((control/'prepare.json').read_text())
    manifest['refinement_config'] = asdict(one)
    manifest['artifacts'] = {str(control/Path(path).name): sha256(control/Path(path).name)
                             for path in manifest['artifacts']}
    atomic_json(manifest, control/'prepare.json')
    class PartialBudget(p.Budget):
        count = 0
        def check(self, work=False):
            if work:
                self.count += 1
                # Six checks complete update 1; update 2 has drawn its random
                # batch by check 9. There was no observation at update 1.
                if self.count == 9:
                    raise Stopped()
            return super().check(work)
    partial = r.train(root, cfg, 'low', 1, PartialBudget(cfg.protocol()))
    baseline = r.train(control, one, 'low', 1, p.Budget(one.protocol()))
    assert partial['completed'] == baseline['completed'] == 1
    a, b = [torch.load(v['resume'], weights_only=False) for v in (partial, baseline)]
    for key in a['model']:
        torch.testing.assert_close(a['model'][key], b['model'][key], atol=0, rtol=0)
    for key in ('rng', 'torch_rng', 'hard'):
        torch.testing.assert_close(a['training'][key], b['training'][key], atol=0, rtol=0)


def test_numerical_gate_rolls_back_and_propagates_review(tmp_path, monkeypatch):
    cfg, root = fixture(tmp_path, monkeypatch)
    original = p.measure
    count = 0
    def measure(*args, **kwargs):
        nonlocal count
        result = original(*args, **kwargs)
        count += 1
        if count == 2:
            result['valid'] = False
        return result
    monkeypatch.setattr(p, 'measure', measure)
    first = r.train(root, cfg, 'low', 1, p.Budget(cfg.protocol()))
    assert first['status'] == 'review' and first['completed'] == 0
    saved = torch.load(first['resume'], weights_only=False)
    initial = torch.load(root/'low_1/initial.pt', weights_only=False)
    for key in saved['model']:
        torch.testing.assert_close(saved['model'][key], initial['model'][key], atol=0, rtol=0)
    second = r.train(root, cfg, 'low', 2, p.Budget(cfg.protocol()))
    assert second['status'] == 'review' and 'skipped' in second['reason']
    assert not (root/'low_2').exists()


def test_changed_config_inputs_duplicate_and_missing_receipts(tmp_path, monkeypatch):
    cfg, root = fixture(tmp_path, monkeypatch)
    with pytest.raises(FileNotFoundError):
        r.train(root, cfg, 'low', 2, p.Budget(cfg.protocol()))
    with pytest.raises(FileNotFoundError):
        r.select(root, cfg, p.Budget(cfg.protocol()))
    other = copy.deepcopy(cfg); other.lr = 2e-6
    with pytest.raises(ValueError, match='configuration changed'):
        r.context(root, other, p.Budget(other.protocol()))
    result = r.summary(root, cfg)
    assert not result['all_jobs_valid'] and not result['sampled_targets_passed']
    r.train(root, cfg, 'low', 1, p.Budget(cfg.protocol()))
    with pytest.raises(ValueError, match='already exists'):
        r.train(root, cfg, 'low', 1, p.Budget(cfg.protocol()))
    (root/'config.yaml').write_text('bad')
    with pytest.raises(ValueError, match='artifact changed'):
        r.context(root, cfg, p.Budget(cfg.protocol()))


@pytest.mark.parametrize('change', [{'batch': 1}, {'replay_fraction': .001}, {'lr': float('nan')},
                                  {'work_minutes': 220.}, {'monitor_uniform': 100000},
                                  {'correction_hidden': []}, {'tail_fraction': 2.}])
def test_strict_config(change):
    with pytest.raises(ValueError):
        r.Config(**change).validate()


def test_stdlib_source_preflight_hash_and_audit_validation(tmp_path):
    # Metadata fixture: deliberately not a numerical audit/convergence test.
    root = tmp_path/'old'; root.mkdir()
    checkpoint = root/'field.pt'; checkpoint.write_bytes(b'fixture')
    oldfile = root/'base.pt'; oldfile.write_bytes(b'base')
    prepared = {'status': 'complete', 'inputs': {'files': {str(oldfile): sha256(oldfile)},
                                               'source_run': '/old', 'base_dir': '/base', 'base_audit': '/audit'}}
    atomic_json(prepared, root/'prepare.json')
    selected = {'checkpoint': str(checkpoint), 'artifacts': {str(checkpoint): sha256(checkpoint)}}
    audits = []
    for i in range(4):
        cloud = root/f'audit_{i}.points.pt'; cloud.write_bytes(b'cloud')
        audit = {'status': 'complete', 'sha256': sha256(checkpoint), 'checkpoint': str(checkpoint),
                 'crosscheck': {'passed': True}, 'growth_passed': True,
                 'combined': {'acs_max': 1e-14, 'tangency_max': 1e-14, 'count': 16896},
                 'groups': {'uniform': {'count': 100000}, 'structured': {'count': 20480}},
                 'artifacts': {str(cloud): sha256(cloud), str(checkpoint): sha256(checkpoint)}}
        atomic_json(audit, root/f'audit_{i}.json'); audits.append(audit)
    comp = {'all_jobs_valid': True, 'full_audit_protocol': True, 'errors': {}, 'selected': selected, 'audits': audits}
    atomic_json(comp, root/'comparison.json')
    assert inputs.verify(root)['checkpoint_sha256'] == sha256(checkpoint)
    (root/'audit_3.points.pt').write_bytes(b'corrupt')
    with pytest.raises(ValueError, match='hash mismatch'):
        inputs.verify(root)


@pytest.mark.parametrize('failure', [None, 'FAIL_CALL', 'EMPTY_CALL'])
def test_launcher_builds_resumable_dependency_chains(shell_env, tmp_path, failure):
    script = REPO/'scratch/launch_S4xS1xS1_CS_refine.sh'
    if not script.exists():
        pytest.skip('scratch intentionally ignored')
    work = tmp_path/'submit'; work.mkdir()
    for name in (script.name, 'submit_S4xS1xS1_CS_refine.sh'):
        shutil.copy(REPO/'scratch'/name, work/name)
    python = Path(shell_env['PATH'].split(':')[0])/'python3'
    python.write_text(f'#!{sys.executable}\n' + '''import json,os,sys
if sys.argv[1].endswith('cs_refine_inputs.py'):
    print(json.dumps({'checkpoint_sha256':'a'*64}))
else:
    os.execv(sys.executable,[sys.executable]+sys.argv[1:])
''')
    env = dict(shell_env)
    if failure:
        env[failure] = '3'
    result = subprocess.run(['bash', str(work/script.name)], cwd=work, env=env, capture_output=True, text=True)
    assert result.returncode == (1 if failure else 0), result.stderr
    submitted = calls(env)
    assert len(submitted) == (3 if failure else 19)
    if not failure:
        assert submitted[13][submitted[13].index('-W')+1] == 'depend=afterany:105.server:109.server:113.server'
        assert submitted[-1][submitted[-1].index('-W')+1] == 'depend=afterany:115.server:116.server:117.server:118.server'
        for i in (1, 5, 9):
            assert submitted[i][submitted[i].index('-W')+1] == 'depend=afterany:101.server'


def test_wrapper_preserves_failure_code_and_arguments(shell_env, tmp_path):
    script = REPO/'scratch/submit_S4xS1xS1_CS_refine.sh'
    if not script.exists():
        pytest.skip('scratch intentionally ignored')
    python = Path(shell_env['PATH'].split(':')[0])/'python3'
    python.write_text('#!/bin/bash\nprintf "%s\\n" "$@"\nexit "${FAKE_EXIT:-0}"\n')
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip()
    env = {**shell_env, 'ROOT': str(tmp_path/'out'), 'CONFIG': '/absolute/config.yaml',
           'REVISION': revision, 'PHASE': 'train', 'BRANCH': 'correction', 'INDEX': '2'}
    result = subprocess.run(['bash', str(script)], env=env, capture_output=True, text=True)
    assert result.returncode == 0 and '--branch\ncorrection\n--index\n2' in result.stdout
    failed = subprocess.run(['bash', str(script)], env={**env, 'FAKE_EXIT': '75'}, capture_output=True)
    assert failed.returncode == 75
