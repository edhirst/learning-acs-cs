"""Synthetic fixtures, real CS derivatives/training; no long production runs."""
import copy
from dataclasses import asdict
import json
from pathlib import Path
import shutil
import subprocess

import pytest
import torch

from extra.acs_training import Stopped, sha256
from extra.cs_matched import Runner, Sampling, comparison, diverse, load_witnesses, read_config
from extra.test_cs_training import setup
from extra.test_acs_training import one_thread
from extra.test_acs_pbs import shell_env, calls

ROOT = Path(__file__).resolve().parents[1]


def fixture(tmp_path, mode='replay'):
    cfg, base, audit, source = setup(tmp_path)
    cfg['matched'] = {'mode': mode, 'sampling': asdict(Sampling(
        target_count=1, capacity=4, refresh_every=1, pool_size=4, starts=2,
        ascent_steps=2, challenge_pool=4, challenge_starts=2, challenge_steps=2))}
    witness = tmp_path / 'witness'
    witness.mkdir()
    from extra.campaign_unit import Geometry
    points = Geometry('s4xs1xs1').sample(6)
    data = {'points': points, 'groups': ['uniform']*2 + ['search_1']*4,
            'metrics': {'nij': torch.arange(6, dtype=torch.float64)},
            'checkpoint_sha256': sha256(source), 'status': 'completed_diagnostic'}
    torch.save(data, witness / 'diagnostic_points.pt')
    report = {'status': 'completed_diagnostic', 'sha256': sha256(source),
              'points_sha256': sha256(witness / 'diagnostic_points.pt')}
    (witness / 'report.json').write_text(json.dumps(report))
    return cfg, base, audit, source, witness


def assert_nested(a, b):
    if isinstance(a, torch.Tensor):
        torch.testing.assert_close(a, b, rtol=0, atol=0)
    elif isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a:
            assert_nested(a[key], b[key])
    elif isinstance(a, (list, tuple)):
        assert len(a) == len(b)
        for aa, bb in zip(a, b):
            assert_nested(aa, bb)
    else:
        assert a == b


def test_both_arms_real_training_common_validation_and_immutable_inputs(tmp_path):
    cfg, base, audit, source, witness = fixture(tmp_path)
    paths = [source, base / 'final.pt', audit, witness / 'report.json', witness / 'diagnostic_points.pt']
    fingerprints = [sha256(p) for p in paths]
    replay = Runner(cfg, base, audit, tmp_path / 'replay', witness, init_from=source)
    uniform_cfg = {**cfg, 'matched': {**cfg['matched'], 'mode': 'uniform'}}
    uniform = Runner(uniform_cfg, base, audit, tmp_path / 'uniform', witness, init_from=source)
    torch.testing.assert_close(replay.validation, uniform.validation, rtol=0, atol=0)
    for runner in (replay, uniform):
        result = runner.run()
        assert result['status'] == 'completed_pending_independent_audit'
        assert result['dense']['acs_floorfree_sup'] < 1e-12
        assert result['dense']['tangency_max'] < 1e-12
        assert result['best_validation']['selection_score'] <= 1
        assert not result['numerically_validated_cs']
        last = torch.load(runner.out / 'state.pt', weights_only=False)
        initial = torch.load(source, weights_only=False)['model']
        assert any(not torch.equal(last['model'][k], initial[k]) for k in initial)
    assert replay.initial == uniform.initial
    torch.testing.assert_close(replay.challenge, uniform.challenge, rtol=0, atol=0)
    assert torch.equal(replay.train_rng.get_state(), uniform.train_rng.get_state())
    assert replay.last_refresh == 2 and uniform.last_refresh == 0
    assert [sha256(p) for p in paths] == fingerprints


@pytest.mark.parametrize('mode', ['uniform', 'replay'])
def test_exact_resume_including_refresh_selection_and_rng(tmp_path, mode):
    cfg, base, audit, source, witness = fixture(tmp_path, mode)
    full = Runner(cfg, base, audit, tmp_path / 'full', witness, init_from=source)
    full.run()
    split = Runner(cfg, base, audit, tmp_path / 'split', witness, init_from=source)
    assert split.run(max_updates=1)['status'] == 'incomplete'
    resumed = Runner(cfg, base, audit, split.out, witness, resume=True)
    resumed.run()
    a, b = [torch.load(p / 'state.pt', weights_only=False) for p in (full.out, split.out)]
    for key in ('model', 'optimizer', 'scheduler', 'train_rng', 'validation', 'best', 'initial', 'training_extension'):
        assert_nested(a[key], b[key])


def test_sampling_uniform_prefix_jitter_and_diversity(tmp_path):
    cfg, base, audit, source, witness = fixture(tmp_path)
    cfg['batch'] = 8
    cfg['matched']['sampling']['target_count'] = 4
    replay = Runner(cfg, base, audit, tmp_path / 'replay', witness, init_from=source)
    uniform_cfg = {**cfg, 'matched': {**cfg['matched'], 'mode': 'uniform'}}
    uniform = Runner(uniform_cfg, base, audit, tmp_path / 'uniform', witness, init_from=source)
    x, y = replay.training_batch(), uniform.training_batch()
    torch.testing.assert_close(x[:4], y[:4], rtol=0, atol=0)
    for block in x.split([5, 2, 2], dim=1):
        torch.testing.assert_close(block.norm(dim=1), torch.ones(8, dtype=torch.float64))
    for exact in x[-2:]:
        assert (replay.replay == exact).all(dim=1).any()
    for jittered in x[-4:-2]:
        assert not (replay.replay == jittered).all(dim=1).any()
    dup = torch.cat((replay.replay, replay.replay))
    chosen = diverse(dup, torch.arange(len(dup)), 100, .05)
    assert len(chosen) == len(replay.replay)


def test_interrupted_refresh_is_transactional_and_initial_stop_resumes(tmp_path, monkeypatch):
    cfg, base, audit, source, witness = fixture(tmp_path)
    runner = Runner(cfg, base, audit, tmp_path / 'run', witness, init_from=source)
    runner.stop_requested = True
    assert runner.run()['status'] == 'incomplete'
    assert not (runner.out / 'final.pt').exists()
    resumed = Runner(cfg, base, audit, runner.out, witness, resume=True)
    resumed.run(max_updates=1)
    before = copy.deepcopy(resumed.training_state())
    original = resumed.check_stop
    checks = [0]
    def stop_during_search():
        checks[0] += 1
        if checks[0] == 4:
            raise Stopped()
        original()
    monkeypatch.setattr(resumed, 'check_stop', stop_during_search)
    with pytest.raises(Stopped):
        resumed.refresh()
    assert_nested(before, resumed.training_state())


def test_tail_regression_and_numerical_failure_cannot_replace_best(tmp_path, monkeypatch):
    cfg, base, audit, source, witness = fixture(tmp_path)
    runner = Runner(cfg, base, audit, tmp_path / 'run', witness, init_from=source)
    runner.observe()
    before = copy.deepcopy(runner.best)
    def bad_tail(points):
        if points is runner.challenge:
            return {**runner.initial['challenge'], 'nij_sup': 2*runner.initial['challenge']['nij_sup']}
        return {**runner.initial, 'nij_rms': .5*runner.initial['nij_rms']}
    monkeypatch.setattr(runner, 'evaluate', bad_tail)
    runner.observe()
    assert_nested(before, runner.best)
    monkeypatch.setattr(runner, 'evaluate', lambda _: {**runner.initial, 'nij_rms': 0, 'tangency_max': 1})
    runner.observe()
    assert_nested(before, runner.best)


def test_witness_provenance_geometry_resume_and_config_rejections(tmp_path):
    cfg, base, audit, source, witness = fixture(tmp_path)
    with pytest.raises(ValueError, match='checksum'):
        Runner(cfg, base, audit, tmp_path / 'wrong', witness, init_from=source, expected_source_sha256='0'*64)
    with pytest.raises(ValueError, match='mismatch'):
        load_witnesses(witness, '0'*64)
    runner = Runner(cfg, base, audit, tmp_path / 'run', witness, init_from=source)
    runner.run(max_updates=1)
    changed = copy.deepcopy(cfg)
    changed['matched']['sampling']['jitter'] = .1
    with pytest.raises(ValueError, match='configuration'):
        Runner(changed, base, audit, runner.out, witness, resume=True)
    data = torch.load(witness / 'diagnostic_points.pt', weights_only=False)
    data['points'][0, 0] += 1
    torch.save(data, witness / 'diagnostic_points.pt')
    with pytest.raises(ValueError, match='mismatch'):
        load_witnesses(witness, sha256(source))
    report = json.loads((witness / 'report.json').read_text())
    report['points_sha256'] = sha256(witness / 'diagnostic_points.pt')
    (witness / 'report.json').write_text(json.dumps(report))
    with pytest.raises(ValueError, match='off the product'):
        load_witnesses(witness, sha256(source))
    for mode in ('uniform', 'replay'):
        production = read_config(ROOT / 'configs/ambient/s4xs1xs1_cs/matched.yaml', mode)
        assert production['steps'] == 2000 and production['batch'] == 32
    with pytest.raises(ValueError):
        Sampling(target_count=32).validate(32)
    with pytest.raises(ValueError):
        Sampling(jitter=float('nan')).validate(32)


def test_summary_missing_and_hash_matched_completed_pair(tmp_path):
    root = tmp_path / 'pair'
    root.mkdir()
    assert not comparison(root)['comparable']
    cfg, base, audit, source, witness = fixture(tmp_path)
    for mode in ('uniform', 'replay'):
        arm_cfg = {**cfg, 'matched': {**cfg['matched'], 'mode': mode}}
        runner = Runner(arm_cfg, base, audit, root / mode, witness, init_from=source)
        runner.run()
        directory = runner.out / 'diagnostic'
        directory.mkdir()
        (directory / 'report.json').write_text(json.dumps({
            'status': 'completed_diagnostic', 'sha256': sha256(runner.out / 'final.pt'),
            'config': {'seed': 123}, 'uniform': {'nij_rms': .2},
            'combined_completed_groups': {'nij_max_found': 1}}))
    result = comparison(root)
    assert result['comparable'] and not result['numerically_validated_cs']
    (root / 'replay/final.pt').write_bytes(b'changed')
    assert not comparison(root)['comparable']


def test_pbs_train_then_diagnose_and_never_diagnose_failed_training(shell_env, tmp_path):
    script = ROOT / 'scratch/submit_S4xS1xS1_CS_matched.sh'
    if not script.exists():
        pytest.skip('scratch intentionally ignored')
    python = Path(shell_env['PATH'].split(':')[0]) / 'python3'
    python.write_text('#!/bin/bash\nprintf "%s\\n" "$@"\nexit "${FAKE_EXIT:-0}"\n')
    out = tmp_path / 'arm output'
    out.mkdir()
    (out / 'final.pt').write_bytes(b'fixture')
    env = {**shell_env, 'BASE_DIR': '/absolute/base', 'BASE_AUDIT': '/absolute/audit',
           'CHECKPOINT': '/absolute/source.pt', 'SOURCE_SHA256': 'a'*64,
           'WITNESS_DIR': '/absolute/witness', 'OUT': str(out), 'MODE': 'replay'}
    result = subprocess.run(['bash', str(script)], env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert '--mode\nreplay' in result.stdout and 'extra/diagnose_cs.py' in result.stdout
    assert 'diagnostic_matched.yaml' in result.stdout
    stopped = subprocess.run(['bash', str(script)], env={**env, 'FAKE_EXIT': '75'}, capture_output=True, text=True)
    assert stopped.returncode == 75 and 'extra/diagnose_cs.py' not in stopped.stdout


@pytest.mark.parametrize('failure', [None, 'FAIL_CALL', 'EMPTY_CALL'])
def test_launcher_dependencies_provenance_and_partial_failure(shell_env, tmp_path, failure):
    script = ROOT / 'scratch/launch_S4xS1xS1_CS_matched.sh'
    if not script.exists():
        pytest.skip('scratch intentionally ignored')
    work = tmp_path / 'submission'
    work.mkdir()
    for name in ('launch_S4xS1xS1_CS_matched.sh', 'submit_S4xS1xS1_CS_matched.sh', 'submit_S4xS1xS1_CS_matched_summary.sh'):
        shutil.copy(ROOT / 'scratch' / name, work / name)
    acs = work / 'outputs/s4xs1xs1_acs/polish3_s307/polish_trilin_uniform/s307'
    cs = work / 'outputs/s4xs1xs1_cs/continue_stages2_s107/run'
    diagnostic = work / 'outputs/s4xs1xs1_cs/diagnostic_stages2_s107'
    for folder in (acs, cs, diagnostic):
        folder.mkdir(parents=True)
    for path in (acs / 'final.pt', acs / 'final.audit.full.seed900001.json', cs / 'final.pt',
                 diagnostic / 'report.json', diagnostic / 'diagnostic_points.pt'):
        path.write_bytes(b'fixture')
    binary = Path(shell_env['PATH'].split(':')[0]) / 'sha256sum'
    binary.write_text('#!/bin/bash\nif [[ "$1" = *diagnostic_points.pt ]]; then\necho 33a263a800a0413132479023029bc0135dc98c17dc0525c60f94b480142b7c87\nelse\necho 886cba9805a2c83e1a7acc35cb14b5e8a0f7d3daa04206706bd862e7bc0ffd11\nfi\n')
    binary.chmod(0o755)
    env = {**shell_env}
    if failure:
        env[failure] = '2'
    result = subprocess.run(['bash', str(work / script.name)], cwd=work, env=env, capture_output=True, text=True)
    submitted = calls(env)
    assert result.returncode == (1 if failure else 0), result.stderr
    assert len(submitted) == (2 if failure else 3)
    if not failure:
        assert submitted[-1][submitted[-1].index('-W')+1] == 'depend=afterany:101.server:102.server'
    else:
        assert all('-W' not in args for args in submitted)
    manifests = list((work / 'outputs/s4xs1xs1_cs').glob('matched_*/submission.tsv'))
    assert len(manifests) == 1 and 'uniform\t101.server' in manifests[0].read_text()
