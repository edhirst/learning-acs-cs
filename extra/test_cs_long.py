"""Small real CS runs plus mocked PBS; no remote training or submissions."""
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

from extra.acs_training import Stopped, sha256
from extra import cs_long, cs_followup_inputs
from extra.cs_long import Controls, Runner, read_config, summary
from extra.test_cs_matched import fixture, assert_nested
from extra.test_cs_replay_followup import input_fixture
from extra.test_acs_training import one_thread
from extra.test_acs_pbs import shell_env, calls

ROOT = Path(__file__).resolve().parents[1]


def setup(tmp_path):
    cfg, base, audit, source, witnesses = fixture(tmp_path)
    cfg['lr'] = 1e-5
    cfg['long_run'] = asdict(Controls(monitor_every=1, monitor_pool=4, monitor_starts=2,
                                    monitor_steps=2, plateau_window=10))
    return cfg, base, audit, source, witnesses


def test_long_real_training_dual_candidates_constant_lr_and_immutable_inputs(tmp_path):
    cfg, base, audit, source, witness = setup(tmp_path)
    paths = [source, base / 'final.pt', audit, witness / 'diagnostic_points.pt']
    hashes = [sha256(p) for p in paths]
    runner = Runner(cfg, base, audit, tmp_path / 'run', witness, init_from=source)
    report = runner.run()
    assert report['status'] == 'completed_pending_independent_audit'
    assert report['ready_for_diagnostics'] and not report['numerically_validated_cs']
    assert runner.optimizer.param_groups[0]['lr'] == 1e-5
    assert len(runner.monitor_points) == 6
    assert runner.best['metrics']['acs_floorfree_sup'] < 1e-12
    assert runner.tail_best['metrics']['tail_max'] <= runner.best['metrics']['tail_max']
    assert (runner.out / 'final_average.pt').exists() and (runner.out / 'final_tail.pt').exists()
    assert [sha256(p) for p in paths] == hashes
    assert not summary(runner.out)['diagnostics_complete']


@pytest.mark.parametrize('interrupt_monitor', [False, True])
def test_exact_resume_including_pending_monitor_observation(tmp_path, monkeypatch, interrupt_monitor):
    cfg, base, audit, source, witness = setup(tmp_path)
    full = Runner(cfg, base, audit, tmp_path / 'full', witness, init_from=source)
    full.run()
    split = Runner(cfg, base, audit, tmp_path / 'split', witness, init_from=source)
    if interrupt_monitor:
        original = split.refresh_monitor
        def interrupt():
            # Raise within ascent, after some values were already computed.
            checks = [0]
            old = split.check_stop
            def check():
                checks[0] += 1
                if checks[0] == 4:
                    raise Stopped()
                old()
            with monkeypatch.context() as patch:
                patch.setattr(split, 'check_stop', check)
                original()
        monkeypatch.setattr(split, 'refresh_monitor', interrupt)
        result = split.run()
        assert result['pending_observation'] and result['steps_done'] == 1
        assert not len(split.monitor_points)
    else:
        result = split.run(max_updates=1)
    assert not result['ready_for_diagnostics']
    assert not (split.out / 'final_average.pt').exists()
    resumed = Runner(cfg, base, audit, split.out, witness, resume=True)
    resumed.run()
    a, b = [torch.load(p / 'state.pt', weights_only=False) for p in (full.out, split.out)]
    for key in ('model', 'optimizer', 'scheduler', 'train_rng', 'validation', 'best', 'initial', 'training_extension'):
        assert_nested(a[key], b[key])


def test_dual_selection_and_growth_plateau_guards(tmp_path, monkeypatch):
    cfg, base, audit, source, witness = setup(tmp_path)
    runner = Runner(cfg, base, audit, tmp_path / 'run', witness, init_from=source)
    runner.observe()
    initial = copy.deepcopy(runner.initial)
    metrics = copy.deepcopy(initial)
    metrics['nij_rms'] *= .9
    metrics['tail_max'] *= 1.01
    monkeypatch.setattr(runner, 'measure', lambda: copy.deepcopy(metrics))
    runner.observe()
    assert runner.best['metrics']['nij_rms'] < initial['nij_rms']
    assert runner.tail_best['metrics']['tail_max'] == initial['tail_max']
    metrics['nij_rms'] = initial['nij_rms'] * 1.01
    metrics['tail_max'] = initial['tail_max'] * .8
    runner.observe()
    assert runner.tail_best['metrics']['tail_max'] < runner.best['metrics']['tail_max']
    metrics['J_frobenius_max'] = initial['J_frobenius_max'] * 2.1
    runner.observe()
    assert runner.terminal_reason == 'growth_guard'
    report = runner.run()
    assert report['status'] == 'stopped_for_review' and report['ready_for_diagnostics']
    assert not report['numerically_validated_cs']

    cfg['long_run']['plateau_window'] = 2
    other = Runner(cfg, base, audit, tmp_path / 'plateau', witness, init_from=source)
    other.observe()
    flat = copy.deepcopy(other.initial)
    monkeypatch.setattr(other, 'measure', lambda: copy.deepcopy(flat))
    monkeypatch.setattr(other, 'refresh_monitor', lambda: setattr(other, 'last_monitor', other.step))
    for step in (1, 2):
        other.step = step
        other.observe()
    assert other.terminal_reason == 'plateau_review'


def test_cli_receipt_chain_noop_and_failed_predecessor(tmp_path, monkeypatch):
    cfg, base, audit, source, witness = setup(tmp_path)
    config = tmp_path / 'config.yaml'
    training = {k: v for k, v in cfg.items() if k not in ('matched', 'long_run')}
    config.write_text(yaml.safe_dump({'training': training, 'sampling': cfg['matched']['sampling'], 'long_run': cfg['long_run']}))
    out = tmp_path / 'run'
    argv = ['cs_long', 'train', '--config', str(config), '--base-from', str(base), '--base-audit', str(audit),
            '--out', str(out), '--witness-dir', str(witness), '--init-from', str(source),
            '--expected-source-sha256', sha256(source), '--max-updates', '1']
    monkeypatch.setattr(sys, 'argv', argv + ['--allocation', '2'])
    with pytest.raises(SystemExit):
        cs_long.main()
    for allocation in (1, 2, 3):
        monkeypatch.setattr(sys, 'argv', argv + ['--allocation', str(allocation)])
        cs_long.main()
        assert (out / f'allocation_{allocation}.json').exists()
    before = sha256(out / 'final_average.pt')
    monkeypatch.setattr(sys, 'argv', argv + ['--allocation', '4'])
    cs_long.main()
    assert sha256(out / 'final_average.pt') == before
    assert (out / 'allocation_4.json').exists()
    # Repeating an older slot cannot silently resume a newer state.
    monkeypatch.setattr(sys, 'argv', argv + ['--allocation', '2'])
    with pytest.raises(SystemExit):
        cs_long.main()


def test_configs_and_s307_source_verification(tmp_path, monkeypatch):
    cfg = read_config(ROOT / 'configs/ambient/s4xs1xs1_cs/long_replay.yaml')
    assert cfg['steps'] == 50000 and cfg['lr'] == 1e-5
    assert cfg['matched']['sampling']['target_count'] == 4
    with pytest.raises(ValueError):
        Controls(monitor_every=251).validate(500)
    source, base = input_fixture(tmp_path, monkeypatch)
    r = json.loads((source / 'result.json').read_text())
    r.update(seed=307, steps_done=5000, steps_planned=5000, selected_step=4750,
             config={'matched': {'sampling': {'target_count': 4}}})
    r['source']['sha256'] = '18abf4ce5e70a3bffe3d9264016cda7f4ff3141de84f38124f2e4e2c01bdb4d4'
    (source / 'result.json').write_text(json.dumps(r))
    assert cs_followup_inputs.verify(source, base, 's307') == sha256(source / 'final.pt')
    r['config']['matched']['sampling']['target_count'] = 8
    (source / 'result.json').write_text(json.dumps(r))
    with pytest.raises(ValueError, match='12.5'):
        cs_followup_inputs.verify(source, base, 's307')


def test_summary_checks_final_hashes_even_after_review_stop(tmp_path):
    report = {'status': 'stopped_for_review', 'terminal_reason': 'plateau_review',
              'steps_done': 10000, 'steps_planned': 50000, 'ready_for_diagnostics': True}
    (tmp_path / 'result.json').write_text(json.dumps(report))
    for name in ('average', 'tail'):
        ck = tmp_path / f'final_{name}.pt'
        ck.write_bytes(name.encode())
        folder = tmp_path / f'diagnostic_{name}'
        folder.mkdir()
        (folder / 'report.json').write_text(json.dumps({'status': 'completed_diagnostic', 'sha256': sha256(ck)}))
    result = summary(tmp_path)
    assert result['diagnostics_complete'] and not result['numerically_validated_cs']
    (tmp_path / 'final_tail.pt').write_bytes(b'changed')
    assert not summary(tmp_path)['diagnostics_complete']


@pytest.mark.parametrize('failure', [None, 'FAIL_CALL', 'EMPTY_CALL'])
def test_long_launcher_chain_and_submission_failure(shell_env, tmp_path, failure):
    script = ROOT / 'scratch/launch_S4xS1xS1_CS_long.sh'
    if not script.exists():
        pytest.skip('scratch deliberately ignored')
    work = tmp_path / 'submission'
    (work / 'outputs/s4xs1xs1_cs').mkdir(parents=True)
    for name in (script.name, 'submit_S4xS1xS1_CS_long.sh'):
        shutil.copy(ROOT / 'scratch' / name, work / name)
    python = Path(shell_env['PATH'].split(':')[0]) / 'python3'
    python.write_text('#!/bin/bash\n[[ "$*" = *"--stage s307"* ]] || exit 2\nprintf "%064d\\n" 0\n')
    env = dict(shell_env)
    if failure:
        env[failure] = '3'
    result = subprocess.run(['bash', str(work / script.name)], cwd=work, env=env, capture_output=True, text=True)
    submitted = calls(env)
    assert result.returncode == (1 if failure else 0), result.stderr
    assert len(submitted) == (3 if failure else 7)
    for index, args in enumerate(submitted):
        if index:
            assert args[args.index('-W')+1] == f'depend=afterany:{100+index}.server'
        variables = args[args.index('-v')+1]
        assert 'replay_followup_s307_YIFBe02B/replay12p5/final.pt' in variables
        if index < 5:
            assert f'ALLOCATION={index+1}' in variables
    if not failure:
        assert 'PHASE=diagnose' in submitted[5][submitted[5].index('-v')+1]
        assert 'PHASE=summary' in submitted[6][submitted[6].index('-v')+1]


def test_pbs_training_exit75_and_both_diagnostic_candidates(shell_env, tmp_path):
    script = ROOT / 'scratch/submit_S4xS1xS1_CS_long.sh'
    if not script.exists():
        pytest.skip('scratch deliberately ignored')
    python = Path(shell_env['PATH'].split(':')[0]) / 'python3'
    python.write_text('#!/bin/bash\nprintf "%s\\n" "$@"\nexit "${FAKE_EXIT:-0}"\n')
    out = tmp_path / 'run'
    out.mkdir()
    for name in ('final_average.pt', 'final_tail.pt', 'allocation_5.json'):
        (out / name).write_bytes(b'fixture')
    env = {**shell_env, 'BASE_DIR': '/base', 'BASE_AUDIT': '/audit', 'CHECKPOINT': '/source',
           'SOURCE_SHA256': 'a'*64, 'WITNESS_DIR': '/witness', 'OUT': str(out), 'PHASE': 'train', 'ALLOCATION': '2'}
    result = subprocess.run(['bash', str(script)], env={**env, 'FAKE_EXIT': '75'}, capture_output=True, text=True)
    assert result.returncode == 75
    assert '--allocation\n2' in result.stdout and '--max-updates\n10000' in result.stdout
    result = subprocess.run(['bash', str(script)], env={**env, 'PHASE': 'diagnose', 'ALLOCATIONS': '5'}, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.count('extra/diagnose_cs.py') == 2
    assert 'final_average.pt' in result.stdout and 'final_tail.pt' in result.stdout
