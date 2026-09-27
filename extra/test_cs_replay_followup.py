"""Follow-up configurations, source-chain preflight, summaries and PBS wiring."""
import json
from pathlib import Path
import shutil
import subprocess

import pytest

from extra import cs_followup_inputs as inputs
from extra.cs_matched import comparison, read_config
from extra.test_acs_pbs import calls, shell_env

ROOT = Path(__file__).resolve().parents[1]


def input_fixture(root, monkeypatch):
    base = root / 'outputs/s4xs1xs1_acs/polish3_s307/polish_trilin_uniform/s307'
    source = root / 'outputs/s4xs1xs1_cs/matched_s207_hs8qhz51/replay'
    base.mkdir(parents=True)
    (source / 'diagnostic').mkdir(parents=True)
    (base / 'final.pt').write_bytes(b'ACS fixture')
    (source / 'final.pt').write_bytes(b'CS fixture')
    (source / 'diagnostic/diagnostic_points.pt').write_bytes(b'witness fixture')
    base_hash = inputs.digest(base / 'final.pt')
    monkeypatch.setattr(inputs, 'ACS_SHA256', base_hash)
    audit = {'case': 's4xs1xs1', 'audit_tier': 'full', 'sha256': base_hash,
             'numerically_validated': True, 'screening_targets_passed': True}
    audit_path = base / 'final.audit.full.seed900001.json'
    audit_path.write_text(json.dumps(audit))
    evidence = {'sha256': base_hash, 'audit_sha256': inputs.digest(audit_path)}
    result = {'case': 's4xs1xs1', 'task': 'cs', 'mode': 'replay', 'seed': 207,
              'steps_done': 2000, 'steps_planned': 2000, 'selected_step': 2000,
              'status': 'completed_pending_independent_audit', 'audited_base': evidence,
              'source': {'sha256': inputs.PARENT_SHA256}}
    diagnostic = {'status': 'completed_diagnostic', 'audited_base': evidence,
                  'sha256': inputs.digest(source / 'final.pt'),
                  'points_sha256': inputs.digest(source / 'diagnostic/diagnostic_points.pt')}
    (source / 'result.json').write_text(json.dumps(result))
    (source / 'diagnostic/report.json').write_text(json.dumps(diagnostic))
    return source, base


def test_followup_configs_differ_only_in_target_count():
    a = read_config(ROOT / 'configs/ambient/s4xs1xs1_cs/replay_followup25.yaml', 'replay')
    b = read_config(ROOT / 'configs/ambient/s4xs1xs1_cs/replay_followup12p5.yaml', 'replay')
    assert a['seed'] == 307 and a['steps'] == 5000 and a['lr'] == 2e-5
    assert a['batch'] == 32
    assert a['matched']['sampling'].pop('target_count') == 8
    assert b['matched']['sampling'].pop('target_count') == 4
    assert a == b


@pytest.mark.parametrize('corruption', ['source', 'points', 'audit', 'parent', 'step', 'status'])
def test_preflight_requires_matching_completed_reviewed_source(tmp_path, monkeypatch, corruption):
    source, base = input_fixture(tmp_path, monkeypatch)
    assert inputs.verify(source, base) == inputs.digest(source / 'final.pt')
    if corruption in ('source', 'points', 'audit'):
        path = {'source': source / 'final.pt', 'points': source / 'diagnostic/diagnostic_points.pt',
                'audit': base / 'final.audit.full.seed900001.json'}[corruption]
        path.write_bytes(path.read_bytes() + b' ')
    else:
        path = source / 'result.json'
        report = json.loads(path.read_text())
        if corruption == 'parent':
            report['source']['sha256'] = '0'*64
        elif corruption == 'step':
            report['selected_step'] = 1750
        else:
            report['status'] = 'incomplete'
        path.write_text(json.dumps(report))
    with pytest.raises(ValueError):
        inputs.verify(source, base)


def test_mixture_summary_accepts_only_intended_config_difference(tmp_path):
    assert not comparison(tmp_path, 'replay-mixtures')['comparable']
    for arm, config in [('replay25', 'replay_followup25.yaml'), ('replay12p5', 'replay_followup12p5.yaml')]:
        folder = tmp_path / arm
        (folder / 'diagnostic').mkdir(parents=True)
        (folder / 'final.pt').write_bytes(arm.encode())
        result = {'config': read_config(ROOT / 'configs/ambient/s4xs1xs1_cs' / config, 'replay'),
                  'status': 'completed_pending_independent_audit', 'source': {'sha256': 'a'*64},
                  'audited_base': {'sha256': 'b'*64}, 'initial_validation': {'nij_rms': .44},
                  'witness_input': {'points_sha256': 'c'*64}}
        (folder / 'result.json').write_text(json.dumps(result))
        diagnostic = {'status': 'completed_diagnostic', 'sha256': inputs.digest(folder / 'final.pt'),
                      'config': {'seed': 941307}}
        (folder / 'diagnostic/report.json').write_text(json.dumps(diagnostic))
    report = comparison(tmp_path, 'replay-mixtures')
    assert report['comparable'] and not report['numerically_validated_cs']
    # Normalization must not replace the actual sampling metadata in the report.
    assert report['arms']['replay25']['training']['config']['matched']['sampling']['target_count'] == 8
    result['config']['lr'] = 1e-4
    (folder / 'result.json').write_text(json.dumps(result))
    assert not comparison(tmp_path, 'replay-mixtures')['comparable']


def test_wrapper_passes_new_configs_and_summary_experiment(shell_env, tmp_path):
    train = ROOT / 'scratch/submit_S4xS1xS1_CS_matched.sh'
    summary = ROOT / 'scratch/submit_S4xS1xS1_CS_matched_summary.sh'
    if not train.exists() or not summary.exists():
        pytest.skip('scratch deliberately ignored')
    python = Path(shell_env['PATH'].split(':')[0]) / 'python3'
    python.write_text('#!/bin/bash\nprintf "%s\\n" "$@"\nexit "${FAKE_EXIT:-0}"\n')
    out = tmp_path / 'run'
    out.mkdir()
    (out / 'final.pt').write_bytes(b'fixture')
    env = {**shell_env, 'BASE_DIR': '/base', 'BASE_AUDIT': '/audit', 'CHECKPOINT': '/source',
           'SOURCE_SHA256': 'a'*64, 'WITNESS_DIR': '/witness', 'OUT': str(out), 'MODE': 'replay',
           'TRAIN_CONFIG': 'configs/ambient/s4xs1xs1_cs/replay_followup12p5.yaml',
           'DIAGNOSTIC_CONFIG': 'configs/ambient/s4xs1xs1_cs/diagnostic_replay_followup.yaml'}
    done = subprocess.run(['bash', str(train)], env=env, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert '--config\n' + env['TRAIN_CONFIG'] in done.stdout
    assert '--config\n' + env['DIAGNOSTIC_CONFIG'] in done.stdout
    env.update(COMPARISON_ROOT=str(out), EXPERIMENT='replay-mixtures')
    done = subprocess.run(['bash', str(summary)], env=env, capture_output=True, text=True)
    assert done.returncode == 0 and '--experiment\nreplay-mixtures' in done.stdout


@pytest.mark.parametrize('failure', [None, 'FAIL_CALL', 'EMPTY_CALL', 'preflight', 'old_wrapper'])
def test_followup_launcher_correct_paths_configs_dependencies_and_failures(shell_env, tmp_path, monkeypatch, failure):
    launcher = ROOT / 'scratch/launch_S4xS1xS1_CS_replay_followup.sh'
    if not launcher.exists():
        pytest.skip('scratch deliberately ignored')
    work = tmp_path / 'submission'
    work.mkdir()
    source, base = input_fixture(work, monkeypatch)
    for name in (launcher.name, 'submit_S4xS1xS1_CS_matched.sh', 'submit_S4xS1xS1_CS_matched_summary.sh'):
        shutil.copy(ROOT / 'scratch' / name, work / name)
    # Mock only the preflight's reviewed ACS constant, not its verification logic.
    python = Path(shell_env['PATH'].split(':')[0]) / 'python3'
    import sys
    python.write_text(f'#!{sys.executable}\nimport sys\nsys.path.insert(0, {str(ROOT)!r})\n'
                      f'from extra import cs_followup_inputs as m\nm.ACS_SHA256 = {inputs.ACS_SHA256!r}\n'
                      'sys.argv.pop(1)\nm.main()\n')
    env = {**shell_env}
    if failure in ('FAIL_CALL', 'EMPTY_CALL'):
        env[failure] = '2'
    elif failure == 'preflight':
        (source / 'final.pt').write_bytes(b'changed')
    elif failure == 'old_wrapper':
        (work / 'submit_S4xS1xS1_CS_matched.sh').write_text('#!/bin/bash\nexit 0\n')
    done = subprocess.run(['bash', str(work / launcher.name)], cwd=work, env=env, capture_output=True, text=True)
    submitted = calls(env)
    if failure in ('preflight', 'old_wrapper'):
        assert done.returncode == 2 and not submitted
        return
    assert done.returncode == (1 if failure else 0), done.stderr
    assert len(submitted) == (2 if failure else 3)
    for index, suffix in enumerate(('25', '12p5')):
        variables = submitted[index][submitted[index].index('-v')+1]
        assert f'CHECKPOINT={source}/final.pt' in variables
        assert f'WITNESS_DIR={source}/diagnostic' in variables
        assert f'TRAIN_CONFIG=configs/ambient/s4xs1xs1_cs/replay_followup{suffix}.yaml' in variables
        assert 'MODE=replay' in variables
    if not failure:
        assert submitted[2][submitted[2].index('-W')+1] == 'depend=afterany:101.server:102.server'
        assert 'EXPERIMENT=replay-mixtures' in submitted[2][submitted[2].index('-v')+1]
    manifests = list((work / 'outputs/s4xs1xs1_cs').glob('replay_followup_*/submission.tsv'))
    assert len(manifests) == 1 and 'replay25\t101.server' in manifests[0].read_text()
