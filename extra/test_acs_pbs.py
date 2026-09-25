"""Read-only summary and mocked PBS submission regressions; no real qsub calls."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

REPO = Path(__file__).resolve().parents[1]


def script(name):
    path = REPO / 'scratch' / name
    if not path.is_file():
        pytest.skip('PBS scripts are deliberately gitignored and copied separately')
    return str(path)


@pytest.fixture
def shell_env(tmp_path):
    work = tmp_path / 'submission directory'
    work.mkdir()
    binaries = tmp_path / 'bin'
    binaries.mkdir()
    for name in ('conda', 'module'):
        stub = binaries / name
        stub.write_text('#!/bin/bash\nexit 0\n')
        stub.chmod(0o755)
    python = binaries / 'python3'
    python.write_text(f'#!/bin/bash\nexec "{sys.executable}" "$@"\n')
    python.chmod(0o755)
    qsub = binaries / 'qsub'
    qsub.write_text(f'#!{sys.executable}\n' + '''import json, os, sys
from pathlib import Path
record = Path(os.environ['QSUB_RECORD'])
calls = [json.loads(line) for line in record.read_text().splitlines()] if record.exists() else []
calls.append(sys.argv[1:])
with record.open('a') as out:
    out.write(json.dumps(sys.argv[1:]) + '\\n')
n = len(calls)
if str(n) == os.environ.get('FAIL_CALL'):
    print('simulated resource limit', file=sys.stderr)
    sys.exit(1)
if str(n) != os.environ.get('EMPTY_CALL'):
    print(f'{100+n}.server')
''')
    qsub.chmod(0o755)
    return {**os.environ, 'PATH': str(binaries) + ':' + os.environ['PATH'],
            'REPO_DIR': str(REPO), 'PBS_O_WORKDIR': str(work), 'PBS_JOBID': '200.server',
            'CAMPAIGN': 'test_campaign', 'SEED': '7',
            'QSUB_RECORD': str(tmp_path / 'qsub.jsonl')}


def calls(env):
    path = Path(env['QSUB_RECORD'])
    return [json.loads(line) for line in path.read_text().splitlines()] if path.exists() else []


def test_launcher_builds_valid_afterany_dependency_and_manifest(shell_env):
    result = subprocess.run(['bash', script('submit_S4xS1xS1_ACS_pilots.sh')],
                            cwd=shell_env['PBS_O_WORKDIR'], env=shell_env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    submitted = calls(shell_env)
    assert len(submitted) == 5
    assert submitted[-1][submitted[-1].index('-W') + 1] == 'depend=afterany:101.server:102.server:103.server:104.server'
    root = Path(shell_env['PBS_O_WORKDIR']) / 'outputs/s4xs1xs1_acs/test_campaign'
    manifest = (root / 'submission.tsv').read_text()
    assert len(manifest.splitlines()) == 6
    assert 'summary\t105.server\t' in manifest
    for args in submitted[:4]:
        assert f'REPO_DIR={REPO}' in args[args.index('-v') + 1]
        assert args[-1].endswith('submit_S4xS1xS1_ACS_v2.sh')
    assert f'SUMMARY_ROOT={root}' in submitted[-1][submitted[-1].index('-v') + 1]


@pytest.mark.parametrize('failure', ['FAIL_CALL', 'EMPTY_CALL'])
def test_launcher_stops_on_failed_or_empty_qsub_without_empty_dependency(shell_env, failure):
    env = {**shell_env, failure: '2'}
    result = subprocess.run(['bash', script('submit_S4xS1xS1_ACS_pilots.sh')],
                            cwd=env['PBS_O_WORKDIR'], env=env, capture_output=True, text=True)
    assert result.returncode == 1
    assert len(calls(env)) == 2
    assert all('-W' not in args for args in calls(env))


def test_launcher_missing_code_fails_before_submitting(shell_env, tmp_path):
    env = {**shell_env, 'REPO_DIR': str(tmp_path / 'old_checkout')}
    result = subprocess.run(['bash', script('submit_S4xS1xS1_ACS_pilots.sh')],
                            cwd=env['PBS_O_WORKDIR'], env=env, capture_output=True, text=True)
    assert result.returncode == 2
    assert 'Commit/add the new files' in result.stderr
    assert not calls(env)


def test_training_explicit_checkout_works_from_unrelated_submission_directory(shell_env):
    # A tiny wall budget tests actual config dispatch without running training.
    env = {**shell_env, 'RECIPE': 'smoke', 'WALL_MINUTES': '0.00000001'}
    result = subprocess.run(['bash', script('submit_S4xS1xS1_ACS_v2.sh')],
                            cwd=env['PBS_O_WORKDIR'], env=env, capture_output=True, text=True)
    assert result.returncode == 75, result.stderr
    path = Path(env['PBS_O_WORKDIR']) / 'outputs/s4xs1xs1_acs/pilots/smoke/s7/200.server/result.json'
    assert json.loads(path.read_text())['status'] == 'incomplete'


def test_training_missing_checkout_explains_paths_and_uncommitted_files(shell_env):
    env = shell_env.copy()
    env.pop('REPO_DIR')
    result = subprocess.run(['bash', script('submit_S4xS1xS1_ACS_v2.sh')],
                            cwd=env['PBS_O_WORKDIR'], env=env, capture_output=True, text=True)
    assert result.returncode == 2
    assert f'PBS_O_WORKDIR={env["PBS_O_WORKDIR"]}' in result.stderr
    assert 'Commit AND push' in result.stderr


def test_summary_job_reports_missing_and_incomplete_pilots(shell_env, tmp_path):
    root = tmp_path / 'results'
    run = root / 'pilot_all2_uniform' / 's7'
    run.mkdir(parents=True)
    record = {'schema_version': 2, 'case': 's4xs1xs1', 'task': 'acs', 'seed': 7,
              'formulation': 'ambient', 'feat': 'all2', 'status': 'incomplete',
              'steps_done': 0, 'counts': {'poly': 0, 'refine': 0, 'lbfgs': 0},
              'stopped_by_wall': True, 'elapsed_min': 0.1}
    (run / 'result.json').write_text(json.dumps(record))
    env = {**shell_env, 'SUMMARY_ROOT': str(root), 'PILOT_SEED': '7'}
    result = subprocess.run(['bash', script('submit_S4xS1xS1_ACS_summary.sh')],
                            env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    report = (root / 'summary.200.server.txt').read_text()
    assert report.count('MISSING RESULT:') == 3
    assert 'status=incomplete' in report
    assert 'wall-stopped=True' in report
    assert not calls(env)


def test_summary_expected_run_warning_when_no_results(tmp_path):
    result = subprocess.run([sys.executable, str(REPO / 'extra/campaign_summary.py'),
                             str(tmp_path), '--all', '--expected-run', 'pilot/s7'],
                            capture_output=True, text=True)
    assert result.returncode == 0
    assert 'MISSING RESULT: pilot/s7' in result.stdout
    assert 'no result.json found' in result.stdout


def test_followup_launcher_three_jobs_checksum_and_summary(shell_env, tmp_path):
    import hashlib
    source = tmp_path / 'audited candidate.pt'
    source.write_bytes(b'checkpoint fixture\n')
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    env = {**shell_env, 'STAGE': 'followup', 'SOURCE_CHECKPOINT': str(source), 'SOURCE_SHA256': digest}
    result = subprocess.run(['bash', script('submit_S4xS1xS1_ACS_pilots.sh')],
                            cwd=env['PBS_O_WORKDIR'], env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    submitted = calls(env)
    assert len(submitted) == 4  # three training jobs plus the held summary
    continuation = submitted[0][submitted[0].index('-v') + 1]
    assert 'RECIPE=polish_trilin_uniform,SEED=107' in continuation
    assert f'INIT_FROM={source},SOURCE_SHA256={digest}' in continuation
    for index, seed in ((1, 8), (2, 9)):
        fresh = submitted[index][submitted[index].index('-v') + 1]
        assert f'RECIPE=pilot_trilin_uniform,SEED={seed}' in fresh
        assert 'INIT_FROM=' not in fresh
    assert submitted[-1][submitted[-1].index('-W') + 1] == 'depend=afterany:101.server:102.server:103.server'
    assert source.read_bytes() == b'checkpoint fixture\n'


def test_followup_hash_mismatch_submits_nothing(shell_env, tmp_path):
    source = tmp_path / 'wrong.pt'
    source.write_bytes(b'changed')
    env = {**shell_env, 'STAGE': 'followup', 'SOURCE_CHECKPOINT': str(source), 'SOURCE_SHA256': '0' * 64}
    result = subprocess.run(['bash', script('submit_S4xS1xS1_ACS_pilots.sh')],
                            cwd=env['PBS_O_WORKDIR'], env=env, capture_output=True, text=True)
    assert result.returncode == 2
    assert 'checksum does not match' in result.stderr
    assert calls(env) == []


def test_polish_wrapper_requires_source_or_resume(shell_env):
    env = {**shell_env, 'RECIPE': 'polish_trilin_uniform', 'SEED': '107'}
    result = subprocess.run(['bash', script('submit_S4xS1xS1_ACS_v2.sh')],
                            cwd=env['PBS_O_WORKDIR'], env=env, capture_output=True, text=True)
    assert result.returncode == 2
    assert 'refusing a fresh MLP start' in result.stderr


def test_summary_reads_followup_manifest_without_assuming_four_pilots(shell_env, tmp_path):
    root = tmp_path / 'followup'
    root.mkdir()
    manifest = 'role\tjob_id\toutput\n'
    for recipe, seed, job in [('polish_trilin_uniform', 107, 101), ('pilot_trilin_uniform', 8, 102),
                              ('pilot_trilin_uniform', 9, 103)]:
        manifest += f'{recipe}\t{job}.server\t{root / recipe / f"s{seed}"}\n'
    manifest += f'summary\t200.server\t{root / "summary.200.server.txt"}\n'
    (root / 'submission.tsv').write_text(manifest)
    env = {**shell_env, 'SUMMARY_ROOT': str(root)}
    result = subprocess.run(['bash', script('submit_S4xS1xS1_ACS_summary.sh')],
                            env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    report = (root / 'summary.200.server.txt').read_text()
    assert report.count('MISSING RESULT:') == 3
    assert 'polish_trilin_uniform/s107' in report
    assert 'pilot_trilin_uniform/s8' in report
    assert 'pilot_trilin_uniform/s9' in report


def test_polish_wrapper_passes_checked_continuation_to_real_runner(shell_env, tmp_path):
    import hashlib
    import torch
    from extra.acs_training import Runner, read_config
    config = read_config(REPO / 'configs/ambient/s4xs1xs1_acs/polish_trilin_uniform.yaml')
    source = Runner(config, tmp_path / 'fixture')
    source.model.refinement_phase(zero_output=False)
    checkpoint = tmp_path / 'candidate.pt'
    torch.save(source.payload(), checkpoint)
    digest = hashlib.sha256(checkpoint.read_bytes()).hexdigest()
    env = {**shell_env, 'RECIPE': 'polish_trilin_uniform', 'SEED': '107',
           'INIT_FROM': str(checkpoint), 'SOURCE_SHA256': digest, 'WALL_MINUTES': '0.00000001'}
    result = subprocess.run(['bash', script('submit_S4xS1xS1_ACS_v2.sh')],
                            cwd=env['PBS_O_WORKDIR'], env=env, capture_output=True, text=True)
    assert result.returncode == 75, result.stderr
    result_path = Path(env['PBS_O_WORKDIR']) / 'outputs/s4xs1xs1_acs/pilots/polish_trilin_uniform/s107/200.server/result.json'
    record = json.loads(result_path.read_text())
    assert record['phase'] == 'refine'
    assert record['source']['sha256'] == digest
    assert hashlib.sha256(checkpoint.read_bytes()).hexdigest() == digest
