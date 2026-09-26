"""Checked CS continuation: test-only ACS fixture, real optimization and AD."""
import copy
import json
from pathlib import Path
import subprocess

import pytest
import torch

from extra.cs_training import Runner
from extra.cs_pilot import checked_base, read_config
from extra.acs_training import sha256
from extra.campaign_unit import Net
from extra.test_acs_recipe_cs_pilot import fixture_base
from extra.test_acs_training import one_thread
from extra.test_acs_pbs import shell_env

ROOT = Path(__file__).resolve().parents[1]


def setup(tmp_path):
    base, audit = fixture_base(tmp_path)
    cfg = read_config(ROOT / 'configs/ambient/s4xs1xs1_cs/continue.yaml')
    cfg.update(steps=3, hidden=[8], batch=2, n_validation=4, n_dense=4,
               eval_chunk=2, eval_every=1, ckpt_every=1, adv_starts=2, adv_steps=2)
    net = Net(9, 9, 2, [8], .01).double()
    source = tmp_path / 'source.pt'
    torch.save({'case': 's4xs1xs1', 'task': 'cs', 'kind': 'ambient', 'model_type': 'mlp',
                'hidden': [8], 'stages': 2, 'model': net.state_dict(),
                'audited_base': checked_base(base, audit)}, source)
    return cfg, base, audit, source


def test_continuation_preserves_source_then_trains_and_selects(tmp_path):
    cfg, base, audit, source = setup(tmp_path)
    original = torch.load(source, weights_only=False)
    digest = sha256(source)
    runner = Runner(cfg, base, audit, tmp_path / 'run', init_from=source,
                    expected_source_sha256=digest)
    for key in original['model']:
        torch.testing.assert_close(runner.model.state_dict()[key], original['model'][key], atol=0, rtol=0)
    result = runner.run()
    assert result['status'] == 'completed_pending_independent_audit'
    assert result['dense']['acs_floorfree_sup'] < 1e-12
    assert result['dense']['conjugation_condition_max'] >= 1
    assert result['best_validation']['nij_rms'] <= result['initial_validation']['nij_rms']
    assert not result['numerically_validated_cs']
    assert sha256(source) == digest
    assert (runner.out / 'final.pt').is_file()
    with pytest.raises(ValueError, match='Already completed'):
        Runner(cfg, base, audit, runner.out, resume=True)


def test_exact_resume_model_optimizer_sampling_and_best(tmp_path):
    cfg, base, audit, source = setup(tmp_path)
    full = Runner(cfg, base, audit, tmp_path / 'full', init_from=source)
    full.run()
    split = Runner(cfg, base, audit, tmp_path / 'split', init_from=source)
    assert split.run(max_updates=1)['status'] == 'incomplete'
    assert not (split.out / 'final.pt').exists()
    resumed = Runner(cfg, base, audit, split.out, resume=True)
    resumed.run()
    for key in full.model.state_dict():
        torch.testing.assert_close(full.model.state_dict()[key], resumed.model.state_dict()[key], atol=0, rtol=0)
    assert full.best['metrics'] == resumed.best['metrics']
    assert torch.equal(full.train_rng.get_state(), resumed.train_rng.get_state())
    a = torch.load(full.out / 'state.pt', weights_only=False)
    b = torch.load(split.out / 'state.pt', weights_only=False)
    assert a['step'] == b['step'] == 3
    assert a['scheduler'] == b['scheduler']
    for key, values in a['optimizer']['state'].items():
        for name, value in values.items():
            torch.testing.assert_close(value, b['optimizer']['state'][key][name], atol=0, rtol=0)


def test_reject_hash_base_architecture_and_changed_resume_config(tmp_path):
    cfg, base, audit, source = setup(tmp_path)
    with pytest.raises(ValueError, match='checksum'):
        Runner(cfg, base, audit, tmp_path / 'wrong', init_from=source, expected_source_sha256='0'*64)
    ck = torch.load(source, weights_only=False)
    for key, value in [('stages', 1), ('audited_base', {'sha256': '0'*64})]:
        bad = tmp_path / f'{key}.pt'
        torch.save({**ck, key: value}, bad)
        with pytest.raises(ValueError):
            Runner(cfg, base, audit, tmp_path / 'bad', init_from=bad)
    runner = Runner(cfg, base, audit, tmp_path / 'run', init_from=source)
    runner.run(max_updates=1)
    with pytest.raises(ValueError, match='configuration'):
        Runner({**cfg, 'steps': 4}, base, audit, runner.out, resume=True)


def test_wall_stop_and_bad_candidate_cannot_replace_best(tmp_path, monkeypatch):
    cfg, base, audit, source = setup(tmp_path)
    runner = Runner({**cfg, 'wall_minutes': 1e-12}, base, audit, tmp_path / 'run', init_from=source)
    result = runner.run()
    assert result['steps_done'] == 0 and result['status'] == 'incomplete'
    assert not (runner.out / 'final.pt').exists()
    before = copy.deepcopy(runner.best)
    with monkeypatch.context() as patch:
        patch.setattr(runner, 'evaluate', lambda _: {**before['metrics'], 'nij_rms': 0, 'conjugation_condition_max': 1e10})
        runner.observe()
    assert runner.best['metrics'] == before['metrics']


def test_pbs_continuation_passes_hash_resume_and_exit_code(shell_env):
    script = ROOT / 'scratch/submit_S4xS1xS1_CS_continue.sh'
    if not script.exists():
        pytest.skip('scratch scripts intentionally gitignored')
    python = Path(shell_env['PATH'].split(':')[0]) / 'python3'
    python.write_text('#!/bin/bash\nprintf "%s\\n" "$@"\nexit "${FAKE_EXIT:-0}"\n')
    env = {**shell_env, 'BASE_DIR': '/absolute/base', 'BASE_AUDIT': '/absolute/audit.json',
           'OUT': '/absolute/new run', 'INIT_FROM': '/absolute/source.pt', 'SOURCE_SHA256': 'a'*64}
    result = subprocess.run(['bash', str(script)], env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert '--expected-source-sha256\n' + 'a'*64 in result.stdout
    assert '--out\n/absolute/new run' in result.stdout
    assert subprocess.run(['bash', str(script)], env={**env, 'FAKE_EXIT': '75'}, capture_output=True).returncode == 75
    env.pop('INIT_FROM'); env.pop('SOURCE_SHA256')
    result = subprocess.run(['bash', str(script)], env={**env, 'RESUME': '1'}, capture_output=True, text=True)
    assert result.returncode == 0 and '--resume' in result.stdout
