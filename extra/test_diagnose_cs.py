"""Read-only diagnostic and saved witnesses; only synthetic test fixtures."""
from pathlib import Path
import subprocess

import pytest
import torch

from extra.diagnose_cs import Config, ascend, diagnose, read_config
from extra.acs_training import sha256
from extra.campaign_unit import Geometry
from extra.test_cs_training import setup
from extra.test_acs_training import one_thread
from extra.test_acs_pbs import shell_env

ROOT = Path(__file__).resolve().parents[1]


def tiny():
    return Config(uniform_count=8, circle_points=2, searches=2, pool_size=8,
                  starts=3, steps=2, chunk=2)


def test_ascent_monotone_and_stays_on_all_three_spheres():
    geo = Geometry('s4xs1xs1')
    starts = geo.sample(4)
    score = lambda x: 2 + x[:, 0] + .1*x[:, 5]
    end, scores = ascend(geo, score, starts, 10, 2)
    assert torch.all(scores >= score(starts))
    torch.testing.assert_close(scores, score(end))
    for block in end.split([5, 2, 2], dim=1):
        torch.testing.assert_close(block.norm(dim=1), torch.ones(4, dtype=torch.float64))


def test_diagnostic_saves_joint_point_metrics_without_mutating_inputs(tmp_path):
    _, base, audit, source = setup(tmp_path)
    paths = [base / 'final.pt', audit, source]
    before = [sha256(p) for p in paths]
    cfg = tiny()
    report = diagnose(cfg, source, base, audit, tmp_path / 'diagnostic', before[-1])
    assert report['status'] == 'completed_diagnostic'
    assert not report['numerically_validated_cs']
    assert report['uniform']['count'] == 8
    assert report['structured']['count'] == 80
    assert len(report['searches']) == 2
    assert all(s['nij_max_found'] >= s['pool_max_found'] - 1e-12 for s in report['searches'])
    points = torch.load(report['points_file'], weights_only=False)
    assert len(points['points']) == 94
    assert [sha256(p) for p in paths] == before
    assert report['combined_completed_groups']['acs_max'] < 1e-12
    for row in report['worst_points']:
        # Construct coordinates in the recorded float64 precision.
        matches = (points['points'] == torch.tensor(row['point'], dtype=torch.float64)).all(dim=1)
        index = matches.nonzero()[0].item()
        for key in ('nij', 'acs', 'tangency', 'G_condition', 'J_frobenius'):
            assert row[key] == points['metrics'][key][index].item()
    with pytest.raises(ValueError, match='empty'):
        diagnose(cfg, source, base, audit, tmp_path / 'diagnostic')


def test_rejects_hash_base_and_marks_interruption(tmp_path):
    _, base, audit, source = setup(tmp_path)
    with pytest.raises(ValueError, match='checksum'):
        diagnose(tiny(), source, base, audit, tmp_path / 'bad', '0'*64)
    ck = torch.load(source, weights_only=False)
    ck['audited_base']['sha256'] = '0'*64
    bad = tmp_path / 'bad_source.pt'
    torch.save(ck, bad)
    with pytest.raises(ValueError, match='different'):
        diagnose(tiny(), bad, base, audit, tmp_path / 'bad')
    report = diagnose(tiny(), source, base, audit, tmp_path / 'stopped', stop=lambda: True)
    assert report['status'] == 'incomplete'
    assert not report['numerically_validated_cs']
    assert 'points_file' not in report


def test_config_and_pbs_arguments(shell_env):
    cfg = read_config(ROOT / 'configs/ambient/s4xs1xs1_cs/diagnostic.yaml')
    assert (cfg.uniform_count, cfg.searches, cfg.starts, cfg.steps) == (16384, 3, 128, 100)
    script = ROOT / 'scratch/submit_S4xS1xS1_CS_diagnostic.sh'
    if not script.exists():
        pytest.skip('scratch scripts intentionally gitignored')
    python = Path(shell_env['PATH'].split(':')[0]) / 'python3'
    python.write_text('#!/bin/bash\nprintf "%s\\n" "$@"\nexit "${FAKE_EXIT:-0}"\n')
    env = {**shell_env, 'BASE_DIR': '/absolute/base', 'BASE_AUDIT': '/absolute/audit.json',
           'CHECKPOINT': '/absolute/CS final.pt', 'SOURCE_SHA256': 'a'*64, 'OUT': '/absolute/new report'}
    result = subprocess.run(['bash', str(script)], env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert '--checkpoint\n/absolute/CS final.pt' in result.stdout
    assert '--out\n/absolute/new report' in result.stdout
    assert subprocess.run(['bash', str(script)], env={**env, 'FAKE_EXIT': '75'}, capture_output=True).returncode == 75
