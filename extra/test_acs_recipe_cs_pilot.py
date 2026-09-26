"""Recipe replay and audited ACS-to-CS handoff; analytic fields are test-only."""
from dataclasses import asdict
import json
from pathlib import Path

import pytest
import torch

from extra.acs_recipe import read_recipe, run_recipe
from extra.acs_training import Config, Runner, read_config, sha256
from extra import cs_pilot
from extra.campaign_unit import Geometry, Net, load_acs_field, make_field, polar_corrected, nij_norm
from extra.test_acs_training import tiny, oracle_coefficients, one_thread

ROOT = Path(__file__).resolve().parents[1]


def test_production_recipe_matches_successful_hyperparameters():
    cfg = read_recipe(ROOT / 'configs/ambient/s4xs1xs1_acs/reproduce.yaml')
    pilot = read_config(ROOT / 'configs/ambient/s4xs1xs1_acs/pilot_trilin_uniform.yaml')
    polish = read_config(ROOT / 'configs/ambient/s4xs1xs1_acs/polish_trilin_uniform.yaml')
    assert asdict(cfg[0]) == asdict(pilot)
    for stage, seed in zip(cfg[1:], [107, 207, 307]):
        expected = {**asdict(polish), 'seed': seed}
        assert asdict(stage) == expected
    assert [c.seed for c in read_recipe(ROOT / 'configs/ambient/s4xs1xs1_acs/reproduce.yaml', 8)] == [8, 108, 208, 308]


def test_recipe_chains_best_checkpoints_and_resumes_completed_stages(tmp_path):
    configs = [tiny(poly_steps=1, refine_steps=1), tiny(seed=107, poly_steps=0, refine_steps=1)]
    report = run_recipe(configs, tmp_path)
    assert report['status'] == 'completed_pending_independent_audit'
    assert not report['numerically_validated']
    source = json.loads((tmp_path / '01_s107/manifest.json').read_text())['source']
    assert source['sha256'] == sha256(tmp_path / '00_s7/final.pt')
    assert source['optimizer_restarted']
    assert run_recipe(configs, tmp_path, resume=True) == report
    with pytest.raises(ValueError, match='empty'):
        run_recipe(configs, tmp_path)
    with pytest.raises(ValueError, match='same recipe'):
        run_recipe([tiny(seed=8)], tmp_path, resume=True)
    with (tmp_path / '00_s7/final.pt').open('ab') as stream:
        stream.write(b'changed')
    with pytest.raises(ValueError, match='changed'):
        run_recipe(configs, tmp_path, resume=True)


def test_recipe_wall_stop_is_not_completion_and_can_resume(tmp_path):
    configs = [tiny(poly_steps=1, refine_steps=1)]
    assert run_recipe(configs, tmp_path, wall_minutes=1e-12)['status'] == 'incomplete'
    assert run_recipe(configs, tmp_path, resume=True)['status'] == 'completed_pending_independent_audit'


def test_recipe_resumes_an_interrupted_stage_exactly(tmp_path, monkeypatch):
    configs = [tiny(poly_steps=2, refine_steps=2), tiny(seed=107, poly_steps=0, refine_steps=2)]
    full = run_recipe(configs, tmp_path / 'full')
    actual_run = Runner.run
    with monkeypatch.context() as patch:
        patch.setattr(Runner, 'run', lambda self: actual_run(self, max_updates=1))
        assert run_recipe(configs, tmp_path / 'split')['status'] == 'incomplete'
    split = run_recipe(configs, tmp_path / 'split', resume=True)
    a = torch.load(full['final_checkpoint'], weights_only=False)
    b = torch.load(split['final_checkpoint'], weights_only=False)
    for key in a['model']:
        torch.testing.assert_close(a['model'][key], b['model'][key], atol=0, rtol=0)


def fixture_base(tmp_path):
    base = tmp_path / 'base'
    runner = Runner(tiny(), base)
    oracle_coefficients(runner.model.poly)
    torch.save(runner.payload(), base / 'final.pt')
    report = tmp_path / 'fixture_audit.json'
    # Synthetic report tests the handoff contract, not a real learned-field audit.
    report.write_text(json.dumps({'case': 's4xs1xs1', 'audit_tier': 'full',
                                 'numerically_validated': True, 'screening_targets_passed': True,
                                 'sha256': sha256(base / 'final.pt')}))
    return base, report


def test_cs_rejects_quick_failed_or_changed_base(tmp_path):
    base, report = fixture_base(tmp_path)
    original = json.loads(report.read_text())
    cs_pilot.checked_base(base, report)
    for key, value in [('audit_tier', 'quick'), ('numerically_validated', False), ('case', 's2'), ('sha256', '0'*64)]:
        report.write_text(json.dumps({**original, key: value}))
        with pytest.raises(ValueError):
            cs_pilot.checked_base(base, report)


@pytest.mark.parametrize('stages', [1, 2])
def test_corrected_frozen_base_identity_initialization_and_gradient(tmp_path, stages):
    base, _ = fixture_base(tmp_path)
    geo = Geometry('s4xs1xs1')
    field = polar_corrected(geo, load_acs_field(geo, base))
    net = Net(9, 9, stages, [8], .5).double()
    with torch.no_grad():
        net.net[-1].weight.zero_()
        net.net[-1].bias.zero_()
    deformed, _ = make_field(geo, 'cs', net, base_fn=field)
    x = geo.sample(2)
    torch.testing.assert_close(deformed(x), field(x), atol=1e-14, rtol=0)
    before = sha256(base / 'final.pt')
    loss = nij_norm(geo, deformed, x).square().mean()
    loss.backward()
    assert all(p.grad is not None and torch.isfinite(p.grad).all() for p in net.parameters())
    assert net.net[-1].weight.grad.norm() > 0
    torch.optim.Adam(net.parameters(), lr=1e-4).step()
    J = deformed(x)
    assert (J @ J + geo.manifold.projector(x)).norm(dim=(1, 2)).max() < 1e-12
    assert sha256(base / 'final.pt') == before


def test_real_cs_pilot_one_step_records_audited_base(tmp_path):
    base, report = fixture_base(tmp_path)
    cfg = cs_pilot.read_config(ROOT / 'configs/ambient/s4xs1xs1_cs/smoke.yaml')
    cfg.update(steps=1, hidden=[8], batch=2, n_initial=2, n_validation=2,
               n_dense=2, n_chart=2, adv_starts=2, adv_steps=2, eval_chunk=2)
    result = cs_pilot.run(cfg, base, report, tmp_path / 'cs')
    assert result['status'] == 'pilot_completed_pending_review'
    assert result['dense']['acs_sup'] < 1e-8
    assert result['dense']['acs_floorfree_sup'] < 1e-12
    assert result['dense']['nij_mean'] > 0
    assert not result['numerically_validated_cs']
    ck = torch.load(tmp_path / 'cs/final.pt', weights_only=False)
    assert ck['audited_base']['sha256'] == sha256(base / 'final.pt')
    with pytest.raises(ValueError, match='new output'):
        cs_pilot.run(cfg, base, report, tmp_path / 'cs')


def test_cs_wall_stop_does_not_publish_final(tmp_path):
    base, report = fixture_base(tmp_path)
    cfg = cs_pilot.read_config(ROOT / 'configs/ambient/s4xs1xs1_cs/smoke.yaml')
    cfg.update(steps=2, hidden=[8], batch=2, n_initial=2, n_validation=2,
               n_dense=2, n_chart=2, adv_starts=2, adv_steps=2, wall_minutes=1e-12)
    result = cs_pilot.run(cfg, base, report, tmp_path / 'cs')
    assert result['status'] == 'incomplete'
    assert not (tmp_path / 'cs/final.pt').exists()
    assert (tmp_path / 'cs/partial.pt').exists()


def test_cs_full_width_base_and_smoke_network_gradient(tmp_path):
    cfg = read_config(ROOT / 'configs/ambient/s4xs1xs1_acs/polish_trilin_uniform.yaml')
    runner = Runner(cfg, tmp_path / 'base')
    oracle_coefficients(runner.model.poly)
    runner.model.refinement_phase(zero_output=True)
    torch.save(runner.payload(), runner.out / 'final.pt')
    geo = Geometry('s4xs1xs1')
    base = polar_corrected(geo, load_acs_field(geo, runner.out))
    cs_cfg = cs_pilot.read_config(ROOT / 'configs/ambient/s4xs1xs1_cs/smoke.yaml')
    model = Net(9, 9, 1, cs_cfg['hidden'], .5).double()
    with torch.no_grad():
        model.net[-1].weight.zero_()
        model.net[-1].bias.zero_()
    field, _ = make_field(geo, 'cs', model, base_fn=base)
    nij_norm(geo, field, geo.sample(cs_cfg['batch'])).square().mean().backward()
    assert all(torch.isfinite(p.grad).all() for p in model.parameters())
