"""Focused regression tests for arbitrary-factor ambient campaigns."""

import itertools

import pytest
import torch

from extra.campaign_unit import (
    CurriculumNet,
    Geometry,
    Net,
    PolyNet,
    acs_residual,
    canonical_J0,
    chart_field,
    chart_protocol,
    load_acs_field,
    make_field,
    polar_corrected,
    product_chart_overlap,
    unit_loss,
)
from src.geometry.ambient import AmbientManifold
from src.geometry.sphere import EvenSphere


DTYPE = torch.float64


def _exact_s4xt2_acs(p):
    """An exact ACS from TS^4 ⊕ R^2 ≅ R^6, used only as a test oracle."""
    batch = len(p)
    x, u, v = p[:, :5], p[:, 5:7], p[:, 7:9]
    P4 = torch.eye(5, dtype=p.dtype, device=p.device).expand(batch, 5, 5)
    P4 = P4 - x.unsqueeze(2) * x.unsqueeze(1)
    tu = torch.stack([-u[:, 1], u[:, 0]], dim=1)
    tv = torch.stack([-v[:, 1], v[:, 0]], dim=1)

    # The first five rows trivialise TS^4 plus the first circle direction:
    # e_a -> (P_x e_a, x_a t_u).  The final row is the second circle tangent.
    rows = [
        torch.cat([P4[:, a, :], x[:, a:a + 1] * tu, p.new_zeros(batch, 2)], dim=1)
        for a in range(5)
    ]
    rows.append(torch.cat([p.new_zeros(batch, 7), tv], dim=1))
    frame = torch.stack(rows, dim=1)
    J0 = canonical_J0(6).to(dtype=p.dtype, device=p.device)
    return frame.transpose(1, 2) @ J0 @ frame


def test_three_factor_geometry_and_projector():
    manifold = AmbientManifold((4, 1, 1))
    x = manifold.sample(12)
    P = manifold.projector(x)
    assert manifold.dim == 6
    assert manifold.ambient_dim == 9
    assert x.shape == (12, 9)
    assert P.shape == (12, 9, 9)
    assert torch.allclose(P @ P, P, atol=1e-12, rtol=0.0)
    assert torch.allclose(P.transpose(1, 2), P, atol=1e-12, rtol=0.0)
    assert torch.allclose(P.diagonal(dim1=1, dim2=2).sum(1), x.new_full((12,), 6.0))


def test_three_factor_pairwise_bilinear_features():
    # Cross-factor counts: 5*2 + 5*2 + 2*2 = 24.
    model = PolyNet(9, 9, out_scale=0.1, feat='bilin', blocks=[5, 2, 2]).to(DTYPE)
    assert model.lin.in_features == 24
    assert model(torch.randn(3, 9, dtype=DTYPE)).shape == (3, 1, 9, 9)


def test_registered_s4xs1xs1_has_exact_test_oracle():
    geo = Geometry('s4xs1xs1', 'ambient')
    x = geo.sample(10)
    assert geo.dims == (4, 1, 1)
    assert geo.d == 6
    assert geo.N == 9
    assert acs_residual(geo, _exact_s4xt2_acs, x).max().item() < 2e-10


def test_three_factor_polar_correction_produces_exact_acs():
    geo = Geometry('s4xs1xs1', 'ambient')

    def perturbed(p):
        scale = 1.0 + 0.05 * torch.sin(p[:, :1])
        return scale.unsqueeze(2) * _exact_s4xt2_acs(p)

    x = geo.sample(32)
    assert acs_residual(geo, perturbed, x).max().item() < 1.0
    corrected = polar_corrected(geo, perturbed)
    assert acs_residual(geo, corrected, x).max().item() < 2e-10


def test_three_factor_chart_transitions_are_consistent():
    geo = Geometry('s4xs1xs1', 'ambient')
    batch = 12
    factor_coords = []
    generator = torch.Generator().manual_seed(123)
    for m in geo.dims:
        x0 = torch.randn(batch, m, dtype=DTYPE, generator=generator)
        x0 = 0.4 * x0 / x0.norm(dim=1, keepdim=True)
        factor_coords.append((x0, EvenSphere.patch_change_coords(x0)))

    patches = list(itertools.product((0, 1), repeat=3))
    J_dict = {}
    for patch in patches:
        xc = torch.cat([factor_coords[a][patch[a]] for a in range(3)], dim=1)
        J_dict[patch] = chart_field(_exact_s4xt2_acs, geo.dims, patch)(xc)

    overlap = product_chart_overlap(J_dict, factor_coords, geo.dims, patches)
    assert overlap.item() < 2e-10


def test_three_factor_full_chart_protocol_runs():
    geo = Geometry('s4xs1xs1', 'ambient')
    chart_nij, overlap = chart_protocol(geo, _exact_s4xt2_acs, n_test=8)
    assert chart_nij >= 0.0
    assert torch.isfinite(torch.tensor(chart_nij))
    assert overlap < 2e-10


def test_three_factor_cs_parametrisation_is_exact_and_differentiable():
    geo = Geometry('s4xs1xs1', 'ambient')
    model = Net(geo.N, geo.N, stages=2, hidden=[8], out_scale=0.01).to(DTYPE)
    J_fn, meta = make_field(geo, 'cs', model, base_fn=_exact_s4xt2_acs)
    x = geo.sample(2)
    loss, _ = unit_loss(geo, 'cs', J_fn, x)
    loss.backward()

    assert meta['exact_acs']
    assert torch.isfinite(loss)
    # acs_residual includes a deliberate sqrt(1e-20) reporting floor.
    assert acs_residual(geo, J_fn, x).max().item() < 2e-10
    grads = [parameter.grad for parameter in model.parameters()]
    assert all(grad is not None and torch.isfinite(grad).all() for grad in grads)


def test_three_factor_acs_checkpoint_round_trip(tmp_path):
    geo = Geometry('s4xs1xs1', 'ambient')
    hidden = [8]
    model = CurriculumNet(
        geo.N, geo.N, hidden, feat='lin', blocks=[m + 1 for m in geo.dims]
    ).to(DTYPE)
    torch.save({
        'model': model.state_dict(),
        'case': geo.case,
        'task': 'acs',
        'hidden': hidden,
        'stages': 1,
        'kind': 'ambient',
        'model_type': 'curriculum',
        'feat': 'lin',
    }, tmp_path / 'final.pt')

    loaded = load_acs_field(geo, str(tmp_path))
    x = geo.sample(4)
    original, _ = make_field(geo, 'acs', model)
    assert torch.allclose(loaded(x), original(x), atol=0.0, rtol=0.0)

    wrong = torch.load(tmp_path / 'final.pt', map_location='cpu', weights_only=False)
    wrong['case'] = 's4xs2'
    torch.save(wrong, tmp_path / 'final.pt')
    with pytest.raises(ValueError, match='expected'):
        load_acs_field(geo, str(tmp_path))

    wrong['case'] = geo.case
    wrong['task'] = 'cs'
    torch.save(wrong, tmp_path / 'final.pt')
    with pytest.raises(ValueError, match="expected 'acs'"):
        load_acs_field(geo, str(tmp_path))
