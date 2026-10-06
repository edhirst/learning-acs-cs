"""Tests for the polynomial free-field CS mode, the LM finisher, the
orientation-controlled ACS attempt and the rational conjugators."""

import pytest
import torch

from extra.campaign_unit import (
    Geometry,
    Net,
    PolyNet,
    acs_residual,
    base_acs,
    cs_nij_norm_fast,
    make_conjugation,
    make_field,
    nij_norm,
    orientation_index,
    pfaffian,
)
from extra.test_ambient_products import _exact_s4xt2_acs
from src.loss.ambient_loss import compute_J_and_DJ
from src.models.poly_field import (
    PolyCache,
    PolyField,
    monomial_keys,
    poly_K_DK,
    poly_residual,
    poly_residual_and_jacobian,
    poly_split_losses,
)
from src.optim.lm import jacobian_fwd, lm_minimise

DTYPE = torch.float64


def _Jstd(n):
    J = torch.zeros(n, n, dtype=DTYPE)
    for k in range(n // 2):
        J[2 * k, 2 * k + 1], J[2 * k + 1, 2 * k] = -1.0, 1.0
    return J


def _hopf_M(p, n):
    """Hopf structure on S^{n-1} x S^1 (test oracle): J dt = xi, J xi = -dt, J|H = i."""
    x, y = p[:, :n], p[:, n:]
    xi, dt = x @ _Jstd(n).T, y @ _Jstd(2).T
    B = len(p)
    top = torch.cat([_Jstd(n).expand(B, n, n), xi.unsqueeze(2) * dt.unsqueeze(1)], dim=2)
    bot = torch.cat([-dt.unsqueeze(2) * xi.unsqueeze(1), torch.zeros(B, 2, 2, dtype=DTYPE)], dim=2)
    return torch.cat([top, bot], dim=1)


def _ce_M(p):
    """Calabi-Eckmann structure on S^3 x S^3 (test oracle)."""
    x1, x2 = p[:, :4], p[:, 4:]
    J0 = _Jstd(4)
    xi1, xi2 = x1 @ J0.T, x2 @ J0.T
    B = len(p)
    top = torch.cat([J0.expand(B, 4, 4), -xi1.unsqueeze(2) * xi2.unsqueeze(1)], dim=2)
    bot = torch.cat([xi2.unsqueeze(2) * xi1.unsqueeze(1), J0.expand(B, 4, 4)], dim=2)
    return torch.cat([top, bot], dim=1)


# ---------------------------------------------------------------------------
# polynomial field
# ---------------------------------------------------------------------------

def test_monomial_classes_nested_and_trilinear():
    k2, k3 = monomial_keys(9, 2), monomial_keys(9, 3)
    assert len(k2) == 55 and set(k2) <= set(k3)
    k2t = monomial_keys(9, 2, blocks=[5, 2, 2], trilin=True)
    assert len(k2t) == 55 + 5 * 2 * 2
    assert set(k2t) <= set(k3)
    with pytest.raises(ValueError):
        monomial_keys(8, 2, blocks=[5, 3], trilin=True)


def test_load_from_preserves_the_field():
    torch.manual_seed(0)
    geo = Geometry('s4xs1xs1', 'ambient')
    small = PolyField(geo.N, 2, init_scale=0.5)
    big = PolyField(geo.N, 2, blocks=[5, 2, 2], trilin=True)
    big.load_from(small)
    x = geo.sample(20)
    assert torch.allclose(small(x), big(x), atol=1e-14, rtol=0)
    bigger = PolyField(geo.N, 3)
    bigger.load_from(big)
    assert torch.allclose(bigger(x), big(x), atol=1e-14, rtol=0)


@pytest.mark.parametrize('case,deg,trilin', [('s3xs1', 2, False), ('s6', 3, False),
                                             ('s4xs1xs1', 2, True)])
def test_explicit_derivative_and_residual_match_autodiff(case, deg, trilin):
    torch.manual_seed(1)
    geo = Geometry(case, 'ambient')
    man = geo.manifold
    f = PolyField(geo.N, deg, init_scale=0.5, blocks=[m + 1 for m in geo.dims], trilin=trilin)
    x = man.sample(24)
    c = PolyCache(man, f, x)
    K, DK = poly_K_DK(f, c)
    K_fn = lambda p: man.projector(p) @ f(p)[:, 0] @ man.projector(p)
    K2, DK2 = compute_J_and_DJ(K_fn, x)
    assert torch.allclose(K, K2, atol=1e-13, rtol=0)
    assert torch.allclose(DK, DK2, atol=1e-12, rtol=0)
    n2, c2 = poly_split_losses(f, c)
    ref_n = (nij_norm(geo, K_fn, x) ** 2).mean()
    ref_c = (acs_residual(geo, K_fn, x) ** 2).mean()
    assert torch.allclose(n2, ref_n, rtol=1e-10, atol=1e-15)
    assert torch.allclose(c2, ref_c, rtol=1e-10, atol=1e-15)
    r = poly_residual(f, c, w_acs=10.0)
    assert torch.allclose((r ** 2).sum(), ref_n + 10.0 * ref_c, rtol=1e-10)


def test_factorised_jacobian_matches_forward_ad():
    torch.manual_seed(2)
    geo = Geometry('s2xs2', 'ambient')
    f = PolyField(geo.N, 2, init_scale=0.5)
    c = PolyCache(geo.manifold, f, geo.manifold.sample(8))
    J1 = jacobian_fwd(lambda th: poly_residual(f, c, th.view_as(f.W)), f.W.detach().reshape(-1), 96)
    r2, J2 = poly_residual_and_jacobian(f, c)
    assert torch.allclose(poly_residual(f, c).detach(), r2, atol=1e-14, rtol=0)
    assert torch.allclose(J1, J2, atol=1e-13, rtol=0)


@pytest.mark.parametrize('case,Mf', [('s3xs1', lambda p: _hopf_M(p, 4)),
                                     ('s5xs1', lambda p: _hopf_M(p, 6)),
                                     ('s3xs3', _ce_M)])
def test_classical_complex_structures_lie_in_the_degree_two_class(case, Mf):
    """Offline oracle check (never used in training): the Hopf and Calabi-Eckmann
    structures are exact integrable ACS and exactly degree-2 ambient fields."""
    geo = Geometry(case, 'ambient')
    man = geo.manifold
    K_fn = lambda p: man.projector(p) @ Mf(p) @ man.projector(p)
    x = man.sample(64)
    assert acs_residual(geo, K_fn, x).max() < 1e-9
    assert nij_norm(geo, K_fn, x).max() < 1e-9
    f = PolyField(geo.N, 2)
    xf = man.sample(600)
    F = f.features(xf)
    W = torch.linalg.lstsq(F, Mf(xf).reshape(len(xf), -1)).solution
    assert (F @ W - Mf(xf).reshape(len(xf), -1)).abs().max() < 1e-12


# ---------------------------------------------------------------------------
# Levenberg-Marquardt
# ---------------------------------------------------------------------------

def test_lm_converges_quadratically_on_a_zero_residual_toy():
    res = lambda th, c: torch.stack([10.0 * (th[1] - th[0] ** 2), 1.0 - th[0]])
    th, log = lm_minimise(res, torch.tensor([-1.2, 1.0], dtype=DTYPE), [None], n_iter=60, lam0=1e-3)
    assert torch.allclose(th, torch.ones(2, dtype=DTYPE), atol=1e-10)
    assert log[-1]['loss'] < 1e-20


def test_lm_polishes_a_perturbed_integrable_structure_to_machine_precision():
    torch.manual_seed(3)
    geo = Geometry('s2xs2', 'ambient')
    man = geo.manifold
    f = PolyField(geo.N, 2)
    xf = man.sample(500)
    F = f.features(xf)
    T = base_acs('s2xs2', xf, man.projector).reshape(len(xf), -1)
    with torch.no_grad():
        f.W.copy_(torch.linalg.lstsq(F, T).solution.T)
        f.W.add_(1e-3 * torch.randn_like(f.W))
    caches = [PolyCache(man, f, man.sample(32)) for _ in range(2)]
    res = lambda W, c: poly_residual(f, c, W)
    jac = lambda W, c: poly_residual_and_jacobian(f, c, W)
    W, log = lm_minimise(res, f.W.detach(), caches, n_iter=20, lam0=1e-3, jac_fn=jac)
    assert log[0]['loss'] > 1e-8 and log[-1]['loss'] < 1e-22
    with torch.no_grad():
        f.W.copy_(W)
    n2, c2 = poly_split_losses(f, PolyCache(man, f, man.sample(64)))
    assert n2 < 1e-18 and c2 < 1e-18                     # held out, not just the fit points


# ---------------------------------------------------------------------------
# orientation index
# ---------------------------------------------------------------------------

def test_pfaffian_squares_to_the_determinant():
    torch.manual_seed(4)
    for d in (2, 4, 6):
        A = torch.randn(7, d, d, dtype=DTYPE)
        A = A - A.transpose(1, 2)
        assert torch.allclose(pfaffian(A) ** 2, torch.linalg.det(A), rtol=1e-10, atol=1e-12)


@pytest.mark.parametrize('case', ['s2xs2', 's6', 's4xs1xs1'])
def test_orientation_index_of_an_exact_acs(case):
    geo = Geometry(case, 'ambient')
    proj = geo.manifold.projector
    J = _exact_s4xt2_acs if case == 's4xs1xs1' else (lambda p: base_acs(case, p, proj))
    x = geo.sample(400)
    o = orientation_index(geo, J, x)
    assert o.abs().min() >= 1.0 - 1e-9                   # |o| >= 1 for every exact ACS
    assert (torch.sign(o) == torch.sign(o[0])).all()     # one orientation everywhere
    o_neg = orientation_index(geo, lambda p: -J(p), x)
    assert torch.allclose(o_neg, (-1.0) ** (geo.d // 2) * o, atol=1e-12)


def test_three_factor_class_contains_the_frame_acs_of_s4xt2():
    """The frame ACS of S4 x T2 (test oracle) is P M P with M = G^T J0 G, where
    frame = G P and G is LINEAR in p; the (x u)(x u) terms of M cancel by the
    antisymmetry of J0, leaving constant, linear, bilinear and three-factor
    trilinear monomials: M lies exactly in the 'all2+trilin' ladder class."""
    geo = Geometry('s4xs1xs1', 'ambient')
    x = geo.sample(600)
    B = len(x)
    xs, u, v = x[:, :5], x[:, 5:7], x[:, 7:9]
    tu = torch.stack([-u[:, 1], u[:, 0]], dim=1)
    tv = torch.stack([-v[:, 1], v[:, 0]], dim=1)
    eye5 = torch.eye(5, dtype=DTYPE)
    rows = [torch.cat([eye5[a].expand(B, 5), xs[:, a:a + 1] * tu, x.new_zeros(B, 2)], dim=1)
            for a in range(5)]
    rows.append(torch.cat([x.new_zeros(B, 7), tv], dim=1))
    G = torch.stack(rows, dim=1)                                          # (B, 6, 9)
    J0 = torch.zeros(6, 6, dtype=DTYPE)
    for k in range(3):
        J0[2 * k, 2 * k + 1], J0[2 * k + 1, 2 * k] = -1.0, 1.0
    M = G.transpose(1, 2) @ J0 @ G
    P = geo.manifold.projector(x)
    assert torch.allclose(P @ M @ P, _exact_s4xt2_acs(x), atol=1e-13, rtol=0)
    pn = PolyNet(geo.N, geo.N, 0.3, feat='all2+trilin', blocks=[5, 2, 2]).to(DTYPE)
    F = pn.features(x)
    W = torch.linalg.lstsq(F, M.reshape(B, -1)).solution
    assert (F @ W - M.reshape(B, -1)).abs().max() < 1e-11
    F2 = PolyNet(geo.N, geo.N, 0.3, feat='all2', blocks=[5, 2, 2]).to(DTYPE).features(x)
    W2 = torch.linalg.lstsq(F2, M.reshape(B, -1)).solution
    assert (F2 @ W2 - M.reshape(B, -1)).abs().max() > 1e-2           # this M needs the trilinear terms


# ---------------------------------------------------------------------------
# rational conjugators and the fast CS loss
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('conj', ['cayley', 'affine'])
def test_rational_conjugator_forward_derivative_on_every_batch_element(conj):
    """Regression: torch.linalg.solve under vmap(jacfwd) returned wrong
    derivatives for batch elements > 0; the conjugators must use inv."""
    torch.manual_seed(5)
    geo = Geometry('s6', 'ambient')
    proj = geo.manifold.projector
    model = Net(geo.N, geo.N, 2, [32, 32], 0.05).to(DTYPE)
    J_fn, _ = make_field(geo, 'cs', model, base_fn=lambda p: base_acs('s6', p, proj), conj=conj)
    x = geo.sample(4)
    _, DJ = compute_J_and_DJ(J_fn, x)
    eps = 1e-6
    for i in range(4):
        for k in (0, 3):
            xp = x[i:i + 1].clone()
            xp[0, k] += eps
            fd = (J_fn(xp) - J_fn(x[i:i + 1]))[0] / eps
            assert (DJ[i, :, :, k] - fd).abs().max() < 1e-5


@pytest.mark.parametrize('conj', ['exp', 'affine'])
def test_fast_cs_loss_matches_reference(conj):
    torch.manual_seed(6)
    geo = Geometry('s6', 'ambient')
    proj = geo.manifold.projector
    base = lambda p: base_acs('s6', p, proj)
    model = Net(geo.N, geo.N, 2, [32, 32], 0.05).to(DTYPE)
    J_fn, _ = make_field(geo, 'cs', model, base_fn=base, conj=conj)
    g_fn = make_conjugation(geo, model, conj)
    x = geo.sample(16)
    ref = nij_norm(geo, J_fn, x).detach()
    fast = cs_nij_norm_fast(geo, g_fn, base, x).detach()
    assert torch.allclose(fast, ref, rtol=1e-10)
    assert acs_residual(geo, J_fn, x).max() < 1e-9       # still an exact ACS


# ---------------------------------------------------------------------------
# exact-ACS measurement modes: polar residual, matrix-free J, polynomial conjugator
# ---------------------------------------------------------------------------

def _s6_admissible_field(deg):
    from src.models.poly_field import fit_field
    geo = Geometry('s6', 'ambient')
    proj = geo.manifold.projector
    f = PolyField(geo.N, deg)
    f.orthonormalize(geo.manifold, n=4000)
    fit_field(f, lambda p: base_acs('s6', p, proj), geo.manifold, n=2000)
    with torch.no_grad():
        f.W.add_(1e-2 * torch.randn_like(f.W))
    return geo, f


def test_orthonormal_basis_removes_sphere_redundancy():
    geo = Geometry('s6', 'ambient')
    f = PolyField(geo.N, 3)
    assert f.orthonormalize(geo.manifold, n=6000) == 1 + 7 + 27 + 77     # harmonics of degree <= 3


def test_polar_residual_and_jacobians_match_reference():
    from extra.campaign_unit import polar_corrected
    from src.models.poly_field import (polar_split, polar_residual, polar_residual_and_jacobian,
                                       polar_factorised)
    torch.manual_seed(7)
    geo, f = _s6_admissible_field(2)
    man = geo.manifold
    x = man.sample(8)
    c = PolyCache(man, f, x)
    K = lambda p: man.projector(p) @ f(p)[:, 0] @ man.projector(p)
    ref = (nij_norm(geo, polar_corrected(geo, K), x) ** 2).mean()
    n2, _ = polar_split(f, c)
    assert torch.allclose(n2, ref, rtol=1e-9)
    r, J = polar_residual_and_jacobian(f, c)
    J2 = jacobian_fwd(lambda th: polar_residual(f, c, th.view_as(f.W)), f.W.detach().reshape(-1), 128)
    assert torch.allclose(J, J2, atol=1e-11, rtol=0)
    fac = polar_factorised(f, c)
    v, u = torch.randn(J.shape[1], dtype=DTYPE), torch.randn(J.shape[0], dtype=DTYPE)
    assert torch.allclose(fac.matvec(v), J @ v, atol=1e-11, rtol=0)
    assert torch.allclose(fac.rmatvec(u), J.T @ u, atol=1e-11, rtol=0)


def test_polynomial_conjugator_residual_is_exact_and_jacobian_matches():
    from src.models.poly_field import ConjCache, conj_residual, conj_residual_and_jacobian
    torch.manual_seed(8)
    geo = Geometry('s6', 'ambient')
    man = geo.manifold
    base = lambda p: base_acs('s6', p, man.projector)
    f = PolyField(geo.N, 1)
    f.orthonormalize(man, n=2000)
    with torch.no_grad():
        f.W.normal_(0, 0.1)
    x = man.sample(6)
    c = ConjCache(man, f, x, base)
    Jf = lambda p: (lambda B: torch.matrix_exp(B) @ base(p) @ torch.matrix_exp(-B))(
        man.projector(p) @ f(p)[:, 0] @ man.projector(p))
    assert acs_residual(geo, Jf, x).max() < 1e-9
    ref = (nij_norm(geo, Jf, x) ** 2).mean()
    assert torch.allclose((conj_residual(f, c) ** 2).sum(), ref, rtol=1e-9)
    r, J = conj_residual_and_jacobian(f, c)
    J2 = jacobian_fwd(lambda th: conj_residual(f, c, th.view_as(f.W)), f.W.detach().reshape(-1), 64)
    assert torch.allclose(J, J2, atol=1e-10, rtol=0)
