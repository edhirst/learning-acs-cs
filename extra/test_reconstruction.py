"""Tests for the period-map PINN and the reconstructed complex structures on
S^6 and S^2 x S^4 (extra/period_map.py, extra/s6_reconstruction.py)."""

import numpy as np
import pytest
import torch

from extra.period_map import (A1, A2, T1, T2, TriangleGroup, collocation, D_function, edge_pairs,
                              functional_residuals, mobius, sample_ford)
from extra.reconstruct_s6 import cusp_check, learn_period_map
from extra.s6_reconstruction import (J_field, J_of_Pi, TorchPeriodMap, blended_field, deck_differential,
                                     mobius_t, nijenhuis)

DTYPE = torch.float64
RHO = np.exp(1j * np.pi / 3)


@pytest.fixture(scope='module')
def solved(tmp_path_factory):
    G, pm = learn_period_map(K=20, c0=-3j, log=lambda *a: None)
    path = tmp_path_factory.mktemp('pm') / 'period_map.npz'
    np.savez(path, a=pm.a, b=pm.b, cb=pm.cb, K=pm.K, q0=pm.q0)
    return G, pm, TorchPeriodMap(str(path))


def test_triangle_group_relations():
    G = TriangleGroup()
    assert all(v < 1e-12 for v in G.check().values())
    mp = np.linalg.matrix_power
    assert np.allclose(mp(T1, 3), np.eye(4)) and np.allclose(mp(T2, 4), np.eye(4))
    N = np.linalg.inv(T1 @ T2) - np.eye(4)
    assert np.allclose(N @ N, 0) and np.linalg.matrix_rank(N) == 2


def test_period_map_reproduces_the_fixed_point_values_without_being_told(solved):
    G, pm, _ = solved
    assert abs(pm.tau(np.array([G.Z1]))[0] - RHO) < 1e-12
    assert abs(pm.tau(np.array([G.Z2]))[0] - 1j) < 1e-12
    assert abs(pm.mu(np.array([G.Z1]))[0] - 1 / (1 + RHO)) < 1e-12        # [AC26, Thm 3.4(ii)]
    assert abs(pm.mu(np.array([G.Z2]))[0] - (1 - 1j) / 2) < 1e-12


def test_functional_equations_hold_on_fresh_gluing_pairs(solved):
    G, pm, _ = solved
    pairs = edge_pairs(G, 800, np.random.default_rng(99))
    r1, r2 = functional_residuals(pm, G, pairs['g1']), functional_residuals(pm, G, pairs['g2'])
    assert max(v for k, v in r1.items() if k.endswith('g1')) < 1e-10
    assert max(v for k, v in r2.items() if k.endswith('g2')) < 1e-10
    assert D_function(pm, sample_ford(G, 2000, np.random.default_rng(3))).max() < 0      # lattices


@pytest.mark.parametrize('lattice', [None, np.diag([1, 1, 1, 0.5])])
def test_reconstructed_structure_is_integrable_and_glues(solved, lattice):
    G, _, pmt = solved
    L = None if lattice is None else torch.tensor(lattice, dtype=DTYPE)
    Ap = {1: A1, 2: A2} if lattice is None else {1: np.linalg.inv(lattice) @ A1 @ lattice,
                                                  2: np.linalg.inv(lattice) @ A2 @ lattice}
    Zf = sample_ford(G, 300, np.random.default_rng(4))
    ab = torch.tensor(np.stack([Zf.real, Zf.imag], 1), dtype=DTYPE)
    J = J_field(pmt, ab, L)
    assert (J @ J + torch.eye(6, dtype=DTYPE)).abs().max() < 1e-12
    assert nijenhuis(lambda p: J_field(pmt, p, L), ab).abs().max() < 1e-12
    pairs = edge_pairs(G, 300, np.random.default_rng(5))
    for j, g in ((1, G.g1), (2, G.g2)):
        Z = torch.tensor(pairs[f'g{j}'], dtype=torch.complex128)
        D = deck_differential(g, Z, Ap[j])
        err = (D @ J_of_Pi(pmt.Pi(Z), L) @ torch.linalg.inv(D) - J_of_Pi(pmt.Pi(mobius_t(g, Z)), L)).abs().max()
        assert err < 1e-10
        Dw = deck_differential(g, Z, np.linalg.inv(Ap[j]))                      # wrong transport fails
        assert (Dw @ J_of_Pi(pmt.Pi(Z), L) @ torch.linalg.inv(Dw) - J_of_Pi(pmt.Pi(mobius_t(g, Z)), L)).abs().max() > 1e-2
        Jb = blended_field(pmt, g, Ap[j], centre=np.median(pairs[f'g{j}'].real), lattice=L, Alat=Ap[j])
        abz = torch.tensor(np.stack([pairs[f'g{j}'].real, pairs[f'g{j}'].imag], 1)[:100], dtype=DTYPE)
        assert nijenhuis(Jb, abz).abs().max() < 1e-9


def test_superlattice_monodromy_matches_guo_fang_lu():
    D = np.diag([1, 1, 1, 0.5])
    A1p, A2p = np.linalg.inv(D) @ A1 @ D, np.linalg.inv(D) @ A2 @ D
    assert np.array_equal(np.round(A1p), [[1, 0, 0, 0], [6, 0, 1, 0], [-6, -1, -1, 0], [-4, 2, 0, 1]])
    assert np.array_equal(np.round(A2p), [[1, 0, 0, 0], [0, 0, -1, 0], [-6, 1, 0, 0], [6, 0, 2, 1]])
    assert np.allclose(A1p, np.round(A1p)) and np.allclose(A2p, np.round(A2p))


def test_cusp_normal_form_satisfies_the_toric_filling_hypothesis(solved):
    _, pm, _ = solved
    c = cusp_check(pm)
    assert c['C_single_valued_err'] < 1e-12 and np.isfinite(c['R_inf']) and c['eps0'] > 0
