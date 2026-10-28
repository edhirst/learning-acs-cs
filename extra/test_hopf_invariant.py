"""Tests for the homotopy class of S^3 x S^1 almost complex structures relative to the Hopf
structure (extra/hopf_invariant.py)."""

import numpy as np
import pytest
import torch

from extra.hopf_invariant import C_MAIN, hopf_invariant, relative_class, sphere_mesh
from src.geometry.hopf_ambient import hopf_ambient_acs


@pytest.fixture(scope='module')
def mesh():
    return sphere_mesh(16)


def hopf_map(q):
    z1, z2 = q[:, 0] + 1j * q[:, 1], q[:, 2] + 1j * q[:, 3]
    w = z1 * np.conj(z2)
    f = np.stack([np.abs(z1) ** 2 - np.abs(z2) ** 2, 2 * w.real, 2 * w.imag], 1)
    return f / np.linalg.norm(f, axis=1, keepdims=True)


def test_mesh_is_a_closed_triangulated_sphere(mesh):
    X, tets = mesh
    assert np.allclose(np.linalg.norm(X, axis=1), 1)
    faces = np.sort(np.concatenate([np.delete(tets, i, axis=1) for i in range(4)]), axis=1)
    _, counts = np.unique(faces, axis=0, return_counts=True)
    assert (counts == 2).all()                                  # every triangle bounds exactly two tets


@pytest.mark.parametrize('pre, expected', [(lambda q: q, 1), (lambda q: q * [1, -1, -1, -1], -1),
                                           (lambda q: np.stack([q[:, 0] ** 2 - (q[:, 1:] ** 2).sum(1),
                                                                *(2 * q[:, 0:1] * q[:, 1:]).T], 1), 2)])
def test_linking_number_of_preimages_is_the_hopf_invariant(mesh, pre, expected):
    X, tets = mesh
    q = pre(X)
    lk = [t['linking'] for t in hopf_invariant(X, hopf_map(q / np.linalg.norm(q, axis=1, keepdims=True)), tets)]
    assert np.allclose(lk, expected, atol=1e-6)


def test_null_homotopic_map_has_invariant_zero(mesh):
    X, tets = mesh
    f = np.c_[X[:, :2], 2 + X[:, 2]]
    f /= np.linalg.norm(f, axis=1, keepdims=True)
    assert all(abs(t['linking']) < 1e-9 for t in hopf_invariant(X, f, tets))


def test_relative_class_of_hopf_structures(mesh):
    Jt = lambda p: hopf_ambient_acs(p, 2, C_MAIN)
    S = torch.diag(torch.tensor([1., -1, -1, -1, 1, -1], dtype=torch.float64))   # (q, theta) -> (conj q, -theta)
    assert relative_class(Jt, Jt, mesh=mesh)[0] == 0
    assert relative_class(lambda p: -Jt(p), Jt, mesh=mesh)[0] == 0
    assert relative_class(lambda p: S @ Jt(p @ S) @ S, Jt, mesh=mesh)[0] == 1
    R = torch.diag(torch.tensor([1., 1, 1, 1, 1, -1], dtype=torch.float64))       # theta -> -theta reverses
    with pytest.raises(ValueError):
        relative_class(lambda p: R @ Jt(p @ R) @ R, Jt, mesh=mesh)
