"""
Registry of the cases of interest for the uniform (chart-free) campaign, and
the explicit classical almost complex structures used as CS-search BASE POINTS
on the non-parallelizable cases (never as targets: they are ACS, not CS,
except on S² and S²×S² where the classical structure happens to be
integrable and the search starts from a random non-integrable deformation
of it anyway).

    case      χ    formulation   ambient embedding          base ACS
    s1xs1     0    frame         ℝ²×ℝ²                      canonical J₀ (frame)
    s3xs1     0    frame         ℝ⁴×ℝ²                      canonical J₀ (frame)
    s3xs3     0    frame         ℝ⁴×ℝ⁴                      canonical J₀ (frame)
    s5xs1     0    frame         ℝ⁶×ℝ²                      canonical J₀ (frame)
    s2        2    ambient       ℝ³                         J v = x × v  (integrable)
    s4        2    ambient       ℝ⁵                         NONE (S⁴ has no ACS)
    s6        2    ambient       ℝ⁷                         octonionic (G₂)
    s2xs2     4    ambient       ℝ³×ℝ³                      product of x × ·
    s4xs2     4    ambient       ℝ⁵×ℝ³                      octonionic-product
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import torch

from .ambient import AmbientManifold
from .octonion_acs import cross7_matrix
from .octonion_product_acs import octonion_product_acs


FRAME_CASES = {
    's1xs1': (1, 1),
    's3xs1': (3, 1),
    's3xs3': (3, 3),
    's5xs1': (5, 1),
}

AMBIENT_CASES = {
    's2':    (2,),
    's4':    (4,),
    's6':    (6,),
    's2xs2': (2, 2),
    's4xs2': (4, 2),
}

ALL_CASES = list(FRAME_CASES) + list(AMBIENT_CASES)


def factor_dims(case: str):
    return FRAME_CASES.get(case) or AMBIENT_CASES[case]


def is_frame_case(case: str) -> bool:
    return case in FRAME_CASES


def ambient_manifold(case: str) -> AmbientManifold:
    return AmbientManifold(AMBIENT_CASES[case])


def _cross3_matrix(x: torch.Tensor) -> torch.Tensor:
    """[x]_× with [x]_× v = x × v, shape (B, 3) → (B, 3, 3)."""
    z = torch.zeros_like(x[:, 0])
    return torch.stack([
        torch.stack([z, -x[:, 2], x[:, 1]], dim=1),
        torch.stack([x[:, 2], z, -x[:, 0]], dim=1),
        torch.stack([-x[:, 1], x[:, 0], z], dim=1),
    ], dim=1)


def base_acs(case: str, p: torch.Tensor, projector) -> torch.Tensor:
    """Explicit tangential ACS at ambient points p for an ambient case.
    Raises for 's4' (no almost complex structure exists on S⁴)."""
    P = projector(p)
    if case == 's2':
        return P @ _cross3_matrix(p) @ P
    if case == 's2xs2':
        A1 = _cross3_matrix(p[:, :3])
        A2 = _cross3_matrix(p[:, 3:])
        top = torch.cat([A1, torch.zeros_like(A1)], dim=2)
        bot = torch.cat([torch.zeros_like(A2), A2], dim=2)
        return P @ torch.cat([top, bot], dim=1) @ P
    if case == 's6':
        return P @ cross7_matrix(p) @ P
    if case == 's4xs2':
        return octonion_product_acs(p, projector)
    if case == 's4':
        raise ValueError("S^4 admits no almost complex structure: no base point exists")
    raise KeyError(case)
