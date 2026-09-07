"""
Nijenhuis tensor of a tangential endomorphism field in the ambient-projector
representation (chart-free; see src/geometry/ambient.py).

For M ⊂ ℝ^n with projector field P(x), take the tangent vector fields
X_i(x) = P(x)e_i (i = 1..n, overcomplete) with smooth ambient extensions, and
a tangential field J = P J̃ P (so J e_i is automatically an ambient extension
of J X_i, since JP = J). With ambient brackets [U,V] = (DV)·U − (DU)·V —
extension-independent on M for fields tangent along M — the (general-endo)
Nijenhuis tensor is

  N_{ij} = [Je_i, Je_j] − J[Je_i, Pe_j] − J[Pe_i, Je_j] + J·J·[Pe_i, Pe_j].

The last term uses J² literally (for exact ACS parametrizations J² = −id_T and
it equals −[Pe_i, Pe_j], which is tangent on M; for soft-constrained fields it
keeps N exact during constraint violation).

Because Σ_i (Pe_i)(Pe_i)ᵀ = P, the per-sample quantity Σ_{a,i,j} (N^a_{ij})²
equals the intrinsic ROUND-METRIC norm² of the induced Nijenhuis tensor (the
induced metric on M is the restricted euclidean one). Validation: the
octonionic ACS on S⁶ must give the G₂-invariant constant ‖N‖ = 19.5959…,
matching the independent chart-coordinate computation.

Layouts: J, P: (B, n, n); DJ, DP: (B, n, n, n) with DJ[b, a, r, c] = ∂_c J_ar.
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import torch


def ambient_brackets(J, DJ, P, DP):
    """The three bracket families, each (B, n, n, n) as [b, a, i, j]."""
    # [Je_i, Je_j]^a = DJ[a,j,c] J[c,i] − DJ[a,i,c] J[c,j]
    JJ = torch.einsum('bajc,bci->baij', DJ, J)
    B1 = JJ - JJ.transpose(2, 3)
    # [Je_i, Pe_j]^a = DP[a,j,c] J[c,i] − DJ[a,i,c] P[c,j]
    B2 = torch.einsum('bajc,bci->baij', DP, J) - \
         torch.einsum('baic,bcj->baij', DJ, P)
    # [Pe_i, Je_j]^a = DJ[a,j,c] P[c,i] − DP[a,i,c] J[c,j]  (= −B2 with i↔j)
    B3 = -B2.transpose(2, 3)
    # [Pe_i, Pe_j]^a = DP[a,j,c] P[c,i] − DP[a,i,c] P[c,j]
    PP = torch.einsum('bajc,bci->baij', DP, P)
    B4 = PP - PP.transpose(2, 3)
    return B1, B2, B3, B4


def ambient_nijenhuis(J, DJ, P, DP):
    """N[b, a, i, j] = N(X_i, X_j)^a — the induced Nijenhuis tensor evaluated
    on the overcomplete tangent system X_i = Pe_i."""
    B1, B2, B3, B4 = ambient_brackets(J, DJ, P, DP)
    Jb = lambda T: torch.einsum('bae,beij->baij', J, T)
    return B1 - Jb(B2) - Jb(B3) + Jb(Jb(B4))


def ambient_nij_norm(J, DJ, P, DP):
    """Per-sample intrinsic (round-metric) Nijenhuis norm, shape (B,)."""
    N = ambient_nijenhuis(J, DJ, P, DP)
    return torch.sqrt((N ** 2).sum(dim=(1, 2, 3)) + 1e-20)


def compute_J_and_DJ(J_fn, x):
    """J (B,n,n) and its ambient Jacobian DJ (B,n,n,n) via forward AD."""
    def single(xs):
        return J_fn(xs.unsqueeze(0)).squeeze(0)
    J = J_fn(x)
    DJ = torch.func.vmap(torch.func.jacfwd(single))(x)
    return J, DJ


def compute_P_and_DP(manifold, x):
    def single(xs):
        return manifold.projector(xs.unsqueeze(0)).squeeze(0)
    P = manifold.projector(x)
    DP = torch.func.vmap(torch.func.jacfwd(single))(x)
    return P, DP
