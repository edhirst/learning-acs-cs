"""
Nijenhuis tensor in a global (non-coordinate) frame — the chart-free
formulation for parallelizable manifolds.

For a frame {E_a} with [E_a, E_b] = c^e_{ab} E_e and an ACS with frame
components J(E_a) = J^e_a E_e, the Nijenhuis tensor
N(X,Y) = [JX,JY] − J[JX,Y] − J[X,JY] − [X,Y] has components

  N^e_{ab} =  J^c_a J^d_b c^e_{cd}
            + J^c_a E_c(J^e_b) − J^d_b E_d(J^e_a)
            − J^e_f ( J^c_a c^f_{cb} − E_b(J^f_a) )
            − J^e_f ( J^d_b c^f_{ad} + E_a(J^f_b) )
            − c^e_{ab}.

Setting c = 0 recovers the coordinate-frame formula used by
src/loss/sphere_loss._nij_tensor (unit-tested). Setting the derivative terms
to zero gives the algebraic Nijenhuis condition on a Lie algebra — the
integrability condition for frame-CONSTANT structures (Samelson).

Layouts:
  J        : (B, d, d)        J[b, e, a] = J^e_a
  dJ_frame : (B, d, d, d)     dJ[b, e, a, c] = E_c(J^e_a)   (derivative last)
  c        : (B, d, d, d)     c[b, e, a, b_] = c^e_{a b_}
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import torch


def nijenhuis_frame_tensor(J: torch.Tensor,
                           dJ_frame: torch.Tensor,
                           c: torch.Tensor) -> torch.Tensor:
    """Assemble N^e_{ab}, shape (B, d, d, d) with N[b, e, a, b_]."""
    # Derivative terms (c = 0 part) — identical structure to the coordinate
    # kernel _nij_tensor: (T − T^t) − (S − S^t), with frame-directional
    # derivatives in place of coordinate partials.
    T = torch.einsum('bca,beBc->beaB', J, dJ_frame)
    S = torch.einsum('bef,bfBa->beaB', J, dJ_frame)
    N = (T - T.transpose(-2, -1)) - (S - S.transpose(-2, -1))
    # Algebraic terms.
    N = N + torch.einsum('bca,bdB,becd->beaB', J, J, c)
    N = N - torch.einsum('bef,bca,bfcB->beaB', J, J, c)
    N = N - torch.einsum('bef,bdB,bfad->beaB', J, J, c)
    N = N - c
    return N


def frame_nij_norm(J: torch.Tensor, dJ_frame: torch.Tensor,
                   c: torch.Tensor) -> torch.Tensor:
    """Per-sample Frobenius norm of the frame-components Nijenhuis tensor.

    For (near-)orthonormal frames this is (near) the metric-natural norm —
    no chart-distortion inflation."""
    N = nijenhuis_frame_tensor(J, dJ_frame, c)
    return torch.sqrt((N ** 2).sum(dim=(1, 2, 3)) + 1e-20)
