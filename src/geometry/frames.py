"""
Global frames on parallelizable sphere products, for the frame-based
(chart-free) Nijenhuis formulation.

Motivation (from the S⁶ complex-structure construction, Alpöge–Claude/Engel
2026, and the session's diagnostics): every known complex structure on the
"known to exist" cases is simple in a suitable GLOBAL frame —

  * S³×S³ ≅ SU(2)×SU(2): left-invariant complex structures exist (Samelson);
    they are CONSTANT matrices in the bi-invariant Lie frame, where the
    integrability condition is purely algebraic.
  * S⁵×S¹ (Hopf manifold (ℂ³∖0)/⟨z ∼ 2z⟩): in the ambient trivialization
    T(S⁵×S¹) ≅ ℝ⁶ via (v, t) ↦ v + t·z, the CE structure is multiplication
    by i — a constant matrix.
  * S³×S¹: both descriptions apply (Lie group SU(2)×U(1), and Hopf S³×S¹).

In chart coordinates these same structures have ‖∂J‖ ~ 34 (measured), which
is what defeats gradient training. In the frame they are ~constant: the
derivative burden is a chart artifact. A global frame also removes patches
entirely — no overlap loss, a single objective.

Frames are represented as ambient vector fields E_a: ℝ^N → ℝ^N (tangent to
the embedded manifold). Structure functions are obtained by autograd from the
ambient bracket [E_a, E_b] = (DE_b)E_a − (DE_a)E_b, expanded in the frame —
no hand-derived constants required (they are verified constant for Lie
frames by the self-tests).

All operations are differentiable torch.
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import torch


# ---------------------------------------------------------------------------
# Quaternion helpers (S³ ⊂ ℍ ≅ ℝ⁴, coordinates (q0, q1, q2, q3))
# ---------------------------------------------------------------------------

def quat_mul(p: torch.Tensor, q: torch.Tensor) -> torch.Tensor:
    """Hamilton product, shapes (..., 4) × (..., 4) → (..., 4)."""
    p0, p1, p2, p3 = p.unbind(-1)
    q0, q1, q2, q3 = q.unbind(-1)
    return torch.stack([
        p0 * q0 - p1 * q1 - p2 * q2 - p3 * q3,
        p0 * q1 + p1 * q0 + p2 * q3 - p3 * q2,
        p0 * q2 - p1 * q3 + p2 * q0 + p3 * q1,
        p0 * q3 + p1 * q2 - p2 * q1 + p3 * q0,
    ], dim=-1)


_QUAT_UNITS = torch.eye(4, dtype=torch.float64)   # 1, i, j, k as rows 0..3


# ---------------------------------------------------------------------------
# Frame definitions (ambient vector fields)
# ---------------------------------------------------------------------------

class ProductFrame:
    """A global frame on S^{m_α} × S^{m_β} embedded in ℝ^{N_α} × ℝ^{N_β}.

    Subclasses define `frame_fields(x) -> (batch, dim, N)`: the dim frame
    vectors' ambient components at each point (rows = frame index).
    """

    dim: int          # manifold dimension (= number of frame fields)
    ambient_dim: int  # N = N_α + N_β

    def frame_fields(self, x: torch.Tensor) -> torch.Tensor:
        raise NotImplementedError

    # -- generic machinery -------------------------------------------------

    def sample(self, n: int, dtype=torch.float64) -> torch.Tensor:
        raise NotImplementedError

    def structure_functions(self, x: torch.Tensor) -> torch.Tensor:
        """c[b, c_, a, b_] = c^{c_}_{a b_} with [E_a, E_b_] = c^{c_}_{a b_} E_{c_}, by autograd:
        upper index FIRST, the layout `src/loss/frame_loss.py` expects. (Before
        2026-10-25 the upper index came last; that agreed only for totally
        antisymmetric structure constants, e.g. the Lie frames of S³×S³ and S³×S¹,
        and gave wrong frame torsion on S⁵×S¹ for non-constant J.)

        Ambient bracket [E_a, E_b] = (DE_b)·E_a − (DE_a)·E_b, where DE is the
        ambient Jacobian of the frame-field map (any smooth ambient extension
        works; the bracket of tangent fields is extension-independent).
        The result is expanded in the frame via least squares (the frame may
        not be orthonormal in ambient coordinates).
        """
        def fields_single(x_single: torch.Tensor) -> torch.Tensor:
            return self.frame_fields(x_single.unsqueeze(0)).squeeze(0)   # (dim, N)

        E = self.frame_fields(x)                                          # (B, dim, N)
        DE = torch.func.vmap(torch.func.jacfwd(fields_single))(x)         # (B, dim, N, N)
        # bracket[b, a, c, :] = DE[b,c,:,:] @ E[b,a,:] − DE[b,a,:,:] @ E[b,c,:]
        DcEa = torch.einsum('bcnk,bak->bacn', DE, E)                      # (DE_c)·E_a
        bracket = DcEa - DcEa.transpose(1, 2)                             # [E_a, E_c] rows? sign:
        # DcEa[b,a,c,:] = (DE_c)E_a ; [E_a,E_c] = (DE_c)E_a − (DE_a)E_c = DcEa[a,c] − DcEa[c,a]
        # Expand in the frame: solve E^T (dim×N)ᵀ ... use lstsq against E rows.
        B, d, N = E.shape
        Emat = E.transpose(1, 2)                                          # (B, N, dim)
        rhs = bracket.reshape(B, d * d, N).transpose(1, 2)                # (B, N, d²)
        sol = torch.linalg.lstsq(Emat, rhs).solution                      # (B, dim, d²)
        c = sol.transpose(1, 2).reshape(B, d, d, d)                       # c[b, a, b_, c_]
        return c.permute(0, 3, 1, 2).contiguous()                         # c[b, c_, a, b_]


class LieFrameS3xS3(ProductFrame):
    """Bi-invariant frame on S³×S³ ⊂ ℝ⁴×ℝ⁴: E_i(q) = q·ê_i per factor.

    Left-invariant right-multiplication fields; the structure functions are
    the su(2)⊕su(2) constants (verified constant by the self-test)."""

    dim = 6
    ambient_dim = 8

    def frame_fields(self, x: torch.Tensor) -> torch.Tensor:
        qa, qb = x[:, :4], x[:, 4:]
        units = _QUAT_UNITS.to(dtype=x.dtype, device=x.device)
        rows = []
        for i in (1, 2, 3):
            e = units[i].expand_as(qa)
            rows.append(torch.cat([quat_mul(qa, e), torch.zeros_like(qb)], dim=1))
        for i in (1, 2, 3):
            e = units[i].expand_as(qb)
            rows.append(torch.cat([torch.zeros_like(qa), quat_mul(qb, e)], dim=1))
        return torch.stack(rows, dim=1)                                   # (B, 6, 8)

    def sample(self, n: int, dtype=torch.float64) -> torch.Tensor:
        qa = torch.randn(n, 4, dtype=dtype)
        qb = torch.randn(n, 4, dtype=dtype)
        qa = qa / qa.norm(dim=1, keepdim=True)
        qb = qb / qb.norm(dim=1, keepdim=True)
        return torch.cat([qa, qb], dim=1)


class LieFrameS3xS1(ProductFrame):
    """Frame on S³×S¹ ⊂ ℝ⁴×ℝ²: su(2) frame on the S³ factor, rotation field
    on the S¹ factor. Structure functions constant."""

    dim = 4
    ambient_dim = 6

    def frame_fields(self, x: torch.Tensor) -> torch.Tensor:
        qa, w = x[:, :4], x[:, 4:]
        units = _QUAT_UNITS.to(dtype=x.dtype, device=x.device)
        rows = []
        for i in (1, 2, 3):
            e = units[i].expand_as(qa)
            rows.append(torch.cat([quat_mul(qa, e), torch.zeros_like(w)], dim=1))
        rot = torch.stack([-w[:, 1], w[:, 0]], dim=1)                     # ∂_θ on S¹
        rows.append(torch.cat([torch.zeros_like(qa), rot], dim=1))
        return torch.stack(rows, dim=1)                                   # (B, 4, 6)

    def sample(self, n: int, dtype=torch.float64) -> torch.Tensor:
        qa = torch.randn(n, 4, dtype=dtype)
        qa = qa / qa.norm(dim=1, keepdim=True)
        th = torch.rand(n, dtype=dtype) * (2 * torch.pi)
        w = torch.stack([torch.cos(th), torch.sin(th)], dim=1)
        return torch.cat([qa, w], dim=1)


class AmbientFrameS5xS1(ProductFrame):
    """Frame on S⁵×S¹ ⊂ ℝ⁶×ℝ² from the Hopf-manifold trivialization
    T(S⁵×S¹) ≅ ℝ⁶: the constant ambient basis e_k of ℝ⁶ ≅ ℂ³ splits at
    z ∈ S⁵ into a tangential part (e_k − ⟨e_k,z⟩z) and a radial coefficient
    ⟨e_k,z⟩ carried by the S¹ (log-radial) direction:

        E_k(z, w) = (e_k − ⟨e_k,z⟩ z,  ⟨e_k,z⟩ · ∂_θ(w)),   k = 1..6.

    In this frame the CE (Hopf) complex structure is multiplication by i on
    ℂ³ — a constant matrix. Structure functions are position-dependent but
    smooth and mild; obtained by autograd."""

    dim = 6
    ambient_dim = 8

    def frame_fields(self, x: torch.Tensor) -> torch.Tensor:
        z, w = x[:, :6], x[:, 6:]
        B = z.shape[0]
        eye = torch.eye(6, dtype=x.dtype, device=x.device)
        coef = z                                                           # ⟨e_k, z⟩ = z_k
        tang = eye.unsqueeze(0) - coef.unsqueeze(2) * z.unsqueeze(1)       # (B, 6, 6)
        rot = torch.stack([-w[:, 1], w[:, 0]], dim=1)                      # (B, 2)
        s1_part = coef.unsqueeze(2) * rot.unsqueeze(1)                     # (B, 6, 2)
        return torch.cat([tang, s1_part], dim=2)                           # (B, 6, 8)

    def sample(self, n: int, dtype=torch.float64) -> torch.Tensor:
        z = torch.randn(n, 6, dtype=dtype)
        z = z / z.norm(dim=1, keepdim=True)
        th = torch.rand(n, dtype=dtype) * (2 * torch.pi)
        w = torch.stack([torch.cos(th), torch.sin(th)], dim=1)
        return torch.cat([z, w], dim=1)


class LieFrameS1xS1(ProductFrame):
    """Frame on the torus S¹×S¹ ⊂ ℝ²×ℝ²: the two rotation fields. Abelian
    Lie group ⇒ structure functions vanish identically."""

    dim = 2
    ambient_dim = 4

    def frame_fields(self, x: torch.Tensor) -> torch.Tensor:
        u, w = x[:, :2], x[:, 2:]
        rot_u = torch.stack([-u[:, 1], u[:, 0]], dim=1)
        rot_w = torch.stack([-w[:, 1], w[:, 0]], dim=1)
        rows = [torch.cat([rot_u, torch.zeros_like(w)], dim=1),
                torch.cat([torch.zeros_like(u), rot_w], dim=1)]
        return torch.stack(rows, dim=1)                                   # (B, 2, 4)

    def sample(self, n: int, dtype=torch.float64) -> torch.Tensor:
        th = torch.rand(n, 2, dtype=dtype) * (2 * torch.pi)
        return torch.cat([torch.cos(th[:, :1]), torch.sin(th[:, :1]),
                          torch.cos(th[:, 1:]), torch.sin(th[:, 1:])], dim=1)


FRAMES = {
    's1xs1': LieFrameS1xS1,
    's3xs3': LieFrameS3xS3,
    's3xs1': LieFrameS3xS1,
    's5xs1': AmbientFrameS5xS1,
}
