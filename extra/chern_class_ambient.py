"""First Chern class of an ambient ACS field evaluated on the S² cycle.

For a tangential ACS J on M ⊂ ℝ^N (J P = P J = J, J² = −P) the (1,0)-bundle
is the image of the complex projector Π = (P + iJ)/2 ⊂ ℂ^N. Its projected
connection has curvature Π dΠ ∧ dΠ, so on a 2-cycle Σ parametrised by (θ, φ):

    ⟨c₁(T^{1,0}M), [Σ]⟩ = (i/2π) ∫_Σ tr(Π [∂_θΠ, ∂_φΠ]) dθ dφ.

c₁ is a homotopy invariant of the ACS, hence constant along the exact-ACS
conjugation flow J = e^B J_base e^{−B}: the CS search can only reach complex
structures in the base point's Chern class. Validation: S² (c₁ = ±2), S²×S²
(c₁[S²×pt] = ±2). Usage: python extra/chern_class_ambient.py
"""
from __future__ import annotations
import os, sys, torch
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from src.geometry.explicit_acs import ambient_manifold, base_acs

torch.set_default_dtype(torch.float64)


def c1_on_s2_cycle(J_fn, embed, n=200):
    """embed(theta, phi) -> ambient point (N,) parametrising the 2-cycle."""
    th = (torch.arange(n) + 0.5) / n * torch.pi
    ph = (torch.arange(2 * n) + 0.5) / (2 * n) * 2 * torch.pi
    TH, PH = torch.meshgrid(th, ph, indexing='ij')
    q = torch.stack([TH.reshape(-1), PH.reshape(-1)], dim=1)

    def PJ_single(t):                     # real-valued: stack(P, J) with P = -J²
        p = embed(t[0], t[1]).unsqueeze(0)
        J = J_fn(p)[0]
        return torch.stack([-J @ J, J])

    PJ = torch.func.vmap(PJ_single)(q)                                   # (B,2,N,N)
    dPJ = torch.func.vmap(torch.func.jacfwd(PJ_single))(q)               # (B,2,N,N,2)
    Pi = torch.complex(PJ[:, 0], PJ[:, 1]) / 2.0
    dPi = torch.complex(dPJ[:, 0], dPJ[:, 1]) / 2.0                      # (B,N,N,2)
    comm = dPi[..., 0] @ dPi[..., 1] - dPi[..., 1] @ dPi[..., 0]
    integrand = (1j * torch.einsum('bij,bji->b', Pi, comm)).real
    cell = (torch.pi / n) * (2 * torch.pi / (2 * n))
    return integrand.sum().item() * cell / (2 * torch.pi)


def sph2(t, p):
    return torch.stack([torch.sin(t) * torch.cos(p), torch.sin(t) * torch.sin(p), torch.cos(t)])


if __name__ == '__main__':
    M2 = ambient_manifold('s2')
    J2 = lambda p: base_acs('s2', p, M2.projector)
    print(f"S²  cross-product ACS:          c1[S²]        = {c1_on_s2_cycle(J2, sph2):+.4f}   (expect ±2)")

    M22 = ambient_manifold('s2xs2')
    J22 = lambda p: base_acs('s2xs2', p, M22.projector)
    y0 = torch.tensor([0.3, -0.5, 0.8]); y0 = y0 / y0.norm()
    print(f"S²×S² product ACS:              c1[S²×pt]     = {c1_on_s2_cycle(J22, lambda t,p: torch.cat([sph2(t,p), y0])):+.4f}   (expect ±2)")
    print(f"S²×S² product ACS:              c1[pt×S²]     = {c1_on_s2_cycle(J22, lambda t,p: torch.cat([y0, sph2(t,p)])):+.4f}   (expect ±2)")

    M42 = ambient_manifold('s4xs2')
    J42 = lambda p: base_acs('s4xs2', p, M42.projector)
    torch.manual_seed(0)
    for k in range(3):
        x0 = torch.randn(5); x0 = x0 / x0.norm()
        c = c1_on_s2_cycle(J42, lambda t, p: torch.cat([x0, sph2(t, p)]))
        print(f"S⁴×S² octonionic-product ACS:   c1[x0×S²]     = {c:+.4f}   (x0 sample {k}; Liang–Liu–Zhang CS has c1 = 0)")
    y0 = torch.tensor([0.0, 0.0, 1.0])
    c = c1_on_s2_cycle(J42, lambda t, p: torch.cat([sph2(t, p), torch.zeros(2), y0]))
    print(f"S⁴×S² octonionic-product ACS:   c1[S²⊂S⁴×pt] = {c:+.4f}   (null-homologous: expect 0)")
