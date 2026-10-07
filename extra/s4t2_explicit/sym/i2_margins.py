"""Numerical margins for the fixed-point condition I2 (docs/s4t2_symbolic_programme.md, section 4.3): on the
linear-blend shell EPS1 <= rad <= EPS2 of each cusp, measure
   C_L  = max |Lambda(vq) - M vq| / |vq|^2          (Lambda = the unblended P_final o moebius o R o h_B, in vq),
   C_L' = max |D Lambda(vq) - M| / |vq|             (Jacobian in vq by autograd through the chart point),
and the resulting lower bound on the smallest singular value of D Lambda_mod from the identity N4:
   sigma_min >= scale - c |D Lambda - M| - |grad c| |Lambda - M vq|,   |grad c| <= 2/(EPS2 - EPS1) = 500.
Also the actual sigma_min of D Lambda_mod (autograd) on the shell.  Sanity for the explicit-constant certificate.
Run from the repo root: python3 extra/s4t2_explicit/sym/i2_margins.py
"""
import math, os, sys
import torch
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
torch.set_default_dtype(torch.float64)
from extra.s4t2_explicit import tcore as T, tlift as L
from torch.func import jvp

EYE4 = torch.eye(4)


def jac(f, p):
    """Batched Jacobian (n, out, 4) of f: R^4 -> R^out by forward-mode jvp along the basis (complex-safe)."""
    return torch.stack([jvp(f, (p,), (EYE4[k].expand_as(p),))[1] for k in range(4)], -1)


torch.manual_seed(5)
n = 4000
for name, C, hb, lam_chart, M, scale, invr in (('cusp 0', T.CUSP0, L.hB_chart0, L.Lam_chart0, L.M0, L.S0, False),
                                               ('cusp inf', T.CUSPINF, L.hB_chartinf, L.Lam_chartinf, L.MINF, L.SINF, True)):
    rad = torch.empty(n).uniform_(L.EPS1, L.EPS2)
    frac = torch.empty(n).uniform_(0, 1)
    a, b = torch.empty(n).uniform_(-math.pi, math.pi), torch.empty(n).uniform_(-math.pi, math.pi)
    x = torch.polar(torch.sqrt(rad * frac / C.a0), a)
    y = torch.polar(torch.sqrt(rad * (1 - frac)), b)

    def lam_unblended(x, y):                       # v of Lambda (no blend) as a function of the chart point
        r = (x * y).abs()
        v, X = L.Lam_from_u(hb(x, y), 1 / r if invr else r)
        return v

    def vq_of(x, y):
        return L.v_q(C, x, y)

    def real4(x, y):
        return torch.stack([x.real, x.imag, y.real, y.imag], -1)

    def from_real4(p):
        return torch.complex(p[..., 0], p[..., 1]), torch.complex(p[..., 2], p[..., 3])

    p = real4(x, y)
    v = lam_unblended(x, y)
    vq = vq_of(x, y)
    err = (v - vq @ M.T).norm(dim=-1)
    nq = vq.norm(dim=-1)
    CL = (err / nq ** 2).max().item()
    # Jacobians in vq: D_p v (3x4) = D_vq v (3x3) . D_p vq (3x4); D_p vq has rank 3 off the origin, so D_vq v = D_p v . pinv(D_p vq)
    Jv = jac(lambda pp: lam_unblended(*from_real4(pp)), p)
    Jq = jac(lambda pp: vq_of(*from_real4(pp)), p)
    DL = Jv @ torch.linalg.pinv(Jq)
    dev = torch.linalg.matrix_norm(DL - M, ord=2)
    CLp = (dev / nq).max().item()
    # the blended map and its smallest singular value in vq
    Jm = jac(lambda pp: lam_chart(*from_real4(pp))[0], p)
    DLm = Jm @ torch.linalg.pinv(Jq)
    smin = torch.linalg.svdvals(DLm)[:, -1]
    bound = scale - dev.max().item() - 500 * err.max().item()
    print(f'{name}: shell rad in [{L.EPS1}, {L.EPS2}], {n} points')
    print(f'   max |Lambda - M vq| / |vq|^2 = {CL:.3f};  max |D Lambda - M|_2 / |vq| = {CLp:.3f};  scale = {scale:.4f}')
    print(f'   max |Lambda - M vq| = {err.max():.2e}, max |D Lambda - M|_2 = {dev.max():.3e}')
    print(f'   certificate bound sigma_min(D Lambda_mod) >= {bound:.4f};  measured min sigma_min = {smin.min():.4f}')
