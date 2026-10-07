"""The complex structure of Y transported to the ROUND S^4 x T^2 by the explicit diffeomorphism F, and its Nijenhuis
tensor computed with the repository's ambient code (src/loss/ambient_loss.py), docs/s4t2_explicit_diffeo.md 12.7.

For a sample p in a holomorphic chart of Y (coordinates R^6, J_Y = J0 standard there):
    A = DF(p) (9 x 6),  J = A J0 A^+  with A^+ = (A^T A)^-1 A^T   (tangential, J^2 = -P on T(S^4 x T^2)),
    DJ[a, r, c] = sum_k (d_k J)[a, r] A^+[k, c],   d_k J = (d_k A) J0 A^+ + A J0 d_k(A^+),
    d_k(A^+) = -A^+ (d_k A) A^+ + (A^T A)^-1 (d_k A)^T (I - A A^+),
which is the ambient derivative of the extension J(q') = J(p + A^+ (q' - F(p))) (exact in tangent directions; the
ambient brackets only use those). D F and D^2 F by nested forward-mode AD (float64). P, DP: AmbientManifold((4,1,1)).
N = ambient_nijenhuis(J, DJ, P, DP); the per-sample norm is the intrinsic round-metric norm (Parseval frame).

Controls: (i) J from two charts at the same point of Y agree (holomorphic gluings, F well defined);
          (ii) a non-integrable almost complex structure J_eps = M J0 M^-1 (M = I + eps sin(p_0) E_{1,3}) pushed through
               the same F gives |N| = O(eps): the pipeline detects non-integrability.
Run from the repo root: python3 extra/s4t2_explicit/nijenhuis_round.py [n_per_region]
"""
from __future__ import annotations
import math, os, sys, time
import torch
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
torch.set_default_dtype(torch.float64)
from torch.func import jvp
from src.geometry.ambient import AmbientManifold
from src.loss.ambient_loss import ambient_brackets, ambient_nijenhuis, compute_P_and_DP
from extra.s4t2_explicit import tcore as T, tfull as Fm, check_full as CF

S4T2 = AmbientManifold((4, 1, 1))
J0 = torch.zeros(6, 6)
for i in range(3):
    J0[2 * i + 1, 2 * i], J0[2 * i, 2 * i + 1] = 1.0, -1.0
EYE6 = torch.eye(6)


def DF_and_D2F(f, p):
    def DF(pp):
        return torch.stack([jvp(f, (pp,), (EYE6[k].expand_as(pp),))[1] for k in range(6)], -1)
    A = DF(p)
    H = torch.stack([jvp(DF, (p,), (EYE6[l].expand_as(p),))[1] for l in range(6)], -1)   # H[b, a, k, l] = d_l d_k F_a
    return f(p), A, H


def transport(A, H, q, Jc=None, dJc=None):
    """Ambient J and DJ from A = DF, H = D^2 F (at q = F(p)) and a chart complex structure Jc (default J0) with chart
    derivative dJc[b, :, :, k] = d_k Jc. A^+ by SVD (error ~ cond(A) eps, not cond^2 as with normal equations);
    I - A A^+ is the exact normal projector I - P(q) of S^4 x T^2."""
    B = A.shape[0]
    Jc = J0.expand(B, 6, 6) if Jc is None else Jc
    Ap = torch.linalg.pinv(A)                                     # (B, 6, 9)
    Gi = Ap @ Ap.transpose(1, 2)                                  # (A^T A)^-1
    J = A @ Jc @ Ap
    Pt = torch.eye(9) - S4T2.projector(q)
    dJ = []
    for k in range(6):
        dA = H[..., k]                                            # d_k A  (B, 9, 6)
        dAp = -Ap @ dA @ Ap + Gi @ dA.transpose(1, 2) @ Pt
        term = dA @ Jc @ Ap + A @ Jc @ dAp
        if dJc is not None:
            term = term + A @ dJc[..., k] @ Ap
        dJ.append(term)
    dJ = torch.stack(dJ, -1)                                      # (B, 9, 9, 6): d_k J
    DJ = torch.einsum('barK,bKc->barc', dJ, Ap)                   # ambient derivative along A^+ (tangent directions)
    return J, DJ


def nij(q, J, DJ):
    P, DP = compute_P_and_DP(S4T2, q)
    N = ambient_nijenhuis(J, DJ, P, DP)
    B1, B2, B3, B4 = ambient_brackets(J, DJ, P, DP)
    Jb = lambda X: torch.einsum('bae,beij->baij', J, X)
    scale = sum(t.flatten(1).norm(dim=-1) for t in (B1, Jb(B2), Jb(B3), Jb(Jb(B4))))
    return N.flatten(1).norm(dim=-1), scale, P


def samples(n):
    g = CF.gen
    out = {}
    z, w, t = CF.mid_sample(n, 0.25, 4.0)
    out['mid (0.25 <= |t| <= 4)'] = (Fm.F_mid_real, CF.real6(z, w, t))
    for nm, C, f in (('c0', T.CUSP0, Fm.F_c0_real), ('ci', T.CUSPINF, Fm.F_ci_real)):
        x, y, w = CF.chart_sample(n, C, rmax=0.25, rmin=1e-8)
        m = len(x) // 10                                          # 10% on the central fibre itself (t = 0 exactly)
        x[:m // 2] = 0.0
        y[m // 2:m] = 0.0
        x[m:m + m // 2] *= 1e-3                                   # and 5% close to the node
        y[m:m + m // 2] *= 1e-3
        out[f'{nm} filling (|t|, |s| <= 0.25, incl. central fibre and node)'] = (f, CF.real6(x, y, w))
    return out


def check_D2F(f, P, h=1e-6):
    """Nested-AD second derivatives vs central differences of the (exact) first derivatives, all (k, l)."""
    worst = 0.0
    for k in range(6):
        d1 = lambda pp: jvp(f, (pp,), (EYE6[k].expand_as(pp),))[1]
        for l in range(k, 6):
            ad = jvp(d1, (P,), (EYE6[l].expand_as(P),))[1]
            fd = (d1(P + h * EYE6[l]) - d1(P - h * EYE6[l])) / (2 * h)
            worst = max(worst, ((ad - fd).norm(dim=-1) / (fd.norm(dim=-1) + d1(P).norm(dim=-1) * 1e-3)).max().item())
    return worst


if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    t0 = time.time()
    for nm, (f, P) in samples(40).items():
        print(f'D^2 F check ({nm}): max relative |nested AD - FD of DF| over all (k, l) = {check_D2F(f, P):.1e}')
    allN = []
    for nm, (f, P) in samples(n).items():
        q, A, H = DF_and_D2F(f, P)
        J, DJ = transport(A, H, q)
        Nn, scale, Pq = nij(q, J, DJ)
        Jsq = (J @ J + Pq).flatten(1).norm(dim=-1).max()
        rel = Nn / scale
        print(f'{nm}: {len(P)} samples')
        print(f'   |J^2 + P| max {Jsq:.1e};  |J| in [{J.flatten(1).norm(dim=-1).min():.2f}, {J.flatten(1).norm(dim=-1).max():.2f}];  '
              f'|DJ| max {DJ.flatten(1).norm(dim=-1).max():.2e}')
        print(f'   |N| (round metric): median {Nn.median():.2e}, max {Nn.max():.2e};  |N| / (sum of bracket-term norms): max {rel.max():.2e}')
        allN.append(rel)
    # control (i): J from the c0 chart and from the mid chart at the same point
    x, y, w = CF.chart_sample(n, T.CUSP0, rmax=0.25, rmin=0.05)
    q1, A1, H1 = DF_and_D2F(Fm.F_c0_real, CF.real6(x, y, w))
    zm, wm, tm, _ = Fm.c0_to_mid(x, y, w)
    q2, A2, H2 = DF_and_D2F(Fm.F_mid_real, CF.real6(zm, wm, tm))
    Ja, DJa = transport(A1, H1, q1)
    Jb, DJb = transport(A2, H2, q2)
    print(f'control (i): same points via the c0 and the mid chart: |q1 - q2| {(q1 - q2).abs().max():.1e}, '
          f'|J1 - J2| / |J| {((Ja - Jb).flatten(1).norm(dim=-1) / Ja.flatten(1).norm(dim=-1)).max():.1e}, '
          f'|DJ1 - DJ2| / |DJ| {((DJa - DJb).flatten(1).norm(dim=-1) / DJa.flatten(1).norm(dim=-1)).max():.1e}')
    # control (ii): a non-integrable almost complex structure through the same F
    f, P = samples(n)['mid (0.25 <= |t| <= 4)']
    q, A, H = DF_and_D2F(f, P)
    E13 = torch.zeros(6, 6); E13[0, 2] = 1.0
    for eps in (1e-1, 1e-3):
        s = torch.sin(P[:, 0])
        M = torch.eye(6) + eps * s[:, None, None] * E13
        Mi = torch.eye(6) - eps * s[:, None, None] * E13                  # E13^2 = 0
        Jc = M @ J0 @ Mi
        dJc = torch.zeros(len(P), 6, 6, 6)
        dM = eps * torch.cos(P[:, 0])[:, None, None] * E13
        dJc[..., 0] = dM @ J0 @ Mi - M @ J0 @ dM                          # d(M^-1) = -dM
        J, DJ = transport(A, H, q, Jc, dJc)
        Nn, scale, _ = nij(q, J, DJ)
        print(f'control (ii): non-integrable J_eps, eps = {eps:.0e}: |N| median {Nn.median():.2e}, '
              f'|N| / scale median {(Nn / scale).median():.2e}')
    print(f'max |N| / scale over all regions (J_Y transported): {torch.cat(allN).max():.2e};  {time.time() - t0:.0f} s')
