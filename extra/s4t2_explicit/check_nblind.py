"""The Nijenhuis check of a transported structure is blind to the transporting map (section 12.12): at sample points q,
replace the 1- and 2-jet of F by those of a RANDOM other local diffeomorphism into S^4 x T^2 (1-jet A M, M random;
tangential second derivatives random symmetric, normal part fixed by the second fundamental form) and compute N with
the same code: N = 0 again. Only non-symmetric 'second derivatives' (not the 2-jet of any map) give N != 0.
Hence N(F_* J_Y) = 0 is the naturality of N, valid for every F; it carries no information about F."""
import os, sys
import torch
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
torch.set_default_dtype(torch.float64)
from extra.s4t2_explicit import nijenhuis_round as NR, check_full as CF, tfull as Fm, tcore as T

g = torch.Generator().manual_seed(1)
x, y, w = CF.chart_sample(400, T.CUSP0, rmax=0.25, rmin=1e-3)
q, A, H = NR.DF_and_D2F(Fm.F_c0_real, CF.real6(x, y, w))
B = len(q)
P = NR.S4T2.projector(q)


def report(nm, A2, H2):
    J, DJ = NR.transport(A2, H2, q)
    Nn, scale, _ = NR.nij(q, J, DJ)
    print(f'{nm:<80s} |N|/scale max {(Nn / scale).max():.1e}  (|J| median {J.flatten(1).norm(dim=1).median():.1f})')


def normal_part(Ab):
    out = torch.zeros(B, 9, 6, 6)
    for a, b in ((0, 5), (5, 7), (7, 9)):
        G = torch.einsum('bik,bil->bkl', Ab[:, a:b], Ab[:, a:b])
        out[:, a:b] = -q[:, a:b, None, None] * G[:, None]
    return out


report('our F (true 1- and 2-jet)', A, H)
AM = A @ torch.randn(B, 6, 6, generator=g)
S = torch.randn(B, 9, 6, 6, generator=g)
S = S + S.transpose(2, 3)
Hgen = torch.einsum('bij,bjkl->bikl', P, S) + normal_part(AM)
report('a random other local diffeomorphism at q (random 1-jet, random symmetric 2-jet)', AM, Hgen)
Hns = Hgen + 0.3 * torch.einsum('bij,bjkl->bikl', P, torch.randn(B, 9, 6, 6, generator=g))
report('the same with NON-symmetric second derivatives (not the 2-jet of any map)', AM, Hns)
