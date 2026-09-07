"""Evaluate an ambient-PINN field under the paper's ORIGINAL chart protocol.

Converts the global ambient field J_amb(p) (optionally polar-corrected to an
exact ACS) to ball-chart components on every patch,

    J_chart(x) = pinv(Dσ(x)) · J_amb(σ(x)) · Dσ(x),

then evaluates exactly what the results table reports: per-patch hard-cutoff
Nijenhuis (CS column) and the cross-patch overlap (ACS column). Because the
ambient field is one global function, the overlap term vanishes identically
for any tangential field; the meaningful ACS quality is the exact-ACS
residual, which the polar correction drives to machine precision.

Run:  python extra/ambient_to_chart_eval.py outputs/ambient/<run>.pt [--polar]
"""
from __future__ import annotations

import os
import sys
import argparse
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.geometry.sphere import EvenSphere, RADIAL_MIDPOINT
from src.geometry.sphere_product import SphereProduct
from src.geometry.ambient import S4xS2, S6
from src.loss.sphere_loss import _nij_tensor, OverlapLoss
from src.loss.product_loss import ProductOverlapLoss
from extra.ambient_pinn_train import AmbientNet, make_J_fn
from extra.ambient_certify import polar_correct_fn
from extra.frame_to_chart_eval import factor_embedding, dJ_autograd

DTYPE = torch.float64
OVERLAP_SIZE = 0.2
R_CUT = RADIAL_MIDPOINT + OVERLAP_SIZE
N_TEST = 4000
CHUNK = 500


def chart_from_ambient(J_amb_fn, dims, charts):
    """dims: factor manifold dims (e.g. (4,2) or (6,)); charts: per-factor
    chart indices. Returns f(x_chart) -> J_chart, differentiable."""
    def emb_single(xc):
        parts, off = [], 0
        for m, ch in zip(dims, charts):
            parts.append(factor_embedding(xc[off:off + m].unsqueeze(0), ch).squeeze(0))
            off += m
        return torch.cat(parts)

    def J_of(x_chart):
        p = torch.func.vmap(emb_single)(x_chart)                      # (B, N)
        D = torch.func.vmap(torch.func.jacfwd(emb_single))(x_chart)   # (B, N, d)
        Ja = J_amb_fn(p)                                              # (B, N, N)
        DtD = D.transpose(1, 2) @ D
        D_pinv = torch.linalg.inv(DtD) @ D.transpose(1, 2)
        return D_pinv @ Ja @ D
    return J_of


def nij_hard_cutoff(J_of, x_chart, dims):
    """Paper-protocol per-patch Nijenhuis: hard cutoff per factor, mean norm."""
    off, mask = 0, torch.ones(len(x_chart), dtype=torch.bool)
    for m in dims:
        mask &= x_chart[:, off:off + m].norm(dim=1) <= R_CUT
        off += m
    xa = x_chart[mask]
    d = sum(dims)
    norms = []
    for s in range(0, len(xa), CHUNK):
        J, dJ = dJ_autograd(J_of, xa[s:s + CHUNK], d)
        N = _nij_tensor(J, dJ)
        norms.append(torch.sqrt((N ** 2).sum(dim=(1, 2, 3)) + 1e-20))
    v = torch.cat(norms)
    return v.mean().item(), v.max().item(), int(mask.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('checkpoint')
    ap.add_argument('--polar', action='store_true', help='apply exact polar correction first (free-K runs)')
    args = ap.parse_args()
    torch.manual_seed(5)

    ck = torch.load(args.checkpoint, map_location='cpu', weights_only=False)
    case = ck['case']
    manifold = S4xS2 if case.startswith('s4s2') else S6
    n = manifold.ambient_dim
    model = AmbientNet(n, n, stages=ck.get('stages', 1), hidden=ck.get('hidden')).to(DTYPE)
    model.load_state_dict(ck['model_state_dict'])
    model.eval()
    J_amb = make_J_fn(case, model, manifold, sign=ck.get('sign', 1.0))
    if args.polar:
        J_amb = polar_correct_fn(J_amb, manifold)

    print(f"CHART-PROTOCOL EVAL of {args.checkpoint} [{case}{' + polar' if args.polar else ''}]")
    if case.startswith('s4s2'):
        dims = (4, 2)
        sp = SphereProduct(4, 2)
        xa0, xb0 = sp.sample_points(N_TEST)
        xa0, xb0 = xa0.to(DTYPE), xb0.to(DTYPE)
        pa = {0: xa0, 1: EvenSphere.patch_change_coords(xa0)}
        pb = {0: xb0, 1: EvenSphere.patch_change_coords(xb0)}
        nijs, J_dict = {}, {}
        for ia, ib in [(0, 0), (0, 1), (1, 0), (1, 1)]:
            J_of = chart_from_ambient(J_amb, dims, (ia, ib))
            xc = torch.cat([pa[ia], pb[ib]], dim=1)
            mean, mx, npts = nij_hard_cutoff(J_of, xc, dims)
            nijs[(ia, ib)] = mean
            print(f"  Nij patch ({ia},{ib}): mean {mean:.4e}  max {mx:.3e}  [{npts} pts]")
            with torch.no_grad():
                J_dict[(ia, ib)] = J_of(xc)
        ov = ProductOverlapLoss(4, 2, OVERLAP_SIZE, OVERLAP_SIZE).to(DTYPE)
        with torch.no_grad():
            ov_loss, _ = ov(J_dict, xa0, xb0, pa[1], pb[1], hard_cutoff=True)
        print(f"  MEAN NIJ (CS column):   {sum(nijs.values())/4:.4e}   [table: 6.28]")
        print(f"  OVERLAP (ACS column):   {ov_loss.item():.4e}   [table: 1.91]")
    else:
        dims = (6,)
        sphere = EvenSphere(dimension=6)
        x1 = sphere.sample_points_patch1(N_TEST).to(DTYPE)
        x2 = EvenSphere.patch_change_coords(x1)
        res = {}
        for patch, xc in [(1, x1), (2, x2)]:
            J_of = chart_from_ambient(J_amb, dims, (patch - 1,))
            mean, mx, npts = nij_hard_cutoff(J_of, xc, dims)
            res[patch] = (J_of, mean)
            print(f"  Nij patch {patch}: mean {mean:.4e}  max {mx:.3e}  [{npts} pts]")
        with torch.no_grad():
            J1 = res[1][0](x1)
            J2 = res[2][0](x2)
            ov_loss, _ = OverlapLoss(overlap_size=OVERLAP_SIZE).to(DTYPE)(J1, J2, x1, x2, hard_cutoff=True)
        print(f"  MEAN NIJ (CS column):   {(res[1][1] + res[2][1]) / 2:.4e}   [table: 9.87; octonionic flat-chart ref 54.0]")
        print(f"  OVERLAP (ACS column):   {ov_loss.item():.4e}   [table: 3.57e-2]")


if __name__ == '__main__':
    main()
