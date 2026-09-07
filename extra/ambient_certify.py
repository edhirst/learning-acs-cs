"""A-posteriori certification of ambient-PINN structures (red-team REC 6).

For a converged run this script:
  1. evaluates SUP norms on a dense sample (MC mean-loss gives L² control only);
  2. for free-K runs (s4s2_*), applies the POLAR CORRECTION
         S = −K² + (I−P)   (SPD iff ‖K²+P‖_∞ < 1 on M;  [S, K] = 0)
         J* = K · S^{−1/2}   ⇒   J*² = −P EXACTLY,
     turning an approximate ACS into an exact one, then recomputes the
     Nijenhuis norm OF THE CORRECTED structure — the honest number;
  3. for s6_cs runs J² = −P holds by construction; it reports sup ‖N‖ and the
     octonionic distance.

Run:  python extra/ambient_certify.py outputs/ambient/<run>.pt [--n 20000]
"""
from __future__ import annotations

import os
import sys
import argparse
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.geometry.ambient import S4xS2, S6
from src.geometry.octonion_acs import cross7_matrix
from src.loss.ambient_loss import (ambient_nij_norm, compute_J_and_DJ,
                                   compute_P_and_DP)
from extra.ambient_pinn_train import AmbientNet, make_J_fn

DTYPE = torch.float64
CHUNK = 1000


def polar_correct_fn(J_fn, manifold):
    """Wrap a free-K field with the exact polar correction J* = K S^{−1/2}."""
    def J_star(x):
        K = J_fn(x)
        P = manifold.projector(x)
        eye = torch.eye(manifold.ambient_dim, dtype=x.dtype).expand_as(P)
        S = -K @ K + (eye - P)
        w, V = torch.linalg.eigh(S)                      # SPD if converged
        S_inv_sqrt = V @ torch.diag_embed(w.clamp(min=1e-12) ** -0.5) @ V.transpose(1, 2)
        return K @ S_inv_sqrt
    return J_star


def dense_eval(J_fn, manifold, n_pts):
    accs, nijs = [], []
    for s in range(0, n_pts, CHUNK):
        x = manifold.sample(min(CHUNK, n_pts - s))
        J, DJ = compute_J_and_DJ(J_fn, x)
        P, DP = compute_P_and_DP(manifold, x)
        with torch.no_grad():
            accs.append(torch.sqrt(((J @ J + P) ** 2).sum(dim=(1, 2)) + 1e-20))
        nijs.append(ambient_nij_norm(J.detach(), DJ.detach(), P, DP))
    acs = torch.cat(accs)
    nij = torch.cat(nijs)
    return acs, nij


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('checkpoint')
    ap.add_argument('--n', type=int, default=20000)
    args = ap.parse_args()
    torch.manual_seed(123)

    ck = torch.load(args.checkpoint, map_location='cpu', weights_only=False)
    case, sign = ck['case'], ck.get('sign', 1.0)
    stages = ck.get('stages', 1)
    manifold = S4xS2 if case.startswith('s4s2') else S6
    n = manifold.ambient_dim
    model = AmbientNet(n, n, stages=stages, hidden=ck.get('hidden')).to(DTYPE)
    model.load_state_dict(ck['model_state_dict'])
    model.eval()
    J_fn = make_J_fn(case, model, manifold, sign=sign)

    print(f"CERTIFY {args.checkpoint}  [{case}, seed {ck.get('seed')}, "
          f"stages {stages}, {args.n} pts]")
    acs, nij = dense_eval(J_fn, manifold, args.n)
    print(f"  raw field:      ||J^2+P||  mean {acs.mean():.4e}  SUP {acs.max():.4e}")
    print(f"                  ||N||      mean {nij.mean():.4e}  SUP {nij.max():.4e}")

    if case in ('s4s2_acs', 's4s2_cs'):
        # Polar correction applies to FREE-K fields only. Conjugation-based
        # fields (s4s2_cs_conj, s6_cs) are exact ACS by construction; applying
        # the eigh-based correction to S = I would give NaN autograd (repeated
        # eigenvalues) in the corrected-field Nijenhuis evaluation.
        if acs.max() >= 1.0:
            print(f"  polar correction NOT admissible (sup ||K^2+P|| = {acs.max():.3f} >= 1);")
            print(f"  no exact structure certified from this run.")
            return
        J_star = polar_correct_fn(J_fn, manifold)
        acs2, nij2 = dense_eval(J_star, manifold, args.n)
        print(f"  CORRECTED (exact ACS by construction; residual = numerics):")
        print(f"                  ||J*^2+P|| mean {acs2.mean():.4e}  SUP {acs2.max():.4e}")
        print(f"                  ||N(J*)||  mean {nij2.mean():.4e}  SUP {nij2.max():.4e}")
        if acs2.max() < 1e-8:
            print(f"  *** EXACT ACS on S4xS2 CERTIFIED (sup residual {acs2.max():.1e}) ***")
    elif case == 's6_cs':
        x = manifold.sample(4000)
        with torch.no_grad():
            J = J_fn(x)
            P = manifold.projector(x)
            Joct = sign * (P @ cross7_matrix(x) @ P)
            d = (J - Joct).norm(dim=(1, 2)) / Joct.norm(dim=(1, 2))
        print(f"  octonionic distance: mean {d.mean():.3f}   "
              f"[octonionic ||N|| = 19.5959; integrable would be ~0]")
    else:  # s4s2_cs_conj: exact ACS by construction; distance to the base point
        from src.geometry.octonion_product_acs import octonion_product_acs
        x = manifold.sample(4000)
        with torch.no_grad():
            J = J_fn(x)
            Jref = octonion_product_acs(x, manifold.projector)
            d = (J - Jref).norm(dim=(1, 2)) / Jref.norm(dim=(1, 2))
        print(f"  distance from octonionic-product ACS base point: mean {d.mean():.3f}   "
              f"[base ||N|| = 15.17; integrable would be ~0]")


if __name__ == '__main__':
    main()
