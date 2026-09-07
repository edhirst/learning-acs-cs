"""S⁶ octonionic-ACS calibration (sphere-side analogue of the CE floor calibration).

Evaluates the paper-protocol losses AT the exact octonionic (G₂) almost
complex structure on S⁶:

  1. ACS floor — the overlap loss at the exact global ACS (expect ~1e-9):
     validates the sphere-side loss stack in d = 6.
  2. Nijenhuis reference — ‖N(J_oct)‖ under the paper protocol: the
     "octonionic basin" reference number. The structure is G₂-homogeneous, so
     the ROUND-METRIC per-sample ‖N‖ should be near-constant across S⁶ — a
     strong self-check of both the structure and the round-metric code path.
     (The flat chart-coordinate norm varies across the chart; both reported.)
  3. Basin diagnostic (optional --checkpoint): distance of a trained S⁶
     model to J_oct (value + derivative) and its Nij vs the reference. If the
     trained field sits at the octonionic reference, free search demonstrably
     converged to the round basin — metrically unrelated to the integrable
     structures of [Alpöge–Claude / Engel 2026], explaining the table's 9.87.

Diagnostics only; never used in training.

Run:  python extra/octonion_calibration.py [--checkpoint path/to/final_model.pt]
"""
from __future__ import annotations

import os
import sys
import argparse
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.geometry.sphere import EvenSphere, RADIAL_MIDPOINT
from src.geometry.octonion_acs import octonion_acs, chart_embedding
from src.loss.sphere_loss import OverlapLoss, _nij_tensor, _metric_norm_sq_12tensor

torch.manual_seed(42)
DTYPE = torch.float64
N_TEST = 10_000
OVERLAP_SIZE = 0.2      # matches configs/{ACS,CS}/config_sphere6.yaml
R_CUT = RADIAL_MIDPOINT + OVERLAP_SIZE
CHUNK = 2048


def compute_dJ(x, patch):
    """dJ[b,i,j,l] = ∂J^i_j/∂x^l via 36 reverse-mode passes (chunked by caller)."""
    xr = x.detach().clone().requires_grad_(True)
    J = octonion_acs(xr, patch=patch)
    B, d, _ = J.shape
    dJ = torch.zeros(B, d, d, d, dtype=x.dtype)
    for i in range(d):
        for j in range(d):
            g, = torch.autograd.grad(J[:, i, j].sum(), xr, retain_graph=True)
            dJ[:, i, j, :] = g
    return J.detach(), dJ


def compute_dJ_fd(x, patch, h=1e-6):
    with torch.no_grad():
        B, d = x.shape
        J = octonion_acs(x, patch=patch)
        dJ = torch.zeros(B, d, d, d, dtype=x.dtype)
        for l in range(d):
            e = torch.zeros(d, dtype=x.dtype); e[l] = h
            dJ[:, :, :, l] = (octonion_acs(x + e, patch=patch)
                              - octonion_acs(x - e, patch=patch)) / (2 * h)
    return J, dJ


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--checkpoint', default=None,
                    help='optional trained S6 checkpoint (final_model.pt) for the basin diagnostic')
    args = ap.parse_args()

    sphere = EvenSphere(dimension=6)
    x1 = sphere.sample_points_patch1(N_TEST).to(DTYPE)
    x2 = EvenSphere.patch_change_coords(x1)

    # ---------------- sanity ----------------
    print("=" * 78)
    print("Octonionic ACS on S^6 — sanity")
    print("=" * 78)
    p1 = chart_embedding(x1, 1)
    p2 = chart_embedding(x2, 2)
    print(f"  same-physical-point embedding agreement: max |σ1(x1) − σ2(φ(x1))| "
          f"= {(p1 - p2).abs().max():.2e}")
    with torch.no_grad():
        J1 = octonion_acs(x1, patch=1)
        J2 = octonion_acs(x2, patch=2)
    I = torch.eye(6, dtype=DTYPE).expand(N_TEST, 6, 6)
    for tag, J in [('patch 1', J1), ('patch 2', J2)]:
        err = (J @ J + I).norm(dim=(1, 2))
        print(f"  J^2+I {tag}: median {err.median():.2e}   max {err.max():.2e}")

    # ---------------- 1. ACS floor (overlap) ----------------
    ov = OverlapLoss(overlap_size=OVERLAP_SIZE).to(DTYPE)
    with torch.no_grad():
        ov_loss, ov_info = ov(J1, J2, x1, x2, hard_cutoff=True)
    print(f"\n  ACS FLOOR (overlap, hard cutoff): {ov_loss.item():.6e}")
    print(f"    [table: trained S^6 ACS = 3.57e-2]")

    # ---------------- 2. Nijenhuis reference ----------------
    print(f"\n  NIJENHUIS REFERENCE (hard cutoff, rho <= {R_CUT:.4f}):")
    # AD-vs-FD cross-check
    _, dJ_ad = compute_dJ(x1[:64], 1)
    _, dJ_fd = compute_dJ_fd(x1[:64], 1)
    rel = (dJ_ad - dJ_fd).norm() / dJ_fd.norm()
    print(f"    dJ cross-check (AD vs FD, 64 pts): rel err {rel.item():.2e}")

    nij_flat = {}
    for patch, xp in [(1, x1), (2, x2)]:
        mask = xp.norm(dim=1) <= R_CUT
        xa = xp[mask]
        flats, rounds = [], []
        for s in range(0, len(xa), CHUNK):
            xc = xa[s:s + CHUNK]
            J, dJ = compute_dJ(xc, patch)
            N = _nij_tensor(J, dJ)
            flats.append(torch.sqrt((N ** 2).sum(dim=(1, 2, 3)) + 1e-20))
            g, g_inv = EvenSphere.round_metric_ball_coords(xc)
            rounds.append(torch.sqrt(_metric_norm_sq_12tensor(N, g, g_inv) + 1e-20))
        flat = torch.cat(flats); rnd = torch.cat(rounds)
        nij_flat[patch] = flat.mean().item()
        print(f"    patch {patch}: flat-chart ||N|| mean {flat.mean():.6e} "
              f"(min {flat.min():.3e}, max {flat.max():.3e})  [{len(xa)} pts]")
        print(f"             round-metric ||N|| mean {rnd.mean():.6e}  "
              f"rel std {(rnd.std() / rnd.mean()).item():.2e}  "
              f"(G2 homogeneity: expect ~0)")
    mean_nij = 0.5 * (nij_flat[1] + nij_flat[2])
    print(f"\n    MEAN NIJ (patches avg): {mean_nij:.6e}")
    print(f"    [table: trained S^6 CS run = 9.87; certified-integrable would be ~0]")

    # ---------------- 3. basin diagnostic ----------------
    print("\n" + "=" * 78)
    print("Basin diagnostic (trained model vs octonionic ACS)")
    print("=" * 78)
    if args.checkpoint is None or not os.path.exists(args.checkpoint or ''):
        print("  no checkpoint supplied — SKIPPED.")
        print("  To run: fetch an S^6 run checkpoint from the HPC, e.g.")
        print("    rsync cenapad:~/homelovelace/NijenhuisPINN/outputs/<s6_run>/checkpoints/final_model.pt .")
        print("    python extra/octonion_calibration.py --checkpoint final_model.pt")
        return

    from src.models import TwoPatchAlmostComplexNet
    ckpt = torch.load(args.checkpoint, map_location='cpu', weights_only=False)
    net_cfg = ckpt['config']['network']
    model = TwoPatchAlmostComplexNet(
        manifold_dim=6,
        hidden_layers=net_cfg['hidden_layers'],
        activation=net_cfg['activation'],
        use_fourier_features=net_cfg['use_fourier_features'],
        num_fourier_features=net_cfg['num_fourier_features'],
        fourier_scale=net_cfg['fourier_scale'],
        initialization=net_cfg['initialization'],
        matrix_exp_terms=net_cfg.get('matrix_exp_terms'),
        use_residual_blocks=net_cfg.get('use_residual_blocks', False),
    ).double()
    state = ckpt.get('ema_state_dict', ckpt['model_state_dict'])
    model.load_state_dict(state)
    model.eval()

    mask = x1.norm(dim=1) <= R_CUT
    xe = x1[mask][:2048]
    with torch.no_grad():
        J_net = model.patch1(xe)
    J_oct = octonion_acs(xe, patch=1).detach()
    rel_val = ((J_net - J_oct).norm(dim=(1, 2)) / J_oct.norm(dim=(1, 2))).mean()
    print(f"  ||J_net − J_oct|| / ||J_oct|| (patch 1, value): {rel_val.item():.2%}")
    print(f"  → small (≲20%): trained field sits in the octonionic basin;")
    print(f"    large (≫50%): it found something else — investigate further.")


if __name__ == '__main__':
    main()
