"""Identify the endpoint of an S³×S¹ run of extra/ambient_supervised_hopf.py (used for the Phase B
analysis in docs/ambient_supervised_hopf_log.md, Stage 4).

Reports, for the field J of RUN_DIR (final.pt if present, else state.pt):
  * intrinsic ‖N‖, ‖J‖_F (‖J‖²_F ≥ d for any ACS, with equality iff J is orthogonal), Dirichlet energy;
  * the closest member of the Hopf family J_c (src/geometry/hopf_ambient.py), c scanned on both signs,
    and the distances to main's J_{1/(2π)} and to the round-orthogonal J_1;
  * whether J is CONSTANT in the global Hopf frame F_k = (e_k − q_k q, q_k ∂_θ): every constant linear
    complex structure K on ℝ⁴ = ℂ² gives, in this frame, an integrable structure (z ∼ e^{2π} z).

Usage:  python extra/ambient_phaseB_identify.py RUN_DIR BASE_DIR [--hidden 128,128,128] [--stages 2]
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import os
import sys
import math
import argparse
import torch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from extra.campaign_unit import Geometry, Net, make_field, nij_norm, dirichlet_energy, DTYPE
from extra.ambient_supervised_hopf import build_base, target_fn, hopf_frame


def rel_rms(A, B):
    return ((A - B).norm(dim=(1, 2)).pow(2).mean() / B.norm(dim=(1, 2)).pow(2).mean()).sqrt().item()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_dir')
    ap.add_argument('base_dir')
    ap.add_argument('--hidden', default='128,128,128')
    ap.add_argument('--stages', type=int, default=2)
    args = ap.parse_args()

    geo, n = Geometry('s3xs1'), 2
    Jt = target_fn(n, 1.0 / (2.0 * math.pi))
    base_fn, _, _ = build_base(geo, n, args.base_dir, 'auto', Jt, 'identify')
    path = os.path.join(args.run_dir, 'final.pt')
    if not os.path.exists(path):
        path = os.path.join(args.run_dir, 'state.pt')
    ck = torch.load(path, map_location='cpu', weights_only=False)
    model = Net(geo.N, geo.N, args.stages, [int(h) for h in args.hidden.split(',')], 0.5).to(DTYPE)
    model.load_state_dict(ck['model'])
    J_fn, _ = make_field(geo, 'cs', model, base_fn=base_fn)

    torch.manual_seed(1)
    x = geo.sample(3000)
    with torch.no_grad():
        J = J_fn(x)
    print(f"{path} (step {ck.get('step', 'final')}): ‖N‖ mean {nij_norm(geo, J_fn, x[:500]).mean():.2e} | "
          f"‖J‖_F mean {J.norm(dim=(1, 2)).mean():.3f} | Dirichlet energy {dirichlet_energy(geo, J_fn, x[:500]).item():.2f}")

    with torch.no_grad():
        best = None
        for sgn in (1.0, -1.0):
            for c in torch.logspace(-1.5, 1.0, 60).tolist():
                d = rel_rms(J, target_fn(n, sgn * c)(x))
                if best is None or d < best[0]:
                    best = (d, sgn * c)
        print(f"closest Hopf J_c': c' = {best[1]:+.4f} (1/c' = {1 / best[1]:+.3f}), rel. RMS distance {best[0]:.3e}")
        for c in (1.0 / (2.0 * math.pi), 1.0):
            print(f"   distance to J_c, c = {c:+.4f}: {rel_rms(J, target_fn(n, c)(x)):.3e}   "
                  f"(Dirichlet energy of J_c: {dirichlet_energy(geo, target_fn(n, c), x[:500]).item():.2f})")

        F = hopf_frame(x, n)
        Kf = F @ J @ F.transpose(1, 2)
        Kbar = Kf.mean(0)
        I4 = torch.eye(4, dtype=DTYPE)
        print(f"frame components: mean relative deviation from their average "
              f"{((Kf - Kbar).norm(dim=(1, 2)) / Kbar.norm()).mean().item():.3e}")
        print(f"   average K̄: ‖K̄² + I‖ {(Kbar @ Kbar + I4).abs().max():.2e}, ‖K̄K̄ᵀ − I‖ {(Kbar @ Kbar.T - I4).abs().max():.2e}")
        print(K_str(Kbar))


def K_str(K):
    return '\n'.join('   ' + '  '.join(f'{v:+.4f}' for v in row) for row in K.tolist())


if __name__ == '__main__':
    main()
