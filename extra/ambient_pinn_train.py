"""Ambient-projector PINN training — chart-free construction of ACS/CS on the
non-parallelizable cases.

Cases:
  s4s2_acs   Build an ACS on S⁴×S² from scratch: K = P K̃_net P, single
             algebraic objective ‖K² + P‖²_F → 0. (Existence known; the old
             chart pipeline scored 1.91 because of cross-chart consistency —
             which the ambient formulation removes by construction.)
  s4s2_cs    OPEN PROBLEM exploration: same field, + ambient Nijenhuis term.
  s6_cs      CS on S⁶ (exists per Alpöge–Claude/Engel 2026): exact-ACS
             parametrization J = e^B (P M_oct P) e^{−B}, B = P B̃_net P —
             J² = −id_T exactly — single Nijenhuis objective. --sign -1 uses
             the opposite-orientation reference −J_oct.

Unbiased: the only training signals are the algebraic ACS condition and the
Nijenhuis residual. The octonionic reference in s6_cs is a base point of the
ACS bundle (like J₀ in the chart parametrization), not a target.

Run:  python extra/ambient_pinn_train.py <case> [--seed 7] [--sign 1]
                                         [--steps 3000] [--save out.pt]
"""
from __future__ import annotations

import os
import sys
import time
import argparse
import torch
import torch.nn as nn

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.geometry.ambient import S4xS2, S6
from src.geometry.octonion_acs import cross7_matrix
from src.geometry.octonion_product_acs import octonion_product_acs
from src.loss.ambient_loss import (ambient_nij_norm, compute_J_and_DJ,
                                   compute_P_and_DP)

DTYPE = torch.float64
HIDDEN = [128, 128, 128]   # default; override with --hidden
BATCH = 192


class AmbientNet(nn.Module):
    """MLP: ambient point → `stages` raw n×n matrix fields."""

    def __init__(self, n_in: int, n_mat: int, out_scale: float = 0.3, stages: int = 1,
                 hidden=None):
        super().__init__()
        layers, prev = [], n_in
        for h in (hidden or HIDDEN):
            layers += [nn.Linear(prev, h), nn.GELU()]
            prev = h
        layers += [nn.Linear(prev, stages * n_mat * n_mat)]
        self.net = nn.Sequential(*layers)
        with torch.no_grad():
            self.net[-1].weight.mul_(out_scale)
            self.net[-1].bias.normal_(0.0, out_scale)
        self.n_mat = n_mat
        self.stages = stages

    def forward(self, x):
        out = self.net(x).reshape(-1, self.stages, self.n_mat, self.n_mat)
        return out


def make_J_fn(case, model, manifold, sign=1.0):
    if case == 's4s2_cs_conj':
        # CS search on S4xS2 from the explicit octonionic-product ACS base point:
        # J = e^{B_k}...e^{B_1} J_ref e^{-B_1}...e^{-B_k}, exact J^2 = -P.
        # Same structure as the S6 search; the base point is an ACS, not a CS.
        def J_fn(x):
            P = manifold.projector(x)
            raw = model(x)
            J = octonion_product_acs(x, manifold.projector)
            for s in range(model.stages):
                Bs = P @ raw[:, s] @ P
                J = torch.matrix_exp(Bs) @ J @ torch.matrix_exp(-Bs)
            return J
    elif case.startswith('s4s2'):
        def J_fn(x):
            P = manifold.projector(x)
            return P @ model(x)[:, 0] @ P
    else:  # s6_cs
        # NOTE (red-team verified): do NOT symmetrize or skew-project the net
        # output — restricting e^B toward SO(7) collapses the search to
        # round-orthogonal J, which LeBrun's theorem excludes for CS on S6.
        # Two nested exponentials (stages=2) make the parametrization reach
        # every smooth same-orientation ACS field (polar canonical form);
        # a single exponential does not (exp: gl -> GL+ is not surjective).
        n = manifold.ambient_dim
        def J_fn(x):
            P = manifold.projector(x)
            raw = model(x)
            Joct = sign * (P @ cross7_matrix(x) @ P)
            J = Joct
            for s in range(model.stages):
                Bs = P @ raw[:, s] @ P
                J = torch.matrix_exp(Bs) @ J @ torch.matrix_exp(-Bs)
            return J
    return J_fn


def acs_residual(J_fn, manifold, x):
    """Per-sample ||K^2 + P||_F (no derivatives needed)."""
    J = J_fn(x)
    P = manifold.projector(x)
    return torch.sqrt(((J @ J + P) ** 2).sum(dim=(1, 2)) + 1e-20)


def losses(case, J_fn, manifold, x, w_c=10.0, power=2.0):
    out = {}
    if case == 's4s2_acs':
        # Algebraic objective only: no DJ needed. `power` > 2 weights the
        # residual tail (localized defect regions that dominate the sup norm
        # while barely moving the mean-squared loss).
        r = acs_residual(J_fn, manifold, x)
        out['total'] = (r ** power).mean()
        out['acs'] = (r ** 2).mean()
        return out
    J, DJ = compute_J_and_DJ(J_fn, x)
    P, DP = compute_P_and_DP(manifold, x)
    if False:
        pass
    elif case == 's4s2_cs':
        acs = ((J @ J + P) ** 2).sum(dim=(1, 2)).mean()
        nij = (ambient_nij_norm(J, DJ, P, DP) ** 2).mean()
        out['acs'], out['nij'] = acs, nij
        out['total'] = w_c * acs + nij
    else:  # s6_cs / s4s2_cs_conj — J²=−P exact by construction
        nij = (ambient_nij_norm(J, DJ, P, DP) ** 2).mean()
        out['nij'] = nij
        out['total'] = nij
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('case', choices=['s4s2_acs', 's4s2_cs', 's4s2_cs_conj', 's6_cs'])
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--sign', type=float, default=1.0)
    ap.add_argument('--steps', type=int, default=3000)
    ap.add_argument('--lbfgs', type=int, default=150)
    ap.add_argument('--stages', type=int, default=1)
    ap.add_argument('--hidden', default=None, help='comma-separated hidden widths, e.g. 256,256,256,256')
    ap.add_argument('--power', type=float, default=2.0, help='ACS residual power (tail weighting)')
    ap.add_argument('--rar', action='store_true', help='residual-adaptive sampling (ACS case)')
    ap.add_argument('--save', default=None)
    ap.add_argument('--init', default=None, help='warm-start from a saved checkpoint')
    args = ap.parse_args()
    torch.manual_seed(args.seed)

    manifold = S4xS2 if args.case.startswith('s4s2') else S6
    n = manifold.ambient_dim
    hidden = [int(h) for h in args.hidden.split(',')] if args.hidden else None
    model = AmbientNet(n, n, stages=args.stages, hidden=hidden).to(DTYPE)
    if args.init:
        model.load_state_dict(torch.load(args.init, map_location='cpu', weights_only=False)['model_state_dict'])
        print(f'warm-started from {args.init}', flush=True)
    J_fn = make_J_fn(args.case, model, manifold, sign=args.sign)
    tag = f"{args.case} seed={args.seed}" + (f" sign={args.sign:+.0f}" if args.case == 's6_cs' else "")

    x0 = manifold.sample(400)
    with torch.no_grad():
        pass
    init = {k: v.item() for k, v in losses(args.case, J_fn, manifold, x0, power=args.power).items()}
    print(f"[{tag}] init: " + "  ".join(f"{k}={v:.4e}" for k, v in init.items()), flush=True)

    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=args.steps, eta_min=1e-5)
    t0 = time.time()
    for step in range(args.steps):
        if args.rar and args.case == 's4s2_acs':
            # Residual-adaptive sampling: half the batch from the worst
            # residuals among 4x candidates, half uniform.
            with torch.no_grad():
                xc = manifold.sample(4 * BATCH)
                r = acs_residual(J_fn, manifold, xc)
                idx = torch.topk(r, BATCH // 2).indices
            x = torch.cat([xc[idx], manifold.sample(BATCH - BATCH // 2)], dim=0)
        elif args.rar:
            # Nijenhuis-residual RAR (attacks the sup tail): score 2x candidates
            # by ||N|| without building a backward graph, keep the worst half.
            xc = manifold.sample(2 * BATCH)
            Jc, DJc = compute_J_and_DJ(J_fn, xc)
            Pc, DPc = compute_P_and_DP(manifold, xc)
            r = ambient_nij_norm(Jc.detach(), DJc.detach(), Pc, DPc)
            idx = torch.topk(r, BATCH // 2).indices
            x = torch.cat([xc[idx], manifold.sample(BATCH - BATCH // 2)], dim=0)
        else:
            x = manifold.sample(BATCH)
        out = losses(args.case, J_fn, manifold, x, power=args.power)
        opt.zero_grad()
        out['total'].backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step()
        sched.step()
        if step % 300 == 0 or step == args.steps - 1:
            print(f"[{tag}] step {step:5d} | " +
                  "  ".join(f"{k}={v.item():.4e}" for k, v in out.items()), flush=True)

    xf = manifold.sample(768)
    lb = torch.optim.LBFGS(model.parameters(), lr=1.0, max_iter=args.lbfgs,
                           history_size=100, line_search_fn='strong_wolfe')

    def closure():
        lb.zero_grad()
        l = losses(args.case, J_fn, manifold, xf, power=args.power)['total']
        l.backward()
        return l
    lb.step(closure)

    # held-out evaluation
    xt = manifold.sample(3000)
    J, DJ = compute_J_and_DJ(J_fn, xt)
    P, DP = compute_P_and_DP(manifold, xt)
    with torch.no_grad():
        acs_res = torch.sqrt(((J @ J + P) ** 2).sum(dim=(1, 2)) + 1e-20)
    nij = ambient_nij_norm(J.detach(), DJ.detach(), P, DP)
    print(f"\n[{tag}] HELD-OUT (3000 pts) after {time.time()-t0:.0f}s:")
    print(f"  ACS residual ||J^2+P||:   mean {acs_res.mean():.6e}  max {acs_res.max():.3e}")
    print(f"  round-metric ||N||:       mean {nij.mean():.6e}  max {nij.max():.3e}")
    if args.case == 's6_cs':
        with torch.no_grad():
            Joct = args.sign * (P @ cross7_matrix(xt) @ P)
            dist = (J - Joct).norm(dim=(1, 2)) / Joct.norm(dim=(1, 2))
        print(f"  distance from octonionic: mean {dist.mean():.3f}  "
              f"[octonionic ||N|| reference: 19.5959]")
    if args.case == 's4s2_cs_conj':
        with torch.no_grad():
            Jref = octonion_product_acs(xt, manifold.projector)
            dist = (J - Jref).norm(dim=(1, 2)) / Jref.norm(dim=(1, 2))
        print(f"  distance from octonionic-product ACS base point: mean {dist.mean():.3f}")
    if False:
        with torch.no_grad():
            Joct = args.sign * (P @ cross7_matrix(xt) @ P)
            dist = (J - Joct).norm(dim=(1, 2)) / Joct.norm(dim=(1, 2))
        print(f"  distance from octonionic: mean {dist.mean():.3f}  "
              f"[octonionic ||N|| reference: 19.5959]")
    if args.save:
        torch.save({'model_state_dict': model.state_dict(), 'case': args.case,
                    'seed': args.seed, 'sign': args.sign, 'stages': args.stages,
                    'hidden': hidden}, args.save)
        print(f"  saved: {args.save}")


if __name__ == '__main__':
    main()
