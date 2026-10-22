"""Reconstruct and certify the complex structures on S^6 and S^2 x S^4.

    python extra/reconstruct_s6.py --out outputs/reconstruction [--K 20]

1. Learn the period map (tau, mu, beta) of [AC26, §3] as exact q-series by
   residual minimisation of its functional equations (extra/period_map.py):
   Levenberg-Marquardt for tau, linear least squares for mu and for beta (up to
   its free constant c0, the family parameter; Im c0 < 0 large enough for
   D < 0).
2. Build the complex structure J on the fibration model (h_Z x R^4/Lambda)/Delta
   (S^6: Lambda; S^2 x S^4: the index-two superlattice Lambda+ of [GFL26]) and
   certify it (extra/s6_reconstruction.py): |N| on the fundamental domain and
   on the multiple-fibre discs (with adversarial sup), deck consistency on the
   gluing pairs, |N| of the blended global field across the identified edges,
   D < 0 (lattices), the functional-equation residuals of the period map.

This is a RECONSTRUCTION: the functional equations, the lattices and the
fibred function class are taken from [AC26] / [GFL26]; no value of the period
map or of J is given to the solver. The cusp filling is the toric quotient of
[AC26, §4] built from the same q-series (standard toric structure: N = 0
identically), the multiple fibres are the log transforms of [AC26, §5] (the
x-invariant J is unchanged by the twisted translations), and X = S^6, S^2 x S^4
diffeomorphically is the topological theorem of those papers.
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import os
import sys
import json
import time
import argparse

import numpy as np
import torch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from extra.period_map import (TriangleGroup, PeriodMap, collocation, collocation_full, solve_tau,
                              solve_mu, solve_beta, sample_ford, edge_pairs, functional_residuals,
                              D_function, mobius, A1, A2)
from extra.s6_reconstruction import (TorchPeriodMap, J_field, J_of_Pi, nijenhuis, mobius_t,
                                     deck_differential, blended_field, adversarial_sup_2d)

DTYPE = torch.float64
C128 = torch.complex128


def learn_period_map(K=20, c0=-3j, log=print):
    G = TriangleGroup()
    Zf = sample_ford(G, 20000, np.random.default_rng(1))
    pm = PeriodMap(K=K, q0=float(np.exp(-2 * np.pi * Zf.imag.min())))
    pm.a[0] = 0.3 + 0.6j                     # tau ~ Z + const at the cusp; no target values
    P1, P2 = collocation_full(G)
    t0 = time.time()
    L = solve_tau(pm, G, P1, P2, iters=150)
    e_mu = solve_mu(pm, G, P1, P2)
    e_b = solve_beta(pm, G, P1, P2, c0=c0)
    log(f"period map: tau loss {L:.1e}, mu fit {e_mu:.1e}, beta fit {e_b:.1e} ({time.time()-t0:.1f}s)")
    return G, pm


def certify(G, pm, pmt, lattice=None, seed=20261002):
    rng = np.random.default_rng(seed)
    rep = {}
    if lattice is None:
        Ap = {1: A1, 2: A2}
    else:
        L = lattice.numpy()
        Ap = {1: np.linalg.inv(L) @ A1 @ L, 2: np.linalg.inv(L) @ A2 @ L}
        rep['monodromy_integral'] = bool(all(np.allclose(A, np.round(A)) for A in Ap.values()))
        rep['A1'], rep['A2'] = np.round(Ap[1]).astype(int).tolist(), np.round(Ap[2]).astype(int).tolist()
    Jf = lambda p: J_field(pmt, p, lattice)
    Zf = sample_ford(G, 4000, rng)
    ab = torch.tensor(np.stack([Zf.real, Zf.imag], 1), dtype=DTYPE)
    N = torch.cat([nijenhuis(Jf, ab[i:i + 500]).detach().flatten(1).norm(dim=1) for i in range(0, len(ab), 500)])
    J = Jf(ab)
    rep['acs_residual_max'] = (J @ J + torch.eye(6, dtype=DTYPE)).abs().max().item()
    rep['N_mean'], rep['N_max'] = N.mean().item(), N.max().item()
    rep['N_adversarial_sup'] = adversarial_sup_2d(lambda p: nijenhuis(Jf, p).flatten(1).norm(dim=1),
                                                  ab[torch.topk(N, 64).indices])
    pairs = edge_pairs(G, 2000, rng)
    P1, P2 = collocation(G, n=128, radius=0.25, rings=4)
    ymin = min(G.Z1.imag, G.Z2.imag) - 0.01
    for j, g, sets in ((1, G.g1, {'pairs': pairs['g1'], 'disc': P1}), (2, G.g2, {'pairs': pairs['g2'], 'disc': P2})):
        for nm, Zs in sets.items():
            Z = torch.tensor(Zs, dtype=C128)
            D = deck_differential(g, Z, Ap[j])
            err = (D @ J_of_Pi(pmt.Pi(Z), lattice) @ torch.linalg.inv(D)
                   - J_of_Pi(pmt.Pi(mobius_t(g, Z)), lattice)).flatten(1).norm(dim=1)
            rep[f'deck_g{j}_{nm}_max'] = err.max().item()
        abd = torch.tensor(np.stack([sets['disc'].real, sets['disc'].imag], 1), dtype=DTYPE)
        rep[f'N_disc_Z{j}_max'] = nijenhuis(Jf, abd).detach().flatten(1).norm(dim=1).max().item()
        Zs = pairs[f'g{j}'][:600]
        Jb = blended_field(pmt, g, Ap[j], centre=np.median(Zs.real), kappa=30.0, lattice=lattice, Alat=Ap[j])
        abz = torch.tensor(np.stack([Zs.real, Zs.imag], 1), dtype=DTYPE)
        Nb = torch.cat([nijenhuis(Jb, abz[i:i + 200]).detach().flatten(1).norm(dim=1) for i in range(0, len(abz), 200)])

        def valid(p, g=g):
            Zp = torch.complex(p[:, 0], p[:, 1])
            return (Zp.imag >= ymin) & (mobius_t(g, Zp).imag >= ymin)
        rep[f'blend_g{j}_N_max'] = Nb.max().item()
        rep[f'blend_g{j}_N_adversarial'] = adversarial_sup_2d(lambda p: nijenhuis(Jb, p).flatten(1).norm(dim=1),
                                                             abz[torch.topk(Nb, 32).indices], steps=40, valid=valid)
    rep['D_max_on_domain'] = float(D_function(pm, Zf).max())
    r1, r2 = functional_residuals(pm, G, pairs['g1']), functional_residuals(pm, G, pairs['g2'])
    rep['functional_residual_max'] = float(max(max(v for k, v in r1.items() if k.endswith('g1')),
                                              max(v for k, v in r2.items() if k.endswith('g2'))))
    rho = np.exp(1j * np.pi / 3)
    rep['tau_Z1_minus_rho'] = float(abs(pm.tau(np.array([G.Z1]))[0] - rho))
    rep['tau_Z2_minus_i'] = float(abs(pm.tau(np.array([G.Z2]))[0] - 1j))
    rep['mu_Z1_minus_1_over_1_plus_rho'] = float(abs(pm.mu(np.array([G.Z1]))[0] - 1 / (1 + rho)))
    rep['mu_Z2_minus_(1-i)/2'] = float(abs(pm.mu(np.array([G.Z2]))[0] - (1 - 1j) / 2))
    keys = [k for k in rep if k.startswith(('N_', 'blend_', 'deck_'))]
    rep['certified_sup'] = max(rep[k] for k in keys if isinstance(rep[k], float))
    return rep


def cusp_check(pm, r=0.3, n=4000, seed=5):
    """Hypothesis of [AC26, Thm 4.5] for the toric filling at the cusp, with the
    learned period map: in the cusp normal form Z = s B0 + C(q), s = Z,
    C = [[6 mu, tau - s], [beta + s, mu]] must be holomorphic at q = 0 (it is a
    q-series) and R = -2 pi Im C bounded on |q| <= r/2; then the Lambda-action on
    Y_eps is free and properly discontinuous for eps <= eps0 = min(r/2, exp(-2 |R|_inf)),
    and N0 = Y_eps / Lambda carries the standard toric structure (N = 0)."""
    rng = np.random.default_rng(seed)
    rad = (r / 2) * np.sqrt(rng.uniform(0, 1, n)); ang = rng.uniform(0, 2 * np.pi, n)
    q = rad * np.exp(1j * ang)
    q = q[np.abs(q) > 1e-12]
    Z = np.log(q) / (2j * np.pi)                     # any branch: C depends on q only
    tau, mu, beta = pm.tau(Z), pm.mu(Z), pm.beta(Z)
    C = np.stack([np.stack([6 * mu, tau - Z], -1), np.stack([beta + Z, mu], -1)], -2)
    R = -2 * np.pi * C.imag
    Rnorm = np.linalg.norm(R, ord=2, axis=(1, 2)).max()
    # single-valuedness in q: C at Z and Z - 1 (the cusp generator) must coincide
    Cs = np.stack([np.stack([6 * pm.mu(Z - 1), pm.tau(Z - 1) - (Z - 1)], -1),
                   np.stack([pm.beta(Z - 1) + (Z - 1), pm.mu(Z - 1)], -1)], -2)
    return {'R_inf': float(Rnorm), 'eps0': float(min(r / 2, np.exp(-2 * Rnorm))),
            'C_single_valued_err': float(np.abs(C - Cs).max()), 'C0': np.round(C[np.abs(q).argmin()], 6).tolist()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--K', type=int, default=20)
    ap.add_argument('--c0', type=float, default=-3.0, help='Im of the free constant of beta (D < 0 needs it negative)')
    ap.add_argument('--target', type=float, default=1e-6)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    torch.set_num_threads(int(os.environ.get('OMP_NUM_THREADS', 4)))
    G, pm = learn_period_map(args.K, c0=1j * args.c0)
    np.savez(os.path.join(args.out, 'period_map.npz'), a=pm.a, b=pm.b, cb=pm.cb, K=pm.K, q0=pm.q0)
    pmt = TorchPeriodMap(os.path.join(args.out, 'period_map.npz'))
    report = {'cusp': cusp_check(pm)}
    report['cusp']['C0'] = str(report['cusp']['C0'])
    print(f"cusp filling [AC26 Thm 4.5]: |R|_inf {report['cusp']['R_inf']:.3f} -> eps0 {report['cusp']['eps0']:.2e} "
          f"(> 0: the toric quotient exists); C single-valued in q to {report['cusp']['C_single_valued_err']:.1e}", flush=True)
    for case, lattice in (('S6', None), ('S2xS4', torch.tensor(np.diag([1, 1, 1, 0.5]), dtype=DTYPE))):
        rep = certify(G, pm, pmt, lattice)
        rep['target'] = args.target
        rep['certified'] = rep['certified_sup'] < args.target and rep['D_max_on_domain'] < 0
        report[case] = rep
        print(f"{case}: certified sup |N| (charts, discs, deck, blended) {rep['certified_sup']:.2e} "
              f"| D max {rep['D_max_on_domain']:.2f} | {'CERTIFIED' if rep['certified'] else 'NOT certified'} "
              f"(target {args.target:g})", flush=True)
    with open(os.path.join(args.out, 'certificate.json'), 'w') as f:
        json.dump(report, f, indent=2)


if __name__ == '__main__':
    main()
