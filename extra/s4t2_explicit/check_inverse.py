"""F^-1 on uniformly distributed points of the round S^4 x T^2: residual |F(F^-1(q)) - q|, the initial-guess quality
(the explicit chain is exact outside the linear-blend shells), and the transported J evaluated there.
Run from the repo root: python3 extra/s4t2_explicit/check_inverse.py [n]"""
from __future__ import annotations
import math, os, sys, time
import torch
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
torch.set_default_dtype(torch.float64)
from extra.s4t2_explicit import tinverse as TI, nijenhuis_round as NR

if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    torch.manual_seed(5)
    q = NR.S4T2.sample(n)
    t0 = time.time()
    tt, lg, x = TI.initial_guess(q)
    guess_err = 0.0
    for key, (mask, P) in TI.to_chart(tt, lg, x).items():
        if mask.any():
            guess_err = max(guess_err, (TI.CHART_F[key](P[mask]) - q[mask]).abs().max().item())
    sol = TI.Finv(q)
    res = torch.cat([v[2] for v in sol.values()])
    print(f'{n} uniform points of S^4 x T^2: charts used ' + ', '.join(f'{k}: {len(v[0])}' for k, v in sol.items()))
    print(f'   explicit chain (before Newton): max |F(p0) - q| = {guess_err:.1e}')
    print(f'   after Newton polish: max |F(p) - q| = {res.max():.1e}, median {res.median():.1e}  ({time.time() - t0:.0f} s)')
    worst = []
    for key, (idx, p, r) in sol.items():
        qq, A, H = NR.DF_and_D2F(TI.CHART_F[key], p)
        J, DJ = NR.transport(A, H, qq)
        Nn, scale, Pq = NR.nij(qq, J, DJ)
        worst.append((Nn / scale).max().item())
        print(f'   {key}: J at F^-1(q): |J^2 + P| max {(J @ J + Pq).flatten(1).norm(dim=-1).max():.1e}; '
              f'|N| / scale max {(Nn / scale).max():.1e}; |J| median {J.flatten(1).norm(dim=-1).median():.1f}')
    ok = res.max() < 1e-10
    print('PASS' if ok else 'FAIL')
