"""Period-map PINN study: convergence in the truncation K and the effect of the
collocation set (fixed-point discs only vs discs + Ford-domain gluing pairs);
grids of tau, mu, beta, D over the fundamental domain for the figures."""
import os, sys, json, time, argparse, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from extra.period_map import *
ap = argparse.ArgumentParser(); ap.add_argument('--out', required=True); a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
G = TriangleGroup()
Zf = sample_ford(G, 20000, np.random.default_rng(1)); qmax = float(np.exp(-2 * np.pi * Zf.imag.min()))
rho = np.exp(1j * np.pi / 3)
pairs = edge_pairs(G, 3000, np.random.default_rng(2))
def held_out(pm):
    w = {}
    for name, Z in pairs.items():
        g = {'g1': G.g1, 'g2': G.g2, 'g1inv': np.linalg.inv(G.g1), 'g2inv': np.linalg.inv(G.g2)}[name]
        W = Z if name in ('g1', 'g2') else mobius(g, Z)
        r = functional_residuals(pm, G, W)
        w[name] = float(max(v for k, v in r.items() if k.endswith('g1' if '1' in name else 'g2')))
    return w
study = []
for colloc in ('discs', 'discs+edges'):
    for K in (8, 12, 16, 20, 24, 28, 32):
        # the same basis scale q0 = qmax (the Ford-domain value, 0.156) for both collocation sets
        # (until 2026-10-25 the discs-only curve used q0 = 0.2: two settings changed at once)
        pm = PeriodMap(K=K, q0=qmax); pm.a[0] = 0.3 + 0.6j
        P1, P2 = (collocation(G, n=64, radius=0.25, rings=3) if colloc == 'discs' else collocation_full(G))
        t0 = time.time(); L = solve_tau(pm, G, P1, P2, iters=150)
        solve_mu(pm, G, P1, P2); solve_beta(pm, G, P1, P2, c0=-3j)
        rec = dict(collocation=colloc, K=K, tau_loss=float(L), seconds=time.time() - t0, held_out=held_out(pm),
                   tau_Z1=float(abs(pm.tau(np.array([G.Z1]))[0] - rho)), tau_Z2=float(abs(pm.tau(np.array([G.Z2]))[0] - 1j)),
                   mu_Z1=float(abs(pm.mu(np.array([G.Z1]))[0] - 1 / (1 + rho))),
                   mu_Z2=float(abs(pm.mu(np.array([G.Z2]))[0] - (1 - 1j) / 2)),
                   coef_a=[float(abs(v)) for v in pm.a], coef_b=[float(abs(v)) for v in pm.b], coef_c=[float(abs(v)) for v in pm.cb])
        study.append(rec); print(colloc, K, f"{max(rec['held_out'].values()):.1e}", flush=True)
        if colloc == 'discs+edges' and K == 20:
            xs, ys = np.linspace(-0.6, 0.6, 241), np.linspace(0.2, 1.6, 281)
            XX, YY = np.meshgrid(xs, ys); ZZ = XX + 1j * YY
            inside = ford_domain_mask(G, ZZ.ravel()).reshape(ZZ.shape)
            np.savez(os.path.join(a.out, 'grids_K20.npz'), X=XX, Y=YY, inside=inside,
                     tau=pm.tau(ZZ.ravel()).reshape(ZZ.shape), mu=pm.mu(ZZ.ravel()).reshape(ZZ.shape),
                     beta=pm.beta(ZZ.ravel()).reshape(ZZ.shape), D=D_function(pm, ZZ.ravel()).reshape(ZZ.shape),
                     Z1=G.Z1, Z2=G.Z2, P1=P1, P2=P2, qmax=qmax)
json.dump({'study': study, 'qmax': qmax, 'Z1': [G.Z1.real, G.Z1.imag], 'Z2': [G.Z2.real, G.Z2.imag]},
          open(os.path.join(a.out, 'done.json'), 'w'), indent=1)
