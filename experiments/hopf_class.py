"""Homotopy class of the S^3 x S^1 fields of the campaign relative to the Hopf target of Table 4
(extra/hopf_invariant.py): the ACS bases, the exact complex structures of the polynomial CS runs
and the conjugation runs. A base with a nonzero class cannot be conjugated to the target.

Every field is taken with the orientation fix of extra/ambient_supervised_hopf.py (identity when
its orientation is the target's, else the reflection theta -> -theta) and classified on two slices
S^3 x {theta}. Writes OUT/done.json.

    python experiments/hopf_class.py --runs outputs/paper_rerun/runs --out outputs/paper_rerun/runs/hopf_class
"""
import os, sys, json, argparse
import torch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extra.campaign_unit import Geometry, load_acs_field, polar_corrected, DTYPE
from extra.ambient_supervised_hopf import target_fn
from extra.hopf_invariant import C_MAIN, sphere_mesh, relative_class


def field(geo, runs, name):
    """J of an ACS unit (polar-corrected base) or of a CS run with an exact field (None otherwise)."""
    run = os.path.join(runs, name)
    if name.startswith('acs_'):
        return polar_corrected(geo, load_acs_field(geo, run))
    from cs_analysis import field_of
    return field_of(geo, run)[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--runs', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--n', type=int, default=32, help='mesh cells per edge of the 4-cube (8 x 6 n^3 tets)')
    ap.add_argument('--thetas', default='0.3,2.5')
    a = ap.parse_args()
    torch.set_num_threads(int(os.environ.get('OMP_NUM_THREADS', 1)))
    os.makedirs(a.out, exist_ok=True)
    names = sorted(n for n in os.listdir(a.runs)
                   if (n.startswith(('acs_s3xs1_', 'cspoly_s3xs1_')) or n in ('csconjexp_s3xs1', 'csconjaff_s3xs1'))
                   and os.path.exists(os.path.join(a.runs, n, 'final.pt')))
    geo = Geometry('s3xs1')
    Jt = target_fn(2, C_MAIN)
    mesh = sphere_mesh(a.n)
    S = torch.ones(geo.N, dtype=DTYPE)
    S[4] = -1.0
    S = torch.diag(S)
    out = {'mesh_n': a.n, 'tets': int(len(mesh[1])), 'units': {}}
    for name in names:
        J = field(geo, a.runs, name)
        if J is None:
            out['units'][name] = {'skipped': 'no exact field'}
            continue
        rec = {'slices': {}}
        for th in [float(t) for t in a.thetas.split(',')]:
            try:
                inv, info = relative_class(J, Jt, theta=th, mesh=mesh)
                rec['orientation_fix'] = 'identity'
            except ValueError:                   # the supervised script's fix for n even: theta -> -theta
                inv, info = relative_class(lambda x: S @ J(x @ S) @ S, Jt, theta=th, mesh=mesh)
                rec['orientation_fix'] = 'reflect'
            rec['slices'][f'{th:g}'] = dict(info, hopf_invariant=inv)
            print(f"{name:24s} theta {th:g}: fix {rec['orientation_fix']:8s} angle to J_t "
                  f"{info['angle_to_target_deg'][0]:5.1f}..{info['angle_to_target_deg'][1]:5.1f} deg, "
                  f"linking {[round(t['linking'], 3) for t in info['trials']]}", flush=True)
        inv = {s['hopf_invariant'] for s in rec['slices'].values()}
        rec['hopf_invariant'] = inv.pop() if len(inv) == 1 else None
        out['units'][name] = rec
    json.dump(out, open(os.path.join(a.out, 'done.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
