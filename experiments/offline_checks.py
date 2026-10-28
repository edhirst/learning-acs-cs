"""Offline checks (never used in training): (i) the Hopf structures on S^3xS^1,
S^5xS^1, the Calabi-Eckmann structure on S^3xS^3 and the product structure on
S^2xS^2 are exact integrable ACS and lie in the degree-2 polynomial class;
(ii) the orientation index o = Pf(E^T (K^T - K) E)/2^{d/2} of every certified
ACS has |o| >= 1 and one sign; (iii) orthogonality defect of the found CS."""
import os, sys, json, argparse, numpy as np, torch
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from extra.campaign_unit import Geometry, acs_residual, nij_norm, load_acs_field, polar_corrected, orientation_index, base_acs
from extra.test_poly_cs import _hopf_M, _ce_M
from src.models.poly_field import PolyField
ap = argparse.ArgumentParser(); ap.add_argument('--out', required=True); ap.add_argument('--acs-root', required=True)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
torch.manual_seed(0)
out = {'classical': {}, 'orientation': {}}
for case, Mf in (('s3xs1', lambda p: _hopf_M(p, 4)), ('s5xs1', lambda p: _hopf_M(p, 6)), ('s3xs3', _ce_M),
                 ('s2xs2', None)):
    geo = Geometry(case, 'ambient'); man = geo.manifold
    if Mf is None:
        Mf = lambda p: base_acs('s2xs2', p, man.projector)
    K = lambda p: man.projector(p) @ Mf(p) @ man.projector(p)
    x = man.sample(400)
    f = PolyField(geo.N, 2); xf = man.sample(3000); F = f.features(xf); T = Mf(xf).reshape(len(xf), -1)
    W = torch.linalg.lstsq(F, T).solution
    out['classical'][case] = {'acs_sup': acs_residual(geo, K, x).max().item(), 'N_sup': nij_norm(geo, K, x).max().item(),
                              'deg2_fit_residual': (F @ W - T).abs().max().item()}
for case in ['s2', 's6', 's1xs1', 's2xs2', 's3xs1', 's3xs3', 's4xs2', 's5xs1', 's4xs1xs1']:
    d = os.path.join(a.acs_root, f'acs_{case}_s7')
    geo = Geometry(case, 'ambient')
    rj = os.path.join(d, 'result.json')
    if not os.path.exists(rj) or not json.load(open(rj)).get('certified_exact_acs'):
        out['orientation'][case] = {'error': 'no certified ACS'}; continue
    try:
        J = polar_corrected(geo, load_acs_field(geo, d))
    except Exception as e:
        out['orientation'][case] = {'error': str(e)}; continue
    o = orientation_index(geo, J, geo.sample(4000)).detach().numpy()
    out['orientation'][case] = {'o_min_abs': float(np.abs(o).min()), 'o_min': float(o.min()), 'o_max': float(o.max()),
                                'one_sign': bool((np.sign(o) == np.sign(o[0])).all()),
                                'hist': np.histogram(o, bins=40)[0].tolist(), 'edges': np.histogram(o, bins=40)[1].tolist()}
    print(case, out['orientation'][case]['o_min_abs'], flush=True)
json.dump(out, open(os.path.join(a.out, 'done.json'), 'w'), indent=1)
