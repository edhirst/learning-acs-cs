"""Every number quoted in the paper "Learning almost complex and complex structures",
read from the run records, in the order of the paper.

    .venv/bin/python experiments/paper_numbers.py [--runs DIR] [--json OUT.json]

Default run directory: $PAPER_OUT/runs (outputs/paper_rerun/runs). The keys name the
paper's section, table or figure; the values are rounded only when printed.
"""
from __future__ import annotations    # PEP 563 — for Python 3.8 HPC compatibility

import os
import sys
import json
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CLASSICAL_CS = ['s2', 's1xs1', 's2xs2', 's3xs1', 's3xs3', 's5xs1']
ACS_CASES = ['s2', 's1xs1', 's4', 's2xs2', 's3xs1', 's6', 's3xs3', 's5xs1', 's4xs2', 's4xs1xs1']
CLASS_LETTER = {'c': 'C', 'lin': 'L', 'bilin': 'B', 'quad_h': 'Q', 'lin+c': 'L+C', 'lin+bilin': 'L+B',
                'all2': 'A2', 'all2+trilin': 'T'}


def collect(mf):
    res, jload, corrected = mf.res, mf.jload, mf.corrected
    out = {}

    def put(key, value, src):
        out[key] = {'value': value, 'source': src}

    # ---- §4.1, Table 2, Figure 1: almost complex searches
    passed, adv_acc = [], []
    for c in ACS_CASES:
        for s in (7, 8):
            name = f'acs_{c}_s{s}'
            r = res(name)
            if not r:
                continue
            lad = [e for e in r.get('ladder', []) if e.get('phase') == 'poly']
            acc = [e for e in lad if e.get('adv_acs_sup', 9) < 1]
            put(f'T2.{c}.s{s}.class', CLASS_LETTER.get(r.get('feat'), r.get('feat')), f'{name}: feat')
            put(f'T2.{c}.s{s}.attempt', acc[-1]['attempt'] if acc else None, f'{name}: accepted ladder attempt')
            put(f'T2.{c}.s{s}.r_adv', r.get('adv_acs_sup'), f'{name}: adv_acs_sup (after refinement)')
            ok = bool(r.get('certified_exact_acs'))
            put(f'T2.{c}.s{s}.pass', ok, f'{name}: certified_exact_acs')
            if ok:
                passed.append(name)
                adv_acc.append(r['adv_acs_sup'])
                put(f'T2.{c}.s{s}.corrected_acs_adv', r['corrected'].get('adv_acs_sup'), f'{name}: corrected.adv_acs_sup')
            lb = [e for e in r.get('lbfgs', []) if e.get('phase') == 'lbfgs']
            if lb:
                put(f'T2.{c}.s{s}.lbfgs_iters_run', sum(e.get('iters_run', 0) for e in lb),
                    f'{name}: L-BFGS iterations actually run')
                put(f'T2.{c}.s{s}.lbfgs_kept', not all(e.get('reverted') for e in lb), f'{name}: L-BFGS kept')
    put('4.1.n_pass', len(passed), 'certified ACS units (of 20)')
    if adv_acc:
        put('4.1.r_adv_range', [min(adv_acc), max(adv_acc)], 'accepted units: adv_acs_sup range')
    for s in (7, 8):
        r = res(f'acs_s4_s{s}')
        if r:
            put(f'4.1.S4.s{s}.mean_acs', r['dense'].get('acs_mean'), f'acs_s4_s{s}: dense.acs_mean')
            put(f'4.1.S4.s{s}.adv_acs', r.get('adv_acs_sup'), f'acs_s4_s{s}: adv_acs_sup')
    r = res('acs_s6_s7')
    if r and r.get('corrected'):
        put('4.1.S6.acs_nij_mean', r['corrected'].get('nij_mean'), 'acs_s6_s7: corrected.nij_mean')

    # ---- §4.2, Table 3, Figure 2: classical complex structures
    deg3, below20, loss_all, nadv_all, mins = [], [], [], [], []
    for c in CLASSICAL_CS:
        for s in (7, 8):
            name = f'cspoly_{c}_s{s}'
            r = res(name)
            if not r:
                continue
            L = [e for e in r.get('ladder', []) if e.get('phase') == 'poly_cs']
            acc = [e for e in L if e.get('accepted')]
            e = acc[-1] if acc else (L[-1] if L else {})
            loss = e.get('final_N2', float('nan')) + 10 * e.get('final_C2', float('nan'))
            put(f'T3.{c}.s{s}.degree', e.get('degree'), f'{name}: accepted rung degree')
            put(f'T3.{c}.s{s}.accepted', bool(acc), f'{name}: ladder accepted')
            put(f'T3.{c}.s{s}.heldout_loss', loss, f'{name}: final_N2 + 10 final_C2 (2000 points, before correction)')
            put(f'T3.{c}.s{s}.N_adv', r.get('adv_nij_sup'), f'{name}: adv_nij_sup (corrected field, 128 x 75)')
            put(f'T3.{c}.s{s}.minutes', r.get('elapsed_min'), f'{name}: elapsed_min')
            if str(e.get('degree')) == '3':
                deg3.append(name)
            if loss <= 1e-20:
                below20.append(name)
            loss_all.append(loss); nadv_all.append(r.get('adv_nij_sup')); mins.append(r.get('elapsed_min'))
    put('4.2.n_runs', len(loss_all), 'classical polynomial runs')
    put('4.2.n_loss_le_1e-20', len(below20), 'runs with held-out loss <= 1e-20')
    put('4.2.degree3_runs', deg3, 'runs accepted at degree 3')
    if nadv_all:
        put('4.2.N_adv_max', max(nadv_all), 'max N_adv over Table 3 (abstract: "below 1e-9")')
        put('4.2.poly_minutes_range_4d6d', [min(mins), max(mins)], 'Table 3 minutes range')
    cm = []
    for c in ('s2xs2', 's3xs1', 's3xs3', 's5xs1'):
        for kind in ('exp', 'aff'):
            name = f'csconj{kind}_{c}'
            r = res(name)
            if r:
                put(f'4.2.conj{kind}.{c}.mean', r['dense']['nij_mean'], f'{name}: dense.nij_mean (6000)')
                put(f'4.2.conj{kind}.{c}.adv', r.get('adv_nij_sup'), f'{name}: adv_nij_sup')
                put(f'4.2.conj{kind}.{c}.minutes', r.get('elapsed_min'), f'{name}: elapsed_min')
                cm.append(r.get('elapsed_min'))
                lb = r.get('lbfgs', [])
                put(f'4.2.conj{kind}.{c}.lbfgs_kept', any(not e.get('reverted') for e in lb), f'{name}: any L-BFGS chunk kept')
                if r.get('conjugator'):
                    put(f'4.2.conj{kind}.{c}.g_sigma_min_adv', r['conjugator'].get('g_sigma_min_adv'),
                        f'{name}: conjugator.g_sigma_min_adv')
    if cm:
        put('4.2.conj_minutes_range', [min(cm), max(cm)], 'classical conjugation runs: minutes range')
    an = jload('cs_analysis', 'done.json') or {}
    po = [v for k, v in an.items() if k.startswith('cspoly_') and 'orth_defect_mean' in v
          and v.get('case') not in ('s2', 's1xs1')]
    if po:
        put('4.2.orth_defect_mean_range', [min(v['orth_defect_mean'] for v in po), max(v['orth_defect_mean'] for v in po)],
            'cs_analysis: mean ||J^T J - I|| over the 4D/6D Table 3 fields')
        put('4.2.singular_value_range', [min(v['sv_min'] for v in po), max(v['sv_max'] for v in po)],
            'cs_analysis: singular values of E^T J E')

    # ---- §4.3, Table 4: supervised benchmark
    for c in ('s3xs1', 's5xs1'):
        val = res(f'hopf_{c}_val_b7')
        if val:
            put(f'T4.{c}.values_only.mean', val['ambient']['nij_mean'], f'hopf_{c}_val_b7: ambient.nij_mean (4000)')
            put(f'T4.{c}.values_only.adv', val.get('adv_nij_sup'), f'hopf_{c}_val_b7: adv_nij_sup')
            put(f'T4.{c}.values_only.rel_val_rms', val['ambient']['rel_val_rms'], f'hopf_{c}_val_b7: rel. value error')
        runs = sorted(n for n in os.listdir(mf.RUNS) if n.startswith(f'hopf_{c}_sob_')) if os.path.isdir(mf.RUNS) else []
        vals = [(n, res(n)) for n in runs]
        vals = [(n, r) for n, r in vals if r]
        if vals:
            put(f'T4.{c}.values_derivs.means', {n: r['ambient']['nij_mean'] for n, r in vals}, 'ambient.nij_mean per run')
            put(f'T4.{c}.values_derivs.rel_val_rms', {n: r['ambient']['rel_val_rms'] for n, r in vals}, 'per run')
            put(f'T4.{c}.values_derivs.base_hopf_invariant', {n: mf.hopf_class_of(r) for n, r in vals},
                'Hopf invariant of the base relative to the target (hopf_class/done.json; 0 = reachable)')
            inc = [r['ambient']['nij_mean'] for n, r in vals if mf.hopf_class_of(r) in (0, None)]
            if inc:
                put(f'T4.{c}.values_derivs.range_in_class', [min(inc), max(inc)], 'bases in the target class')
        for kind in ('exp', 'aff'):
            r = res(f'csconj{kind}_{c}')
            if r:
                put(f'T4.{c}.unsup_{kind}.mean', r['dense']['nij_mean'], f'csconj{kind}_{c}: dense.nij_mean')

    # ---- §4.4, Table 5, Figure 3: hard cases
    for c in ('s6', 's4xs2', 's4xs1xs1'):
        b = res(f'acs_{c}_s7')
        if b and b.get('corrected'):
            put(f'T5.{c}.acs_start.mean', b['corrected']['nij_mean'], f'acs_{c}_s7: corrected.nij_mean')
            put(f'T5.{c}.acs_start.adv', b.get('adv_nij_sup'), f'acs_{c}_s7: adv_nij_sup')
        for name, lab in ((f'csconjexp_{c}', 'exp'), (f'csconjaff_{c}', 'affine'), (f'csconjtail_{c}', 'tail'),
                          (f'csconjtail_{c}_cont', 'tail_cont')):
            r = res(name)
            if r:
                put(f'T5.{c}.{lab}.mean', r['dense']['nij_mean'], f'{name}: dense.nij_mean')
                put(f'T5.{c}.{lab}.max', r['dense']['nij_sup'], f'{name}: dense.nij_sup (6000)')
                put(f'T5.{c}.{lab}.adv', r.get('adv_nij_sup'), f'{name}: adv_nij_sup')
                if r.get('conjugator'):
                    put(f'T5.{c}.{lab}.g_sigma_min_adv', r['conjugator'].get('g_sigma_min_adv'), f'{name}: conjugator')
                    put(f'T5.{c}.{lab}.g_sigma_min_dense', r['conjugator'].get('g_sigma_min_dense'), f'{name}: conjugator')
        r = res(f'conjpoly_{c}')
        if r:
            put(f'T5.{c}.polygen.degree', r.get('degree'), f'conjpoly_{c}: retained degree')
            put(f'T5.{c}.polygen.mean', r['dense']['nij_mean'], f'conjpoly_{c}: dense.nij_mean')
            put(f'T5.{c}.polygen.adv', r.get('adv_nij_sup'), f'conjpoly_{c}: adv_nij_sup')
            put(f'T5.{c}.polygen.ladder', [(e['degree'], e['N_mean'], e['adv_N_sup']) for e in mf.conjpoly_table(f'conjpoly_{c}')],
                'per degree: 2000-point mean, 64 x 60 ascent')
        r = res(f'cspoly_{c}_s7')
        if r:
            L = [e for e in r.get('ladder', []) if e.get('phase') == 'poly_cs']
            put(f'4.4.cspoly.{c}.returned_degree', r.get('degree'), f'cspoly_{c}_s7: degree of the returned field')
            put(f'4.4.cspoly.{c}.loss_per_degree', {str(e['degree']): e['final_N2'] + 10 * e['final_C2'] for e in L},
                'held-out loss per degree')
            put(f'4.4.cspoly.{c}.r_adv', r.get('adv_acs_sup'), f'cspoly_{c}_s7: adv_acs_sup of the returned field')
    for name, q in (('csconjexp_s6', '0.05'), ('csconjexp_s4xs2', '0.05'), ('csconjexp_s5xs1', '0.01'),
                    ('csconjaff_s5xs1', '0.01'), ('csconjaff_s4xs2', '0.01'), ('csconjtail_s6', '0.05'),
                    ('csconjtail_s4xs2', '0.05'), ('csconjtail_s4xs1xs1_cont', '0.05')):
        v = an.get(name, {}).get('share_top', {})
        if q in v:
            put(f'4.4.share.{name}.worst{q}', v[q], f'cs_analysis: share of E||N||^2 in the worst {q}')

    # ---- §5: reconstruction
    cert = jload('reconstruction', 'certificate.json')
    if cert:
        for m in ('S6', 'S2xS4'):
            c = cert[m]
            put(f'T6.{m}.acs_residual_max', c.get('acs_residual_max'), f'reconstruction: {m}.acs_residual_max')
            put(f'T6.{m}.N_adversarial_sup', c.get('N_adversarial_sup'), f'reconstruction: {m}.N_adversarial_sup')
            put(f'T6.{m}.N_disc_max', max(c.get('N_disc_Z1_max', 0), c.get('N_disc_Z2_max', 0)), 'disc maxima')
            put(f'T6.{m}.functional_residual_max', c.get('functional_residual_max'), 'certificate functional residual')
            put(f'T6.{m}.deck_max', max(v for k, v in c.items() if k.startswith('deck_') and isinstance(v, float)),
                'largest deck discrepancy')
            put(f'T6.{m}.blend_N_max', max(v for k, v in c.items() if k.startswith('blend_') and k.endswith('_max')),
                'blended-overlap N')
            put(f'T6.{m}.D_max', c.get('D_max_on_domain'), 'max of D on the domain (must be < 0)')
    st = jload('reconstruction_K', 'done.json')
    if st:
        put('5.2.q0', st.get('qmax'), 'reconstruction_K: qmax (basis scale q0)')
        for rec in st['study']:
            put(f"F4.{rec['collocation']}.K{rec['K']}", max(rec['held_out'].values()), 'held-out functional residual')
    s4 = jload('s4t2_reconstruction', 'certificate.json')
    if s4:
        for nm in ('paper', 'learned_best'):
            put(f'T6.S4xT2.{nm}.det_normalised_min', s4[nm]['det_normalised_min'], f'{nm}: smallest normalised det')
            put(f'T6.S4xT2.{nm}.own_annulus.N_max', s4[nm]['N_max'], f'{nm}: own annulus')
            put(f'T6.S4xT2.{nm}.own_annulus.chart_inf', s4[nm]['chart_infinity_consistency_max'], f'{nm}: own annulus')
            put(f'T6.S4xT2.{nm}.filling_eps', [s4[nm]['filling_eps_0'], s4[nm]['filling_eps_inf']], 'filling radii')
            ca = (s4.get('common_annulus') or {}).get(nm)
            if ca:
                put(f'T6.S4xT2.{nm}.common_annulus.N_max', ca['N_max'], f"common annulus eps {ca['annulus_eps']:.3g}")
                put(f'T6.S4xT2.{nm}.common_annulus.chart_inf', ca['chart_infinity_consistency_max'], 'common annulus')
        ls = s4['learned_summary']
        put('5.4.starts_admissible', [ls['admissible'], ls['starts']], 'admissible starts / starts')
        put('5.4.margins', sorted(set(round(m, 4) for m in ls['margins'])), 'final feasible margins')

    # ---- Appendix B
    oc = jload('offline_checks', 'done.json')
    if oc:
        put('B.offline_checks', oc, 'offline_checks: done.json')
    hc = jload('hopf_class', 'done.json')
    if hc:
        put('4.3.hopf_invariant', {n: u.get('hopf_invariant') for n, u in hc['units'].items()},
            'hopf_class: Hopf invariant relative to the Hopf target, S^3 x S^1 fields')
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--runs', default=None, help='run directory (default $PAPER_OUT/runs)')
    ap.add_argument('--json', default=None)
    a = ap.parse_args()
    if a.runs:
        os.environ['PAPER_RUNS'] = os.path.abspath(a.runs)
    import make_figures as mf
    out = collect(mf)
    for k, v in out.items():
        val = v['value']
        if isinstance(val, float):
            val = f'{val:.4g}'
        print(f'{k:48s} {val}   [{v["source"]}]')
    if a.json:
        json.dump(out, open(a.json, 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
