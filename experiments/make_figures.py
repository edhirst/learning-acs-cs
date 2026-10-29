"""Figures (PDF) and LaTeX tables of the report, from the rerun data.

    .venv/bin/python experiments/make_figures.py [--only fig_acs_ladder,tab_summary]

Reads $PAPER_OUT/runs/<job>/{result.json, progress.jsonl, done.json,
certificate.json} (default PAPER_OUT = outputs/paper_rerun); writes
$PAPER_OUT/figures/*.pdf, $PAPER_OUT/tables/*.tex and $PAPER_OUT/tables/numbers.tex
(macros for every number quoted in the text).
Missing runs are skipped, so the script can be run while the campaign is
still going. Colours: the validated categorical palette (fixed order) and the
one-hue blue ramp of the dataviz reference instance.
"""
from __future__ import annotations

import os
import re
import sys
import json
import math
import argparse
import traceback

import numpy as np
import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
OUT = os.path.abspath(os.environ.get('PAPER_OUT', os.path.join(ROOT, 'outputs', 'paper_rerun')))
RUNS = os.environ.get('PAPER_RUNS', os.path.join(OUT, 'runs'))
FIG = os.environ.get('PAPER_FIG_DIR', os.path.join(OUT, 'figures'))
TAB = os.environ.get('PAPER_TAB_DIR', os.path.join(OUT, 'tables'))

# ---- palette (validated: scripts/validate_palette.py, light mode, adjacent pairs;
#      first three slots also pass all-pairs, used for scatter plots)
PAL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = PAL
RAMP = ['#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7', '#3987e5', '#2a78d6',
        '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b']          # sequential blue 100 -> 700
ORDINAL = ['#86b6ef', '#5598e7', '#2a78d6', '#1c5cab', '#104281']     # validated --ordinal
INK, INK2, MUTED, GRID = '#0b0b0b', '#52514e', '#898781', '#e1e0d9'
SEQ = LinearSegmentedColormap.from_list('blue_seq', RAMP)

CASES = ['s2', 's1xs1', 's4', 's2xs2', 's3xs1', 's6', 's3xs3', 's5xs1', 's4xs2', 's4xs1xs1']
CLASSICAL = ['s2xs2', 's3xs1', 's3xs3', 's5xs1']
HARD = ['s6', 's4xs2', 's4xs1xs1']
TEX = {'s2': r'$S^2$', 's4': r'$S^4$', 's6': r'$S^6$', 's1xs1': r'$S^1\times S^1$', 's2xs2': r'$S^2\times S^2$',
       's3xs1': r'$S^3\times S^1$', 's3xs3': r'$S^3\times S^3$', 's4xs2': r'$S^4\times S^2$',
       's5xs1': r'$S^5\times S^1$', 's4xs1xs1': r'$S^4\times S^1\times S^1$'}
DIM = {'s2': 2, 's1xs1': 2, 's4': 4, 's2xs2': 4, 's3xs1': 4, 's6': 6, 's3xs3': 6, 's5xs1': 6, 's4xs2': 6,
       's4xs1xs1': 6}
NFAC = {'s2': 1, 's4': 1, 's6': 1, 's4xs1xs1': 3}
LADDER = ['c', 'lin', 'bilin', 'quad_h', 'lin+c', 'lin+bilin', 'all2', 'all2+trilin']
LADDER_SHORT = {'c': 'c', 'lin': 'lin', 'bilin': 'bilin', 'quad_h': 'quad', 'lin+c': 'lin+c',
                'lin+bilin': 'lin+bi', 'all2': 'all2', 'all2+trilin': '+trilin'}
ATTEMPT_MARK = {'plain': 'o', 'dirichlet': 's', 'orient': '^'}


def style():
    mpl.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic',
        'mathtext.bf': 'Arial:bold', 'mathtext.sf': 'Arial', 'mathtext.fallback': 'cm',
        'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
        'legend.fontsize': 7, 'axes.titlelocation': 'left', 'axes.titlecolor': INK, 'axes.titlepad': 5,
        'axes.edgecolor': MUTED, 'axes.linewidth': 0.6, 'axes.labelcolor': INK,
        'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'axes.axisbelow': True,
        'grid.color': GRID, 'grid.linewidth': 0.5, 'grid.linestyle': '-',
        'xtick.color': MUTED, 'ytick.color': MUTED, 'xtick.labelcolor': INK2, 'ytick.labelcolor': INK2,
        'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.minor.width': 0.4, 'ytick.minor.width': 0.4,
        'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.minor.size': 1.5, 'ytick.minor.size': 1.5,
        'lines.linewidth': 1.4, 'lines.solid_capstyle': 'round', 'lines.solid_joinstyle': 'round',
        'lines.dash_capstyle': 'round', 'lines.markersize': 4.5,
        'legend.frameon': False, 'legend.handlelength': 1.6, 'legend.borderaxespad': 0.3, 'legend.labelcolor': INK,
        'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
        'pdf.fonttype': 42, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.03,
    })


# ---------------------------------------------------------------------------
# data access
# ---------------------------------------------------------------------------

def jload(*parts):
    p = os.path.join(RUNS, *parts)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def res(name):
    return jload(name, 'result.json')


def prog(name):
    p = os.path.join(RUNS, name, 'progress.jsonl')
    if not os.path.exists(p):
        return []
    out = []
    with open(p) as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    out.append(json.loads(line))
                except json.JSONDecodeError:
                    pass
    return out


def corrected(r):
    """(mean |N|, max |N|) of the exact field reported by a CS unit."""
    if r is None:
        return None, None
    c = r.get('corrected') if not r.get('exact_acs_by_construction') else r.get('dense')
    c = c or r.get('dense') or {}
    return c.get('nij_mean'), c.get('nij_sup')


def best_unbiased(case):
    """Best (smallest adversarial sup |N|) unbiased CS run on a case; (name, result)."""
    best = None
    for name in sorted(os.listdir(RUNS)) if os.path.isdir(RUNS) else []:
        if not name.startswith(('cspoly_', 'csconj', 'conjpoly_', 'lp_', 'polar_', 'atlas_')):
            continue
        r = res(name)
        if r is None or r.get('case') != case or r.get('status') == 'no_certified_base':
            continue
        a = r.get('adv_nij_sup')
        if a is None or not np.isfinite(a):
            continue
        if r.get('cs_mode', 'conj') in ('poly', 'polar', 'atlas') and not r.get('certified_exact_acs', True):
            continue
        if best is None or a < best[1]['adv_nij_sup']:
            best = (name, r)
    return best


# ---------------------------------------------------------------------------
# formatting
# ---------------------------------------------------------------------------

def sci(x, d=1, tex=True):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return '--'
    x = float(x)
    if x == 0:
        return '$0$'
    e = int(math.floor(math.log10(abs(x))))
    if -2 <= e <= 2:
        s = f'{x:.{max(0, d + 1 - e - 1) if e < 0 else max(0, d - e) if e < d else 0}f}'
        if e >= 0:
            s = f'{x:.{max(0, d - e)}f}' if abs(x) < 100 else f'{x:.0f}'
        else:
            s = f'{x:.{d - e}f}'
        return f'${s}$' if tex else s
    m = x / 10 ** e
    if round(m, d) >= 10:
        m, e = m / 10, e + 1
    s = f'{m:.{d}f}\\times10^{{{e}}}'
    return f'${s}$' if tex else s


def mt(x, d=1):
    """Matplotlib mathtext version of sci()."""
    return sci(x, d).replace('\\times', r'\times')


def macro(name, value):
    v = str(value)
    if len(v) > 1 and v.startswith('$') and v.endswith('$') and v.count('$') == 2:
        v = '\\ensuremath{' + v[1:-1] + '}'          # usable in and out of math mode
    return f'\\newcommand{{\\{name}}}{{{v}}}\n'


NUMBERS = {}


def num(name, value):
    if not re.fullmatch(r'[A-Za-z]+', name):
        raise ValueError(name)
    NUMBERS[name] = value


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, name + '.pdf'))
    if os.environ.get('FIG_PNG_DIR'):                     # previews for visual inspection
        fig.savefig(os.path.join(os.environ['FIG_PNG_DIR'], name + '.png'), dpi=220)
    plt.close(fig)
    print('wrote', name, flush=True)


def write_table(name, body):
    os.makedirs(TAB, exist_ok=True)
    with open(os.path.join(TAB, name + '.tex'), 'w') as f:
        f.write(body)
    print('wrote table', name, flush=True)


def logaxis(ax, axis='y'):
    (ax.set_yscale if axis == 'y' else ax.set_xscale)('log')
    loc = mpl.ticker.LogLocator(base=10, numticks=12)
    (ax.yaxis if axis == 'y' else ax.xaxis).set_major_locator(loc)
    (ax.yaxis if axis == 'y' else ax.xaxis).set_minor_locator(mpl.ticker.NullLocator())


def decade_ticks(ax, lo, hi, axis='y', every=1):
    ticks = [10.0 ** k for k in range(lo, hi + 1, every)]
    (ax.set_yticks if axis == 'y' else ax.set_xticks)(ticks)
    (ax.set_yticklabels if axis == 'y' else ax.set_xticklabels)([f'$10^{{{k}}}$' for k in range(lo, hi + 1, every)])


def hline(ax, y, text, color=MUTED, where='right', va='bottom', dy=1.15, fs=6.5):
    ax.axhline(y, color=color, lw=0.8, zorder=1)
    x = ax.get_xlim()[1 if where == 'right' else 0]
    ax.text(x, y * (dy if va == 'bottom' else 1 / dy), text, color=INK2, fontsize=fs,
            ha='right' if where == 'right' else 'left', va=va)


def vline(ax, x, text, color=MUTED, fs=6.5, y=None, ha='left'):
    ax.axvline(x, color=color, lw=0.8, zorder=1)
    if text:
        yy = ax.get_ylim()[1] if y is None else y
        ax.text(x, yy, ' ' + text, color=INK2, fontsize=fs, ha=ha, va='top')


def mark(ax, x, y, color, marker='o', ms=4.5, **kw):
    ax.plot(x, y, marker=marker, ms=ms, color=color, mec='white', mew=0.7, ls='none', **kw)


# ---------------------------------------------------------------------------
# Figure 1: overview
# ---------------------------------------------------------------------------

def fig_overview():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.3, 3.4), sharey=True, gridspec_kw={'wspace': 0.08})
    ys = np.arange(len(CASES))[::-1]
    for y, c in zip(ys, CASES):
        for s, col, off in ((7, BLUE, 0.12), (8, ORANGE, -0.12)):
            r = res(f'acs_{c}_s{s}')
            if r and r.get('adv_acs_sup') is not None:
                mark(a1, r['adv_acs_sup'], y + off, col)
            if r and r.get('corrected'):
                mark(a1, r['corrected']['adv_acs_sup'], y + off, AQUA, marker='D', ms=3.8)
    a1.set_xscale('log'); a1.set_xlim(3e-11, 12)
    decade_ticks(a1, -10, 1, 'x', every=2)
    a1.axvline(1.0, color=MUTED, lw=0.8, zorder=1)
    a1.text(0.85, -0.62, 'degeneracy bound', fontsize=6.0, color=INK2, va='bottom', ha='right',
            bbox=dict(facecolor='white', edgecolor='none', pad=0.4))
    a1.set_ylim(-0.95, len(CASES) - 0.3)
    a1.set_yticks(ys); a1.set_yticklabels([TEX[c] for c in CASES])
    a1.set_xlabel(r'adversarial sup of $\|K^2+P\|$')
    a1.set_title('Almost complex structure (ACS)')
    a1.grid(axis='y', visible=False)
    h = [Line2D([], [], marker='o', ls='none', color=BLUE, mec='white', label='raw field, seed 7'),
         Line2D([], [], marker='o', ls='none', color=ORANGE, mec='white', label='raw field, seed 8'),
         Line2D([], [], marker='D', ls='none', color=AQUA, mec='white', ms=3.8, label='after polar correction')]
    a1.legend(handles=h, loc='upper left', bbox_to_anchor=(-0.02, -0.2), ncol=2, fontsize=6.3)
    cert = jload('reconstruction', 'certificate.json')
    for y, c in zip(ys, CASES):
        if c == 's4':
            a2.text(2e-12, y, 'refused: no certified ACS (none exists)', fontsize=6.5, color=INK2, va='center',
                    bbox=dict(facecolor='white', edgecolor='none', pad=0.6))
            continue
        b = best_unbiased(c)
        if b is not None:
            mark(a2, b[1]['adv_nij_sup'], y + 0.1, BLUE)
        if cert and c in ('s6', 's4xs2'):
            key = 'S6' if c == 's6' else 'S2xS4'
            mark(a2, cert[key]['certified_sup'], y - 0.1, ORANGE, marker='s', ms=4.2)
        c4 = jload('s4t2_reconstruction', 'certificate.json')
        if c4 and c == 's4xs1xs1':
            mark(a2, c4['paper']['certified_sup'], y - 0.1, ORANGE, marker='s', ms=4.2)
        inv = [jload(n, 'done.json') for n in ('s4t2inv_s7', 's4t2inv_s8')]
        advs = [e['N_adv'] for d in inv if d for e in d.get('cs', []) if 'N_adv' in e]
        if advs and c == 's4xs1xs1':
            mark(a2, min(advs), y + 0.3, AQUA, marker='D', ms=3.8)
    a2.set_xscale('log'); a2.set_xlim(1e-12, 1e3)
    decade_ticks(a2, -12, 3, 'x', every=3)
    a2.axvline(1e-6, color=MUTED, lw=0.8, zorder=1)
    a2.text(1.4e-6, -0.55, 'target $10^{-6}$', fontsize=6.0, color=INK2, va='bottom')
    a2.set_xlabel(r'adversarial sup of $\|N_J\|$ on an exact ACS')
    a2.set_title('Complex structure (CS)')
    a2.grid(axis='y', visible=False)
    h = [Line2D([], [], marker='o', ls='none', color=BLUE, mec='white', label='unsupervised, unbiased PINN (best run)'),
         Line2D([], [], marker='s', ls='none', color=ORANGE, mec='white', ms=4.2,
                label='numerical reconstruction (fibration model)')]
    if jload('s4t2inv_s7', 'done.json'):
        h.append(Line2D([], [], marker='D', ls='none', color=AQUA, mec='white', ms=3.8,
                        label=r'unsupervised, $T^2$-invariant ansatz'))
    a2.legend(handles=h, loc='upper left', bbox_to_anchor=(-0.02, -0.2), ncol=1, fontsize=6.3)
    save(fig, 'fig_overview')


# ---------------------------------------------------------------------------
# ACS figures
# ---------------------------------------------------------------------------

def ladder_entries(r):
    return [e for e in (r or {}).get('ladder', []) if e.get('phase') == 'poly']


def fig_acs_ladder():
    fig, axes = plt.subplots(2, 5, figsize=(6.3, 3.5), sharey=True, gridspec_kw={'hspace': 0.75, 'wspace': 0.1})
    for ax, c in zip(axes.flat, CASES):
        nf = NFAC.get(c, 2)
        classes = [f for f in LADDER if ('bilin' not in f or nf >= 2) and ('trilin' not in f or nf >= 3)]
        xi = {f: i for i, f in enumerate(classes)}
        for s, col, off in ((7, BLUE, -0.13), (8, ORANGE, 0.13)):
            r = res(f'acs_{c}_s{s}')
            for e in ladder_entries(r):
                k = e.get('attempt') or ('dirichlet' if e.get('lam0', 0) > 0 else 'plain')
                mark(ax, xi[e['feat']] + off, e['adv_acs_sup'], col, marker=ATTEMPT_MARK.get(k, 'o'), ms=4)
        # down to 1e-11: since the polynomial attempts are purely polynomial (2026-10-25), an exact
        # rung (S^3 x S^1 seed 7, the Dirichlet rungs of S^2 / S^1 x S^1) reaches the 1e-10 floor
        ax.set_yscale('log'); ax.set_ylim(3e-11, 6)
        ax.axhline(1.0, color=MUTED, lw=0.8, zorder=1)
        ax.set_xticks(range(len(classes)))
        ax.set_xticklabels([LADDER_SHORT[f] for f in classes], rotation=90, fontsize=6)
        ax.set_xlim(-0.6, len(classes) - 0.4)
        ax.set_title(TEX[c], fontsize=7.5)
        ax.grid(axis='x', visible=False)
    decade_ticks(axes[0, 0], -10, 0, 'y', every=2)
    for ax in axes[:, 0]:
        ax.set_ylabel(r'adv. sup $\|K^2+P\|$', fontsize=7)
    h = [Line2D([], [], marker='o', ls='none', color=BLUE, mec='white', label='seed 7'),
         Line2D([], [], marker='o', ls='none', color=ORANGE, mec='white', label='seed 8'),
         Line2D([], [], marker='o', ls='none', color=MUTED, mec='white', label='plain'),
         Line2D([], [], marker='s', ls='none', color=MUTED, mec='white', label='Dirichlet homotopy'),
         Line2D([], [], marker='^', ls='none', color=MUTED, mec='white', label='orientation penalty')]
    fig.legend(handles=h, loc='upper center', ncol=5, bbox_to_anchor=(0.5, 1.03), fontsize=6.5)
    save(fig, 'fig_acs_ladder')


def fig_acs_refine():
    fig, axes = plt.subplots(2, 5, figsize=(6.3, 3.2), sharey=True, sharex=True,
                             gridspec_kw={'hspace': 0.45, 'wspace': 0.1})
    for ax, c in zip(axes.flat, CASES):
        r = res(f'acs_{c}_s7')
        p = [q for q in prog(f'acs_{c}_s7') if 'step' in q]
        if p:
            st = [q['step'] for q in p]
            ax.plot(st, [q['acs_mean'] for q in p], color=BLUE, lw=1.2)
            ax.plot(st, [q['acs_sup'] for q in p], color=ORANGE, lw=1.2)
        if r and r.get('dense'):
            x_end = 10800
            mark(ax, x_end, r['dense']['acs_mean'], BLUE, marker='o', ms=3.5)
            mark(ax, x_end, r['dense']['acs_sup'], ORANGE, marker='o', ms=3.5)
            if r.get('corrected'):
                ax.text(0.97, 0.06, 'polar: ' + mt(r['corrected']['acs_sup']), transform=ax.transAxes,
                        fontsize=5.8, color=INK2, ha='right')
            else:
                ax.text(0.97, 0.06, 'not admissible', transform=ax.transAxes, fontsize=5.8, color=INK2, ha='right')
        ax.set_yscale('log'); ax.set_ylim(1e-6, 5)
        ax.axhline(1.0, color=MUTED, lw=0.8, zorder=1)
        ax.set_title(TEX[c], fontsize=7.5)
        ax.set_xticks([0, 5000, 10000]); ax.set_xticklabels(['0', '5k', '10k'])
    decade_ticks(axes[0, 0], -6, 0, 'y', every=2)
    for ax in axes[1]:
        ax.set_xlabel('Adam step', fontsize=7)
    for ax in axes[:, 0]:
        ax.set_ylabel(r'$\|K^2+P\|$ held out', fontsize=7)
    h = [Line2D([], [], color=BLUE, label='mean over 1500 fresh points'),
         Line2D([], [], color=ORANGE, label='max over 1500 fresh points'),
         Line2D([], [], marker='o', ls='none', color=MUTED, mec='white', ms=3.5, label='after L-BFGS (6000 points)')]
    fig.legend(handles=h, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.02), fontsize=6.5)
    save(fig, 'fig_acs_refine')


def fig_s4():
    fig, ax = plt.subplots(figsize=(3.3, 2.6))
    for grp, col, lab in ((['s4'], ORANGE, r'$S^4$ (no ACS exists)'),
                          ([c for c in CASES if c != 's4'], BLUE, 'the nine other cases')):
        xs, ys = [], []
        for c in grp:
            for s in (7, 8):
                r = res(f'acs_{c}_s{s}')
                for e in ladder_entries(r):
                    xs.append(e['acs_mean']); ys.append(e['adv_acs_sup'])
                if r and r.get('dense'):
                    xs.append(r['dense']['acs_mean']); ys.append(r['adv_acs_sup'])
        mark(ax, xs, ys, col, ms=4, label=lab, alpha=0.9)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlim(1e-6, 5); ax.set_ylim(1e-5, 6)
    lim = np.array([1e-6, 5])
    ax.plot(lim, lim, color=GRID, lw=0.8, zorder=0)
    ax.axhline(1.0, color=MUTED, lw=0.8, zorder=1)
    ax.text(1.3e-6, 1.25, r'degeneracy bound: $\|K^2+P\|\geq 1$ at a real eigenvalue', fontsize=6, color=INK2)
    decade_ticks(ax, -6, 0, 'x', every=2); decade_ticks(ax, -5, 0, 'y')
    ax.set_xlabel(r'held-out mean of $\|K^2+P\|$')
    ax.set_ylabel(r'adversarial sup of $\|K^2+P\|$')
    ax.legend(loc='lower right')
    save(fig, 'fig_s4')


def fig_orientation():
    d = jload('offline_checks', 'done.json')
    if not d:
        return
    o = d['orientation']
    cases = [c for c in CASES if c in o and 'error' not in o[c]]
    fig, ax = plt.subplots(figsize=(3.3, 2.4))
    for y, c in enumerate(cases[::-1]):
        e = np.array(o[c]['edges']); h = np.array(o[c]['hist'])
        sgn = 1.0 if o[c]['o_max'] > 0 else -1.0
        lo, hi = sorted([sgn * o[c]['o_min'], sgn * o[c]['o_max']])
        mids = sgn * 0.5 * (e[1:] + e[:-1])
        cdf = np.cumsum(h) / h.sum()
        order = np.argsort(mids)
        med = np.interp(0.5, np.cumsum(h[order]) / h.sum(), mids[order])
        ax.plot([lo, hi], [y, y], color=BLUE, lw=3.0, solid_capstyle='round')
        mark(ax, med, y, BLUE, ms=5)
        ax.text(hi + 0.01, y, f'{lo:.3f}–{hi:.3f}', fontsize=6, color=INK2, va='center')
    ax.axvline(1.0, color=MUTED, lw=0.8)
    ax.text(1.0, len(cases) - 0.35, r' $|o|\geq 1$ (Lemma 4.3)', fontsize=6, color=INK2, va='bottom')
    ax.set_yticks(range(len(cases))); ax.set_yticklabels([TEX[c] for c in cases[::-1]])
    ax.set_ylim(-0.6, len(cases) - 0.1)
    ax.set_xlabel(r'orientation index $s\,o(p)$ of the certified exact ACS (range, median)')
    ax.grid(axis='y', visible=False)
    xmax = max(max(abs(o[c]['o_min']), abs(o[c]['o_max'])) for c in cases)
    ax.set_xlim(0.95, xmax + 0.12 * (xmax - 0.95) + 0.08)
    save(fig, 'fig_orientation')


# ---------------------------------------------------------------------------
# CS figures
# ---------------------------------------------------------------------------

def poly_curve(name, w=10.0):
    """(minutes, held-out loss, degree) along the Adam warm-up and the LM rounds."""
    out = []
    for q in prog(name):
        if q.get('phase') in ('poly_adam', 'poly_lm'):
            out.append((q['elapsed_min'], q['N2'] + w * q['C2'], str(q['degree']), q['phase']))
    return out


def fig_cs_poly():
    cases = CLASSICAL + HARD
    fig, axes = plt.subplots(2, 4, figsize=(6.3, 3.4), gridspec_kw={'hspace': 0.55, 'wspace': 0.32})
    for ax, c in zip(axes.flat, cases):
        for s, col in ((7, BLUE), (8, ORANGE)):
            cv = poly_curve(f'cspoly_{c}_s{s}')
            if not cv:
                continue
            t = [v[0] for v in cv]; L = [max(v[1], 1e-30) for v in cv]
            ax.plot(t, L, color=col, lw=1.2, marker='o', ms=2.6, mec='white', mew=0.4)
            for i in range(1, len(cv)):
                if cv[i][2] != cv[i - 1][2]:
                    ax.axvline(t[i], color=GRID, lw=0.8, zorder=0)
                    ax.text(t[i], 1e-2 if c in CLASSICAL else 2e3, f' deg {cv[i][2]}', fontsize=5.8, color=INK2, va='top')
        ax.set_yscale('log')
        if c in CLASSICAL:
            ax.set_ylim(1e-26, 1e4); decade_ticks(ax, -24, 0, 'y', every=6)
            ax.axhline(1e-20, color=MUTED, lw=0.8, zorder=1)
        else:
            ax.set_ylim(0.3, 3e4); decade_ticks(ax, 0, 4, 'y', every=1)
        ax.set_title(TEX[c], fontsize=7.5)
        ax.set_xlabel('wall clock [min]', fontsize=6.5)
    axes.flat[-1].axis('off')
    axes[0, 0].set_ylabel(r'held-out $\mathbb{E}\|N_K\|^2+10\,\mathbb{E}\|K^2+P\|^2$', fontsize=6.5)
    axes[1, 0].set_ylabel(r'held-out loss', fontsize=6.5)
    axes[0, 0].text(0.02, 1e-20 * 3, 'acceptance $10^{-20}$', fontsize=5.8, color=INK2,
                    transform=mpl.transforms.blended_transform_factory(axes[0, 0].transAxes, axes[0, 0].transData))
    h = [Line2D([], [], color=BLUE, marker='o', ms=2.6, mec='white', label='seed 7'),
         Line2D([], [], color=ORANGE, marker='o', ms=2.6, mec='white', label='seed 8')]
    axes.flat[-1].legend(handles=h, loc='center left', fontsize=7,
                         title='Adam warm-up, then one point\nper Levenberg–Marquardt round',
                         title_fontsize=6.5, alignment='left')
    save(fig, 'fig_cs_poly')


def conj_curve(name):
    return [(q['elapsed_min'], q['nij_mean'], q['nij_sup'], q['step']) for q in prog(name) if 'step' in q]


def fig_cs_conj_vs_poly():
    fig, axes = plt.subplots(1, 4, figsize=(6.3, 2.2), sharey=True, gridspec_kw={'wspace': 0.12})
    for ax, c in zip(axes, CLASSICAL):
        for name, col, lab in ((f'csconjexp_{c}', BLUE, 'conjugation, exp, Adam + L-BFGS'),
                               (f'csconjaff_{c}', ORANGE, 'conjugation, affine, Adam + L-BFGS')):
            cv = conj_curve(name)
            if cv:
                ax.plot([v[0] for v in cv], [v[1] for v in cv], color=col, lw=1.2, label=lab)
                r = res(name)
                if r and r.get('dense'):
                    mark(ax, r['elapsed_min'], r['dense']['nij_mean'], col, ms=3.5)
        cv = poly_curve(f'cspoly_{c}_s7', w=0.0)
        if cv:
            ax.plot([v[0] for v in cv], [math.sqrt(max(v[1], 1e-40)) for v in cv], color=AQUA, lw=1.2,
                    label='polynomial field, Levenberg–Marquardt')
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_ylim(1e-16, 1e3); decade_ticks(ax, -15, 2, 'y', every=3)
        ax.set_xlim(0.05, 200)
        ax.set_xticks([0.1, 1, 10, 100]); ax.set_xticklabels(['0.1', '1', '10', '100'])
        ax.set_title(TEX[c], fontsize=7.5)
        ax.set_xlabel('wall clock [min]', fontsize=6.5)
    axes[0].set_ylabel(r'held-out $\|N_J\|$ (mean; RMS for LM)', fontsize=6.5)
    hs, ls = axes[0].get_legend_handles_labels()
    if not hs:
        hs, ls = axes[1].get_legend_handles_labels()
    fig.legend(hs, ls, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.08), fontsize=6.5)
    save(fig, 'fig_cs_conj_vs_poly')


def fig_hard_conj():
    fig, axes = plt.subplots(2, 3, figsize=(6.3, 3.4), sharex='col', gridspec_kw={'hspace': 0.18, 'wspace': 0.25})
    runs = (('csconjexp_{}', BLUE, 'exp conjugation, mean-square'),
            ('csconjaff_{}', ORANGE, 'affine conjugation, mean-square'),
            ('csconjtail_{}', AQUA, 'exp conjugation, tail-aware sampling'))
    for j, c in enumerate(HARD):
        for pat, col, lab in runs:
            name = pat.format(c)
            cv = conj_curve(name)
            if not cv:
                continue
            cont = conj_curve(name + '_cont')
            if cont:                                         # warm-started continuation of the same run
                off = cv[-1][3]
                cv = cv + [(v[0], v[1], v[2], v[3] + off) for v in cont]
                if res(name + '_cont'):
                    name = name + '_cont'
                axes[0, j].axvline(off, color=GRID, lw=0.8, zorder=0)
                axes[0, j].text(off, 2e3, ' continuation', fontsize=5.8, color=INK2, va='top')
            st = [v[3] for v in cv]
            axes[0, j].plot(st, [v[1] for v in cv], color=col, lw=1.2, label=lab)
            axes[1, j].plot(st, [v[2] for v in cv], color=col, lw=1.2)
            r = res(name)
            if r and r.get('adv_nij_sup') is not None:
                a = r['adv_nij_sup']
                if a > 3e3:
                    mark(axes[1, j], st[-1], 3e3, col, marker='^', ms=5)
                    axes[1, j].text(st[-1], 1.5e3, mt(a, 0) + ' ', fontsize=5.8, color=INK2, ha='right', va='top')
                else:
                    mark(axes[1, j], st[-1], a, col, marker='D', ms=4)
        axes[0, j].set_title(TEX[c], fontsize=7.5)
        for i in range(2):
            axes[i, j].set_yscale('log')
        axes[0, j].set_ylim(0.05, 3e3); axes[1, j].set_ylim(0.3, 5e3)
        decade_ticks(axes[0, j], -1, 3, 'y'); decade_ticks(axes[1, j], 0, 3, 'y')
        axes[1, j].set_xlabel('Adam step', fontsize=6.5)
    axes[0, 0].set_ylabel(r'mean $\|N_J\|$ (held out)', fontsize=6.5)
    axes[1, 0].set_ylabel(r'max $\|N_J\|$ (held out)', fontsize=6.5)
    hs, ls = [], []
    for ax in axes[0]:
        for h_, l_ in zip(*ax.get_legend_handles_labels()):
            if l_ not in ls:
                hs.append(h_); ls.append(l_)
    hs.append(Line2D([], [], marker='D', ls='none', color=MUTED, mec='white', ms=4)); ls.append('final adversarial sup')
    fig.legend(hs, ls, loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.06), fontsize=6.5)
    save(fig, 'fig_hard_conj')


def conjpoly_table(name, final=False):
    """Ladder entries of a polynomial-generator unit (per degree: 2000-point mean and max,
    64 x 60 ascent). With `final`, the retained degree carries the unit's final report
    instead (6000 dense points, 128 x 75 ascent), as quoted in the paper's Table 5."""
    r = res(name)
    if not r:
        return []
    L = [dict(e) for e in r.get('ladder', []) if e.get('phase') == 'conjpoly_ladder']
    if final and r.get('dense'):
        for e in L:
            if str(e.get('degree')) == str(r.get('degree')):
                e.update(N_mean=r['dense']['nij_mean'], N_sup_dense=r['dense']['nij_sup'],
                         adv_N_sup=r.get('adv_nij_sup'), final_report=True)
    return L


def fig_conjpoly():
    fig, axes = plt.subplots(1, 2, figsize=(6.3, 2.4), sharey=True, gridspec_kw={'wspace': 0.1})
    for ax, c in zip(axes, ('s6', 's4xs2')):
        base = res(f'acs_{c}_s7')
        L = conjpoly_table(f'conjpoly_{c}')
        if not base or not L:
            continue
        deg = [0] + [e['degree'] for e in L]
        mean = [base['corrected']['nij_mean']] + [e['N_mean'] for e in L]
        dmax = [base['corrected']['nij_sup']] + [e['N_sup_dense'] for e in L]
        adv = [base['adv_nij_sup']] + [e['adv_N_sup'] for e in L]
        for ys, col, lab in ((mean, BLUE, 'mean (2000 points)'), (dmax, ORANGE, 'max (2000 points)'),
                             (adv, AQUA, 'adversarial sup')):
            ax.plot(deg, ys, color=col, lw=1.2, marker='o', ms=4, mec='white', mew=0.7, label=lab)
        if c == 's6':
            Lp = conjpoly_table('lp_s6')
            if Lp:
                e = Lp[-1]
                for v, col in ((e['N_mean'], BLUE), (e['N_sup_dense'], ORANGE), (e['adv_N_sup'], AQUA)):
                    mark(ax, e['degree'] + 0.12, v, col, marker='s', ms=4)
                ax.text(e['degree'] + 0.2, e['adv_N_sup'], r'$L^8$', fontsize=6, color=INK2, va='center')
        ax.set_yscale('log'); ax.set_ylim(0.3, 60)
        ax.set_yticks([0.3, 1, 3, 10, 30]); ax.set_yticklabels(['0.3', '1', '3', '10', '30'])
        ax.set_xticks(deg); ax.set_xticklabels(['own ACS' if k == 0 else str(k) for k in deg])
        ax.set_xlim(-0.3, max(deg) + 0.5)
        ax.set_xlabel(r'degree of the generator $B$ in $J=e^{B}J_b e^{-B}$', fontsize=6.5)
        ax.set_title(TEX[c], fontsize=7.5)
    axes[0].set_ylabel(r'$\|N_J\|$', fontsize=7)
    hs, ls = axes[0].get_legend_handles_labels()
    hs.append(Line2D([], [], marker='s', ls='none', color=MUTED, mec='white', ms=4)); ls.append(r'degree 3, $L^p$ objective (IRLS)')
    fig.legend(hs, ls, loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.07), fontsize=6.5)
    save(fig, 'fig_conjpoly')


def fig_conjpoly_rounds():
    fig, axes = plt.subplots(1, 2, figsize=(6.3, 2.3), gridspec_kw={'wspace': 0.25})
    off = 0
    for name, col, lab in (('conjpoly_s6', BLUE, r'conjugation ladder, $L^2$'),
                           ('lp_s6', ORANGE, r'degree 3, $L^p$ ramp $p: 2\to 8$')):
        recs = [q for q in prog(name) if q.get('phase') == 'conjpoly_lm']
        if not recs:
            continue
        x = np.arange(off + 1, off + len(recs) + 1)
        axes[0].plot(x, [q['N2'] for q in recs], color=col, lw=1.2, marker='o', ms=2.5, mec='white', mew=0.4, label=lab)
        axes[1].plot(x, [q.get('N_sup_sampled', np.nan) for q in recs], color=col, lw=1.2, marker='o', ms=2.5, mec='white',
                     mew=0.4, label=lab)
        degs = [q['degree'] for q in recs]
        for i in range(len(recs)):
            if i == 0 or degs[i] != degs[i - 1]:
                for ax in axes:
                    ax.axvline(x[i] - 0.5, color=GRID, lw=0.8, zorder=0)
                axes[0].text(x[i] - 0.3, 600, f'deg {degs[i]}', fontsize=5.8, color=INK2, va='top')
        off += len(recs)
    for ax in axes:
        ax.set_yscale('log'); ax.set_xlabel('Levenberg–Marquardt round', fontsize=6.5)
    axes[0].set_ylim(0.5, 800); axes[1].set_ylim(1, 100)
    axes[0].set_ylabel(r'held-out $\mathbb{E}\|N_J\|^2$', fontsize=6.5)
    axes[1].set_ylabel(r'sampled sup $\|N_J\|$ (2000 points)', fontsize=6.5)
    axes[0].set_title(r'$S^6$: mean square', fontsize=7.5); axes[1].set_title(r'$S^6$: sup', fontsize=7.5)
    axes[1].legend(loc='upper right', fontsize=6.5)
    save(fig, 'fig_conjpoly_rounds')


def fig_polar():
    recs = [q for q in prog('polar_s6') if q.get('phase') == 'polar_lm']
    if not recs:
        return
    fig, ax = plt.subplots(figsize=(3.3, 2.3))
    x = np.arange(1, len(recs) + 1)
    y = np.array([q['N2'] for q in recs], dtype=float)
    rb = np.array([q['rolled_back'] for q in recs])
    ax.plot(x, y, color=BLUE, lw=1.2, label=r'held-out $\mathbb{E}\|N(J^*)\|^2$ after the round')
    mark(ax, x[~rb], y[~rb], BLUE, ms=4)
    if rb.any():
        ax.plot(x[rb], y[rb], ls='none', marker='X', ms=5, color=ORANGE, mec='white', mew=0.5,
                label='round rolled back (admissibility lost)')
    ent = [e for e in res('polar_s6').get('ladder', [])] if res('polar_s6') else []
    if ent and ent[0].get('init_N2'):
        hline(ax, ent[0]['init_N2'], 'start: own ACS projected to degree 2', dy=1.03)
    ax.set_yscale('log')
    ax.set_xlabel('Levenberg–Marquardt round', fontsize=7)
    ax.set_title(r'$S^6$, polar mode, degree 2', fontsize=7.5)
    ax.legend(loc='lower left', fontsize=6.3)
    save(fig, 'fig_polar')


def fig_atlas():
    r = res('atlas_s2xs2')
    recs = [q for q in prog('atlas_s2xs2') if q.get('phase') == 'atlas_lm']
    if not recs:
        return
    fig, ax = plt.subplots(figsize=(3.3, 2.3))
    x = np.arange(1, len(recs) + 1)
    rb = np.array([q['rolled_back'] for q in recs])
    ax.plot(x, [q['fit_mismatch'] for q in recs], color=BLUE, lw=1.2, label='gluing mismatch, collocation pairs')
    ax.plot(x, [q['held_out'] for q in recs], color=ORANGE, lw=1.2, label='gluing mismatch, fresh pairs')
    if rb.any():
        ax.plot(x[rb], np.array([q['held_out'] for q in recs])[rb], ls='none', marker='X', ms=5, color=ORANGE,
                mec='white', mew=0.5, label='round rolled back')
    ax.set_yscale('log')
    ax.set_xlabel('Levenberg–Marquardt round (CG)', fontsize=7)
    ax.set_title(r'$S^2\times S^2$, integrable atlas, 48 charts', fontsize=7.5)
    ax.legend(loc='upper right', fontsize=6.3)
    save(fig, 'fig_atlas')


def fig_fibre():
    d = jload('fibre_topology_s6', 'done.json') or jload('fibre_topology_s6', 'log.json')
    if not d:
        return
    R = d['rounds']
    fig, axes = plt.subplots(1, 2, figsize=(6.3, 2.2), gridspec_kw={'wspace': 0.28})
    x = [q['round'] for q in R]
    axes[0].plot(x, [q['res_mean'] for q in R], color=BLUE, lw=1.2, marker='o', ms=2.5, mec='white', mew=0.4, label='mean')
    axes[0].plot(x, [q['res_max'] for q in R], color=ORANGE, lw=1.2, marker='o', ms=2.5, mec='white', mew=0.4, label='max')
    axes[1].plot(x, [q['N_mean'] for q in R], color=BLUE, lw=1.2, marker='o', ms=2.5, mec='white', mew=0.4, label='mean')
    axes[1].plot(x, [q['N_max'] for q in R], color=ORANGE, lw=1.2, marker='o', ms=2.5, mec='white', mew=0.4, label='max')
    for ax in axes:
        ax.set_yscale('log'); ax.set_xlabel('Levenberg–Marquardt round', fontsize=6.5)
    axes[0].set_ylabel(r'$\bar\partial$-integrability residual', fontsize=6.5)
    axes[1].set_ylabel(r'$\|N_J\|$ (coordinate norm)', fontsize=6.5)
    axes[0].set_title('scale-invariant period residual', fontsize=7.5)
    axes[1].set_title('Nijenhuis tensor of the fibred structure', fontsize=7.5)
    axes[1].legend(loc='upper right', fontsize=6.5)
    save(fig, 'fig_fibre')


def fig_defect():
    d = jload('cs_analysis', 'done.json')
    if not d:
        return
    fig, ax = plt.subplots(figsize=(3.3, 2.5))
    sel = [('acs_s6_s7', BLUE, r'$S^6$, own exact ACS'),
           ('csconjexp_s6', ORANGE, r'$S^6$, conjugation, mean-square'),
           ('csconjtail_s6', AQUA, r'$S^6$, conjugation, tail-aware'),
           ('csconjexp_s4xs2', YELLOW, r'$S^4\times S^2$, conjugation, mean-square'),
           ('csconjexp_s5xs1', MAGENTA, r'$S^5\times S^1$, conjugation, mean-square'),
           ('csconjtail_s4xs1xs1_cont', GREEN, r'$S^4\times S^1\times S^1$, tail-aware, 5500 steps')]
    for name, col, lab in sel:
        if name not in d or 'share_top' not in d[name]:
            continue
        st = d[name]['share_top']
        q = np.array([float(k) for k in st]); s = np.array(list(st.values()))
        o = np.argsort(q)
        ax.plot(q[o], s[o], color=col, lw=1.2, marker='o', ms=2.8, mec='white', mew=0.4, label=lab)
    qq = np.logspace(-3, 0, 50)
    ax.plot(qq, qq, color=MUTED, lw=0.8, label='uniform defect')
    ax.set_xscale('log'); ax.set_xlim(1e-3, 1); ax.set_ylim(0, 1.02)
    ax.set_xticks([1e-3, 1e-2, 1e-1, 1]); ax.set_xticklabels(['0.1%', '1%', '10%', '100%'])
    ax.set_xlabel('fraction of points (worst first)', fontsize=7)
    ax.set_ylabel(r'share of $\mathbb{E}\|N_J\|^2$', fontsize=7)
    # below the axes: inside, the legend covered the S^5 x S^1 curve, which runs along the top
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize=5.6, columnspacing=1.2)
    save(fig, 'fig_defect')


# ---------------------------------------------------------------------------
# reconstruction figures
# ---------------------------------------------------------------------------

def fig_recon_domain():
    path = os.path.join(RUNS, 'reconstruction_K', 'grids_K20.npz')
    if not os.path.exists(path):
        return
    g = np.load(path)
    X, Y, inside = g['X'], g['Y'], g['inside']
    fig, axes = plt.subplots(1, 3, figsize=(6.3, 2.35), gridspec_kw={'wspace': 0.32})
    ax = axes[0]
    ax.contourf(X, Y, inside.astype(float), levels=[0.5, 1.5], colors=['#eef4fc'])
    ax.contour(X, Y, inside.astype(float), levels=[0.5], colors=[MUTED], linewidths=0.7)
    P1, P2 = g['P1'], g['P2']
    for P, col, lab in ((P1, BLUE, r'$g_1$-equations'), (P2, ORANGE, r'$g_2$-equations')):
        inF = (P.real > X.min()) & (P.real < X.max()) & (P.imag > Y.min()) & (P.imag < Y.max())
        ax.plot(P.real[inF], P.imag[inF], ls='none', marker='o', ms=0.5, mew=0, color=col, alpha=0.45, label=lab)
    for Zc, nm in ((complex(g['Z1']), r'$Z_1$'), (complex(g['Z2']), r'$Z_2$')):
        mark(ax, Zc.real, Zc.imag, INK, ms=4.5)
        ax.text(Zc.real + 0.04, Zc.imag + 0.04, nm, fontsize=7, color=INK)
    ax.set_xlim(X.min(), X.max()); ax.set_ylim(Y.min(), Y.max())
    ax.set_aspect('equal')
    ax.set_xlabel(r'Re $Z$', fontsize=7); ax.set_ylabel(r'Im $Z$', fontsize=7)
    ax.set_title('Ford domain, collocation', fontsize=7.5)
    lg = ax.legend(loc='upper right', fontsize=5.8, markerscale=4, handletextpad=0.2)
    ax.grid(False)
    for ax, field, title, lab in ((axes[1], g['tau'].imag, r'Im $\tau(Z)$', r'Im $\tau$'),
                                  (axes[2], -g['D'], r'$-D(Z)$, $D=\mathrm{Im}\,\beta-6(\mathrm{Im}\,\mu)^2/\mathrm{Im}\,\tau$',
                                   r'$-D$')):
        F = np.where(inside, field, np.nan)
        pc = ax.pcolormesh(X, Y, F, cmap=SEQ, shading='auto', rasterized=True)
        ax.contour(X, Y, inside.astype(float), levels=[0.5], colors=[MUTED], linewidths=0.7)
        ax.set_aspect('equal'); ax.grid(False)
        ax.set_xlim(X.min(), X.max()); ax.set_ylim(Y.min(), Y.max())
        ax.set_xlabel(r'Re $Z$', fontsize=7)
        ax.set_title(title, fontsize=7.5)
        cb = fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.03)
        cb.outline.set_visible(False); cb.ax.tick_params(labelsize=6, color=MUTED, labelcolor=INK2)
    save(fig, 'fig_recon_domain')


def fig_recon_convergence():
    d = jload('reconstruction_K', 'done.json')
    if not d:
        return
    st = d['study']
    fig, axes = plt.subplots(1, 2, figsize=(6.3, 2.3), gridspec_kw={'wspace': 0.28})
    for colloc, col, lab in (('discs', BLUE, 'fixed-point discs only'),
                             ('discs+edges', ORANGE, 'discs + Ford-domain gluing pairs')):
        S = sorted([s for s in st if s['collocation'] == colloc], key=lambda s: s['K'])
        if not S:
            continue
        axes[0].plot([s['K'] for s in S], [max(s['held_out'].values()) for s in S], color=col, lw=1.2,
                     marker='o', ms=3.5, mec='white', mew=0.6, label=lab)
    axes[0].set_yscale('log'); axes[0].set_xlabel('truncation order $K$ of the $q$-series', fontsize=7)
    axes[0].set_ylabel('max functional-equation residual\non fresh gluing pairs', fontsize=6.5)
    axes[0].legend(loc='upper right', fontsize=6.3)
    S = sorted([s for s in st if s['collocation'] == 'discs+edges'], key=lambda s: s['K'])
    if S:
        top = next((r for r in S if r['K'] == 20), S[-1])       # the production truncation
        n = np.arange(len(top['coef_a']))
        for key, col, lab in (('coef_a', BLUE, r'$\tau-Z$'), ('coef_b', ORANGE, r'$\mu$'), ('coef_c', AQUA, r'$\beta+Z$')):
            v = np.maximum(np.array(top[key]), 1e-300)
            axes[1].plot(n, v, color=col, lw=1.1, marker='o', ms=2.6, mec='white', mew=0.4, label=lab)
        axes[1].set_yscale('log'); axes[1].set_ylim(1e-18, 1e2)
        decade_ticks(axes[1], -18, 2, 'y', every=4)
        axes[1].set_xlabel(r'coefficient index $n$ (basis $(q/q_0)^n$)', fontsize=7)
        axes[1].set_ylabel('coefficient modulus', fontsize=6.5)
        axes[1].set_title(f"$K={top['K']}$, discs + gluing pairs", fontsize=7.5)
        axes[1].axhline(1e-16, color=MUTED, lw=0.8)
        q0 = d['qmax']
        axes[1].plot(n, 10 * q0 ** n, color=MUTED, lw=0.8, label=r'$10\,q_0^{\,n}$')
        axes[1].legend(loc='upper right', fontsize=6.3)
    save(fig, 'fig_recon_convergence')


CERT_ROWS = [('functional_residual_max', 'functional equations, gluing pairs'),
             ('acs_residual_max', r'$\|J^2+I\|$'),
             ('N_max', r'$\|N_J\|$, 4000 domain points (max)'),
             ('N_adversarial_sup', r'$\|N_J\|$, domain, adversarial'),
             ('N_disc_Z1_max', r'$\|N_J\|$, disc around $Z_1$'),
             ('N_disc_Z2_max', r'$\|N_J\|$, disc around $Z_2$'),
             ('deck_g1_pairs_max', r'deck consistency $g_1$, gluing pairs'),
             ('deck_g2_pairs_max', r'deck consistency $g_2$, gluing pairs'),
             ('deck_g1_disc_max', r'deck consistency $g_1$, disc'),
             ('deck_g2_disc_max', r'deck consistency $g_2$, disc'),
             ('blend_g1_N_adversarial', r'$\|N\|$ blended across $g_1$-edge, adv.'),
             ('blend_g2_N_adversarial', r'$\|N\|$ blended across $g_2$-edge, adv.')]


def fig_recon_cert():
    cert = jload('reconstruction', 'certificate.json')
    if not cert:
        return
    fig, ax = plt.subplots(figsize=(4.6, 2.9))
    rows = [r for r in CERT_ROWS if r[0] in cert['S6']]
    ys = np.arange(len(rows))[::-1]
    for key, col, off, lab in (('S6', BLUE, 0.12, r'$S^6$'), ('S2xS4', ORANGE, -0.12, r'$S^2\times S^4$')):
        mark(ax, [max(cert[key][r[0]], 1e-17) for r in rows], ys + off, col, ms=4.5, label=lab)
    ax.set_yticks(ys); ax.set_yticklabels([r[1] for r in rows], fontsize=6.5)
    ax.set_xscale('log'); ax.set_xlim(1e-17, 1e-3); decade_ticks(ax, -16, -4, 'x', every=2)
    ax.axvline(1e-6, color=MUTED, lw=0.8); ax.text(1.2e-6, ys[0] + 0.5, 'target', fontsize=6, color=INK2, va='bottom')
    ax.axvline(2.2e-16, color=GRID, lw=0.8)
    ax.text(2.6e-16, ys[-1] - 0.75, 'double-precision $\\epsilon$', fontsize=5.8, color=INK2)
    ax.set_ylim(-1.0, len(rows) - 0.2)
    ax.grid(axis='y', visible=False)
    ax.set_xlabel('certificate quantity (fresh points; coordinate norms)', fontsize=7)
    ax.legend(loc='lower right', fontsize=6.5)
    save(fig, 'fig_recon_cert')


# ---------------------------------------------------------------------------
# tables
# ---------------------------------------------------------------------------

def tab_acs():
    lines = []
    for c in CASES:
        for s in (7, 8):
            r = res(f'acs_{c}_s{s}')
            if not r:
                continue
            L = ladder_entries(r)
            acc = [e for e in L if e['adv_acs_sup'] < 1.0]
            e = acc[0] if acc else (min(L, key=lambda e: e['adv_acs_sup']) if L else None)
            att = (e.get('attempt') or 'plain') if e else '--'
            att = {'plain': 'plain', 'dirichlet': 'Dirichlet', 'orient': 'orientation'}.get(att, att)
            cls = e['feat'].replace('_', r'\_') if e else '--'
            corr = r.get('corrected') or {}
            ok = r.get('certified_exact_acs')
            lines.append(' & '.join([
                TEX[c] if s == 7 else '', str(s), f'\\texttt{{{cls}}}' if acc else f'(\\texttt{{{cls}}})',
                att if acc else '--', str(len(L)), sci(e['adv_acs_sup']) if e else '--',
                sci(r['dense']['acs_mean']), sci(r['dense']['acs_sup']), sci(r.get('adv_acs_sup')),
                sci(corr.get('adv_acs_sup')) if corr else '--',
                r'\cmark' if ok else r'\xmark', f"{r['elapsed_min']:.0f}"]) + r' \\')
        lines.append(r'\addlinespace[2pt]')
    body = '\n'.join(lines[:-1])
    write_table('tab_acs', body + '\n')


def tab_cs_poly():
    lines = []
    for c in CASES:
        for s in (7, 8):
            name = f'cspoly_{c}_s{s}'
            r = res(name)
            if not r:
                continue
            if r.get('status') == 'no_certified_base':
                lines.append(f"{TEX[c] if s == 7 else ''} & {s} & \\multicolumn{{10}}{{l}}{{refused: own ACS not "
                             f"certified (adversarial sup {sci(r['base']['adv_acs_sup_raw'], 2)})}} \\\\")
                continue
            L = [e for e in r.get('ladder', []) if e.get('phase') == 'poly_cs']
            acc = [e for e in L if e.get('accepted')]
            ret = [e for e in L if str(e.get('degree')) == str(r.get('degree'))]   # the returned field
            last = acc[-1] if acc else (ret[-1] if ret else (L[-1] if L else {}))
            rounds = sum(len(e.get('lm', [])) for e in L)
            mean, mx = corrected(r)
            ok = r.get('certified_exact_acs')
            score = (last.get('final_N2', float('nan')) + 10 * last.get('final_C2', float('nan'))) if last else None
            an = jload('cs_analysis', 'done.json') or {}
            orth = an.get(name, {}).get('orth_defect_max')
            lines.append(' & '.join([
                TEX[c] if s == 7 else '', str(s), str(last.get('degree', '--')) + ('' if acc else '*'),
                f"{r['n_params']:,}".replace(',', '{,}'), str(rounds), sci(score),
                sci(r.get('adv_acs_sup')), sci(mean) if ok else 'n/a', sci(mx) if ok else 'n/a',
                sci(r.get('adv_nij_sup')) if ok else 'n/a',
                sci(orth, 2) if (orth is not None and ok) else '--', f"{r['elapsed_min']:.1f}"]) + r' \\')
        lines.append(r'\addlinespace[2pt]')
    write_table('tab_cs_poly', '\n'.join(lines[:-1]) + '\n')


def tab_conj():
    lines = []
    for name, lab in [(f'csconjexp_{c}', 'exp') for c in CLASSICAL + HARD] + \
                     [(f'csconjaff_{c}', 'affine') for c in CLASSICAL + HARD] + \
                     [(f'csconjtail_{c}', 'exp, tail') for c in HARD] + \
                     [(f'csconjtail_{c}_cont', 'exp, tail, continued') for c in HARD]:
        r = res(name)
        if not r:
            continue
        c = r['case']
        lines.append(' & '.join([TEX[c], lab, f"{r['steps_done']:,}".replace(',', '{,}'),
                                 sci(r.get('init', {}).get('nij')), sci(r['dense']['nij_mean']),
                                 sci(r['dense']['nij_sup']), sci(r.get('adv_nij_sup')),
                                 f"{r['elapsed_min']:.0f}"]) + r' \\')
    write_table('tab_conj', '\n'.join(lines) + '\n')


NAMB = {'s2': 3, 's4': 5, 's6': 7, 's1xs1': 4, 's2xs2': 6, 's3xs1': 6, 's3xs3': 8, 's4xs2': 8, 's5xs1': 8,
        's4xs1xs1': 9}


def mlp_params(r):
    """Parameter count of the conjugation MLP of a CS unit (GELU MLP, `stages` N x N outputs)."""
    n, h, st = NAMB[r['case']], r['hidden'], r['stages']
    dims = [n] + list(h) + [st * n * n]
    return sum(a * b + b for a, b in zip(dims[:-1], dims[1:]))


def tab_hard():
    lines = []
    for c in HARD:
        base = res(f'acs_{c}_s7')
        if base and base.get('corrected'):
            lines.append(' & '.join([TEX[c], 'own certified ACS $J^*$ (no CS search)', '--',
                                     sci(base['corrected']['nij_mean']), sci(base['corrected']['nij_sup']),
                                     sci(base.get('adv_nij_sup'))]) + r' \\')
        r = res(f'cspoly_{c}_s7')
        if r:
            L = [e for e in r.get('ladder', []) if e.get('phase') == 'poly_cs']
            sc = min((e.get('final_N2', np.inf) + 10 * e.get('final_C2', np.inf) for e in L), default=None)
            lines.append(' & '.join(['', f"free polynomial field, degree {r['degree']} (LM)",
                                     f"{r['n_params']:,}".replace(',', '{,}'),
                                     r'\multicolumn{3}{l}{not admissible: adv.\ sup $\|K^2+P\|=' +
                                     sci(r.get('adv_acs_sup'), 2).strip('$') + '$; loss ' + sci(sc, 2) + '}']) + r' \\')
        for name, lab in ((f'csconjexp_{c}', 'MLP conjugation, exp, mean-square'),
                          (f'csconjaff_{c}', 'MLP conjugation, affine, mean-square'),
                          (f'csconjtail_{c}', 'MLP conjugation, exp, tail-aware'),
                          (f'csconjtail_{c}_cont', 'the same, continued (+4000 steps)')):
            r = res(name)
            if r:
                lines.append(' & '.join(['', lab, f"{mlp_params(r):,}".replace(',', '{,}'),
                                         sci(r['dense']['nij_mean']), sci(r['dense']['nij_sup']),
                                         sci(r.get('adv_nij_sup'))]) + r' \\')
        for name, lab in ((f'polar_{c}', 'polar mode, degree 2'),):
            r = res(name)
            if r:
                L = [e for e in r.get('ladder', []) if e.get('phase') == 'polar_ladder']
                for e in L:
                    lines.append(' & '.join(['', f"polar mode, degree {e['degree']} (LM, barrier)",
                                             f"{e['n_params']:,}".replace(',', '{,}'), sci(e.get('N_mean')),
                                             sci(e.get('N_sup_dense')), sci(e.get('adv_N_sup'))]) + r' \\')
        for name, lab in ((f'conjpoly_{c}', ''), (f'lp_{c}', r', $L^8$')):
            for e in conjpoly_table(name, final=True):
                lines.append(' & '.join(['', f"polynomial generator, degree {e['degree']}{lab} (LM)",
                                         f"{e['n_params']:,}".replace(',', '{,}'), sci(e['N_mean']),
                                         sci(e['N_sup_dense']), sci(e['adv_N_sup'])]) + r' \\')
        lines.append(r'\addlinespace[3pt]')
    write_table('tab_hard', '\n'.join(lines[:-1]) + '\n')


def tab_recon():
    cert = jload('reconstruction', 'certificate.json')
    if not cert:
        return
    lines = []
    for key, lab in CERT_ROWS + [('D_max_on_domain', r'$\max D$ on the domain (lattice condition $D<0$)'),
                                 ('tau_Z1_minus_rho', r'$|\tau(Z_1)-e^{i\pi/3}|$ (not imposed)'),
                                 ('tau_Z2_minus_i', r'$|\tau(Z_2)-i|$ (not imposed)'),
                                 ('mu_Z1_minus_1_over_1_plus_rho', r'$|\mu(Z_1)-(1+e^{i\pi/3})^{-1}|$ (not imposed)'),
                                 ('mu_Z2_minus_(1-i)/2', r'$|\mu(Z_2)-(1-i)/2|$ (not imposed)'),
                                 ('certified_sup', r'\textbf{certified sup (max of the rows marked $\ast$)}')]:
        if key not in cert['S6']:
            continue
        star = r'$^\ast$' if key.startswith(('N_', 'blend_', 'deck_')) else ''
        f = (lambda v: f'${v:.2f}$') if key == 'D_max_on_domain' else sci
        lines.append(f"{lab}{star} & {f(cert['S6'][key])} & {f(cert['S2xS4'][key])} \\\\")
    cu = cert['cusp']
    lines.append(r'\addlinespace[2pt]')
    lines.append(f"cusp: $C(q)$ single-valued in $q$ & \\multicolumn{{2}}{{c}}{{{sci(cu['C_single_valued_err'])}}} \\\\")
    lines.append(f"cusp: $\\|R\\|_\\infty$, $\\epsilon_0=\\min(r/2, e^{{-2\\|R\\|_\\infty}})$ & "
                 f"\\multicolumn{{2}}{{c}}{{${cu['R_inf']:.2f}$, {sci(cu['eps0'])}}} \\\\")
    write_table('tab_recon', '\n'.join(lines) + '\n')


def tab_classical():
    d = jload('offline_checks', 'done.json')
    if not d:
        return
    names = {'s2xs2': r'product structure on $S^2\times S^2$', 's3xs1': r'Hopf structure on $S^3\times S^1$',
             's5xs1': r'Hopf structure on $S^5\times S^1$', 's3xs3': r'Calabi--Eckmann structure on $S^3\times S^3$'}
    lines = []
    for c in ('s2xs2', 's3xs1', 's5xs1', 's3xs3'):
        e = d['classical'].get(c)
        if e:
            lines.append(f"{names[c]} & {sci(e['acs_sup'])} & {sci(e['N_sup'])} & {sci(e['deg2_fit_residual'])} \\\\")
    write_table('tab_classical', '\n'.join(lines) + '\n')


BEFORE_CS = {'s2': 5.27e-4, 's4': 8.50, 's6': 9.87, 's1xs1': 1.90e-4, 's2xs2': 4.04e-2, 's3xs1': 2.00e-1,
             's3xs3': 7.63, 's4xs2': 6.28, 's5xs1': 6.42}


def tab_chart():
    """Earlier chart-based campaign (mean coordinate |N| over stereographic charts) vs the best
    unbiased run of this report, evaluated with the same chart protocol."""
    lines = []
    for c in CASES:
        b = best_unbiased(c)
        now = '--' if c == 's4' else (sci(b[1].get('chart_cs')) if b else '--')
        how = 'refused' if c == 's4' else (b[0].split('_')[0] if b else '--')
        lines.append(f"{TEX[c]} & {sci(BEFORE_CS.get(c))} & {now} & \\texttt{{{how}}} \\\\")
    write_table('tab_chart', '\n'.join(lines) + '\n')


def tab_summary():
    known = {'s2': ('yes', 'yes'), 's1xs1': ('yes', 'yes'), 's4': ('no', 'no'), 's2xs2': ('yes', 'yes'),
             's3xs1': ('yes', 'yes'), 's6': ('yes', r'yes \cite{AC26,Engel26}'), 's3xs3': ('yes', 'yes'),
             's5xs1': ('yes', 'yes'), 's4xs2': ('yes', r'yes \cite{LLZ26,GFL26}'), 's4xs1xs1': ('yes', r'yes \cite{Gao:2026tyy}')}
    cert = jload('reconstruction', 'certificate.json')
    lines = []
    for c in CASES:
        acs = [res(f'acs_{c}_s{s}') for s in (7, 8)]
        ok = [bool(r and r.get('certified_exact_acs')) for r in acs]
        acsv = ', '.join(sci(r.get('adv_acs_sup')) for r in acs if r)
        if c == 's4':
            acs_txt = r'refused (\xmark\xmark)'
        else:
            acs_txt = ('\\cmark' * sum(ok) + '\\xmark' * (2 - sum(ok))) if any(acs) else '--'
        b = best_unbiased(c)
        if c == 's4':
            cs_txt, cs_val, how = 'refused', '--', 'no certified ACS'
        elif b is None:
            cs_txt, cs_val, how = '--', '--', '--'
        else:
            mode = b[1].get('cs_mode', 'conj')
            lab = {'poly': 'polynomial field + LM', 'conj': 'MLP conjugation', 'conjpoly': 'polynomial generator',
                   'polar': 'polar mode', 'atlas': 'atlas'}.get(mode, mode)
            cs_txt, cs_val = lab, sci(b[1]['adv_nij_sup'])
            how = 'unbiased PINN' if b[1]['adv_nij_sup'] < 1e-6 else 'not found'
            if DIM[c] == 2:
                how = r'trivial ($N\equiv0$ in dim.\ 2)'
        rec = ''
        if cert and c in ('s6', 's4xs2'):
            rec = sci(cert['S6' if c == 's6' else 'S2xS4']['certified_sup'])
            how = 'reconstruction'
        c4 = jload('s4t2_reconstruction', 'certificate.json')
        if c4 and c == 's4xs1xs1':
            rec = sci(c4['paper']['certified_sup'])
            how = 'reconstruction'
        lines.append(' & '.join([TEX[c], str(DIM[c]), known[c][0], acs_txt, acsv, known[c][1], cs_txt, cs_val,
                                 rec or '--', how]) + r' \\')
    write_table('tab_summary', '\n'.join(lines) + '\n')


def tab_compute():
    """Wall-clock per job: both scheduler logs (the jobs kept from the interrupted
    first attempt and the rerun), the continuation job and the S^4 x T^2 runs."""
    rows = {}
    # scheduler_manual.log: jobs run by hand, in the scheduler's format (2026-10-25: extra S^3 x S^1 bases)
    logs = [os.path.join(OUT, n) for n in ('scheduler_run1.log', 'scheduler.log', 'scheduler_cont.log',
                                           'scheduler_manual.log')]
    for log in logs:
        if not os.path.exists(log):
            continue
        for line in open(log):
            m = re.search(r'DONE (\S+) \(([\d.]+) min', line)
            if m and os.path.exists(os.path.join(RUNS, m.group(1))):
                rows[m.group(1)] = float(m.group(2))
    for name, fn in (('s4t2_reconstruction', 'certificate.json'), ('s4t2inv_s7', 'done.json'),
                     ('s4t2inv_s8', 'done.json'), ('s4t2invconj_s7', 'done.json')):
        d = jload(name, fn)
        if d and d.get('elapsed_min'):
            rows[name] = float(d['elapsed_min'])
    if not rows:
        return
    groups = {}
    for name, mins in rows.items():
        g = name.split('_')[0]
        if name.endswith('_cont') and g == 'csconjtail':
            g = 'cont'
        if name == 'hopf_class':
            g = 'hopfclass'
        groups.setdefault(g, []).append(mins)
    lab = {'acs': 'ACS units (10 cases $\\times$ 2 seeds)', 'cspoly': 'CS, polynomial mode + LM',
           'csconjexp': 'CS, MLP conjugation, exp', 'csconjaff': 'CS, MLP conjugation, affine',
           'csconjtail': 'CS, MLP conjugation, tail-aware', 'cont': 'CS, tail-aware continuation ($S^4\\times T^2$)',
           'polar': 'CS, polar mode', 'conjpoly': 'CS, polynomial generator',
           'lp': r'CS, $L^p$ objective', 'atlas': 'CS, integrable atlas', 'fibre': 'fibration, topology only',
           'reconstruction': 'reconstruction $S^6$, $S^2\\times S^4$ (incl.\\ $K$ study)',
           'offline': 'offline checks', 'cs': 'CS analysis (post hoc)',
           's4t2': '$S^4\\times T^2$ reconstruction (20 starts)', 's4t2inv': '$S^4\\times T^2$, $T^2$-invariant polynomial',
           's4t2invconj': '$S^4\\times T^2$, $T^2$-invariant MLP conjugation (incl.\\ continuation)',
           's4t2invcp': '$S^4\\times T^2$, $T^2$-invariant polynomial generator',
           'hopf': 'supervised Hopf benchmark (Table 4)', 'hopfclass': '$S^3\\times S^1$ homotopy classes (post hoc)'}
    order = ['acs', 'cspoly', 'csconjexp', 'csconjaff', 'csconjtail', 'cont', 'polar', 'conjpoly', 'lp', 'atlas',
             'fibre', 'reconstruction', 'offline', 'cs', 's4t2', 's4t2inv', 's4t2invconj', 's4t2invcp', 'hopf',
             'hopfclass']
    def fmt(x):
        # the reconstruction jobs take seconds: "<1" rather than a rounded 0
        return '$<1$' if x < 1 else f'{x:.0f}'

    lines = []
    for g in sorted(groups, key=lambda g: order.index(g) if g in order else 99):
        v = groups[g]
        lines.append(f"{lab.get(g, g)} & {len(v)} & {fmt(sum(v))} & {fmt(max(v))} \\\\")
    lines.append(r'\midrule')
    tot = sum(rows.values())
    lines.append(f"total & {len(rows)} & {tot:.0f} & -- \\\\")
    write_table('tab_compute', '\n'.join(lines) + '\n')
    num('totalJobs', str(len(rows)))
    num('totalCPUHours', f'{tot / 60:.0f}')


def numbers():
    """Macros for the numbers quoted in the running text."""
    for c in CASES:
        key = re.sub(r'\d', lambda m: 'abcdefghij'[int(m.group())], c).replace('x', 'X')
        for s in (7, 8):
            r = res(f'acs_{c}_s{s}')
            if r:
                num(f'acsAdv{key}s{"g" if s == 7 else "h"}', sci(r.get('adv_acs_sup')))
                L = ladder_entries(r)
                acc = [e for e in L if e['adv_acs_sup'] < 1.0]
                if acc:
                    num(f'acsClass{key}s{"g" if s == 7 else "h"}', acc[0]['feat'].replace('_', r'\_'))
                    num(f'acsAttempt{key}s{"g" if s == 7 else "h"}', acc[0].get('attempt', 'plain'))
                if r.get('corrected'):
                    num(f'acsNijMean{key}s{"g" if s == 7 else "h"}', sci(r['corrected']['nij_mean'], 2))
                    num(f'acsNijAdv{key}s{"g" if s == 7 else "h"}', sci(r.get('adv_nij_sup'), 2))
                if c == 's4':
                    num(f'acsMinMean{key}s{"g" if s == 7 else "h"}', sci(min(e['acs_mean'] for e in L)) if L else '--')
                    num(f'acsMinAdv{key}s{"g" if s == 7 else "h"}', sci(min(e['adv_acs_sup'] for e in L), 2) if L else '--')
            r = res(f'cspoly_{c}_s{s}')
            if r and r.get('certified_exact_acs'):
                num(f'csAdv{key}s{"g" if s == 7 else "h"}', sci(r.get('adv_nij_sup')))
                num(f'csMin{key}s{"g" if s == 7 else "h"}', f"{r['elapsed_min']:.1f}")
        b = best_unbiased(c)
        if b:
            num(f'bestAdv{key}', sci(b[1]['adv_nij_sup']))
            m, x = corrected(b[1])
            num(f'bestMean{key}', sci(m))
    cert = jload('reconstruction', 'certificate.json')
    if cert:
        for k, lab in (('S6', 'Six'), ('S2xS4', 'TwoFour')):
            num(f'recSup{lab}', sci(cert[k]['certified_sup']))
            num(f'recNadv{lab}', sci(cert[k]['N_adversarial_sup']))
            num(f'recFun{lab}', sci(cert[k]['functional_residual_max']))
            num(f'recDmax{lab}', f"${cert[k]['D_max_on_domain']:.2f}$")
        num('recEpsZero', sci(cert['cusp']['eps0']))
        num('recRinf', f"${cert['cusp']['R_inf']:.2f}$")
    an = jload('cs_analysis', 'done.json') or {}
    for name, lab in (('csconjexp_s6', 'SixMS'), ('csconjtail_s6', 'SixTail'), ('conjpoly_s6', 'SixCP'),
                      ('csconjexp_s4xs2', 'FourTwoMS'), ('csconjaff_s4xs2', 'FourTwoAff'),
                      ('csconjexp_s5xs1', 'FiveOneMS'), ('csconjtail_s4xs1xs1_cont', 'FourOneOneTail'),
                      ('csconjtail_s4xs2', 'FourTwoTail')):
        if name in an and 'share_top' in an[name]:
            num(f'share{lab}', f"{100 * an[name]['share_top']['0.05']:.0f}")
            num(f'shareOne{lab}', f"{100 * an[name]['share_top']['0.01']:.0f}")
    po = [v for k, v in an.items() if k.startswith('cspoly_') and 'orth_defect_max' in v and v['case'] not in ('s2', 's1xs1')]
    if po:
        num('orthMin', f"{min(v['orth_defect_mean'] for v in po):.2f}")
        num('orthMax', f"{max(v['orth_defect_mean'] for v in po):.1f}")
        num('svMin', f"{min(v['sv_min'] for v in po):.2f}")
        num('svMax', f"{max(v['sv_max'] for v in po):.2f}")
    os.makedirs(TAB, exist_ok=True)
    with open(os.path.join(TAB, 'numbers.tex'), 'w') as f:
        f.write('% generated by experiments/make_figures.py -- do not edit\n')
        for k in sorted(NUMBERS):
            f.write(macro(k, NUMBERS[k]))
    print(f'wrote {len(NUMBERS)} number macros', flush=True)


# ---------------------------------------------------------------------------
# S^4 x S^1 x S^1 (Section 9): reconstruction and T^2-invariant search
# ---------------------------------------------------------------------------

def fig_s4t2_recon():
    path = os.path.join(RUNS, 's4t2_reconstruction', 'grids.npz')
    if not os.path.exists(path):
        return
    g = np.load(path)
    fig, axes = plt.subplots(1, 3, figsize=(6.3, 2.3), gridspec_kw={'wspace': 0.62})
    ax = axes[0]
    pc = ax.pcolormesh(g['U'], g['TH'], g['det_paper'], cmap=SEQ, shading='auto', rasterized=True)
    cb = fig.colorbar(pc, ax=ax, fraction=0.05, pad=0.03)
    cb.outline.set_visible(False); cb.ax.tick_params(labelsize=6, color=MUTED, labelcolor=INK2)
    ax.set_xlabel(r'$\log|t|$', fontsize=7); ax.set_ylabel(r'$\arg t$', fontsize=7)
    ax.set_yticks([-np.pi, 0, np.pi]); ax.set_yticklabels([r'$-\pi$', '0', r'$\pi$'])
    ax.set_title(r'normalised $\det_{\mathbb{R}}$, supplied periods', fontsize=7.5)
    ax.grid(False)
    ax = axes[1]
    for key, col, lab in (('det_paper', BLUE, 'supplied data'), ('det_learned', ORANGE, 'learned data')):
        if key in g:
            ax.plot(g['U'][:, 0], g[key].min(axis=1), color=col, lw=1.2, label=lab)
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_xlabel(r'$\log|t|$', fontsize=7); ax.set_ylabel(r'min over $\arg t$', fontsize=6.5)
    ax.set_title('non-degeneracy over $\\mathbb{C}^*$', fontsize=7.5)
    ax.legend(loc='lower center', fontsize=6)
    ax = axes[2]
    r = g['ratio_r']
    ax.axhspan(0.5, 1.5, color='#eef4fc', zorder=0)
    for key, col, lab in (('ratio_paper_0', BLUE, 'supplied, cusp 0'), ('ratio_paper_inf', BLUE, r'supplied, cusp $\infty$'),
                          ('ratio_learned_0', ORANGE, 'learned, cusp 0'), ('ratio_learned_inf', ORANGE, r'learned, cusp $\infty$')):
        if key in g:
            ls = '-' if key.endswith('_0') else (0, (3, 1.5))
            ax.plot(r, g[key][:, 1], color=col, lw=1.1, ls=ls, label=lab)
    ax.set_xscale('log'); ax.set_ylim(0.9, 2.0)
    ax.set_xlabel(r'$|t|$ (resp. $|s|$)', fontsize=7); ax.set_ylabel(r'$\max\,P/\log|t|$', fontsize=6.5)
    ax.set_title('filling condition', fontsize=7.5)
    ax.legend(loc='upper left', fontsize=5.6)
    save(fig, 'fig_s4t2_recon')


def tab_s4t2_recon():
    c = jload('s4t2_reconstruction', 'certificate.json')
    if not c:
        return
    rows = [('chirality_0', r'$\omega(\bar e_2,\bar e_3)$ at the cusp $t=0$'),
            ('chirality_inf', r'$\omega(\bar e_2,\bar e_3)$ at the cusp $t=\infty$'),
            ('det_normalised_min', r'$\min_{\mathbb{C}^*}$ normalised $\det_\mathbb{R}$ (adversarial)'),
            ('filling_eps_0', r'filling radius $\varepsilon$ at $0$ (Step 0 and (H2))'),
            ('filling_eps_inf', r'filling radius $\varepsilon$ at $\infty$'),
            ('acs_residual_max', r'$\|J^2+I\|$'),
            ('N_max', r'$\|N_J\|$ on the regular part$^\ast$'),
            ('N_disc_0_max', r'$\|N_J\|$ on the filling disc at $0$ ($u=\log t$)$^\ast$'),
            ('N_disc_inf_max', r'$\|N_J\|$ on the filling disc at $\infty$$^\ast$'),
            ('monodromy_consistency_max', r'monodromy $e_4\mapsto e_4+e_1$, consistency$^\ast$'),
            ('monodromy_wrong_transport', r'\quad the same with the wrong transport'),
            ('chart_infinity_consistency_max', r'chart change $\Theta$ at $\infty$, consistency$^\ast$'),
            ('log_transform_consistency_max', r'logarithmic transform $\Phi$, consistency$^\ast$'),
            ('certified_sup', r'\textbf{certified sup (starred rows)}')]
    cols = [('paper', 'supplied'), ('learned_best', 'learned')]
    lines = []
    for k, lab in rows:
        vals = []
        for key, _ in cols:
            v = c.get(key, {}).get(k)
            if v is None:
                vals.append('--')
            elif k.startswith(('chirality', 'det_')):
                vals.append(f'${v:+.3f}$')
            elif k.startswith('filling'):
                vals.append(f'${v:.3g}$')
            else:
                vals.append(sci(v))
        lines.append(f"{lab} & {' & '.join(vals)} \\\\")
    write_table('tab_s4t2_recon', '\n'.join(lines) + '\n')
    ls = c.get('learned_summary', {})
    num('sftAdmissible', str(ls.get('admissible', '--')))
    num('sftStarts', str(ls.get('starts', '--')))
    num('sftSupPaper', sci(c['paper']['certified_sup']))
    num('sftDetPaper', f"${c['paper']['det_normalised_min']:.2f}$")
    if 'learned_best' in c:
        num('sftSupLearned', sci(c['learned_best']['certified_sup']))
        num('sftDetLearned', f"${c['learned_best']['det_normalised_min']:.2f}$")
    learned = jload('s4t2_reconstruction', 'learned.json') or []
    ok = [r for r in learned if r.get('feasible') and r.get('certificate', {}).get('admissible')]
    if ok:
        num('sftNullDefect', f"{min(r['null_conic_defect'] for r in ok):.2f}--{max(r['null_conic_defect'] for r in ok):.2f}")
        pair = sum(1 for r in ok if r['relation_c2_plus_a2'] < 0.1 and r['relation_c3_minus_a3'] < 0.1)
        num('sftPairCount', {0: 'none', 1: 'one', 2: 'two', 3: 'three'}.get(pair, str(pair)))
        num('sftSupLearnedMax', sci(max(r['certificate']['certified_sup'] for r in ok)))


def fig_s4t2_invariant():
    runs = [(n, jload(n, 'done.json')) for n in ('s4t2inv_s7', 's4t2inv_s8')]
    runs = [(n, d) for n, d in runs if d and d.get('cs')]
    if not runs:
        return
    fig, axes = plt.subplots(1, 2, figsize=(6.3, 2.3), gridspec_kw={'wspace': 0.28})
    for (name, d), col in zip(runs, (BLUE, ORANGE)):
        off = 0
        for e in d['cs']:
            x = np.arange(off + 1, off + len(e['lm']) + 1)
            lab = f"seed {d['seed']}" if off == 0 else None
            axes[0].plot(x, [q['N2'] for q in e['lm']], color=col, lw=1.2, marker='o', ms=2.5, mec='white', mew=0.4, label=lab)
            axes[1].plot(x, [q['C2'] for q in e['lm']], color=col, lw=1.2, marker='o', ms=2.5, mec='white', mew=0.4)
            if len(runs) == 1 or col == BLUE:
                for ax in axes:
                    ax.axvline(off + 0.5, color=GRID, lw=0.8, zorder=0)
                axes[0].text(off + 0.7, 0.98, f"deg {e['degree']}", transform=mpl.transforms.blended_transform_factory(
                    axes[0].transData, axes[0].transAxes), fontsize=5.8, color=INK2, va='top')
            off += len(e['lm'])
    for ax in axes:
        ax.set_yscale('log'); ax.set_xlabel('Levenberg–Marquardt round', fontsize=6.5)
        ax.grid(axis='x', visible=False)
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
    axes[0].set_ylabel(r'held-out $\mathbb{E}\|N_K\|^2$', fontsize=6.5)
    axes[1].set_ylabel(r'held-out $\mathbb{E}\|K^2+P\|^2$', fontsize=6.5)
    axes[0].set_title(r'$T^2$-invariant $S^4\times T^2$: Nijenhuis term', fontsize=7.5)
    axes[1].set_title('algebraic term', fontsize=7.5)
    if len(runs) > 1:
        axes[0].legend(loc='upper right', fontsize=6.3)
    save(fig, 'fig_s4t2_invariant')


def tab_s4t2_invariant():
    lines = []
    for name in ('s4t2inv_s7', 's4t2inv_s8'):
        d = jload(name, 'done.json')
        if not d:
            continue
        a = d.get('acs', {})
        acc = [e for e in d.get('acs_ladder', []) if e['adv_acs_sup'] < 1]
        lab = f"{acc[0]['degree']}, {acc[0]['attempt']}" if acc else '--'
        lines.append(f"{d['seed']} & own ACS (degree {lab}) & {acc[0]['n_params'] if acc else '--'} & -- & "
                     f"-- & {sci(a.get('N_mean'))} & {sci(a.get('N_max'))} & {sci(a.get('N_adv'))} \\\\")
        for e in d.get('cs', []):
            last = e['lm'][-1] if e['lm'] else {}
            loss = (last.get('N2', np.nan) + 10 * last.get('C2', np.nan)) if last else None
            adm = e.get('adv_acs_sup')
            lines.append(f" & degree {e['degree']} & {e['n_params']:,} & {len(e['lm'])} & {sci(loss)} & "
                         .replace(',', '{,}') +
                         (f"{sci(e.get('N_mean'))} & {sci(e.get('N_max'))} & {sci(e.get('N_adv'))} \\\\"
                          if 'N_mean' in e else f"\\multicolumn{{3}}{{l}}{{not admissible (adv.\\ $\\|K^2+P\\|$ {sci(adm).strip('$')})}} \\\\"))
        lines.append(r'\addlinespace[2pt]')
    dc = jload('s4t2invconj_s7', 'done.json')
    dcc = jload('s4t2invconj_s7_cont', 'done.json')
    npar = 5 * 128 + 128 + 2 * (128 * 128 + 128) + 128 * 98 + 98
    for dd, lab in ((dc, None), (dcc, 'cont')):
        if not dd:
            continue
        what = (f"MLP conjugation, tail-aware Adam ({dd['steps']:,} steps)" if lab is None
                else f"the same, continued (+{dd['steps']:,} steps)")
        lines.append(f"{dd['seed']} & {what} & {npar:,} & -- & -- & ".replace(',', '{,}') +
                     f"{sci(dd['N_mean'])} & {sci(dd['N_max'])} & {sci(dd['N_adv'])} \\\\")
    if dc or dcc:
        lines.append(r'\addlinespace[2pt]')
    dcp_ = jload('s4t2invcp_s7', 'done.json')
    if dcp_:
        for e in dcp_.get('conjpoly', []):
            lines.append(f"{dcp_['seed']} & polynomial generator, degree {e['degree']} (LM) & {e['n_params']:,} & {len(e['lm'])} & "
                         .replace(',', '{,}') + f"{sci(e['lm'][-1]['N2']) if e['lm'] else '--'} & {sci(e['N_mean'])} & "
                         f"{sci(e['N_max'])} & {sci(e['N_adv'])} \\\\")
        lines.append(r'\addlinespace[2pt]')
    if lines:
        write_table('tab_s4t2_invariant', '\n'.join(lines[:-1]) + '\n')
        best = None
        for name in ('s4t2inv_s7', 's4t2inv_s8'):
            d = jload(name, 'done.json') or {}
            for e in d.get('cs', []):
                if 'N_adv' in e and (best is None or e['N_adv'] < best[0]):
                    best = (e['N_adv'], e['N_mean'], e['degree'])
        for dd in (dc, dcc):
            if dd and (best is None or dd['N_adv'] < best[0]):
                best = (dd['N_adv'], dd['N_mean'], 'MLP')
        words = {1: 'One', 2: 'Two', 3: 'Three', 4: 'Four'}
        di = jload('s4t2inv_s7', 'done.json') or {}
        if di.get('acs'):
            num('sftInvAcsMean', sci(di['acs']['N_mean'])); num('sftInvAcsAdv', sci(di['acs']['N_adv']))
        for e in di.get('cs', []):
            w = words.get(e['degree'])
            if w and 'N_adv' in e:
                num(f'sftInvPolyMean{w}', sci(e['N_mean'])); num(f'sftInvPolyAdv{w}', sci(e['N_adv']))
        for dd, w in ((dc, 'Conj'), (dcc, 'Cont')):
            if dd:
                num(f'sftInv{w}Mean', sci(dd['N_mean'])); num(f'sftInv{w}Adv', sci(dd['N_adv']))
                num(f'sftInv{w}Max', sci(dd['N_max']))
        dcp = jload('s4t2invcp_s7', 'done.json') or {}
        for e in dcp.get('conjpoly', []):
            w = words.get(e['degree'])
            num(f'sftInvCPMean{w}', sci(e['N_mean'])); num(f'sftInvCPAdv{w}', sci(e['N_adv']))
        if dcp.get('conjpoly'):
            bcp = min(dcp['conjpoly'], key=lambda e: e['N_adv'])
            num('sftInvCPBestAdv', sci(bcp['N_adv'])); num('sftInvCPBestDeg', str(bcp['degree']))
            num('sftInvCPBestMean', sci(bcp['N_mean']))
        if best:
            num('sftInvBestAdv', sci(best[0])); num('sftInvBestMean', sci(best[1])); num('sftInvBestDeg', f'degree-{best[2]} class' if isinstance(best[2], int) else 'MLP conjugation')


def fig_s4t2_search():
    """All unbiased attempts on S^4 x T^2: training curves of the conjugation runs and the
    final adversarial suprema of every route."""
    fig = plt.figure(figsize=(6.3, 4.3))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], hspace=0.55, wspace=0.3)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])]
    cv = conj_curve('csconjtail_s4xs1xs1')
    cont = conj_curve('csconjtail_s4xs1xs1_cont')
    if cv:
        off = cv[-1][3]
        full = cv + [(v[0], v[1], v[2], v[3] + off) for v in cont]
        axes[0].plot([v[3] for v in full], [v[1] for v in full], color=BLUE, lw=1.2, label='full model')
        axes[1].plot([v[3] for v in full], [v[2] for v in full], color=BLUE, lw=1.2)
    d = jload('s4t2invconj_s7', 'done.json')
    h = d['hist'] if d else []
    if not h:                                   # running: parse the log
        pth = os.path.join(RUNS, 's4t2invconj_s7', 'log.txt')
        if os.path.exists(pth):
            for line in open(pth):
                m = re.search(r'step\s+(\d+) \| batch \S+ \| held-out mean (\S+) max (\S+)', line)
                if m:
                    h.append({'step': int(m.group(1)), 'N_mean': float(m.group(2)), 'N_max': float(m.group(3))})
    dcont = jload('s4t2invconj_s7_cont', 'done.json')
    hc = dcont['hist'] if dcont else []
    if not hc:
        pth = os.path.join(RUNS, 's4t2invconj_s7_cont', 'log.txt')
        if os.path.exists(pth):
            for line in open(pth):
                m = re.search(r'step\s+(\d+) \| batch \S+ \| held-out mean (\S+) max (\S+)', line)
                if m:
                    hc.append({'step': int(m.group(1)), 'N_mean': float(m.group(2)), 'N_max': float(m.group(3))})
    if h and hc:
        off = h[-1]['step']
        h = h + [dict(q, step=q['step'] + off) for q in hc if q['step'] > 1]
        for ax in axes[:2]:
            ax.axvline(off, color=GRID, lw=0.8, zorder=0)
    if h:
        axes[0].plot([q['step'] for q in h], [q['N_mean'] for q in h], color=ORANGE, lw=1.2, label=r'$T^2$-invariant')
        axes[1].plot([q['step'] for q in h], [q['N_max'] for q in h], color=ORANGE, lw=1.2)
    for ax, lab in ((axes[0], r'held-out mean $\|N_J\|$'), (axes[1], r'held-out max $\|N_J\|$')):
        ax.set_yscale('log'); ax.set_xlabel('Adam step (tail-aware)', fontsize=6.5); ax.set_ylabel(lab, fontsize=6.5)
    axes[0].set_title('MLP conjugation', fontsize=7.5); axes[1].set_title('MLP conjugation', fontsize=7.5)
    axes[0].legend(loc='upper right', fontsize=6.2)
    # final adversarial suprema of every unbiased route
    rows = []
    b = res('acs_s4xs1xs1_s7')
    if b:
        rows.append(('own ACS (full model)', b.get('adv_nij_sup'), BLUE))
    for name, lab in (('csconjexp_s4xs1xs1', 'MLP conj., mean-square'), ('csconjtail_s4xs1xs1', 'MLP conj., tail, 1500 st.'),
                      ('csconjtail_s4xs1xs1_cont', 'MLP conj., tail, 5500 st.'), ('cspoly_s4xs1xs1_s7', 'free poly 2t (LM)')):
        r = res(name)
        if r and r.get('adv_nij_sup') is not None and (r.get('cs_mode') != 'poly' or r.get('certified_exact_acs')):
            rows.append((lab, r['adv_nij_sup'], BLUE))
    di = jload('s4t2inv_s7', 'done.json')
    if di:
        if di.get('acs'):
            rows.append((r'own ACS ($T^2$-inv.)', di['acs']['N_adv'], ORANGE))
        for e in di.get('cs', []):
            if 'N_adv' in e:
                rows.append((f"free poly deg {e['degree']} (LM)", e['N_adv'], ORANGE))
    if d:
        rows.append((f"MLP conj., tail, {d['steps']} st.", d['N_adv'], ORANGE))
    if dcont:
        rows.append((f"MLP conj., tail, {d['steps'] + dcont['steps'] if d else dcont['steps']} st.", dcont['N_adv'], ORANGE))
    dcp = jload('s4t2invcp_s7', 'done.json')
    if dcp:
        for e in dcp.get('conjpoly', []):
            rows.append((f"poly generator deg {e['degree']} (LM)", e['N_adv'], ORANGE))
    ax = axes[2]
    ys = np.arange(len(rows))[::-1]
    for y, (lab, v, col) in zip(ys, rows):
        mark(ax, v, y, col, ms=4.2)
    ax.set_yticks(ys); ax.set_yticklabels([r[0] for r in rows], fontsize=6.3)
    ax.set_xscale('log'); ax.grid(axis='y', visible=False)
    ax.set_xlabel(r'adversarial sup of $\|N_J\|$ (exact ACS)', fontsize=6.8)
    ax.set_title('every unbiased route on $S^4\\times T^2$', fontsize=7.5)
    hh = [Line2D([], [], marker='o', ls='none', color=BLUE, mec='white', label='full model $S^4\\times S^1\\times S^1$'),
          Line2D([], [], marker='o', ls='none', color=ORANGE, mec='white', label='$T^2$-invariant reduced model')]
    ax.legend(handles=hh, loc='upper center', fontsize=6.3)
    save(fig, 'fig_s4t2_search')


def hopf_class_of(r):
    """Hopf invariant, relative to the target, of the base of a supervised S^3 x S^1 run
    (runs/hopf_class/done.json); 0 for S^5 x S^1 (one homotopy class per orientation); None if unknown."""
    if r.get('case') != 's3xs1':
        return 0
    if 'hopf_invariant' in r.get('base', {}):            # recorded by the supervised script since 2026-10-25
        return r['base']['hopf_invariant']
    hc = jload('hopf_class', 'done.json')
    unit = os.path.basename(r.get('base', {}).get('dir', ''))
    return (hc or {}).get('units', {}).get(unit, {}).get('hopf_invariant')


def tab_supervised():
    """Paper Table 4: mean intrinsic Nijenhuis norm of supervised fits of the Hopf structure
    (z -> e z, extra/ambient_supervised_hopf.py; 4000 fresh points) against the unsupervised
    conjugation units of Figure 2 (6000 fresh points). Ranges are over the listed runs.
    On S^3 x S^1 a base outside the homotopy class of the target (Hopf invariant != 0,
    experiments/hopf_class.py) cannot be conjugated to it: such runs get their own row."""
    lines = []
    for c in ('s3xs1', 's5xs1'):
        rows = []
        val = res(f'hopf_{c}_val_b7')
        sob = [res(n) for n in sorted(os.listdir(RUNS)) if n.startswith(f'hopf_{c}_sob_')] \
            if os.path.isdir(RUNS) else []
        sob = [r for r in sob if r]
        if val:
            rows.append(('supervised, values only', sci(val['ambient']['nij_mean'], 2),
                         sci(val.get('adv_nij_sup'), 2), '1'))
        groups = {True: [], False: []}
        for r in sob:
            groups[hopf_class_of(r) in (0, None)].append(r)
        for ok, lab in ((True, 'supervised, values and derivatives'),
                        (False, 'supervised, values and derivatives, base in another class')):
            if groups[ok]:
                m = [r['ambient']['nij_mean'] for r in groups[ok]]
                a = [r.get('adv_nij_sup') for r in groups[ok]]
                rng = sci(min(m), 2) + ('--' + sci(max(m), 2) if len(m) > 1 else '')
                rows.append((lab, rng, sci(max(a), 2), str(len(m))))
        for name, lab in ((f'csconjexp_{c}', 'unsupervised, exponential conjugation'),
                          (f'csconjaff_{c}', 'unsupervised, affine conjugation')):
            r = res(name)
            if r:
                rows.append((lab, sci(r['dense']['nij_mean'], 2), sci(r.get('adv_nij_sup'), 2), '1'))
        for k, (lab, mean, adv, nrun) in enumerate(rows):
            lines.append(' & '.join([TEX[c] if k == 0 else '', lab, mean, adv, nrun]) + r' \\')
        lines.append(r'\addlinespace[2pt]')
    if lines:
        write_table('tab_supervised', '\n'.join(lines[:-1]) + '\n')


ALL = [fig_overview, fig_acs_ladder, fig_acs_refine, fig_s4, fig_orientation, fig_cs_poly, fig_cs_conj_vs_poly,
       fig_hard_conj, fig_conjpoly, fig_conjpoly_rounds, fig_polar, fig_atlas, fig_fibre, fig_defect,
       fig_recon_domain, fig_recon_convergence, fig_recon_cert,
       tab_acs, tab_cs_poly, tab_conj, tab_hard, tab_recon, tab_classical, tab_summary, tab_compute, tab_chart,
       fig_s4t2_recon, tab_s4t2_recon, fig_s4t2_invariant, tab_s4t2_invariant, fig_s4t2_search, tab_supervised,
       numbers]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', default=None)
    a = ap.parse_args()
    style()
    todo = ALL if not a.only else [f for f in ALL if f.__name__ in a.only.split(',')]
    for f in todo:
        try:
            f()
        except Exception:
            print(f'FAILED {f.__name__}', flush=True)
            traceback.print_exc()
    if not a.only:
        with open(os.path.join(OUT, 'figures_done.txt'), 'w') as fh:
            fh.write('ok\n')


if __name__ == '__main__':
    main()
