"""E4: meQTL detection layer on the realistic-theta variants of E2b and E3.

  python scripts/python/e4.py     # -> csv/e4/*.csv, pdf/e4/e4_detection.pdf

Each segregating variant snapshot (csv/e2b/variants.csv, csv/e3/variants.csv) has a true derived-allele frequency p
and a per-allele effect beta = dR on regional methylation load. Phenotype: the individual's methylation load (sum
over its two haplotypes) plus N(0, SIGMA^2) noise from measurement and biology. A variant is detected with
  power = P(chi2_1(lambda) > c),  lambda = n 2p(1-p) beta^2 / SIGMA^2,  c = chi2_1 quantile at P_THRESH,
for MAF >= MAF_MIN. Detected sets are power-weighted (expected values; no extra sampling noise).
SIGMA = 3: a methylated CpG-loss variant (|beta| ~ 0.95) at MAF 0.1 explains ~2% of variance.
Classes: CpG-SNP meQTL (changes R, not occupancy) and motif-SNP meQTL (changes occupancy and R through protection).
Statistics, per model x n x class:
  n_detected    expected detections per locus snapshot
  maf_detected  power-weighted mean MAF of detected variants; rel_null = ratio to the neutral model under the same
                detection (the ascertainment-matched test)
  corr_beta_maf power-weighted Pearson correlation of |beta| and MAF among detected variants (the naive
                "effect-frequency" relationship; nonzero under neutrality = detection bias)
  frac_increase fraction of detected CpG-SNP meQTLs whose derived allele increases methylation
SEs are bootstrap over replicates.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2, ncx2

PROJ = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
CSV, PDF = f'{PROJ}/csv/e4', f'{PROJ}/pdf/e4'
SIGMA, P_THRESH, MAF_MIN = 3.0, 1e-6, 0.01
NS = [500, 5000, 30000]
CRIT = chi2.isf(P_THRESH, 1)
SOURCES = {'e2b': ['M0', 'P_k10', 'C_k2', 'C_k10_Rmatch', 'C_k10'], 'e3': ['M0', 'Pcre_k10', 'MS_k10', 'CP_k10', 'CP_k40']}


def load():
    out = []
    for src, conds in SOURCES.items():
        v = pd.read_csv(f'{PROJ}/csv/{src}/variants.csv', usecols=['daf', 'dR', 'dO', 'cond', 'rep'])
        v = v[v.cond.isin(conds)].copy()
        v['src'] = src
        out.append(v)
    v = pd.concat(out, ignore_index=True)
    v['maf'] = np.minimum(v.daf, 1 - v.daf)
    v['abs_beta'] = v.dR.abs()
    v['cls'] = np.where((v.dO.abs() < 1e-9) & (v.abs_beta >= 0.05), 'CpG_SNP',
                        np.where((v.dO.abs() >= 0.05) & (v.abs_beta >= 0.05), 'motif_SNP', 'none'))
    return v[(v.cls != 'none') & (v.maf >= MAF_MIN)]


def power(maf, beta, n):
    lam = n * 2 * maf * (1 - maf) * beta ** 2 / SIGMA ** 2
    return ncx2.sf(CRIT, 1, lam)


def stats(g, w):
    if w.sum() <= 0:
        return dict(maf=np.nan, corr=np.nan, inc=np.nan)
    maf = np.average(g.maf, weights=w)
    x, y = g.abs_beta.values, g.maf.values
    mx, my = np.average(x, weights=w), np.average(y, weights=w)
    cov = np.average((x - mx) * (y - my), weights=w)
    sx, sy = np.sqrt(np.average((x - mx) ** 2, weights=w)), np.sqrt(np.average((y - my) ** 2, weights=w))
    return dict(maf=maf, corr=cov / (sx * sy) if sx > 0 and sy > 0 else np.nan,
                inc=np.average(g.dR.values > 0, weights=w))


def main():
    os.makedirs(CSV, exist_ok=True)
    os.makedirs(PDF, exist_ok=True)
    v = load()
    nsnap = {'e2b': 1000, 'e3': 1000}          # variant snapshots per run (MEAS / VEVERY = 200000 / 200)
    rows = []
    rng = np.random.default_rng(0)
    for n in NS:
        v[f'w{n}'] = power(v.maf.values, v.abs_beta.values, n)
    for (src, cond, cls), g in v.groupby(['src', 'cond', 'cls']):
        reps = g.rep.unique()
        null = v[(v.src == src) & (v.cond == 'M0') & (v.cls == cls)]
        for n in NS:
            w = g[f'w{n}'].values
            st = stats(g, w)
            sn = stats(null, null[f'w{n}'].values) if len(null) else dict(maf=np.nan)
            boots = []
            byrep = {r: d for r, d in g.groupby('rep')}
            nullrep = {r: d for r, d in null.groupby('rep')} if len(null) else {}
            for _ in range(300):
                gg = pd.concat([byrep[r] for r in rng.choice(reps, len(reps))])
                b = stats(gg, gg[f'w{n}'].values)
                if nullrep:
                    nr = list(nullrep)
                    nn = pd.concat([nullrep[r] for r in rng.choice(nr, len(nr))])
                    b['rel'] = b['maf'] / stats(nn, nn[f'w{n}'].values)['maf']
                boots.append(b)
            bt = pd.DataFrame(boots)
            rows.append(dict(src=src, cond=cond, cls=cls, n=n, reps=len(reps),
                             n_detected_per_snapshot=w.sum() / (len(reps) * nsnap[src]),
                             maf_detected=st['maf'], maf_detected_se=bt['maf'].std(),
                             maf_rel_null=st['maf'] / sn['maf'] if np.isfinite(sn['maf']) else np.nan,
                             maf_rel_null_se=bt['rel'].std() if 'rel' in bt else np.nan,
                             corr_beta_maf=st['corr'], corr_se=bt['corr'].std(),
                             frac_increase=st['inc'], frac_increase_se=bt['inc'].std(),
                             maf_all_segregating=g.maf.mean()))
    tab = pd.DataFrame(rows)
    tab.to_csv(f'{CSV}/detection.csv', index=False)
    with pd.option_context('display.width', 250, 'display.max_rows', 200):
        print(tab[['src', 'cond', 'cls', 'n', 'n_detected_per_snapshot', 'maf_detected', 'maf_rel_null',
                   'maf_rel_null_se', 'corr_beta_maf', 'frac_increase']].to_string(index=False,
                                                                                    float_format=lambda x: f'{x:.3g}'))
    plot(tab)


def plot(tab):
    INK, MUTED, GRID = '#0b0b0b', '#8a8984', '#e4e3df'
    COL = {500: '#9ec5f4', 5000: '#2a78d6', 30000: '#0d3a73'}
    plt.rcParams.update({'font.size': 8, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                         'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6})
    order = [('e2b', 'M0'), ('e2b', 'P_k10'), ('e2b', 'C_k2'), ('e2b', 'C_k10_Rmatch'), ('e2b', 'C_k10'),
             ('e3', 'Pcre_k10'), ('e3', 'MS_k10'), ('e3', 'CP_k10'), ('e3', 'CP_k40')]
    lab = ['neutral', 'P', 'C κ=2', 'C κ=10\nR-match', 'C κ=10', 'P (CRE)', 'meth-sens\nTF', 'TF+meth\nκ=10',
           'TF+meth\nκ=40']
    fig, ax = plt.subplots(1, 3, figsize=(14, 3.8), constrained_layout=True)
    for a, (col, se, yl, t, cls) in zip(ax, [
            ('maf_rel_null', 'maf_rel_null_se', 'mean MAF of detected meQTLs ÷ neutral (same detection)',
             'A  CpG-SNP meQTLs, ascertainment-matched', 'CpG_SNP'),
            ('maf_rel_null', 'maf_rel_null_se', 'mean MAF of detected meQTLs ÷ neutral (same detection)',
             'B  motif-SNP meQTLs (few variants without a TF)', 'motif_SNP'),
            ('corr_beta_maf', 'corr_se', 'corr(|β|, MAF) among detected', 'C  naive effect–frequency correlation',
             'CpG_SNP')]):
        for j, n in enumerate(NS):
            ys, es = [], []
            for (src, cond) in order:
                r = tab[(tab.src == src) & (tab.cond == cond) & (tab.cls == cls) & (tab.n == n)]
                ys.append(r[col].iloc[0] if len(r) else np.nan)
                es.append(r[se].iloc[0] if len(r) else np.nan)
            x = np.arange(len(order)) + (j - 1) * 0.22
            ys, es = np.asarray(ys, float), 2 * np.asarray(es, float)
            if cls == 'motif_SNP':                   # log scale: drop error bars that cross zero (few variants)
                es = np.where(ys - es > 0, es, np.nan)
            a.errorbar(x, ys, yerr=es, fmt='o', ms=4.5, color=COL[n], mec='white', elinewidth=1,
                       label=f'n = {n:,}')
        a.axhline(1 if col == 'maf_rel_null' else 0, color=MUTED, lw=0.8, ls=':')
        if cls == 'motif_SNP' and col == 'maf_rel_null':
            a.set_yscale('log')
        a.axvline(4.5, color=GRID, lw=1)
        a.set_xticks(range(len(order)), lab, fontsize=6.5)
        a.set(ylabel=yl, title=t)
    ax[0].legend(frameon=False)
    fig.savefig(f'{PDF}/e4_detection.pdf')
    fig.savefig(f'{PDF}/e4_detection.png', dpi=150)


if __name__ == '__main__':
    main()
