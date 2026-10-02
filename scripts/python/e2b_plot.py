"""E2b figure: E2 statistics at theta = 0.024 (E2) vs realistic theta = 0.001 (E2b). -> pdf/e2b/e2b_theta.pdf
Colour = theta (grey: E2, blue: E2b); x = condition."""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJ = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
OUT = f'{PROJ}/pdf/e2b'
HI, LO, INK, MUTED, GRID = '#8a8984', '#2a78d6', '#0b0b0b', '#52514e', '#e4e3df'
CONDS = ['M0', 'P_k2', 'P_k10', 'C_k2', 'C_k10_Rmatch', 'C_k10']
SHORT = ['neutral', 'P κ=2', 'P κ=10', 'C κ=2', 'C κ=10\nR-matched', 'C κ=10']
plt.rcParams.update({'font.size': 8, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6})


def main():
    os.makedirs(OUT, exist_ok=True)
    hi = pd.read_csv(f'{PROJ}/csv/e2/summary.csv').set_index('cond').loc[CONDS]
    lo = pd.read_csv(f'{PROJ}/csv/e2b/summary.csv').set_index('cond').loc[CONDS]
    tv = pd.read_csv(f'{PROJ}/csv/e2/turnover.csv')
    tv16 = tv[tv.lag == 16].groupby('cond')[['cpg_turnover', 'dR']].mean().loc[CONDS]   # 8000 gen = 0.5/u
    x = np.arange(len(CONDS))
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6), constrained_layout=True)

    a = ax[0]
    a.errorbar(x - 0.12, hi.PD_methCpG_loss_rel_M0, yerr=2 * hi.PD_se, fmt='o', ms=5, color=HI, mec='white',
               label='θ = 0.024 (E2)')
    a.errorbar(x + 0.12, lo.PD_methCpG_loss_rel_M0, yerr=2 * lo.PD_se, fmt='o', ms=5, color=LO, mec='white',
               label='θ = 0.001 (E2b, ~human)')
    a.axhline(1, color=MUTED, lw=0.8, ls=':')
    a.set(yscale='log', ylabel='P/D, loss of methylated CpGs (÷ neutral)', title='A  P/D')
    a.legend(frameon=False)

    a = ax[1]
    m0h, m0l = hi.loc['M0', 'maf_R_only'], lo.loc['M0', 'daf_R_only']
    a.errorbar(x - 0.12, hi.maf_R_only / m0h, yerr=2 * hi.maf_R_only_se / m0h, fmt='o', ms=5, color=HI, mec='white',
               label='E2: mean MAF')
    a.errorbar(x + 0.12, lo.daf_R_only / m0l, yerr=2 * lo.daf_R_only_se / m0l, fmt='o', ms=5, color=LO, mec='white',
               label='E2b: mean DAF')
    a.axhline(1, color=MUTED, lw=0.8, ls=':')
    a.set(ylabel='frequency of R-only variants (÷ neutral)', title='B  allele frequency, variants changing R only',
          ylim=(0, 1.2))
    a.legend(frameon=False)

    a = ax[2]
    for c, sh, i in zip(CONDS, SHORT, x):
        a.plot([tv16.loc[c, 'cpg_turnover'], lo.loc[c, 'cpg_turnover_full']], [tv16.loc[c, 'dR'], lo.loc[c, 'dR_full']],
               color=MUTED, lw=0.8)
        a.plot(tv16.loc[c, 'cpg_turnover'], tv16.loc[c, 'dR'], 'o', ms=5, color=HI, mec='white')
        a.plot(lo.loc[c, 'cpg_turnover_full'], lo.loc[c, 'dR_full'], 'o', ms=5, color=LO, mec='white')
        a.annotate(sh.replace('\n', ' '), (lo.loc[c, 'cpg_turnover_full'], lo.loc[c, 'dR_full']), xytext=(4, -9),
                   textcoords='offset points', fontsize=6.5, color=MUTED)
    a.set(xlabel='CpG positional turnover over 0.5/u', ylabel='|Δ methylation load R| over 0.5/u',
          title='C  turnover needs high θ when steps are strong')
    for a in ax[:2]:
        a.set_xticks(x, SHORT)
    fig.savefig(f'{OUT}/e2b_theta.pdf')
    fig.savefig(f'{OUT}/e2b_theta.png', dpi=150)


if __name__ == '__main__':
    main()
