"""E3 figure: combined (CP) and methylation-sensitive (MS) models at realistic theta. -> pdf/e3/e3_combined.pdf"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import e3  # noqa: E402

PROJ = e3.PROJ
OUT = f'{PROJ}/pdf/e3'
BLUE, ORANGE, AQUA, INK, MUTED, GRID = '#2a78d6', '#eb6834', '#1baf7a', '#0b0b0b', '#8a8984', '#e4e3df'
plt.rcParams.update({'font.size': 8, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6})
ORDER = ['M0', 'Pcre_k10', 'MS_k10', 'CP_k10', 'CP_k40']
SHORT = {'M0': 'neutral', 'Pcre_k10': 'P (CRE)', 'MS_k10': 'meth-sensitive\nTF', 'CP_k10': 'TF + meth\nκ=10',
         'CP_k40': 'TF + meth\nκ=40', 'P_k10': 'P (E2b)', 'C_k10_Rmatch': 'C (E2b)'}
COL = {'M0': MUTED, 'Pcre_k10': ORANGE, 'MS_k10': AQUA, 'CP_k10': BLUE, 'CP_k40': BLUE, 'P_k10': ORANGE,
       'C_k10_Rmatch': BLUE}


def main():
    os.makedirs(OUT, exist_ok=True)
    runs = pd.read_csv(f'{e3.CSV}/runs.csv')
    tab = pd.read_csv(f'{e3.CSV}/summary.csv').set_index('cond')
    ref = pd.read_csv(f'{PROJ}/csv/e2b/summary.csv').set_index('cond')
    fig, ax = plt.subplots(1, 4, figsize=(15, 3.8), constrained_layout=True)

    a = ax[0]
    R = np.linspace(0.5, 5, 100)
    es = e3.CONDS['CP_k10'][3]
    a.plot(R, np.minimum(es * np.exp(R / e3.e2.R0), 1), color=INK, lw=1, label='optimum O·e^(−R/R₀) = E*')
    for c, mk in (('CP_k10', 'o'), ('CP_k40', 's')):
        g = runs[runs.cond == c]
        a.scatter(g.R_mean, g.O_mean + (0.015 if c == 'CP_k40' else -0.015), s=16, marker=mk, alpha=0.7,
                  facecolor=BLUE if c == 'CP_k10' else 'white', edgecolor=BLUE, lw=0.8, label=SHORT[c].replace('\n', ' '))
    for c in ('P_k10', 'C_k10_Rmatch'):
        a.scatter(ref.loc[c, 'R_mean'], ref.loc[c, 'O_mean'], s=50, marker='*', color=COL[c], label=SHORT[c] + ' mean')
    a.set(xlabel='methylation load R (run mean)', ylabel='TF occupancy O (run mean)', ylim=(-0.03, 1.03),
          title='A  TF + methylation: where on the trade-off?')
    a.legend(frameon=False, fontsize=6.5)

    a = ax[1]
    conds = ORDER + ['P_k10', 'C_k10_Rmatch']
    for i, c in enumerate(conds):
        t = tab.loc[c] if c in tab.index else ref.loc[c]
        a.errorbar(i, t.PD_methCpG_loss_rel_M0, yerr=2 * t.PD_se, fmt='o', ms=5, color=COL[c], mec='white')
    a.axhline(1, color=MUTED, lw=0.8, ls=':')
    a.axvline(len(ORDER) - 0.5, color=GRID, lw=1)
    a.set_xticks(range(len(conds)), [SHORT[c] for c in conds], fontsize=6.5, rotation=45, ha='right')
    a.set(yscale='log', ylabel='P/D, loss of methylated CpGs outside motif (÷ neutral)', title='B  P/D')

    a = ax[2]
    classes = [('R_only', 'changes R only'), ('motif', 'motif (non-CpG)'), ('motif_CpG', 'motif CpG')]
    for k, (cls, nm) in enumerate(classes):
        for i, c in enumerate(ORDER[1:]):
            y = tab.loc[c, f'daf_{cls}_rel_M0']
            if np.isfinite(y):
                a.plot(k + (i - 1.5) * 0.16, y, 'o', ms=5, color=COL[c], mec=COL[c], mew=1.2,
                       mfc=COL[c] if c != 'CP_k40' else 'white', label=SHORT[c].replace('\n', ' ') if k == 0 else None)
    a.axhline(1, color=MUTED, lw=0.8, ls=':')
    a.set_xticks(range(3), [nm for _, nm in classes])
    a.set(yscale='log', ylabel='mean DAF ÷ neutral (same class)', title='C  allele frequency by effect class')
    a.legend(frameon=False, fontsize=6.5)

    a = ax[3]
    for k, (cls, nm) in enumerate(classes):
        for i, c in enumerate(ORDER[1:]):
            y = tab.loc[c, f'2Ns_het_{cls}']
            if np.isfinite(y) and y > 0:
                a.plot(k + (i - 1.5) * 0.16, y, 'o', ms=5, color=COL[c], mec=COL[c], mew=1.2,
                       mfc=COL[c] if c != 'CP_k40' else 'white')
    a.axhline(1, color=MUTED, lw=0.8, ls=':')
    a.set_xticks(range(3), [nm for _, nm in classes])
    a.set(yscale='log', ylabel='realised 2Ns (heterozygote, at optimum)', title='D  selection per variant class')
    fig.savefig(f'{OUT}/e3_combined.pdf')
    fig.savefig(f'{OUT}/e3_combined.png', dpi=150)


if __name__ == '__main__':
    main()
