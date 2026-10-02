"""E2 figure: causal (C) vs passenger (P) methylation under stabilising selection. -> pdf/e2/e2_causal_passenger.pdf
Colour = model (C blue, P orange, M0 grey); marker = selection strength (open kappa = 2, filled kappa = 10)."""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))

PROJ = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
CSV, PDF = f'{PROJ}/csv/e2', f'{PROJ}/pdf/e2'
BLUE, ORANGE, INK, MUTED, GRID = '#2a78d6', '#eb6834', '#0b0b0b', '#8a8984', '#e4e3df'
STYLE = {'M0': (MUTED, 'o', True), 'P_k2': (ORANGE, 'o', False), 'P_k10': (ORANGE, 'o', True),
         'C_k2': (BLUE, 's', False), 'C_k10': (BLUE, 's', True), 'C_k10_Rmatch': (BLUE, 'D', True)}
LABEL = {'M0': 'neutral', 'P_k2': 'passenger κ=2', 'P_k10': 'passenger κ=10', 'C_k2': 'causal κ=2',
         'C_k10': 'causal κ=10', 'C_k10_Rmatch': 'causal κ=10, R-matched'}
plt.rcParams.update({'font.size': 8, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'lines.linewidth': 1.4})


def mk(cond):
    c, m, filled = STYLE[cond]
    return dict(color=c, marker=m, mfc=c if filled else 'white', mec=c, ms=5, mew=1.2)


def main():
    os.makedirs(PDF, exist_ok=True)
    tv = pd.read_csv(f'{CSV}/turnover.csv')
    tab = pd.read_csv(f'{CSV}/summary.csv').set_index('cond')
    conds = list(STYLE)
    fig, ax = plt.subplots(1, 4, figsize=(14.5, 3.7), constrained_layout=True)

    for a, (xc, yc, xl, yl, t) in zip(ax[:2], [
            ('cpg_turnover', 'dR', 'CpG positional turnover (1 − Jaccard)', '|Δ methylation load R|',
             'A  methylation load vs CpG turnover'),
            ('motif_turnover', 'dO', 'motif positional turnover', '|Δ occupancy O|',
             'B  occupancy vs motif turnover')]):
        for c in conds:
            g = tv[tv.cond == c].groupby('lag')[[xc, yc]].mean()
            a.plot(g[xc], g[yc], ls='-', label=LABEL[c], **mk(c))
        a.set(xlabel=xl, ylabel=yl, title=t)
    ax[0].legend(frameon=False, fontsize=7, title='lags 500–32,000 gen', title_fontsize=7)

    a = ax[2]
    classes = [('R_only', 'changes R only'), ('O', 'changes O'), ('null', 'neither')]
    for k, (cls, nm) in enumerate(classes):
        for i, c in enumerate(conds):
            x = k + (i - 2.5) * 0.12
            y, e = tab.loc[c, f'maf_{cls}'], tab.loc[c, f'maf_{cls}_se']
            if np.isfinite(y):
                a.errorbar(x, y, yerr=2 * e if np.isfinite(e) else None, ls='none', elinewidth=1,
                           label=LABEL[c] if k == 0 else None, **mk(c))
    a.set_xticks(range(3), [nm for _, nm in classes])
    a.set(ylabel='mean minor-allele frequency', title='C  allele frequency by effect class')

    a = ax[3]
    for i, c in enumerate(conds):
        y, e = tab.loc[c, 'PD_methCpG_loss_rel_M0'], tab.loc[c, 'PD_se']
        a.errorbar(i, y, yerr=2 * e, ls='none', elinewidth=1, **mk(c))
        s = tab.loc[c, '2Ns_het_R_only']
        a.annotate(f'2Ns={s:.1f}', (i, y), xytext=(4, 6), textcoords='offset points', fontsize=6.5, color=MUTED)
    a.axhline(1, color=MUTED, lw=0.8, ls=':')
    short = {'M0': 'neutral', 'P_k2': 'P κ=2', 'P_k10': 'P κ=10', 'C_k2': 'C κ=2', 'C_k10': 'C κ=10',
             'C_k10_Rmatch': 'C κ=10\nR-matched'}
    a.set_xticks(range(len(conds)), [short[c] for c in conds])
    a.set(ylabel='P/D, loss of methylated CpGs (÷ neutral)', title='D  polymorphism / divergence')
    fig.savefig(f'{PDF}/e2_causal_passenger.pdf')
    fig.savefig(f'{PDF}/e2_causal_passenger.png', dpi=150)


if __name__ == '__main__':
    main()
