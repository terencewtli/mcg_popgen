"""E1 figure: passenger-methylation footprint. SLiM = colored marks (mean over replicates, ±2 SE); oracle = line.
-> pdf/e1/e1_passenger.pdf (+ .png)"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJ = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
CSV, PDF = f'{PROJ}/csv/e1', f'{PROJ}/pdf/e1'
SERIES = ['#2a78d6', '#eb6834', '#1baf7a']
INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e4e3df'
plt.rcParams.update({'font.size': 8, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'lines.linewidth': 1.4})


def profile_panel(ax, prof, conds, col, ylabel, title, labels):
    for c, cond in zip(SERIES, conds):
        s = prof[(prof.cond == cond) & (prof.source == 'slim')].groupby('bin')[col]
        o = prof[(prof.cond == cond) & (prof.source == 'oracle')].set_index('bin')[col]
        x = s.mean().index + 5
        ax.plot(o.index + 5, o.values, color=c, lw=1.2, alpha=0.9)
        ax.errorbar(x, s.mean(), yerr=2 * s.sem(), fmt='o', ms=4, color=c, mec='white', mew=0.6, elinewidth=1,
                    label=labels[cond])
    ax.set(xlabel='distance from motif edge (bp)', ylabel=ylabel, title=title)


def main():
    os.makedirs(PDF, exist_ok=True)
    prof = pd.read_csv(f'{CSV}/profiles.csv')
    sm = pd.read_csv(f'{CSV}/summary.csv').set_index('cond')
    lab = {'M0_ctx': 'neutral, context-only', 'P_ctx': 'TF selected, context-only', 'P_meth': 'TF selected, meth-gated',
           'M0_meth': 'neutral, meth-gated', 'M0_meth_fb': 'neutral, gated + density fb',
           'P_meth_fb': 'TF selected, gated + density fb'}
    fig, ax = plt.subplots(2, 3, figsize=(12.5, 7), constrained_layout=True)

    profile_panel(ax[0, 0], prof, ['M0_ctx', 'P_ctx', 'P_meth'], 'density', 'CpGs per 100 bp',
                  'A  CpG retention around a selected motif', lab)
    ax[0, 0].legend(frameon=False, title='dots: SLiM ±2SE; lines: oracle', title_fontsize=7)
    profile_panel(ax[0, 1], prof, ['M0_meth', 'P_meth', 'P_meth_fb'], 'meth', 'mean methylation of CpGs',
                  'B  methylation (never enters fitness)', lab)
    ax[0, 1].legend(frameon=False)
    profile_panel(ax[0, 2], prof, ['M0_meth', 'M0_meth_fb', 'P_meth_fb'], 'density', 'CpGs per 100 bp',
                  'C  density feedback, with and without a TF', lab)
    ax[0, 2].legend(frameon=False)

    def contrast_panel(a, conds, xt, title):
        x = np.arange(len(conds))
        y, e = sm.loc[conds, 'contrast'], sm.loc[conds, 'contrast_se']
        a.errorbar(x, y, yerr=2 * e, fmt='o', ms=6, color=SERIES[0], mec='white', label='SLiM')
        a.plot(x, sm.loc[conds, 'oracle_contrast'], '_', ms=18, mew=2, color=INK, label='oracle (motif frozen)')
        a.axhline(1, color=MUTED, lw=0.8, ls=':')
        a.set_xticks(x, xt)
        a.set(ylabel='CpG density, <25 bp / ≥75 bp', title=title)

    contrast_panel(ax[1, 0], ['P_meth_lam10', 'P_meth', 'P_meth_lam50'], ['λ=10', 'λ=25', 'λ=50'],
                   'D  island contrast vs protection range λ')
    ax[1, 0].legend(frameon=False)
    contrast_panel(ax[1, 1], ['P_meth_r5', 'P_meth', 'P_meth_r50'], ['r=5', 'r=20', 'r=50'],
                   'E  island contrast vs CpG hypermutability r')
    q = ['P_meth_Q1', 'P_meth', 'P_meth_Q4']
    contrast_panel(ax[1, 2], q, [f'Q={k}\nN={int(sm.loc[c, "N"])}' for k, c in zip([1, 2, 4], q)],
                   'F  rescaling check (2Ns, θ, t·u fixed)')
    fig.savefig(f'{PDF}/e1_passenger.pdf')
    fig.savefig(f'{PDF}/e1_passenger.png', dpi=150)


if __name__ == '__main__':
    main()
