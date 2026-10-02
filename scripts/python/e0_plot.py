"""E0 validation figure: observed (SLiM, colored marks) vs expected (oracle/analytic, dark line). -> pdf/e0/e0_validation.pdf"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJ = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
CSV, PDF = f'{PROJ}/csv/e0', f'{PROJ}/pdf/e0'
SERIES = ['#2a78d6', '#eb6834', '#1baf7a']          # categorical slots 1-3 (validated all-pairs)
EXPECT, MUTED, GRID = '#0b0b0b', '#52514e', '#e4e3df'

plt.rcParams.update({'font.size': 8, 'axes.edgecolor': MUTED, 'axes.labelcolor': EXPECT, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'lines.linewidth': 2})


def main():
    os.makedirs(PDF, exist_ok=True)
    fig, ax = plt.subplots(1, 4, figsize=(13, 3.2), constrained_layout=True)

    # A: CpG trajectories, mean over replicates, SLiM vs oracle
    tr = pd.read_csv(f'{CSV}/traj.csv')
    labels = {'ctx_r20': 'context r=20', 'ctx_r10': 'context r=10', 'gate_r20_g0.5': 'r=20, gate 0.5'}
    for c, (cond, lab) in zip(SERIES, labels.items()):
        g = tr[tr.cond == cond].groupby('gen')
        m, se = g.mean_cpg.mean(), g.mean_cpg.std() / np.sqrt(g.size())
        ax[0].plot(m.index, g.oracle_mean.mean(), color=EXPECT, lw=1.2, zorder=1)
        ax[0].errorbar(m.index, m, yerr=2 * se, fmt='o', ms=4, color=c, mec='white', mew=0.8, elinewidth=1,
                       label=lab, zorder=2)
    ax[0].set(xlabel='generation', ylabel='CpGs per 1 kb haplotype', title='A  CpG decay: SLiM vs oracle')
    ax[0].legend(frameon=False, title='dots: SLiM ±2SE; line: oracle', title_fontsize=7)

    # B: pi across rescaling Q
    pi = pd.read_csv(f'{CSV}/pi.csv')
    theta = 0.01
    exp_pi = theta / (1 + 4 * theta / 3)
    for k, (q, g) in enumerate(pi.groupby('Q')):
        x = np.full(len(g), k) + np.linspace(-0.15, 0.15, len(g))
        ax[1].scatter(x, g.pi * 1e3, s=10, color=SERIES[0], alpha=0.5, lw=0)
        se = g.pi.std() / np.sqrt(len(g))
        ax[1].errorbar(k, g.pi.mean() * 1e3, yerr=2 * se * 1e3, fmt='s', ms=6, color=SERIES[0], mec='white', zorder=3)
    ax[1].axhline(exp_pi * 1e3, color=EXPECT, lw=1.2)
    ax[1].text(2.4, exp_pi * 1e3, ' JC\n expectation', va='center', fontsize=7, color=MUTED)
    ax[1].set_xticks(range(3), [f'Q={q}\nN={n}' for q, n in pi.groupby('Q').N.first().items()])
    ax[1].set(ylabel='π per kb', title='B  neutral π, θ=0.01, rescaled', xlim=(-0.5, 2.9))

    # C: folded SFS (Q=1), observed vs expected
    sfs = pd.read_csv(f'{CSV}/sfs.csv')
    for c, (q, g) in zip(SERIES, sfs.groupby('Q')):
        ax[2].plot(g.i, g.observed, 'o', ms=4, color=c, mec='white', mew=0.6, label=f'Q={q}')
    g = sfs[sfs.Q == 1]
    ax[2].plot(g.i, g.expected, color=EXPECT, lw=1.2, label='1/i + 1/(n−i)')
    ax[2].set(yscale='log', xlabel='minor allele count (n=50)', ylabel='proportion of SNPs', title='C  folded SFS')
    ax[2].legend(frameon=False)

    # D: fixation probability relative to neutral, vs 2Ns
    fx = pd.read_csv(f'{CSV}/fix.csv')
    xs = np.linspace(-2.5, 10.5, 200)
    n0 = 100
    ku = [(-np.expm1(-x / (2 * n0)) / -np.expm1(-x) if x != 0 else 1 / (2 * n0)) * 2 * n0 for x in xs]
    ax[3].plot(xs, ku, color=EXPECT, lw=1.2, label='Kimura (N=100)')
    for c, (n, g) in zip(SERIES, fx.groupby('N', sort=False)):
        se = np.sqrt(g.p_hat * (1 - g.p_hat) / g.K)
        ax[3].errorbar(g.two_Ns, g.p_hat * 2 * n, yerr=2 * se * 2 * n, fmt='o', ms=5, color=c, mec='white',
                       label=f'SLiM N={n}')
    ax[3].set(yscale='log', xlabel='2Ns', ylabel='P(fix) × 2N', title='D  fixation, 20,000 trials each')
    ax[3].legend(frameon=False)

    fig.savefig(f'{PDF}/e0_validation.pdf')
    fig.savefig(f'{PDF}/e0_validation.png', dpi=150)


if __name__ == '__main__':
    main()
