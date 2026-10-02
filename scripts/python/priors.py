"""P0: primate-scale priors for CpG evolution under different selection / mutation models.

  python scripts/python/priors.py run     # SLiM check of the PRF polymorphism density -> sim/p0/
  python scripts/python/priors.py calc    # analytic priors + PRF check -> csv/p0/, pdf/p0/priors.pdf

Per-CpG two-state model (theory.pdf §9) with origin-fixation weights, the weak-mutation limit validated in E0:
  loss rate  = alpha * l(m) * phi(-S),  l(m) = 2 (1 + (r - 1) m) + 4       (m = GERMLINE methylation)
  gain rate  = alpha * gamma * phi(+S), gamma = 6/15                       (matches E0: 0.90 vs 0.94 CpG/100 bp)
  phi(S) = S / (1 - e^-S): substitution rate relative to neutral for population-scaled advantage S = 2Ns
  stationary f = gamma e^S / (gamma e^S + l(m))  (Sella-Hirsh: neutral x e^{S * state})
Time is measured in units of neutral pairwise divergence D at non-CpG sites: alpha * T_pair = D / 3.
D values are approximate genome-wide figures: human-chimp ~0.012, human-macaque ~0.065. S > 0 means a CpG is
favoured. That covers selection on the CpG itself and also gBGC, which favours GC alleles with B = 4 Ne b, often
quoted ~0.3-1. The two can be told apart because gBGC also acts on non-CpG W<->S sites.
"""
import glob
import os
import subprocess
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.integrate import quad

PROJ = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SLIM = os.path.expanduser('~/miniconda3/envs/slim4/bin/slim')
SIM, CSV, PDF = f'{PROJ}/sim/p0', f'{PROJ}/csv/p0', f'{PROJ}/pdf/p0'
GAMMA = 6 / 15
PAIRS = {'human-chimp': 0.012, 'human-macaque': 0.065}
PRF_N, PRF_L, PRF_MU, PRF_NSAMP, PRF_SIGMAS, PRF_REPS = 500, 100000, 1e-7, 20, [0, -1, -2, -4, -8], 6


def phi(S):
    S = np.asarray(S, float)
    with np.errstate(divide='ignore', invalid='ignore'):
        out = S / -np.expm1(-S)
    return np.where(np.abs(S) < 1e-9, 1.0, out)


def loss_mult(m, r):
    return 2 * (1 + (r - 1) * m) + 4


def stationary_f(m, S, r):
    return GAMMA * np.exp(S) / (GAMMA * np.exp(S) + loss_mult(m, r))


def oe(m, S, r):
    return stationary_f(m, S, r) * 16            # uniform base composition reference


def p_present_given_present(m, S, r, D):
    """P(CpG at the orthologous position in species B | CpG in species A), stationary reversible two-state chain."""
    a = loss_mult(m, r) * phi(-S)
    b = GAMMA * phi(S)
    f = b / (a + b)
    t = D / 3.0                                    # alpha * T_pair
    return f + (1 - f) * np.exp(-(a + b) * t)


def prf_density(x, sigma):
    """Sawyer-Hartl density (per unit theta) of derived-allele frequency for new mutations, sigma = 2Ns."""
    if abs(sigma) < 1e-9:
        return 1.0 / x
    return -np.expm1(-sigma * (1 - x)) / (-np.expm1(-sigma) * x * (1 - x))


def prf_seg(sigma, n):
    return quad(lambda x: prf_density(x, sigma) * (1 - x ** n - (1 - x) ** n), 0, 1, limit=200)[0]


def prf_singletons(sigma, n):
    return quad(lambda x: prf_density(x, sigma) * n * x * (1 - x) ** (n - 1), 0, 1, limit=200)[0]


# ---------------------------------------------------------------------------------------------------------------
def prf_task(job):
    sigma, rep = job
    out = f'{SIM}/prf_sigma{sigma}_rep{rep}.tsv'
    if os.path.exists(out):
        return
    s = sigma / (2 * PRF_N)
    r = subprocess.run([SLIM, '-s', str(500 + 10 * rep - sigma), '-d', f'N={PRF_N}', '-d', f'L={PRF_L}',
                        '-d', f'MU={PRF_MU}', '-d', f'S={s}', '-d', f'NSAMP={PRF_NSAMP}', '-d', 'NSAMPLES=10',
                        '-d', f"OUT='{out}'", f'{PROJ}/scripts/slim/p0_prf.slim'], capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(r.stdout[-1500:] + r.stderr[-1500:])


def run():
    os.makedirs(SIM, exist_ok=True)
    with Pool(9) as pool:
        pool.map(prf_task, [(s, k) for s in PRF_SIGMAS for k in range(PRF_REPS)], chunksize=1)


def calc():
    os.makedirs(CSV, exist_ok=True)
    os.makedirs(PDF, exist_ok=True)
    ms = np.linspace(0, 1, 101)
    Ss = np.linspace(0, 4, 81)
    grid = pd.DataFrame([dict(r=r, m=m, S=S, f=stationary_f(m, S, r), oe=oe(m, S, r))
                         for r in (10, 15, 20) for m in ms for S in Ss])
    grid.to_csv(f'{CSV}/oe_grid.csv', index=False)
    cons = pd.DataFrame([dict(pair=p, D=D, r=15, m=m, S=S, p_cons=p_present_given_present(m, S, 15, D))
                         for p, D in PAIRS.items() for m in (0.0, 0.5, 1.0) for S in Ss])
    cons.to_csv(f'{CSV}/conservation.csv', index=False)
    n = PRF_NSAMP
    pd_rows = []
    for S in np.linspace(0, 8, 33):
        sig = -S                                     # CpG-loss allele is deleterious when a CpG is favoured
        P = prf_seg(sig, n) / prf_seg(0, n)
        Dv = float(phi(sig))
        pd_rows.append(dict(S=S, poly_rel=P, div_rel=Dv, pd_ratio=P / Dv,
                            singleton_frac=prf_singletons(sig, n) / prf_seg(sig, n)))
    pdr = pd.DataFrame(pd_rows)
    pdr.to_csv(f'{CSV}/pd_ratio.csv', index=False)
    # PRF check against SLiM
    theta_L = 4 * PRF_N * PRF_MU * PRF_L
    chk = []
    for sig in PRF_SIGMAS:
        files = sorted(glob.glob(f'{SIM}/prf_sigma{sig}_rep*.tsv'))
        v = np.concatenate([pd.read_csv(f, sep='\t').seg.values for f in files])
        reps = np.array([pd.read_csv(f, sep='\t').seg.mean() for f in files])   # replicate means (samples within a run are correlated)
        e = theta_L * prf_seg(sig, n)
        chk.append(dict(sigma=sig, n_samples=len(v), observed=reps.mean(), se=reps.std(ddof=1) / np.sqrt(len(reps)),
                        expected=e, z=(reps.mean() - e) / (reps.std(ddof=1) / np.sqrt(len(reps)))))
    chk = pd.DataFrame(chk)
    chk.to_csv(f'{CSV}/prf_check.csv', index=False)
    # key numbers
    key = []
    for r in (10, 15, 20):
        key.append(dict(quantity=f'o/e, m=1, S=0, r={r}', value=oe(1, 0, r)))
        key.append(dict(quantity=f'o/e, m=0, S=0, r={r}', value=oe(0, 0, r)))
        # S needed on a fully methylated CpG to reach CGI-like o/e 0.6
        target = 0.6 / 16
        key.append(dict(quantity=f'S for o/e=0.6 at m=1, r={r}',
                        value=np.log(target * loss_mult(1, r) / (GAMMA * (1 - target)))))
        key.append(dict(quantity=f'm for o/e=0.6 at S=0, r={r}',
                        value=((GAMMA * (1 - target) / target - 4) / 2 - 1) / (r - 1)))
    for p, D in PAIRS.items():
        for m in (0.0, 1.0):
            for S in (0, 1, 2, 4):
                key.append(dict(quantity=f'P(CpG in B | CpG in A), {p}, m={m:g}, S={S}, r=15',
                                value=p_present_given_present(m, S, 15, D)))
    pd.DataFrame(key).to_csv(f'{CSV}/key_numbers.csv', index=False)
    with pd.option_context('display.width', 200, 'display.max_colwidth', 80):
        print(chk.to_string(index=False, float_format=lambda x: f'{x:.3g}'))
        print(pd.DataFrame(key).to_string(index=False, float_format=lambda x: f'{x:.3g}'))
        print(pdr.iloc[::4].to_string(index=False, float_format=lambda x: f'{x:.3g}'))
    plot(grid, cons, pdr, chk)


def plot(grid, cons, pdr, chk):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    SER = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']
    INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e4e3df'
    plt.rcParams.update({'font.size': 8, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                         'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'lines.linewidth': 2})
    fig, ax = plt.subplots(1, 4, figsize=(14, 3.6), constrained_layout=True)

    a = ax[0]
    g = grid[grid.S == 0]
    lo, hi = g[g.r == 20].set_index('m').oe, g[g.r == 10].set_index('m').oe
    a.fill_between(lo.index, lo.values, hi.values, color=SER[0], alpha=0.15, lw=0, label='S=0, r = 10–20')
    for c, S in zip(SER, (0, 1, 2, 3)):
        h = grid[(grid.r == 15) & np.isclose(grid.S, S)]
        a.plot(h.m, h.oe, color=c, label=f'S = 2Ns = {S}')
    for y, t in ((0.22, 'bulk genome ~0.2–0.25'), (0.6, 'CGI threshold 0.6')):
        a.axhline(y, color=MUTED, lw=0.8, ls=':')
        a.text(0.99, y, t, ha='right', va='bottom', fontsize=7, color=MUTED)
    a.set(xlabel='germline methylation m', ylabel='CpG o/e at equilibrium', yscale='log',
          title='A  CpG o/e: mutation (m) vs selection (S)')
    a.legend(frameon=False, fontsize=7)

    a = ax[1]
    h = grid[grid.r == 15].pivot(index='S', columns='m', values='oe')
    im = a.contourf(h.columns, h.index, np.log10(h.values), levels=20, cmap='Blues')
    cs = a.contour(h.columns, h.index, h.values, levels=[0.22, 0.4, 0.6, 0.8], colors=INK, linewidths=1)
    a.clabel(cs, fmt='o/e=%.2g', fontsize=7)
    a.set(xlabel='germline methylation m', ylabel='S = 2Ns favouring the CpG',
          title='B  same o/e, different causes (r=15)')

    a = ax[2]
    for c, m in zip(SER, (1.0, 0.0)):
        for ls, (p, D) in zip(('-', '--'), PAIRS.items()):
            h = cons[(cons.pair == p) & (cons.m == m)]
            a.plot(h.S, h.p_cons, color=c, ls=ls, label=f'm={m:g}, {p}')
    a.set(xlabel='S = 2Ns favouring the CpG', ylabel='P(CpG in B | CpG in A)', ylim=(0, 1.02),
          title='C  orthologous CpG conservation')
    a.legend(frameon=False, fontsize=7)

    a = ax[3]
    a.plot(pdr.S, pdr.poly_rel, color=SER[0], label='polymorphism (PRF)')
    a.plot(pdr.S, pdr.div_rel, color=SER[1], label='divergence (Kimura)')
    a.plot(pdr.S, pdr.pd_ratio, color=INK, lw=1.2, label='P/D ratio')
    neutral = chk[chk.sigma == 0].observed.iloc[0]
    a.errorbar(-chk.sigma, chk.observed / neutral, yerr=2 * chk.se / neutral, fmt='o', ms=5, color=SER[0],
               mec='white', label='SLiM polymorphism ±2SE')
    a.axhline(1, color=MUTED, lw=0.8, ls=':')
    a.set(yscale='log', xlabel='S = 2Ns against CpG loss', ylabel='relative to neutral (n = 20)',
          title='D  rate change vs selection: P/D')
    a.legend(frameon=False, fontsize=7)
    fig.savefig(f'{PDF}/priors.pdf')
    fig.savefig(f'{PDF}/priors.png', dpi=150)


if __name__ == '__main__':
    {'run': run, 'calc': calc}[sys.argv[1]]()
