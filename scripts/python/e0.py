"""E0 validation: SLiM nucleotide/WF machinery vs analytic or oracle expectations. No methylation, no TF.

  python scripts/python/e0.py run       # SLiM + oracle runs -> sim/e0/ (raw) and csv/e0/ (summaries)
  python scripts/python/e0.py check     # recompute summaries from sim/e0/ and print the pass/fail table

Checks (each reported as z = (observed - expected) / SE; |z| < 3 is a pass):
  traj : CpG count trajectory under CpG hypermutation, SLiM population vs single-lineage oracle, per replicate,
         for context-only hypermutation (r=20, r=10) and a constant mutation() gate (r=20, gate=0.5 == r=10).
  pi   : neutral pairwise diversity vs JC finite-sites expectation pi = theta / (1 + 4 theta / 3), theta = 4 N u,
         across rescaling factors Q (N/Q, u*Q, r*Q, generations/Q).
  sfs  : folded SFS proportions vs the neutral 1/i expectation (pooled over replicates, per Q).
  fix  : fixation probability of a single new codominant mutation vs Kimura u = (1 - e^-s) / (1 - e^-2Ns).
"""
import os
import subprocess
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import oracle  # noqa: E402

PROJ = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SLIM = os.path.expanduser('~/miniconda3/envs/slim4/bin/slim')
SIM, CSV = f'{PROJ}/sim/e0', f'{PROJ}/csv/e0'
NPROC = 9

# trajectory test
TRAJ = dict(N=100, L=1000, ALPHA=1e-5, GENS=20000, NCHK=20, NSAMP=2, RREC=0.0)
TRAJ_CONDS = {'ctx_r20': (20, 1.0), 'ctx_r10': (10, 1.0), 'gate_r20_g0.5': (20, 0.5)}
TRAJ_REPS, ORACLE_LINEAGES = 30, 20
# pi / SFS test
THETA, PI_N0, PI_L, PI_RHO, PI_NSAMP, PI_REPS = 0.01, 500, 2000, 0.4, 50, 20
QS = [1, 2, 5]
# fixation test
FIX_TWO_NS = [-2, 0, 1, 2, 5, 10]
FIX_N, FIX_K = 100, 20000
FIX_Q2 = [(50, 0), (50, 5)]


def slim(script, seed, params):
    args = [SLIM, '-s', str(seed)]
    for k, v in params.items():
        args += ['-d', f"{k}='{v}'" if isinstance(v, str) else f'{k}={v}']
    args.append(f'{PROJ}/scripts/slim/{script}')
    r = subprocess.run(args, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stdout[-2000:] + r.stderr[-2000:])


def traj_task(job):
    cond, rep = job
    rcpg, gate = TRAJ_CONDS[cond]
    out = f'{SIM}/traj_{cond}_rep{rep}'
    slim('e0_neutral.slim', 1000 + rep, dict(TRAJ, RCPG=rcpg, GATE=gate, OUT=out))
    t = pd.read_csv(out + '.traj.tsv', sep='\t')
    anc = open(out + '.anc.txt').read()
    m, se = oracle.expected_cpg_trajectory(anc, TRAJ['ALPHA'], rcpg, gate, t.gen.values, ORACLE_LINEAGES, 5000 + rep)
    t['oracle_mean'], t['oracle_se'], t['cond'], t['rep'] = m, se, cond, rep
    return t


def pi_task(job):
    q, rep = job
    n = PI_N0 // q
    alpha = THETA / (12 * n)                       # u = 3 alpha per site; theta = 4 N u
    out = f'{SIM}/pi_Q{q}_rep{rep}'
    slim('e0_neutral.slim', 2000 + 100 * q + rep,
         dict(N=n, L=PI_L, ALPHA=alpha, RCPG=1, GATE=1.0, RREC=PI_RHO / (4 * n), GENS=10 * n, NCHK=1,
              NSAMP=PI_NSAMP, OUT=out))
    return q, rep, out


def fix_task(job):
    n, two_ns = job
    s = two_ns / (2 * n)
    out = f'{SIM}/fix_N{n}_2Ns{two_ns}.tsv'
    slim('e0_fixation.slim', 3000 + n + two_ns, dict(N=n, K=FIX_K, S=s, H=0.5, OUT=out))
    return out


def seq_summaries(path):
    """pi per site and folded SFS counts (minor-allele count 1..n/2) over biallelic sites of a SLiM sample file."""
    seqs = np.array([list(l.strip()) for l in open(path) if l.strip()])
    n, L = seqs.shape
    pi_sum, sfs = 0.0, np.zeros(n // 2 + 1, int)
    for j in range(L):
        _, c = np.unique(seqs[:, j], return_counts=True)
        pi_sum += (n * n - np.sum(c * c)) / (n * (n - 1))
        if len(c) == 2:
            sfs[min(c)] += 1
    return pi_sum / L, sfs


def kimura(n, two_ns):
    if two_ns == 0:
        return 1 / (2 * n)
    s = two_ns / (2 * n)
    return -np.expm1(-s) / -np.expm1(-two_ns)


def run():
    os.makedirs(SIM, exist_ok=True)
    os.makedirs(CSV, exist_ok=True)
    with Pool(NPROC) as pool:
        fix_jobs = [(FIX_N, x) for x in FIX_TWO_NS] + FIX_Q2
        fix_out = pool.map_async(fix_task, fix_jobs)
        pi_out = pool.map_async(pi_task, [(q, r) for q in QS for r in range(PI_REPS)])
        traj = pool.map(traj_task, [(c, r) for c in TRAJ_CONDS for r in range(TRAJ_REPS)])
        fix_out.get()
        pi_out.get()
    pd.concat(traj).to_csv(f'{CSV}/traj.csv', index=False)


def check():
    rows = []
    # trajectories: per-replicate SLiM - oracle difference, averaged over replicates, at each checkpoint
    tr = pd.read_csv(f'{CSV}/traj.csv')
    tr['d'] = tr.mean_cpg - tr.oracle_mean
    for cond, g in tr.groupby('cond'):
        a = g.groupby('gen').d.agg(['mean', 'std', 'count'])
        a = a[a.index > 1]
        z = a['mean'] / (a['std'] / np.sqrt(a['count']))
        worst = z.abs().idxmax()
        end = a.index.max()
        rows.append(dict(check='traj', case=cond, observed=g[g.gen == end].mean_cpg.mean(),
                         expected=g[g.gen == end].oracle_mean.mean(), z=z[end], note=f'max|z| over 20 times = {abs(z[worst]):.2f} at gen {worst}'))
    # gate vs context-only at the matched effective rate, SLiM vs SLiM
    e = tr[tr.gen == tr.gen.max()].groupby('cond').mean_cpg.agg(['mean', 'std', 'count'])
    a, b = e.loc['gate_r20_g0.5'], e.loc['ctx_r10']
    rows.append(dict(check='traj', case='gate_r20_g0.5 vs ctx_r10 (SLiM)', observed=a['mean'], expected=b['mean'],
                     z=(a['mean'] - b['mean']) / np.sqrt(a['std'] ** 2 / a['count'] + b['std'] ** 2 / b['count']),
                     note='different ancestral seqs; compares endpoint means'))
    # pi and SFS
    pis, sfs_rows = [], []
    for q in QS:
        n = PI_N0 // q
        vals, sfs = [], 0
        for r in range(PI_REPS):
            p, s = seq_summaries(f'{SIM}/pi_Q{q}_rep{r}.sample.txt')
            vals.append(p)
            sfs = sfs + s
            pis.append(dict(Q=q, N=n, rep=r, pi=p))
        vals = np.array(vals)
        exp_pi = THETA / (1 + 4 * THETA / 3)
        rows.append(dict(check='pi', case=f'Q={q} (N={n})', observed=vals.mean(), expected=exp_pi,
                         z=(vals.mean() - exp_pi) / (vals.std(ddof=1) / np.sqrt(len(vals))), note=f'{PI_REPS} reps'))
        k = np.arange(1, PI_NSAMP // 2 + 1)
        w = 1 / k + 1 / (PI_NSAMP - k)
        w[k == PI_NSAMP - k] = 1 / k[k == PI_NSAMP - k]       # i = n/2 counted once
        exp_p = w / w.sum()
        obs = sfs[1:]
        chi2 = np.sum((obs - obs.sum() * exp_p) ** 2 / (obs.sum() * exp_p))
        dof = len(k) - 1
        rows.append(dict(check='sfs', case=f'Q={q} (N={n})', observed=obs[0] / obs.sum(), expected=exp_p[0],
                         z=(chi2 - dof) / np.sqrt(2 * dof), note=f'chi2={chi2:.1f}, dof={dof}; obs/exp = singleton share'))
        for i, (o, ep) in enumerate(zip(obs, exp_p), 1):
            sfs_rows.append(dict(Q=q, i=i, observed=o / obs.sum(), expected=ep, count=o))
    pd.DataFrame(pis).to_csv(f'{CSV}/pi.csv', index=False)
    pd.DataFrame(sfs_rows).to_csv(f'{CSV}/sfs.csv', index=False)
    # fixation
    fx = []
    for n, two_ns in [(FIX_N, x) for x in FIX_TWO_NS] + FIX_Q2:
        f = pd.read_csv(f'{SIM}/fix_N{n}_2Ns{two_ns}.tsv', sep='\t').iloc[0]
        u = kimura(n, two_ns)
        p_hat = f.nfix / f.K
        z = (p_hat - u) / np.sqrt(u * (1 - u) / f.K)
        fx.append(dict(N=n, two_Ns=two_ns, K=f.K, nfix=f.nfix, p_hat=p_hat, kimura=u, z=z))
        rows.append(dict(check='fix', case=f'N={n}, 2Ns={two_ns}', observed=p_hat, expected=u, z=z,
                         note=f'{int(f.nfix)}/{int(f.K)} fixed'))
    pd.DataFrame(fx).to_csv(f'{CSV}/fix.csv', index=False)
    out = pd.DataFrame(rows)
    out['pass'] = out.z.abs() < 3
    out.to_csv(f'{CSV}/e0_checks.csv', index=False)
    with pd.option_context('display.width', 200, 'display.max_colwidth', 60):
        print(out.to_string(index=False, float_format=lambda x: f'{x:.4g}'))


if __name__ == '__main__':
    {'run': run, 'check': check}[sys.argv[1]]()
