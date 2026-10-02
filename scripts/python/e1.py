"""E1: how much CpG "conservation" does passenger methylation buy for free?

  python scripts/python/e1.py run      # SLiM + oracle -> sim/e1/ (raw), csv/e1/ (per-run profiles)
  python scripts/python/e1.py check    # parity, SLiM vs oracle, rescaling, summaries -> csv/e1/*.csv
  python scripts/python/e1_plot.py     # -> pdf/e1/

No condition puts methylation (or any CpG) into fitness. Selection, where present, acts on TF occupancy only.
Conditions (Q = 2 baseline: N = 100, alpha = 2e-5 per target, r = 20, burn-in 5/u, measurement 2/u, u = 3 alpha):
  M0_ctx      neutral, context-only CpG hypermutation
  M0_meth     neutral, methylation-gated (motif decays, so TF protection disappears)
  M0_meth_fb  neutral, gated, + CpG-density feedback on methylation (A1 = 0.5): can islands self-sustain?
  P_ctx       selection on occupancy (2Ns = 10 for one heterozygous mismatch), context-only
  P_meth      same selection, gated: the passenger-methylation footprint
  P_meth_fb   same, + density feedback
  sweeps on P_meth: LAMBDA in {10, 50}, r in {5, 50}; rescaling: P_meth at Q = 1, 4
Oracle (scripts/python/oracle.py e1_lineage) is exact for M0_* and freezes the motif for P_*.
"""
import glob
import os
import subprocess
import sys
import zlib
from dataclasses import replace
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import model  # noqa: E402
import oracle  # noqa: E402

PROJ = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SLIM = os.path.expanduser('~/miniconda3/envs/slim4/bin/slim')
QUICK = os.environ.get('E1_QUICK') == '1'          # pipeline test: 1/10 run length, few reps, separate dirs
SIM, CSV = (f'{PROJ}/sim/e1_quick', f'{PROJ}/sim/e1_quick/csv') if QUICK else (f'{PROJ}/sim/e1', f'{PROJ}/csv/e1')
NPROC = 9
TWO_NS = 10.0
REPS_MAIN, REPS_SWEEP, ORACLE_LINEAGES = (3, 2, 20) if QUICK else (40, 20, 200)
BIN = 10


def base_params(q=2):
    """Q-rescaled parameters. Q = 2 is the baseline; Q = 1 is N = 200."""
    n = 200 // q
    alpha = 1e-5 * q
    u = 3 * alpha
    f = 10 if QUICK else 1
    return model.Params(N=n, ALPHA=alpha, BURN=int(round(5 / u / f)), MEAS=int(round(2 / u / f)), EVERY=1000 // q // f)


def conditions():
    b = base_params(2)
    sel = model.sel_for_two_ns(TWO_NS, b.N, b)
    c = {
        'M0_ctx': replace(b, MODE=0),
        'M0_meth': replace(b, MODE=1),
        'M0_meth_fb': replace(b, MODE=1, A1=0.5),
        'P_ctx': replace(b, MODE=0, SEL=sel),
        'P_meth': replace(b, MODE=1, SEL=sel),
        'P_meth_fb': replace(b, MODE=1, SEL=sel, A1=0.5),
        'P_meth_lam10': replace(b, MODE=1, SEL=sel, LAMBDA=10.0),
        'P_meth_lam50': replace(b, MODE=1, SEL=sel, LAMBDA=50.0),
        'P_meth_r5': replace(b, MODE=1, SEL=sel, RCPG=5.0),
        'P_meth_r50': replace(b, MODE=1, SEL=sel, RCPG=50.0),
    }
    for q in (1, 4):
        bq = base_params(q)
        c[f'P_meth_Q{q}'] = replace(bq, MODE=1, SEL=model.sel_for_two_ns(TWO_NS, bq.N, bq))
    return c


def nreps(cond):
    return REPS_MAIN if cond in ('M0_ctx', 'M0_meth', 'M0_meth_fb', 'P_ctx', 'P_meth', 'P_meth_fb') else REPS_SWEEP


def slim_task(job):
    cond, rep = job
    P = conditions()[cond]
    out = f'{SIM}/{cond}_rep{rep}'
    if os.path.exists(out + '.parity.tsv'):
        return
    args = [SLIM, '-s', str(10000 + 100 * rep + zlib.crc32(cond.encode()) % 97)]
    for k, v in P.slim_args().items():
        args += ['-d', f'{k}={float(v)}' if isinstance(v, float) else f'{k}={v}']
    args += ['-d', f"OUT='{out}'", f'{PROJ}/scripts/slim/e1_passenger.slim']
    r = subprocess.run(args, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(cond + r.stdout[-1500:] + r.stderr[-1500:])


def oracle_task(job):
    cond, chunk, n = job
    out = f'{SIM}/oracle_{cond}_chunk{chunk}.npz'
    if os.path.exists(out):
        return
    P = conditions()[cond]
    rng = np.random.default_rng(777 + 1000 * chunk + zlib.crc32(cond.encode()) % 997)
    times = np.arange(P.BURN + P.EVERY, P.BURN + P.MEAS + 1, P.EVERY, dtype=float)
    tot = None
    for _ in range(n):
        anc = rng.integers(0, 4, P.L)
        anc[P.mstart:P.mstart + 8] = model.MOTIF
        a = oracle.e1_lineage(anc, P, times, rng, freeze_motif=cond.startswith('P_'))
        tot = a if tot is None else {k: tot[k] + a[k] for k in a}
    np.savez(out, **tot)


def run():
    os.makedirs(SIM, exist_ok=True)
    os.makedirs(CSV, exist_ok=True)
    conds = conditions()
    # longest jobs first
    order = sorted(conds, key=lambda c: -(conds[c].N * (conds[c].BURN + conds[c].MEAS) * (1 + 6 * (conds[c].SEL > 0))))
    slim_jobs = [(c, r) for c in order for r in range(nreps(c))]
    oracle_jobs = [(c, k, ORACLE_LINEAGES // 10) for c in conds for k in range(10)]
    with Pool(NPROC) as pool:
        a = pool.map_async(oracle_task, oracle_jobs)
        pool.map(slim_task, slim_jobs, chunksize=1)
        a.get()


# ---------------------------------------------------------------------------------------------------------------
def binned(prof, P):
    """prof: dict/DataFrame with per-position cpg, c, g, msum and scalar n. Returns per-distance-bin stats."""
    L = P.L
    pos = np.arange(L - 1)
    d = model.dist_to_motif(pos, P)
    keep = (pos < P.mstart - 1) | (pos > P.mstart + 7)          # CpGs overlapping the motif excluded
    n = prof['n']
    cpg = np.asarray(prof['cpg'])[:L - 1]
    exp_cpg = np.asarray(prof['c'])[:L - 1] * np.asarray(prof['g'])[1:] / n
    msum = np.asarray(prof['msum'])[:L - 1]
    b = (d // BIN).astype(int)
    rows = []
    for k in np.unique(b[keep]):
        sel = keep & (b == k)
        rows.append(dict(bin=k * BIN, density=cpg[sel].sum() / (n * sel.sum()) * 100,       # CpGs per 100 bp
                         oe=cpg[sel].sum() / max(exp_cpg[sel].sum(), 1e-12),
                         meth=msum[sel].sum() / cpg[sel].sum() if cpg[sel].sum() else np.nan))
    return pd.DataFrame(rows)


def island_stats(prof, P):
    pos = np.arange(P.L - 1)
    d = model.dist_to_motif(pos, P)
    keep = (pos < P.mstart - 1) | (pos > P.mstart + 7)
    cpg = np.asarray(prof['cpg'])[:P.L - 1]
    near, far = keep & (d < 25), keep & (d >= 75)
    dn, df = cpg[near].sum() / near.sum(), cpg[far].sum() / far.sum()
    return dict(density_near=dn / prof['n'] * 100, density_far=df / prof['n'] * 100,
                contrast=dn / df if df > 0 else np.nan)


def load_slim(cond, rep):
    f = f'{SIM}/{cond}_rep{rep}'
    t = pd.read_csv(f + '.prof.tsv', sep='\t')
    s = pd.read_csv(f + '.scal.tsv', sep='\t').iloc[0]
    prof = dict(cpg=t.cpg.values, c=t.c.values, g=t.g.values, msum=t.msum.values, n=int(s.n))
    return prof, s


def parity(cond, rep, P):
    f = f'{SIM}/{cond}_rep{rep}'
    seqs = [np.array(['ACGT'.index(ch) for ch in l.strip()]) for l in open(f + '.seqs.txt') if l.strip()]
    par = pd.read_csv(f + '.parity.tsv', sep='\t')
    worst = 0.0
    for h, g in par.groupby('hap'):
        nuc = seqs[h]
        p = model.cpg_pos(nuc)
        assert np.array_equal(p, g.pos.values), f'CpG positions differ {cond} rep{rep} hap{h}'
        worst = max(worst, np.max(np.abs(model.meth_at(nuc, p, P) - g.m.values)), abs(model.occ(nuc, P) - g.occ.iloc[0]))
    return worst


def boot_ratio_se(a, b, nboot=2000, seed=0):
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(a), (nboot, len(a)))
    return np.std(a[idx].mean(1) / b[idx].mean(1), ddof=1)


def check():
    conds = conditions()
    prof_rows, summ, par_rows = [], [], []
    for cond, P in conds.items():
        reps = [r for r in range(nreps(cond)) if os.path.exists(f'{SIM}/{cond}_rep{r}.scal.tsv')]
        isl, occs, mism = [], [], []
        for r in reps:
            prof, s = load_slim(cond, r)
            b = binned(prof, P)
            b['cond'], b['source'], b['rep'] = cond, 'slim', r
            prof_rows.append(b)
            isl.append(island_stats(prof, P))
            occs.append(s.osum / s.n)
            mism.append(1 - s.s8 / s.n)
            par_rows.append(dict(cond=cond, rep=r, max_abs_diff=parity(cond, r, P)))
        tot = None
        for f in sorted(glob.glob(f'{SIM}/oracle_{cond}_chunk*.npz')):
            z = dict(np.load(f))
            tot = z if tot is None else {k: tot[k] + z[k] for k in z}
        tot['n'] = int(tot['n'])
        b = binned(tot, P)
        b['cond'], b['source'], b['rep'] = cond, 'oracle', -1
        prof_rows.append(b)
        io = island_stats(tot, P)
        isl = pd.DataFrame(isl)
        summ.append(dict(cond=cond, N=P.N, alpha=P.ALPHA, r=P.RCPG, mode=P.MODE, A1=P.A1, LAMBDA=P.LAMBDA,
                         two_Ns_8_7=model.selection_table(P)['8/7'] if P.SEL > 0 else 0.0, reps=len(reps),
                         occ=np.mean(occs), frac_motif_imperfect=np.mean(mism),
                         density_near=isl.density_near.mean(), density_near_se=isl.density_near.sem(),
                         density_far=isl.density_far.mean(), density_far_se=isl.density_far.sem(),
                         # ratio of pooled means (same estimator as the oracle); SE by bootstrap over replicates
                         contrast=isl.density_near.mean() / isl.density_far.mean(),
                         contrast_se=boot_ratio_se(isl.density_near.values, isl.density_far.values),
                         oracle_density_near=io['density_near'], oracle_density_far=io['density_far'],
                         oracle_contrast=io['contrast'], oracle_lineages=ORACLE_LINEAGES))
    pd.concat(prof_rows).to_csv(f'{CSV}/profiles.csv', index=False)
    sm = pd.DataFrame(summ)
    # SLiM vs oracle, near and far density; the oracle SE is small (200 lineages x 66 samples) and ignored
    sm['z_near_vs_oracle'] = (sm.density_near - sm.oracle_density_near) / sm.density_near_se
    sm['z_far_vs_oracle'] = (sm.density_far - sm.oracle_density_far) / sm.density_far_se
    sm.to_csv(f'{CSV}/summary.csv', index=False)
    par = pd.DataFrame(par_rows)
    par.to_csv(f'{CSV}/parity.csv', index=False)
    q = sm[sm.cond.isin(['P_meth_Q1', 'P_meth', 'P_meth_Q4'])].copy()
    ref = q[q.cond == 'P_meth'].iloc[0]
    q['z_contrast_vs_Q2'] = (q.contrast - ref.contrast) / np.sqrt(q.contrast_se ** 2 + ref.contrast_se ** 2)
    q['z_near_vs_Q2'] = (q.density_near - ref.density_near) / np.sqrt(q.density_near_se ** 2 + ref.density_near_se ** 2)
    q.to_csv(f'{CSV}/qcheck.csv', index=False)
    print(f'parity: max |m_SLiM - m_python| over {len(par)} runs = {par.max_abs_diff.max():.2e}')
    cols = ['cond', 'reps', 'occ', 'frac_motif_imperfect', 'density_near', 'density_far', 'contrast', 'contrast_se',
            'oracle_contrast', 'z_near_vs_oracle', 'z_far_vs_oracle']
    with pd.option_context('display.width', 250):
        print(sm[cols].to_string(index=False, float_format=lambda x: f'{x:.3g}'))
        print(q[['cond', 'N', 'contrast', 'contrast_se', 'density_near', 'z_contrast_vs_Q2', 'z_near_vs_Q2']]
              .to_string(index=False, float_format=lambda x: f'{x:.3g}'))


if __name__ == '__main__':
    {'run': run, 'check': check}[sys.argv[1]]()
