"""E2b: E2 at realistic theta (4N u = 0.001 per site, about human; CpG recurrence 4N r alpha ~ 0.007).

  python scripts/python/e2b.py run     # -> sim/e2b/
  python scripts/python/e2b.py check   # -> csv/e2b/*.csv (compared with E2 at theta = 0.024)

Same models, maps, kappa and optima as E2 (scripts/python/e2.py). Composition equilibrates on a 1/u timescale
(~400k generations here), so each run starts from the matching E2 run's final consensus (near that model's
equilibrium composition). Then 20N burn-in for polymorphism and 0.5/u generations of measurement. Composition drift
is checked as first-half vs second-half R. Variants are characterised in SLiM (derived allele known), so frequencies
are derived-allele frequencies (DAF).
"""
import os
import subprocess
import sys
import zlib
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import e2  # noqa: E402
import model  # noqa: E402

PROJ = e2.PROJ
SIM, CSV = f'{PROJ}/sim/e2b', f'{PROJ}/csv/e2b'
REPS, NBG, VEVERY, EVERY = 40, 10, 200, 500
ALPHA = 0.001 / (4 * 100 * 3)                         # theta = 4 N (3 alpha) = 0.001
BURN, MEAS = 2000, 200000
CONDS = list(e2.STAGE1) + list(e2.STAGE2_KAPPA)


def task(job):
    cond, rep = job
    out = f'{SIM}/{cond}_rep{rep}'
    if os.path.exists(out + '.done'):
        return
    src = rep % e2.REPS                                # E2 had 30 reps; reuse ancestors cyclically
    t = pd.read_csv(f'{e2.SIM}/{cond}_rep{src}.ts.tsv', sep='\t')
    fa = out + '.anc.fa'
    open(fa, 'w').write('>anc\n' + t.consensus.iloc[-1] + '\n')
    gm, k, es = e2.cond_spec(cond)
    P = e2.BASE
    d = dict(N=P.N, L=P.L, ALPHA=ALPHA, RCPG=P.RCPG, SEL=0, A0=P.A0, A1=P.A1, A2=P.A2, W=P.W, LAMBDA=P.LAMBDA,
             BETA=P.BETA, S0=P.S0, GMODE=gm, KAPPA=k, ESTAR=es, R0=e2.R0, BURN=BURN, MEAS=MEAS, EVERY=EVERY,
             VEVERY=VEVERY, NBG=NBG)
    args = [e2.SLIM, '-s', str(40000 + 100 * rep + zlib.crc32(cond.encode()) % 97)]
    for key, v in d.items():
        args += ['-d', f'{key}={float(v)}' if isinstance(v, float) else f'{key}={v}']
    args += ['-d', f"ANCFILE='{fa}'", '-d', f"OUT='{out}'", f'{PROJ}/scripts/slim/e2b_lowtheta.slim']
    r = subprocess.run(args, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(cond + r.stdout[-1500:] + r.stderr[-1500:])
    open(out + '.done', 'w').close()


def run():
    os.makedirs(SIM, exist_ok=True)
    os.makedirs(CSV, exist_ok=True)
    with Pool(e2.NPROC) as pool:
        pool.map(task, [(c, r) for c in CONDS for r in range(REPS)], chunksize=1)


def classify(v):
    return np.where(v.dO.abs() >= 0.05, 'O',
                    np.where((v.dO.abs() < 1e-9) & (v.dR.abs() >= 0.2), 'R_only',
                             np.where((v.dO.abs() < 1e-9) & (v.dR.abs() < 0.02), 'null', 'other')))


def run_stats(job):
    cond, rep = job
    P = e2.BASE
    ts = pd.read_csv(f'{SIM}/{cond}_rep{rep}.ts.tsv', sep='\t')
    ser = []
    for c in ts.consensus:
        x = e2.nuc_arr(c)
        p = model.cpg_pos(x)
        m = model.meth_at(x, p, P)
        ser.append(dict(cpg=set(p.tolist()), R=float(m.sum()), O=float(model.occ(x, P)),
                        meth_cpg={int(i) for i, mi in zip(p, m) if mi > 0.5 and e2.outside_motif(i, P)}))
    half = len(ser) // 2
    lost = sum(len(ser[t]['meth_cpg'] - ser[t + 1]['cpg']) for t in range(len(ser) - 1))
    at_risk = sum(len(ser[t]['meth_cpg']) for t in range(len(ser) - 1))
    a, b = ser[0], ser[-1]
    v = pd.read_csv(f'{SIM}/{cond}_rep{rep}.var.tsv', sep='\t')
    nsamp = MEAS // VEVERY
    v['cls'] = classify(v) if len(v) else []
    v['cond'], v['rep'] = cond, rep
    mean_meth_cpg = np.mean([len(s['meth_cpg']) for s in ser])
    poly = ((v.cls == 'R_only') & (v.cpgloss > 0.5)).sum() / nsamp / max(mean_meth_cpg, 1e-9) if len(v) else 0.0
    summ = dict(cond=cond, rep=rep, R_mean=np.mean([s['R'] for s in ser]),
                R_first_half=np.mean([s['R'] for s in ser[:half]]), R_second_half=np.mean([s['R'] for s in ser[half:]]),
                R_sd_time=np.std([s['R'] for s in ser]), O_mean=np.mean([s['O'] for s in ser]),
                cpg_turnover_full=1 - len(a['cpg'] & b['cpg']) / len(a['cpg'] | b['cpg']) if (a['cpg'] | b['cpg']) else 0,
                dR_full=abs(a['R'] - b['R']), loss_rate=lost / max(at_risk, 1) / EVERY, poly_loss_per_cpg=poly,
                nseg_per_sample=len(v) / nsamp)
    return summ, v


def check():
    os.makedirs(CSV, exist_ok=True)
    with Pool(e2.NPROC) as pool:
        res = pool.map(run_stats, [(c, r) for c in CONDS for r in range(REPS)], chunksize=2)
    summ = pd.DataFrame([r[0] for r in res])
    var = pd.concat([r[1] for r in res if len(r[1])])
    summ.to_csv(f'{CSV}/runs.csv', index=False)
    var.to_csv(f'{CSV}/variants.csv', index=False)
    m0 = summ[summ.cond == 'M0']
    rows = []
    rng = np.random.default_rng(0)
    for c in CONDS:
        s = summ[summ.cond == c]
        vv = var[var.cond == c]
        pdr = (s.poly_loss_per_cpg.mean() / s.loss_rate.mean()) / (m0.poly_loss_per_cpg.mean() / m0.loss_rate.mean())
        boots = []
        for _ in range(1000):
            i, j = rng.integers(0, len(s), len(s)), rng.integers(0, len(m0), len(m0))
            den = s.loss_rate.values[i].mean() * m0.poly_loss_per_cpg.values[j].mean()
            if den > 0:
                boots.append(s.poly_loss_per_cpg.values[i].mean() * m0.loss_rate.values[j].mean() / den)
        row = dict(cond=c, reps=len(s), R_mean=s.R_mean.mean(), R_drift=(s.R_second_half - s.R_first_half).mean(),
                   R_drift_se=(s.R_second_half - s.R_first_half).sem(), R_sd_time=s.R_sd_time.mean(),
                   O_mean=s.O_mean.mean(), cpg_turnover_full=s.cpg_turnover_full.mean(), dR_full=s.dR_full.mean(),
                   nseg_per_sample=s.nseg_per_sample.mean(), PD_methCpG_loss_rel_M0=pdr, PD_se=np.std(boots, ddof=1))
        for cls in ('R_only', 'O', 'null'):
            vc = vv[vv.cls == cls]
            per_rep = vc.groupby('rep').daf.mean()
            row[f'daf_{cls}'] = vc.daf.mean()
            row[f'daf_{cls}_se'] = e2.boot_se(per_rep.values) if len(per_rep) > 1 else np.nan
            row[f'rare_frac_{cls}'] = (vc.daf < 0.05).mean() if len(vc) else np.nan
            row[f'n_rows_{cls}'] = len(vc)
        rows.append(row)
    tab = pd.DataFrame(rows)
    tab.to_csv(f'{CSV}/summary.csv', index=False)
    hi = pd.read_csv(f'{e2.CSV}/summary.csv').set_index('cond')
    comp = pd.DataFrame(dict(PD_theta0024=hi.PD_methCpG_loss_rel_M0, PD_theta0001=tab.set_index('cond').PD_methCpG_loss_rel_M0,
                             PD_theta0001_se=tab.set_index('cond').PD_se,
                             R_theta0024=hi.R_mean, R_theta0001=tab.set_index('cond').R_mean))
    comp.to_csv(f'{CSV}/compare_theta.csv')
    with pd.option_context('display.width', 250, 'display.max_columns', 40):
        print(tab.T.to_string(float_format=lambda x: f'{x:.3g}'))
        print(comp.to_string(float_format=lambda x: f'{x:.3g}'))


if __name__ == '__main__':
    {'run': run, 'check': check}[sys.argv[1]]()
