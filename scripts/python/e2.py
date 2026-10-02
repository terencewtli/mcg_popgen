"""E2: causal (C) vs passenger (P) methylation under stabilising selection at an intermediate optimum.

  python scripts/python/e2.py run      # stage 1 (M0, P) -> calibrate C's optimum -> stage 2 (C); sim/e2/
  python scripts/python/e2.py check    # statistics -> csv/e2/*.csv
  python scripts/python/e2_plot.py     # -> pdf/e2/

Shared maps (model.py); methylation-gated mutation (r = 20); Q = 2 (N = 100). Readouts:
  P: g = O,               optimum E* = 0.5 (6 of 8 motif positions matching; 28 optimal motif configurations)
  C: g = exp(-R / R0),    R = sum of m_i over the locus, R0 = 3; E* = mean over P-equilibrium haplotypes of
                          exp(-R / R0), so C is optimised for the methylation load P actually has (matching)
  M0: neutral.  KAPPA in {2, 10}: w = exp(-KAPPA (E - E*)^2).
Statistics: (1) trait stability vs sequence turnover on the consensus sequence over time; (2) minor-allele
frequency of segregating variants by effect class (R-only, O-affecting, null); (3) P/D for loss of methylated CpGs.
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

PROJ = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SLIM = os.path.expanduser('~/miniconda3/envs/slim4/bin/slim')
QUICK = os.environ.get('E2_QUICK') == '1'           # pipeline test: 1/10 run length, 2 reps, separate dirs
SIM, CSV = (f'{PROJ}/sim/e2_quick', f'{PROJ}/sim/e2_quick/csv') if QUICK else (f'{PROJ}/sim/e2', f'{PROJ}/csv/e2')
NPROC, R0, NPOP = 9, 3.0, 4
REPS, EVERY = (2, 50) if QUICK else (30, 500)
LAGS = [1, 2, 4, 8, 16, 32, 64]                     # in snapshots of EVERY ticks
NHAP_EFFECT = 20                                    # haplotypes used to average a variant's effect
BASE = replace(model.Params(N=100, ALPHA=2e-5, RCPG=20.0, MODE=1, A1=0.0), BURN=8300 if QUICK else 83000, MEAS=3400 if QUICK else 34000, EVERY=EVERY)
STAGE1 = {'M0': (0, 0.0, 0.0), 'P_k2': (1, 2.0, 0.5), 'P_k10': (1, 10.0, 0.5)}
STAGE2_KAPPA = {'C_k2': 2.0, 'C_k10': 10.0, 'C_k10_Rmatch': 10.0}
# C_k2 / C_k10: E* = mean of exp(-R/R0) over P haplotypes (matched on expression). That under-matches R when C narrows
# the R distribution (Jensen), so C_k10_Rmatch sets E* = exp(-mean_R_P / R0) (matched on mean methylation load).


def calibration():
    f = f'{CSV}/calibration.csv'
    return pd.read_csv(f).set_index('cond') if os.path.exists(f) else None


def cond_spec(cond):
    """(gmode, kappa, estar)"""
    if cond in STAGE1:
        return STAGE1[cond]
    cal = calibration()
    k = STAGE2_KAPPA[cond]
    col = 'estar_Rmatch' if cond.endswith('Rmatch') else 'estar_for_C'
    return 2, k, float(cal.loc[f'P_k{int(k)}', col])


def slim_task(job):
    cond, rep = job
    out = f'{SIM}/{cond}_rep{rep}'
    if os.path.exists(out + f'.pop{NPOP - 1}.txt'):
        return
    gm, k, es = cond_spec(cond)
    P = BASE
    args = [SLIM, '-s', str(20000 + 100 * rep + zlib.crc32(cond.encode()) % 97)]
    d = dict(N=P.N, L=P.L, ALPHA=P.ALPHA, RCPG=P.RCPG, SEL=0, A0=P.A0, A1=P.A1, A2=P.A2, W=P.W, LAMBDA=P.LAMBDA,
             BETA=P.BETA, S0=P.S0, GMODE=gm, KAPPA=k, ESTAR=es, R0=R0, BURN=P.BURN, MEAS=P.MEAS, EVERY=P.EVERY,
             NPOP=NPOP)
    for key, v in d.items():
        args += ['-d', f'{key}={float(v)}' if isinstance(v, float) else f'{key}={v}']
    args += ['-d', f"OUT='{out}'", f'{PROJ}/scripts/slim/e2_stab.slim']
    r = subprocess.run(args, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(cond + r.stdout[-1500:] + r.stderr[-1500:])


def to_int(seq):
    return np.frombuffer(seq.strip().encode(), dtype=np.uint8).copy()


_LUT = np.zeros(256, dtype=np.int64)
for _i, _c in enumerate(b'ACGT'):
    _LUT[_c] = _i


def nuc_arr(seq):
    return _LUT[to_int(seq)]


def load_pop(cond, rep, k):
    return np.array([nuc_arr(l) for l in open(f'{SIM}/{cond}_rep{rep}.pop{k}.txt') if l.strip()])


def run():
    os.makedirs(SIM, exist_ok=True)
    os.makedirs(CSV, exist_ok=True)
    with Pool(NPROC) as pool:
        pool.map(slim_task, [(c, r) for c in STAGE1 for r in range(REPS)], chunksize=1)
    # calibrate C's optimum on the P equilibrium: mean of exp(-R/R0) over all P haplotypes in all snapshots
    rows = []
    for cond in ('P_k2', 'P_k10'):
        Rs = np.concatenate([[model.meth_load(x, BASE) for x in load_pop(cond, r, k)]
                             for r in range(REPS) for k in range(NPOP)])
        rows.append(dict(cond=cond, mean_R=Rs.mean(), sd_R=Rs.std(), estar_for_C=np.exp(-Rs / R0).mean(),
                         estar_Rmatch=np.exp(-Rs.mean() / R0)))
    pd.DataFrame(rows).to_csv(f'{CSV}/calibration.csv', index=False)
    print(pd.DataFrame(rows))
    with Pool(NPROC) as pool:
        pool.map(slim_task, [(c, r) for c in STAGE2_KAPPA for r in range(REPS)], chunksize=1)


# ---------------------------------------------------------------------------------------------------------------
OUTSIDE = None


def outside_motif(p, P):
    return (p < P.mstart - 1) | (p > P.mstart + 7)


def consensus_series(cond, rep, P):
    ts = pd.read_csv(f'{SIM}/{cond}_rep{rep}.ts.tsv', sep='\t')
    out = []
    for _, row in ts.iterrows():
        x = nuc_arr(row.consensus)
        p = model.cpg_pos(x)
        m = model.meth_at(x, p, P)
        out.append(dict(gen=row.gen, mean_E=row.mean_E, mean_w=row.mean_w, x=x, cpg=set(p.tolist()),
                        meth_cpg={int(i) for i, mi in zip(p, m) if mi > 0.5 and outside_motif(i, P)},
                        R=float(m.sum()), O=float(model.occ(x, P)),
                        match=(x[P.mstart:P.mstart + 8] == model.MOTIF)))
    return out


def turnover(series):
    rows = []
    for lag in LAGS:
        a = [(series[t], series[t + lag]) for t in range(len(series) - lag)]
        if not a:
            continue
        jac = [1 - len(u['cpg'] & v['cpg']) / len(u['cpg'] | v['cpg']) if (u['cpg'] | v['cpg']) else 0.0 for u, v in a]
        rows.append(dict(lag=lag, cpg_turnover=np.mean(jac),
                         motif_turnover=np.mean([np.mean(u['match'] != v['match']) for u, v in a]),
                         dR=np.mean([abs(u['R'] - v['R']) for u, v in a]),
                         dO=np.mean([abs(u['O'] - v['O']) for u, v in a])))
    return rows


def variant_rows(pop, P, gmode, kappa, rng):
    rows = []
    n = len(pop)
    for j in range(P.L):
        col = pop[:, j]
        cnt = np.bincount(col, minlength=4)
        if (cnt > 0).sum() < 2:
            continue
        order = np.argsort(-cnt)
        a, b = order[0], order[1]
        carriers = np.flatnonzero(col == a)
        hs = pop[rng.choice(carriers, min(NHAP_EFFECT, len(carriers)), replace=False)]
        dR, dO, dgP, dgC, loss = [], [], [], [], []
        for x in hs:
            y = x.copy()
            y[j] = b
            Rx, Ry = model.meth_load(x, P), model.meth_load(y, P)
            Ox, Oy = model.occ(x, P), model.occ(y, P)
            dR.append(Ry - Rx)
            dO.append(Oy - Ox)
            dgC.append(np.exp(-Ry / R0) - np.exp(-Rx / R0))
            pos_x, pos_y = set(model.cpg_pos(x).tolist()), set(model.cpg_pos(y).tolist())
            loss.append(len(pos_x - pos_y) > 0 and len(pos_y - pos_x) == 0)
        dR, dO = np.mean(dR), np.mean(dO)
        dg = np.mean(dO) if gmode == 1 else np.mean(dgC)
        cls = 'O' if abs(dO) >= 0.05 else ('R_only' if (abs(dO) < 1e-9 and abs(dR) >= 0.2) else
                                          ('null' if (abs(dO) < 1e-9 and abs(dR) < 0.02) else 'other'))
        rows.append(dict(site=j, maf=cnt[b] / n, dR=dR, dO=dO, cls=cls, cpg_loss=bool(np.mean(loss) > 0.5),
                         two_Ns_het=2 * P.N * kappa * (dg / 2) ** 2 if gmode else 0.0,
                         dg_C=np.mean(dgC), dg_P=dO))
    return rows


def run_stats(job):
    cond, rep = job
    P = BASE
    gm, k, es = cond_spec(cond)
    rng = np.random.default_rng(rep)
    ser = consensus_series(cond, rep, P)
    tv = pd.DataFrame(turnover(ser))
    tv['cond'], tv['rep'] = cond, rep
    # divergence-like rate: methylated CpGs (outside motif) present at t and absent at t+1 snapshot
    lost = sum(len(ser[t]['meth_cpg'] - ser[t + 1]['cpg']) for t in range(len(ser) - 1))
    at_risk = sum(len(ser[t]['meth_cpg']) for t in range(len(ser) - 1))
    var, par, poly_loss, meth_cpg_snap = [], [], 0, 0
    ts = pd.read_csv(f'{SIM}/{cond}_rep{rep}.ts.tsv', sep='\t').set_index('gen')
    popt = P.BURN + np.round(np.arange(1, NPOP + 1) * P.MEAS / NPOP).astype(int)
    for kk in range(NPOP):
        pop = load_pop(cond, rep, kk)
        vr = pd.DataFrame(variant_rows(pop, P, gm, k, rng))
        vr['cond'], vr['rep'], vr['snap'] = cond, rep, kk
        var.append(vr)
        # methylated-CpG-loss polymorphism, per methylated CpG of the snapshot consensus (outside motif)
        cons = np.array([np.bincount(pop[:, j], minlength=4).argmax() for j in range(P.L)])
        p = model.cpg_pos(cons)
        mc = p[(model.meth_at(cons, p, P) > 0.5) & outside_motif(p, P)]
        meth_cpg_snap += len(mc)
        sites = set(mc.tolist()) | set((mc + 1).tolist())
        poly_loss += int(((vr.cls == 'R_only') & vr.cpg_loss & vr.site.isin(sites)).sum()) if len(vr) else 0
        if gm:
            py = np.mean([model.readout(x, P, gm, R0) for x in pop])
            par.append(dict(cond=cond, rep=rep, snap=kk, slim_mean_E=ts.loc[popt[kk], 'mean_E'], python_mean_E=py))
    summ = dict(cond=cond, rep=rep,
                R_mean=np.mean([s['R'] for s in ser]), R_sd_time=np.std([s['R'] for s in ser]),
                O_mean=np.mean([s['O'] for s in ser]), O_sd_time=np.std([s['O'] for s in ser]),
                ncpg_mean=np.mean([len(s['cpg']) for s in ser]), mean_w=np.mean([s['mean_w'] for s in ser]),
                mean_E=np.mean([s['mean_E'] for s in ser]),
                loss_rate=lost / max(at_risk, 1) / EVERY, poly_loss_per_cpg=poly_loss / max(meth_cpg_snap, 1))
    return summ, tv, pd.concat(var), par


def boot_se(x, nb=2000, seed=0):
    x = np.asarray(x)
    rng = np.random.default_rng(seed)
    return np.std(x[rng.integers(0, len(x), (nb, len(x)))].mean(1), ddof=1)


def check():
    conds = list(STAGE1) + list(STAGE2_KAPPA)
    with Pool(NPROC) as pool:
        res = pool.map(run_stats, [(c, r) for c in conds for r in range(REPS)], chunksize=2)
    summ = pd.DataFrame([r[0] for r in res])
    tv = pd.concat([r[1] for r in res])
    var = pd.concat([r[2] for r in res])
    par = pd.DataFrame([p for r in res for p in r[3]])
    summ.to_csv(f'{CSV}/runs.csv', index=False)
    tv.to_csv(f'{CSV}/turnover.csv', index=False)
    var.to_csv(f'{CSV}/variants.csv', index=False)
    par['abs_diff'] = (par.slim_mean_E - par.python_mean_E).abs()
    par.to_csv(f'{CSV}/parity.csv', index=False)
    print(f'parity: max |mean_E SLiM - Python| = {par.abs_diff.max():.2e} over {len(par)} snapshots')
    print(calibration())

    # per-condition table
    m0 = summ[summ.cond == 'M0']
    pd_m0 = m0.poly_loss_per_cpg.mean() / m0.loss_rate.mean()
    rows = []
    for c in conds:
        s = summ[summ.cond == c]
        v = var[var.cond == c]
        t64 = tv[(tv.cond == c) & (tv.lag == 64)]
        pdr = s.poly_loss_per_cpg.mean() / s.loss_rate.mean() / pd_m0
        rng = np.random.default_rng(1)
        boots = []
        for _ in range(1000):
            i = rng.integers(0, len(s), len(s))
            j = rng.integers(0, len(m0), len(m0))
            den = s.loss_rate.values[i].mean() * m0.poly_loss_per_cpg.values[j].mean()
            if den > 0:
                boots.append(s.poly_loss_per_cpg.values[i].mean() * m0.loss_rate.values[j].mean() / den)
        row = dict(cond=c, reps=len(s), R_mean=s.R_mean.mean(), R_sd_time=s.R_sd_time.mean(),
                   O_mean=s.O_mean.mean(), O_sd_time=s.O_sd_time.mean(), ncpg=s.ncpg_mean.mean(),
                   mean_E=s.mean_E.mean(), mean_w=s.mean_w.mean(),
                   cpg_turnover_2u=t64.cpg_turnover.mean(), motif_turnover_2u=t64.motif_turnover.mean(),
                   dR_2u=t64.dR.mean(), dO_2u=t64.dO.mean(),
                   PD_methCpG_loss_rel_M0=pdr, PD_se=np.std(boots, ddof=1))
        for cls in ('R_only', 'O', 'null'):
            vc = v[v.cls == cls]
            row[f'maf_{cls}'] = vc.maf.mean()
            row[f'maf_{cls}_se'] = boot_se(vc.groupby('rep').maf.mean().values) if len(vc) else np.nan
            row[f'nseg_{cls}'] = len(vc) / (len(s) * NPOP)
            row[f'2Ns_het_{cls}'] = vc.two_Ns_het.mean()
        rows.append(row)
    tab = pd.DataFrame(rows)
    tab.to_csv(f'{CSV}/summary.csv', index=False)
    with pd.option_context('display.width', 250, 'display.max_columns', 40):
        print(tab.T.to_string(float_format=lambda x: f'{x:.3g}'))


if __name__ == '__main__':
    {'run': run, 'check': check}[sys.argv[1]]()
