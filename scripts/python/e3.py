"""E3: TF and methylation both matter (realistic theta, as E2b).

  python scripts/python/e3.py run      # -> sim/e3/
  python scripts/python/e3.py check    # -> csv/e3/*.csv
  python scripts/python/e3_plot.py     # -> pdf/e3/

Conditions (methylation-gated mutation, r = 20; N = 100; high-theta burn-in 5/u at alpha = 2e-5, then realistic
theta 4Nu = 0.001 for 20N + 0.5/u generations; see scripts/slim/e3_combined.slim):
  M0        neutral (same protocol)
  CP_k10/40 g = O exp(-R/R0), motif TGACTCAT, E* = 0.5 exp(-mean R_P / R0) (the E2 passenger equilibrium is one
            optimal point among many O-R trade-offs)
  Pcre_k10  g = O, motif CRE TGACGTCA with its CpG required for binding, E* = O(7 matches)   [control for MS]
  MS_k10    g = O (1 - m_c), same motif, E* = O7 (1 - m_c(O7))   (methylation-sensitive TF)
The CRE optima sit at 7 of 8 matches so that a homozygote can reach them. At 6 matches the discrete MS readout makes
the 6/7 heterozygote fittest (overdominance), which would contaminate the frequency statistics.
"""
import os
import subprocess
import sys
import zlib
from dataclasses import replace
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import e2  # noqa: E402
import e2b  # noqa: E402
import model  # noqa: E402

PROJ = e2.PROJ
SIM, CSV = f'{PROJ}/sim/e3', f'{PROJ}/csv/e3'
REPS, NBG, VEVERY, EVERY = 40, 10, 200, 500
ALPHA_HI, ALPHA_LO = 2e-5, e2b.ALPHA
BURN_HI, BURN_LO, MEAS = 83000, 2000, 200000
P0 = model.Params()


def _o(k):
    return float(model.occ_from_matches(k, P0))


R_P = float(pd.read_csv(f'{e2.CSV}/calibration.csv').set_index('cond').loc['P_k10', 'mean_R'])
O7 = _o(7)
MC7 = 1 / (1 + np.exp(-(P0.A0 - P0.A2 * O7)))          # motif-CpG methylation at 7 matches (d = 0)
CONDS = {  # name: (GMODE, MOTIFID, KAPPA, ESTAR)
    'M0': (1, 0, 0.0, 0.5),
    'CP_k10': (3, 0, 10.0, 0.5 * np.exp(-R_P / e2.R0)),
    'CP_k40': (3, 0, 40.0, 0.5 * np.exp(-R_P / e2.R0)),
    'Pcre_k10': (1, 1, 10.0, O7),
    'MS_k10': (4, 1, 10.0, O7 * (1 - MC7)),
}


def task(job):
    cond, rep = job
    out = f'{SIM}/{cond}_rep{rep}'
    if os.path.exists(out + '.done'):
        return
    gm, mid, k, es = CONDS[cond]
    P = P0
    d = dict(N=P.N, L=P.L, ALPHA_HI=ALPHA_HI, ALPHA_LO=ALPHA_LO, RCPG=P.RCPG, A0=P.A0, A1=P.A1, A2=P.A2, W=P.W,
             LAMBDA=P.LAMBDA, BETA=P.BETA, S0=P.S0, GMODE=gm, MOTIFID=mid, KAPPA=k, ESTAR=float(es), R0=e2.R0,
             BURN_HI=BURN_HI, BURN_LO=BURN_LO, MEAS=MEAS, EVERY=EVERY, VEVERY=VEVERY, NBG=NBG)
    args = [e2.SLIM, '-s', str(60000 + 100 * rep + zlib.crc32(cond.encode()) % 97)]
    for key, v in d.items():
        args += ['-d', f'{key}={float(v)}' if isinstance(v, float) else f'{key}={v}']
    args += ['-d', f"OUT='{out}'", f'{PROJ}/scripts/slim/e3_combined.slim']
    r = subprocess.run(args, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(cond + r.stdout[-1500:] + r.stderr[-1500:])
    open(out + '.done', 'w').close()


def run():
    os.makedirs(SIM, exist_ok=True)
    os.makedirs(CSV, exist_ok=True)
    order = sorted(CONDS, key=lambda c: CONDS[c][2] == 0)          # selected (slower) first
    with Pool(e2.NPROC) as pool:
        pool.map(task, [(c, r) for c in order for r in range(REPS)], chunksize=1)


def classify(v, mid, P):
    mcpg = v.pos.isin([P.mstart + 3, P.mstart + 4]) if mid == 1 else np.zeros(len(v), bool)
    return np.where(mcpg, 'motif_CpG',
                    np.where(v.dO.abs() >= 0.05, 'motif',
                             np.where((v.dO.abs() < 1e-9) & (v.dR.abs() >= 0.2), 'R_only',
                                      np.where((v.dO.abs() < 1e-9) & (v.dR.abs() < 0.02), 'null', 'other'))))


def run_stats(job):
    cond, rep = job
    gm, mid, k, es = CONDS[cond]
    P = P0
    motif = model.MOTIFS[mid]
    ts = pd.read_csv(f'{SIM}/{cond}_rep{rep}.ts.tsv', sep='\t')
    ser, par = [], 0.0
    for c, g in zip(ts.consensus, ts.cons_g):
        x = e2.nuc_arr(c)
        p = model.cpg_pos(x)
        m = model.meth_at_m(x, p, P, motif)
        par = max(par, abs(model.readout_e3(x, P, gm, mid, e2.R0) - g))
        ser.append(dict(cpg=set(p.tolist()), R=float(m.sum()), O=model.occ_m(x, P, motif),
                        match=(x[P.mstart:P.mstart + 8] == motif),
                        meth_cpg={int(i) for i, mi in zip(p, m) if mi > 0.5 and e2.outside_motif(i, P)}))
    lost = sum(len(ser[t]['meth_cpg'] - ser[t + 1]['cpg']) for t in range(len(ser) - 1))
    at_risk = sum(len(ser[t]['meth_cpg']) for t in range(len(ser) - 1))
    a, b = ser[0], ser[-1]
    v = pd.read_csv(f'{SIM}/{cond}_rep{rep}.var.tsv', sep='\t')
    v['cls'] = classify(v, mid, P) if len(v) else []
    v['cond'], v['rep'] = cond, rep
    v['two_Ns_het'] = 2 * P.N * k * (v.dg / 2) ** 2 if len(v) else []
    nsamp = MEAS // VEVERY
    mean_meth_cpg = np.mean([len(s['meth_cpg']) for s in ser])
    poly = ((v.cls == 'R_only') & (v.cpgloss > 0.5)).sum() / nsamp / max(mean_meth_cpg, 1e-9) if len(v) else 0.0
    Os, Rs = np.array([s['O'] for s in ser]), np.array([s['R'] for s in ser])
    summ = dict(cond=cond, rep=rep, parity=par, mean_E=ts.mean_E.mean(), mean_w=ts.mean_w.mean(),
                O_mean=Os.mean(), R_mean=Rs.mean(), R_sd_time=Rs.std(), O_sd_time=Os.std(),
                corr_O_R=np.corrcoef(Os, Rs)[0, 1] if Os.std() > 0 and Rs.std() > 0 else np.nan,
                motif_matches=np.mean([s['match'].sum() for s in ser]),
                motif_turnover_full=np.mean(a['match'] != b['match']),
                motif_cpg_present=np.mean([(P.mstart + 3) in s['cpg'] for s in ser]),
                cpg_turnover_full=1 - len(a['cpg'] & b['cpg']) / len(a['cpg'] | b['cpg']) if (a['cpg'] | b['cpg']) else 0,
                dR_full=abs(a['R'] - b['R']), loss_rate=lost / max(at_risk, 1) / EVERY, poly_loss_per_cpg=poly)
    return summ, v


def check():
    os.makedirs(CSV, exist_ok=True)
    with Pool(e2.NPROC) as pool:
        res = pool.map(run_stats, [(c, r) for c in CONDS for r in range(REPS)], chunksize=2)
    summ = pd.DataFrame([r[0] for r in res])
    var = pd.concat([r[1] for r in res if len(r[1])])
    summ.to_csv(f'{CSV}/runs.csv', index=False)
    var.to_csv(f'{CSV}/variants.csv', index=False)
    print(f'parity: max |readout(consensus) SLiM - Python| = {summ.parity.max():.2e}')
    m0 = summ[summ.cond == 'M0']
    vm0 = var[var.cond == 'M0']
    rng = np.random.default_rng(0)
    rows = []
    for c in CONDS:
        s, vv = summ[summ.cond == c], var[var.cond == c]
        pdr = (s.poly_loss_per_cpg.mean() / s.loss_rate.mean()) / (m0.poly_loss_per_cpg.mean() / m0.loss_rate.mean())
        boots = []
        for _ in range(1000):
            i, j = rng.integers(0, len(s), len(s)), rng.integers(0, len(m0), len(m0))
            den = s.loss_rate.values[i].mean() * m0.poly_loss_per_cpg.values[j].mean()
            if den > 0:
                boots.append(s.poly_loss_per_cpg.values[i].mean() * m0.loss_rate.values[j].mean() / den)
        row = dict(cond=c, reps=len(s), estar=CONDS[c][3], mean_E=s.mean_E.mean(), mean_w=s.mean_w.mean(),
                   O_mean=s.O_mean.mean(), R_mean=s.R_mean.mean(), O_sd_time=s.O_sd_time.mean(),
                   R_sd_time=s.R_sd_time.mean(), corr_O_R=s.corr_O_R.mean(), motif_matches=s.motif_matches.mean(),
                   motif_cpg_present=s.motif_cpg_present.mean(), motif_turnover_full=s.motif_turnover_full.mean(),
                   cpg_turnover_full=s.cpg_turnover_full.mean(), dR_full=s.dR_full.mean(),
                   PD_methCpG_loss_rel_M0=pdr, PD_se=np.std(boots, ddof=1))
        for cls in ('R_only', 'motif', 'motif_CpG', 'null'):
            vc = vv[vv.cls == cls]
            base = vm0[vm0.cls == ('motif' if cls == 'motif_CpG' else cls)].daf.mean()
            per_rep = vc.groupby('rep').daf.mean()
            row[f'daf_{cls}'] = vc.daf.mean()
            row[f'daf_{cls}_se'] = e2.boot_se(per_rep.values) if len(per_rep) > 1 else np.nan
            row[f'daf_{cls}_rel_M0'] = vc.daf.mean() / base if len(vc) else np.nan
            row[f'n_rows_{cls}'] = len(vc)
            row[f'2Ns_het_{cls}'] = vc.two_Ns_het.mean() if len(vc) else np.nan
        rows.append(row)
    tab = pd.DataFrame(rows)
    tab.to_csv(f'{CSV}/summary.csv', index=False)
    with pd.option_context('display.width', 250, 'display.max_columns', 60):
        print(tab.set_index('cond').T.to_string(float_format=lambda x: f'{x:.3g}'))


if __name__ == '__main__':
    {'run': run, 'check': check}[sys.argv[1]]()
