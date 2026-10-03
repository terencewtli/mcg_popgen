"""R2c: germline-matched comparator (R1) and the R2 discriminator (P/D vs somatic methylation), genome-wide.

CpG side: R2b tables (ancestral CpGs with events; R1e definitions), CpGs in strong HOMER hits (R1g margin >= 2)
excluded (F6). Non-CpG side: R1f counts run with --human-sperm / --loyfer (S>W changes).

Germline bins for both sides use the human regional sperm m (+-500 bp, cov >= 5). The CpG side additionally
requires the chimp site m in the same bin ('conserved'), as in R1e.

Tables (P/D with 95% CIs from a joint 1 Mb block bootstrap, so CpG and non-CpG numbers are resampled together):
  germ_matched   per germline bin x element: CpG P/D relative to its background, S>W relative to its background,
                 and their ratio (excess) = CpG_rel / SW_rel. Excess > 1: CpG loss more constrained than other S>W
                 changes in the same element and germline context.
  soma_<var>     per germline bin x element group x somatic bin (<var> = soma_med / soma_min, regional): CpG P/D,
                 S>W P/D and ratio = CpG / SW. Causal methylation predicts the ratio rising with somatic m at fixed
                 germline m; passenger predicts it flat. 'trend' rows give ratio(top bin) / ratio(bottom bin).
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from R1e_pd_tables import events  # noqa: E402
from R2b_cpg_soma import GERM_BINS, GERM_LABELS, SOMA_BINS, SOMA_LABELS  # noqa: E402


def load_cpg(soma_glob: str, motif_dir: str, strong: float, min_cov: int) -> pd.DataFrame:
    parts = []
    for f in sorted(glob.glob(soma_glob)):
        x = pd.read_csv(f, sep='\t', dtype={'snv_pos': str, 'snv_kind': str, 'cons_C': str, 'cons_G': str})
        chrom = f.split('/')[-1].split('.')[0]
        mo = pd.read_csv(f'{motif_dir}/{chrom}.motif.tsv.gz', sep='\t', usecols=['pos', 'best_margin'])
        x = x.merge(mo, on='pos', how='left')
        x = x[~(x.best_margin >= strong)]
        x['chrom'] = chrom
        parts.append(x)
    d = pd.concat(parts, ignore_index=True)
    d = d[d.state != 'other']
    d = pd.concat([d, events(d)], axis=1)
    d = d[(d.P + d.D) > 0]
    d['block'] = d.chrom + ':' + (d.pos // 1_000_000).astype(str)
    h_bin = pd.cut(d.h_reg, GERM_BINS, labels=GERM_LABELS).astype(str)
    c_bin = np.where(d.c_cov >= min_cov, pd.cut(d.c_m, GERM_BINS, labels=GERM_LABELS).astype(str), 'nan')
    d['germ'] = np.where(h_bin == c_bin, h_bin, 'discordant')
    d['element'] = np.where(d.cons != '.', 'coding', np.where(d.cgi == 1, 'CGI', np.where(d.ccre != 'none', d.ccre, 'background')))
    d['soma_med'] = pd.cut(d.soma_med_reg, SOMA_BINS, labels=SOMA_LABELS).astype(str)
    d['soma_min'] = pd.cut(d.soma_min_reg, SOMA_BINS, labels=SOMA_LABELS).astype(str)
    return d


def block_counts(df: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    return df.groupby(keys + ['block'])[['P', 'D']].sum()


def boot(cpg: pd.DataFrame, sw: pd.DataFrame, stat, n_boot: int, rng: np.random.Generator) -> tuple[float, float, float]:
    """stat(cpg_block_counts, sw_block_counts) -> float; joint resampling of 1 Mb blocks."""
    blocks = np.union1d(cpg.index.get_level_values('block').unique(), sw.index.get_level_values('block').unique())
    est = stat(cpg, sw)
    cb = cpg.reset_index().pivot_table(index='block', columns=[c for c in cpg.index.names if c != 'block'],
                                       values=['P', 'D'], aggfunc='sum', fill_value=0).reindex(blocks, fill_value=0)
    sb = sw.reset_index().pivot_table(index='block', columns=[c for c in sw.index.names if c != 'block'],
                                      values=['P', 'D'], aggfunc='sum', fill_value=0).reindex(blocks, fill_value=0)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(blocks), len(blocks))
        c = cb.iloc[idx].sum().unstack(0)
        s = sb.iloc[idx].sum().unstack(0)
        vals.append(stat(c, s))
    lo, hi = np.nanpercentile(vals, [2.5, 97.5])
    return est, lo, hi


def pdr(t: pd.DataFrame, key) -> float:
    try:
        r = t.loc[key]
    except KeyError:
        return np.nan
    if isinstance(r, pd.DataFrame):
        r = r.sum()
    return r.P / r.D if r.D > 0 else np.nan


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--soma-glob', required=True)
    ap.add_argument('--motif-dir', required=True)
    ap.add_argument('--sw-glob', required=True, help='R1f counts with germ / soma columns')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--strong', type=float, default=2.0)
    ap.add_argument('--min-cov', type=int, default=5)
    ap.add_argument('--n-boot', type=int, default=200)
    a = ap.parse_args()
    rng = np.random.default_rng(1)

    cpg = load_cpg(a.soma_glob, a.motif_dir, a.strong, a.min_cov)
    sw = pd.concat([pd.read_csv(f, sep='\t') for f in sorted(glob.glob(a.sw_glob))], ignore_index=True)
    sw = sw[sw.direction == 'S>W']
    print(f'CpG events: P {cpg.P.sum():,} D {cpg.D.sum():,}; non-CpG S>W: P {sw.P.sum():,} D {sw.D.sum():,}', flush=True)

    # germline-matched comparator
    rows = []
    c_bc, s_bc = block_counts(cpg, ['germ', 'element']), block_counts(sw, ['germ', 'element'])
    for g in GERM_LABELS:
        for el in sorted(set(cpg.element) - {'background'}):
            def stat(c, s, g=g, el=el):
                c, s = c.groupby(level=[0, 1]).sum() if 'block' in c.index.names else c, \
                       s.groupby(level=[0, 1]).sum() if 'block' in s.index.names else s
                return (pdr(c, (g, el)) / pdr(c, (g, 'background'))) / (pdr(s, (g, el)) / pdr(s, (g, 'background')))
            cs = cpg[(cpg.germ == g) & (cpg.element == el)]
            ss = sw[(sw.germ == g) & (sw.element == el)]
            if cs.D.sum() < 20 or ss.D.sum() < 20:
                continue
            est, lo, hi = boot(c_bc.loc[[g]], s_bc.loc[[g]], stat, a.n_boot, rng)
            cb, sb = cpg[(cpg.germ == g) & (cpg.element == 'background')], sw[(sw.germ == g) & (sw.element == 'background')]
            rows.append(dict(germ=g, element=el, cpg_n_events=len(cs), cpg_PD=cs.P.sum() / cs.D.sum(),
                             cpg_rel=(cs.P.sum() / cs.D.sum()) / (cb.P.sum() / cb.D.sum()),
                             sw_PD=ss.P.sum() / ss.D.sum(), sw_rel=(ss.P.sum() / ss.D.sum()) / (sb.P.sum() / sb.D.sum()),
                             excess=est, excess_lo=lo, excess_hi=hi))
    t = pd.DataFrame(rows)
    t.to_csv(f'{a.outdir}/R2c_germ_matched.tsv', sep='\t', index=False, float_format='%.4g')
    print('\n== germline-matched comparator (excess = CpG_rel / SW_rel)\n' + t.to_string(index=False, float_format=lambda v: f'{v:.3g}'), flush=True)

    # R2: ratio CpG / SW across somatic bins, at fixed germline bin, non-coding only
    groups = {'noncoding_all': lambda e: e != 'coding', 'elements': lambda e: ~e.isin(['coding', 'background']),
              'background': lambda e: e == 'background', 'CGI': lambda e: e == 'CGI'}
    for var in ['soma_med', 'soma_min']:
        rows = []
        for g in GERM_LABELS:
            for gname, f in groups.items():
                cs, ss = cpg[(cpg.germ == g) & f(cpg.element)], sw[(sw.germ == g) & f(sw.element)]
                c_b, s_b = block_counts(cs, [var]), block_counts(ss, [var])
                for b in SOMA_LABELS:
                    c1, s1 = cs[cs[var] == b], ss[ss[var] == b]
                    if c1.D.sum() < 20 or s1.D.sum() < 20:
                        continue
                    def stat(c, s, b=b):
                        c = c.groupby(level=0).sum() if 'block' in c.index.names else c
                        s = s.groupby(level=0).sum() if 'block' in s.index.names else s
                        return pdr(c, b) / pdr(s, b)
                    est, lo, hi = boot(c_b, s_b, stat, a.n_boot, rng)
                    rows.append(dict(germ=g, group=gname, bin=b, cpg_P=c1.P.sum(), cpg_D=c1.D.sum(),
                                     cpg_PD=c1.P.sum() / c1.D.sum(), sw_PD=s1.P.sum() / s1.D.sum(),
                                     ratio=est, ratio_lo=lo, ratio_hi=hi))
                present = [b for b in SOMA_LABELS if ((cs[var] == b).sum() and cs[cs[var] == b].D.sum() >= 20
                                                     and ss[ss[var] == b].D.sum() >= 20)]
                if len(present) >= 2:
                    lo_b, hi_b = present[0], present[-1]
                    def stat(c, s, lo_b=lo_b, hi_b=hi_b):
                        c = c.groupby(level=0).sum() if 'block' in c.index.names else c
                        s = s.groupby(level=0).sum() if 'block' in s.index.names else s
                        return (pdr(c, hi_b) / pdr(s, hi_b)) / (pdr(c, lo_b) / pdr(s, lo_b))
                    est, lo, hi = boot(c_b, s_b, stat, a.n_boot, rng)
                    rows.append(dict(germ=g, group=gname, bin=f'trend {hi_b} / {lo_b}', ratio=est, ratio_lo=lo, ratio_hi=hi))
        t = pd.DataFrame(rows)
        t.to_csv(f'{a.outdir}/R2c_{var}.tsv', sep='\t', index=False, float_format='%.4g')
        print(f'\n== R2 by {var} (ratio = CpG P/D / non-CpG S>W P/D)\n' + t.to_string(index=False, float_format=lambda v: f'{v:.3g}'), flush=True)


if __name__ == '__main__':
    main()
