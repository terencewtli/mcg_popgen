"""R1e: CpG-loss P/D tables (DESIGN_PLAN R1) from the R1d annotated ancestral CpGs (all chromosomes).

Events, per ancestral CpG (CG in chimp and macaque):
  P  a segregating CpG-loss transition: retained site with C>T / G>A (0 < derived count < AN), or a human-lineage
     loss whose ancestral allele still segregates
  D  a fixed human-lineage loss by transition (TG / CA in hg38, no restoring SNV)
  Ps P restricted to derived count >= 2 (non-singletons)

Germline stratum (per site): chimp sperm methylation at the site (c_m, cov >= 5), which exists for P and D sites
alike (every ancestral CpG is still CG in chimp). 'Conserved' additionally requires the human regional sperm
methylation (h_reg, CpGs within +-500 bp) in the same bin; a human-specific change in germline methylation would
otherwise make P/D != 1 without selection (REVIEW 3, F3).

Element strata: cCRE class, CGI, coding consequence at the event position (positive control: missense / stop vs
synonymous), and background (no cCRE, no CGI, not coding). phyloP is not used to define background (it includes
human, so it is circular with D).

Report: P/D with 95% CI from a 1 Mb block bootstrap, and P/D relative to background in the same germline bin.
P/D is computed on counts, so opportunity (number of sites) cancels.
"""
from __future__ import annotations

import argparse
import glob

import numpy as np
import pandas as pd

M_BINS = [-0.01, 0.2, 0.8, 1.01]
M_LABELS = ['low', 'mid', 'high']


def events(d: pd.DataFrame) -> pd.DataFrame:
    ts = d.snv_kind == 'ts'
    seg = (d.dac > 0) & (d.dac < d.an)
    lost_restoring = ((d.state == 'lost_C') & (d.snv_pos == 'C') | (d.state == 'lost_G') & (d.snv_pos == 'G')) & ts
    p = ((d.state == 'retained') & ts & seg) | lost_restoring
    dd = d.state.isin(['lost_C', 'lost_G']) & ~lost_restoring
    ev_pos = np.where(d.state == 'lost_C', 'C', np.where(d.state == 'lost_G', 'G', d.snv_pos))
    cons = np.where(ev_pos == 'C', d.cons_C, np.where(ev_pos == 'G', d.cons_G, '.'))
    return pd.DataFrame({'P': p.astype(int), 'Ps': (p & (d.dac >= 2)).astype(int), 'D': dd.astype(int),
                         'cons': cons}, index=d.index)


def boot_ci(g: pd.DataFrame, n_boot: int, rng: np.random.Generator) -> tuple[float, float]:
    blk = g.groupby('block')[['P', 'D']].sum()
    if len(blk) < 10 or blk.D.sum() == 0:
        return np.nan, np.nan
    idx = rng.integers(0, len(blk), size=(n_boot, len(blk)))
    P, D = blk.P.values[idx].sum(1), blk.D.values[idx].sum(1)
    r = P / np.where(D > 0, D, np.nan)
    return tuple(np.nanpercentile(r, [2.5, 97.5]))


def table(d: pd.DataFrame, by: list[str], n_boot: int, rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    for key, g in d.groupby(by):
        key = key if isinstance(key, tuple) else (key,)
        P, Ps, D = g.P.sum(), g.Ps.sum(), g.D.sum()
        lo, hi = boot_ci(g, n_boot, rng)
        rows.append(dict(zip(by, key), n_sites=len(g), P=P, Ps=Ps, D=D,
                         PD=P / D if D else np.nan, PsD=Ps / D if D else np.nan, PD_lo=lo, PD_hi=hi))
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--annot-glob', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--n-boot', type=int, default=200)
    ap.add_argument('--min-cov', type=int, default=5)
    ap.add_argument('--motif-dir', help='R1g per-chromosome motif tables; needed for --motif-mode other than all')
    ap.add_argument('--motif-mode', default='all', choices=['all', 'exclude_any', 'exclude_strong', 'only_strong'],
                    help='exclude_any: drop CpGs in any HOMER hit; exclude_strong: drop hits with margin >= '
                         '--strong; only_strong: keep only those (the F6 "TF-selected" reference)')
    ap.add_argument('--strong', type=float, default=2.0)
    a = ap.parse_args()
    rng = np.random.default_rng(1)

    parts = []
    for f in sorted(glob.glob(a.annot_glob)):
        x = pd.read_csv(f, sep='\t', dtype={'snv_pos': str, 'snv_kind': str, 'cons_C': str, 'cons_G': str})
        x['chrom'] = f.split('/')[-1].split('.')[0]
        if a.motif_mode != 'all':
            mo = pd.read_csv(f'{a.motif_dir}/{x.chrom.iloc[0]}.motif.tsv.gz', sep='\t')
            x = x.merge(mo[['pos', 'n_hits', 'best_margin']], on='pos', how='left')
            strong = x.best_margin >= a.strong
            keep = {'exclude_any': x.n_hits == 0, 'exclude_strong': ~strong, 'only_strong': strong}[a.motif_mode]
            x = x[keep.values]
        parts.append(x)
        print(f'{f}: {len(x):,}', flush=True)
    d = pd.concat(parts, ignore_index=True)
    d = d[d.state != 'other']
    d = pd.concat([d, events(d)], axis=1)
    d = d[(d.P + d.D) > 0]  # only sites with an event enter P/D
    d['block'] = d.chrom + ':' + (d.pos // 1_000_000).astype(str)
    print(f'events: P {d.P.sum():,} (non-singleton {d.Ps.sum():,}), D {d.D.sum():,}', flush=True)

    d['germ'] = np.where(d.c_cov >= a.min_cov, pd.cut(d.c_m, M_BINS, labels=M_LABELS).astype(str), 'nocov')
    h_bin = pd.cut(d.h_reg, M_BINS, labels=M_LABELS).astype(str)
    d['germ_cons'] = np.where((d.germ != 'nocov') & (h_bin == d.germ), d.germ, 'discordant_or_nocov')
    d['coding'] = np.where(d.cons == '.', 'noncoding', d.cons)
    d['element'] = np.where(d.coding != 'noncoding', 'coding',
                   np.where(d.cgi == 1, 'CGI',
                   np.where(d.ccre != 'none', d.ccre, 'background')))
    # phyloP is NOT used: phyloP100way includes human, so sites with a human-lineage substitution score lower, and a
    # phyloP < 0 'background' is enriched for D (P/D biased down; the phyloP >= 0 remainder biased up)

    out = {
        'overall': table(d.assign(all='all'), ['all'], a.n_boot, rng),
        'germ': table(d, ['germ'], a.n_boot, rng),
        'germ_cons': table(d, ['germ_cons'], a.n_boot, rng),
        'element': table(d, ['element'], a.n_boot, rng),
        'germ_cons_x_element': table(d[d.germ_cons != 'discordant_or_nocov'], ['germ_cons', 'element'], a.n_boot, rng),
        'coding_control': table(d[d.element == 'coding'], ['coding'], a.n_boot, rng),
        'coding_control_x_germ': table(d[(d.element == 'coding') & (d.germ != 'nocov')], ['germ', 'coding'],
                                       a.n_boot, rng),
    }
    # relative to background in the same germline bin
    t = out['germ_cons_x_element']
    bg = t[t.element == 'background'].set_index('germ_cons').PD
    t['PD_rel_bg'] = t.PD / t.germ_cons.map(bg)
    for k, t in out.items():
        tag = '' if a.motif_mode == 'all' else f'.{a.motif_mode}'
        t.to_csv(f'{a.outdir}/R1e_pd_{k}{tag}.tsv', sep='\t', index=False, float_format='%.4g')
        print(f'\n== {k}\n{t.to_string(index=False, float_format=lambda v: f"{v:.3g}")}', flush=True)


if __name__ == '__main__':
    main()
