"""R2d: are the sperm-low / oocyte-high ("low|high") CpGs with excess CpG-loss constraint imprinted loci?

Inputs: R2b tables run with --oocyte (soma_oo), R1g motif tables (strong hits excluded, as in R2c), GENCODE genes,
and geneimprint's human table (tsv/R2/geneimprint_human_2026-10-03.tsv; 'Imprinted' status only, 131 genes).

Germline bins as in R2c: conserved sperm bin (chimp site m and human regional m agree) | MII oocyte regional bin.
Outputs (tsv/R2/):
  R2d_pd_imprint.tsv     P/D by germline bin x imprinted locus (gene body +-100 kb) x element group x soma_med bin
  R2d_lowhigh_regions.tsv  low|high element CpGs with an event, clustered (<= 2 kb gaps): counts, P/D, mean somatic
                         m, nearest gene, imprinted flag
  R2d_lowhigh_soma_profile.tsv  low|high element CpGs by somatic profile: imprint-like (constitutive
                         hemimethylation: soma_min_reg >= 0.25 and soma_med_reg 0.3-0.7) vs other
"""
from __future__ import annotations

import argparse
import glob
import gzip
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from R1e_pd_tables import events  # noqa: E402
from R2b_cpg_soma import GERM_BINS, GERM_LABELS, SOMA_BINS, SOMA_LABELS  # noqa: E402


def genes(gtf: str) -> pd.DataFrame:
    rows = []
    with gzip.open(gtf, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.split('\t')
            if f[2] != 'gene' or not f[0].startswith('chr') or '_' in f[0]:
                continue
            name = f[8].split('gene_name "')[1].split('"')[0]
            gtype = f[8].split('gene_type "')[1].split('"')[0]
            rows.append((f[0], int(f[3]) - 1, int(f[4]), name, gtype))
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'gene', 'gene_type'])


def nearest(g: pd.DataFrame, chrom: str, pos: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = g[g.chrom == chrom]
    if len(x) == 0:
        return np.full(len(pos), '.'), np.full(len(pos), np.inf)
    st, en = x.start.values, x.end.values
    d = np.where(pos[:, None] < st[None, :], st[None, :] - pos[:, None],
                 np.where(pos[:, None] > en[None, :], pos[:, None] - en[None, :], 0))
    k = d.argmin(1)
    return x.gene.values[k], d[np.arange(len(pos)), k]


def pdr(x: pd.DataFrame) -> float:
    return x.P.sum() / x.D.sum() if x.D.sum() > 0 else np.nan


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--soma-glob', required=True)
    ap.add_argument('--motif-dir', required=True)
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--imprint', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--flank', type=int, default=100_000)
    ap.add_argument('--strong', type=float, default=2.0)
    a = ap.parse_args()

    imp = pd.read_csv(a.imprint, sep='\t')
    imp_genes = set(imp.loc[imp.Status == 'Imprinted', 'Gene'])
    g = genes(a.gtf)
    gi = g[g.gene.isin(imp_genes)]
    print(f'imprinted genes: {len(imp_genes)} listed, {gi.gene.nunique()} found in GENCODE', flush=True)
    pc = g[g.gene_type == 'protein_coding']

    parts = []
    for f in sorted(glob.glob(a.soma_glob)):
        x = pd.read_csv(f, sep='\t', dtype={'snv_pos': str, 'snv_kind': str, 'cons_C': str, 'cons_G': str})
        chrom = f.split('/')[-1].split('.')[0]
        mo = pd.read_csv(f'{a.motif_dir}/{chrom}.motif.tsv.gz', sep='\t', usecols=['pos', 'best_margin'])
        x = x.merge(mo, on='pos', how='left')
        x = x[~(x.best_margin >= a.strong) & (x.state != 'other')]
        x = pd.concat([x, events(x)], axis=1)
        x = x[(x.P + x.D) > 0].copy()
        x['chrom'] = chrom
        xi = gi[gi.chrom == chrom]
        p = x.pos.values
        near = np.zeros(len(x), dtype=bool)
        for s, e in zip(xi.start.values, xi.end.values):
            near |= (p >= s - a.flank) & (p <= e + a.flank)
        x['imprinted_locus'] = near
        parts.append(x)
    d = pd.concat(parts, ignore_index=True)
    h_bin = pd.cut(d.h_reg, GERM_BINS, labels=GERM_LABELS).astype(str)
    c_bin = np.where(d.c_cov >= 5, pd.cut(d.c_m, GERM_BINS, labels=GERM_LABELS).astype(str), 'nan')
    sperm = np.where(h_bin == c_bin, h_bin, 'discordant')
    d['germ'] = pd.Series(sperm, index=d.index) + '|' + pd.cut(d.oo_reg, GERM_BINS, labels=GERM_LABELS).astype(str)
    d['element'] = np.where(d.cons != '.', 'coding', np.where(d.cgi == 1, 'CGI', np.where(d.ccre != 'none', d.ccre, 'background')))
    d['egroup'] = np.where(d.element == 'coding', 'coding', np.where(d.element == 'background', 'background', 'elements'))
    d['soma_med'] = pd.cut(d.soma_med_reg, SOMA_BINS, labels=SOMA_LABELS).astype(str)
    print(f'events: {len(d):,}; imprinted-locus events {d.imprinted_locus.sum():,}', flush=True)

    rows = []
    for key, x in d[d.germ.isin(['low|high', 'low|low', 'high|high', 'high|low'])].groupby(
            ['germ', 'imprinted_locus', 'egroup', 'soma_med']):
        if x.D.sum() >= 5:
            rows.append(dict(zip(['germ', 'imprinted_locus', 'egroup', 'soma_med'], key), n=len(x), P=x.P.sum(), D=x.D.sum(), PD=pdr(x)))
    t = pd.DataFrame(rows)
    t.to_csv(f'{a.outdir}/R2d_pd_imprint.tsv', sep='\t', index=False, float_format='%.4g')
    print('\n== P/D by germline x imprinted locus x element group x somatic bin\n' + t.to_string(index=False, float_format=lambda v: f'{v:.3g}'), flush=True)

    lh = d[(d.germ == 'low|high') & (d.egroup == 'elements')].sort_values(['chrom', 'pos']).copy()
    lh['imprint_like'] = (lh.soma_min_reg >= 0.25) & lh.soma_med_reg.between(0.3, 0.7)
    prof = lh.groupby(['imprint_like', 'imprinted_locus']).apply(lambda x: pd.Series({'n': len(x), 'P': x.P.sum(), 'D': x.D.sum(), 'PD': pdr(x)}))
    prof.to_csv(f'{a.outdir}/R2d_lowhigh_soma_profile.tsv', sep='\t', float_format='%.4g')
    print('\n== low|high element CpGs by somatic profile\n' + prof.to_string(float_format=lambda v: f'{v:.3g}'), flush=True)
    print('\nsomatic median m of low|high element CpGs: ' + lh.soma_med.value_counts().sort_index().to_string().replace('\n', '; '), flush=True)

    # clusters
    new = (lh.chrom != lh.chrom.shift()) | (lh.pos - lh.pos.shift() > 2000)
    lh['cluster'] = new.cumsum()
    cl = lh.groupby('cluster').agg(chrom=('chrom', 'first'), start=('pos', 'min'), end=('pos', 'max'), n=('pos', 'size'),
                                   P=('P', 'sum'), D=('D', 'sum'), soma_med=('soma_med_reg', 'mean'),
                                   soma_min=('soma_min_reg', 'mean'), oo=('oo_reg', 'mean'), cgi=('cgi', 'max'),
                                   imprinted_locus=('imprinted_locus', 'max'))
    gene, dist = [], []
    for chrom, x in cl.groupby('chrom'):
        gg, dd = nearest(pc, chrom, ((x.start + x.end) // 2).values)
        gene += list(zip(x.index, gg, dd))
    gd = pd.DataFrame(gene, columns=['cluster', 'nearest_gene', 'dist']).set_index('cluster')
    cl = cl.join(gd)
    cl['imprinted_gene_nearest'] = cl.nearest_gene.isin(imp_genes)
    cl = cl.sort_values('n', ascending=False)
    cl.to_csv(f'{a.outdir}/R2d_lowhigh_regions.tsv', sep='\t', float_format='%.3f')
    print(f'\n== low|high element clusters: {len(cl):,}; with an imprinted gene within {a.flank // 1000} kb: '
          f'{cl.imprinted_locus.sum()}; events in them {lh.imprinted_locus.sum()} of {len(lh)}', flush=True)
    print(cl.head(25).to_string(float_format=lambda v: f'{v:.2f}'), flush=True)


if __name__ == '__main__':
    main()
