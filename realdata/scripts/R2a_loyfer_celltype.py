"""R2a: Loyfer 2023 atlas (253 hg38 .beta files) -> per-CpG methylation by cell-type group, one chromosome.

.beta format (wgbstools): uint8 pairs (methylated, covered) per CpG, in wgbstools hg38 CpG order. That order is all
'CG' dinucleotides on chr1-22, X, Y, M of UCSC hg38 (29,401,795; checked against the file size). Chromosome order is
not documented in the files, so each task tests the candidate orders and keeps the one where CpG-island CpGs are
least methylated relative to the rest (a misaligned offset smears the CGI dip away). The contrast for each order is
logged.

Samples: the sorted-cell atlas only (titles ending in a '-Z...' donor ID); the cfDNA_* and CNVS-* samples in the same
GEO series are plasma / CNV-study samples, not cell types. Groups: the title without the donor suffix
(e.g. 'Saphenous-Vein-Endothelium-Z000000RM' -> 'Saphenous-Vein-Endothelium'), i.e. tissue-specific cell types. Counts are pooled within a group; m = meth / cov where cov >= min_cov.

Outputs:
  <out>.npz      pos (0-based C), groups, m (n x G float16, NaN below min_cov), cov (n x G uint16)
  <out>.tsv.gz   pos, m_max, m_med, m_min, n_groups (groups with cov >= min_cov), g_max, g_min
"""
from __future__ import annotations

import argparse
import gzip
import re

import numpy as np
import pandas as pd
import pysam

ORDERS = {
    'natural_XYM': [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM'],
    'M_first': ['chrM'] + [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY'],
    'natural_MXY': [f'chr{i}' for i in range(1, 23)] + ['chrM', 'chrX', 'chrY'],
}


def sample_groups(matrix: str) -> pd.DataFrame:
    rows = {}
    with gzip.open(matrix, 'rt') as fh:
        for line in fh:
            if line.startswith('!Sample_title') or line.startswith('!Sample_geo_accession'):
                f = [x.strip('"') for x in line.rstrip('\n').split('\t')]
                rows[f[0]] = f[1:]
    s = pd.DataFrame({'gsm': rows['!Sample_geo_accession'], 'title': rows['!Sample_title']})
    s = s[s.title.str.contains(r'-Z[0-9A-Z]+$')]  # sorted-cell atlas only; drops cfDNA_* and CNVS-* (plasma / CNV study)
    s['group'] = [re.sub(r'-Z[0-9A-Z]+$', '', t) for t in s.title]
    return s


def cgi_contrast(beta: np.ndarray, is_cgi: np.ndarray) -> float:
    cov = beta[:, 1].astype(float)
    m = np.where(cov >= 5, beta[:, 0] / np.maximum(cov, 1), np.nan)
    return float(np.nanmean(m[~is_cgi]) - np.nanmean(m[is_cgi]))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--chrom', required=True)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--beta-dir', required=True)
    ap.add_argument('--urls', required=True, help='tsv/R0/loyfer_hg38_beta_urls.tsv (gsm, file, size, url)')
    ap.add_argument('--matrix', required=True, help='GSE186458 series matrix')
    ap.add_argument('--cgi', required=True)
    ap.add_argument('--min-cov', type=int, default=5)
    ap.add_argument('--out', required=True, help='prefix')
    a = ap.parse_args()

    fa = pysam.FastaFile(a.fasta)
    seqs = {c: fa.fetch(c).upper() for c in ORDERS['natural_XYM']}
    pos = {c: np.array([m.start() for m in re.finditer('CG', s)], dtype=np.int64) for c, s in seqs.items()}
    n_tot = sum(len(p) for p in pos.values())
    print(f'{a.chrom}: {len(pos[a.chrom]):,} CpGs; genome {n_tot:,}', flush=True)
    p = pos[a.chrom]

    files = pd.read_csv(a.urls, sep='\t', header=None, names=['gsm', 'file', 'size', 'url'])
    meta = sample_groups(a.matrix).merge(files[['gsm', 'file']], on='gsm')
    print(f'{len(meta)} samples in {meta.group.nunique()} groups', flush=True)

    cgi = pd.read_csv(a.cgi, sep='\t', header=None, usecols=[1, 2, 3], names=['chrom', 'start', 'end'])
    cgi = cgi[cgi.chrom == a.chrom]
    d = np.zeros(len(seqs[a.chrom]) + 1, dtype=np.int32)
    np.add.at(d, cgi.start.values, 1)
    np.add.at(d, cgi.end.values, -1)
    is_cgi = (np.cumsum(d[:-1]) > 0)[p]

    # pick the chromosome order with the strongest CGI contrast, using a blood sample (deep, CGIs unmethylated)
    test = meta[meta.group.str.contains('Blood')].file.iloc[0] if meta.group.str.contains('Blood').any() else meta.file.iloc[0]
    best, scores = None, {}
    for name, order in ORDERS.items():
        off = sum(len(pos[c]) for c in order[:order.index(a.chrom)])
        mm = np.memmap(f'{a.beta_dir}/{test}', dtype=np.uint8, mode='r', offset=2 * off, shape=(len(p), 2))
        scores[name] = cgi_contrast(np.asarray(mm), is_cgi)
    best = max(scores, key=scores.get)
    print('CGI contrast (non-CGI m - CGI m) by order: ' + ', '.join(f'{k} {v:.3f}' for k, v in scores.items())
          + f' -> {best}', flush=True)
    order = ORDERS[best]
    off = sum(len(pos[c]) for c in order[:order.index(a.chrom)])

    groups = sorted(meta.group.unique())
    meth = np.zeros((len(p), len(groups)), dtype=np.uint32)
    cov = np.zeros((len(p), len(groups)), dtype=np.uint32)
    for j, g in enumerate(groups):
        for f in meta.file[meta.group == g]:
            mm = np.memmap(f'{a.beta_dir}/{f}', dtype=np.uint8, mode='r', offset=2 * off, shape=(len(p), 2))
            meth[:, j] += mm[:, 0]
            cov[:, j] += mm[:, 1]
    with np.errstate(invalid='ignore', divide='ignore'):
        m = np.where(cov >= a.min_cov, meth / cov, np.nan).astype(np.float32)
    ok = ~np.isnan(m)
    n_ok = ok.sum(1)
    print(f'groups covered per CpG: median {np.median(n_ok):.0f} of {len(groups)}', flush=True)

    np.savez_compressed(f'{a.out}.npz', pos=p, groups=np.array(groups), m=m.astype(np.float16),
                        cov=np.minimum(cov, 65535).astype(np.uint16))
    with np.errstate(all='ignore'):
        mx, md, mn = np.nanmax(m, 1), np.nanmedian(m, 1), np.nanmin(m, 1)
    m_fill_hi = np.where(ok, m, -1)
    m_fill_lo = np.where(ok, m, 2)
    gnames = np.array(groups)
    summ = pd.DataFrame({'pos': p, 'm_max': mx, 'm_med': md, 'm_min': mn, 'n_groups': n_ok,
                         'g_max': np.where(n_ok > 0, gnames[m_fill_hi.argmax(1)], '.'),
                         'g_min': np.where(n_ok > 0, gnames[m_fill_lo.argmin(1)], '.')})
    summ.to_csv(f'{a.out}.tsv.gz', sep='\t', index=False, float_format='%.3f', compression='gzip')
    print(f'CGI median of m_med {np.nanmedian(md[is_cgi]):.3f}, non-CGI {np.nanmedian(md[~is_cgi]):.3f}', flush=True)
    print('done', flush=True)


if __name__ == '__main__':
    main()
