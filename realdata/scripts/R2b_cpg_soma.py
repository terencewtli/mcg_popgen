"""R2b: regional germline and somatic methylation for ancestral CpGs (R1d) on one chromosome.

At a fixed human loss the hg38 site is no longer a CpG, so neither human sperm nor the Loyfer atlas measures it. Both
are therefore taken regionally, identically for P and D sites: the mean over the other CpGs within +-flank bp.
  soma_med_reg   mean of Loyfer per-CpG median-across-cell-types m (CpGs covered in >= min_groups groups)
  soma_min_reg   mean of per-CpG minimum across cell types (low = demethylated in some cell type, i.e. active there)
  soma_med_site, soma_min_site   the site's own values (retained CpGs only; NaN at losses)
regional() is shared with R1f, so the non-CpG comparator gets the same definitions.
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

SOMA_BINS = [-0.01, 0.2, 0.5, 0.8, 1.01]
SOMA_LABELS = ['0-0.2', '0.2-0.5', '0.5-0.8', '0.8-1']
GERM_BINS = [-0.01, 0.2, 0.8, 1.01]
GERM_LABELS = ['low', 'mid', 'high']


def regional(p: np.ndarray, v: np.ndarray, q: np.ndarray, flank: int = 500) -> np.ndarray:
    """Mean of v at sorted positions p within +-flank of each query q, excluding a value exactly at q."""
    ok = ~np.isnan(v)
    p, v = p[ok], v[ok]
    cs = np.concatenate([[0.0], np.cumsum(v)])
    lo = np.searchsorted(p, q - flank, side='left')
    hi = np.searchsorted(p, q + flank, side='right')
    tot, n = cs[hi] - cs[lo], (hi - lo).astype(float)
    j = np.searchsorted(p, q)
    jc = np.minimum(j, len(p) - 1)
    at = (j < len(p)) & (p[jc] == q)
    tot = tot - np.where(at, v[jc], 0)
    n = n - at
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(n > 0, tot / n, np.nan)


def loyfer_arrays(path: str, min_groups: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    s = pd.read_csv(path, sep='\t', usecols=['pos', 'm_med', 'm_min', 'n_groups'])
    s.loc[s.n_groups < min_groups, ['m_med', 'm_min']] = np.nan
    return s.pos.values, s.m_med.values.astype(float), s.m_min.values.astype(float)


def human_sperm_arrays(path: str, chrom: str, seq: str, min_cov: int = 5) -> tuple[np.ndarray, np.ndarray]:
    s = pd.read_csv(path, sep='\t', header=None, names=['chrom', 'start', 'end', 'm', 'cov'])
    s = s[(s.chrom == chrom) & (s['cov'] >= min_cov)]
    st = s.start.values
    on_g = np.array([seq[x] == 'G' and x > 0 and seq[x - 1] == 'C' for x in st], dtype=bool)
    s = pd.DataFrame({'p': np.where(on_g, st - 1, st), 'm': s.m.values}).groupby('p').m.mean()
    return s.index.values, s.values.astype(float)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--chrom', required=True)
    ap.add_argument('--annot', required=True)
    ap.add_argument('--loyfer', required=True, help='R2a chrN.tsv.gz')
    ap.add_argument('--min-groups', type=int, default=40)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    d = pd.read_csv(a.annot, sep='\t', dtype={'snv_pos': str, 'snv_kind': str, 'cons_C': str, 'cons_G': str})
    p0 = d.pos.values - 1
    lp, med, mn = loyfer_arrays(a.loyfer, a.min_groups)
    d['soma_med_reg'] = regional(lp, med, p0)
    d['soma_min_reg'] = regional(lp, mn, p0)
    site = pd.DataFrame({'med': med, 'mn': mn}, index=lp)
    d['soma_med_site'] = site.med.reindex(p0).values
    d['soma_min_site'] = site.mn.reindex(p0).values
    print(f'{a.chrom}: {len(d):,} CpGs; regional somatic m available for {np.mean(~np.isnan(d.soma_med_reg)):.3f}; '
          f'corr(site, regional) med {pd.Series(d.soma_med_site).corr(pd.Series(d.soma_med_reg)):.2f}', flush=True)
    d.to_csv(a.out, sep='\t', index=False, compression='gzip', float_format='%.4g')
    print('done', flush=True)


if __name__ == '__main__':
    main()
