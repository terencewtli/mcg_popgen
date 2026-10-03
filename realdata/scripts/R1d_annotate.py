"""R1d: annotate R1b ancestral CpGs on one chromosome.

Added columns:
  h_m h_cov   human sperm methylation / coverage (Molaro 2011, lifted to hg38; R1c)
  c_m c_cov   chimp sperm methylation / coverage at the orthologous position (panTro2 -> hg19 -> hg38)
  h_reg h_reg_n, c_reg c_reg_n   mean sperm methylation of the other CpGs (cov >= 5) within +-500 bp, and their
              number. At fixed human losses the site has no human CpG, so h_m is missing; h_reg is the human
              germline state there.
  ccre        ENCODE V3 cCRE class (PLS > pELS > dELS > CTCF-only > DNase-H3K4me3 > none)
  ctcf        1 if the cCRE is CTCF-bound
  cgi         1 inside a UCSC CpG island
  phylop      mean phyloP100way over the two bases
  cons_C      coding consequence of C>T at the C on the ancestral codon (syn / mis / stop / "." if not coding),
              Ensembl canonical CDS; fixed losses are scored on the ancestral CpG
  cons_G      same for G>A at the G
  tx gene     canonical transcript (no version) and gene name, '.' if not coding

Sperm values are matched to the C of the dinucleotide; a lifted value that lands on the G (minus-strand chain) is
moved to the C.
"""
from __future__ import annotations

import argparse
import gzip

import numpy as np
import pandas as pd
import pyBigWig
import pysam

CCRE_ORDER = ['PLS', 'pELS', 'dELS', 'CTCF-only', 'DNase-H3K4me3']
SEVERITY = {'.': 0, 'syn': 1, 'mis': 2, 'stop': 3}
COMP = str.maketrans('ACGT', 'TGCA')
BASES = 'TCAG'
AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES)
         for k, c in enumerate(BASES)}


def sperm_at(path: str, chrom: str, pos0: np.ndarray, seq: str) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    s = pd.read_csv(path, sep='\t', header=None, names=['chrom', 'start', 'end', 'm', 'cov'])
    s = s[s.chrom == chrom]
    st = s.start.values
    on_g = np.array([seq[p] == 'G' and p > 0 and seq[p - 1] == 'C' for p in st], dtype=bool)
    st = np.where(on_g, st - 1, st)
    s = pd.DataFrame({'p': st, 'm': s.m.values, 'cov': s['cov'].values}).groupby('p').agg({'m': 'mean', 'cov': 'sum'})
    m = s.m.reindex(pos0).values
    cov = s['cov'].reindex(pos0).fillna(0).astype(int).values
    return m, cov, s


def regional(s: pd.DataFrame, pos0: np.ndarray, flank: int = 500, min_cov: int = 5) -> tuple[np.ndarray, np.ndarray]:
    """Mean sperm methylation of the other CpGs (cov >= min_cov) within +-flank bp, and how many there are."""
    s = s[s['cov'] >= min_cov]
    p = s.index.values
    cm = np.concatenate([[0.0], np.cumsum(s.m.values)])
    lo = np.searchsorted(p, pos0 - flank, side='left')
    hi = np.searchsorted(p, pos0 + flank, side='right')
    tot, n = cm[hi] - cm[lo], hi - lo
    self_i = np.searchsorted(p, pos0)
    is_self = (self_i < len(p)) & (p[np.minimum(self_i, len(p) - 1)] == pos0)
    tot = tot - np.where(is_self, s.m.values[np.minimum(self_i, len(p) - 1)], 0)
    n = n - is_self
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(n > 0, tot / n, np.nan), n


def interval_mask(starts: np.ndarray, ends: np.ndarray, L: int) -> np.ndarray:
    d = np.zeros(L + 1, dtype=np.int32)
    np.add.at(d, starts, 1)
    np.add.at(d, ends, -1)
    return np.cumsum(d[:-1]) > 0


def coding(gtf: str, chrom: str, seq: str, want: set[int]) -> dict[int, tuple[str, str, str]]:
    """0-based genomic position -> (consequence of the CpG transition at that base, tx, gene)."""
    tx: dict[str, list] = {}
    with gzip.open(gtf, 'rt') as fh:
        for line in fh:
            if not line.startswith(chrom + '\t'):
                continue
            f = line.rstrip('\n').split('\t')
            if f[2] != 'CDS' or 'Ensembl_canonical' not in f[8] or 'gene_type "protein_coding"' not in f[8]:
                continue
            if 'start_NF' in f[8] or 'end_NF' in f[8]:
                continue
            tid = f[8].split('transcript_id "')[1].split('"')[0].split('.')[0]
            gene = f[8].split('gene_name "')[1].split('"')[0]
            tx.setdefault(tid, [gene, f[6], []])[2].append((int(f[3]) - 1, int(f[4])))
    out: dict[int, tuple[str, str, str]] = {}
    for tid, (gene, strand, segs) in tx.items():
        segs.sort()
        g = np.concatenate([np.arange(s, e) for s, e in segs])
        if strand == '-':
            g = g[::-1]
        cds = ''.join(seq[p] for p in g)
        if strand == '-':
            cds = cds.translate(COMP)
        n = len(cds) - len(cds) % 3
        for k in np.where(np.isin(g[:n], list(want)))[0]:
            p = int(g[k])
            b = seq[p]
            mut = {'C': 'T', 'G': 'A'}.get(b)
            if mut is None:
                continue
            if strand == '-':
                mut = mut.translate(COMP)
            c0 = k - k % 3
            ref = cds[c0:c0 + 3]
            alt = ref[:k % 3] + mut + ref[k % 3 + 1:]
            ra, aa = CODON.get(ref), CODON.get(alt)
            if ra is None or aa is None:
                continue
            cons = 'syn' if ra == aa else ('stop' if aa == '*' else 'mis')
            if SEVERITY[cons] > SEVERITY[out.get(p, ('.',))[0]]:
                out[p] = (cons, tid, gene)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--chrom', required=True)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--anc', required=True)
    ap.add_argument('--human-sperm', required=True)
    ap.add_argument('--chimp-sperm', required=True)
    ap.add_argument('--ccre', required=True)
    ap.add_argument('--cgi', required=True)
    ap.add_argument('--phylop', required=True)
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    seq = pysam.FastaFile(a.fasta).fetch(a.chrom).upper()
    L = len(seq)
    d = pd.read_csv(a.anc, sep='\t')
    p0 = d.pos.values - 1
    print(f'{a.chrom}: {len(d):,} ancestral CpGs', flush=True)

    d['h_m'], d['h_cov'], hs = sperm_at(a.human_sperm, a.chrom, p0, seq)
    d['c_m'], d['c_cov'], cs = sperm_at(a.chimp_sperm, a.chrom, p0, seq)
    d['h_reg'], d['h_reg_n'] = regional(hs, p0)
    d['c_reg'], d['c_reg_n'] = regional(cs, p0)
    print(f'  sperm: human {np.mean(d.h_cov > 0):.2f}, chimp {np.mean(d.c_cov > 0):.2f} of sites covered', flush=True)

    cc = pd.read_csv(a.ccre, sep='\t', header=None, usecols=[0, 1, 2, 5], names=['chrom', 'start', 'end', 'cls'])
    cc = cc[cc.chrom == a.chrom]
    code = np.zeros(L, dtype=np.int8)
    for r, c in enumerate(CCRE_ORDER[::-1], start=1):  # higher code = higher priority, written last
        sub = cc[cc.cls.str.split(',').str[0] == c]
        code[interval_mask(sub.start.values, sub.end.values, L)] = r
    names = np.array(['none'] + CCRE_ORDER[::-1])
    d['ccre'] = names[code[p0]]
    ctcf = cc[cc.cls.str.contains('CTCF')]
    d['ctcf'] = interval_mask(ctcf.start.values, ctcf.end.values, L)[p0].astype(int)

    cgi = pd.read_csv(a.cgi, sep='\t', header=None, usecols=[1, 2, 3], names=['chrom', 'start', 'end'])
    cgi = cgi[cgi.chrom == a.chrom]
    d['cgi'] = interval_mask(cgi.start.values, cgi.end.values, L)[p0].astype(int)

    bw = pyBigWig.open(a.phylop)
    v = np.array(bw.values(a.chrom, 0, L), dtype=np.float32)
    d['phylop'] = np.nanmean(np.vstack([v[p0], v[p0 + 1]]), axis=0)
    del v

    # consequences are for the ancestral -> derived change, so restore the ancestral CpG at fixed losses first
    anc_arr = np.frombuffer(seq.encode(), dtype=np.uint8).copy()
    anc_arr[p0], anc_arr[p0 + 1] = ord('C'), ord('G')
    cod = coding(a.gtf, a.chrom, anc_arr.tobytes().decode(), set(p0.tolist()) | set((p0 + 1).tolist()))
    d['cons_C'] = [cod.get(p, ('.',))[0] for p in p0]
    d['cons_G'] = [cod.get(p + 1, ('.',))[0] for p in p0]
    d['tx'] = [cod.get(p, cod.get(p + 1, ('.', '.', '.')))[1] for p in p0]
    d['gene'] = [cod.get(p, cod.get(p + 1, ('.', '.', '.')))[2] for p in p0]
    print('  ccre: ' + ', '.join(f'{k} {v:,}' for k, v in d.ccre.value_counts().items()), flush=True)
    print(f'  cgi {d.cgi.sum():,}; coding C: ' + ', '.join(f'{k} {v:,}' for k, v in d.cons_C.value_counts().items()),
          flush=True)
    d.to_csv(a.out, sep='\t', index=False, compression='gzip', float_format='%.4g')
    print('done', flush=True)


if __name__ == '__main__':
    main()
