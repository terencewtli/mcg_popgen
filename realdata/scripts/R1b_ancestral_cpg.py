"""R1b: ancestral CpGs on one hg38 chromosome, their fate on the human lineage, and their 1000G polymorphism.

Ancestral CpG: the hg38 dinucleotide (i, i+1) aligns to CG in both chimp (panTro6) and macaque (rheMac10) in the
UCSC reciprocal-best axt files, as two adjacent alignment columns with no insertion between them.

Human state in hg38:
  CG               retained
  TG / CA          lost by a transition at the C / G (human-lineage substitution, D)
  anything else    lost otherwise (transversion or double change; kept, labelled 'other')

Polymorphism, joined from R1a (2,504 unrelated 1000G samples):
  retained site, SNV C>T at i or G>A at i+1        polymorphic loss; derived count = AC
  lost site (TG / CA), SNV restoring the C / G      ancestral allele still segregating; derived count = AN - AC.
                                                    These are P, not D; D counts only fixed losses.
Transversions at CpGs are recorded in snv_kind but are not the CpG-loss class used for P/D.

Output (one row per ancestral CpG): pos (1-based, the C), hg38 dinucleotide, state, snv_pos ('C', 'G' or '.'),
snv_kind (ts / tv / '.'), dac (derived allele count), an.
"""
from __future__ import annotations

import argparse
import gzip

import numpy as np
import pandas as pd
import pysam

GAP = ord('-')


def load_axt(path: str, chrom: str, length: int) -> tuple[np.ndarray, np.ndarray]:
    """Aligned query base per hg38 position (0 = unaligned) and a break flag after each position."""
    base = np.zeros(length, dtype=np.uint8)
    brk = np.zeros(length, dtype=bool)
    n_blocks = 0
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            f = line.split()
            tseq, qseq = next(fh).rstrip('\n'), next(fh).rstrip('\n')
            if f[1] != chrom:
                continue
            n_blocks += 1
            t = np.frombuffer(tseq.upper().encode(), dtype=np.uint8)
            q = np.frombuffer(qseq.upper().encode(), dtype=np.uint8)
            is_t = t != GAP
            hpos = int(f[2]) - 1 + np.cumsum(is_t) - 1  # 0-based hg38 position of each column (or of the previous base)
            cols = np.where(is_t)[0]
            qb = q[cols].copy()
            qb[qb == GAP] = 0
            base[hpos[cols]] = qb
            ins = np.where(~is_t)[0]  # query insertion: no adjacency across it
            ins = ins[hpos[ins] >= 0]
            brk[hpos[ins]] = True
            brk[hpos[cols[-1]]] = True  # block end
    print(f'  {path.split("/")[-1]}: {n_blocks:,} blocks on {chrom}, {np.count_nonzero(base):,} aligned bases', flush=True)
    return base, brk


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--chrom', required=True)
    ap.add_argument('--chimp-axt', required=True)
    ap.add_argument('--macaque-axt', required=True)
    ap.add_argument('--snv', required=True, help='R1a table for this chromosome')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    seq = np.frombuffer(pysam.FastaFile(a.fasta).fetch(a.chrom).upper().encode(), dtype=np.uint8)
    L = len(seq)
    print(f'{a.chrom}: {L:,} bp', flush=True)
    pt, pt_brk = load_axt(a.chimp_axt, a.chrom, L)
    rm, rm_brk = load_axt(a.macaque_axt, a.chrom, L)

    C, G, T, A = ord('C'), ord('G'), ord('T'), ord('A')
    i = np.arange(L - 1)
    anc = (pt[:-1] == C) & (pt[1:] == G) & (rm[:-1] == C) & (rm[1:] == G) & ~pt_brk[:-1] & ~rm_brk[:-1]
    i = i[anc]
    h0, h1 = seq[i], seq[i + 1]
    state = np.full(len(i), 'other', dtype=object)
    state[(h0 == C) & (h1 == G)] = 'retained'
    state[(h0 == T) & (h1 == G)] = 'lost_C'
    state[(h0 == C) & (h1 == A)] = 'lost_G'
    state[(h0 == ord('N')) | (h1 == ord('N'))] = 'N'
    keep = state != 'N'
    i, state, h0, h1 = i[keep], state[keep], h0[keep], h1[keep]
    out = pd.DataFrame({'pos': i + 1, 'dinuc': [chr(x) + chr(y) for x, y in zip(h0, h1)], 'state': state})
    print(f'  ancestral CpGs: {len(out):,}; ' + ', '.join(f'{k} {v:,}' for k, v in out.state.value_counts().items()),
          flush=True)

    # polymorphism: SNVs at the C (pos) or the G (pos + 1)
    snv = pd.read_csv(a.snv, sep='\t', usecols=['pos', 'ref', 'alt', 'ac', 'an'])
    snv = snv.drop_duplicates('pos', keep=False)  # biallelic table; drop the rare positions listed twice
    snv = snv.set_index('pos')
    out['snv_pos'], out['snv_kind'], out['dac'], out['an'] = '.', '.', 0, 0
    for which, off, anc_base in (('C', 0, 'C'), ('G', 1, 'G')):
        p = out.pos + off
        hit = p.isin(snv.index)
        s = snv.loc[p[hit]]
        s.index = out.index[hit]
        ref_anc = s.ref == anc_base  # hg38 carries the ancestral base, ALT is derived
        ts_from_anc = ref_anc & (s.alt == {'C': 'T', 'G': 'A'}[which])
        ts_to_anc = (s.ref == {'C': 'T', 'G': 'A'}[which]) & (s.alt == anc_base)
        dac = np.where(ref_anc, s.ac, s.an - s.ac)
        kind = np.where(ts_from_anc | ts_to_anc, 'ts', 'tv')
        # when both positions carry an SNV, keep the one with the larger derived count
        upd = s.index[(dac > out.loc[s.index, 'dac'].values)]
        sel = s.index.get_indexer(upd)
        out.loc[upd, 'snv_pos'] = which
        out.loc[upd, 'snv_kind'] = kind[sel]
        out.loc[upd, 'dac'] = dac[sel]
        out.loc[upd, 'an'] = s.an.values[sel]
    print('  SNVs at ancestral CpGs: ' + ', '.join(
        f'{k} {v:,}' for k, v in out[out.snv_pos != '.'].groupby(['state', 'snv_kind']).size().items()), flush=True)
    out.to_csv(a.out, sep='\t', index=False, compression='gzip')
    print('done', flush=True)


if __name__ == '__main__':
    main()
