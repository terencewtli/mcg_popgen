"""R1g: TF-motif overlap of ancestral CpGs (for the F6 motif exclusion in R2), one chromosome.

Two modes:
  fasta     write the chromosome with every ancestral CpG restored (C at pos, G at pos + 1; R1b), plus a bed of
            the CpGs (2 bp). Scanning the ancestral sequence means a CpG lost on the human lineage is not
            also lost from its motif, which would bias P/D.
  collapse  read `bedtools intersect -wa -wb` lines (HOMER hit bed + CpG bed) on stdin and write one row per CpG:
            n_hits (full motifs; half-sites are dropped upstream), best_motif, best_margin (log-odds score minus
            HOMER's threshold for that motif name; slightly negative where a name occurs twice in the library with
            different thresholds), best_motif_has_cg (the motif consensus
            contains CG).
HOMER known vertebrate motifs (472, with thresholds) are used; JASPAR is not needed for this step.
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd
import pysam


def thresholds(motif_file: str) -> dict[str, tuple[float, bool]]:
    out = {}
    with open(motif_file) as fh:
        for line in fh:
            if line.startswith('>'):
                f = line[1:].rstrip('\n').split('\t')
                name = f[1].split('/')[0]
                out[name] = (float(f[2]), 'CG' in f[0].upper())
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['fasta', 'collapse'])
    ap.add_argument('--chrom', required=True)
    ap.add_argument('--anc', required=True)
    ap.add_argument('--fasta')
    ap.add_argument('--out-fa')
    ap.add_argument('--out-bed')
    ap.add_argument('--motifs')
    ap.add_argument('--out')
    a = ap.parse_args()

    d = pd.read_csv(a.anc, sep='\t', usecols=['pos'])
    p0 = d.pos.values - 1
    if a.mode == 'fasta':
        s = np.frombuffer(pysam.FastaFile(a.fasta).fetch(a.chrom).upper().encode(), dtype=np.uint8).copy()
        s[p0], s[p0 + 1] = ord('C'), ord('G')
        seq = s.tobytes().decode()
        with open(a.out_fa, 'w') as fh:
            fh.write(f'>{a.chrom}\n')
            for i in range(0, len(seq), 100):
                fh.write(seq[i:i + 100] + '\n')
        pd.DataFrame({'c': a.chrom, 's': p0, 'e': p0 + 2}).to_csv(a.out_bed, sep='\t', header=False, index=False)
        print(f'{a.chrom}: {len(p0):,} ancestral CpGs restored and written', flush=True)
        return

    thr = thresholds(a.motifs)
    best: dict[int, list] = {}
    n = 0
    for line in sys.stdin:
        f = line.rstrip('\n').split('\t')
        # hit: chrom start end name score strand | cpg: chrom start end
        name, score, cpg = f[3].split('/')[0], float(f[4]), int(f[7])
        t, has_cg = thr.get(name, (np.nan, False))
        margin = score - t
        b = best.get(cpg)
        if b is None:
            best[cpg] = [1, name, margin, has_cg]
        else:
            b[0] += 1
            if margin > b[2]:
                b[1], b[2], b[3] = name, margin, has_cg
        n += 1
    out = pd.DataFrame({'pos': p0 + 1})
    hit = pd.DataFrame.from_dict(best, orient='index', columns=['n_hits', 'best_motif', 'best_margin', 'best_motif_has_cg'])
    hit.index = hit.index + 1
    out = out.merge(hit, left_on='pos', right_index=True, how='left')
    out['n_hits'] = out.n_hits.fillna(0).astype(int)
    out['best_motif'] = out.best_motif.fillna('.')
    out['best_motif_has_cg'] = out.best_motif_has_cg.fillna(False).astype(int)
    out.to_csv(a.out, sep='\t', index=False, compression='gzip', float_format='%.3f')
    print(f'{a.chrom}: {n:,} hit-CpG overlaps; {np.mean(out.n_hits > 0):.3f} of CpGs in >= 1 motif; '
          f'margin >= 2: {np.mean(out.best_margin >= 2):.3f}', flush=True)


if __name__ == '__main__':
    main()
