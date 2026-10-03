"""R1a: annotate every biallelic 1000G SNV on one chromosome with its hg38 context and CpG / W-S class.

Counts are from the 2,504 unrelated samples (sum of the per-superpopulation *_unrel INFO fields).
Input is the INFO-only stream from bcftools (see scripts/qsub/R1a_snv_context.sh). Polarisation (ancestral CpG,
human-lineage losses) is done later in R1b with the chimp / macaque alignments; here everything is relative to hg38.

Output columns:
  pos      1-based hg38 position
  ref alt  hg38 REF / ALT
  ac an    allele count / number in the unrelated samples
  ctx      hg38 trinucleotide centred on the SNV
  cpg      'C' if REF is the C of a CpG, 'G' if it is the G of a CpG, '.' otherwise
  klass    cpg_loss_ts (CpG C>T or G>A), cpg_loss_tv (other change at a CpG), cpg_gain (ALT creates a CpG),
           noncpg (none of these)
  ws       W>S, S>W, W>W or S>S (for the gBGC control)
"""
from __future__ import annotations

import argparse
import gzip
import sys

import pysam

STRONG = set('GC')


def classify(seq: str, i: int, ref: str, alt: str) -> tuple[str, str, str, str]:
    prev, nxt = seq[i - 1] if i > 0 else 'N', seq[i + 1] if i + 1 < len(seq) else 'N'
    ctx = prev + ref + nxt
    if ref == 'C' and nxt == 'G':
        cpg = 'C'
    elif ref == 'G' and prev == 'C':
        cpg = 'G'
    else:
        cpg = '.'
    if cpg != '.':
        transition = (cpg == 'C' and alt == 'T') or (cpg == 'G' and alt == 'A')
        klass = 'cpg_loss_ts' if transition else 'cpg_loss_tv'
    elif (alt == 'C' and nxt == 'G') or (alt == 'G' and prev == 'C'):
        klass = 'cpg_gain'
    else:
        klass = 'noncpg'
    ws = ('S' if ref in STRONG else 'W') + '>' + ('S' if alt in STRONG else 'W')
    return ctx, cpg, klass, ws


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--chrom', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    seq = pysam.FastaFile(a.fasta).fetch(a.chrom).upper()
    print(f'{a.chrom}: {len(seq):,} bp loaded', flush=True)

    n = n_mismatch = n_mono = 0
    with gzip.open(a.out, 'wt') as out:
        out.write('pos\tref\talt\tac\tan\tctx\tcpg\tklass\tws\n')
        # stdin: POS REF ALT AC_<pop>_unrel x5 AN_<pop>_unrel x5
        for line in sys.stdin:
            f = line.rstrip('\n').split('\t')
            pos, ref, alt = int(f[0]), f[1].upper(), f[2].upper()
            ac = sum(int(x) for x in f[3:8] if x != '.')
            an = sum(int(x) for x in f[8:13] if x != '.')
            n += 1
            if n % 1_000_000 == 0:
                print(f'{n:,} SNVs', flush=True)
            if seq[pos - 1] != ref:
                n_mismatch += 1
                continue
            if ac == 0 or ac == an:
                n_mono += 1  # private to the related samples, or fixed in the unrelated set
                continue
            ctx, cpg, klass, ws = classify(seq, pos - 1, ref, alt)
            out.write(f'{pos}\t{ref}\t{alt}\t{ac}\t{an}\t{ctx}\t{cpg}\t{klass}\t{ws}\n')
    print(f'done: {n:,} SNVs read, {n_mismatch:,} REF mismatches, {n_mono:,} monomorphic in unrelated', flush=True)


if __name__ == '__main__':
    main()
