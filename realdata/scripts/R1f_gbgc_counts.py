"""R1f: gBGC control. Non-CpG P and D counts by substitution direction (W>S, S>W, W>W, S>S) and element, one chromosome.

CpG loss is S>W, and gBGC raises P/D for S>W changes (W>S alleles are favoured in transmission, so S>W alleles
segregate but fix less often). R1 therefore compares CpG-loss P/D with non-CpG S>W P/D in the same element class.

Sites: hg38 positions where chimp and macaque carry the same base (that base is ancestral) at the site and both
neighbours, with the axt alignments as in R1b. Excluded: sites in an ancestral CpG (C before G, or G after C), and
events whose derived allele would create a CpG (derived C before an ancestral G, or derived G after an ancestral C).
  D  hg38 base differs from the ancestral base, with no SNV restoring it
  P  1000G SNV with one allele ancestral (derived count = AC if REF is ancestral, AN - AC otherwise); 0 < dac < AN
Elements: coding (Ensembl-canonical CDS) > CGI > cCRE class (PLS > pELS > dELS > CTCF-only > DNase-H3K4me3) >
background. phyloP is not used (it includes human, so it is circular with D).
Output: counts by element x direction x 1 Mb block (P, Ps = P with derived count >= 2, D); with --human-sperm /
--loyfer also by regional germline and somatic m bins (R2b definitions; the matched comparator for R1e / R2).
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from R1b_ancestral_cpg import load_axt  # noqa: E402
from R1d_annotate import CCRE_ORDER, interval_mask  # noqa: E402
from R2b_cpg_soma import GERM_BINS, GERM_LABELS, SOMA_BINS, SOMA_LABELS, human_sperm_arrays, loyfer_arrays, regional  # noqa: E402

STRONG = np.zeros(256, dtype=bool)
STRONG[[ord('C'), ord('G')]] = True


def cds_mask(gtf: str, chrom: str, L: int) -> np.ndarray:
    st, en = [], []
    with gzip.open(gtf, 'rt') as fh:
        for line in fh:
            if not line.startswith(chrom + '\t'):
                continue
            f = line.split('\t')
            if f[2] == 'CDS' and 'Ensembl_canonical' in f[8] and 'gene_type "protein_coding"' in f[8]:
                st.append(int(f[3]) - 1)
                en.append(int(f[4]))
    return interval_mask(np.array(st, dtype=np.int64), np.array(en, dtype=np.int64), L)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--chrom', required=True)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--chimp-axt', required=True)
    ap.add_argument('--macaque-axt', required=True)
    ap.add_argument('--snv', required=True)
    ap.add_argument('--ccre', required=True)
    ap.add_argument('--cgi', required=True)
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--human-sperm', help='R1c human.hg38.bed; adds germ = regional human sperm m bin (+-500 bp)')
    ap.add_argument('--oocyte', help='R0c MII oocyte bed; adds oo = regional oocyte m bin (+-500 bp)')
    ap.add_argument('--loyfer', help='R2a chrN.tsv.gz; adds soma_med / soma_min = regional somatic m bins (+-500 bp)')
    a = ap.parse_args()

    h = np.frombuffer(pysam.FastaFile(a.fasta).fetch(a.chrom).upper().encode(), dtype=np.uint8)
    L = len(h)
    pt, _ = load_axt(a.chimp_axt, a.chrom, L)
    rm, _ = load_axt(a.macaque_axt, a.chrom, L)
    acgt = np.isin(pt, [ord(x) for x in 'ACGT'])
    anc = np.where(acgt & (pt == rm), pt, 0).astype(np.uint8)  # 0 = unknown
    C, G = ord('C'), ord('G')
    prev = np.concatenate([[0], anc[:-1]])
    nxt = np.concatenate([anc[1:], [0]])
    ok = (anc > 0) & (prev > 0) & (nxt > 0)
    ok &= ~((anc == C) & (nxt == G)) & ~((anc == G) & (prev == C))
    print(f'{a.chrom}: {ok.sum():,} usable non-CpG ancestral sites', flush=True)

    # polymorphism
    snv = pd.read_csv(a.snv, sep='\t', usecols=['pos', 'ref', 'alt', 'ac', 'an'])
    snv = snv.drop_duplicates('pos', keep=False)
    i = snv.pos.values - 1
    keep = ok[i]
    snv, i = snv[keep], i[keep]
    a_b = anc[i]
    ref = np.array([ord(x) for x in snv.ref], dtype=np.uint8)
    alt = np.array([ord(x) for x in snv.alt], dtype=np.uint8)
    ref_anc, alt_anc = ref == a_b, alt == a_b
    sel = ref_anc | alt_anc
    der = np.where(ref_anc, alt, ref)
    dac = np.where(ref_anc, snv.ac.values, snv.an.values - snv.ac.values)
    seg = sel & (dac > 0) & (dac < snv.an.values)
    restoring = np.zeros(L, dtype=bool)
    restoring[i[alt_anc]] = True  # hg38 carries a derived base whose ancestral allele still segregates
    P = pd.DataFrame({'i': i[seg], 'anc': a_b[seg], 'der': der[seg], 'dac': dac[seg]})

    # fixed differences
    di = np.where(ok & (h != anc) & np.isin(h, [ord(x) for x in 'ACGT']) & ~restoring)[0]
    D = pd.DataFrame({'i': di, 'anc': anc[di], 'der': h[di], 'dac': -1})
    ev = pd.concat([P.assign(kind='P'), D.assign(kind='D')], ignore_index=True)
    # derived allele must not create a CpG
    ev = ev[~(((ev.der == C) & (nxt[ev.i] == G)) | ((ev.der == G) & (prev[ev.i] == C)))]
    ev['direction'] = (np.where(STRONG[ev.anc], 'S', 'W').astype(object) + '>'
                       + np.where(STRONG[ev.der], 'S', 'W').astype(object))

    # elements
    cc = pd.read_csv(a.ccre, sep='\t', header=None, usecols=[0, 1, 2, 5], names=['chrom', 'start', 'end', 'cls'])
    cc = cc[cc.chrom == a.chrom]
    code = np.zeros(L, dtype=np.int8)
    for r, c in enumerate(CCRE_ORDER[::-1], start=1):
        sub = cc[cc.cls.str.split(',').str[0] == c]
        code[interval_mask(sub.start.values, sub.end.values, L)] = r
    names = np.array(['none'] + CCRE_ORDER[::-1])
    cgi = pd.read_csv(a.cgi, sep='\t', header=None, usecols=[1, 2, 3], names=['chrom', 'start', 'end'])
    cgi = cgi[cgi.chrom == a.chrom]
    in_cgi = interval_mask(cgi.start.values, cgi.end.values, L)
    el = names[code[ev.i.values]].astype(object)
    el = np.where(el == 'none', 'background', el)
    el = np.where(in_cgi[ev.i.values], 'CGI', el)
    el = np.where(cds_mask(a.gtf, a.chrom, L)[ev.i.values], 'coding', el)
    ev['element'] = el
    ev['block'] = a.chrom + ':' + (ev.i // 1_000_000).astype(str)
    ev['P'] = (ev.kind == 'P').astype(int)
    ev['Ps'] = ((ev.kind == 'P') & (ev.dac >= 2)).astype(int)
    ev['D'] = (ev.kind == 'D').astype(int)
    keys = ['element', 'direction', 'block']
    if a.human_sperm:
        hp, hm = human_sperm_arrays(a.human_sperm, a.chrom, h.tobytes().decode())
        ev['germ'] = pd.cut(regional(hp, hm, ev.i.values), GERM_BINS, labels=GERM_LABELS).astype(str)
        keys.append('germ')
    if a.oocyte:
        op, om = human_sperm_arrays(a.oocyte, a.chrom, h.tobytes().decode())
        ev['oo'] = pd.cut(regional(op, om, ev.i.values), GERM_BINS, labels=GERM_LABELS).astype(str)
        keys.append('oo')
    if a.loyfer:
        lp, med, mn = loyfer_arrays(a.loyfer, 40)
        ev['soma_med'] = pd.cut(regional(lp, med, ev.i.values), SOMA_BINS, labels=SOMA_LABELS).astype(str)
        ev['soma_min'] = pd.cut(regional(lp, mn, ev.i.values), SOMA_BINS, labels=SOMA_LABELS).astype(str)
        keys += ['soma_med', 'soma_min']
    out = ev.groupby(keys)[['P', 'Ps', 'D']].sum().reset_index()
    out.to_csv(a.out, sep='\t', index=False, compression='gzip')
    tot = ev.groupby('direction')[['P', 'D']].sum()
    print((tot.assign(PD=tot.P / tot.D)).to_string(), flush=True)
    print('done', flush=True)


if __name__ == '__main__':
    main()
