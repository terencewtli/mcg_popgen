# realdata/ (cluster-owned)

Real-data arm of mcg_popgen, owned by the cluster session (`md/ENDPOINT.md` repo rules: the laptop never edits this
directory; the cluster works in this clone and pulls before committing).

- `md/REVIEW_2026-10-02.md`: math and logic review of the simulation phase (corrections to F2 and the real-N job)
- `md/DESIGN_PLAN.md`: real-data design (R0 data, R1 CpG-loss P/D, R2 P/D vs somatic methylation, R3 o/e +
  conservation, R4 HPRC2), compute and time estimates, open decisions

Large data (methylomes, VCFs, alignments) stays in `/u/project/cluo_scratch/terencew/claude/mcg_popgen/` and is never
committed. Record paths and md5s here as they are downloaded.

## Data on disk (R0, 2026-10-03)

Paths are set in `config/paths.sh`. Downloads go to `/u/project/cluo_scratch/terencew/claude/mcg_popgen/R0/`, and md5s are written to `R0/MD5SUMS.txt` by R0a.

| data | path under `R0/` | source |
|---|---|---|
| UCSC hg38 fasta (all contigs; indexed as `hg38.fa`) | `ucsc/hg38.fa.gz` | hgdownload bigZips (md5 checked) |
| hg38 × rheMac10 reciprocal-best axt | `ucsc/hg38.rheMac10.rbest.axt.gz` | hgdownload vsRheMac10 (md5 checked) |
| CpG islands, phyloP100way | `ucsc/cpgIslandExt.txt.gz`, `ucsc/hg38.phyloP100way.bw` | hgdownload |
| liftOver binary + chains (hg18→hg38; panTro2→hg19→hg38) | `../tools/liftOver`, `liftover/` | hgdownload |
| Human (hg18) + chimp (panTro2) sperm WGBS: CpG methylation, coverage, HMRs | `sperm/GSE30340_*` | GEO GSE30340 (Molaro 2011) |
| Loyfer 2023 atlas, 253 hg38 `.beta` (wgbstools CpG order, 29,401,795 CpGs) + series matrix | `loyfer/` | GEO GSE186458; URL list `tsv/R0/loyfer_hg38_beta_urls.tsv` |
| JASPAR 2024 CORE vertebrates non-redundant (MEME) | `jaspar/` | jaspar.elixir.no |
| gnomAD v4.1 constraint metrics | `gnomad/` | gnomAD GCS bucket |
| 1000G chr1 (re-download; the local copy fails `gzip -t`; same md5 as the earlier clean copy in `latent_genos/reference/1000G_30x/`) | `kg/` | EBI 20220422_3202_phased_SNV_INDEL_SV |

Already on disk (not copied): 1000G NYGC 30x per-chromosome VCFs, hg38 × panTro6 rbest axt and ENCODE V3 cCREs (paths in `config/paths.sh`).

## Outputs

| step | output | status |
|---|---|---|
| R1a | `R1/snv/chrN.snv.tsv.gz` (biallelic SNVs, 2,504 unrelated: AC / AN, hg38 trinucleotide, CpG class, W/S) + `chrN.summary.tsv` | chr2–22 done 2026-10-03; chr1 rerun on the R0b copy |
| R1b | `R1/anc/chrN.anc.tsv.gz` (ancestral CpGs: hg38 state, SNV, derived count) | chr2–22 done; chr1 queued |
| R1c | `R1/sperm/{human,chimp}.hg38.bed` (sperm m, cov on hg38) | running (chimp in 40 chunks) |
| R1d | `R1/annot/chrN.annot.tsv.gz` (+ sperm, cCRE, CGI, phyloP, coding consequence) | queued |
| R1e | `tsv/R1/R1e_pd_*.tsv` (P/D tables, small; committed) | queued |
