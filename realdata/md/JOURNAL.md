# Journal (cluster session)

## 2026-10-03: R0 downloads and R1a started (defaults accepted)

- The user accepted the DESIGN_PLAN defaults: sperm as the germline methylome (oocyte as a later sensitivity check, not
  downloaded yet), Loyfer atlas as the somatic readout, and the R2 simulation on the laptop.
- GSE30340 has both species: human sperm (hg18) and chimp sperm (panTro2), each with CpG methylation, coverage and HMRs.
  The conserved-germline-state filter in R1 is therefore possible. The chimp data reach hg38 through panTro2→hg19→hg38
  chains.
- R0a (job 15025255): all downloads in one job; see `README.md`.
- R1a (`scripts/R1a_snv_context.py`): 21 of 22 chromosomes done. chr22: 925,730 biallelic SNVs, 0 REF mismatches,
  AN = 5,008, 199,083 CpG transitions (~22%).
- **The local 1000G chr1 VCF is corrupt.** `gzip -t` gives a length error, and bcftools stops at chr1:190,673,509.
  The local file has exactly the size of the EBI file. R0b (job 15025308) re-downloads it to `R0/kg/`, checks it and
  byte-compares the two. R1a chr1 (job 15025309) is held on it and uses the new copy automatically.
- Job-script gotchas on this cluster: conda activate breaks under `set -u` (wrap it in `set +u`), and `module` needs
  `source /u/local/Modules/default/init/bash` in batch jobs.
- R0b result: the fresh chr1 passes `gzip -t`. It differs from the local copy in 199 bytes, starting at byte
  1,762,378,770, so the damage is local and on disk. md5 589762a4… matches the clean copy latent_genos made on
  2026-09-30 (`latent_genos/reference/1000G_30x/`, A04h). This was already known and written up in latent_genos and
  pool_design; the re-download duplicates it. Now also documented at the source:
  `demux_benchmark/pool_design/vcf/1000G/README_CHR1_CORRUPT.md`. The merged VCFs there (`1000G.merged*.vcf.gz`) have
  no chr1 past 190.67 Mb either.

## 2026-10-03: R1b–R1e (ancestral CpGs → P/D)

- **R1b** (`scripts/R1b_ancestral_cpg.py`): an ancestral CpG is CG in both panTro6 and rheMac10, as adjacent axt
  columns. Its hg38 state is retained / lost_C (TG) / lost_G (CA) / other. Joined to R1a: polymorphic losses at
  retained sites; at lost sites, the ancestral allele still segregating. chr22: 223,101 ancestral CpGs, 7.1% lost on
  the human lineage by transition (≈ 0.6% divergence × ~12 CpG hypermutability), and 24% of retained sites
  polymorphic for a loss allele at n = 5,008.
- **R1c**: human sperm hg18→hg38 took 3 min. Chimp panTro2→hg19 in a single pass was on course for ~9 h, so it now runs
  as 40 chunks. The UCSC liftOver binary needs glibc ≥ 2.25, so the conda `babel` env's liftOver is used
  (`config/paths.sh`).
- **Germline state at fixed losses.** A fixed human loss has no human CpG, so human sperm methylation cannot be
  measured at the site. Chimp keeps the CpG at every ancestral site, so the per-site germline stratum is **chimp sperm
  m** (cov ≥ 5). The human-chimp conservation filter uses the human **regional** sperm m (other CpGs within ±500 bp,
  cov ≥ 5) in the same bin.
- **R1d** (`scripts/R1d_annotate.py`): sperm m (site + regional, both species), cCRE V3 class, CGI, phyloP100way,
  coding consequence on the Ensembl-canonical CDS. Bug found in testing: consequences were computed on hg38, where
  fixed losses are already T / A, so coding D was 0. They are now computed on the ancestral codon (CpG restored).
- **R1e** (`scripts/R1e_pd_tables.py`): P/D by germline bin, element, germline × element (relative to the phyloP < 0
  background in the same bin), and the coding control. 1 Mb block bootstrap.
- **chr22 test (code check; not stratified by germline, chimp sperm not lifted yet):**
  - **The positive control works.** Missense CpG loss P/D 14.1 (11.9–16.2) vs synonymous 5.6 (4.8–6.5);
    stop 29 (n = 30).
  - Elements: background 3.3, dELS 4.4, pELS 5.0, PLS 5.8, CGI 7.2, phyloP ≥ 0 non-element 10.1.
  - Do not read the element ordering yet. Without germline stratification, CGI / PLS P/D can reflect a lineage change
    in germline methylation as well as selection (REVIEW §3, F3).
- Queued: R1c chimp chunks (15025528) → merge (15025529) → R1d all chromosomes (15025530; also held on R1b chr1) →
  R1e genome-wide tables to `tsv/R1/` (15025531).
- Not yet: the TF-motif exclusion (JASPAR scan; needed before R2), the gBGC W↔S control (needs ancestral states at
  non-CpG sites), the gnomAD LOEUF split of the coding control, and oocyte methylation.
