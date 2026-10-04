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

## 2026-10-03 ~11:00: genome-wide R1e; R1f, R1g and R2a started

- **R1c:** human 28.15M / 28.16M CpGs lifted to hg38 (99.97%), chimp 25.38M / 26.60M (95.4%).
- **R1e genome-wide** (`tsv/R1/R1e_pd_*.tsv`): 2.80M ancestral CpGs with an event; P 2,207,489, D 590,760; overall P/D
  3.74.
  - **Positive control:** missense CpG loss P/D 14.6 (14.1–15.1) vs synonymous 4.98 (4.90–5.08). Stop-gains: 1,278 P
    vs 1 D. The machinery detects selection; a null in R2 will be meaningful.
  - **Conserved germline-high CpGs** (chimp sperm m > 0.8 and human regional m > 0.8; 1.26M sites), P/D relative to
    the phyloP < 0 background in the same bin: CGI 2.54, PLS 1.70, pELS 1.37, dELS 1.31, DNase-H3K4me3 1.29,
    CTCF-only 1.17. Coding 2.50 and conserved non-element sequence (phyloP ≥ 0) 2.22.
  - **Germline-low CpGs:** CGI 1.79, PLS 1.24, dELS 1.32, pELS 1.14.
  - **Reading:** elements carry excess P/D at fixed germline state, but conserved non-element sequence shows the same
    excess. Selection on sequence (TF motifs; F6) predicts this as well as causal methylation does. Not
    discriminating until (i) motif CpGs are separated (R1g), (ii) it is compared with non-CpG S>W P/D in the same
    elements (R1f), and (iii) R2 tests the slope against somatic m at fixed germline state.
  - **Caveat:** CpG and non-CpG P/D are not directly comparable (recurrent mutation at hypermutable CpGs;
    chr22 non-CpG S>W 4.4 vs CpG background 3.45). Compare element / background ratios within each class.
  - Sites without chimp coverage (nocov, 312k) have P/D 2.25, lower than every covered bin; check what they are
    (repeats? poorly aligned?) before using them.
- **R2a** (`scripts/R2a_loyfer_celltype.py`): Loyfer betas → per-CpG methylation for 82 tissue-specific cell types
  (207 sorted samples; cfDNA_* and CNVS-* samples excluded). The CpG order was checked from data: the CGI contrast is
  0.51 for chr1–22, X, Y, M vs 0.008 with chrM first. CGI median m 0.03, elsewhere 0.90. Done, all 22 chromosomes
  (`$scratch/R2/loyfer/chrN.{npz,tsv.gz}`).
- **R1g** (running): HOMER known motifs (472) on the ancestral-CpG-restored sequence. 2 Mb test: 60% of CpGs lie in ≥ 1
  hit; 17% have a hit ≥ 2 log-odds above threshold. Per-CpG n_hits / best motif / margin, so R2 can choose the cutoff.
- **R1f** (running): non-CpG P / D by direction × element. chr22 test: S>W 4.40, W>S 3.46, S>S 4.09, W>W 4.40, the
  expected gBGC pattern.

## 2026-10-03 ~11:30: phyloP background was circular; corrected comparison with non-CpG S>W

- **Bug in the 11:00 entry: the background was biased.** It was non-element sites with phyloP100way < 0. phyloP
  includes human, so a site with a human-lineage substitution scores lower, and phyloP < 0 enriches for D. Background
  P/D was biased down and the phyloP ≥ 0 remainder up (non-CpG S>W 2.56 vs 13.4). The `PD_rel_bg` values in the 11:00
  entry are inflated. Background is now every non-element, non-coding site, and phyloP is not used anywhere in R1e
  or R1f. R1f also labels Ensembl-canonical CDS as 'coding' (first priority), as R1e does. Tables rerun;
  `tsv/R1/R1e_pd_*.tsv` and `tsv/R1/R1f_gbgc_pd.tsv` are the corrected versions.
- **Corrected** (P/D relative to the non-element background of the same class; CpG split by conserved germline m):

| element | CpG loss, germline low | CpG loss, germline high | non-CpG S>W | non-CpG W>S |
|---|---|---|---|---|
| CGI | 1.53 | 2.34 (5.6k sites) | 1.50 | 0.93 |
| PLS | 1.06 | 1.56 (4.5k) | 1.30 | 1.12 |
| pELS | 0.98 | 1.26 | 1.17 | 1.00 |
| dELS | 1.13 | 1.20 | 1.14 | 0.97 |
| CTCF-only | (n = 11) | 1.07 | 1.10 | 0.95 |
| coding | 1.79 | 2.29 | 1.61 | 1.55 |

- **Reading (provisional; motif exclusion still pending):**
  - In enhancers and CTCF sites, and in germline-unmethylated CGIs / promoters (the usual CGI), CpG-loss P/D is no
    higher than for any S>W change in the same element. Part of the element S>W excess is gBGC: S>W / W>S is 1.33–1.86
    in elements vs 1.15 in background. So per-CpG selection on methylation is not detectable there, which fits the
    passenger prediction (F4) at per-CpG 2Ns ≳ 1.
  - The exception is **germline-methylated CGIs and promoters** (about 10k CpGs): 2.34 vs 1.50 and 1.56 vs 1.30. This
    is the only CpG-specific excess. Candidate explanations:
    - real selection on CpGs that are methylated in the germline and sit in CGIs / promoters (imprinted or
      germline-specific promoters?);
    - this CGI class is more constrained in general. The non-CpG comparator is not stratified by regional germline m
      yet; **next fix:** stratify the R1f non-CpG sites by human regional sperm m;
    - motif CpGs (R1g).
- Also next: characterise the nocov stratum (P/D 2.25), and the R2 slope against somatic m (R2a done).

## 2026-10-03 ~12:30: R1g done; motif-excluded R1e

- **R1g:** all 22 chromosomes done. 58.7–59.2% of ancestral CpGs lie in at least one HOMER known-motif hit (472
  motifs, half-sites dropped); 16.3–16.7% are in a strong hit (margin ≥ 2). Tables: `$scratch/R1/motif/chrN.motif.tsv.gz`.
- **R1e by motif mode** (`tsv/R1/R1e_pd_*.{exclude_any,exclude_strong,only_strong}.tsv`), P/D relative to background:

| conserved germline m | element | all | no motif hit | no strong hit | strong hit only |
|---|---|---|---|---|---|
| high | CGI | 2.34 | 2.30 | 2.33 | 2.37 |
| high | PLS | 1.56 | 1.55 | 1.53 | 1.69 |
| high | pELS | 1.26 | 1.22 | 1.24 | 1.35 |
| high | dELS | 1.20 | 1.17 | 1.19 | 1.26 |
| high | CTCF-only | 1.07 | 1.00 | 1.07 | 1.11 |
| low | CGI | 1.53 | 1.51 | 1.53 | 1.54 |
| low | PLS | 1.06 | 1.00 | 1.04 | 1.10 |

- **Reading:**
  - Motif overlap does not explain the germline-methylated CGI / promoter excess. It survives excluding every motif CpG.
  - Motif CpGs are only slightly more constrained than non-motif CpGs (e.g. high PLS 1.69 vs 1.55). HOMER hits are
    mostly not bound sites (they cover 59% of CpGs), so this is a weak version of the F6 test. Use ChIP / footprint
    support before claiming motif CpGs behave as "TF-selected".
  - The non-CpG S>W comparator (R1f) is not motif-filtered.
- **Still open, in order:**
  1. Stratify the non-CpG comparator by regional germline m (the remaining confounder for the germline-methylated
     CGI excess).
  2. Find out what the ~10k germline-methylated CGI / PLS CpGs are: imprinted DMRs, germline-specific promoters,
     intragenic CGIs?
  3. Characterise the nocov stratum.
  4. Run R2, the slope against somatic m from Loyfer at fixed germline m.

## 2026-10-03 ~15:30: germline-matched comparator and R2 (P/D vs somatic m)

Setup (`scripts/R2b_cpg_soma.py`, `R2c_pd_soma.py`; R1f rerun with `--human-sperm --loyfer`; tables in `tsv/R2/`):
- **Regional measures.** At fixed human losses there is no human CpG, so germline and somatic m are measured
  regionally: the mean of the other CpGs within ±500 bp, for P and D sites alike. Regional vs site somatic median m
  correlate at r = 0.93 (chr22).
- **CpG side:** conserved germline bin (chimp site m and human regional m agree); CpGs in strong HOMER hits excluded.
- **Non-CpG side:** S>W changes, human regional germline bin.
- Joint 1 Mb block bootstrap, 200 replicates.

**Germline-matched comparator** (excess = CpG P/D relative to its background ÷ S>W P/D relative to its background,
both in the same germline bin):

| germline | CGI | PLS | pELS | dELS | CTCF-only | coding |
|---|---|---|---|---|---|---|
| high | **1.95 (1.77–2.20)** | 1.22 (1.12–1.35) | 1.08 (1.03–1.13) | 1.05 (1.03–1.07) | 0.98 (0.91–1.08) | 1.44 (1.39–1.49) |
| low | 1.06 (1.01–1.11) | 0.83 (0.76–0.93) | 0.86 (0.82–0.92) | 0.95 (0.88–1.03) | – | 1.08 (1.00–1.19) |

- The germline-methylated CGI excess **survives the matched comparator**.
- Deficits as low as 0.83 also appear (germline-low PLS / pELS). Read ±15% as the scale of residual level differences
  between CpG and non-CpG P/D, not as signal. Only the germline-high CGI excess (and PLS, marginally) clearly exceeds it.
- Coding 1.44 is composition: CpG transitions in CDS are more often missense than the average S>W change. It is not a
  methylation effect.

**R2: ratio = CpG P/D ÷ S>W P/D, by regional somatic m, at fixed germline bin (non-coding):**
- **Germline-high, the bulk (≈ 1.0M CpG events):** flat.
  - soma_min: 0.92 / 0.89 / 0.94 / 0.91 across bins; top / bottom 0.98 (0.91–1.06). Elements only: 1.00 (0.95–1.06).
  - soma_med: 0.97 → 0.89 → 0.92 over 0.2–1. The top / bottom 0.73 is driven by the small 0–0.2 bin (196 D).
  - **No rise in CpG-loss constraint with somatic methylation**, the passenger prediction. Sensitivity bound: a rise
    of more than ~6% from somatically unmethylated to methylated is excluded (soma_min, elements).
- **Germline-low CpGs: a rise.** soma_med top / bottom 1.13 (1.08–1.19), elements 1.27 (1.16–1.39); soma_min
  elements 1.23 (1.07–1.44). These CpGs are unmethylated in sperm but methylated in soma. A candidate, but the obvious
  confounder is germline mutation rate: CpGs unmethylated in sperm and methylated in oocyte or early embryo mutate
  faster than the sperm stratum implies. **Needs the oocyte / embryo methylome (Okae 2014; Guo 2014) as a rate
  covariate before reading.**
- **Germline-high CGI:** soma_med 1.39 (0.5–0.8) → 1.83 (0.8–1), trend 1.31 (0.92–1.71); soma_min flat (1.15,
  0.85–1.59). The CGI excess does not clearly scale with somatic m.

**Overall:**
- For the large majority of CpGs (germline-methylated, regulatory or not), CpG loss is no more constrained than other
  S>W changes in matched context, and constraint does not track somatic methylation. This fits methylation as a
  passenger at per-CpG 2Ns ≳ 1 (the simulation's F4 / R2 outcome 1).
- Two candidate exceptions remain:
  1. ~5k germline-methylated CGI CpGs: excess ~1.9, not dose-dependent. What are they (imprinted gDMRs?
     germline-specific promoters)?
  2. germline-unmethylated, soma-methylated CpGs: rising. The oocyte methylome must come first.
- Next: download Okae 2014 oocyte / blastocyst WGBS (JGA / DDBJ?) or Guo 2014 (GSE49828); annotate the germline-high
  CGIs with imprinting control regions; translate the bound into per-CpG S with a stabilising-selection simulation
  (laptop).

## 2026-10-03 ~21:00: oocyte methylome (R0c) and R2 stratified by sperm × oocyte

- **R0c** (`scripts/qsub/R0c_oocyte.sh`): Zhu 2018 Nat Genet (GSE81233) single-cell WGBS, hg19 → hg38. All 63 files
  match the GEO sizes. Cells pooled per group:

| group | cells | CpGs on hg38 | cov ≥ 5 | mean m |
|---|---|---|---|---|
| MII oocyte | 36 | 25.99M | 19.16M | 0.57 |
| GV oocyte | 8 | 18.63M | 5.70M | 0.42 |
| ICM | 19 | 26.43M | 20.76M | 0.25 |

  - Strand handling checked: 99.2% of lifted positions sit on the C of an hg38 CpG.
  - MII is used as the oocyte stratum (regional ±500 bp, cov ≥ 5, same definition as sperm, on both the CpG and
    non-CpG sides).
- **R2 with sperm | oocyte germline bins** (`tsv/R2/R2c_*_oo.tsv`; R2b `--oocyte`, R1f `--oocyte`,
  R2c `--with-oocyte`).
- **Germline-matched comparator:** the sperm-methylated CGI excess does not depend on oocyte state:
  - high|high 1.91 (1.62–2.20);
  - high|mid 2.05 (1.60–2.68);
  - high|low 2.06 (1.52–2.84);
  - sperm-low CGI: 1.08 (low|low), 1.15 (low|mid), 1.47 (low|high, n = 663 events).
- **Germline-methylated (sperm-high), the bulk:** flat in every oocyte stratum.
  - high|high soma_med trend 0.96 (0.87–1.06), elements 0.94 (0.83–1.07); soma_min 0.98 (0.91–1.08).
  - high|low: 0.90 (0.72–1.17) and 1.05 (0.97–1.14).
  - The passenger reading is robust to oocyte methylation.
- **Sperm-low CpGs:**
  - **low|low** (unmethylated in both gametes): the rise **persists but weakens**.
    - soma_med trend 1.12 (1.03–1.24), elements 1.17 (1.00–1.39); sperm-only was 1.13 / 1.27.
    - It is a top-bin effect: the ratio is about 1.0 for soma m 0–0.8 and 1.11–1.15 at 0.8–1.
    - soma_min is not significant (1.05, 0.94–1.17).
  - **low|mid:** similar (1.17, 1.05–1.30).
  - **low|high** (unmethylated in sperm, methylated in oocyte, the pattern of maternal gDMRs / imprinting control
    regions): elements 0.8–1 vs 0.5–0.8 = 1.86 (1.39–2.58). Few events; candidate imprinted loci.
- **Reading:**
  - Oocyte methylation does not explain the sperm-low rise.
  - The remaining rate route is the **post-implantation epiblast / early PGC window**. The germline is re-methylated
    to a soma-like state there before PGC erasure, so germline mutations arising in that window see soma-like
    methylation, and "somatically methylated" partly means "methylated in part of the germline cycle".
  - Testing that needs human epiblast or early PGC methylomes (e.g. Guo 2015 Cell, Tang 2015 Cell; check access),
    or a mutation-rate model fitted to de novo mutations (e.g. DNM rates at CpGs by soma m at fixed gamete m).
  - Until then, the sperm-low rise is not evidence for causal methylation.
