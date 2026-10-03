# Real-data design plan (cluster arm)

Written 2026-10-02 after reviewing the simulation phase (`realdata/md/REVIEW_2026-10-02.md`). Nothing run yet.
It builds on SUMMARY R1–R7 and REAL_DATA analyses A–D, reordered by robustness and cost.

## The question, restated for real data

Is DNA methylation at regulatory CpGs under selection **as methylation** (causal: fitness reads the methylation
level), or is it a **passenger** of TF binding (fitness reads occupancy; methylation follows)?

The simulations say the cleanest discriminator is CpG-loss variants **outside TF motifs**:
- **passenger:** neutral (P/D = 1, F4);
- **causal:** depleted (P/D > 1, F5);
- **motif variants:** depleted either way (F6), so they cannot discriminate.

## Design principle: avoid detection bias by not using detection

E4 shows meQTL-based tests need an ascertainment-matched null that is hard to build in real data. But the effect of
a **CpG-destroying** allele on methylation is known *a priori*, without genotypes: it removes that CpG's methylation.
Its size is the CpG's methylation level (plus neighbour effects). So the causal model predicts a dose–response that
can be tested on **all** CpG-loss variants, with methylation measured in reference methylomes, and no meQTL calling.

## Analyses, in order

### R0. Data assembly (1–2 days; mostly public data, partly already on disk)

| data | use | source / status |
|---|---|---|
| Human polymorphism | P (CpG-loss SNVs, SFS) | 1000G NYGC 30x, 3,202 samples; **on disk** (`/u/project/cluo/terencew/demux_benchmark/pool_design/vcf/1000G/by_chrom/`; its chr1 is corrupt past 190.67 Mb, re-downloaded by R0b). Use the 2,504 unrelated samples. At n ≈ 5,000 haplotypes, methylated-CpG transitions are far from saturation (θ_CpG ≈ 4×10⁻³ → ~4% of methylated CpGs polymorphic), so recurrence is a small correction, not a blocker. |
| Human–chimp–macaque alignment | D on the human lineage; polarisation | hg38×panTro6 reciprocal-best axt **on disk** (mcg_evo `reference/ucsc/`); add hg38×rheMac10 rbest axt from UCSC (~1.5 GB) |
| Germline methylation | stratify the mutation rate | human and chimp sperm WGBS (Molaro 2011, GEO GSE30340: verify); oocyte / early-embryo WGBS (Okae 2014; Guo 2014: verify access) |
| Somatic methylation | the functional readout | human WGBS atlas by cell type (Loyfer 2023, GEO GSE186458: verify), ENCODE / Roadmap tissue WGBS |
| Regulatory annotation | element classes | ENCODE V3 cCREs (**on disk** in mcg_evo); CGIs (UCSC) |
| TF motifs | exclude / separate motif CpGs | JASPAR 2024 + FIMO (or HOMER known motifs) on hg38 |
| Positive-control sites | prove the P/D machinery detects selection | coding CpGs: CpG-loss alleles that are missense in constrained genes (gnomAD LOEUF) vs synonymous |

### R1. Analysis B core: CpG-loss P/D by germline methylation × element class (≈ 1 week)

- **Unit.** Ancestral CpGs (CpG in chimp and macaque, aligned in all three species).
  - **D:** human-lineage loss (TpG/CpA in hg38).
  - **P:** a segregating CpG-destroying allele in 1000G at an ancestral CpG that is CpG in hg38. Use derived
    allele count ≥ 1, split singleton vs non-singleton.
- **Strata.**
  - Germline methylation of the CpG: sperm m in bins (< 0.2, 0.2–0.8, > 0.8).
  - **Require conserved germline state in human and chimp sperm.** This is the main confounder missing from the
    simulation docs (REVIEW §3, F3 row): a lineage-specific change in germline methylation makes P/D ≠ 1 with no
    selection.
  - Crossed with element class: PLS, pELS, dELS, CTCF-only, CGI, and non-cCRE background.
- **Neutral comparator.** The same germline stratum in non-cCRE, non-conserved background (phyloP < 0), as the
  denominator. Report (P/D in stratum) ÷ (P/D in matched background).
- **gBGC control.** Non-CpG W→S and S→W P/D in the same elements (gBGC raises S→W P/D, and CpG loss is S→W).
- **Positive control.** Coding CpGs, missense vs synonymous CpG-loss. Missense P/D must come out > 1, otherwise
  the pipeline cannot see selection and a null in R2 means nothing.
- **Reading.**
  - P/D = 1 in elements: no per-CpG selection detectable.
  - P/D > 1: selection on CpGs (or something correlated). Do **not** convert to S with the P0 formula; under
    stabilising selection it overestimates S by orders of magnitude (E2b: 8 vs 158).
- **Predicted magnitudes** at this n (REVIEW §1 table): ~1.6 / 2.8 / 11 for S = 1 / 2 / 4 under directional
  selection; less under stabilising selection.
- **Power.** Human-lineage CpG losses number in the hundreds of thousands genome-wide and tens of thousands within
  cCREs. Polymorphic CpG-loss alleles number in the millions. A P/D shift of 1.2 per stratum is detectable;
  confounders, not counts, are the limit.

### R2. The discriminator: P/D vs somatic methylation at fixed germline methylation (≈ 1 week after R1)

- **Logic.**
  - Within one germline-m stratum the mutation rate is fixed.
  - Under causal methylation, the fitness effect of losing a CpG scales with that CpG's **somatic** methylation
    (its contribution to the methylation load), so P/D rises with somatic m.
  - Under passenger methylation, P/D is flat in somatic m.
- **Design.**
  - Within germline stratum × element class, bin CpGs by somatic m (cell-type maximum and median from the WGBS
    atlas, and in the cell type where the element is active).
  - Fit P/D ~ somatic m (bootstrap over elements).
  - Exclude motif-overlapping CpGs (F6). Analyse them separately as the "TF selected" reference, which should be
    depleted regardless of somatic m.
- **Confounders and controls.**
  - Somatic and germline m correlate, so stratify germline m finely (or model the rate continuously).
  - Somatic m correlates with element type and GC, so carry the gBGC (W↔S) control into each bin.
  - Somatically hypomethylated CpGs cluster in TF-bound regions, so control motif density.
- **Simulation support (laptop or cluster, small).** Turn the predicted curve into numbers: extend P0 / E2b so
  per-CpG S scales with the CpG's somatic contribution under C (stabilising, κ grid), against flat under P. This
  is a new, cheap simulation, and it is the one most directly tied to a real-data figure.
- **Possible outcomes:**
  1. Flat everywhere, with the positive control working: no evidence for causal methylation at per-CpG 2Ns ≳ 1.
     State the sensitivity bound from the simulation.
  2. Rising in a specific element class (e.g. CGI promoters, CTCF sites): candidate causal-methylation elements.
  3. Rising everywhere, including background: points to a confounder (gBGC or rate misspecification). The controls
     are designed to catch this.

### R3. Secondary: joint o/e + conservation (≈ 3 days; uses R0 data)

- REVIEW §2: at matched CpG o/e, selection-maintained CpGs turn over 1.5–1.9× faster (less conserved to macaque)
  than mutation-maintained ones.
- **Test.** Within germline-unmethylated elements (where passive retention explains o/e), compare macaque
  conservation with the mutation-only expectation. Estimate the gain rate γ empirically per context from
  background.
- **Read only alongside R1.** Its identification needs γ, r and α, so it is a consistency check, not a primary test.

### R4. Optional, later: haplotype-resolved methylation (HPRC2)

- Fit the methylation map (protection width λ, CpG-density feedback a₁) from per-haplotype 5mC (REAL_DATA §4.1). This
  replaces the hand-set map before any further simulation.
- Measure allele-specific effects of CpG-SNPs directly, including **neighbour** effects (the a₁ term), which R2
  treats as zero.

### Deprioritised

- **E4-style detected-meQTL MAF test:** keep only as a check (REVIEW §3). Its real-data null is hard to build.
- **Analysis D** (null for the Ma et al. sign test): wait. Our reproduction found the sign-test results and the
  promoter table not reproducible (mcg_evo `docs/REPRODUCTION.md` R04b / R05).
- **Real-N SLiM check:** re-scope to polymorphism snapshots for C κ=10 only, and run in the background (REVIEW §4).

## Compute and time

| step | compute | wall / analyst time |
|---|---|---|
| R0 downloads + liftover / indexing | < 20 CPU-h; ~50–100 GB scratch (methylomes are the bulk) | 1–2 days |
| R1 site classification + P/D tables | ~20–50 CPU-h (genome-wide per-site joins; chromosome job array) | ~1 week |
| R2 + its small simulation | < 20 CPU-h | ~1 week |
| R3 | < 10 CPU-h | ~3 days |
| Real-N SLiM (background) | ~100–300 CPU-h as an array (polymorphism-only) | unattended |

Total: about 3–4 weeks to a first answer on R1–R2, mostly analyst time; compute is modest.

## Repo layout (cluster-owned, per `md/ENDPOINT.md`)

- `realdata/README.md`: data paths and md5s (large data stays in `/u/project/cluo_scratch/...`, never committed).
- `realdata/scripts/`: `R0a_*` downloads, `R1a_classify_cpgs`, `R1b_pd_tables`, `R2a_pd_vs_soma`, `R3a_oe_conservation`.
- `realdata/md/`: this plan, the review, and `JOURNAL.md` for the cluster session.
- `realdata/tsv/`, `realdata/pdf/`: small tables and figures.

## Open decisions for the user

1. **Which germline methylome is authoritative?**
   - sperm only (the paternal-age-dominated mutation process makes sperm the main source);
   - or sperm + oocyte (CpG deamination in the female germline is a minority contribution).

   Proposed default: sperm, with oocyte as a sensitivity check.
2. **Somatic readout.** A cell-type-resolved atlas (Loyfer 2023) or bulk tissues? The atlas is better for R2, because
   causal methylation should matter in the cell type where the element acts.
3. **R2's simulation.** Run it on the laptop (it owns `scripts/`, per repo rules), with the cluster providing
   observed bins.
