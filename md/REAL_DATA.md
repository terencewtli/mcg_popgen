# From toy simulations to real data

Written 2026-10-01. Status: ideas and design, nothing run on real data. Citations were verified against PubMed or
journal pages by a search agent unless marked (unverified). Full list in `md/LITERATURE.md`.

## 1. What is "toy" about E0–E2, and what is not

**Toy (the biology):** one 300-bp locus, one 8-bp motif, hand-set occupancy and methylation curves, no trans
factors, no methylation → TF arrow, no recombination.

**Not toy (the population-genetic regime):** results depend on 2Ns, mutation-rate ratios and time measured in
substitutions. Q-rescaling preserves all three, and that invariance was checked (E0; E1 at Q = 1, 2, 4). N = 100 is a
computational convenience, not an approximation.

**The one inflated regime parameter:** θ per site is 0.024, about 20× human. This speeds up polymorphism but
increases recurrent CpG mutation (4N·rα ≈ 0.16 vs ≈ 0.005 in humans). Divergence and footprint results (E1, P0
divergence) do not depend on θ. Polymorphism statistics (E2 MAF, P/D) do, and should be re-run at realistic θ before
any comparison with real data.

## 2. Will the effects change at real scale?

| effect | expected at real scale | why |
|---|---|---|
| Passive CpG retention (E1) | **similar magnitude (a few-fold)** | Equilibrium composition depends on r and on the germline protection width, not on N or θ. Real r ≈ 10–15 (lower than E1's 20, so a smaller footprint); real low-methylated regions are hundreds of bp (wider than λ = 25, so a larger footprint). CGI o/e 0.6–1 vs genome ~0.2 is consistent with this (Cohen 2011). |
| Selection signals (P0 P/D, E2 MAF) | **same per-site 2Ns, far more power** | Millions of CpGs and thousands of meQTLs make P/D shifts of 1.2–1.5 detectable. |
| Interpretability | **worse** | The real genome mixes passenger and causal elements; genome-wide averages blur them. Stratify by element class. |

**At scale the binding constraint is confounders, not power.** Each of these produces shifts of the same size as the
signal:
1. **Recurrent CpG mutation.** CpG polymorphism saturates in samples of order 10⁵. Use moderate n or model recurrence.
2. **gBGC.** B = 4Ne·b ≈ 0.3–1 behaves like S favouring CpGs. Control with non-CpG W↔S sites at the same loci.
3. **Demography.** Expansion inflates rare variants. Compare against a matched neutral class rather than the
   standard neutral model.
4. **LD.** Lead meQTL SNPs are often not causal, which dilutes effect–frequency relationships. Fine-map, or restrict to
   variants with direct mechanisms (CpG-SNPs).
5. **Ascertainment (the most dangerous).** An meQTL is detected only if 2p(1−p)β² exceeds the power threshold, so
   discovered rare variants have large effects *with no selection* (the "trumpet" from detection alone; cf. Souaiaia
   et al. 2026 for complex traits). Any effect–MAF test must model detection power.

**What "scaling up" should mean** (rather than a longer sequence):
1. Many loci with heterogeneous parameters (element types, germline methylation, CpG density). Cheap with
   origin-fixation for divergence.
2. Realistic θ and human demography for polymorphism: an msprime/stdpopsim neutral background with SLiM for selected
   loci.
3. **A measurement layer**: bisulfite/array sampling noise, realistic sample sizes, meQTL calling with real power. This
   reproduces ascertainment and turns the simulations into a null for real meQTL analyses. It matters most.
4. A map fitted to data rather than set by hand (§4).

## 3. Four real-data analyses the model sets up

### A. Passive retention: germline methylation, not somatic activity, predicts CpG excess (from E1)
- **Prediction:** CpG o/e and cross-species CpG conservation at regulatory elements track **germline** hypomethylation
  (sperm, oocyte, early embryo). Elements hypomethylated only in soma should show no passive excess. Any excess there
  needs selection or gBGC.
- **Data:** germline methylomes (Molaro 2011 Cell, human/chimp sperm; Okae 2014 PLoS Genet, WGBS of oocytes, sperm and blastocysts; Guo 2014 Nature, single-cell RRBS of gametes and early embryos); tissue methylomes;
  ENCODE cCREs; human–chimp–macaque alignments (or Zoonomia/Cactus).
- **Readout:** CpG o/e and P(CpG conserved to macaque) against germline m, stratified by somatic-only vs
  germline-active elements.

### B. McDonald–Kreitman-type P/D test for CpG loss (from P0)
- **Prediction:** a change in mutation rate (germline m) leaves CpG-loss P/D at 1; selection with S = 1, 2, 4 gives
  P/D ≈ 1.5, 2.4, 7.8. o/e and divergence alone lie on a ridge in (m, S), which P/D breaks.
- **Data:** CpG-loss polymorphism from 1000G/gnomAD at moderate n (or with recurrence modelled); CpG-loss divergence on
  the human lineage (human–chimp–macaque); methylation-aware CpG mutation rates (Chandra & Gao 2026; gnomAD's
  methylation-adjusted mutability (Karczewski 2020: CpG transitions in 3 methylation bins, but from **somatic**
  methylation averaged over 37 tissues, not germline. That mismatch is exactly what analysis A tests)).
- **Design:** stratify by germline methylation class × element class. The neutral comparator is the same germline-m
  class in non-functional sequence. The gBGC control is non-CpG W→S vs S→W in the same elements.

### C. meQTL effect size vs allele frequency, by mechanism class (from E2)
- **Existing method:** for complex traits and expression, β² ∝ [2p(1−p)]^S estimates negative selection (BayesS,
  Zeng 2018 Nat Genet; α-model, Schoech 2019 Nat Commun; for eQTLs, Glassberg 2019 Genetics: constraint is real but
  weak).
- **For meQTLs specifically:**
  - Hatcher, Hemani, …, Min (bioRxiv 2021, 10.1101/2021.11.25.469994; no journal version found) estimated a
    BayesS-type S across DNAm sites and found both negative and positive selection. This is the closest prior work:
    **pooled meQTL effect–MAF selection inference is not new.**
  - GoDMC itself (Min 2021 Nat Genet, 32,851 people) has no effect–MAF analysis.
  - Shang 2023 Nat Commun (GENOA; first author unverified) reports an effect–MAF correlation of about −0.36,
    descriptively.
- **What is new here: a prediction by mechanism class**, plus a simulated null that includes recurrence and
  ascertainment.
  - **CpG-SNPs** (create or destroy a CpG; direct sequence → methylation effect): effect–MAF coupling is expected
    under C and not under P.
  - **Motif-breaking SNPs:** coupled under both models, so not diagnostic alone.
  - **The test is the difference in S between the two classes, not a shared S.**
- **Required care:** CpG-SNPs are hypermutable, so their MAF baseline is inflated by recurrence and must come from a
  mutation-matched neutral class. Detection power must be modelled (§2.5). Array meQTLs (450k/EPIC, e.g. GoDMC) cover
  CpGs biased toward promoters and CGIs.
- **Data:** GoDMC-scale array meQTLs; WGBS meQTL panels; long-read haplotype-resolved methylation (§4).

### D. A mutation-aware null for the Ma/Fraser sign test (from E1 + P0)
- **Issue:** their test assumes human- and chimp-biased changes are a fair coin, corrected by a global background
  fraction. CpG losses outnumber gains, and losses happen mostly at methylated CpGs. The coin is therefore biased by
  mutation, and the bias differs between gene sets with different CpG content and germline methylation (e.g. CGI vs
  non-CGI promoters).
- **Analysis:** simulate gene-set-level cis-DMR direction counts under neutral, methylation-gated CpG mutation with each
  set's own CpG density and germline methylation, and see how often the sign test fires.
- **Data:** their hybrid ASM DMRs (bioRxiv 10.64898/2026.01.20.700710) plus primate alignments.

## 4. Haplotype-resolved long-read methylation is the natural data for this model

The model's unit is **per-haplotype methylation as a function of that haplotype's own sequence (cis)**. Phased
long-read 5mC measures exactly that, and it is free of array CpG bias. Sources:
- **HPRC release 2** (232 individuals). HiFi and ONT data carry 5mC calls, and UCSC's HPRC2 hub has per-CpG bigWigs for
  each haplotype assembly. That is haplotype-resolved in practice, though "phased" is not stated explicitly.
- **T2T ape genomes** (Yoo 2025 Nature). Long-read 5mC was used to classify promoters by whether their methylation is
  conserved across species, in fibroblast and lymphoblastoid lines. That gives somatic-cell comparative methylation on
  complete ape assemblies.
- Allele-specific methylation (ASM) pipelines on any phased long-read data.
1. **Fit the map.** Estimate (a₀, a₁, a₂, λ, W) from haplotype-resolved methylation:
   - λ from the width of LMRs around TF footprints;
   - a₁, W from how far CpG-SNPs shift neighbouring CpGs. Ma et al. Fig. 3C reports effects out to ~50 bp, which is an
     independent check.
   The simulations then run on a data-calibrated map instead of a hand-set one.
2. **Allele-specific effect sizes for analysis C,** with frequencies from gnomAD/1000G and mechanism classes from
   sequence (CpG-SNP vs motif SNP).
3. **Germline vs soma.** Somatic long-read methylomes give the functional readout; germline data give the mutational
   one. The two disagree exactly where analysis A predicts an informative contrast.

## 5. What would make this more than a toy

1. Re-run E2 statistics at realistic θ (and with recurrence matched to the human sample size of interest).
2. Add the measurement and ascertainment layer, and show the effect–MAF test behaves under the null (M0 and P).
3. Fit the map to haplotype-resolved data (§4.1).
4. Then run analysis B, which is the cheapest real-data test (public polymorphism + alignments + germline methylomes)
   and needs no new methylation data.
