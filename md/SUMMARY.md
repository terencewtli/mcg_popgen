# Summary: what the simulations say about causal vs passenger methylation

Simulation phase closed 2026-10-01. Figure-led version: `pdf/summary.pdf`. Model: `pdf/theory.pdf`. Checklist and
cluster hand-off: `md/ENDPOINT.md`. Real-data designs A–D: `md/REAL_DATA.md`.

## Question

> Which population-genetic signatures distinguish **causal** methylation (fitness reads methylation) from
> **passenger** methylation (methylation is set by TF binding, and fitness reads binding)? At what selection strength,
> and which survive realistic θ and meQTL detection?

## Model in one paragraph

- **Population and locus.** A Wright–Fisher population (SLiM 5.2) carries a 300-bp locus with one TF motif. Occupancy
  O and per-CpG methylation m are deterministic functions of each haplotype's own sequence.
  - Methylation is regenerated every generation (germline reprogramming).
  - A bound TF protects nearby CpGs, giving low-methylated regions.
- **Mutation.** Methylated CpGs mutate up to r = 20× faster (deamination).
- **Fitness.** Fitness reads O (passenger), methylation load R = Σm (causal), both, or O·(1 − m_motifCpG)
  (methylation-sensitive TF). Selection is stabilising at an intermediate optimum.
- **Scaling.** Everything is reported as 2Ns and θ. Rescaled N = 100 is exact for these quantities (checked), and the
  headline results are at human-like θ = 0.001.

**Validation.** All of these were checked:
- the simulator against analytic theory and an exact single-lineage oracle (E0: 18/18 checks);
- Python–SLiM parity for every map (≤ 5×10⁻⁷);
- rescaling invariance with selection (E1);
- the Sawyer–Hartl polymorphism density (P0).

## Findings

| # | finding | evidence |
|---|---|---|
| F1 | **Passive retention.** With no selection on any CpG, a selected TF leaves a 2.6× CpG enrichment within 25 bp (3.6× within 10 bp), but only when CpG mutation is methylation-gated. The footprint scales with protection width and r. | E1 |
| F2 | **Composition cannot separate mutation from selection.** CGI-like o/e 0.6 arises from germline m ≈ 0.15 with no selection, or full methylation with per-CpG S ≈ 1.2. Divergence data lie on the same ridge. Genome-wide CpG depletion needs no selection (o/e 0.14–0.26 at r = 20–10). | P0 |
| F3 | **CpG-loss P/D breaks the ridge.** A mutation-rate change leaves P/D = 1; selection raises it. At human-like θ: 2.1 for 2Ns ≈ 1.6 (P0 predicted 2.0), 8–14 for 2Ns ≈ 8–11. P0's directional formula is an upper bound under stabilising selection. | P0, E2b, E3 |
| F4 | **Passenger methylation is invisible to frequency statistics.** CpG-only variants under P, P (CRE) and the methylation-sensitive TF have P/D 0.98–1.09, DAF 0.97–1.00 of neutral, and ascertainment-matched meQTL MAF 1.00 ± 0.05, at both θ values and all sample sizes. This is the most robust result. | E2, E2b, E3, E4 |
| F5 | **Causal methylation depletes CpG-only variants, alone or alongside a selected TF.** DAF 27–78% of neutral; detected-meQTL MAF 0.45–0.90 of matched neutral at n = 30,000 when per-CpG 2Ns ≳ 2, with power rising with n. The signal sits mostly below MAF 0.01. | E2b, E3, E4 |
| F6 | **Motif conservation and motif-SNP meQTL depletion mark TF selection, not causal methylation.** Motif-SNP meQTL MAF is 0.04–0.17 of neutral in every model with a selected TF. | E3, E4 |
| F7 | **The naive effect–frequency correlation is biased.** Among detected meQTLs, corr(\|β\|, MAF) = −0.09 to −0.13 under passenger (detection bias), as negative as under causal. | E4 |
| F8 | **Methylation can be conserved while CpG positions turn over, but only under weak per-CpG selection** (2Ns ≈ 1.6). At 2Ns ≈ 8–11 and human-like θ, positions mostly freeze. The full turnover seen at high θ in E2 needed co-segregating compensatory mutations. | E2, E2b |
| F9 | **Alternative TF/methylation optima.** When both TF and methylation matter, lineages commit to different (motif strength, methylation load) combinations with the same expression and do not switch. The proportions depend on the burn-in; that two such states exist does not. | E3 |
| F10 | **Methylation sensitivity of a bound TF is undetectable.** A strongly bound TF keeps its own CpG unmethylated, so MS ≈ passenger control on every statistic. | E3 |

**Corrected along the way:** E2's "full CpG turnover at κ = 10" (a high-θ effect; F8), and E2's first matching
(on expression, not R; Jensen's inequality). Both are documented in place.

## Predictions for real data

| # | prediction | test | real-data design |
|---|---|---|---|
| R1 | CpG excess and cross-species CpG conservation at regulatory elements track **germline** hypomethylation. Soma-only hypomethylated elements show no passive excess; any excess there implies selection or gBGC. | F1, F2 | A |
| R2 | Stratified by germline methylation, CpG-loss **P/D = 1 where only the mutation rate differs**, and P/D > 1 where CpGs are selected. Non-CpG W↔S sites at the same loci control for gBGC. | F2, F3 | B |
| R3 | **CpG-SNP meQTLs outside motifs:** detected-MAF ratio = 1 against an ascertainment-matched neutral under passenger methylation, and < 1 (decreasing with n) under causal. Motif-SNP meQTLs are depleted either way. | F4–F6 | C |
| R4 | Pooled meQTL effect–MAF correlations without a matched null **do not** indicate selection. | F7 | C |
| R5 | Elements under weak per-CpG constraint show CpG positional turnover with conserved methylation across species; strongly constrained elements conserve CpG positions. | F8 | A (comparative) |
| R6 | Orthologous elements with matched expression can differ in both motif strength and methylation, so **methylation divergence is not evidence of expression divergence.** This is relevant to hybrid cis-DMR interpretation. | F9 | D |
| R7 | Mutation-driven direction bias in CpG gains vs losses differs between gene sets by CpG content and germline methylation. The Ma et al. sign test needs that null. *Not yet simulated.* | F1–F3 | D |

## Limits

- **Toy biology:** one locus, a hand-set map, cis only, no trans factors, and no feedback between methylation and
  binding.
- **Frequency resolution:** at N = 100 the minimum MAF is 0.005, so the rare tail where selection acts is barely
  represented. Power in E4 is probably underestimated.
- **E4 detection model:** analytic power, no LD, and a regional phenotype rather than single-CpG arrays.
- **Composition starts** come from high-θ runs (E2b/E3); no composition drift was detected.
- **Open item (E1):** at r = 50, the footprint is about 2 SE below the oracle.

## Hand-off: first cluster jobs

1. **Real-N check** (ENDPOINT item 8): C κ=10 and P κ=10 at N = 10⁴ with human-scale rates, no rescaling. Expect E2b's
   P/D and DAF to hold, and E4's power to rise because the rare tail is resolved. This is a SLiM job array, with SLiM
   from conda.
2. **Real-data analysis B** (CpG-loss P/D; public data only), or **map fitting on HPRC2 haplotype-resolved
   methylation** (`md/REAL_DATA.md` §4.1), which also supplies effect sizes for analysis C.
3. Follow the repo rules in `md/ENDPOINT.md`: the cluster owns `realdata/`, works directly in a clone, and pulls
   before committing.
