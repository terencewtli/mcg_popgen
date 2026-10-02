# P0: priors for CpG evolution across primates under different models

Code: `scripts/python/priors.py`, `scripts/slim/p0_prf.slim`. Output: `csv/p0/`, `pdf/p0/priors.pdf`. Run 2026-10-01.

## The problem, stated precisely

At each orthologous position across species there are three observables. First, whether a CpG is present in each
species. Second, its methylation in each species and tissue. Third, CpG-destroying polymorphism within each species.
A faithful population-genetic model for one CpG has two knobs:

- **Mutation:** the loss rate is α·ℓ(m), with ℓ(m) = 2[1 + (r−1)m] + 4. Here m is **germline** methylation, because
  deamination has to happen in the germline to be inherited. r ≈ 10–20.
- **Selection:** S = 2Ns favouring the CpG. Origin-fixation multiplies loss substitutions by φ(−S) and gains by φ(+S),
  where φ(S) = S/(1 − e^−S).

Does origin-fixation suffice, or is the full Wright–Fisher simulation needed? For divergence between species,
origin-fixation suffices. Per-site 4Neμ is about 10⁻³ (about 10⁻² at methylated CpGs), well inside the weak-mutation
regime, and E0 validated Kimura's formula in SLiM. Forward WF simulation is still needed for four things:
polymorphism, linked selection, ancestral polymorphism (a substantial share of human–chimp divergence) and
non-equilibrium demography.

All priors below are in units of neutral non-CpG divergence D, so Ne, μ and generation time drop out.
D values are approximate genome-wide figures: human–chimp ≈ 0.012, human–macaque ≈ 0.065.

## Results

**A/B. CpG o/e cannot tell mutation from selection.** With no selection, fully methylated CpGs equilibrate at o/e =
0.14–0.26 for r = 20–10. That matches the bulk genome (~0.2–0.25), so no selection is needed to explain genome-wide
depletion. A CGI-like o/e of 0.6 is reached by either of two routes:
- germline methylation m ≈ 0.15 with **no selection** (r = 15); or
- full methylation with **S ≈ 1.2** per CpG (0.85–1.46 for r = 10–20).

These lie on one ridge (panel B), and divergence or composition data alone cannot separate them. Cohen et al. 2011's
"minimal selection" conclusion is identified only because germline methylation is supplied as data. Note also that
S ≥ 2 on methylated CpGs overshoots real CGIs (o/e > 1). Any selection on CpGs per se must be weak.

**C. Orthologous CpG conservation (r = 15).**

| | human–chimp | human–macaque |
|---|---|---|
| methylated, S = 0 | 0.87 | 0.48 |
| methylated, S = 1 / 2 / 4 | 0.92 / 0.96 / 0.99 | 0.65 / 0.80 / 0.95 |
| unmethylated, S = 0 | 0.98 | 0.88 |

The human–chimp distance is nearly saturated and carries little information; the macaque distance is where models
separate. The same ridge appears here too: an unmethylated neutral CpG (0.88 at human–macaque) looks like a methylated
CpG with S ≈ 2.

**D. The way out: polymorphism versus divergence.** A change in mutation rate (m) scales polymorphism and divergence
equally, so the P/D ratio for CpG-loss alleles stays at 1. Selection removes alleles before fixation much more
efficiently than it removes them from polymorphism:

| S | polymorphism | divergence | P/D | singleton fraction (n = 20) |
|---|---|---|---|---|
| 0 | 1 | 1 | 1 | 0.28 |
| 1 | 0.87 | 0.58 | 1.5 | 0.31 |
| 2 | 0.76 | 0.31 | 2.4 | 0.35 |
| 4 | 0.58 | 0.075 | 7.8 | 0.43 |

The polymorphism curve (Sawyer–Hartl Poisson random field) was checked against SLiM at 2Ns = 0, −1, −2, −4, −8: all
|z| ≤ 1.3 (neutral case with 16 replicates). So the statistic that separates "unmethylated" from "selected" is a
McDonald–Kreitman-style contrast: CpG-loss P/D, stratified by germline methylation, against a matched neutral class.
The SFS shift is weak in comparison (0.28 → 0.35 at S = 2).

## Discussion points this sets up

1. **Germline and soma are different arrows.** The mutational footprint (E1) is written by *germline* or
   early-embryo methylation, while function reads *somatic* methylation.
   - Passive retention predicts CpG excess only at elements that are unmethylated in the germline (e.g. TFs active in
     germline or early embryo).
   - A soma-only enhancer with CpG excess cannot be passive. That leaves selection or gBGC.
   - This is testable now with sperm/oocyte methylomes (Molaro 2011; Qu 2018) and tissue methylomes.
2. **gBGC looks like S > 0 on CpGs.** B = 4Ne·b is often quoted at ~0.3–1 genome-wide (approximate; higher in hotspots),
   the same scale as the S ≈ 1 needed for CGIs. It acts on every W↔S site, so non-CpG W→S/S→W rates at the same loci
   are the control.
3. **Where causal vs passenger enters.** These priors only ask whether CpGs carry S > 0 once germline m is accounted
   for. In the passenger model, S = 0 at every CpG (selection sits on motif bases). In the causal model, S > 0 at CpGs,
   but under a degenerate map it is spread thinly across many of them. So the order is:
   - P0: is there any S at CpGs?
   - E2: if so, is it shaped like selection on a methylation *level* (degenerate, compensatory turnover, an
     intermediate optimum) or on specific CpGs (e.g. CpG-containing or methyl-sensitive motifs)?
4. **Caveats.**
   - o/e is computed against uniform composition; GC-rich CGIs shift the reference.
   - Sites are treated as independent (no neighbour effects).
   - Very large human samples saturate CpG polymorphism through recurrent mutation, so P must come from modest n or a
     model of recurrence.
   - Demography distorts the SFS, but much of that cancels against a matched neutral CpG class.
