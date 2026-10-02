# Compact analysis plan: causal vs passenger methylation in a 200-bp SLiM locus

Written 2026-10-01 after reading CHATGPT.md, then CLAUDE_REVISIONS.md. Scope: simulations only, MacBook only, one locus.

## 0. The framing in one paragraph

Mammalian germline reprogramming means methylation is rebuilt from sequence every generation. Under that assumption, methylation
is a **molecular phenotype of the haplotype**, not a second inheritance channel. "Selection on methylation" is then selection on
sequence, and *causal* vs *passenger* differ only in **which mutations carry fitness effects**. A mutation tells the two apart
only if it changes methylation without changing TF occupancy (or the reverse). So the question that can actually be answered is:

> For which mutation classes, and in which parameter regime (map degeneracy × N·s × CpG hypermutability), do causal and passenger
> methylation leave different footprints in (a) substitution patterns and (b) the allele frequencies of meQTL-like variants?

Most of this can be predicted on paper. The simulations check the prediction and show how large the effects are. The non-obvious
parts are the two places where intuition tends to fail:
1. **Methylation-gated mutation** (methylated CpGs deaminate, unmethylated ones don't) lets a passenger model build a CpG-rich,
   unmethylated island around a bound motif with **zero selection on CpGs or methylation**. That footprint looks like selection.
2. **A degenerate map** (methylation set by CpG density over a window) spreads selection on methylation across many sites, each
   with a tiny effect. Every per-site test then reads "neutral" even though the aggregate trait is strongly conserved and CpGs
   turn over in a compensatory way.

## 1. Generative model (exact, shared by every experiment)

**Sequence.** Haplotype x ∈ {A,C,G,T}^L, L = 200, no recombination. One 8-bp motif window W in the centre.

**TF occupancy** (per haplotype): O(x) = σ(β·(S(x) − S₀)), where S is the PWM match score (start with a match count, 0–8).

**Methylation** (per CpG i, per haplotype, deterministic expectation over cells):
m_i(x) = σ(a₀ − a₁·ρ_i(x) − a₂·O(x)·K(d_i))
- ρ_i = CpG count within ±w bp. The **degeneracy knob** is a₁ and w: a₁ = 0 gives per-CpG methylation (not degenerate), while a
  large a₁ with w ≈ 50 gives a region-level CGI-like switch (very degenerate).
- K(d_i) = protection kernel by distance to the motif (Stadler 2011 LMR logic). Setting a₂ = 0 removes TF → methylation coupling.
- Not modelled at first: stochastic epialleles, transgenerational inheritance, methylation → TF binding (methyl-sensitive TFs such
  as NRF1, Domcke 2015). Methylation → TF binding is the one deferral that matters, because it makes "passenger" impure. It is
  queued as model C′.

**Expression / fitness.** Diploid, codominant: E = ½[g(h₁) + g(h₂)], w = exp(−(E − θ)² / 2V_s).
- **P (passenger):** g = O. Methylation has no causal effect.
- **C (causal):** g = 1 − m̄_W, the mean methylation of CpGs within the promoter window. The TF matters *only* through methylation.
- **0 (neutral):** w ≡ 1.

P and C share **identical** occupancy and methylation maps and differ only in which molecular state the fitness function reads.
That is the cleanest causal/passenger contrast available. CLAUDE_REVISIONS.md's Model 2 (direct sequence selection) is just P with a₂ = 0,
so it is merged into P.

**Mutation.** K80 background rate μ, with CpG transitions multiplied by r. A switch selects between:
- `ctx`: r constant (context-only hypermutability, as in a trinucleotide matrix)
- `meth`: r·m_i(x) (methylation-gated, evaluated on the mutating haplotype)
- `off`: r = 1

gBGC is excluded in v1 because there is no recombination and therefore no gene conversion. That exclusion is explicit and
deliberate.

**Reported scale for every run:** 4Nμ per bp, r, and 2N·s_eff for each single-mutation class, computed from Δw of the
heterozygote at the optimum genotype. Raw s is never reported alone.

**SLiM implementation** (prototype runs, `scripts/slim/smoke_m4_methmut.slim`): a nucleotide model with a trinucleotide matrix
sets CpG transitions at the maximum rate r·μ. A `mutation()` callback accepts each proposed CpG transition with probability
m_i(haplotype), which gives the `meth` mode exactly. A `fitnessEffect()` callback computes E from both haplosomes. N = 1000,
L = 200, 10,000 generations: **12.7 s, 11 MB, one core**. With 10 cores that is about 45 replicates per minute.

## 2. Experiments (each is under ~1 h on the laptop)

**E0: validation, no methylation (do first).** Check neutral π against 4Nμ (finite-sites JC), the SFS against 1/i, single-site
fixation probability against Kimura under SLiM's 1+hs convention, and invariance of all of these to the rescaling factor Q (N/Q,
μ·Q, generations/Q). Then add `ctx` CpG hypermutation and compare equilibrium CpG o/e with a ~50-line Python origin-fixation
Gillespie on the same matrix. That is the only part of CLAUDE_REVISIONS.md's "Engine A" I would build, as an oracle rather than a second
engine.

**E1: the passenger footprint, no selection on methylation.** Compare 0+`ctx`, 0+`meth`, P+`ctx` and P+`meth`. Readouts are CpG
o/e and m̄ as functions of distance to the motif, at stationarity (10N burn-in).
Prediction: only P+`meth` produces a local CpG-rich, hypomethylated island. Its size should scale with the protection kernel and
with r. The quantity of interest is **how much "CpG conservation" is free**: the motif-proximal CpG substitution rate relative to
flanks, with s = 0 on every CpG. This is the mechanism behind Cohen, Kenigsberg & Tanay 2011, rebuilt as a forward model.

**Status (2026-10-01).** E0 is done (`md/E0.md`): all checks pass. E1 is done (`md/E1.md`): a 2.6× CpG footprint
from gating alone. P0 is done (`md/PRIORS.md`): o/e and divergence alone cannot separate germline hypomethylation from weak selection (S ≈ 1); the CpG-loss P/D ratio can. E2 is next.

**E2: causal vs passenger, matched (the core experiment).** *Revised 2026-10-01: use an intermediate optimum*
*(e.g. E\* = 0.5 for promoter methylation, and the equivalent for occupancy). That makes selection truly stabilising,*
*so many genotypes sit at the optimum and compensatory turnover can happen. An optimum at the trait's maximum, as in*
*E1, is directional with diminishing returns.* Turn `meth` on and tune P and C on a small grid until mean methylation
and the CpG o/e profile match. Then test which statistics still separate them:
- **Constraint map:** substitution rate per site class (motif / motif-proximal CpG / distal CpG / non-CpG proximal), divided by
  the flank rate. Prediction: distal-but-in-window CpGs are constrained in C and free in P, while motif sites are constrained in
  both.
- **meQTL × allele-frequency test:** for each segregating variant, take Δm (its methylation effect) and ΔO, then regress derived
  allele frequency on |Δm| for variants with ΔO ≈ 0. Prediction: the slope is ~0 in P and negative in C. This is the statistic
  that carries over to real data (meQTL effect size vs DAF/SFS). The simulation shows the confounder it must survive: CpG-loss
  alleles are hypermutable, so recurrent mutation pushes their frequency up regardless of selection.
- **Degeneracy sweep:** run (a₁, w) × N·s and plot where each statistic stops separating C from P. Analytic expectation: per-site
  identifiability needs 4N·s·|Δm_i| ≳ 1, and near saturation of a degenerate switch |Δm_i| → 0. **This phase diagram is the main
  figure.**

**E3 (only if E2 is interesting): two-population split.** Run 4N–20N generations after the split and plot methylation divergence
against sequence divergence for 0, P and C (CHATGPT Q3/Q7: is methylation conserved while sequence diverges?).

**Rescaling for all runs:** N = 1000, 4Nμ ≈ 0.004 per bp (μ = 1e-6), r ∈ {1, 10, 30}. The smoke test's μ = 2.5e-5 is for timing
only. Q-invariance is checked in E0 before anything is interpreted.

## 3. What I'd change from CLAUDE_REVISIONS.md

- Keep its core points: regenerate methylation from sequence, include the passenger null, report N·s, rescale and burn in, and
  treat null results as results.
- Drop the two-engine architecture, the five-species phylogeny and the Eidos↔Python parity suite for now. They are too heavy for
  "think small". A tiny Python oracle for E0 is enough. In practice, parity can be tested by computing the maps on a fixed sequence
  in both.
- Drop the "SLiM can't do methylation-dependent mutation" worry. The rejection-sampling `mutation()` callback is exact and cheap,
  because it runs only on proposed mutations.
- Reframe C vs P around **fitness reading different molecular states from shared maps**, rather than around different maps.

## 4. Honest assessment

It is a good learning vehicle and a weak paper on its own. The headline answer ("distinguishable iff decoupling mutations exist
with 4Ns|Δm| > 1") falls out of the setup, and Cohen 2011 plus the TFBS-turnover literature (Berg/Willmann/Lässig 2004;
Mustonen & Lässig 2005, 2008) cover adjacent ground. The value is in two places. E1 quantifies how much methylation-associated
sequence conservation is free under passive retention. E2 gives a meQTL × DAF test with known power and known confounders, which
connects directly to real QTL data later.
