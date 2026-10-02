# Revisions to the project spec: methylation evolution generative model

This document amends the original project spec ("A population-genetic generative model of sequence-encoded DNA methylation evolution"). Read both. **Where they conflict, this document takes precedence.** Where this document is silent, the original spec holds.

The original spec's strengths should be kept: explicit competing models, staged development, exposing assumptions rather than hiding them, and treating identifiability as a result rather than a given. The revisions below sharpen the central question, add a missing null model, fix the computational plan, and narrow the first deliverable.

---

## 1. Revised central question

The original spec lists ten questions. The first version of the project should answer **one**:

> Under a degenerate sequence → methylation map (many sequences produce the same methylation state), with or without methylation-dependent mutation, can comparative sequence + methylation data distinguish **selection acting through methylation** from **selection acting on TF binding, with methylation as a passive readout**?

Questions 3 and 7 of the original spec (methylation conserved while sequence diverges; different sequences converging on the same methylation state) are the scientific core of this question. Treat the other questions as secondary and only pursue them if they fall out naturally.

### Why this framing

If methylation is regenerated from sequence each generation (the mammalian case, given germline reprogramming), then "selection through methylation" is formally selection on a sequence-encoded trait. Models 1 and 2 in the original spec therefore differ **only in the shape of the genotype → phenotype map**: which mutations have fitness effects, and how they interact. The informative signal will come from redundancy, compensatory substitutions, and site turnover under a degenerate map. This mirrors the transcription-factor binding-site turnover literature, which should be used as a template (see Section 6).

---

## 2. Required addition: the passenger-methylation null (Model 4)

The original spec omits the most important confound. Add it as a first-class model, implemented alongside Models 0–3.

### Model 4: Methylation as a passenger of TF binding

```
sequence
   ↓
TF binding (motif occupancy)
   ↓            ↘
expression     methylation   (binding protects CpGs from methylation)
   ↓
fitness
```

Methylation is determined by sequence **and** by TF occupancy, but has **no causal effect** on expression or fitness.

### Why it matters

Much mammalian methylation is a consequence of regulatory activity rather than a cause of it. TF binding protects nearby CpGs from methylation (e.g., low-methylated regions; Stadler et al. 2011), and transcription itself shapes methylation patterns. Model 4 can produce:

- strong methylation conservation across lineages,
- genotype → methylation correlations,
- apparent constraint on methylation-relevant sites,

all without methylation contributing to fitness. This is a harder null to reject than Model 2 (direct sequence selection) and is the one reviewers would raise first.

### Implementation requirements

- TF occupancy should be a function of motif match score (e.g., a logistic function of a position-weight-matrix score), kept as simple as the methylation model.
- Methylation probability at a CpG should decrease with nearby TF occupancy.
- The fitness function must take expression as input, and expression must depend on TF occupancy only, never on methylation.
- Models 1 and 4 should be parameterized so they can be **matched on observable methylation conservation**. The inference question is whether anything else distinguishes them once that is matched.

### Updated model set

| Model | Sequence → methylation | Methylation → fitness | Methylation → mutation | TF binding → fitness |
|---|---|---|---|---|
| 0 (neutral) | yes | no | no | no |
| 1 (methylation-mediated selection) | yes | yes (via expression) | optional | no |
| 2 (direct sequence selection) | yes | no | optional | sequence features act directly |
| 3 (methylation-dependent mutation) | yes | no | yes | no |
| 4 (passenger methylation) | yes, modulated by TF occupancy | no | optional | yes (via expression) |

Methylation-dependent mutation should be a switch that can be combined with any model, not only a standalone model.

---

## 3. Prior work to benchmark against

Before building, read and summarize the following, and state explicitly how this project differs from each:

1. **Cohen, Kenigsberg & Tanay (2011), *Cell*. "Primate CpG islands are maintained by heterogeneous evolutionary regimes involving minimal selection."** This paper modeled CpG island dynamics with methylation-dependent deamination and concluded that little selection is needed to maintain them. It is essentially the original spec's Null B, already tested on primate data. The project should build on it and identify what it could not distinguish, not rediscover it. A useful early sanity check is reproducing its qualitative conclusion under Model 3.
2. **Sella & Hirsh (2005), *PNAS*.** Stationary distribution of genotypes under weak-mutation origin-fixation dynamics. This is the basis for the analytic engine in Section 4.
3. **Berg, Willmann & Lässig (2004), *BMC Evolutionary Biology*; Mustonen & Lässig (2005, 2008).** Evolutionary models of TF binding-site turnover and fitness landscapes over binding energy. These are the closest template for modeling a degenerate sequence → phenotype map and site turnover.

If any citation details are wrong, flag it rather than silently correcting or inventing.

GC-biased gene conversion is a known alternative explanation for CpG/GC patterns. Include it as an optional mutation-process switch, or justify explicitly why it is excluded.

---

## 4. Revised architecture: two engines, not one

The original spec routes everything through SLiM. That is inefficient for the comparative question, which concerns substitutions between species rather than polymorphism within one population.

### Engine A: origin-fixation model (build first)

For a short sequence (start with ~20–200 bp, a handful of CpGs and one motif) under weak mutation (Nμ ≪ 1):

- Represent evolution as a continuous-time Markov chain over sequences.
- Transition rate from sequence *g* to single-mutant neighbor *g′* = (population mutation supply) × (mutation rate for that change, which may depend on methylation state of *g*) × (fixation probability given the fitness difference).
- Use Kimura's diffusion fixation probability. **Derive and document the exact form for the chosen ploidy, dominance, and Ne convention**, rather than copying a formula, and verify it numerically against SLiM (Section 7).
- Under reversible mutation, Sella & Hirsh give the stationary distribution in closed form (stationary probability ∝ neutral stationary probability × fitness raised to a power of order the effective population size; derive the exact exponent for the chosen convention). Methylation-dependent mutation breaks reversibility; solve the chain numerically instead.
- Run the chain along a phylogeny (start with a two-species split, then a five-species primate-like tree) to generate tip sequences and tip methylation states.

This engine makes Questions 3, 4, and 7 approachable semi-analytically and makes large parameter sweeps cheap.

Sequence space grows as 4^L, so either keep L small enough for exact computation or simulate trajectories on the chain (Gillespie) for longer sequences. State which is used.

### Engine B: SLiM forward simulation (build second)

Use SLiM only for questions that need within-population polymorphism:

- site-frequency spectra,
- meQTL-like genotype → methylation correlations,
- linked selection, if recombination is added,
- checking where the weak-mutation assumption of Engine A breaks down.

### Shared components

The sequence → methylation map, TF-occupancy map, expression function, and fitness function must be **one shared Python module** used by both engines, so that any difference between engines reflects population dynamics rather than implementation drift. SLiM will need these reimplemented in Eidos; write a test asserting the Eidos and Python versions agree on a fixed set of sequences.

---

## 5. Computational corrections

### The original benchmark simulates almost nothing

With N = 1,000 diploids, L = 5 kb, and μ ≈ 10⁻⁸ per site per generation, the new-mutation supply is 2NLμ ≈ 0.1 per generation. Over 2,000 generations (0.5N) the population barely leaves its starting sequence and is nowhere near mutation–selection–drift equilibrium.

### Required fixes for Engine B

- **Rescale.** Multiply μ, s, and recombination rate by a factor Q while dividing N and generation count by Q, holding Nμ, Ns, and Nr fixed. Document Q for every run. Verify that summary statistics are approximately invariant to Q over a reasonable range (Section 7).
- **Burn in.** Run at least ~10N generations (in rescaled units) before recording anything, or start from an msprime-generated equilibrium population. Report the burn-in used.
- **Divergence needs a split.** Comparative statistics require at least two populations descended from a common ancestor and run for times comparable to the divergence of interest. A single population cannot produce between-species divergence.

### SLiM implementation snag

SLiM's nucleotide models support sequence-context mutation matrices (e.g., trinucleotide), but a mutation rate that depends on **methylation state**, where methylation depends on CpG density over a window or on TF occupancy, is not expressible as a fixed context matrix. It will need a `mutation()` callback or a custom approach. Prototype this on a tiny case and measure runtime before committing to it. If it is prohibitively slow, report that and propose alternatives rather than quietly simplifying the model.

---

## 6. Modeling assumptions to state explicitly

In addition to the list in Section 11 of the original spec, document and justify:

1. **How degenerate the sequence → methylation map is.** For example, a CpG-density threshold model is highly degenerate; a per-CpG model is not. This choice largely determines whether Questions 3 and 7 have interesting answers, so treat degeneracy as a swept parameter, not a fixed choice.
2. **The direction of causality between methylation and TF binding** in each model.
3. **Whether fitness is defined on expression only**, and what selection strength (Ns) each parameter setting corresponds to. Report Ns for every run; raw s alone is not interpretable.
4. **Ploidy, dominance, and Ne conventions** used in fixation probabilities, consistent across both engines.
5. **What "methylation conservation" means operationally** (e.g., correlation of mean methylation across orthologous CpGs, or agreement in binary methylated/unmethylated calls) and how alignment of non-conserved CpGs is handled.

---

## 7. Validation tests

### Engine A

- Neutral case: stationary distribution matches the mutation-only expectation.
- Reversible mutation with selection: numerical stationary distribution matches the Sella & Hirsh closed form.
- Fixation probabilities match Kimura's formula for single-site cases.

### Engine B

- Neutral case: expected nucleotide diversity ≈ 4Nμ (in the infinite-sites limit), neutral SFS shape, neutral fixation probability ≈ 1/(2N).
- Selected single site: fixation probability matches Kimura.
- Rescaling invariance: key summary statistics stable across a range of Q.

### Cross-engine

- In the weak-mutation regime, substitution rates and tip distributions from Engine B agree with Engine A within sampling error. Report where they diverge as Nμ increases.

### Shared module

- Eidos and Python implementations of the methylation, occupancy, expression, and fitness functions agree on a fixed test set.

---

## 8. Revised first deliverable

Replace Section 14 of the original spec with the following, in order:

1. A written summary (one page or less) of the three prior-work items in Section 3 and how this project differs.
2. A written statement of the exact generative assumptions for Models 0–4, including the answers to Section 6 above, **before writing simulation code**. Stop and ask for review at this point.
3. The shared Python module for sequence → methylation, TF occupancy, expression, and fitness, with unit tests.
4. Engine A on a two-species split, with its validation tests passing.
5. A first identifiability experiment using Engine A only: simulate Models 1 and 4 matched on methylation conservation, and test whether any candidate statistic separates them. A null result is a valid result; report it as such.
6. Only after the above: Engine B, its validation tests, and a runtime/memory benchmark on the local machine.

---

## 9. Explicitly deferred

Do not build these in the first version, but keep the architecture compatible with them:

- **Trans effects.** Everything above is cis-only. A natural later extension is a trans factor that co-evolves with sequence, such as a KRAB-zinc-finger protein targeting transposable elements for methylation (cf. Jacobs et al. 2014). This is where cis/trans decomposition from interspecies hybrid data could eventually connect to the model.
- Real comparative methylome data and all the issues listed in Section 6 of the original spec.
- Recombination, unless a linked-selection question requires it.
- Realistic (e.g., neural-network) sequence → methylation models.

---

## 10. Working principles

- Prefer the smallest model that can answer the central question. Add complexity only when a specific question requires it, and say which question.
- Never tune a model or statistic until it "works." If a statistic is uninformative, that is a result.
- When an implementation constraint forces a modeling compromise, stop and report it rather than making the change silently.
