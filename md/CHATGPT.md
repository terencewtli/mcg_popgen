# Project: A population-genetic generative model of sequence-encoded DNA methylation evolution

## Goal

Build a small, interpretable population-genetic simulation to study the evolutionary dynamics of:

**DNA sequence → methylation → regulatory phenotype → fitness → sequence evolution**

The initial goal is NOT to model real genomes or reproduce real methylomes. The goal is to establish a minimal generative model and determine what evolutionary signatures arise when methylation is:

1. purely a sequence-dependent molecular state,
2. phenotypically relevant but not itself directly inherited,
3. subject to selection through its effect on phenotype/fitness.

Eventually, the model could be connected to comparative methylome data across species, but the first stage should be entirely theoretical/simulation-based.

---

# 1. Core generative model

Consider a short regulatory sequence (initially ~1–10 kb) in a diploid population.

Each individual has:

### A. DNA sequence

A nucleotide sequence containing:

* CpGs
* potentially methylation-sensitive CpGs
* optional sequence motifs
* optional regulatory features

Sequence evolves through:

* point mutation
* optionally recombination
* optionally methylation-dependent mutation

Start with point mutation only unless recombination becomes necessary.

---

### B. Methylation state

Methylation at each relevant CpG/site is generated probabilistically from local sequence context.

Conceptually:

```
sequence → P(methylation)
```

Possible determinants:

* CpG presence/absence
* local CpG density
* short sequence motifs
* GC content
* optional TF-like motifs

Initially use a deliberately simple function rather than a realistic neural network.

For example:

```
logit(P(M_i = 1)) =
    intercept
    + β_CpG * local_CpG_density
    + β_motif * motif_score
    + noise
```

The exact parameterization should be modular so that alternative sequence→methylation models can be tested.

---

### C. Regulatory phenotype

Methylation is translated into a phenotype such as gene expression.

For example:

```
expression = f(methylation)
```

Start with a simple monotonic relationship.

Later consider:

* threshold effects
* nonlinear effects
* multiple CpGs
* methylation at enhancers vs promoters
* opposing methylation effects at different sites

---

### D. Fitness

Fitness depends on the regulatory phenotype.

The initial model should include at least:

### Stabilizing selection

```
w = exp[-s * (expression - optimum)^2]
```

This is useful because it creates a biologically interpretable setting in which multiple molecular configurations could potentially produce similar fitness.

Eventually consider:

* directional selection
* stabilizing selection
* weak selection
* no selection

---

### E. Population reproduction

Use a Wright–Fisher-style diploid population.

Each generation:

1. calculate methylation
2. calculate phenotype
3. calculate fitness
4. sample parents according to fitness
5. generate offspring
6. introduce mutations
7. repeat

Initially use a single population with constant population size.

---

# 2. Crucial biological distinction

The model must distinguish between:

## Model 0: Neutral

Sequence evolves neutrally.

Methylation is a sequence-dependent molecular consequence but has no effect on fitness.

```
sequence → methylation
```

but NOT:

```
methylation → fitness
```

---

## Model 1: Sequence-mediated selection

Sequence affects methylation.

Methylation affects expression.

Expression affects fitness.

```
sequence
   ↓
methylation
   ↓
expression
   ↓
 fitness
```

Selection therefore acts indirectly on sequence through its methylation consequences.

---

## Model 2: Direct sequence selection

Some sequence features affect fitness directly, independently of methylation.

```
sequence ─────────→ fitness
   ↓
methylation
   ↓
expression
   ↓
 fitness
```

This is an important confounding/null model.

---

## Model 3: Methylation-dependent mutation

Methylation changes the mutation process.

For example:

```
methylated CpG → elevated C→T mutation probability
```

This creates a feedback loop:

```
sequence
   ↓
methylation
   ↓
mutation rate
   ↓
sequence evolution
```

This model is especially important because methylation-associated sequence evolution does NOT necessarily imply selection.

---

# 3. Essential scientific questions

The simulation should be designed around these questions rather than around producing a particular result.

### Question 1

If methylation is completely neutral, what sequence/methylation evolutionary patterns arise purely from:

* mutation
* drift
* methylation-dependent mutation?

### Question 2

If methylation affects fitness, what signatures distinguish this from neutral sequence evolution?

### Question 3

Can selection maintain a methylation state without maintaining the exact underlying DNA sequence?

In other words:

> How conserved can methylation be while sequence diverges?

### Question 4

Conversely, can sequence be highly conserved while methylation evolves?

### Question 5

When does selection preserve:

* individual CpGs?
* CpG density?
* transcription-factor motifs?
* broader sequence architecture?
* methylation itself?

### Question 6

How much sequence divergence can occur before methylation becomes evolutionarily unstable?

### Question 7

Can different DNA sequences converge on similar methylation states?

This is especially important because it would imply that comparative methylation conservation does not necessarily imply sequence conservation.

### Question 8

Can present-day sequence + methylation data distinguish:

* neutral sequence evolution
* selection acting directly on sequence
* selection acting through methylation?

### Question 9

How does the answer depend on:

* selection strength
* population size
* mutation rate
* methylation-dependent mutation
* strength of sequence→methylation coupling
* number of methylation sites
* phenotype optimum?

### Question 10

What statistics calculated from extant populations/species are actually informative about the underlying evolutionary mechanism?

This last question should be treated as an inference/identifiability problem, not assumed in advance.

---

# 4. Key conceptual experiment

A central experiment should compare:

### Null A

Neutral sequence + sequence-dependent methylation

### Null B

Neutral sequence + methylation-dependent mutation

### Alternative A

Sequence → methylation → phenotype → fitness

### Alternative B

Direct sequence → fitness + sequence → methylation

Then ask whether the resulting present-day data can distinguish these models.

The simulator should generate both:

* the full evolutionary history
* present-day observable data

This lets us deliberately hide the true generating model and test inference.

---

# 5. What to measure

At minimum record:

### Sequence statistics

* allele frequencies
* nucleotide diversity
* sequence divergence
* CpG density
* motif frequencies
* site-frequency spectrum
* haplotype diversity

### Methylation statistics

* mean methylation
* methylation variance
* methylation divergence between populations
* methylation divergence between homologous loci
* methylation conservation despite sequence divergence

### Joint sequence/methylation statistics

* correlation between genotype and methylation
* correlation between local sequence distance and methylation distance
* CpG genotype → methylation effect
* motif genotype → methylation effect
* sequence distance required to produce methylation divergence

### Evolutionary statistics

* fixation probability
* allele-frequency trajectories
* lineage-specific substitutions
* methylation-state transitions
* fitness trajectories

Do not assume any particular statistic is optimal. The purpose of the simulation is partly to discover which observable statistics retain information about the generating process.

---

# 6. Comparative-evolution extension

Eventually the model could be used to ask:

> Given sequence and methylation data from multiple species, which evolutionary mechanisms could plausibly generate the observed joint patterns?

Potential empirical data:

* human
* chimpanzee
* gorilla
* orangutan
* macaque
* other primates
* eventually broader mammals

But DO NOT build the first version around real comparative data.

First establish whether the theoretical model produces identifiable signatures.

Real comparative methylome data introduces additional issues:

* tissue/cell-type mismatch
* developmental stage
* technical batch effects
* methylome measurement error
* alignment uncertainty
* incomplete orthology
* demographic history
* ancestral methylation uncertainty

These should be added only after the core model works.

---

# 7. Software architecture

Prefer established population-genetic software rather than implementing population genetics from scratch.

Recommended stack:

### SLiM

Use SLiM for forward-time population-genetic simulation.

Responsibilities:

* population
* reproduction
* mutation
* selection
* allele-frequency evolution
* optional recombination
* tree-sequence recording

### msprime

Potentially use msprime to generate an initial population with realistic neutral genetic variation before beginning the forward simulation.

Do NOT require this for the first prototype.

### tskit

Use tree sequences for compact ancestry/history storage.

Do not save every individual's genome at every generation.

### Python

Use Python for:

* parameter management
* simulation orchestration
* summary-statistic calculation
* parameter sweeps
* plotting
* inference/identifiability experiments

### R

Optional. Do not introduce R unless it is clearly useful.

---

# 8. MacBook-first computational constraints

The first implementation MUST be designed to run comfortably on a MacBook Pro.

Start with approximately:

```
population size: 500–2,000 diploid individuals
sequence length: 1–10 kb
generations: 1,000–10,000
replicates: 50–500
```

The first benchmark should be approximately:

```
N = 1,000
sequence = 5 kb
generations = 2,000
replicates = 100
```

Do not optimize prematurely. First measure runtime and memory.

The initial goal is to have simulations that run in minutes-to-hours, not days.

---

# 9. How to keep the model computationally small

Use the following simplifications initially:

### 1. Tiny genome

Use 1–10 kb rather than chromosome-scale sequence.

### 2. Small population

N = 500–2,000 initially.

### 3. Short evolutionary history

1,000–10,000 generations initially.

### 4. Few functional elements

Start with perhaps 5–20 methylation-relevant sites rather than millions of CpGs.

### 5. No recombination initially

Add recombination only if linkage/haplotype evolution becomes scientifically important.

### 6. No realistic methylome model initially

Use a simple interpretable sequence→methylation function.

### 7. No real comparative data initially

Simulated data only.

### 8. Development mode

During model development:

```
simulation → summary statistics → discard simulation
```

Do not save every run.

### 9. Production mode

Once the model is stable:

```
simulation → tree sequence → downstream analysis
```

Use tree-sequence recording rather than storing every generation.

### 10. Parallelize independent replicates

Parameter combinations and replicates are embarrassingly parallel.

Do not build distributed infrastructure initially. Python multiprocessing/job arrays should be sufficient.

---

# 10. Development strategy

Implement this in stages.

## Stage 1 — population-genetic sanity check

Build a minimal Wright–Fisher/SLiM simulation with:

* mutation
* drift
* selection
* no methylation

Validate against known population-genetic expectations.

The point is to ensure the simulator is correct before adding biological complexity.

---

## Stage 2 — sequence → methylation

Add a simple sequence-dependent methylation model.

No selection through methylation yet.

Demonstrate that changing sequence parameters produces predictable methylation changes.

---

## Stage 3 — methylation → phenotype → fitness

Add expression/phenotype and fitness.

Compare neutral and selected models.

---

## Stage 4 — methylation-dependent mutation

Add elevated mutation probability at methylated CpGs.

Quantify how much apparent methylation-associated sequence conservation/divergence can arise without selection.

---

## Stage 5 — inference experiment

Generate data under one model.

Pretend the generating model is unknown.

Ask whether summary statistics from present-day sequence + methylation can identify the generating mechanism.

This is likely the most scientifically important stage.

---

# 11. Important modeling questions to expose rather than silently decide

Before implementing the full model, explicitly list and discuss:

1. Is methylation inherited directly, or regenerated from sequence in each generation?
2. Does methylation affect fitness directly or through expression?
3. Is methylation deterministic or stochastic conditional on sequence?
4. Is methylation itself heritable?
5. Does methylation affect mutation rate?
6. Does methylation affect recombination?
7. Does sequence affect methylation locally or through long-range features?
8. Is fitness stabilizing, directional, or context-dependent?
9. Are there multiple equivalent methylation states with similar fitness?
10. How many methylation sites should jointly determine phenotype?
11. Is the fitness effect additive across sites?
12. Should there be epistasis?
13. How should methylation noise affect fitness?
14. What aspects of the model are identifiable from extant data?

Do not hide these assumptions in implementation.

---

# 12. Desired outputs

The project should eventually produce:

### Figure 1

Generative model:

```
DNA sequence
     ↓
methylation
     ↓
 expression
     ↓
  fitness
     ↓
reproduction
     ↓
sequence evolution
```

with the methylation-dependent mutation feedback loop shown separately.

### Figure 2

Example evolutionary trajectories under neutral vs selected models.

### Figure 3

Sequence divergence vs methylation divergence.

### Figure 4

Methylation conservation despite sequence divergence.

### Figure 5

Ability of observable statistics to distinguish the different generating models.

### Figure 6

Parameter sweep showing where evolutionary signatures are identifiable.

These are targets, not requirements. If the simulations reveal that a proposed statistic is uninformative, that should be treated as a scientific result rather than engineered away.

---

# 13. Important scientific principle

Do not begin by assuming that selection is the explanation.

The central challenge is to distinguish:

```
mutation
drift
methylation-dependent mutation
demographic history
direct sequence selection
selection mediated through methylation
```

The model should therefore be capable of producing convincing methylation/sequence evolutionary patterns **without selection**.

Only then should we ask whether selection generates additional, identifiable signatures.

---

# 14. First implementation request

Start by creating a minimal, well-tested prototype rather than the complete system.

Deliver:

1. a SLiM model for a small diploid population,
2. a simple sequence-dependent methylation model,
3. optional methylation-dependent mutation,
4. a simple phenotype/fitness model,
5. Python scripts to run replicate simulations,
6. summary statistics,
7. basic plots,
8. runtime/memory benchmark on the local machine,
9. tests validating the neutral model against expected population-genetic behavior.

Before adding complexity, explain:

* the exact generative assumptions,
* which parameters are biologically meaningful,
* which are merely simulation conveniences,
* what is expected to be identifiable,
* and what alternative model specifications would substantially change the conclusions.

Favor **small, interpretable, falsifiable models over biological realism** in the first version.

