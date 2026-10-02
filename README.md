# mcg_popgen

Small forward population-genetic simulations (SLiM 5.2) of sequence-encoded DNA methylation: when does methylation leave a
footprint of being *causal* for fitness, as opposed to a *passenger* of TF binding? Simulation only; it runs on a laptop.

- `md/CHATGPT.md`, `md/CLAUDE_REVISIONS.md`: original spec and revisions
- `md/PLAN.md`: current compact analysis plan
- `tex/theory.tex` → `pdf/theory.pdf`: full generative model and validation theory
- Results, in order:
  - `md/E0.md`: simulator validation (18/18 checks pass)
  - `md/E1.md`: passenger-methylation CpG footprint (2.6×, with no selection on CpGs)
  - `md/PRIORS.md`: primate-scale priors; P/D separates germline hypomethylation from selection
  - `md/E2.md`: causal vs passenger under stabilising selection (methylation conserved through full CpG turnover)
- `md/LITERATURE.md`: precedent map
- `md/REAL_DATA.md`: what's toy vs not, scaling, and four real-data analyses the model sets up
- Figures in `pdf/`, summary tables in `csv/`. Raw simulation output (`sim/`) is not mirrored; regenerate it with
  `scripts/python/{e0,e1,priors,e2}.py run`.
- `scripts/slim/`: SLiM models; `scripts/sync_to_github.sh` mirrors the working dir here

SLiM: `~/miniconda3/envs/slim4/bin/slim` (5.2, passes `-testSLiM`/`-testEidos`).
