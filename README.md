# mcg_popgen

Small forward population-genetic simulations (SLiM 5.2) of sequence-encoded DNA methylation: when does methylation leave a
footprint of being *causal* for fitness, as opposed to a *passenger* of TF binding? Simulation only; it runs on a laptop.

- `md/CHATGPT.md`, `md/CLAUDE_REVISIONS.md`: original spec and revisions
- `md/PLAN.md`: current compact analysis plan
- `md/ENDPOINT.md`: endpoint checklist for the simulation phase; cluster hand-off and two-account repo rules
- `tex/theory.tex` → `pdf/theory.pdf`: full generative model and validation theory
- Results, in order:
  - `md/E0.md`: simulator validation (18/18 checks pass)
  - `md/E1.md`: passenger-methylation CpG footprint (2.6×, with no selection on CpGs)
  - `md/PRIORS.md`: primate-scale priors; P/D separates germline hypomethylation from selection
  - `md/E2.md`: causal vs passenger under stabilising selection (methylation conserved through full CpG turnover)
  - `md/E2b.md`: E2 at realistic θ; passenger still neutral, causal signals stronger, and turnover with conserved R only at weak per-CpG selection
  - `md/E3.md`: TF + methylation both causal (both signatures; alternative optima across lineages); methylation-sensitive TF indistinguishable from passenger
  - `md/E4.md`: meQTL detection layer; ascertainment-matched CpG-SNP MAF test works, naive effect–MAF correlation does not
- `md/LITERATURE.md`: precedent map
- `md/REAL_DATA.md`: what's toy vs not, scaling, and four real-data analyses the model sets up
- Figures in `pdf/`, summary tables in `csv/`. Raw simulation output (`sim/`) is not mirrored; regenerate it with
  `scripts/python/{e0,e1,priors,e2,e2b,e3}.py run` and `scripts/python/e4.py`.
- `scripts/slim/`: SLiM models; `scripts/sync_to_github.sh` mirrors the working dir here

SLiM: `~/miniconda3/envs/slim4/bin/slim` (5.2, passes `-testSLiM`/`-testEidos`).
