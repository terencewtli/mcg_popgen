# mcg_popgen

Small forward population-genetic simulations (SLiM 5.2) of sequence-encoded DNA methylation: when does methylation leave a
footprint of being *causal* for fitness, as opposed to a *passenger* of TF binding? Simulation only; it runs on a laptop.

- `md/CHATGPT.md`, `md/CLAUDE_REVISIONS.md`: original spec and revisions
- `md/PLAN.md`: current compact analysis plan
- `tex/theory.tex` → `pdf/theory.pdf`: full generative model and validation theory
- `scripts/slim/`: SLiM models; `scripts/sync_to_github.sh` mirrors the working dir here

SLiM: `~/miniconda3/envs/slim4/bin/slim` (5.2, passes `-testSLiM`/`-testEidos`).
