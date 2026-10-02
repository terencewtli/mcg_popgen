# E2b: E2 at realistic θ

Code: `scripts/slim/e2b_lowtheta.slim`, `scripts/python/{e2b,e2b_plot}.py`. Results: `csv/e2b/`, `pdf/e2b/e2b_theta.pdf`.
Run 2026-10-01: 240 runs, 16 min on 9 cores.

## Design

- **Same models, maps, κ and optima as E2;** only α changes. θ = 4N·3α = 0.001 per site, about human, versus 0.024 in
  E2. CpG recurrence 4N·rα ≈ 0.007, close to the human ~0.005, versus 0.16 in E2.
- **Starting point.** Composition equilibrates on a 1/u ≈ 400k-generation timescale, so each run starts from the
  matching E2 run's final consensus. That is near that model's equilibrium.
- **Run length.** 20N burn-in for polymorphism, then 0.5/u = 200k generations of measurement. 40 replicates.
- **Variants are characterised inside SLiM,** so the derived allele is known and frequencies are DAF, not MAF. Each
  variant's effect is averaged over 10 non-carrier backgrounds.
- **No composition drift:** mean R over the second half minus the first half is within ±1 SE for every condition
  except P κ=2, which is −0.21 ± 0.18.

## Results

| | M0 | P κ2 | P κ10 | C κ2 | C κ10 R-matched | C κ10 |
|---|---|---|---|---|---|---|
| mean R | 2.72 | 2.67 | 2.77 | 2.77 | 2.82 | 1.91 |
| P/D, loss of methylated CpGs (÷ M0), E2b | 1.00 ± 0.05 | 1.00 ± 0.04 | 1.04 ± 0.05 | **2.09 ± 0.10** | **8.3 ± 0.8** | **11 ± 4** |
| same, E2 (θ = 0.024) | 1.00 | 0.98 | 1.03 | 1.21 | 1.89 | 2.58 |
| P0 directional prediction at realised 2Ns | | | | 2.0 (2Ns 1.6) | 158 (2Ns 8.1) | 2,400 (2Ns 11.4) |
| DAF of R-only variants (÷ M0) | 1 | 0.99 | 0.97 | **0.78** | **0.69** | **0.31** |
| CpG positional turnover over 0.5/u, E2b | 0.98 | 0.97 | 0.96 | 0.98 | **0.51** | **0.16** |
| same window, E2 | 0.99 | 0.98 | 0.98 | 0.98 | 0.91 | 0.78 |
| \|ΔR\| over 0.5/u, E2b | 1.86 | 1.39 | 1.58 | 0.22 | 0.05 | 0.001 |

1. **Passenger stays invisible at realistic θ.** P/D (1.00, 1.04) and R-only DAF (0.97–0.99 of neutral) are neutral.
   This is the most robust result in the project so far: it holds at both θ values and both κ values.
2. **Causal signals are stronger at realistic θ, not weaker.** P/D is 2.1 (κ=2) to 8–11 (κ=10), versus 1.2–2.6 at high
   θ. R-only DAF falls by 22–69%.
   - At high θ, recurrent CpG mutation and many co-segregating variants at a non-recombining locus blunt selection.
   - At realistic θ, the weak-selection case matches P0's origin-fixation prediction (2.09 vs 2.0 at 2Ns = 1.6).
   - Strong-step cases stay far below P0's directional numbers (8 vs 158), because of stabilising selection and
     compensation (`theory.pdf` §11).
3. **Correction to E2: "methylation conserved through complete CpG turnover" depends on θ when steps are strong.**
   Per-step 2Ns ≈ 1.6 (C κ=2): full turnover with R conserved, at both θ. **This part of the E2 claim holds.** Per-step
   2Ns ≈ 8–11: turnover over 0.5/u drops from 0.78–0.91 (E2) to 0.16–0.51 at realistic θ, because CpG positions
   mostly freeze. The likely mechanism:
   - Compensating a lost CpG means gaining one elsewhere.
   - At high θ the gain often arises while the loss is still segregating (stochastic tunnelling).
   - At realistic θ the loss usually has to fix first, which is a deleterious intermediate costing 2Ns ≈ 8–11, so
     valley crossing is rare.
   - The mechanism is inferred from the θ dependence; the tunnelling rate was not measured directly.
   **Revised statement:** methylation conserved while CpG positions turn over requires a degenerate map with weak
   per-CpG selection (2Ns ≲ 2). Strong per-CpG selection conserves both.

## Implications

- **For E4 and real data:** polymorphism statistics from E2 were conservative. At human-like θ the causal signal in
  P/D and DAF is larger. Detection bias (E4) is now the main threat, not θ.
- **For CHATGPT Q3** ("how conserved can methylation be while sequence diverges?"): it can be fully conserved through
  full CpG turnover only in the weak per-site regime. That is the same regime where per-site tests have the least
  power (P/D ≈ 2, DAF −22%). Sequence-level conservation scores cannot see the trait-level conservation, but
  frequency statistics can, weakly.
- **Comparison caveat:** E2 reports MAF and E2b reports DAF. Panel B therefore normalises each to its own neutral
  baseline.
