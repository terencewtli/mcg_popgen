# Endpoint for the simulation phase, and hand-off to real data

Agreed 2026-10-01. The aim is a finished, validated simulation unit within about two days. After that the project moves
to the cluster and works exclusively on real data.

## The question the simulation phase answers

> Which population-genetic signatures can distinguish causal from passenger methylation, at what selection strength,
> and which of them survive realistic θ and meQTL detection bias?

**Deliverable:** a short list of concrete, falsifiable predictions for real data. Each one links to the simulation that
supports it.

## Checklist

| # | item | status | write-up |
|---|---|---|---|
| 1 | E0 simulator validation | done (18/18 checks) | `md/E0.md` |
| 2 | E1 passive CpG footprint (methylation never in fitness) | done (2.6× within 25 bp) | `md/E1.md` |
| 3 | P0 primate-scale priors (mutation vs selection on CpGs) | done (P/D breaks the o/e ridge) | `md/PRIORS.md` |
| 4 | E2 causal vs passenger, stabilising selection | done | `md/E2.md` |
| 5 | **E2b: E2 at realistic θ** (4Nu = 0.001, CpG recurrence ≈ human) | done: passenger still neutral; causal P/D 2–11; CpG turnover with conserved R only at weak per-CpG 2Ns | `md/E2b.md` |
| 6 | **E3: combined model**: TF and methylation both causal (g = O·exp(−R/R₀)), plus a methylation-sensitive TF variant (m lowers O, cf. NRF1) | next | `md/E3.md` |
| 7 | **E4: meQTL detection layer**: analytic power applied to simulated variants (detected if n·2p(1−p)β²/σ² passes the threshold; n = 500, 5,000, 30,000). Does the CpG-SNP vs motif-SNP difference in effect–frequency survive detection bias under M0, P and C? | | `md/E4.md` |
| 8 | One real-N check (N = 10⁴, no rescaling) for the headline condition: overnight on the laptop, or the first cluster job | | in `md/E2b.md` or `md/E4.md` |
| 9 | **Consolidated summary**: `md/SUMMARY.md` (numbered predictions → supporting experiment → real-data analysis A–D in `md/REAL_DATA.md`) plus a figure-led PDF | | `md/SUMMARY.md` |

**Deliberately out of scope**, so the phase actually ends:
- large phase-diagram sweeps (two κ values per model bracket the boundary);
- trans factors;
- bistable CpG islands;
- the E1 r = 50 open item, unless a headline claim depends on it.

## What it is worth (honest view)

- **As a unit it is calibration and a null model for real data.** On its own it is probably a methods note, not a
  paper. Pooled meQTL selection inference already exists (Hatcher 2021 bioRxiv).
- **Most likely to hold up as findings:**
  - methylation conserved through complete CpG turnover (degenerate map + weak steps);
  - passenger methylation being invisible to frequency tests;
  - the germline-vs-soma mismatch in gnomAD-style CpG mutability.
- **The cluster phase decides whether these appear in real data.**

## Cluster hand-off and two-account working convention

**Split:**
- **Laptop session: simulations.** E2b–E4, plus large sweeps sent to the cluster as SLiM job arrays if needed. SLiM
  installs through conda with no admin rights.
- **Cluster session: real data.** First choice is analysis B (CpG-loss P/D: polymorphism + primate alignments +
  germline methylomes; public data only). Alternative: fitting the methylation map to HPRC2 haplotype-resolved
  methylation (`md/REAL_DATA.md` §4.1).

**Repository rules** (one repo, `terencewtli/mcg_popgen`, two writers):
1. **The laptop owns:** `md/` (except `md/realdata/`), `scripts/`, `tex/`, `csv/e*`, `csv/p0`, `pdf/e*`, `pdf/p0`,
   `pdf/theory.pdf`. These are written by `scripts/sync_to_github.sh` from the laptop working directory. Never edit them
   on the cluster.
2. **The cluster owns:** `realdata/` (code, small tables, figures, and its notes in `realdata/md/`). The cluster works
   directly in a clone of the repo and does not run the sync script.
3. **Before committing, both sides pull.** The sync script now runs `git pull --ff-only` first.
4. The sync script's `--delete` rsyncs only touch laptop-owned directories, so a laptop sync never removes `realdata/`.
5. **Large data** (methylomes, VCFs, alignments) stays on the cluster and is never committed. Record paths and md5s in
   `realdata/README.md`.

**First thing a cluster session should read:** this file, `md/SUMMARY.md` (once written), `md/REAL_DATA.md` and
`pdf/theory.pdf`.
