# Precedent: population-genetic models applied to DNA methylation (lit check 2026-10-01)

Citations were checked against PubMed or bioRxiv metadata by a search agent. Flagged entries were not fully verified.
Gaffney 2012 was my own mis-suggestion: it is an eQTL paper, and the methylation companion is Bell 2011.

## Where this project sits

| Strand | Has a sequence→methylation map? | Has selection/drift? | Key papers |
|---|---|---|---|
| Neutral CpG decay gated by methylation | implicit (methylated vs not) | mutation + drift (+BGC) | Sved & Bird 1990 PNAS; Cohen, Kenigsberg & Tanay 2011 Cell; Molaro 2011 Cell; Chandra & Gao 2026 PLoS Genet (gnomAD methylated vs unmethylated CpG rates by context, the best source for r) |
| Epiallele population genetics (plants) | no: methylation is its own mutating allele | yes, SFS-based | Charlesworth & Jain 2014 Genetics; van der Graaf 2015 PNAS; Wang & Fan 2014 GBE (D^m); Muyle, Ross-Ibarra, Seymour & Gaut 2021 Genetics (4Nₑs ≈ 1.4 on gene-body methylation); Shahryary 2020 (AlphaBeta) |
| Comparative / phylogenetic methylation | no | OU / CTMC / sign test | Vilgalys 2019 MBE (OU, Papio); Hernando-Herraez 2013, 2015; Qu 2018 Genome Res; Capra & Kostka 2014; Ma, Starr, Gokhman & Fraser 2026 bioRxiv 10.64898/2026.01.20.700710 |
| Experimental sequence determinants | yes (empirical) | no | Lienert 2011 Nat Genet; Krebs 2014 and Wachter 2014 eLife; Stadler 2011 Nature (TF → LMR, i.e. passenger); Domcke 2015 Nature (methylation blocks NRF1, i.e. reverse arrow); Long 2013 eLife |
| meQTL causality | statistical | mostly no | Banovich 2014 PLoS Genet (meQTLs often downstream of TF binding); Hawe 2022 Nat Genet; Husquin 2018 Genome Biol. **No SFS or origin-fixation fit to human meQTL alleles found.** |
| TF-binding-site population genetics | occupancy map, no methylation | full mutation–selection–drift | Berg, Willmann & Lässig 2004; Mustonen & Lässig 2005; Mustonen 2008 PNAS; Tuğrul 2015 PLoS Genet; Lynch & Hagner 2015 PNAS |
| Forward simulation with sequence-computed methylation | | | **none found.** Nearest: Branciamore 2014 PNAS (TE silencing); Cohen 2011's simulated alignments |

## Synthesis

This project fuses two existing lines: Cohen 2011's methylation-gated mutation and Mustonen & Lässig's binding-energy
fitness landscapes. The Schübeler-lab experiments are the empirical justification for the map. The pieces are old; the
combination appears unexplored. Two untested questions follow directly:
1. Does neutral methylation-gated CpG decay generate apparent lineage- or pathway-specific sign-test signals (Ma et al.)?
2. Can passenger and causal methylation be separated by meQTL Δm × allele-frequency statistics?

Flags: Vidalis 2016 was checked via the publisher, not PubMed. Ayres 2025 J Theor Biol has no DOI captured. Si 2023 and
Zhang 2025 are preprints only.
