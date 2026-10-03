# Literature search 2026-10-03: CpG / CGI population-genetic theory and methylation-aware binding-site evolution

Scope: is there prior rigorous population-genetic theory or simulation of CpG / CpG-island evolution under mutation,
selection and drift, and has any Lässig-style binding-site model (Berg/Willmann/Lässig 2004; Mustonen & Lässig 2005;
Mustonen 2008; Tuğrul 2015) been extended to DNA methylation? This extends `md/LITERATURE.md` (2026-10-01), which is not edited here.

Verification: every entry in the tables was checked against PubMed E-utilities (title, first author, year, journal,
PMID, DOI) or, for preprints, the bioRxiv API. Abstracts were read for every entry marked with a summary. Entries already
in `md/LITERATURE.md` are re-listed only where they are needed for the comparison and are marked (L).

## Verdict

**The combination looks novel.** Searches across PubMed, bioRxiv and the web found no model in which (i) methylation is
computed from each haplotype's sequence (CpG density plus TF-occupancy protection), (ii) CpG mutation is gated by that
computed methylation, and (iii) fitness acts on methylation or on occupancy inside a forward mutation–selection–drift
simulation. No Lässig-style binding-energy model has been extended to methylation. The papers we found either build that
framework with sequence-independent mutation (Berg 2004; Mustonen & Lässig 2005; Mustonen 2008; Tuğrul 2015;
Haldane 2014; Stewart & Plotkin 2012; Lynch & Hagner 2015), or treat CpG deamination only as a neutral accelerator of motif
*creation* (Behrens & Vingron 2010; Zemojtel 2011). Each piece exists separately:
- **Methylation-gated neutral CpG decay** has explicit models: Sved & Bird 1990 (analytic equilibrium); Arndt 2003 and
  Hwang & Green 2004 (neighbour-dependent substitution); Cohen, Kenigsberg & Tanay 2011 (evolutionary regimes with gBGC,
  with population data arguing for minimal selection); Mugal 2015 (methylation-aware equilibrium GC).
- **Selection on methylated CpGs** has been inferred from frequency data: Agarwal & Przeworski 2021 (mutation saturation);
  Ying & Huttley 2011; Harpak 2016 (recurrent mutation shapes the SFS); Si, Kang & Zöllner 2023 preprint.
- **Selection on methylation as a trait** appears only as an effect–MAF fit (Hatcher 2021 preprint) or as plant epiallele
  SFS tests (L).
- **Mechanistic sequence→methylation maps** exist with no evolution: Haerter 2014 / Lövkvist 2016 (CpG-density
  cooperativity); Schübeler-lab experiments (L); Long 2016.
- **Forward simulations of methylation**: individual-based simulations treat methylation as a sequence-independent
  epiallele or quantitative-genetic layer (Branciamore 2014; Geoghegan & Spencer 2012; Furrow & Feldman 2014;
  Klironomos 2013; Ayres 2025). The only simulator we found that runs on genealogies (MethEvolSIM, extending Grosser &
  Metzler 2020) evolves methylation states, not sequence.

The project's "passive retention" result (F1) is foreshadowed verbally by Cohen 2011 and by Long 2016 ("protection is
DNA-encoded"). The ape-scale observation that germline hypomethylation preserves CpG-rich sequence is now in Son et al.
2026 (bioRxiv, T2T apes). The causal-vs-passenger P/D, DAF and meQTL contrast has no precedent we could find. Hartl 2019
(CpG density raises promoter activity independent of DNMTs) is worth citing: it is direct evidence for a third fitness
route, CpG→function without methylation, which the model does not include.

**Five closest precedents**
1. **Cohen, Kenigsberg & Tanay 2011 Cell** (PMID 21620139): explicit evolutionary models of CGI maintenance through
   low-methylation-gated deamination plus gBGC, with population data arguing for minimal selection. It is the neutral
   half of this project, but it has no sequence→methylation map, no TF and no forward simulation.
2. **Mustonen & Lässig 2005 PNAS / Mustonen 2008 PNAS / Tuğrul 2015 PLoS Genet** (binding-energy fitness landscapes with
   full mutation–selection–drift). This is the selection framework the project borrows; none of these papers include
   methylation or CpG-gated mutation.
3. **Agarwal & Przeworski 2021 eLife** (PMID 34806592): treats germline-methylated CpGs as a mutation-saturation
   experiment and infers strong selection from missing polymorphism. It is the closest frequency-based selection
   inference that conditions on methylation-gated mutation. It covers coding sites only and does not address causal vs
   passenger.
4. **Zemojtel 2011 GBE** (PMID 22016335) **with Behrens & Vingron 2010 J Comput Biol** (PMID 21128851): CpG deamination
   creates TFBSs efficiently, and the waiting-time model has neighbour-dependent mutation. These are the only
   binding-site evolution models we found with methylation-type mutation, and they are neutral, with no selection or
   methylation map.
5. **Hatcher et al. 2021 bioRxiv** (10.1101/2021.11.25.469994): a BayesS-type joint effect–MAF model for DNAm sites. It
   is the only selection inference on methylation phenotypes in humans. It pools across sites, has no mechanism
   stratification and no matched null, and we found no journal version (bioRxiv API, 2026-10-03).

Runner-up: **Son et al. 2026 bioRxiv** (germline hypomethylation preserves CpG reservoirs in T2T apes). It is the most
recent empirical support for F1/R1, but it is descriptive, with no model.

## Strand 1: CpG / CGI origin, maintenance, decay

| Paper | IDs | What it does | Seq→meth map? | Selection / drift / sim? |
|---|---|---|---|---|
| Bird AP 1980 Nucleic Acids Res, "DNA methylation and the frequency of CpG in animal DNA" | PMID 6253938; 10.1093/nar/8.7.1499 | Original link between genome methylation and CpG deficit. | no | no |
| Bird AP 1986 Nature, "CpG-rich islands and the function of DNA methylation" | PMID 2423876; 10.1038/321209a0 | Classic proposal that CGIs are unmethylated footprints. | no | no |
| Sved J, Bird A 1990 PNAS, "The expected equilibrium of the CpG dinucleotide in vertebrate genomes under a mutation model" (L) | PMID 2352943; 10.1073/pnas.87.12.4692 | Analytic mutation-only equilibrium CpG o/e given methylation-driven deamination. | methylated fraction as parameter | mutation only |
| Arndt PF, Burge CB, Hwa T 2003 J Comput Biol, "DNA sequence evolution with neighbor-dependent mutation" | PMID 12935330; 10.1089/10665270360688039 | Exactly solvable neighbour-dependent (CpG) substitution model; equilibrium dinucleotide frequencies. | no | substitution process only |
| Hwang DG, Green P 2004 PNAS, "Bayesian MCMC sequence analysis reveals varying neutral substitution patterns in mammalian evolution" | PMID 15292512; 10.1073/pnas.0404142101 | Context-dependent (CpG) neutral substitution inference across mammals. | no | neutral |
| Fryxell KJ, Moon WJ 2005 MBE, "CpG mutation rates in the human genome are highly dependent on local GC content" | PMID 15537806; 10.1093/molbev/msi043 | CpG transition rate depends on local GC (duplex stability). | no | no |
| Zhao Z, Jiang C 2007 MBE, "Methylation-dependent transition rates are dependent on local sequence lengths and genomic regions" | PMID 17056644; 10.1093/molbev/msl156 | Empirical CpG transition-rate heterogeneity. | no | no |
| Elango N, Yi SV 2008 MBE, "DNA methylation and structural and functional bimodality of vertebrate promoters" | PMID 18469331; 10.1093/molbev/msn110 | Links promoter CpG bimodality to germline methylation across vertebrates. | no | no |
| Walser JC, Ponger L, Furano AV 2008 Genome Res, "CpG dinucleotides and the mutation rate of non-CpG DNA" | PMID 18550801; 10.1101/gr.076455.108 | CpG density modulates non-CpG mutation; relevant to composition-dependent rates. | no | no |
| Cohen NM, Kenigsberg E, Tanay A 2011 Cell (L) | PMID 21620139; 10.1016/j.cell.2011.04.024 | **Explicit evolutionary models**: low methylation → slow deamination maintains CGIs; gBGC stabilises methylated CGIs; population data → minimal selection. | methylation as input (measured) | mutation + BGC; selection tested, found minimal |
| Mugal CF, Ellegren H 2011 Genome Biol, "Substitution rate variation at human CpG sites correlates with non-CpG divergence, methylation level and GC content" | PMID 21696599; 10.1186/gb-2011-12-6-r58 | CpG substitution rate vs sperm methylation and GC. | no | no |
| Mugal CF, Arndt PF, Holm L, Ellegren H 2015 G3, "Evolutionary consequences of DNA methylation on the GC content in vertebrate genomes" | PMID 25591920; 10.1534/g3.114.015545 | Neighbour-dependent substitution model with sperm methylation maps; methylation shifts equilibrium GC by about 15%. | methylation as input | substitution/equilibrium; no selection |
| Berglund J, Quilez J, Arndt PF, Webster MT 2014 GBE, "Germline methylation patterns determine the distribution of recombination events in the dog genome" | PMID 25527838; 10.1093/gbe/evu282 | Disentangles CpG mutability and gBGC at dog CGIs; gBGC creates new CGIs absent PRDM9. | no | gBGC + mutation (substitution model) |
| Glémin S et al. 2015 Genome Res, "Quantification of GC-biased gene conversion in the human genome" | PMID 25995268; 10.1101/gr.185488.114 | SFS-based gBGC strength estimates in humans. Use it to set gBGC controls (R2). | no | SFS-based |
| Bergman J, Schierup MH 2021 Genetics, "Population dynamics of GC-changing mutations in humans and great apes" | PMID 34081117; 10.1093/genetics/iyab083 | CpG transitions most affected by gBGC in human/chimp; reduced in gorilla/orangutan. Directly relevant to R2's W↔S control. | no | gBGC from segregating frequencies |
| Miyahara H et al. 2015 BMC Genomics, "Factors to preserve CpG-rich sequences in methylated CpG islands" | PMID 25879481; 10.1186/s12864-015-1286-x | Partitions maintenance of methylated CGIs into sperm hypomethylation vs CpG selection vs BGC; mostly germline hypomethylation. | no | descriptive, divergence |
| Branciamore S, Chen ZX, Riggs AD, Rodin SN 2010 PNAS, "CpG island clusters and pro-epigenetic selection for CpGs in protein-coding exons of HOX..." | PMID 20716685; 10.1073/pnas.1010506107 | Argues for selection preserving silent CpGs in developmental TF genes ("pro-epigenetic selection"). | no | inferred selection, no model |
| Long HK ... Klose RJ 2013 eLife, "Epigenetic conservation at gene regulatory elements revealed by non-methylated DNA profiling in seven vertebrates" (L) | PMID 23467541; 10.7554/eLife.00348 | Non-methylated islands across vertebrates. | empirical | no |
| Long HK, King HW ... Klose RJ 2016 NAR, "Protection of CpG islands from DNA methylation is DNA-encoded and evolutionarily conserved" | PMID 27084945; 10.1093/nar/gkw258 | Human chromosome in mouse: CGI protection is sequence-encoded; distal elements depend on species TF occupancy. It is direct support for the map. | yes (empirical) | no |
| Rademacher K ... Horsthemke B 2014 GBE, "Evolutionary origin and methylation status of human intronic CpG islands that are not present in mouse" | PMID 24923327; 10.1093/gbe/evu125 | Lineage-specific CGIs (some retrotransposed), sequence-dependent methylation fate. | qualitative | no |
| Kocher AA ... Noonan JP 2024 Genome Biol, "CpG island turnover events predict evolutionary changes in enhancer activity" | PMID 38872220; 10.1186/s13059-024-03300-z | CGI turnover across 9 mammals tracks species-specific enhancer activity. Empirical counterpart to F8/R5. | no | no (comparative) |
| You JS, Pierce S, Liang G, Jones PA 2025 PNAS, "Roles of transposable elements and DNA methylation in the formation of CpG islands and CpG-depleted regulatory elements" | PMID 41134632; 10.1073/pnas.2502963122 | Verbal model: CGIs persist because TE insertion is counterselected; ancestral CpG-rich genome eroded by methylation. | no | verbal |
| Son DR, Loh EY, Jeong H, Ma J, Eichler EE, Yi SV 2026 bioRxiv, "Germline hypomethylation shapes dynamic CpG reservoirs in ape genomes" | PMID 42368016; 10.64898/2026.06.15.732472 | T2T apes + long/short-read germline methylomes: sperm methylation drives CpG erosion; hypomethylated satellites/SDs act as CpG reservoirs. Most recent empirical support for F1/R1. | no | descriptive |
| Hartl D, Krebs AR ... Schübeler D 2019 Genome Res, "CG dinucleotides enhance promoter activity independent of DNA methylation" | PMID 30709850; 10.1101/gr.241653.118 | CpG density raises promoter output with TF motifs, independent of DNMTs. **Third fitness route (CpG itself, not methylation)**; the model omits it. | yes (empirical) | no |

## Strand 2: TF-binding-site evolution and methylation

| Paper | IDs | What it does | Methylation? | Selection / drift / sim? |
|---|---|---|---|---|
| Berg J, Willmann S, Lässig M 2004 BMC Evol Biol, "Adaptive evolution of transcription factor binding sites" (L) | PMID 15511291; 10.1186/1471-2148-4-42 | Binding-energy fitness landscape + mutation–selection–drift; fast adaptive site formation. | no | yes, full pop-gen |
| Mustonen V, Lässig M 2005 PNAS, "Evolutionary population genetics of promoters" (L) | PMID 16236723; 10.1073/pnas.0505537102 | Infers selection on sites from binding-energy model. | no | yes |
| Mustonen V, Kinney J, Callan CG, Lässig M 2008 PNAS, "Energy-dependent fitness" (L) | PMID 18723669; 10.1073/pnas.0805909105 | Yeast TFBS fitness vs binding energy. | no | yes |
| Haldane A, Manhart M, Morozov AV 2014 PLoS Comput Biol, "Biophysical fitness landscapes for transcription factor binding sites" | PMID 25010228; 10.1371/journal.pcbi.1003683 | Infers fitness(binding) under a monomorphic-population model; epistasis common. | no | yes (SSWM) |
| Tuğrul M, Paixão T, Barton NH, Tkačik G 2015 PLoS Genet, "Dynamics of Transcription Factor Binding Site Evolution" (L) | PMID 26545200; 10.1371/journal.pgen.1005639 | Waiting times for site emergence under selection; site length/population size dependence. | no | yes, simulation + theory |
| Stewart AJ, Hannenhalli S, Plotkin JB 2012 Genetics, "Why transcription factor binding sites are ten nucleotides long" | PMID 22887818; 10.1534/genetics.112.143370 | Pop-gen of site length/specificity. | no | yes |
| Lynch M, Hagner K 2015 PNAS (L) | PMID 25535374; 10.1073/pnas.1421641112 | Drift-barrier meandering of binding interfaces. | no | yes |
| Behrens S, Vingron M 2010 J Comput Biol, "Studying the evolution of promoter sequences: a waiting time problem" | PMID 21128851; 10.1089/cmb.2010.0084 | Waiting time for k-mer emergence under i.i.d. vs neighbour-dependent (CpG) mutation; CpG deamination accelerates site creation. | CpG-mutation only | neutral, no selection |
| Zemojtel T ... Vingron M 2011 GBE, "CpG deamination creates transcription factor-binding sites with high efficiency" | PMID 22016335; 10.1093/gbe/evr107 | Empirical: deamination creates TFBSs (p53, Oct4, Myc) efficiently. | CpG-mutation only | no |
| Hernando-Herraez I ... Marques-Bonet T 2015 NAR, "The interplay between DNA methylation and sequence divergence in recent human evolution" | PMID 26170231; 10.1093/nar/gkv693 | Human-specific substitutions in TFBSs within human DMRs; methylation gains coupled to local sequence change. | measured | comparative |
| Sahm A, Koch P, Horvath S, Hoffmann S 2021 MBE, "An analysis of methylome evolution in primates" | PMID 34175932; 10.1093/molbev/msab189 | Phylogenetic models on methylation; conservation correlates with TF-binding density. | no | phylogenetic |
| Kumar P ... Singh U 2025 Transcription, "CGGBP1 from higher amniotes restricts cytosine methylation and drives a GC-bias in transcription factor-binding sites at repressed promoters" | PMID 40740140; 10.1080/21541264.2025.2533598 | A non-TF protector keeps GC-rich TFBSs unmethylated; links protection to GC retention across >100 species. Same "protection → retention" logic as F1. | yes (empirical) | comparative |
| Yin Y ... Taipale J 2017 Science, "Impact of cytosine methylation on DNA binding specificities of human transcription factors" | PMID 28473536; 10.1126/science.aaj2239 | Methyl-sensitive/-preferring TFs (basis for the MS-TF fitness model). | biochemical | no |
| Domcke S ... Schübeler D 2015 Nature (L) | PMID 26675734; 10.1038/nature16462 | Methylation blocks NRF1 (reverse arrow). | yes | no |
| Stadler MB ... Schübeler D 2011 Nature (L) | PMID 22170606; 10.1038/nature10716 | TF binding creates LMRs (passenger arrow). | yes | no |

No paper was found that puts methylation, methylation-gated mutation, or methyl-sensitive binding into a Lässig-type
fitness landscape.

## Strand 3: population-genetic inference on CpGs / methylation

| Paper | IDs | What it does | Relation |
|---|---|---|---|
| Xia J, Han L, Zhao Z 2012 BMC Genomics, "Investigating the relationship of DNA methylation with mutation rate and allele frequency in the human genome" | PMID 23281708; 10.1186/1471-2164-13-S8-S7 | SNP density and allele frequency vs hESC methylation level per CpG. | Descriptive; methylation level correlates with allele frequency, with no selection model. |
| Ying H, Huttley G 2011 GBE, "Exploiting CpG hypermutability to identify phenotypically significant variation within human protein-coding genes" | PMID 21398426; 10.1093/gbe/evr021 | Uses CpG hypermutability as a baseline for selection inference (coding). | P/D logic at CpGs, coding only. |
| Cocozza S, Akhtar MM, Miele G, Monticelli A 2011 PLoS One, "CpG islands undermethylation in human genomic regions under selective pressure" | PMID 21829712; 10.1371/journal.pone.0023156 | CGIs in selected regions are less methylated and less polymorphic. | Correlational; confounds mutation and selection, the ridge F2 warns about. |
| Scala G ... Cocozza S 2014 PLoS One | PMID 25474578; 10.1371/journal.pone.0114432 | SNP frequency classes around TSSs inside vs outside CGIs; purifying selection + BGC. | Frequency-stratified, no methylation-gated null. |
| Harpak A, Bhaskar A, Pritchard JK 2016 PLoS Genet, "Mutation rate variation is a primary determinant of the distribution of allele frequencies in humans" | PMID 27977673; 10.1371/journal.pgen.1006489 | Recurrent mutation at hypermutable (CpG) sites distorts SFS. | The CpG DAF/SFS null must include recurrence (relevant to real-N runs and to R2/R3). |
| Agarwal I, Przeworski M 2021 eLife, "Mutation saturation for fitness effects at human CpG sites" | PMID 34806592; 10.7554/eLife.71513 | Germline-methylated CpGs are about 99% polymorphic for synonymous C>T at n = 390k; absence of polymorphism implies strong selection. | Closest frequency-based selection inference conditioned on methylation-gated mutation (coding). |
| Agarwal I, Fuller ZL, Myers SR, Przeworski M 2023 eLife | PMID 36648429; 10.7554/eLife.83172 | Mutation–selection-balance posterior on hs for LoF. | Method template for per-element selection from CpG frequency. |
| Dukler N, Mughal MR, Ramani R, Huang YF, Siepel A 2022 Nat Commun, "Extreme purifying selection against point mutations in the human genome" | PMID 35879308; 10.1038/s41467-022-31872-6 | ExtRaINSIGHT: rare-variant depletion after mutation-rate control; little ultraselection in TFBSs. | Noncoding selection baseline. |
| Arbiza L ... Siepel A 2013 Nat Genet, "Genome-wide inference of natural selection on human transcription factor binding sites" | PMID 23749186; 10.1038/ng.2658 | INSIGHT P/D-type inference on TFBSs. | P/D framework for R2; no methylation stratification. |
| Carlson J ... Zöllner S 2018 Nat Commun, "Extremely rare variants reveal patterns of germline mutation rate heterogeneity in humans" | PMID 30218074; 10.1038/s41467-018-05936-5 | Mutation-rate model from ERVs including methylation covariates. | Mutation-model input. |
| Seplyarskiy VB ... Sunyaev SR 2021 Science, "Population sequencing data reveal a compendium of mutational processes in the human germ line" | PMID 34385354; 10.1126/science.aba7408 | Decomposes germline mutation processes from population data. | Mutation-model input. |
| Chen S ... Karczewski KJ 2024 Nature, "A genomic mutational constraint map using variation in 76,156 human genomes" (Gnocchi) | PMID 38057664; 10.1038/s41586-023-06045-0 | Noncoding constraint with context- and methylation-adjusted mutation model. | Constraint baseline; not mechanism-stratified. |
| Chandra S, Gao Z 2026 PLoS Genet (L) | PMID 42224293; 10.1371/journal.pgen.1011957 | Methylated vs unmethylated CpG mutation rates by 4/6-mer context from gnomAD with recurrence. | Best source for r. |
| Si Y, Kang HM, Zöllner S 2023 bioRxiv, "Germline CpG methylation signatures in the human population inferred from genetic polymorphism" | 10.1101/2023.03.24.534151 (v2 2024-01-03; no published version per bioRxiv API) | HMM infers germline methylation from TOPMed polymorphism; hypermethylated monomorphic CpGs flagged as likely constrained. | Note circularity for R2: it infers methylation *from* polymorphism. |
| Hatcher C ... Min JL 2021 bioRxiv (L) | 10.1101/2021.11.25.469994 | Effect–MAF selection (S) for 2000 DNAm sites, ALSPAC. | Closest meQTL-selection precedent. |
| Gutierrez-Arcelus M ... Dermitzakis ET 2013 eLife, "Passive and active DNA methylation and the interplay with genetic variation in gene regulation" | PMID 23755361; 10.7554/eLife.00523 | Argues that inter-individual methylation is often a passive consequence of TF binding/expression. | Empirical "passenger" framing (cite with Banovich 2014). |
| Do C ... Tycko B 2016 AJHG, "Mechanisms and disease associations of haplotype-dependent allele-specific DNA methylation" | PMID 27153397; 10.1016/j.ajhg.2016.03.027 | hap-ASM explained largely by TF-binding-site disruption. | Supports passenger arrow; mechanism for R3. |
| Onuchic V ... Milosavljevic A 2018 Science, "Allele-specific epigenome maps reveal sequence-dependent stochastic switching at regulatory loci" | PMID 30139913; 10.1126/science.aar3146 | Allele-specific methylation linked to TF motif disruption and stochastic switching. | Map-fitting target. |
| Charlesworth B, Jain K 2014 Genetics (L) | PMID 25230951; 10.1534/genetics.114.167973 | Selection–drift with reversible, arbitrarily high mutation (epialleles). | Theory for high-rate reversible states; methylation not sequence-encoded. |
| Wang J, Fan C 2014 GBE (L) | PMID 25539727; 10.1093/gbe/evu271 | SMP-SFS neutrality test for methylation. | Plant epiallele SFS. |
| Muyle A, Ross-Ibarra J, Seymour DK, Gaut BS 2021 Genetics (L) | PMID 33871638; 10.1093/genetics/iyab061 | Selection on gene-body methylation from epiallele SFS. | Plant; not sequence-encoded. |
| Costa CE ... Lea AJ 2025 Mol Ecol, "Genetic architecture of immune cell DNA methylation in the rhesus macaque" | PMID 39582237; 10.1111/mec.17576 | meQTL in a natural primate population; enriched at methyl-sensitive TF sites. | Possible R3 dataset; no selection analysis. |

## Strand 4: simulations / models where methylation is computed or selected

| Paper | IDs | What it does | Seq→meth? | Selection / sim? |
|---|---|---|---|---|
| Haerter JO, Lövkvist C, Dodd IB, Sneppen K 2014 NAR, "Collaboration between CpG sites is needed for stable somatic inheritance of DNA methylation states" | PMID 24288373; 10.1093/nar/gkt1235 | Stochastic cooperative methylation model; CpG density → bistable CGI states. | **yes (mechanistic, from CpG density)** | none (somatic, no evolution) |
| Lövkvist C, Dodd IB, Sneppen K, Haerter JO 2016 NAR, "DNA methylation in human epigenomes depends on local topology of CpG sites" | PMID 26932361; 10.1093/nar/gkw124 | Same model fitted to human methylomes: CpG spacing predicts methylation. | **yes** | none. Candidate drop-in replacement for the hand-set map. |
| Sormani G et al. 2016 Mol Biosyst, "Stabilization of epigenetic states of CpG islands by local cooperation" | PMID 26923344; 10.1039/c6mb00044d | Related CpG-cooperation model. | yes | none |
| Grosser K, Metzler D 2020 BMC Bioinformatics, "Modeling methylation dynamics with simultaneous changes in CpG islands" | PMID 32183713; 10.1186/s12859-020-3438-5 | IWE-SSE: CTMC for methylation on cell/species trees with island-wide events. (MethEvolSIM, CRAN 2024, extends it to coalescent trees; no paper found.) | no | neutral epimutation along trees |
| Capra JA, Kostka D 2014 Bioinformatics (L) | PMID 25161227; 10.1093/bioinformatics/btu445 | Phylogenetic methylation models. | no | neutral |
| Branciamore S, Rodin AS, Riggs AD, Rodin SN 2014 PNAS (L) | PMID 24733912; 10.1073/pnas.1402585111 | Drift analysis + simulation of stochastic epigenetic marks after gene duplication. | no | selection + drift, epigenetic state not sequence-encoded |
| Geoghegan JL, Spencer HG 2012 Theor Popul Biol, "Population-epigenetic models of selection" | PMID 21855559; 10.1016/j.tpb.2011.08.001 | Epiallele selection models. | no | selection, analytic |
| Furrow RE, Feldman MW 2014 Evolution, "Genetic variation and the evolution of epigenetic regulation" | PMID 24588347; 10.1111/evo.12225 | Modifier models for epigenetic regulation. | no | yes |
| Klironomos FD, Berg J, Collins S 2013 BioEssays, "How epigenetic mutations can affect genetic evolution: model and mechanism" | PMID 23580343; 10.1002/bies.201200169 | Simulation: fast epimutations alter adaptive paths. | no | yes |
| Kronholm I, Collins S 2016 Mol Ecol, "Epigenetic mutations can both help and hinder adaptive evolution" | PMID 26139359; 10.1111/mec.13296 | Simulations of epimutation in adaptation. | no | yes |
| Ayres L, Bovenhuis H, Calus MPL 2025 J Theor Biol, "A single-locus quantitative genetic model incorporating DNA methylation" | PMID 40189137; 10.1016/j.jtbi.2025.112110 | Fisher decomposition with a methylation layer (basic vs expressed genetic value). | methylation as allele-modulated quantity | quantitative genetics, no drift sim |
| López-Catalina A et al. 2025 J Anim Breed Genet | PMID 39868874; 10.1111/jbg.12925 | Breeding-value model with simulated methylation effects. | no | no pop-gen |
| Mueller SA, Merondun J, Lečić S, Wolf JBW 2025 Nat Commun, "Epigenetic variation in light of population genetic practice" | PMID 39863592; 10.1038/s41467-025-55989-6 | Perspective on integrating epigenetic data into pop-gen (genetically determined vs autonomous). | — | review |

No forward simulation (SLiM or other) was found in which methylation is computed from the evolving sequence and feeds
back on CpG mutation or fitness.

## Strand 5: recent long-read / T2T methylation and mutation models

| Paper | IDs | Relevance |
|---|---|---|
| Gershman A ... Timp W 2022 Science, "Epigenetic patterns in a complete human genome" | PMID 35357915; 10.1126/science.abj5089 | T2T-CHM13 5mC baseline. |
| Kolmogorov M ... 2023 Nat Methods, "Scalable Nanopore sequencing of human genomes provides a comprehensive view of haplotype-resolved variation and methylation" | PMID 37710018; 10.1038/s41592-023-01993-x | Haplotype-resolved methylation at scale. |
| Gustafson JA ... 2024 Genome Res, "High-coverage nanopore sequencing of samples from the 1000 Genomes Project..." | PMID 39358015; 10.1101/gr.279273.124 | 1KG-ONT methylation for meQTL-like analyses (C). |
| Porubsky D ... Eichler EE 2025 Nature, "Human de novo mutation rates from a four-generation pedigree reference" | PMID 40269156; 10.1038/s41586-025-08922-2 | Pedigree DNMs on T2T assemblies; CpG rate calibration. |
| Lucas JK ... 2026 bioRxiv, "HPRC2: A human pangenome reference with near-complete coverage of common genetic variation" | PMID 42539208; 10.64898/2026.07.21.739710 | Citable HPRC2 reference for the map-fitting plan. |
| Son DR ... Yi SV 2026 bioRxiv (above) | 10.64898/2026.06.15.732472 | Ape T2T + germline methylomes; CpG erosion vs hypomethylated reservoirs. |
| Tian Y ... 2026 HGG Adv, "Fine mapping regulatory variants by characterizing native CpG methylation with nanopore long-read sequencing" | PMID 41109954; 10.1016/j.xhgg.2025.100532 | Long-read allele-specific methylation for variant fine-mapping. |
| Spisak N, de Manuel M, Milligan W, Sella G, Przeworski M 2024 PLoS Biol, "The clock-like accumulation of germline and somatic mutations can arise from the interplay of DNA damage and repair" | PMID 38885262; 10.1371/journal.pbio.3002678 | Damage/repair model incl. methyl-CpG; mutation-rate theory. |
| de Manuel M, Przeworski M, Spisak N, Stolyarova A 2025/2026 PLoS Biol essay, "What sets the mutation rate of a cell type in an animal species?" | bioRxiv 10.64898/2025.12.19.695482; published 10.1371/journal.pbio.3003799 (per bioRxiv API; PubMed record not checked) | Perspective; context only. |
| Ma Z, Starr AL, Gokhman D, Fraser H 2026 bioRxiv (L) | 10.64898/2026.01.20.700710 (v updated 2026-09-29) | Cis/trans methylation divergence and selection; R7 target. |

## Unverified leads (not confirmed; do not cite yet)

- Matsuo K, Clay O, Takahashi T, Silke J, Schaffner W 1993 Gene, "Evidence for erosion of mouse CpG islands during
  mammalian evolution" (title from memory; not found in PubMed search).
- Jiang C, Han L, Su B, Li WH, Zhao Z 2007, "Features and trend of loss of promoter-associated CpG islands in the human
  and mouse genomes" (MBE?), not located.
- Arndt PF, Petrov DA, Hwa T 2003 MBE, "Distinct changes of genomic biases in nucleotide substitution at the time of
  mammalian radiation" (not retrieved in this session's author search).
- MethEvolSIM (Castillo Vicente S, Metzler D; CRAN v0.1.1, 2024): software verified on CRAN only, no paper found.
- A journal version of Hatcher 2021: none found (bioRxiv API reports no published version).
- Whether gnomAD v4 / Gnocchi's methylation stratification uses germline vs somatic methylation: not re-checked here
  (LITERATURE.md notes Karczewski 2020 used somatic).

## Suggested additions to md/LITERATURE.md (for the laptop session)

1. **Strand table, "Neutral CpG decay gated by methylation" row**: add Arndt, Burge & Hwa 2003 J Comput Biol (PMID
   12935330) and Mugal 2015 G3 (PMID 25591920) as the explicit substitution-model treatments. Add Bergman & Schierup 2021
   Genetics (PMID 34081117) and Berglund 2014 GBE (PMID 25527838) for gBGC at CpGs, which feed R2's W↔S control.
2. **"TF-binding-site population genetics" row**: add Behrens & Vingron 2010 (PMID 21128851) and Zemojtel 2011 (PMID
   22016335) as the only binding-site models with CpG-type mutation (neutral). Add Haldane 2014 (PMID 25010228) and
   Stewart 2012 (PMID 22887818). State explicitly: "no extension to methylation found (search 2026-10-03)".
3. **New row "Selection inferred at methylated CpGs"**: Agarwal & Przeworski 2021 eLife (PMID 34806592), Harpak 2016
   (PMID 27977673; recurrent mutation and SFS), Ying & Huttley 2011 (PMID 21398426), Si/Kang/Zöllner 2023 bioRxiv
   (circularity caveat), Cocozza 2011 (PMID 21829712; the confounded design F2 warns about).
4. **"Experimental sequence determinants" row**: add Long 2016 NAR (PMID 27084945; protection is DNA-encoded) and Hartl
   2019 Genome Res (PMID 30709850; CpG density acts independent of DNMTs). Add a "Limits" note that the model omits a
   CpG→function route that bypasses methylation.
5. **"Forward simulation" row**: add Haerter 2014 / Lövkvist 2016 NAR (PMIDs 24288373, 26932361) as mechanistic
   sequence→methylation maps with no evolution, a candidate replacement for the hand-set map. Add Grosser & Metzler 2020
   (PMID 32183713) as tree-based epimutation simulation. Keep "none found" for the combination.
6. **meQTL causality row**: add Gutierrez-Arcelus 2013 eLife (PMID 23755361), Do 2016 AJHG (PMID 27153397) and Onuchic
   2018 Science (PMID 30139913) as empirical "passenger" evidence. Add Costa 2025 Mol Ecol (PMID 39582237) as a primate
   meQTL dataset.
7. **Real-data additions**: Son 2026 bioRxiv (10.64898/2026.06.15.732472) for R1; Kocher 2024 Genome Biol (PMID 38872220)
   for R5; HPRC2 preprint Lucas 2026 (10.64898/2026.07.21.739710); Porubsky 2025 Nature (PMID 40269156); Gustafson 2024
   Genome Res (PMID 39358015); Chen 2024 Nature Gnocchi (PMID 38057664).
8. **Synthesis wording**: "Cohen 2011 + Mustonen & Lässig" remains the correct two-line fusion. Name Agarwal & Przeworski
   2021 as the frequency-side precedent and Behrens/Zemojtel as the only binding-site + CpG-mutation work.
