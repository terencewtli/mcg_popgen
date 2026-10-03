#!/bin/bash
#$ -cwd
#$ -o /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/logs/R1b_ancestral_cpg.$JOB_ID.$TASK_ID
#$ -j y
#$ -N R1b_ancestral_cpg
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
#$ -t 1-22
#$ -tc 22

# R1b: ancestral CpGs (CG in chimp and macaque) per hg38 chromosome, human-lineage fate, and 1000G polymorphism.
set -euo pipefail
source /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/github/mcg_popgen/realdata/config/paths.sh
source ~/project-cluo/miniconda3/etc/profile.d/conda.sh
set +u; conda activate allcools; set -u

# ID=1
ID=$SGE_TASK_ID
chrom=chr$ID
mkdir -p $r1/anc
echo "$chrom start $(date)"
time python $repodir/scripts/R1b_ancestral_cpg.py --fasta $hg38_fa --chrom $chrom \
    --chimp-axt $pantro6_axt --macaque-axt $rhemac10_axt \
    --snv $r1/snv/$chrom.snv.tsv.gz --out $r1/anc/$chrom.anc.tsv.gz.tmp
mv $r1/anc/$chrom.anc.tsv.gz.tmp $r1/anc/$chrom.anc.tsv.gz
echo "$chrom end $(date)"
