#!/bin/bash
#$ -cwd
#$ -o /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/logs/R2d_imprinting.$JOB_ID
#$ -j y
#$ -N R2d_imprinting
#$ -l h_data=8G,h_rt=2:00:00
#$ -pe shared 4

# R2d: sperm-low / oocyte-high CpGs vs imprinted loci (geneimprint 'Imprinted' genes) and imprint-like somatic profiles.
set -euo pipefail
source /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/github/mcg_popgen/realdata/config/paths.sh
source ~/project-cluo/miniconda3/etc/profile.d/conda.sh
set +u; conda activate allcools; set -u
echo "start $(date)"
time python $repodir/scripts/R2d_imprinting.py --soma-glob "$scratch/R2/soma_oo/chr*.soma.tsv.gz" --motif-dir $r1/motif \
    --gtf $gencode_gtf --imprint $repodir/tsv/R2/geneimprint_human_2026-10-03.tsv --outdir $repodir/tsv/R2
echo "end $(date)"
