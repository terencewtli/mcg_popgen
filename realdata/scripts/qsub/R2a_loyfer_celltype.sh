#!/bin/bash
#$ -cwd
#$ -o /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/logs/R2a_loyfer_celltype.$JOB_ID.$TASK_ID
#$ -j y
#$ -N R2a_loyfer_celltype
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
#$ -t 1-22
#$ -tc 20

# R2a: Loyfer atlas betas -> per-CpG methylation by tissue-specific cell type (82 groups, 207 sorted samples), per chromosome.
set -euo pipefail
source /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/github/mcg_popgen/realdata/config/paths.sh
source ~/project-cluo/miniconda3/etc/profile.d/conda.sh
set +u; conda activate allcools; set -u

# ID=1
ID=$SGE_TASK_ID
chrom=chr$ID
mkdir -p $scratch/R2/loyfer
echo "$chrom start $(date)"
time python $repodir/scripts/R2a_loyfer_celltype.py --chrom $chrom --fasta $hg38_fa --beta-dir $dl/loyfer \
    --urls $repodir/tsv/R0/loyfer_hg38_beta_urls.tsv --matrix $dl/loyfer/GSE186458_series_matrix.txt.gz \
    --cgi $dl/ucsc/cpgIslandExt.txt.gz --out $scratch/R2/loyfer/$chrom
echo "$chrom end $(date)"
