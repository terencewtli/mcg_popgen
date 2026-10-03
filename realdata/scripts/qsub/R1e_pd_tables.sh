#!/bin/bash
#$ -cwd
#$ -o /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/logs/R1e_pd_tables.$JOB_ID
#$ -j y
#$ -N R1e_pd_tables
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4

# Usage: qsub [-v MODES='exclude_any exclude_strong only_strong'] R1e_pd_tables.sh
# R1e: genome-wide CpG-loss P/D tables (germline x element strata, coding positive control) from the R1d annotations.
set -euo pipefail
source /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/github/mcg_popgen/realdata/config/paths.sh
source ~/project-cluo/miniconda3/etc/profile.d/conda.sh
set +u; conda activate allcools; set -u
mkdir -p $repodir/tsv/R1
echo "start $(date)"
# MODES: all (no motif filter), and the R1g motif filters (exclude_any / exclude_strong / only_strong)
for mode in ${MODES:-all}; do
    echo "== motif mode: $mode"
    time python $repodir/scripts/R1e_pd_tables.py --annot-glob "$r1/annot/chr*.annot.tsv.gz" --outdir $repodir/tsv/R1 \
        --motif-dir $r1/motif --motif-mode $mode
done
echo "end $(date)"
