#!/bin/bash
#$ -cwd
#$ -o /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/logs/R1f_gbgc.$JOB_ID.$TASK_ID
#$ -j y
#$ -N R1f_gbgc
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
#$ -tc 22

# R1f: gBGC control. MODE=counts (-t 1-22): non-CpG P / D counts per chromosome; MODE=tables (-t 1): genome-wide P/D.
# Usage: j=$(qsub -terse -v MODE=counts -t 1-22 R1f_gbgc.sh | cut -d. -f1); qsub -v MODE=tables -t 1 -hold_jid $j R1f_gbgc.sh
set -euo pipefail
source /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/github/mcg_popgen/realdata/config/paths.sh
source ~/project-cluo/miniconda3/etc/profile.d/conda.sh
set +u; conda activate allcools; set -u
MODE=${MODE:-counts}

# ID=1
ID=$SGE_TASK_ID
mkdir -p $r1/gbgc $repodir/tsv/R1
echo "$MODE $ID start $(date)"
if [ "$MODE" = counts ]; then
    chrom=chr$ID
    time python $repodir/scripts/R1f_gbgc_counts.py --chrom $chrom --fasta $hg38_fa \
        --chimp-axt $pantro6_axt --macaque-axt $rhemac10_axt --snv $r1/snv/$chrom.snv.tsv.gz \
        --ccre $ccre_dir/GRCh38-cCREs.V3.bed --cgi $dl/ucsc/cpgIslandExt.txt.gz --gtf $gencode_gtf \
        --out $r1/gbgc/$chrom.gbgc.tsv.gz.tmp
    mv $r1/gbgc/$chrom.gbgc.tsv.gz.tmp $r1/gbgc/$chrom.gbgc.tsv.gz
else
    n=$(ls $r1/gbgc/chr*.gbgc.tsv.gz | wc -l)
    [ "$n" -eq 22 ] || { echo "only $n of 22 chromosomes"; exit 1; }
    time python $repodir/scripts/R1f_gbgc_tables.py --counts-glob "$r1/gbgc/chr*.gbgc.tsv.gz" --out $repodir/tsv/R1/R1f_gbgc_pd.tsv
fi
echo "$MODE $ID end $(date)"
