#!/bin/bash
#$ -cwd
#$ -o /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/logs/R2b_soma.$JOB_ID.$TASK_ID
#$ -j y
#$ -N R2b_soma
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
#$ -tc 22

# R2b / R2c. MODE=chrom (-t 1-22): regional germline + somatic m for ancestral CpGs (R2b) and the matching non-CpG
# counts (R1f with --human-sperm / --loyfer). MODE=tables (-t 1): R2c genome-wide tables to tsv/R2.
# Usage: j=$(qsub -terse -v MODE=chrom -t 1-22 R2b_soma.sh | cut -d. -f1); qsub -v MODE=tables -t 1 -hold_jid $j R2b_soma.sh
set -euo pipefail
source /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/github/mcg_popgen/realdata/config/paths.sh
source ~/project-cluo/miniconda3/etc/profile.d/conda.sh
set +u; conda activate allcools; set -u
MODE=${MODE:-chrom}
r2=$scratch/R2
mkdir -p $r2/soma $r2/gbgc_strat $repodir/tsv/R2

# ID=1
ID=$SGE_TASK_ID
echo "$MODE $ID start $(date)"
if [ "$MODE" = chrom ]; then
    chrom=chr$ID
    time python $repodir/scripts/R2b_cpg_soma.py --chrom $chrom --annot $r1/annot/$chrom.annot.tsv.gz \
        --loyfer $r2/loyfer/$chrom.tsv.gz --out $r2/soma/$chrom.soma.tsv.gz.tmp
    mv $r2/soma/$chrom.soma.tsv.gz.tmp $r2/soma/$chrom.soma.tsv.gz
    time python $repodir/scripts/R1f_gbgc_counts.py --chrom $chrom --fasta $hg38_fa \
        --chimp-axt $pantro6_axt --macaque-axt $rhemac10_axt --snv $r1/snv/$chrom.snv.tsv.gz \
        --ccre $ccre_dir/GRCh38-cCREs.V3.bed --cgi $dl/ucsc/cpgIslandExt.txt.gz --gtf $gencode_gtf \
        --human-sperm $r1/sperm/human.hg38.bed --loyfer $r2/loyfer/$chrom.tsv.gz \
        --out $r2/gbgc_strat/$chrom.gbgc.tsv.gz.tmp
    mv $r2/gbgc_strat/$chrom.gbgc.tsv.gz.tmp $r2/gbgc_strat/$chrom.gbgc.tsv.gz
else
    n=$(ls $r2/soma/chr*.soma.tsv.gz $r2/gbgc_strat/chr*.gbgc.tsv.gz | wc -l)
    [ "$n" -eq 44 ] || { echo "only $n of 44 per-chromosome files"; exit 1; }
    time python $repodir/scripts/R2c_pd_soma.py --soma-glob "$r2/soma/chr*.soma.tsv.gz" --motif-dir $r1/motif \
        --sw-glob "$r2/gbgc_strat/chr*.gbgc.tsv.gz" --outdir $repodir/tsv/R2
fi
echo "$MODE $ID end $(date)"
