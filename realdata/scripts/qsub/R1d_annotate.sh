#!/bin/bash
#$ -cwd
#$ -o /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/logs/R1d_annotate.$JOB_ID.$TASK_ID
#$ -j y
#$ -N R1d_annotate
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
#$ -t 1-22
#$ -tc 22

# R1d: annotate ancestral CpGs (R1b) with human / chimp sperm methylation (R1c), cCRE, CGI, phyloP and coding consequence.
set -euo pipefail
source /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/github/mcg_popgen/realdata/config/paths.sh
source ~/project-cluo/miniconda3/etc/profile.d/conda.sh
set +u; conda activate allcools; set -u

# ID=1
ID=$SGE_TASK_ID
chrom=chr$ID
mkdir -p $r1/annot
echo "$chrom start $(date)"
time python $repodir/scripts/R1d_annotate.py --chrom $chrom --fasta $hg38_fa \
    --anc $r1/anc/$chrom.anc.tsv.gz \
    --human-sperm $r1/sperm/human.hg38.bed --chimp-sperm $r1/sperm/chimp.hg38.bed \
    --ccre $ccre_dir/GRCh38-cCREs.V3.bed --cgi $dl/ucsc/cpgIslandExt.txt.gz \
    --phylop $dl/ucsc/hg38.phyloP100way.bw --gtf $gencode_gtf \
    --out $r1/annot/$chrom.annot.tsv.gz.tmp
mv $r1/annot/$chrom.annot.tsv.gz.tmp $r1/annot/$chrom.annot.tsv.gz
echo "$chrom end $(date)"
