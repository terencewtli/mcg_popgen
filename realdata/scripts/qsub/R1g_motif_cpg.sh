#!/bin/bash
#$ -cwd
#$ -o /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/logs/R1g_motif_cpg.$JOB_ID.$TASK_ID
#$ -j y
#$ -N R1g_motif_cpg
#$ -l h_data=8G,h_rt=8:00:00
#$ -pe shared 4
#$ -t 1-22
#$ -tc 22

# R1g: HOMER known vertebrate motifs scanned on the ancestral-CpG-restored chromosome; per-CpG motif overlap.
# Half-site motifs are dropped (they cover most of the genome). Only the CpG overlaps are kept, never the full hit bed.
set -euo pipefail
source /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/github/mcg_popgen/realdata/config/paths.sh
source ~/project-cluo/miniconda3/etc/profile.d/conda.sh
set +u; conda activate allcools; set -u
export PATH=/u/project/cluo/terencew/programs/homer/bin:$HOME/bin:$PATH
motifs=/u/project/cluo/terencew/programs/homer/data/knownTFs/vertebrates/known.motifs

# ID=1
ID=$SGE_TASK_ID
chrom=chr$ID
tmp=$scratch/tmp/R1g_$chrom; mkdir -p $tmp $r1/motif
echo "$chrom start $(date)"
time python $repodir/scripts/R1g_motif_cpg.py fasta --chrom $chrom --anc $r1/anc/$chrom.anc.tsv.gz \
    --fasta $hg38_fa --out-fa $tmp/anc.fa --out-bed $tmp/cpg.bed
time (scanMotifGenomeWide.pl $motifs $tmp/anc.fa -bed -keepAll -p ${NSLOTS:-4} 2> $tmp/scan.log \
    | awk -F'\t' 'tolower($4) !~ /half/' \
    | bedtools intersect -a stdin -b $tmp/cpg.bed -wa -wb \
    | python $repodir/scripts/R1g_motif_cpg.py collapse --chrom $chrom --anc $r1/anc/$chrom.anc.tsv.gz \
        --motifs $motifs --out $r1/motif/$chrom.motif.tsv.gz.tmp)
mv $r1/motif/$chrom.motif.tsv.gz.tmp $r1/motif/$chrom.motif.tsv.gz
tail -n 3 $tmp/scan.log
echo "$chrom end $(date)"
