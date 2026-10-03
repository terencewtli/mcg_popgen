#!/bin/bash
#$ -cwd
#$ -o /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/logs/R1a_snv_context.$JOB_ID.$TASK_ID
#$ -j y
#$ -N R1a_snv_context
#$ -l h_data=8G,h_rt=8:00:00
#$ -pe shared 4
#$ -t 1-22
#$ -tc 22

# R1a: per-chromosome table of biallelic 1000G SNVs (2,504 unrelated) with hg38 context and CpG / W-S class.
set -euo pipefail
source /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/github/mcg_popgen/realdata/config/paths.sh
source ~/project-cluo/miniconda3/etc/profile.d/conda.sh
set +u; conda activate allcools; set -u
source /u/local/Modules/default/init/bash
module load bcftools

# ID=1
ID=$SGE_TASK_ID
chrom=chr$ID
vcf=$kg_vcf_dir/1000G.$chrom.vcf.gz
[ -s $dl/kg/1000G.$chrom.vcf.gz ] && vcf=$dl/kg/1000G.$chrom.vcf.gz   # re-downloaded copy (R0b; local chr1 is corrupt)
echo "vcf: $vcf"
fasta=/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosome.fa   # autosomes, chr-prefixed; same sequence as UCSC hg38
mkdir -p $r1/snv
out=$r1/snv/$chrom.snv.tsv.gz
echo "$chrom start $(date)"

fmt='%POS\t%REF\t%ALT\t%AC_AFR_unrel\t%AC_AMR_unrel\t%AC_EAS_unrel\t%AC_EUR_unrel\t%AC_SAS_unrel\t%AN_AFR_unrel\t%AN_AMR_unrel\t%AN_EAS_unrel\t%AN_EUR_unrel\t%AN_SAS_unrel\n'
time (bcftools view -G -v snps -m2 -M2 --threads 3 $vcf \
    | bcftools query -f "$fmt" \
    | python $repodir/scripts/R1a_snv_context.py --fasta $fasta --chrom $chrom --out $out.tmp)
mv $out.tmp $out

# per-chromosome summary: class x ws counts, and singletons
time zcat $out | awk -F'\t' 'NR>1{k=$8"\t"$9; n[k]++; if($4==1||$4==$5-1) s[k]++} END{for(k in n) print k"\t"n[k]"\t"s[k]+0}' \
    | sort > $r1/snv/$chrom.summary.tsv
cat $r1/snv/$chrom.summary.tsv
echo "$chrom end $(date)"
