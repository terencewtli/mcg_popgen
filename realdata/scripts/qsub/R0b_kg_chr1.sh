#!/bin/bash
#$ -cwd
#$ -o /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/logs/R0b_kg_chr1.$JOB_ID
#$ -j y
#$ -N R0b_kg_chr1
#$ -l h_data=8G,h_rt=8:00:00
#$ -pe shared 4

# R0b: the on-disk 1000G chr1 VCF (demux_benchmark copy) fails gzip -t (BGZF error at chr1:190.67 Mb).
# Re-download the same NYGC 2022 phased panel file into scratch, test it, and compare with the local copy.
set -uo pipefail
source /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/github/mcg_popgen/realdata/config/paths.sh
mkdir -p $dl/kg
b=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV
f=1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
echo "start $(date)"
time wget -q -c -O $dl/kg/$f $b/$f
time wget -q -c -O $dl/kg/$f.tbi $b/$f.tbi
echo "== gzip -t on the fresh download"
time gzip -t $dl/kg/$f && echo "fresh: OK" || echo "fresh: CORRUPT (upstream file is broken too)"
echo "== byte comparison with the local copy (first differences)"
time cmp -l $kg_vcf_dir/1000G.chr1.vcf.gz $dl/kg/$f | head -5
echo "differing bytes: $(cmp -l $kg_vcf_dir/1000G.chr1.vcf.gz $dl/kg/$f | wc -l)"
ln -sf $dl/kg/$f $dl/kg/1000G.chr1.vcf.gz
ln -sf $dl/kg/$f.tbi $dl/kg/1000G.chr1.vcf.gz.tbi
md5sum $dl/kg/$f
echo "end $(date)"
