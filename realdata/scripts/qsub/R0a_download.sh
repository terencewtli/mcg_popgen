#!/bin/bash
#$ -cwd
#$ -o /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/logs/R0a_download.$JOB_ID
#$ -j y
#$ -N R0a_download
#$ -l h_data=8G,h_rt=24:00:00
#$ -pe shared 4

# R0a: public downloads for the real-data arm (DESIGN_PLAN R0). Re-runnable: wget -c resumes, finished files are skipped.
set -uo pipefail
source /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/github/mcg_popgen/realdata/config/paths.sh
mkdir -p $dl/{ucsc,liftover,sperm,loyfer,jaspar,gnomad} $tools
echo "start $(date)"

get() {  # url outdir
    local f=$2/$(basename $1)
    if [ -s $f.done ]; then echo "skip $f"; return; fi
    wget -q -c -O $f $1 && echo ok > $f.done && echo "got $f $(stat -c %s $f)" || echo "FAILED $1"
}

ucsc=https://hgdownload.soe.ucsc.edu

echo "== tools"
get $ucsc/admin/exe/linux.x86_64/liftOver $tools; chmod +x $tools/liftOver

echo "== UCSC hg38"
get $ucsc/goldenPath/hg38/bigZips/hg38.fa.gz $dl/ucsc
wget -q -O $dl/ucsc/hg38.bigZips.md5sum.txt $ucsc/goldenPath/hg38/bigZips/md5sum.txt
get $ucsc/goldenPath/hg38/vsRheMac10/reciprocalBest/axtRBestNet/hg38.rheMac10.rbest.axt.gz $dl/ucsc
wget -q -O $dl/ucsc/rheMac10.axt.md5sum.txt $ucsc/goldenPath/hg38/vsRheMac10/reciprocalBest/axtRBestNet/md5sum.txt
get $ucsc/goldenPath/hg38/database/cpgIslandExt.txt.gz $dl/ucsc
get $ucsc/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw $dl/ucsc

echo "== liftOver chains (human sperm hg18 -> hg38; chimp sperm panTro2 -> hg19 -> hg38)"
get $ucsc/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz $dl/liftover
get $ucsc/goldenPath/panTro2/liftOver/panTro2ToHg19.over.chain.gz $dl/liftover
get $ucsc/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz $dl/liftover

echo "== sperm WGBS, Molaro 2011 (GSE30340)"
geo=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30340/suppl
for f in human_sperm_CpG_methylation_hg18.bedgraph human_sperm_CpG_coverage_hg18.bedgraph human_sperm_hmr_hg18.bed \
         chimp_sperm_CpG_methylation_panTro2.bedgraph chimp_sperm_CpG_coverage_panTro2.bedgraph chimp_sperm_hmr_panTro2.bed; do
    get $geo/GSE30340_$f.gz $dl/sperm
done
get https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30340/matrix/GSE30340-GPL9115_series_matrix.txt.gz $dl/sperm
get https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30340/matrix/GSE30340-GPL9378_series_matrix.txt.gz $dl/sperm

echo "== JASPAR 2024, gnomAD v4.1 constraint"
get https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt $dl/jaspar
get https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv $dl/gnomad

echo "== Loyfer 2023 atlas (GSE186458), 253 hg38 beta files, ~15 GB"
get https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186458/matrix/GSE186458_series_matrix.txt.gz $dl/loyfer
export -f get; export dl
time cut -f4 $repodir/tsv/R0/loyfer_hg38_beta_urls.tsv | xargs -P 4 -I{} bash -c 'get {} $dl/loyfer'

echo "== checks"
# Loyfer: every file has the GEO size
awk -F'\t' -v d=$dl/loyfer '{cmd="stat -c %s "d"/"$2" 2>/dev/null"; cmd | getline s; close(cmd); if (s!=$3) {bad++; print "size mismatch", $2, s, $3}} END{print "loyfer size mismatches:", bad+0}' $repodir/tsv/R0/loyfer_hg38_beta_urls.tsv
# UCSC md5s
(cd $dl/ucsc && grep ' hg38.fa.gz$' hg38.bigZips.md5sum.txt | md5sum -c -)
(cd $dl/ucsc && grep 'hg38.rheMac10.rbest.axt.gz' rheMac10.axt.md5sum.txt | md5sum -c -)

echo "== hg38 fasta: unzip + index"
source /u/local/Modules/default/init/bash
module load samtools
if [ ! -s $hg38_fa.fai ]; then
    time (zcat $dl/ucsc/hg38.fa.gz > $hg38_fa && samtools faidx $hg38_fa)
fi

echo "== md5 manifest"
time (cd $dl && find . -type f ! -name '*.done' ! -name '*.fai' ! -name 'hg38.fa' | sort | xargs md5sum > $dl/MD5SUMS.txt)
du -sh $dl/*
echo "end $(date)"
