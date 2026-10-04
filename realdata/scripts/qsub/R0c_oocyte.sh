#!/bin/bash
#$ -cwd
#$ -o /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/logs/R0c_oocyte.$JOB_ID
#$ -j y
#$ -N R0c_oocyte
#$ -l h_data=8G,h_rt=24:00:00
#$ -pe shared 4

# R0c: human oocyte (MII 36 cells, GV 8) and ICM (19) single-cell WGBS from Zhu 2018 Nat Genet (GSE81233; hg19,
# bismark, all C contexts per cell). Download, keep CpGs, merge strands onto the C, pool cells per group, lift to hg38.
# Output: $r1/oocyte/<group>.hg38.bed (chrom start end m cov; same format as R1c, so R2b / R1f can reuse it).
# Re-runnable: finished downloads and pooled groups are skipped.
set -uo pipefail
source /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/github/mcg_popgen/realdata/config/paths.sh
urls=$repodir/tsv/R0/zhu2018_oocyte_icm_urls.tsv   # group gsm file size url
d=$dl/zhu2018; out=$r1/oocyte; mkdir -p $d $out
echo "start $(date)"

echo "== download ($(wc -l < $urls) files, ~32 GB)"
time (cut -f5 $urls | xargs -P 4 -I{} sh -c 'f='$d'/$(basename {}); [ -s $f.done ] || (wget -q -c -O $f {} && echo ok > $f.done) || echo "FAILED {}"')
awk -F'\t' -v d=$d '{cmd="stat -c %s "d"/"$3" 2>/dev/null"; cmd | getline s; close(cmd); if (s!=$4) {bad++; print "size mismatch", $3, s, $4}} END{print "size mismatches:", bad+0}' $urls

for g in MII GV ICM; do
    [ -s $out/$g.hg38.bed ] && { echo "skip $g"; continue; }
    echo "== pool $g ($(awk -v g=$g '$1==g' $urls | wc -l) cells)"
    # CpG rows only; minus-strand calls (the G of the CpG, 1-based pos) moved to the C: 0-based C = pos - 2 for '-', pos - 1 for '+'
    time (awk -v g=$g -F'\t' '$1==g{print $3}' $urls | while read f; do zcat $d/$f; done \
        | awk -F'\t' 'BEGIN{OFS="\t"} $10=="CpG"{p=($4=="-") ? $2-2 : $2-1; print $1,p,$6,$5}' \
        | sort -k1,1 -k2,2n -S 16G -T $scratch/tmp \
        | awk -F'\t' 'BEGIN{OFS="\t"} {k=$1"\t"$2; if (k!=pk && NR>1) print pc,pp,pp+1,pm":"pt; if (k!=pk) {pm=0; pt=0} pm+=$3; pt+=$4; pk=k; pc=$1; pp=$2} END{print pc,pp,pp+1,pm":"pt}' \
        > $out/$g.hg19.counts.bed)
    # name = met:total -> liftOver keeps it; then m = met / total
    time $liftover_bin -minMatch=0.9 $out/$g.hg19.counts.bed $dl/liftover/hg19ToHg38.over.chain.gz $out/$g.hg38.raw.bed $out/$g.unmapped.bed
    awk 'BEGIN{OFS="\t"} {split($4,a,":"); if (a[2]>0) print $1,$2,$3,a[1]/a[2],a[2]}' $out/$g.hg38.raw.bed | sort -k1,1 -k2,2n -S 8G > $out/$g.hg38.bed
    echo "$g: $(wc -l < $out/$g.hg19.counts.bed) CpGs hg19, $(wc -l < $out/$g.hg38.bed) on hg38; cov>=5: $(awk '$5>=5' $out/$g.hg38.bed | wc -l)"
done
echo "end $(date)"
