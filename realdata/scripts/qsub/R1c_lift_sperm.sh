#!/bin/bash
#$ -cwd
#$ -o /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/logs/R1c_lift_sperm.$JOB_ID.$TASK_ID
#$ -j y
#$ -N R1c_lift_sperm
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
#$ -tc 40

# R1c: human (hg18) and chimp (panTro2) sperm CpG methylation (Molaro 2011) lifted to hg38.
# Chimp goes through the human genome (panTro2 -> hg19 -> hg38), so each lifted chimp CpG sits on its human ortholog.
# Output bed (0-based, 1 bp, sorted): chrom start end meth cov. The lifted base may be the C or the G of the hg38
# dinucleotide (minus-strand chains); R1d resolves it against the sequence.
# Modes (qsub -v MODE=...):
#   prep   -t 1      build both input beds, lift human (3 min), split the chimp bed into 40 chunks
#   chimp  -t 1-40   lift one chimp chunk (single-pass panTro2 -> hg19 took ~9 h, so it is split)
#   merge  -t 1      concatenate and sort
# Usage: j1=$(qsub -terse -v MODE=prep -t 1 R1c_lift_sperm.sh | cut -d. -f1)
#        j2=$(qsub -terse -v MODE=chimp -t 1-40 -hold_jid $j1 R1c_lift_sperm.sh | cut -d. -f1)
#        qsub -v MODE=merge -t 1 -hold_jid $j2 R1c_lift_sperm.sh
set -euo pipefail
source /u/project/cluo/terencew/claude/project_ideas/mcg_popgen/github/mcg_popgen/realdata/config/paths.sh
out=$r1/sperm; mkdir -p $out/chunks
lo=$liftover_bin; ch=$dl/liftover; sp=$dl/sperm
MODE=${MODE:-prep}
# ID=1
ID=$SGE_TASK_ID
echo "$MODE $ID start $(date)"

mkbed() {  # meth.bedgraph.gz cov.bedgraph.gz -> bed with name = meth:cov
    paste <(zcat $1) <(zcat $2) | awk -F'\t' 'BEGIN{OFS="\t"} {if($1!=$5||$2!=$6) {print "coordinate mismatch at line "NR > "/dev/stderr"; exit 1} print $1,$2,$3,$4":"$8}'
}
fmt() {  # liftOver bed (name = meth:cov) -> chrom start end meth cov
    awk 'BEGIN{OFS="\t"} {split($4,a,":"); print $1,$2,$3,a[1],a[2]}'
}

if [ "$MODE" = prep ]; then
    if [ ! -s $out/human.hg38.raw.bed ]; then
        time mkbed $sp/GSE30340_human_sperm_CpG_methylation_hg18.bedgraph.gz $sp/GSE30340_human_sperm_CpG_coverage_hg18.bedgraph.gz > $out/human.hg18.bed
        time $lo -minMatch=0.9 $out/human.hg18.bed $ch/hg18ToHg38.over.chain.gz $out/human.hg38.raw.bed $out/human.unmapped.bed
    fi
    if [ ! -s $out/chimp.panTro2.bed ]; then
        time mkbed $sp/GSE30340_chimp_sperm_CpG_methylation_panTro2.bedgraph.gz $sp/GSE30340_chimp_sperm_CpG_coverage_panTro2.bedgraph.gz > $out/chimp.panTro2.bed
    fi
    time split -n l/40 -d -a 2 $out/chimp.panTro2.bed $out/chunks/chimp.panTro2.
    ls $out/chunks | head -3
elif [ "$MODE" = chimp ]; then
    k=$(printf '%02d' $((ID - 1)))
    c=$out/chunks/chimp.panTro2.$k
    time $lo -minMatch=0.5 $c $ch/panTro2ToHg19.over.chain.gz $c.hg19 $c.unmapped_hg19
    time $lo -minMatch=0.9 $c.hg19 $ch/hg19ToHg38.over.chain.gz $c.hg38.tmp $c.unmapped_hg38
    mv $c.hg38.tmp $c.hg38
    echo "chunk $k: $(wc -l < $c) in, $(wc -l < $c.hg38) lifted"
elif [ "$MODE" = merge ]; then
    n=$(ls $out/chunks/chimp.panTro2.??.hg38 | wc -l)
    [ "$n" -eq 40 ] || { echo "only $n of 40 chimp chunks lifted"; exit 1; }
    time (fmt < $out/human.hg38.raw.bed | sort -k1,1 -k2,2n -S 8G > $out/human.hg38.bed)
    time (cat $out/chunks/chimp.panTro2.??.hg38 | fmt | sort -k1,1 -k2,2n -S 8G > $out/chimp.hg38.bed)
    for s in human chimp; do echo "$s: $(wc -l < $out/$s.hg38.bed) lifted to hg38"; done
    echo "inputs: human $(wc -l < $out/human.hg18.bed), chimp $(wc -l < $out/chimp.panTro2.bed)"
fi
echo "$MODE $ID end $(date)"
