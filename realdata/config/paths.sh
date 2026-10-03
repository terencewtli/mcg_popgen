# Paths for the real-data arm. Source this from every script; nothing else hardcodes paths.
projdir=/u/project/cluo/terencew/claude/project_ideas/mcg_popgen
repodir=$projdir/github/mcg_popgen/realdata
scratch=/u/project/cluo_scratch/terencew/claude/mcg_popgen
logdir=$projdir/logs

# R0 downloads (large, scratch only)
dl=$scratch/R0
tools=$scratch/tools
liftover_bin=/u/project/cluo/terencew/miniconda3/envs/babel/bin/liftOver   # the UCSC binary in $tools needs glibc >= 2.25 (CentOS 7 has 2.17)

# On disk already
kg_vcf_dir=/u/project/cluo/terencew/demux_benchmark/pool_design/vcf/1000G/by_chrom   # 1000G.chr{1..22}.vcf.gz, NYGC 30x, 3,202 samples
pantro6_axt=/u/project/cluo/terencew/claude/project_ideas/mcg_evo/reference/ucsc/hg38.panTro6.rbest.axt.gz
ccre_dir=/u/project/cluo/terencew/claude/project_ideas/mcg_evo/reference/encode_ccre
gencode_gtf=/u/project/cluo/terencew/reference/hg38_igvf/gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz

# Downloaded in R0a (UCSC hg38, chr-prefixed, all chromosomes; matches the wgbstools CpG order of the Loyfer betas)
hg38_fa=$dl/ucsc/hg38.fa
rhemac10_axt=$dl/ucsc/hg38.rheMac10.rbest.axt.gz

# R1 outputs
r1=$scratch/R1
