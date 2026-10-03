# realdata/ (cluster-owned)

Real-data arm of mcg_popgen, owned by the cluster session (`md/ENDPOINT.md` repo rules: the laptop never edits this
directory; the cluster works in this clone and pulls before committing).

- `md/REVIEW_2026-10-02.md`: math and logic review of the simulation phase (corrections to F2 and the real-N job)
- `md/DESIGN_PLAN.md`: real-data design (R0 data, R1 CpG-loss P/D, R2 P/D vs somatic methylation, R3 o/e +
  conservation, R4 HPRC2), compute and time estimates, open decisions

Large data (methylomes, VCFs, alignments) stays in `/u/project/cluo_scratch/terencew/claude/mcg_popgen/` and is never
committed. Record paths and md5s here as they are downloaded.
