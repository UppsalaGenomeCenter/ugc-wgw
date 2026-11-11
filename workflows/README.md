# workflows/

The ten entrypoints `ugc_wgw_<stage>.wdl` (singleton, upstream, cohort_qc,
cohort_call, downstream, cohort_merge, cohort_freq, assembly, assembly_call,
cohort_pangenome) are derived from the real signatures in
`vendor/hifi-human-wgs-wdl/workflows/` at the tag in `upstream.lock`. Stage model:
`docs/DESIGN.md` §5–§6; contracts: `docs/ENTRYPOINTS.md`. `ugc_wgw_upstream.wdl`,
`ugc_wgw_downstream.wdl` and `ugc_wgw_singleton.wdl` forward every upstream output by name;
regenerate their output blocks from the vendored files on a sync rather than
editing them by hand.

- `ugc_wgw/cohort/` — ours: regions, concat, gvcf_slice, glnexus_scatter, split_by_sample, svx, sv_merge_bcftools, trgt_merge_lps, freq
- `ugc_wgw/qc/` — ours: somalier (the extract, the one-sample relate and the sex resolution the upstream override calls; ped, relate, ancestry and the summary tables of `ugc_wgw_cohort_qc`)
- `ugc_wgw/mito/` — ours: collapse (mitorsaw's per-haplotype chrM genotypes to diploid, the haplogroups with haplogrep3; called by the upstream override and, for members from before 0.15, by `ugc_wgw_cohort_merge`), merge (the cohort's mt sets and the gathered haplogroups)
- `ugc_wgw/methyl/` — ours: report (`methbat report` per sample, called by `ugc_wgw_singleton` and `ugc_wgw_downstream`, and by `ugc_wgw_cohort_methyl` for members whose run predates it), cohort (`methbat build`, the scoring against the background with the gene annotation and the events, `compare`, `joint-segment`, `deconvolve`, the summary tables)
- `ugc_wgw/tr/` — ours: cohort (`ugc_wgw_cohort_tr`'s tasks: the per-shard allele-length statistics and expansion-outlier scan over the merged TRGT VCF, the screen of the named disease loci against STRchive's thresholds, the gather)
- `ugc_wgw/intake/` — ours: checks (the input checks of an aligned BAM and pre-called VCFs against the installed reference, the VCF sample column renamed; called by `ugc_wgw_intake` and, for an aligned BAM given instead of reads, by the upstream override)
- `ugc_wgw/assembly/` — ports (rows in `DERIVED_FILES.md`): hifiasm, gfatools and yak from HiFi-human-assembly-WDL v1.0.2; align_asm, that repository's align_hifiasm rewritten for our minimap2 2.31 image (the paftools call is retired); flagger, HMM-Flagger v1.2.0's commands as three tasks around the upstream images. Ours: chrom_coverage, dipcall, pav, vcf_finish, masking, concordance
- `ugc_wgw/pangenome/` — ours: cactus_split and cactus_chrom (Minigraph-Cactus v3.3.0's own commands), names, vcf_post, truvari
- `ugc_wgw/fail.wdl` — ugc_wgw_fail: a task that fails with a message (WDL 1.0 has no assert)
- `ugc_wgw/provenance.wdl` — ugc_wgw_manifest_write
- `overrides/` — glnexus, pbmm2, sawfish (the image as a task input; 2.3.0_build1 by default since 0.17.0), upstream: copied+modified upstream files; every file listed in `DERIVED_FILES.md`

Header for every authored file:

```wdl
# ugc-pacbio-wgw — <one-line purpose>
# origin: <none | vendor/<name>/<path> @ <tag>>
# see docs/DESIGN.md §<n>
```
