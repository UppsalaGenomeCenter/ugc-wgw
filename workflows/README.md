# workflows/

The nine entrypoints `ugc_wgw_<stage>.wdl` (singleton, upstream, cohort_call,
downstream, cohort_merge, cohort_freq, assembly, assembly_call,
cohort_pangenome) are derived from the real signatures in
`vendor/hifi-human-wgs-wdl/workflows/` at the tag in `upstream.lock`. Stage model:
`docs/DESIGN.md` §5–§6; contracts: `docs/ENTRYPOINTS.md`. `ugc_wgw_upstream.wdl`,
`ugc_wgw_downstream.wdl` and `ugc_wgw_singleton.wdl` forward every upstream output by name;
regenerate their output blocks from the vendored files on a sync rather than
editing them by hand.

- `ugc_wgw/cohort/` — ours: regions, concat, gvcf_slice, glnexus_scatter, split_by_sample, svx, sv_merge_bcftools, trgt_merge_lps, freq
- `ugc_wgw/assembly/` — ports (rows in `DERIVED_FILES.md`): hifiasm, gfatools and yak from HiFi-human-assembly-WDL v1.0.2; align_asm, that repository's align_hifiasm rewritten for our minimap2 2.31 image (the paftools call is retired); flagger, HMM-Flagger v1.2.0's commands as three tasks around the upstream images. Ours: chrom_coverage, dipcall, pav, vcf_finish, masking, concordance
- `ugc_wgw/pangenome/` — ours: cactus_split and cactus_chrom (Minigraph-Cactus v3.3.0's own commands), names, vcf_post, truvari
- `ugc_wgw/fail.wdl` — ugc_wgw_fail: a task that fails with a message (WDL 1.0 has no assert)
- `ugc_wgw/provenance.wdl` — ugc_wgw_manifest_write
- `overrides/` — glnexus, pbmm2, upstream: copied+modified upstream files; every file listed in `DERIVED_FILES.md`

Header for every authored file:

```wdl
# ugc-pacbio-wgw — <one-line purpose>
# origin: <none | vendor/<name>/<path> @ <tag>>
# see docs/DESIGN.md §<n>
```
