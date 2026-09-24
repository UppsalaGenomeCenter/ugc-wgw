# workflows/

Entrypoints `ugc_wgw_<stage>.wdl` are derived from the real signatures in
`vendor/hifi-human-wgs-wdl/workflows/` at the tag in `upstream.lock`. Stage model:
`docs/DESIGN.md` §5–§6; contracts: `docs/ENTRYPOINTS.md`. `ugc_wgw_upstream.wdl`,
`ugc_wgw_downstream.wdl` and `ugc_wgw_singleton.wdl` forward every upstream output by name;
regenerate their output blocks from the vendored files on a sync rather than
editing them by hand.

- `ugc_wgw/cohort/` — regions, concat, gvcf_slice, glnexus_scatter, split_by_sample, svx, sv_merge_bcftools, trgt_merge_lps, freq
- `ugc_wgw/assembly/` — hifiasm, gfatools, align_hifiasm (minimap2 + paftools), yak, bcftools_zip_stats: ported from the assembly repo into v3 conventions (rows in `DERIVED_FILES.md`)
- `ugc_wgw/provenance.wdl` — ugc_wgw_manifest_write
- `overrides/` — copied+modified upstream tasks; every file listed in `DERIVED_FILES.md`

Header for every authored file:

```wdl
# ugc-pacbio-wgw — <one-line purpose>
# origin: <none | vendor/<name>/<path> @ <tag>>
# see docs/DESIGN.md §<n>
```
