# Derived files

Files under `workflows/overrides/` that were copied from `vendor/` and modified,
and files under `workflows/ugc_wgw/assembly/` that were ported from PacBio's
HiFi-human-assembly-WDL at v1.0.2 (commit `025b3910`, wdl-common `e37b3274`;
vendored for the port in session 6 and dropped on 2026-09-24, since it was
never imported and upstream has not moved past that commit). This table is the
complete surface area of an upstream sync. Keep it short.

| path | origin path | origin tag | why | last reconciled at tag |
|---|---|---|---|---|
| `workflows/ugc_wgw/assembly/hifiasm.wdl` | HiFi-human-assembly-WDL `workflows/assemble_genome/assemble_genome.wdl` (hifiasm_assemble) | v1.0.2 | ported to v3 conventions; our hifiasm 0.25.0 image (`containers/hifiasm`); unitig globs widened to `*.[pr]_utg`; asserts one graph per haplotype | v1.0.2 |
| `workflows/ugc_wgw/assembly/gfatools.wdl` | `…/assemble_genome.wdl` (gfa2fa) | v1.0.2 | ported; gfatools 0.5_34e0fcf_build2 image | v1.0.2 |
| `workflows/ugc_wgw/assembly/align_hifiasm.wdl` | `…/assemble_genome.wdl` (align_hifiasm, paftools) | v1.0.2 | ported; `set -euo pipefail`, explicit haplotype input and output names, paftools 4 CPU / 32 GB with `sort -S --parallel` | v1.0.2 |
| `workflows/ugc_wgw/assembly/yak.wdl` | HiFi-human-assembly-WDL `workflows/de_novo_assembly_trio/de_novo_assembly_trio.wdl` (yak_count) | v1.0.2 | ported; yak 0.1_build3 image | v1.0.2 |
| `workflows/ugc_wgw/assembly/bcftools_zip_stats.wdl` | HiFi-human-assembly-WDL wdl-common `wdl/tasks/zip_index_vcf.wdl` + `bcftools_stats.wdl` | v1.0.2 | folded into one task on pb_wdl_base | v1.0.2 |
| `workflows/overrides/glnexus.wdl` | `vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/tasks/glnexus.wdl` | v4.0.0 | one shard per task: `regions_bed` required, `--list` input file instead of the gVCF symlinks v4 makes for cloud engines, `ulimit -Sn`, scratch on `$TMPDIR`, output named by shard; `threads`/`mem_gb = 60` inputs and the GLnexus config block as upstream | v4.0.0 |
