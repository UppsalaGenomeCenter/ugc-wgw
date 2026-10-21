# Derived files

Files under `workflows/overrides/` that were copied from `vendor/` and modified,
and the files under `workflows/ugc_wgw/assembly/` that were ported from another
project: PacBio's HiFi-human-assembly-WDL at v1.0.2 (commit `025b3910`;
vendored for the port in session 6 and dropped on 2026-09-24, since it was
never imported and upstream has not moved past that commit) for `hifiasm.wdl`,
`gfatools.wdl`, `yak.wdl` and `align_asm.wdl`, and HMM-Flagger v1.2.0
(`mobinasri/flagger`) for `flagger.wdl`. The other files there
(`chrom_coverage.wdl`, `concordance.wdl`, `dipcall.wdl`, `masking.wdl`,
`pav.wdl`, `vcf_finish.wdl`) and everything under `workflows/ugc_wgw/pangenome/`
are ours (`origin: none`) and are not listed; the third-party images they run
(dipcall, PAV, HMM-Flagger, secphase, Minigraph-Cactus, Truvari) are pinned by
digest in `image_manifest.ugc-wgw.txt`. This table is the complete surface area
of an upstream sync. Keep it short.

| path | origin path | origin tag | why | last reconciled at tag |
|---|---|---|---|---|
| `workflows/ugc_wgw/assembly/hifiasm.wdl` | HiFi-human-assembly-WDL `workflows/assemble_genome/assemble_genome.wdl` (hifiasm_assemble) | v1.0.2 | ported to v3 conventions; our hifiasm 0.25.0 image (`containers/hifiasm`); unitig globs widened to `*.[pr]_utg`; asserts one graph per haplotype | v1.0.2 |
| `workflows/ugc_wgw/assembly/gfatools.wdl` | `…/assemble_genome.wdl` (gfa2fa) | v1.0.2 | ported; gfatools 0.5_34e0fcf_build2 image | v1.0.2 |
| `workflows/ugc_wgw/assembly/align_asm.wdl` | the ported `align_hifiasm.wdl` (HiFi-human-assembly-WDL `…/assemble_genome.wdl`, align_hifiasm) | v1.0.2 | rewritten for 0.11.0: our minimap2 2.31 image, one run to SAM, `samtools sort` and `paftools.js sam2paf` (BAM + PAF); paftools call retired with the stage's calling | v1.0.2 |
| `workflows/ugc_wgw/assembly/flagger.wdl` | HMM-Flagger v1.2.0 (mobinasri/flagger) `wdls/workflows/hmm_flagger_end_to_end_with_mapping.wdl`, `tasks/alignment/long_read_aligner.wdl`, `tasks/coverage/bam_coverage.wdl`, `tasks/hmm_flagger/hmm_flagger.wdl`, `ext/secphase` correct_bam | v1.2.0 | the commands ported into three tasks around the upstream images (minimap2 ours): read alignment, correct_bam, bam2cov + hmm_flagger with the whole-assembly annotation built as Flagger's task does; not vendored (Cromwell runtime keys, submodule imports) | v1.2.0 |
| `workflows/ugc_wgw/assembly/yak.wdl` | HiFi-human-assembly-WDL `workflows/de_novo_assembly_trio/de_novo_assembly_trio.wdl` (yak_count) | v1.0.2 | ported; yak 0.1_build3 image | v1.0.2 |
| `workflows/overrides/glnexus.wdl` | `vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/tasks/glnexus.wdl` | v4.0.0 | one shard per task: `regions_bed` required, `--list` input file instead of the gVCF symlinks v4 makes for cloud engines, `ulimit -Sn`, scratch on `$TMPDIR`, output named by shard; `threads`/`mem_gb = 60` inputs and the GLnexus config block as upstream | v4.0.0 |
| `workflows/overrides/pbmm2.wdl` | `vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/workflows/pbmm2/pbmm2.wdl` (workflow block only; the tasks are imported from the vendored file) | v4.0.0 | optional `pbmm2_index` input: the prebuilt reference index from the reference map; `create_pbmm2_index` runs only when it is absent | v4.0.0 |
| `workflows/overrides/upstream.wdl` | `vendor/hifi-human-wgs-wdl/workflows/upstream/upstream.wdl` | v4.0.0 | imports `overrides/pbmm2.wdl` for the HiFi alignment and the vendored pbmm2 tasks for the fail-reads path; optional `pbmm2_index` input passed to both, so the 5.8 GB index is not rebuilt once per input BAM; otherwise call for call upstream | v4.0.0 |

