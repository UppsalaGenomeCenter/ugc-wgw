# `cohort_call`

Joint genotyping of the members' gVCFs with GLnexus, per region shard, and
each member's slice.

| | |
|---|---|
| Subject | cohort |
| Modes | joint |
| All stages | [08 Stages](../08-stages.md) |

GLnexus over the cohort's gVCFs, scattered by region; a per-sample split; an
optional multi-sample sawfish call. Mirrors upstream's `joint.wdl`.

```mermaid
flowchart TD
  regions[/"scatter_regions.GRCh38_GIABv3.tsv"/]:::data --> prep["ugc_wgw_regions_prepare: one BED per shard"]
  subgraph slice["per sample"]
    gv[/"gVCF"/]:::data --> sl["bcftools view -r, all shards in one task"]
  end
  prep --> sl
  subgraph shard["per shard, 26 with the shipped regions"]
    gl["glnexus_cli --bed, every sample's slice"]
  end
  sl --> gl
  prep --> gl
  gl --> concat["bcftools concat"] --> cvcf[/"cohort.joint.GRCh38_GIABv3.small_variants.vcf.gz"/]:::data
  cvcf --> split["bcftools +split per contig, then one concat per sample"] --> per[/"per-sample small-variant VCFs, cohort order"/]:::result
  subgraph sv["if run_sawfish_joint_call, default false"]
    tars[/"discover tars, aligned BAMs"/]:::data --> sc["sawfish call, multi-sample"] --> ssplit["split by sample"] --> psv[/"per-sample SV VCFs"/]:::result
  end
  classDef default fill:none,stroke:#8a8f95
  classDef data fill:#5d6875,stroke:#5d6875,color:#ffffff
  classDef result fill:#2f7d4a,stroke:#2f7d4a,color:#ffffff
  classDef box fill:#1f7a70,fill-opacity:0.10,stroke:#1f7a70,stroke-opacity:0.45
  class slice,shard,sv box
```

Why the slicing: `glnexus_cli` ignores gVCF indices and reads every gVCF end
to end for each region, so 26 shards over 1000 samples would read 26 000
whole gVCFs. One `bcftools view` per sample cuts all 26 slices index-based,
and each GLnexus shard reads only its slices.

| Driver fills | You may override |
|---|---|
| `cohort_id`, `sample_ids`, `gvcfs`, `gvcf_indices`; with `run_sawfish_joint_call`: `discover_tars`, `aligned_bams`, `aligned_bam_indices` | `split_keep_homref` (true), `run_sawfish_joint_call` (false), `glnexus_threads` (32), `glnexus_mem_gb` (32 + 0.1 N), `split_mem_gb` (8 + 0.02 N), `scatter_regions_file`, `skip_contigs` (`["chrM"]`) |

GLnexus skips the `chrM` shard (`skip_contigs`): the cohort's mitochondrial
genotypes come from mitorsaw through `cohort_merge`, not from DeepVariant's
diploid calls on a genome that is not diploid, so the cohort small-variant
VCF carries no chrM record since 0.15.

Resources per shard: GLnexus `glnexus_threads` CPUs and `glnexus_mem_gb` GB,
which is also its `--mem-gbytes` budget; the database lives on the node's
local scratch when the site set `TASK_SCRATCH_ENV` (chapter 04), else in
the task's work directory.
The split runs one bgzip writer per sample. Outputs: cohort VCF, per-sample
split VCFs in cohort order, optional joint SV set (chapter 07).
