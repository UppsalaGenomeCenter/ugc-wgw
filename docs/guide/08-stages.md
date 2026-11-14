# Stages

For each stage: what it runs, which inputs the driver fills and which you may
override through `stage_inputs`, the switches, the resources, and where the
outputs are described. The diagrams are simplified; the exact call graphs are
in [workflow-graphs.md](workflow-graphs.md). Upstream documents its tools in
`vendor/hifi-human-wgs-wdl/docs/` (`pbmm2.md`, `deepvariant.md`, `trgt.md`,
`kivvi.md`, `methbat.md`, `gpu.md`).

Every stage takes `ugc_wgw_version` (echoed back and checked), `ref_map_file`
(the one rendered map, chapter 04), `backend` (`HPC`) and `preemptible` from
the driver; `cohort_merge`, `assembly` and `assembly_call` also take
`ugc_wgw_container_registry` for the images this repository builds.
`ugc-wgw stage-inputs --stage <stage>` prints these tables from the
installed WDL: every input with its type, default, description, who fills it
and the project's override; `--nested` adds the inputs of the tasks inside
(chapter 05).

Per-sample resource requests are upstream's own and unchanged: every task's
request is in the generated [inventory](task-resources.md), and chapter 12
shows how a site caps or overrides them. The numbers below are the ones
ugc-pacbio-wgw sets. The N-scaled memory formulas are first estimates until
the cohort smoke test replaces them.

## The stages

One page per stage, in the order the modes run them:

| Stage | Subject | Modes | What it does |
|---|---|---|---|
| [`singleton`](stages/singleton.md) | per sample | standalone | Everything for one sample in one run: alignment, small variants, SVs, phasing, repeats, methylation, the targeted callers. |
| [`upstream`](stages/upstream.md) | per sample | joint | The first half of a sample in joint mode: alignment, the small-variant gVCF, SV discovery, the targeted callers. |
| [`cohort_qc`](stages/cohort_qc.md) | cohort, optional | standalone, joint | Identity: relatedness of every pair, the pedigree, sex and ancestry with somalier. |
| [`cohort_call`](stages/cohort_call.md) | cohort | joint | Joint genotyping of the members' gVCFs with GLnexus, per region shard, and each member's slice. |
| [`downstream`](stages/downstream.md) | per sample | joint | The second half of a sample in joint mode: phasing, repeats, methylation, pharmacogenes on the cohort-consistent calls. |
| [`intake`](stages/intake.md) | per sample | standalone, joint | For a sample that arrives with its own alignment and calls: the checks and what phasing needs. |
| [`cohort_merge`](stages/cohort_merge.md) | cohort | standalone, joint | The cohort's SV, tandem-repeat and mitochondrial call sets, and in standalone mode the joint small-variant set. |
| [`cohort_freq`](stages/cohort_freq.md) | cohort | standalone, joint | Allele frequencies of the cohort's call sets with sex strata. |
| [`cohort_tr`](stages/cohort_tr.md) | cohort, optional | standalone, joint | Tandem repeats across the cohort: the disease-locus screen, expansion outliers, the allele-length resource. |
| [`cohort_methyl`](stages/cohort_methyl.md) | cohort | standalone, joint | Methylation across the cohort: the background, each member's events, comparisons, imprinting, cell-type composition. |
| [`cohort_mqtl`](stages/cohort_mqtl.md) | cohort, optional | standalone, joint | Structural variants and tandem repeats tested against the methylation of nearby CpG islands; the island-by-member methylation matrix. |
| [`assembly`](stages/assembly.md) | per sample | assembly | A phased de novo assembly with hifiasm and its quality control. |
| [`assembly_call`](stages/assembly_call.md) | per sample | assembly | Variants from the assembly against the reference with dipcall and PAV. |
| [`cohort_pangenome`](stages/cohort_pangenome.md) | cohort | assembly | A Minigraph-Cactus pangenome graph of the members' assemblies and its variants. |

```
standalone:  singleton -> cohort_qc -> cohort_merge -> cohort_freq -> cohort_tr -> cohort_methyl -> cohort_mqtl
joint:       upstream -> cohort_qc -> cohort_call -> downstream -> cohort_merge -> cohort_freq -> cohort_tr -> cohort_methyl -> cohort_mqtl
assembly:    assembly -> assembly_call -> cohort_pangenome
```

`cohort_qc`, `cohort_tr` and `cohort_mqtl` are optional: no stage waits
for them. A sample that enters as an aligned BAM starts at `upstream`
without alignment; one that enters with its calls runs `intake` in the
place of `singleton` or `upstream` (chapter 05).

Each page has the same parts: what the stage runs (a simplified diagram
and its text version), the inputs the driver fills and the ones you may
override, the resources, the outputs, and for the cohort analyses how to
read the result. Chapter 14 walks through the pages those results end up
on, chapter 13 lists the QC flags.
