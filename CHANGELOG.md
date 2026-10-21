# Changelog

Release notes of the public copy. Each version corresponds to a verified
offline bundle; the upstream tag it is built on is named first.

## 0.12.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.5.0,
unchanged. The release is about the assembly track's reporting.

- **The analysis summaries read every assembly-track output.** The sample
  page's Assembly section shows contiguity with the NGx curve and hifiasm's
  low-quality regions, Flagger's composition overall and per haplotype, a
  sequence QC table (GC, N bases, telomeric contig ends) with the contigs
  carrying the most non-Hap bases, trio QC with the contigs carrying the
  most switches, the chromosomes (assembled fraction per haplotype and a
  track of `reliable.bed` against `dip.bed`), the karyotype's evidence, both
  callers with Ts/Tv, het/hom, spectra, the PASS-versus-masked split and
  PAV's structural variants by size, the reliable region and the concordance
  in full. The cohort page draws the assembly metrics as distributions with
  outliers and threshold failures, lists one row per member, counts
  karyotypes, and for the pangenome shows the reference length and
  haplotype sequence per chromosome, the graph's structure (reference gaps,
  node degrees, nesting levels, vcfwave types and allele lengths), truvari's
  genotype concordance and the members' genotypes in the graph. Digest
  schema 2: every per-sample digest rebuilds once.
- New QC: `asm_concordance_low` fails a sample whose assembly disagrees with
  its read-based SNVs (`asm_genotype_concordance_min`, 0.90: a sample swap
  scores 0.5–0.7), `asm_karyotype_unknown` warns, and
  `pg_member_call_rate_min` (0.90) lists members with little assembled
  sequence in the graph on the cohort page.
- New small tasks on existing outputs, none of the expensive tasks touched:
  `assembly` writes a per-contig QC table (`contig_qc`: GC, N bases and
  runs, telomeric bases at each end, Flagger's bases per label; seqtk in the
  Flagger image), `assembly_call` the callers' spectra (`call_spectra`),
  `cohort_pangenome` the members' genotype statistics (`member_stats`) and,
  per chromosome, the join's path statistics and reference gaps, the bases
  per genome and the node degree histogram; its frequency summary counts
  records by nesting level, vcfwave type and allele length.
- The image manifest drops the retired `align_hifiasm` image (32 images in
  the bundle); the smoke harness drives `assembly_call` and `cohort_pangenome`
  as single stages; the masking task's pipeline passes shellcheck 0.9.0 (the
  0.11.0 tag's CI failure).

## 0.11.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.5.0,
unchanged. There was no 0.10.0. This release is the assembly track: the
`assembly` mode becomes `assembly → assembly_call → cohort_pangenome`. The tree
at tag v0.11.0 fails its own CI on a shellcheck 0.9.0 finding in a task
command; public commit 50a73a0 on main carries the one-line fix.

- **Breaking, `assembly`.** The stage ends at an output contract
  (`assembler_used`, `assembler_version`, `haplotype_labels`, `hap_fastas`,
  statistics, graphs); `haplotypes`, `zipped_assembly_fastas`, `asm_bams`,
  `merged_asm_bam` and the `paftools_*` outputs are gone; the paftools call is
  gone with them, while the `align_hifiasm` image, now unused, still ships in
  this release's manifest.
- **Expensive.** hifiasm runs with `--telo-m CCCTAA` (`hifiasm_telomere_motif`,
  empty = off) and writes its version file, so every sample's assembly is
  recomputed once on this version; no later release touches the hifiasm
  command. HMM-Flagger v1.2.0 (own tasks around the published images) and,
  for trio-binned samples, `yak trioeval` give per-assembly QC.
- New stage `assembly_call` (per sample): minimap2 2.31 aligns each
  haplotype to a decoy-free GRCh38 (BAM and PAF), the stage reports assembled
  bases per chromosome, infers the karyotype, and calls with dipcall v0.3 and
  PAV 3. Primary VCFs carry the sample's name, unphased genotypes unless
  trio-binned, and `FILTER=ASM_UNRELIABLE` outside the Flagger-reliable
  region (added, never dropped); with a read-based VCF a concordance table.
- New stage `cohort_pangenome` (per cohort): a Minigraph-Cactus v3.3.0 graph
  of the members' assemblies, built per chromosome behind one genome-wide
  split (`chromosomes` under `stage_inputs` limits the build); the cohort
  VCFs (raw and vcfwave) under the sample IDs with non-trio members
  unphased, allele frequencies with karyotype strata, GBZ and GFA per
  chromosome, optional `truvari bench` against the cohort's read-based SV
  merge (`pangenome_sv_cohort`, a profile key). The graph's paths of members
  that were not trio-binned are not haplotypes: no LD, imputation or
  PanGenie-style use without statistical phasing (guide chapter 02).
- References: the `ugc-wgw-asm` 0.1.0 lock entry (the decoy-free contig list,
  the chromosome list, dipcall's PAR BED; the installer derives the FASTA and
  index from the data container); the map gains `asm_fasta`,
  `asm_fasta_index`, `asm_par_bed`, `asm_graph_chromosomes`.
- Images: a new `minimap2` image of ours (minimap2 2.31, samtools 1.24,
  paftools.js) and the published dipcall, HMM-Flagger, secphase, PAV,
  Minigraph-Cactus and Truvari images pinned by digest.
- Summary: assembly QC thresholds and flags, an "Assembly QC" block per
  sample, QC columns in the cohort table, `summary --pass-list` for freezing
  a QC-gated pangenome cohort, and `summary --mode assembly` with the cohort
  page's pangenome section.

## 0.9.4

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.5.0,
unchanged. The only workflow change is a parameter description.

- The `parabricks` profile asks for 4 GPUs of any type (`gpu_type` empty)
  and gives `run_parabricks_deepvariant` 256 GiB through the `mem_gb` stage
  input of `singleton` and `upstream`, after measuring the GPU shapes on
  five samples: four GPUs per job with an empty type use every GPU node,
  six GPUs take as long as four, three 2-GPU jobs give a 6-GPU node its
  best throughput, and the calls do not depend on the GPU count or model.
  Site copies of the profile under `<prefix>/profiles/` are never
  overwritten by the installer and need the same edit. Chapter 12 and the
  HPC checklist carry the numbers.
- The `gpuType` descriptions of `ugc_wgw_singleton` and `ugc_wgw_upstream`
  say the plugin applies the type at submit; they still described the
  behaviour before 0.9.3.

## 0.9.3

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.5.0. No
workflow changed since 0.9.1.

- The GPU type is applied by the plugin at submit instead of being passed
  as the workflow's `gpuType` input. That input reached every task through
  `RuntimeAttributes`, so the type was part of every call-cache key and
  moving a project between GPU node types re-aligned every sample. The
  driver now hands the project's `gpu_type` to the plugin per run
  (`MINIWDL__UGC_WGW__GPU_TYPE`, miniwdl's override of `[ugc_wgw]
  gpu_type`); an empty type requests any GPU of the partition, so one
  project can use nodes of different GPU types; alignments are shared across
  flavours. The run manifest records `gpu_type`; `--dry-run` prints the
  environment line. Projects created before this release align once more.
- `parabricks_gpus` must be even: upstream's task runs `pbrun deepvariant
  --run-partition`, which refuses an odd GPU count; `init`, `config.json`
  and profiles refuse it before a job is submitted. Chapter 12 explains the
  shapes a 6-GPU node allows and how a project shapes its Parabricks jobs
  through stage inputs.

## 0.9.2

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.4.0,
unchanged since 0.9.0. No workflow changed since 0.9.1.

- The run report's Provenance names the container backend, the node-local
  scratch (`scratch_env` and its mount, or `off`) and the resource policy
  file of the install's `miniwdl.cfg`, so a report says whether tasks had
  the node's local disk.
- `tests/bundle.sh` runs without Apptainer again: its check of the
  `ugc_slurm_singularity` backend resolves the class through miniwdl's
  plugin registry instead of constructing it, whose initialisation probes
  `apptainer --version`. This is what failed the CI of 0.9.0 and 0.9.1
  after every other step had passed.

## 0.9.1

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.4.0,
unchanged since 0.9.0.

- The pbmm2 reference index is built once at install instead of once per
  input BAM. Upstream's alignment subworkflow and its fail-reads path build
  the same 5.8 GB index for every BAM; the call cache catches it only after
  the first build, so a campaign's first wave starts hundreds of identical
  builds. `install-bundle.sh` builds it once with the bundled pbmm2 image
  (about a minute, 12 GB of memory, verified against the md5 in
  `references.lock`; it says how to build it in a job when the login node
  cannot), the reference map gains `pbmm2_index`, and two derived files,
  `workflows/overrides/pbmm2.wdl` and `workflows/overrides/upstream.wdl`,
  pass it to every alignment in the `singleton` and `upstream` stages;
  without an index they behave as upstream. Alignment inputs change, so a
  sample replayed into a project on this version aligns once more.

## 0.9.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.4.0.

- Node-local scratch for tasks. The plugin registers a container backend,
  `ugc_slurm_singularity`, which the rendered `miniwdl.cfg` selects: it is
  miniwdl-slurm's `slurm_singularity`, and when the new site key
  `TASK_SCRATCH_ENV` names the node's per-job scratch variable (`TMPDIR`,
  `SNIC_TMP`) it binds that directory at `/mnt/scratch` inside every task's
  container and sets `TMPDIR` to it, so pbmm2's sort spill files, GLnexus's
  database and other TMPDIR-aware temporaries leave the shared file system.
  A job whose variable is unset, unwritable or tmpfs runs as before, and
  one `[ugc-wgw] scratch:` line in the job's `slurm-<jobid>.out` says which.
  The installer validates the key and prints the setting; `ugc-wgw
  resources` shows `scratch=` in its summary line. Until now `TMPDIR` was
  never set inside a container, so GLnexus's database had always been on
  the shared work directory; the guide is corrected.

## 0.8.2

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.3.0,
unchanged since 0.8.0.

- Failure classification reads the failed task's `stdout.txt` as well
  (after the SLURM log and `stderr.txt`) and counts a CUDA host-memory
  failure (`cudaErrorMemoryAllocation`, `CUDA_ERROR_OUT_OF_MEMORY`) as kind
  `resource`: Parabricks logs to stdout, leaves stderr empty and exits 255
  without a SLURM out-of-memory state when a pinned allocation does not
  fit.
- The worked resource policy gives `run_parabricks_deepvariant` 256G: the
  declared 192 GiB held at 19–25x and failed on deeper samples. Chapters 09
  and 12 explain where Parabricks' log is and why SLURM shows nothing.

## 0.8.1

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.3.0,
unchanged since 0.8.0.

- `samples add` and `samples check` show the progress of the input check
  (`inspecting read files N/M (p%), elapsed, about T left`, rewritten in
  place on a terminal, one line every ten seconds otherwise) and end with
  `inspected N read file(s) in T`; `-v` logs each sample's yield.
- Guide: the examples show the current version and the installer's current
  report; the quick start describes the input check, `--profile cpu` and
  `ugc-wgw usage`; the command reference lists exit code 2.

## 0.8.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.3.0.

- Input BAM check: `ugc-wgw samples add` opens every read file with a
  stdlib BGZF reader. An empty, truncated or read-less BAM is refused
  (`--drop-empty` registers the sample without its empty or read-less
  files; `--no-inspect` keeps only the existence check); aligned inputs and
  missing read groups warn; the read count and bases are estimated from
  the first records (exact from a `.pbi`) and recorded with the sample;
  `input_thresholds` (`file_reads_min`, `sample_gbases_min`; a config and
  profile key) warn about thin files and samples. New `ugc-wgw samples
  check [ID...] [--stored]` lists size, reads, bases, movie and problems
  per file.
- `submit` checks the raw read files just before each run that names them
  (present, same size as registered, BGZF end marker) and fails the attempt
  as an input error without a SLURM job when they do not hold; `--dry-run`
  prints the result.
- Outputs across file systems: a hardlink miniwdl cannot make across file
  systems (a cached output of a run whose results live elsewhere, an input
  passed through) becomes a symlink with a warning instead of
  `Invalid cross-device link`; the installer writes `ugc_wgw_miniwdl.pth`
  into the venv for it and warns when the prefix, references and call cache
  are on different file systems, as `ugc-wgw init` does for the results.

## 0.7.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is unchanged.

- Site profiles: `ugc-wgw init --profile NAME|FILE` applies a JSON profile
  of the project keys that depend on the cluster and the sample set
  (`stage_inputs`, the DeepVariant flavour and GPU choice, `max_inflight`,
  prices, QC thresholds) from `<prefix>/profiles/`, which the installer
  fills from `backends/hpc/profiles/` once and never overwrites, or from
  the shipped examples; flags win over the profile; `config.json` records
  the profile's name, checksum and the keys it set, `ugc-wgw stage-inputs`
  marks its overrides and the report names it.
- Shipped profiles `cpu` and `parabricks`: no alignment chunking, pbmm2 on
  24 threads so two alignments share a 48-core node, 20 runs in flight.

## 0.6.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is unchanged.

- `ugc-wgw usage` and a "Resource usage" section in the run report:
  core-hours (allocated and requested), CPU-hours used with the efficiency,
  GPU-hours, memory GB-hours, queue wait, jobs, disk usage (`--sizes`) and a
  cost estimate from prices you supply (`--price KEY=VALUE`, `config.json`
  `prices`; `--basis allocated|requested`), per stage, task and subject,
  with the mean per sample (clean and as run).
- The driver reads `sacct` for every run's jobs when the run finishes and
  keeps it as `accounting.json` next to the manifest (compute time without
  the queue wait, the wait itself, CPUs, `TotalCPU`, peak memory, GPUs,
  disk); `usage --collect` reads it later for older runs; attempts without
  accounting are estimated from the workflow log and labelled so.
- The smoke harness's `sacct` step also prints the driver's per-task table.

## 0.5.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is unchanged.

- `ugc-wgw stage-inputs`: every input of every stage with its type, default,
  description, who fills it (the driver, a `config.json` key, or nobody, so
  it is yours under `stage_inputs`) and the project's current override; it
  warns about overrides `submit` would refuse. `--nested` adds the
  call-qualified inputs of the tasks inside each entrypoint (thread counts,
  memory, tool options), listed by miniwdl itself, and checks dotted
  overrides against them.
- Upstream's task thread counts are settable as such nested keys
  (`upstream.pbmm2.pbmm2_align_wgs.threads`), which size both the SLURM
  request and the command. `ugc-wgw resources` says so when a policy row
  sets `cpu` for a task that interpolates its thread count, and names the
  entrypoint input where one exists (`hifiasm_threads`); the task inventory
  records those inputs.
- Guide: `examples/resources.tsv`, a worked site policy for 48-core, 384 GB
  nodes with its reasoning in chapter 12; every entrypoint input now carries
  a `parameter_meta` description.
- Smoke harness: `run.sh init --hifiasm-threads N`; the drivers it tees run
  with `--color never`.
- A non-object `stage_inputs.<stage>` in `config.json` is a clean error.

## 0.4.1

Patch release: the shell scripts pass shellcheck 0.9.0, which the CI runs.
No workflow or driver change; a 0.4.0 bundle stays valid.

## 0.4.0

First public release. Built on PacBio HiFi-human-WGS-WDL v4.0.0
(`15e82cb9`) and its reference data container `GRCh38_GIABv3`.

- Seven entrypoints: `singleton`, `upstream`, `cohort_call`, `downstream`,
  `cohort_merge`, `cohort_freq`, `assembly`; standalone, joint and assembly
  modes; cohort stages scattered over 26 region shards at every cohort size.
- The `ugc-wgw` driver: sample sheet and frozen cohorts, plan and submit with
  a concurrency cap, automatic retries of transient failures, orphan
  cancellation, SQLite state, per-run provenance manifests, `status`,
  `progress`, `logs`, `report` and `summary`.
- Offline bundle build and install with digest-pinned containers, engine
  wheels and the reference data container; site caps and a per-task resource
  policy applied by a miniwdl plugin; GPU DeepVariant and Parabricks as a
  project switch.
- Per-sample and per-cohort HTML analysis summaries with QC flags; a
  self-contained run report.
- `scripts/manifest-to-samples.py` turns a per-file delivery manifest into
  the sample sheet.
