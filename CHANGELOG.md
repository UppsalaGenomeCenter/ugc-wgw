# Changelog

Release notes of the public copy. Each version corresponds to a verified
offline bundle; the upstream tag it is built on is named first.

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
