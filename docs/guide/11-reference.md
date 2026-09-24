# Reference

## `ugc-wgw` command line

Global options come before the command: `--project DIR` (default `.`) selects
the project, `-v` (`--verbose`) streams INFO log lines to stderr as well as to
the log file. `--color auto|always|never` colours the stderr copy of the log by
meaning (green `run.success`, red failures, cyan submissions, dimmed waiting
lines, magenta driver actions, bold `subject=`/`stage=` values, red failure
details); `auto` (the default, or the `UGC_WGW_COLOR` variable) colours only
when stderr is a terminal and `NO_COLOR` is unset. The file log
`.ugc-wgw/logs/ugc-wgw.log` is never coloured.

### `ugc-wgw init <project_dir>`

| Flag | Default | Meaning |
|---|---|---|
| `--install DIR` | | Install version directory (`<prefix>/versions/<v>` or `<prefix>/current`); derives the code, miniwdl, config and venv paths. |
| `--code DIR`, `--miniwdl EXE`, `--cfg FILE` | | The paths individually; either `--install` or all three. |
| `--venv DIR` | | Engine venv, for version probing only. |
| `--results DIR` | project dir | Results root. |
| `--ref-map FILE` | required | The rendered reference map (`ref_map=` in the installer's report). |
| `--registry STR` | `ghcr.io/uppsalagenomecenter` | Registry prefix of the ugc-built images. |
| `--max-inflight N` | `4` | Default concurrency for `submit`. |
| `--poll-interval SEC` | `30` | Default poll period. |
| `--no-assembly-parents` | off | Assembly never trio-bins. |
| `--deepvariant cpu\|gpu\|parabricks` | `cpu` | Small-variant caller of `singleton`/`upstream` (chapter 12). |
| `--gpu-type TYPE` | empty | SLURM gres type of the GPU tasks; empty = any GPU. |
| `--parabricks-gpus N` | `4` | GPUs per Parabricks task. |

### `ugc-wgw samples add <tsv>`

| Flag | Meaning |
|---|---|
| `--no-check` | Do not require read paths to exist. |
| `--replace` | Overwrite samples that are already registered. |

### `ugc-wgw samples remove <ID>...`

| Flag | Meaning |
|---|---|
| `--force` | Also delete the samples' run rows. |

Refused for a member of any cohort, for a sample with an active run, and,
without `--force`, for a sample with recorded runs. Results on disk are never
touched. Every id is checked before anything is deleted.

### `ugc-wgw samples list`

| Flag | Default | Meaning |
|---|---|---|
| `--mode M` | `standalone` | Which sample stages become columns. |
| `--stage S` | | Only this stage. |
| `--status a,b` | | Keep rows whose latest status is one of these. |
| `--any-version` | off | Latest run at any version. |
| `--json` | off | JSON instead of the table. |

### `ugc-wgw cohort freeze <cohort_id> --samples <file>`

Freezes the ordered ID list in the file; every ID must be registered. No
flags besides `--samples`.

### `ugc-wgw submit` and `ugc-wgw retry`

| Flag | Default | Meaning |
|---|---|---|
| `--mode M` | `standalone` | `standalone`, `joint` or `assembly`. |
| `--stage S` | all stages of the mode | One stage. Required for `retry`. |
| `--samples ID ...` / `--samples-file FILE` | all | Sample selection (mutually exclusive). |
| `--cohort ID` | | Cohort context; required from the first stage in `joint`. |
| `--max-inflight N` | config | Concurrent miniwdl runs. |
| `--poll-interval SEC` | config | Seconds between polls. |
| `--any-version` | off | Accept prerequisites from any version. |
| `--dry-run` | off | Print the first wave and the plan; change nothing. |
| `--takeover` | off | Take the project lease from a driver on another host that is known to be dead. |
| `--max-attempts N` | `3` | `retry` only: refuse a new attempt at or past this count. |

### `ugc-wgw status`

| Flag | Default | Meaning |
|---|---|---|
| `--mode M` | `standalone` | Stage columns. |
| `--cohort ID` | | Members of this cohort plus the cohort row. |
| `--failed` | off | Only rows with a failed or cancelled cell, plus `kind` and `message` columns. |
| `--any-version` | off | Latest run at any version. |
| `--json` | off | JSON; rows carry `kind` and `message`. |

### `ugc-wgw logs <subject_id> --stage S`

| Flag | Default | Meaning |
|---|---|---|
| `--attempt N` | latest | Which attempt. |
| `--tail [N]` | all; `50` when given without a number | Last N lines. |
| `-f`, `--follow` | off | Keep printing until the run ends, then print its final status, class, kind and message. |
| `--json` | off | Show `workflow.log.json` (JSON lines) instead of the plain `workflow.log`. |

### `ugc-wgw progress`

| Flag | Default | Meaning |
|---|---|---|
| `--mode M`, `--stage S`, `--samples`, `--samples-file`, `--cohort ID`, `--any-version` | as for `submit` | The selection to summarise. |
| `--json` | off | The summary as JSON (what the `submit.progress` event carries). |

### `ugc-wgw report`

| Flag | Default | Meaning |
|---|---|---|
| `--mode M` | all | Only runs of this mode. |
| `--cohort ID` | all | Only the cohort and its members. |
| `--any-version` | current | Runs of every ugc-wgw version. |
| `--sizes` | off | Measure each run's `out/` (walks the results tree). |
| `--out FILE` | `<results>/reports/ugc-wgw-report-<stamp>.html` | Output file; the path is printed. |

### `ugc-wgw summary`

| Flag | Default | Meaning |
|---|---|---|
| `--samples ID...` / `--samples-file F` | every registered sample | Samples to summarise. |
| `--cohort ID` | none | Also write the cohort page; covers the cohort's members. |
| `--mode M` | `standalone` | `standalone` (singleton) or `joint` (upstream + downstream). |
| `--any-version` | current | Runs of every ugc-wgw version. |
| `--out-dir DIR` | `<results>/reports/summary` | Where `<id>.summary.html` and `<id>.summary.json` go. |
| `--force` | off | Rebuild the per-sample digests even when their runs are unchanged. |
| `--threshold KEY=VALUE` | config or defaults | QC threshold override; repeatable. Keys: `depth_mean_min`, `mapped_read_percent_min`, `read_quality_median_min`. |
| `--jobs N` | `1` | Build digests in N processes. |

Exit 1 when nothing could be written (no sample has a successful run of the
mode), on an unknown sample, cohort or threshold key.

### `ugc-wgw resources`

| Flag | Default | Meaning |
|---|---|---|
| `--stage S` | all | Only tasks of this stage. |
| `--changed` | off | Only tasks a site cap or a policy row changes. |
| `--json` | off | Rows as JSON. |

Prints, per task, the declared request, the effective cores, memory, wall
time, partition and constraint with `(clamp)` or `(policy)` marks, the GPU
request of the two GPU tasks (`gpus`: `[type:]N`), the rows that applied
and the inventory's notes; the summary line names the project's
`deepvariant` flavour. Warns about rows that match no task, a missing
policy file, a partition defeated by `[slurm] extra_args`, an engine
without the plugin, and, when a GPU flavour is selected, `run_options`
without `--nv` or no GPU partition. Exit 1 on a malformed policy (chapter
12).

### `ugc-wgw inputs <subject_id> --stage S`

| Flag | Default | Meaning |
|---|---|---|
| `--mode M` | `standalone` | Affects the wiring of `downstream`, `cohort_merge` and `cohort_freq`. |
| `--cohort ID` | | Cohort context. A cohort subject without `--cohort` is looked up by its own ID. |
| `--any-version` | off | Read earlier outputs at any version. |

### Exit codes

`0` success; `1` an error was printed or a run failed in this session;
`130` stopped by SIGINT or SIGTERM.

## Files of a project

| Path | Content |
|---|---|
| `<project>/.ugc-wgw/config.json` | Configuration; keys below. |
| `<project>/.ugc-wgw/state.sqlite` | Samples, cohorts, runs, events. Never edit. |
| `<project>/.ugc-wgw/lock` | The lease of the running `submit` or `retry`: `pid= host= since= last_seen=`, flock-ed on its host, `last_seen` refreshed every poll. |
| `<project>/.ugc-wgw/logs/ugc-wgw.log` | Human log, UTC. |
| `<project>/.ugc-wgw/logs/events.jsonl` | One JSON event per line. |

## `config.json` keys

| Key | Default | Meaning |
|---|---|---|
| `schema` | `1` | Config schema version. |
| `code_dir` | | Installed code directory; `VERSION` is read from it. |
| `miniwdl` | | miniwdl executable. |
| `miniwdl_cfg` | | Config passed as `--cfg` to every run. |
| `venv_dir` | `null` | Engine venv. |
| `results_dir` | project dir | Results root. |
| `ref_map_file` | | Every stage. |
| `ugc_wgw_container_registry` | `ghcr.io/uppsalagenomecenter` | `cohort_merge`, `assembly`. |
| `backend` | `"HPC"` | Every stage. |
| `preemptible` | `true` | Every stage; no effect on the HPC. |
| `max_inflight` | `4` | |
| `poll_interval` | `30.0` | |
| `assembly_use_parents` | `true` | |
| `auto_retry_max` | `2` | Automatic re-attempts of transient failures while the attempt number is at most this; `0` disables. |
| `backoff_seconds` | `300` | Wait before the first automatic re-attempt; doubles per attempt, capped at 3600. |
| `cancel_orphans` | `true` | Run `scancel` on the SLURM jobs of runs settled as `driver_lost` or killed after the stop grace. |
| `lease_seconds` | `900` | A lease older than this is expired; a driver on another host may then take the project. |
| `progress_interval` | `300` | Seconds between progress summaries in the log while `submit` runs (also logged when the counts change); `0` disables. |
| `stage_inputs` | `{}` | `{"<stage>": {"<input>": value}}`; `ugc_wgw_<stage>` also accepted and wins. |
| `summary_thresholds` | `{}` | `ugc-wgw summary` QC thresholds, e.g. `{"depth_mean_min": 25}` (chapter 07). |
| `project_url` | the public repository | Repository link printed in the analysis summaries. |
| `deepvariant` | `"cpu"` | `cpu`, `gpu` or `parabricks`: the small-variant caller of `singleton` and `upstream` (chapter 12). |
| `gpu_type` | `""` | SLURM gres type behind `--gres gpu:<type>:N`; empty = `gpu:N`. |
| `parabricks_gpus` | `4` | GPUs per Parabricks task. |

## Sample sheet

Tab-separated, header row required. Columns `sample_id` (required,
`[A-Za-z0-9._-]+`), `hifi_reads` (required, comma-separated uBAM paths),
`sex` (`MALE`, `FEMALE`, blank), `fail_reads`, `father_id`, `mother_id`;
other columns are kept as metadata. Example: `examples/samples.tsv`.
`scripts/manifest-to-samples.py` writes one from a per-file delivery
manifest with a PLINK-style pedigree (chapter 05).

## Cohort list

One sample ID per line; `#` comments and blank lines ignored; only the first
token of a line is read; order is significant. Example: `examples/cohort.txt`.

## Site file

`KEY=VALUE` lines for `install-bundle.sh --site`: `SLURM_PARTITION`,
`SLURM_ACCOUNT`, `SLURM_EXTRA_ARGS` (never the partition or account),
`SLURM_PARTITION_GPU`, `SLURM_ACCOUNT_GPU`, `SINGULARITY_NV` (0),
`TASK_CPU_MAX` (0), `TASK_MEMORY_MAX` (0), `TASK_RESOURCES`
(`<prefix>/resources.tsv`), `TASK_CONCURRENCY` (50), `TASK_TIME_MINUTES`
(4320), `CALL_CACHE_DIR` (`<prefix>/call_cache`). Example:
`examples/site.cfg`; details in chapter 04.

## Resource policy

`<prefix>/resources.tsv` (the site file's `TASK_RESOURCES`), tab-separated
with the header `task cpu memory time partition constraint`; `#` lines
ignored; `-` or empty keeps the declared value. `task` is a name from
`backends/hpc/resources.declared.tsv` or a glob; every matching row applies
in file order, later cells win. `memory` needs a unit, binary whatever the
spelling (`256G` = 256 GiB); `time` is `sbatch` syntax (minutes, `HH:MM:SS`,
`D-HH:MM:SS`). Extra columns `stages`, `source`, `notes` are ignored. Rows
are applied after `cpu_max`/`memory_max` and are not capped by them.
Chapter 12.

## `install-bundle.sh`

| Flag | Meaning |
|---|---|
| `--bundle TAR` | Required. |
| `--prefix DIR` | Required. Install root. |
| `--references DIR` | Default `<prefix>/references`. |
| `--site FILE` | Site values. |
| `--python EXE` | Default `python3`; must match the bundle's Python minor version. |
| `--verify-only` | Verify, check the workflows, install nothing. |
| `--activate` | Repoint `<prefix>/current`. |
| `--force` | Overwrite an existing `versions/<v>`. |
| `-h`, `--help` | Print the header comment of the script. |

## `miniwdl.cfg` keys an operator may touch

`[scheduler] task_concurrency`; `[task_runtime] cpu_max`, `memory_max`
(bytes, `0` = none); `[task_runtime] defaults` (`time_minutes`,
`maxRetries`, `slurm_partition`, `slurm_account`, `slurm_partition_gpu`,
`slurm_account_gpu`; one JSON line); `[slurm] extra_args`;
`[singularity] run_options` (`--nv` for the GPU tasks, from
`SINGULARITY_NV`); `[call_cache] dir`; `[ugc_wgw] resources`. Everything else
is required as rendered (chapter 04).

## Attempt directory

`inputs.json`, `miniwdl.stdout`, `miniwdl.stderr`, `workflow.log` (plain),
`workflow.log.json`, `run.json`, `outputs.json`, `error.json`, `out/`,
`call-<task>/` (with `slurm_singularity.log.txt`), `run_manifest.json`;
`current` in the stage directory points at the latest finished attempt
(chapter 07).

## Run states, error classes, events

States: `pending`, `submitted`, `running` (active); `success`, `failed`,
`cancelled` (terminal). Driver error classes: `driver_lost`,
`version_mismatch`, `NoResult`, `launch_error`, `killed`; miniwdl classes as
reported in `error.json` (`CommandFailed` and others). Failure kinds:
`transient`, `resource`, `input`, `tool`, `version`, `cancelled`, `unknown`
(chapter 09). Events: `sample.added`, `sample.warning`, `sample.removed`,
`cohort.frozen`, `submit.start`, `submit.stop`, `run.created`,
`run.submitted`, `run.running`, `run.blocked`, `run.auto_retry`,
`run.terminating`, `run.success`, `run.failed`, `run.cancelled`,
`run.reconciled`, `run.scancel`, `submit.progress`. miniwdl log lines the
report reads: `runtime.cpu adjusted to host limit` (a site cap) and
`ugc-wgw resource policy applied` (a policy row).

## Blocked reasons

| Text | Meaning |
|---|---|
| ``failed (attempt N, <class>, <kind>); run `ugc-wgw retry --stage S` `` | A terminal failure of a kind the driver does not retry on its own; `, auto-retries exhausted` when it was transient but past `auto_retry_max`. |
| `backing off until <ts> after transient <class> (attempt N)` | The driver will re-attempt at that time. |
| `failed after N attempts (max M)` | The retry cap was reached. |
| `waiting for <stage> of <ids> (+N more)` | An earlier sample stage has not succeeded for these samples. |
| `waiting for <stage> of cohort <id>` | An earlier cohort stage has not succeeded. |
| `<stage> needs --cohort <id>` | A cohort stage was selected without a cohort. |

## Environment assumptions

`python3` 3.9 or newer on `PATH` for the driver (`bin/ugc-wgw` uses
`/usr/bin/env python3`); `apptainer` on `PATH`; SLURM `sbatch` (and `scancel`,
optional) on the host running the driver; the results root, `<prefix>` and the
reference directory readable from compute nodes; `$TMPDIR` on compute nodes
large enough for GLnexus; one driver per project, whose `.ugc-wgw/state.sqlite`
is not on plain NFS.
