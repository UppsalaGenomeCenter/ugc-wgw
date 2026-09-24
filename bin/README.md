# bin/ — the `ugc-wgw` driver

`bin/ugc-wgw` sequences samples and cohorts through the stages of a mode, generates
each run's inputs JSON from the state DB, launches `miniwdl run` with a concurrency
cap, and writes `run_manifest.json` when a run finishes. Specification:
`docs/DESIGN.md` §8 (responsibilities, CLI, state DB) and §9 (results layout,
provenance). Python 3.9+, stdlib only; the package lives in `bin/ugc_wgw/`.
The operator's walkthrough is `docs/guide/` (chapters 05, 06 and 09); this
page is the reference.

## Layout of a project

```
<project>/.ugc-wgw/config.json     paths and defaults (code dir, miniwdl, cfg, ref maps, registry, stage_inputs,
                               auto_retry_max, backoff_seconds, cancel_orphans, lease_seconds, progress_interval)
<project>/.ugc-wgw/state.sqlite    samples, cohorts, runs, events (never edit by hand; schema migrated on open)
<project>/.ugc-wgw/lock            flock + lease line: pid= host= since= last_seen= (heartbeat every poll)
<project>/.ugc-wgw/logs/ugc-wgw.log    human log;  events.jsonl: one JSON event per line
<results>/samples/<id>/<ugc_wgw_version>/<stage>/attempt-<n>/   miniwdl run dir: inputs.json, run.json,
                                                             outputs.json, out/, workflow.log (plain),
                                                             workflow.log.json, run_manifest.json,
                                                             call-<task>/slurm_singularity.log.txt (SLURM job id)
<results>/samples/<id>/<ugc_wgw_version>/<stage>/current -> attempt-<n>
<results>/cohorts/<id>/<ugc_wgw_version>/<stage>/attempt-<n>/
```

Every attempt gets a fresh directory because miniwdl refuses to finish in a
directory that already holds an `out/`. `<stage>/current/run_manifest.json` is
the canonical provenance path of a subject and stage.

## Commands

```
ugc-wgw init <project_dir> --install <root>/current --ref-map F [--deepvariant cpu|gpu|parabricks --gpu-type T]
ugc-wgw init <project_dir> --code <repo> --miniwdl <exe> --cfg <miniwdl.cfg> ...     # dev machine, no install prefix
ugc-wgw samples add <samples.tsv> [--no-check] [--replace]
ugc-wgw samples remove <id>... [--force]                # refused for cohort members; --force deletes run rows too
ugc-wgw samples list [--mode M] [--stage S] [--status failed,running] [--json]
ugc-wgw cohort freeze <cohort_id> --samples <ids.txt>
ugc-wgw submit --mode {standalone,joint,assembly} [--stage S] [--samples ID... | --samples-file F | --cohort ID]
           [--max-inflight N] [--poll-interval SEC] [--any-version] [--dry-run] [--takeover]
ugc-wgw retry  --stage S [--mode M] [--samples ... | --cohort ID] [--max-attempts N] [--dry-run] [--takeover]
ugc-wgw status [--mode M] [--cohort ID] [--failed] [--json]   # --failed/--json add the failure kind and message
ugc-wgw logs   <subject_id> --stage S [--tail N] [--attempt N] [--follow] [--json]
ugc-wgw inputs <subject_id> --stage S [--mode M] [--cohort ID]
ugc-wgw progress [--mode M] [--stage S] [--samples ... | --cohort ID] [--json]  # done/total per stage, ETA
ugc-wgw report [--mode M] [--cohort ID] [--any-version] [--sizes] [--out FILE]   # self-contained HTML run summary
ugc-wgw summary [--samples ... | --cohort ID] [--mode M] [--any-version] [--out-dir DIR] [--force]
            [--threshold KEY=VALUE ...] [--jobs N]  # per-sample and cohort analysis pages (what was found, QC flags)
ugc-wgw resources [--stage S] [--changed] [--json]     # what each task asks SLURM for: declared → site caps → policy
```

`--project DIR` (default `.`) selects the project; `-v` streams INFO logs to stderr;
`--color auto|always|never` (default `auto`, or `UGC_WGW_COLOR`) colours that stderr copy
by meaning on a terminal (the file log stays plain; `NO_COLOR` is honoured).
`ugc-wgw report` is the technical view of a campaign (runs, timings, failures);
`ugc-wgw summary` is the analysis view (coverage, variants, phasing, repeats,
methylation, targeted callers, software versions) rendered from per-sample
JSON digests under `<results>/reports/summary/` (guide chapter 07).

## Samples TSV

Header row required. Columns: `sample_id` (`[A-Za-z0-9._-]+`, becomes a directory
name), `sex` (`MALE`, `FEMALE` or blank; advisory, upstream infers the sex),
`hifi_reads` (one or more uBAM paths, comma-separated), `fail_reads` (optional,
same form), `father_id`, `mother_id` (optional). Extra columns are kept as metadata.
Paths must exist unless `--no-check`.

## Modes and gating

`standalone` = `singleton` → `cohort_merge` → `cohort_freq`; `joint` =
`upstream` → `cohort_call` → `downstream` → `cohort_merge` → `cohort_freq`
(needs `--cohort`); `assembly` = the single per-sample stage `assembly` (no
cohort). `cohort_freq` turns the cohort's joint small-variant VCF and merged
SV VCF into sites-only frequency files; its inputs builder takes each
member's sex from the sheet, else from the sex upstream inferred. In assembly mode a sample whose
`father_id` and `mother_id` are both registered samples with reads is
trio-binned with the parents' reads; otherwise it runs in sample mode with a
warning. `assembly_use_parents = false` in `config.json` (or `ugc-wgw init
--no-assembly-parents`) forces sample mode. A stage is runnable for a
subject when every earlier stage of the mode succeeded for it (for a cohort: for
every frozen member) **at the current ugc-wgw version**; `--any-version` accepts
older successes and records their version in the manifest. A failed or
cancelled run blocks its subject until `ugc-wgw retry` starts a new attempt, with
one exception: a failure of kind `transient` (see below) is re-attempted by the
submit loop itself after a backoff of `backoff_seconds × 2^(attempt-1)`
(capped at one hour) while the attempt number is at most `auto_retry_max`
(default 2). Terminal runs are never modified.

## Failure kinds

Every finished attempt is classified from miniwdl's error document and the
failed task's `slurm_singularity.log.txt` and `stderr.txt`
(`ugc_wgw/failures.py`):

| Kind | Rule | Driver's reaction |
|---|---|---|
| `transient` | `driver_lost`, `Interrupted`, `Terminated`, `launch_error`, `killed`, `NoResult`; SLURM preemption / node-failure text | automatic re-attempt after backoff |
| `resource` | exit status 137 or 253; `oom_kill`, `OUT_OF_MEMORY`, `DUE TO TIME LIMIT`, `TIMEOUT`, or an sbatch refusal (`Requested node configuration is not available`, `CPU count per node can not be satisfied`) in the task's logs | blocked until `ugc-wgw retry` (cap or change the request first: site limits, a policy row, or `stage_inputs`) |
| `input` | `InputError`, `DownloadFailed`; `No such file`, `does not exist`, `Permission denied` from a failed command | blocked until the input is fixed and `ugc-wgw retry` |
| `tool` | any other `CommandFailed` / `OutputError` / class | blocked until `ugc-wgw retry` |
| `version` | `version_mismatch` | blocked; fix the project's install or override |
| `cancelled` | the run was stopped by the driver | blocked until `ugc-wgw retry` |
| `unknown` | no error class at all | blocked until `ugc-wgw retry` |

The kind, a short message and the evidence line are on the run row
(`error_kind`, `error_message`), in `run_manifest.json` (`error`), in the
`run.failed` event and in `ugc-wgw status --failed` / `--json`.

## Task resources

What a task asks SLURM for is decided in three layers, none of them in the
driver: the WDL declares it, miniwdl caps it (`[task_runtime] cpu_max`,
`memory_max` from the site file), and the `ugc_wgw_resources` task plugin
(`plugins/ugc_wgw_miniwdl`, shipped in the bundle) applies the site policy
`<prefix>/resources.tsv` per task (cores, memory, wall time, partition,
constraint). `ugc-wgw resources` prints the result per task from the same
inventory (`backends/hpc/resources.declared.tsv`), the same rendered
`miniwdl.cfg` and the same parser, so it shows what sbatch will be asked
for; `submit` and `retry` refuse a malformed policy before launching
anything, warn when a row matches no task or a partition row is defeated
by `[slurm] extra_args`, and record the policy's path and sha256 in
`run_manifest.json` (`resource_policy`). The report tabulates the caps
and the policy changes miniwdl logged. Guide: chapter 12.

GPU DeepVariant is a project switch: `ugc-wgw init --deepvariant gpu|parabricks
[--gpu-type a100] [--parabricks-gpus 4]` (or the `deepvariant`, `gpu_type`,
`parabricks_gpus` keys of `config.json`) sets `use_gpu`,
`use_parabricks_deepvariant`, `gpuType` and Parabricks' GPU count for
`singleton` and `upstream`; the install needs `SINGULARITY_NV=1` and
`SLURM_PARTITION_GPU` in its site file, and `ugc-wgw resources` shows the
resulting gres per task and warns when either is missing (guide chapter 12).

## How submit runs

`ugc-wgw submit` is a foreground loop: it takes `.ugc-wgw/lock`, reconciles runs left
active by a previous driver, launches up to `--max-inflight` `miniwdl run`
processes, polls them, finalizes each finished run (manifest, DB, `current`
link), refreshes the lease, and launches what became runnable, until nothing
is runnable or in flight; while a subject is backing off after a transient
failure it waits for that moment instead of exiting. Run it under `tmux` or
as a SLURM job. Ctrl-C or SIGTERM forwards SIGTERM to the children (miniwdl
cancels its SLURM steps), waits up to two minutes, marks them `cancelled`
and, for a child that had to be SIGKILLed, runs `scancel` on its job ids;
exit status 130. Exit 1 if any subject's latest attempt failed this session
(a subject recovered by an automatic retry does not count).

One driver per project. `.ugc-wgw/lock` is a flock plus a lease line
(`pid= host= since= last_seen=`) refreshed every poll. A `submit` on another
host refuses while the lease is younger than `lease_seconds` (900) and takes
over an expired one; `--takeover` takes a fresh one when the holder is known
to be dead, after which that holder stops itself on its next heartbeat. The
state database is SQLite in WAL mode: keep it on the host that runs the
driver or on a filesystem with working POSIX locks, never on plain NFS.

Total SLURM pressure is `max_inflight × [scheduler] task_concurrency` of the
miniwdl config.

Reconciliation after a crash: runs the DB still calls active are settled from
their run directory (`run.json`, `outputs.json`, `error.json`); if there is no
result and the recorded driver PID is dead on this host, or the run belongs to
another host whose lease has expired (or this process now holds the lease),
the run becomes `failed/driver_lost`, its SLURM job ids are cancelled, and the
next `submit` re-attempts it automatically (kind `transient`). Runs of another
host with a fresh lease are left alone with a warning. Every command that
reads run state reconciles first.

`--dry-run` prints the resource summary (caps and policy), then, for the
first wave, the attempt directory, the generated inputs JSON and the exact
`miniwdl run` command line, plus the blocked and active lists, and changes
nothing.

## Tests

```
python3 -m unittest discover -s tests/driver -t . -v
```

The suite runs the driver against `tests/driver/fake_miniwdl.py`, which
reproduces miniwdl's run-directory conventions (input validation against the real
WDL input blocks, `run.json`, `outputs.json`, `out/`, `error.json`, exit codes).

## Not implemented

A per-task container mapping in `run_manifest.json` (today: the union of the
image manifests); driving one project from several hosts at once (by design,
see the lease); `samples remove` of a cohort member (cohorts are immutable);
combined modes that run assembly beside the variant stages in one session
(deferred, `docs/DESIGN.md` §5.1; use a second project).
