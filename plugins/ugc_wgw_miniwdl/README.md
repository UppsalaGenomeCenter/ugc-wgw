# ugc-wgw-miniwdl

A miniwdl task plugin. It reads the site's per-task resource policy
(`<prefix>/resources.tsv`, named by `[ugc_wgw] resources` in the rendered
`miniwdl.cfg`) and rewrites a task's evaluated `cpu`, `memory`,
`time_minutes`, `slurm_partition` and `slurm_constraint` before
miniwdl-slurm turns them into `sbatch` arguments.

Why a plugin: every upstream task hard-codes its cores and memory as
private declarations, so no inputs file can change them, and editing 37
vendored tasks would fight every upstream sync. miniwdl's task-plugin
hook runs after the runtime block is evaluated and clamped
(`[task_runtime] cpu_max`, `memory_max`) and before the container is
submitted, which is exactly the point where a site policy belongs.

Layout:

| File | Role |
|---|---|
| `ugc_wgw_miniwdl/policy.py` | The TSV format: parsing, units, matching, application. Stdlib only; the driver imports it too (`ugc-wgw resources`). |
| `ugc_wgw_miniwdl/resources.py` | The coroutine miniwdl calls per task (`miniwdl.plugin.task` entry point `ugc_wgw_resources`). |
| `ugc_wgw_miniwdl/scratch.py`, `ugc_wgw_miniwdl/scripts/exec_scratch.sh` | The `ugc_slurm_singularity` container backend (`miniwdl.plugin.container_backend`): miniwdl-slurm's `slurm_singularity` with the job's node-local scratch bound into every task's container and `TMPDIR` pointing at it, when `[ugc_wgw] scratch_env` names the job variable. |
| `ugc_wgw_miniwdl/crossdev.py` | Hardlink, else symlink: wraps miniwdl's `symlink_force` so an output on another file system becomes a symlink instead of `Invalid cross-device link`. Imported by `resources.py` and by the `ugc_wgw_miniwdl.pth` the installer writes into the venv. |
| `pyproject.toml` | Package metadata and the entry point. |

Build and ship: `scripts/make-bundle.sh` builds the wheel on the dev
machine (`pip wheel --no-deps`, network for the setuptools build backend)
and lists it in the wheelhouse's `requirements.txt`;
`scripts/install-bundle.sh` installs it into the version's venv offline
and checks the entry point. Bump `version` in `pyproject.toml` and
`__init__.py` whenever the plugin changes.

Development: `.venv/bin/pip install -e plugins/ugc_wgw_miniwdl`, then
`.venv/bin/miniwdl --version` lists `ugc_wgw_resources`. Tests:
`python3 -m unittest discover -s tests/plugin -t . -v` (no miniwdl
needed).

Format, precedence and examples: `docs/guide/12-resources.md`.

## Outputs across file systems

miniwdl's `[file_io] output_hardlinks = true` (our configuration, so that
`delete_work = success` can reclaim work directories) links every output
into the run's `out/` with `os.link`. A target on another file system, the
cached output of an earlier run whose results live elsewhere or an input a
workflow passes through, fails with `OSError: [Errno 18] Invalid
cross-device link`. `crossdev.install()` wraps `WDL._util.symlink_force`:
a hardlink that fails with EXDEV becomes a symlink, logged as a WARNING
(`ugc-wgw cross-device output: symlinked instead of hardlinked`) on a child
of miniwdl's logger, so it appears in the run's stderr and workflow log.
Work directories always share the file system of their `out/`, so only
links to files that outlive the run change, and those are what miniwdl
makes anyway when hardlinks are off. The installer writes
`ugc_wgw_miniwdl.pth` into the venv so the wrap is in place before a
cached workflow's outputs are linked, which happens before any plugin
loads.

## Node-local scratch

Apptainer's `--containall` clears the environment, so inside a task's
container `TMPDIR` and `SLURM_JOB_ID` are unset, and miniwdl binds the
container's `/tmp` to a directory under the task's run directory, on the
shared file system; miniwdl-slurm passes the container command to `sbatch`
as literal script arguments, so a `$TMPDIR` in `run_options` is never
expanded. Only a script running on the node can find the job's scratch
directory. `scratch.py` registers `ugc_slurm_singularity`, miniwdl-slurm's
backend with one override: when `[ugc_wgw] scratch_env` names a job
variable (`TMPDIR`, `SNIC_TMP`), the trailing `exec_script.sh` of the
`sbatch` arguments becomes `scripts/exec_scratch.sh <var> <mount> --`.
On the node that script inserts `--bind <dir>:<mount> --env TMPDIR=<mount>`
after the `exec` subcommand when the variable names a writable directory
that is not tmpfs, prints one `[ugc-wgw] scratch:` line to stderr (the
job's `slurm-<id>.out` in the task directory) and execs the command. pbmm2's
sort spill files, GLnexus's database (our override) and other
TMPDIR-aware temporaries then stay on the node's disk, which the scheduler
removes with the job. With `scratch_env` empty the backend is byte for byte
`slurm_singularity`. The log file the driver reads keeps its name
(`cli_name` is inherited).

## The GPU type

`[ugc_wgw] gpu_type` names the SLURM gres type (`l40s`, `a100`) that every
task with a `gpuCount` is submitted with; empty means any GPU of the
partition (`--gres gpu:N`). The driver passes the project's `gpu_type`
per run through miniwdl's environment override
(`MINIWDL__UGC_WGW__GPU_TYPE`), which outranks the configuration file.
Why not a workflow input: upstream's `RuntimeAttributes.gpuType` is an
input of every task, so a type given that way is part of every task's
call-cache key and switching between GPU nodes would re-align every
sample. A non-empty `gpuType` the workflow set itself (an explicit stage
input) wins over the configuration. The `ugc-wgw gpu request` NOTICE says
where the type came from (`gpu_type_from`).

