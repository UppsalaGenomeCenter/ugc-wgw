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
