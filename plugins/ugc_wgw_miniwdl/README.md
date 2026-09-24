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
