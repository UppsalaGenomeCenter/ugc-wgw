#!/usr/bin/env bash
# Static check of every entrypoint. Run before committing anything under workflows/.
set -euo pipefail
cd "$(dirname "$0")/.."
status=0
shopt -s nullglob
entrypoints=(workflows/ugc_wgw_*.wdl)
if [ ${#entrypoints[@]} -eq 0 ]; then
  echo "no entrypoints yet (workflows/ugc_wgw_*.wdl); nothing to check"
  exit 0
fi
for wf in "${entrypoints[@]}"; do
  echo "== miniwdl check $wf"
  miniwdl check --strict "$wf" || status=1
done
# Guard: docs/guide/workflow-graphs.md is generated from the entrypoints and must be current.
echo "== wdl-graph drift check"
if ! python3 scripts/wdl-graph.py --check; then
  echo "ERROR: docs/guide/workflow-graphs.md is stale or scripts/wdl-graph.py failed;" \
       "run: python3 scripts/wdl-graph.py --write (needs miniwdl importable: activate .venv)" >&2
  status=1
fi
# Guard: the task resource inventory (backends/hpc/resources.declared.tsv, docs/guide/task-resources.md) is
# generated from the same tree and must be current.
echo "== wdl-resources drift check"
if ! python3 scripts/wdl-resources.py --check; then
  echo "ERROR: backends/hpc/resources.declared.tsv or docs/guide/task-resources.md is stale or" \
       "scripts/wdl-resources.py failed; run: python3 scripts/wdl-resources.py --write" >&2
  status=1
fi
# Guard: nothing under vendor/ may be modified relative to HEAD.
if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  if ! git diff --quiet HEAD -- vendor/ 2>/dev/null; then
    echo "ERROR: uncommitted modifications under vendor/ (read-only; use scripts/vendor-upstream.py)" >&2
    status=1
  fi
fi
exit $status
