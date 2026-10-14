#!/usr/bin/env bash
# ugc-wgw: run a task's container with the job's node-local scratch bound in and TMPDIR pointing at it.
#
# Usage: exec_scratch.sh <env var> <mount point> -- <apptainer> [--verbose] exec <options...> <image> <command...>
#
# miniwdl-slurm submits every task as `sbatch ... exec_script.sh <container command>`, and its script is
# `exec "$@"`. The ugc_slurm_singularity backend (ugc_wgw_miniwdl/scratch.py) substitutes this one when
# [ugc_wgw] scratch_env names a job environment variable. Running on the node, where that variable is known,
# it inserts `--bind <dir>:<mount>` and `--env TMPDIR=<mount>` after the `exec` subcommand when the variable
# points at a writable directory that is not tmpfs (memory); tools that honour TMPDIR then keep their
# temporaries on the node's disk. The directory is the job's own and the scheduler removes it when the job
# ends. One line on stderr, which SLURM writes to slurm-<jobid>.out in the task directory, says what happened.
set -euo pipefail

usage() { echo "usage: exec_scratch.sh <env var> <mount point> -- <container command...>" >&2; exit 2; }

[ $# -ge 4 ] || usage
var=$1
mount=$2
shift 2
[ "$1" = "--" ] || usage
shift

scratch=${!var:-}
note=""
if [ -z "$scratch" ]; then
  note="scratch: none ($var is unset)"
elif [ ! -d "$scratch" ] || [ ! -w "$scratch" ]; then
  note="scratch: none ($var=$scratch is not a writable directory)"
elif [ "$(stat -f -c %T "$scratch" 2>/dev/null || echo unknown)" = "tmpfs" ]; then
  note="scratch: none ($var=$scratch is tmpfs, memory)"
fi
if [ -n "$note" ]; then
  echo "[ugc-wgw] $note" >&2
  exec "$@"
fi

# insert after the `exec` subcommand; the executable may be followed by --verbose first
args=()
inserted=0
for a in "$@"; do
  args+=("$a")
  if [ "$inserted" -eq 0 ] && [ "$a" = "exec" ]; then
    args+=("--bind" "$scratch:$mount" "--env" "TMPDIR=$mount")
    inserted=1
  fi
done
if [ "$inserted" -eq 0 ]; then
  echo "[ugc-wgw] scratch: none (no exec subcommand in the container command)" >&2
  exec "$@"
fi
echo "[ugc-wgw] scratch: $scratch -> $mount (TMPDIR)" >&2
exec "${args[@]}"
