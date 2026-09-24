#!/usr/bin/env bash
#SBATCH --job-name ugc-cohort2026
#SBATCH --partition core
#SBATCH --account <account>
#SBATCH --cpus-per-task 2
#SBATCH --mem 8G
#SBATCH --time 7-00:00:00
#SBATCH --output ugc-submit-%j.log
# Runs the driver as a long-lived job on a partition that allows jobs to submit jobs.
# The driver itself is light: each miniwdl process it launches is one python process
# that waits on sbatch --wait; the heavy work happens in the tasks' own jobs.
set -euo pipefail
# module load python/3.12.7          # the interpreter the bundle's wheels and venv were made for
export PATH=/proj/ugc/current/code/bin:$PATH
cd /proj/ugc/projects/cohort2026
# scancel of this job sends SIGTERM: the driver forwards it to its miniwdl children,
# which cancel their SLURM steps, and marks the runs cancelled (exit 130).
exec ugc-wgw -v submit --mode standalone --cohort C1 --max-inflight 40 --poll-interval 60
