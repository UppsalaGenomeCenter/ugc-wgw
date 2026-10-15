#!/usr/bin/env bash
# Bundle round-trip: make-bundle.sh (--images data when the dev machine's SIF cache already holds the
# reference data container, so the installer's copy-and-verify step runs; --images none otherwise)
# -> install-bundle.sh --verify-only and a real install into a scratch prefix, then `ugc-wgw init --install`.
# Needs network for pip download (dev machine or CI). Usage: tests/bundle.sh [<scratch dir>]
# A scratch dir it made itself (about 8 GB) is removed after a passing run and kept, with its path printed, after
# a failing one or when UGC_WGW_KEEP_SCRATCH=1; a scratch dir given as the argument is always kept.
set -euo pipefail
cd "$(dirname "$0")/.."
keep=${UGC_WGW_KEEP_SCRATCH:-}
if [ $# -ge 1 ]; then
  scratch=$1; keep=1
else
  scratch=$(mktemp -d "${TMPDIR:-/tmp}/ugc-bundle-test.XXXXXX")
fi
cleanup() {
  rc=$?
  if [ "$rc" -eq 0 ] && [ -z "$keep" ]; then
    rm -rf "$scratch"
  else
    echo "== scratch kept at $scratch" >&2
  fi
}
trap cleanup EXIT
echo "== bundle test in $scratch"
images=none; bundle_args=()
if ls bundle/out/sif-cache/*workflow-data-container-hifi-human-wgs-wdl-grch38_giabv3*.sif >/dev/null 2>&1; then
  images=data; bundle_args=(--sif-cache "$(pwd)/bundle/out/sif-cache")
  # with the pbmm2 image cached too, the installer's index build runs (about a minute, 12 GB of memory)
  if ls bundle/out/sif-cache/*pbmm2@sha256_86e39aa6*.sif >/dev/null 2>&1; then images=install; fi
fi
echo "== images mode $images"
tarball=$(scripts/make-bundle.sh --version HEAD --out "$scratch/out" --images "$images" "${bundle_args[@]}")
[ -f "$tarball" ] && [ -f "$tarball.sha256" ]
scripts/install-bundle.sh --bundle "$tarball" --prefix "$scratch/prefix" --verify-only 2>&1 | tee "$scratch/verify.log" >&2
grep -q 'plugin: ugc-wgw-miniwdl' "$scratch/verify.log"
grep -q 'reference data container:' "$scratch/verify.log"
# a site file that keeps --partition in SLURM_EXTRA_ARGS next to SLURM_PARTITION must be refused before anything is installed
printf 'SLURM_PARTITION=core\nSLURM_EXTRA_ARGS=--partition core\n' > "$scratch/bad-site.cfg"
if scripts/install-bundle.sh --bundle "$tarball" --prefix "$scratch/prefix-bad" --site "$scratch/bad-site.cfg" >/dev/null 2>&1; then
  echo "installer accepted --partition in SLURM_EXTRA_ARGS next to SLURM_PARTITION" >&2; exit 1
fi
[ ! -d "$scratch/prefix-bad/versions" ]
printf 'SLURM_PARTITION=core\nSLURM_ACCOUNT=acc-1\nSLURM_EXTRA_ARGS=--qos short\nTASK_CPU_MAX=48\nTASK_MEMORY_MAX=360G\nTASK_RESOURCES=<prefix>/resources.tsv\nSLURM_PARTITION_GPU=gpu\nSINGULARITY_NV=1\nTASK_SCRATCH_ENV=TMPDIR\n' > "$scratch/site.cfg"
report=$(scripts/install-bundle.sh --bundle "$tarball" --prefix "$scratch/prefix" --site "$scratch/site.cfg" --activate)
install_dir=$(echo "$report" | sed -n 's/^install_dir=//p')
[ -x "$install_dir/venv/bin/miniwdl" ]
echo "$report" | grep -q '^profiles='
[ -f "$scratch/prefix/profiles/cpu.json" ] && [ -f "$scratch/prefix/profiles/parabricks.json" ]
[ -f "$install_dir/miniwdl.cfg" ] && ! grep -q '{{' "$install_dir/miniwdl.cfg"
[ "$(readlink "$scratch/prefix/current")" = "versions/$(cat VERSION)" ]
grep -q 'command_shell = /bin/bash' "$install_dir/miniwdl.cfg"
grep -q '^cpu_max = 48$' "$install_dir/miniwdl.cfg"
grep -q '^memory_max = 386547056640$' "$install_dir/miniwdl.cfg"   # 360 GiB in bytes
grep -q '^defaults = {"docker": "ubuntu:20.04", "maxRetries": 1, "time_minutes": 4320, "slurm_partition": "core", "slurm_account": "acc-1", "slurm_partition_gpu": "gpu"}$' "$install_dir/miniwdl.cfg"
grep -q '^run_options = \["--containall", "--no-mount", "hostfs", "--nv"\]$' "$install_dir/miniwdl.cfg"
grep -q '^extra_args = "--qos short"$' "$install_dir/miniwdl.cfg"
grep -q "^resources = $scratch/prefix/resources.tsv$" "$install_dir/miniwdl.cfg"
grep -q '^container_backend = ugc_slurm_singularity$' "$install_dir/miniwdl.cfg"
grep -q '^scratch_env = TMPDIR$' "$install_dir/miniwdl.cfg"
grep -q '^scratch_mount = /mnt/scratch$' "$install_dir/miniwdl.cfg"
[ -f "$scratch/prefix/resources.tsv" ] && grep -q '^task	cpu	memory	time	partition	constraint$' "$scratch/prefix/resources.tsv"
echo "$report" | grep -q '^plugin=ugc-wgw-miniwdl '
# an invalid scratch variable name is refused before anything is installed
# shellcheck disable=SC2016  # a literal dollar on purpose: variable names are not expanded by the installer
echo 'TASK_SCRATCH_ENV=$TMPDIR' > "$scratch/bad-scratch.cfg"
if scripts/install-bundle.sh --bundle "$tarball" --prefix "$scratch/prefix-bad2" --site "$scratch/bad-scratch.cfg" >/dev/null 2>&1; then
  echo "installer accepted TASK_SCRATCH_ENV=\$TMPDIR" >&2; exit 1
fi
[ ! -d "$scratch/prefix-bad2/versions" ]
# the node-local scratch backend: registered, selected by the rendered cfg, its exec script shipped in the wheel and
# working: a fake apptainer shows the bind and TMPDIR inserted after `exec`
"$install_dir/venv/bin/python" - "$install_dir/miniwdl.cfg" "$scratch" <<'PYEOF'
import logging, os, pathlib, stat, subprocess, sys
from WDL.runtime.config import Loader
from WDL.runtime import task_container
from ugc_wgw_miniwdl import scratch
cfg = Loader(logging.getLogger("bundle-test"), filenames=[sys.argv[1]])
assert cfg["scheduler"]["container_backend"] == "ugc_slurm_singularity"
host = pathlib.Path(sys.argv[2]) / "backend-host"; host.mkdir()
backend = task_container.new(cfg, logging.getLogger("bundle-test"), "call-x-1", str(host))
assert type(backend).__name__ == "SlurmSingularityScratch", type(backend)
args = backend._slurm_invocation()
assert args[-4:] == [scratch.EXEC_SCRIPT, "TMPDIR", "/mnt/scratch", "--"], args[-4:]
assert os.path.isfile(scratch.EXEC_SCRIPT) and "site-packages" in scratch.EXEC_SCRIPT, scratch.EXEC_SCRIPT
fake = pathlib.Path(sys.argv[2]) / "fake-apptainer"
fake.write_text('#!/usr/bin/env bash\nfor a in "$@"; do printf "%s\\n" "$a"; done\n'); fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
job_scratch = pathlib.Path(sys.argv[2]) / "job-scratch"; job_scratch.mkdir()
res = subprocess.run(["bash", scratch.EXEC_SCRIPT, "TMPDIR", "/mnt/scratch", "--", str(fake), "exec", "--pwd", "/w", "img", "cmd"],
                     env={"PATH": os.environ["PATH"], "TMPDIR": str(job_scratch)}, capture_output=True, text=True)
assert res.returncode == 0, res.stderr
assert res.stdout.splitlines()[:5] == ["exec", "--bind", f"{job_scratch}:/mnt/scratch", "--env", "TMPDIR=/mnt/scratch"], res.stdout   # the fake prints its arguments
assert "scratch:" in res.stderr, res.stderr
PYEOF
# the hardlink-else-symlink wrap is in place at interpreter start (ugc_wgw_miniwdl.pth), before any plugin loads
ls "$install_dir"/venv/lib/python*/site-packages/ugc_wgw_miniwdl.pth >/dev/null
"$install_dir/venv/bin/python" -c 'import sys, WDL._util as u; sys.exit(0 if getattr(u.symlink_force, "ugc_wgw_crossdev", False) else 1)'
# the rendered cfg parses as miniwdl config and the plugin sees the policy path through it
"$install_dir/venv/bin/python" - "$install_dir/miniwdl.cfg" "$scratch/prefix/resources.tsv" <<'PYEOF'
import logging, sys
from WDL.runtime.config import Loader
from ugc_wgw_miniwdl import resources
cfg = Loader(logging.getLogger("t"), filenames=[sys.argv[1]])
assert resources.policy_path(cfg) == sys.argv[2], resources.policy_path(cfg)
assert cfg["task_runtime"].get_int("cpu_max") == 48
assert cfg["task_runtime"].get_dict("defaults")["slurm_partition"] == "core"
assert cfg["task_runtime"].get_dict("defaults")["slurm_partition_gpu"] == "gpu"
assert cfg["singularity"].get_list("run_options") == ["--containall", "--no-mount", "hostfs", "--nv"]
PYEOF
ref_map=$install_dir/references/ugc_wgw_ref_map.GRCh38_GIABv3.tsv
grep -q "^name	GRCh38_GIABv3$" "$ref_map"
grep -q "^scatter_regions	$scratch/prefix/references/ugc-wgw-extras-0.3.0/scatter_regions.GRCh38_GIABv3.tsv$" "$ref_map"
mmi=$scratch/prefix/references/ugc-wgw-extras-0.3.0/human_GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.mmi
grep -q "^pbmm2_index	$mmi$" "$ref_map"
grep -q "^run_starphase	true$" "$ref_map"
[ -f "$scratch/prefix/references/ugc-wgw-extras-0.3.0/scatter_regions.GRCh38_GIABv3.tsv" ]
echo "$report" | grep -q "^ref_map=$ref_map$"
if [ "$images" = data ] || [ "$images" = install ]; then
  # the tree was copied out of the container and every file verified against references.lock
  tree=$scratch/prefix/references/hifi-wdl-resources-v4.0.0-GRCh38_GIABv3
  [ -f "$tree/manifest.json" ] && [ -f "$tree/GRCh38_GIABv3/trgt/adotto_strchive_20250827.hg38.bed.gz" ]
  echo "$report" | grep -q "^references=$tree$"
  if [ "$images" = install ]; then
    # the reference's pbmm2 index was built with the bundled image and matches references.lock
    [ -f "$mmi" ]
    [ "$(md5sum "$mmi" | cut -d' ' -f1)" = "$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["bundles"]["ugc-wgw-extras"]["built"][0]["md5"])' references.lock)" ]
    while read -r f; do [ -f "$f" ] || { echo "map names a missing file: $f" >&2; exit 1; }; done < <(cut -f2 "$ref_map" | grep '^/')
  else
    [ ! -f "$mmi" ]   # no pbmm2 image in the bundle: the installer says how to build it and the map still names it
  fi
  # a second install must reuse the tree (no copy), verify it and keep the index
  scripts/install-bundle.sh --bundle "$tarball" --prefix "$scratch/prefix" --site "$scratch/site.cfg" --force > "$scratch/reinstall.log" 2>&1
  grep -q 'reference tree present' "$scratch/reinstall.log"
  grep -q 'profile kept' "$scratch/reinstall.log"
  if [ "$images" = install ]; then grep -q 'pbmm2 index present and verified' "$scratch/reinstall.log"; else grep -q 'no pbmm2 index built' "$scratch/reinstall.log"; fi
fi
if grep -rq 'v3p1p0\|tertiary' "$install_dir/references" "$install_dir/inputs"; then echo "v3 map names survive in the install" >&2; exit 1; fi
# WDL read_map() needs exactly two fields per line: no comments, no blank lines, in any rendered map
for m in "$install_dir"/references/*.tsv; do
  if grep -qE '^#|^$' "$m" || awk -F'\t' 'NF != 2 { bad = 1 } END { exit bad }' "$m"; then :; else echo "rendered map is not a two-column TSV: $m" >&2; exit 1; fi
done
"$install_dir/code/bin/ugc-wgw" init "$scratch/proj" --install "$scratch/prefix/current" --ref-map "$ref_map"
grep -q '"miniwdl_cfg"' "$scratch/proj/.ugc-wgw/config.json"
# the rendered inputs templates point at files the installer really wrote
for f in $(grep -ho '"[^"]*/references/[^"]*"' "$install_dir"/inputs/*.json | tr -d '"' | sort -u); do
  [ -f "$f" ] || { echo "rendered template points at a missing file: $f" >&2; exit 1; }
done
echo "== bundle test passed ($(du -h "$tarball" | cut -f1) bundle)"
