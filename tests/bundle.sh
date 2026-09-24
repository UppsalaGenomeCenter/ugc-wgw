#!/usr/bin/env bash
# Bundle round-trip: make-bundle.sh (--images data when the dev machine's SIF cache already holds the
# reference data container, so the installer's copy-and-verify step runs; --images none otherwise)
# -> install-bundle.sh --verify-only and a real install into a scratch prefix, then `ugc-wgw init --install`.
# Needs network for pip download (dev machine or CI). Usage: tests/bundle.sh [<scratch dir>]
set -euo pipefail
cd "$(dirname "$0")/.."
scratch=${1:-$(mktemp -d "${TMPDIR:-/tmp}/ugc-bundle-test.XXXXXX")}
echo "== bundle test in $scratch"
images=none; bundle_args=()
if ls bundle/out/sif-cache/*workflow-data-container-hifi-human-wgs-wdl-grch38_giabv3*.sif >/dev/null 2>&1; then
  images=data; bundle_args=(--sif-cache "$(pwd)/bundle/out/sif-cache")
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
printf 'SLURM_PARTITION=core\nSLURM_ACCOUNT=acc-1\nSLURM_EXTRA_ARGS=--qos short\nTASK_CPU_MAX=48\nTASK_MEMORY_MAX=360G\nTASK_RESOURCES=<prefix>/resources.tsv\nSLURM_PARTITION_GPU=gpu\nSINGULARITY_NV=1\n' > "$scratch/site.cfg"
report=$(scripts/install-bundle.sh --bundle "$tarball" --prefix "$scratch/prefix" --site "$scratch/site.cfg" --activate)
install_dir=$(echo "$report" | sed -n 's/^install_dir=//p')
[ -x "$install_dir/venv/bin/miniwdl" ]
[ -f "$install_dir/miniwdl.cfg" ] && ! grep -q '{{' "$install_dir/miniwdl.cfg"
[ "$(readlink "$scratch/prefix/current")" = "versions/$(cat VERSION)" ]
grep -q 'command_shell = /bin/bash' "$install_dir/miniwdl.cfg"
grep -q '^cpu_max = 48$' "$install_dir/miniwdl.cfg"
grep -q '^memory_max = 386547056640$' "$install_dir/miniwdl.cfg"   # 360 GiB in bytes
grep -q '^defaults = {"docker": "ubuntu:20.04", "maxRetries": 1, "time_minutes": 4320, "slurm_partition": "core", "slurm_account": "acc-1", "slurm_partition_gpu": "gpu"}$' "$install_dir/miniwdl.cfg"
grep -q '^run_options = \["--containall", "--no-mount", "hostfs", "--nv"\]$' "$install_dir/miniwdl.cfg"
grep -q '^extra_args = "--qos short"$' "$install_dir/miniwdl.cfg"
grep -q "^resources = $scratch/prefix/resources.tsv$" "$install_dir/miniwdl.cfg"
[ -f "$scratch/prefix/resources.tsv" ] && grep -q '^task	cpu	memory	time	partition	constraint$' "$scratch/prefix/resources.tsv"
echo "$report" | grep -q '^plugin=ugc-wgw-miniwdl '
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
grep -q "^scatter_regions	$scratch/prefix/references/ugc-wgw-extras-0.2.0/scatter_regions.GRCh38_GIABv3.tsv$" "$ref_map"
grep -q "^run_starphase	true$" "$ref_map"
[ -f "$scratch/prefix/references/ugc-wgw-extras-0.2.0/scatter_regions.GRCh38_GIABv3.tsv" ]
echo "$report" | grep -q "^ref_map=$ref_map$"
if [ "$images" = data ]; then
  # the tree was copied out of the container and every file verified against references.lock
  tree=$scratch/prefix/references/hifi-wdl-resources-v4.0.0-GRCh38_GIABv3
  [ -f "$tree/manifest.json" ] && [ -f "$tree/GRCh38_GIABv3/trgt/adotto_strchive_20250827.hg38.bed.gz" ]
  echo "$report" | grep -q "^references=$tree$"
  while read -r f; do [ -f "$f" ] || { echo "map names a missing file: $f" >&2; exit 1; }; done < <(cut -f2 "$ref_map" | grep '^/')
  # a second install must reuse the tree (no copy) and still verify it
  scripts/install-bundle.sh --bundle "$tarball" --prefix "$scratch/prefix" --site "$scratch/site.cfg" --force > "$scratch/reinstall.log" 2>&1
  grep -q 'reference tree present' "$scratch/reinstall.log"
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
