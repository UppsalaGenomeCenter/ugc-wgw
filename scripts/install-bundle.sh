#!/usr/bin/env bash
# Install a code bundle on the isolated HPC. See docs/DESIGN.md §14.
#
# Usage: install-bundle.sh --bundle <ugc-pacbio-wgw-X.Y.Z.tar> --prefix <root>
#            [--references <dir>] [--site <site.cfg>] [--python <exe>] [--activate] [--verify-only] [--force]
#
# Steps:
#   1. verify the tar's .sha256 (if present), unpack into <root>/.staging/, verify every manifest.json entry
#   2. move to <root>/versions/<version>/{code,sif,wheels}; python3 -m venv --without-pip venv;
#      bootstrap pip from the bundled wheel; pip install --no-index --find-links wheels -r requirements.txt
#   3. render code/backends/hpc/miniwdl.cfg.template -> <root>/versions/<version>/miniwdl.cfg
#      ({{SIF_CACHE_DIR}}, {{CALL_CACHE_DIR}}, {{TASK_CONCURRENCY}}, {{SLURM_EXTRA_ARGS}},
#      {{TASK_RUNTIME_DEFAULTS}} from TASK_TIME_MINUTES + SLURM_PARTITION + SLURM_ACCOUNT (+ the _GPU pair),
#      {{SINGULARITY_RUN_OPTIONS}} from SINGULARITY_NV,
#      {{TASK_CPU_MAX}}, {{TASK_MEMORY_MAX}}, {{TASK_RESOURCES}}; site values from --site KEY=VALUE
#      lines, see backends/hpc/site.cfg.example); create <root>/resources.tsv from the example when
#      absent, validate it when present
#   4. references: copy the build's tree and manifest.json out of the bundled reference data
#      container SIF into <references>/<install_subdir>/ (--references, default <root>/references)
#      unless present, verify every md5 in references.lock, install ugc-wgw-extras (scatter regions),
#      render references/ugc_wgw_ref_map.<build>.template.tsv into <root>/versions/<version>/references/
#      and check it against the manifest, render the inputs templates into .../inputs/
#   5. --verify-only: step 1, a temporary venv, the ugc_wgw_resources plugin entry point, `miniwdl check
#      --strict` on every entrypoint, whether the reference data container is in the bundle; nothing installed
#   6. --activate: <root>/current -> versions/<version> (never done implicitly)
#   7. print the paths the driver needs and the `ugc-wgw init` command line
# Needs only bash, coreutils, tar, python3 (matching the bundle's python_version), apptainer for step 4. No network.
set -euo pipefail

bundle=""; prefix=""; references=""; site=""; python_exe=""; activate=0; verify_only=0; force=0
while [ $# -gt 0 ]; do
  case "$1" in
    --bundle) bundle="$2"; shift 2 ;;
    --prefix) prefix="$2"; shift 2 ;;
    --references) references="$2"; shift 2 ;;
    --site) site="$2"; shift 2 ;;
    --python) python_exe="$2"; shift 2 ;;
    --activate) activate=1; shift ;;
    --verify-only) verify_only=1; shift ;;
    --force) force=1; shift ;;
    -h|--help) sed -n '2,19p' "$0"; exit 0 ;;
    *) echo "unknown option $1" >&2; exit 2 ;;
  esac
done
if [ -z "$bundle" ] || [ -z "$prefix" ]; then echo "--bundle <tar> and --prefix <root> are required" >&2; exit 2; fi
[ -f "$bundle" ] || { echo "bundle not found: $bundle" >&2; exit 1; }
python_exe=${python_exe:-python3}
command -v "$python_exe" >/dev/null || { echo "python not found: $python_exe" >&2; exit 1; }

log() { echo "[install-bundle] $*" >&2; }
die() { echo "[install-bundle] error: $*" >&2; exit 1; }

mkdir -p "$prefix"
prefix=$(cd "$prefix" && pwd)
references=${references:-$prefix/references}

# ---- 1. verify and unpack ----------------------------------------------
if [ -f "$bundle.sha256" ]; then
  log "1. verifying $(basename "$bundle").sha256"
  (cd "$(dirname "$bundle")" && sha256sum --check --status "$(basename "$bundle").sha256") || die "bundle checksum mismatch"
else
  log "1. no .sha256 next to the bundle; relying on manifest.json only"
fi
staging="$prefix/.staging"
rm -rf "$staging"
mkdir -p "$staging"
tar -C "$staging" -xf "$bundle"
name=$(ls "$staging")
if [ -z "$name" ] || [ ! -d "$staging/$name" ]; then die "unexpected bundle layout"; fi
src="$staging/$name"
[ -f "$src/manifest.json" ] || die "manifest.json missing"
log "   verifying every file in manifest.json"
"$python_exe" - "$src" <<'PYEOF'
import hashlib, json, os, sys
src = sys.argv[1]
doc = json.load(open(os.path.join(src, "manifest.json")))
bad = []
seen = set()
for f in doc["files"]:
    p = os.path.join(src, f["path"])
    seen.add(f["path"])
    if not os.path.isfile(p):
        bad.append(f"missing: {f['path']}"); continue
    if os.path.getsize(p) != f["size"]:
        bad.append(f"size: {f['path']}"); continue
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    if h.hexdigest() != f["sha256"]:
        bad.append(f"sha256: {f['path']}")
for root, _, names in os.walk(src):
    for fn in names:
        rel = os.path.relpath(os.path.join(root, fn), src)
        if rel != "manifest.json" and rel not in seen:
            bad.append(f"not in manifest: {rel}")
want = doc["python_version"]
have = f"{sys.version_info[0]}.{sys.version_info[1]}"
if want != have:
    bad.append(f"python {have} does not match the bundle's wheels for python {want} (rebuild with make-bundle.sh --python-version {have})")
if bad:
    print("\n".join("[install-bundle]   " + b for b in bad), file=sys.stderr)
    sys.exit(1)
print(f"[install-bundle]   {len(doc['files'])} files OK; {len(doc.get('images', []))} images; version {doc['ugc_pacbio_wgw']['version']}", file=sys.stderr)
PYEOF
version=$("$python_exe" -c 'import json,sys; print(json.load(open(sys.argv[1]))["ugc_pacbio_wgw"]["version"])' "$src/manifest.json")

make_venv() {  # $1 = venv dir, $2 = wheels dir
  "$python_exe" -m venv --without-pip "$1"
  pipwhl=$(find "$2" -maxdepth 1 -name 'pip-*.whl' | head -1)
  [ -n "$pipwhl" ] || die "no pip wheel in the bundle"
  "$1/bin/python" "$pipwhl/pip" install --quiet --no-index --find-links "$2" pip
  "$1/bin/python" -m pip install --quiet --no-index --find-links "$2" -r "$2/requirements.txt"
}

plugin_version() {  # $1 = venv dir; prints the ugc-wgw-miniwdl version, fails when the task plugin is not registered
  "$1/bin/python" - <<'PYEOF'
import importlib.metadata as m, sys
eps = m.entry_points()
eps = eps.select(group="miniwdl.plugin.task") if hasattr(eps, "select") else eps.get("miniwdl.plugin.task", [])
hit = [e for e in eps if e.name == "ugc_wgw_resources"]
if not hit or hit[0].value != "ugc_wgw_miniwdl.resources:task" or not callable(hit[0].load()):
    sys.exit(1)
print(m.version("ugc-wgw-miniwdl"))
PYEOF
}

# ---- 5. verify-only ------------------------------------------------------
if [ "$verify_only" -eq 1 ]; then
  log "5. --verify-only: temporary venv and miniwdl check --strict"
  make_venv "$staging/venv" "$src/wheels"
  status=0
  for wf in "$src"/code/workflows/ugc_wgw_*.wdl; do
    if "$staging/venv/bin/miniwdl" check --strict "$wf" >/dev/null 2>&1; then
      log "   ok   $(basename "$wf")"
    else
      log "   FAIL $(basename "$wf")"; status=1
    fi
  done
  log "   engine: $("$staging/venv/bin/miniwdl" --version | head -1)"
  if "$python_exe" -c 'import json,sys; m=json.load(open(sys.argv[1])); r=json.load(open(sys.argv[2]))["bundles"]["hifi-wdl-resources"]; sys.exit(0 if any(i.get("image")==r["source"] for i in m.get("images",[])) else 1)' "$src/manifest.json" "$src/code/references.lock"; then
    log "   reference data container: in the bundle"
  else
    log "   reference data container: not in the bundle (images mode $("$python_exe" -c 'import json,sys; print(json.load(open(sys.argv[1])).get("images_mode"))' "$src/manifest.json")); the install needs an existing reference tree"
  fi
  if pv=$(plugin_version "$staging/venv"); then
    log "   plugin: ugc-wgw-miniwdl $pv (task plugin ugc_wgw_resources registered)"
  else
    log "   FAIL ugc_wgw_resources task plugin not registered in the venv"; status=1
  fi
  rm -rf "$staging"
  [ "$status" -eq 0 ] || die "verify-only failed"
  log "verify-only passed for $name"
  exit 0
fi

# ---- site values (read before anything is installed, so a bad site file costs nothing) ------
TASK_CONCURRENCY=50; CALL_CACHE_DIR="$prefix/call_cache"; SLURM_EXTRA_ARGS=""; TASK_TIME_MINUTES=4320
SLURM_PARTITION=""; SLURM_ACCOUNT=""; TASK_CPU_MAX=0; TASK_MEMORY_MAX=0; TASK_RESOURCES="$prefix/resources.tsv"
SLURM_PARTITION_GPU=""; SLURM_ACCOUNT_GPU=""; SINGULARITY_NV=0
if [ -n "$site" ]; then
  [ -f "$site" ] || die "site file not found: $site"
  while IFS='=' read -r key value; do
    key=$(echo "$key" | tr -d '[:space:]')
    case "$key" in
      ''|'#'*) ;;
      TASK_CONCURRENCY|CALL_CACHE_DIR|SLURM_EXTRA_ARGS|TASK_TIME_MINUTES|SLURM_PARTITION|SLURM_ACCOUNT|TASK_CPU_MAX|TASK_MEMORY_MAX|TASK_RESOURCES|SLURM_PARTITION_GPU|SLURM_ACCOUNT_GPU|SINGULARITY_NV)
        printf -v "$key" '%s' "${value//<prefix>/$prefix}" ;;
      *) log "   warning: unknown site key $key ignored" ;;
    esac
  done < "$site"
fi
if { [ -n "$SLURM_PARTITION" ] || [ -n "$SLURM_ACCOUNT" ]; } && echo " $SLURM_EXTRA_ARGS" | grep -qE ' (-p|--partition|-A|--account)([ =]|$)'; then
  die "SLURM_EXTRA_ARGS carries --partition/--account while SLURM_PARTITION/SLURM_ACCOUNT are set; sbatch keeps the last one it sees, so a per-task partition would never apply. Move them out of SLURM_EXTRA_ARGS."
fi
case "$TASK_CPU_MAX" in ''|*[!0-9-]*) die "TASK_CPU_MAX must be an integer (0 = no cap), got '$TASK_CPU_MAX'" ;; esac
case "$SINGULARITY_NV" in 0|1) ;; *) die "SINGULARITY_NV must be 0 or 1, got '$SINGULARITY_NV'" ;; esac
case "$TASK_TIME_MINUTES" in ''|*[!0-9]*) die "TASK_TIME_MINUTES must be a positive integer, got '$TASK_TIME_MINUTES'" ;; esac

# ---- 2. place and venv ---------------------------------------------------
dest="$prefix/versions/$version"
if [ -e "$dest" ]; then
  [ "$force" -eq 1 ] || die "$dest exists (use --force to replace)"
  rm -rf "$dest"
fi
mkdir -p "$prefix/versions"
mv "$src" "$dest"
rm -rf "$staging"
log "2. venv from bundled wheels"
make_venv "$dest/venv" "$dest/wheels"
plugin=$(plugin_version "$dest/venv") || die "ugc_wgw_resources task plugin not registered in the venv (wheel ugc_wgw_miniwdl missing from the bundle?)"
log "   $("$dest/venv/bin/miniwdl" --version | head -1); miniwdl-slurm $("$dest/venv/bin/python" -c 'import importlib.metadata as m; print(m.version("miniwdl-slurm"))'); ugc-wgw-miniwdl $plugin"

# ---- 3. miniwdl.cfg and the resource policy ---------------------------------
log "3. miniwdl.cfg"
mkdir -p "$CALL_CACHE_DIR"
"$python_exe" - "$dest/code/backends/hpc/miniwdl.cfg.template" "$dest/miniwdl.cfg" \
  "SIF_CACHE_DIR=$dest/sif" "CALL_CACHE_DIR=$CALL_CACHE_DIR" "TASK_CONCURRENCY=$TASK_CONCURRENCY" \
  "SLURM_EXTRA_ARGS=$SLURM_EXTRA_ARGS" "TASK_TIME_MINUTES=$TASK_TIME_MINUTES" \
  "SLURM_PARTITION=$SLURM_PARTITION" "SLURM_ACCOUNT=$SLURM_ACCOUNT" "TASK_CPU_MAX=$TASK_CPU_MAX" \
  "TASK_MEMORY_MAX=$TASK_MEMORY_MAX" "TASK_RESOURCES=$TASK_RESOURCES" \
  "SLURM_PARTITION_GPU=$SLURM_PARTITION_GPU" "SLURM_ACCOUNT_GPU=$SLURM_ACCOUNT_GPU" "SINGULARITY_NV=$SINGULARITY_NV" <<'PYEOF'
import json, os, sys
tpl, out, *pairs = sys.argv[1:]
kv = dict(pair.split("=", 1) for pair in pairs)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(tpl)), "..", "..", "plugins", "ugc_wgw_miniwdl"))
from ugc_wgw_miniwdl.policy import parse_memory  # binary units, as the policy file and the WDL use them
defaults = {"docker": "ubuntu:20.04", "maxRetries": 1, "time_minutes": int(kv.pop("TASK_TIME_MINUTES"))}
for key, site_key in (("slurm_partition", "SLURM_PARTITION"), ("slurm_account", "SLURM_ACCOUNT"),
                      ("slurm_partition_gpu", "SLURM_PARTITION_GPU"), ("slurm_account_gpu", "SLURM_ACCOUNT_GPU")):
    value = kv.pop(site_key).strip()
    if value:
        defaults[key] = value   # the *_gpu keys are what miniwdl-slurm uses for tasks with a gpuCount
kv["TASK_RUNTIME_DEFAULTS"] = json.dumps(defaults)  # one line: configparser continuation rules
run_options = ["--containall", "--no-mount", "hostfs"]
if kv.pop("SINGULARITY_NV") == "1":
    run_options.append("--nv")   # bind the host's NVIDIA driver into every container; a no-op on nodes without one
kv["SINGULARITY_RUN_OPTIONS"] = json.dumps(run_options)
mm = kv["TASK_MEMORY_MAX"].strip()
try:
    kv["TASK_MEMORY_MAX"] = mm if mm in ("0", "-1") else str(parse_memory(mm))
except ValueError as exc:
    sys.exit(f"TASK_MEMORY_MAX: {exc}")
text = open(tpl).read()
for k, v in kv.items():
    text = text.replace("{{" + k + "}}", v)
left = [l for l in text.splitlines() if "{{" in l and not l.lstrip().startswith("#")]
if left:
    sys.exit("unrendered placeholders: " + "; ".join(left))
open(out, "w").write(text)
PYEOF
if [ -e "$TASK_RESOURCES" ]; then
  "$python_exe" - "$dest/code/plugins/ugc_wgw_miniwdl" "$TASK_RESOURCES" <<'PYEOF' || die "the resource policy is malformed (fix it, then rerun with --force)"
import sys
sys.path.insert(0, sys.argv[1])
from ugc_wgw_miniwdl.policy import PolicyError, parse_policy
try:
    rules = parse_policy(sys.argv[2])
except PolicyError as exc:
    print(f"[install-bundle]   {exc}", file=sys.stderr)
    sys.exit(1)
print(f"[install-bundle]   resource policy {sys.argv[2]}: {len(rules)} row(s)", file=sys.stderr)
PYEOF
else
  mkdir -p "$(dirname "$TASK_RESOURCES")"
  cp "$dest/code/backends/hpc/resources.tsv.example" "$TASK_RESOURCES"
  log "   resource policy created from the example (no active rows): $TASK_RESOURCES"
fi

# ---- 4. references and inputs templates ---------------------------------
log "4. references under $references"
mkdir -p "$dest/references" "$dest/inputs"
# references.lock names the build, the container paths, the install subdirs and every file's md5;
# manifest.json names the SIF of the reference data container when the bundle carries it
ref_build=""; ref_subdir=""; ref_sif=""; ref_data_path=""; ref_manifest_path=""; extras_subdir=""; extras_files=""
eval "$("$python_exe" - "$dest/code/references.lock" "$dest/manifest.json" <<'PYEOF'
import json, shlex, sys
b = json.load(open(sys.argv[1]))["bundles"]
r, x = b["hifi-wdl-resources"], b["ugc-wgw-extras"]
m = json.load(open(sys.argv[2]))
sif = next((i["sif"] for i in m.get("images", []) if i.get("image") == r["source"]), "")
for k, v in {"ref_build": r["build"], "ref_subdir": r["install_subdir"], "ref_sif": sif,
             "ref_data_path": r["container_paths"]["data"], "ref_manifest_path": r["container_paths"]["manifest"],
             "extras_subdir": x["install_subdir"], "extras_files": " ".join(f["name"] for f in x["files"])}.items():
    print(f"{k}={shlex.quote(v)}")
PYEOF
)"
tree=$references/$ref_subdir
if [ -f "$tree/manifest.json" ]; then
  log "   reference tree present: $tree"
elif [ -n "$ref_sif" ] && [ -f "$dest/sif/$ref_sif" ]; then
  command -v apptainer >/dev/null || die "apptainer is needed to copy the reference data out of $ref_sif"
  { mkdir -p "$tree" 2>/dev/null && [ -w "$tree" ]; } || die "cannot create $tree (use --references <writable dir>)"
  log "   copying $ref_data_path and $ref_manifest_path out of the reference data container ($ref_sif; a few minutes)"
  apptainer exec --bind "$tree:/mnt/ugc-references" "$dest/sif/$ref_sif" \
    sh -c "cp -r '$ref_data_path' /mnt/ugc-references/ && cp '$ref_manifest_path' /mnt/ugc-references/" \
    || die "copying the reference data out of the container failed"
else
  log "   warning: $tree is absent and the bundle carries no reference data container (make-bundle.sh --images all or data);"
  log "   copy the container's $ref_data_path and $ref_manifest_path there yourself (docs/guide/04-install.md)"
fi
if [ -f "$tree/manifest.json" ]; then
  "$python_exe" - "$dest/code/references.lock" "$tree" <<'PYEOF' || die "reference tree $tree does not match references.lock"
import hashlib, json, os, sys
r = json.load(open(sys.argv[1]))["bundles"]["hifi-wdl-resources"]
bad = []
for f in r["files"]:
    p = os.path.join(sys.argv[2], f["name"])
    if not os.path.isfile(p):
        bad.append("missing: " + f["name"]); continue
    if os.path.getsize(p) != f["bytes"]:
        bad.append("size: " + f["name"]); continue
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    if h.hexdigest() != f["md5"]:
        bad.append("md5: " + f["name"])
if bad:
    print("\n".join("[install-bundle]   " + b for b in bad), file=sys.stderr)
    sys.exit(1)
print(f"[install-bundle]   {len(r['files'])} reference files verified (md5) under {sys.argv[2]}", file=sys.stderr)
PYEOF
fi
if mkdir -p "$references/$extras_subdir" 2>/dev/null && [ -w "$references/$extras_subdir" ]; then
  for f in $extras_files; do cp "$dest/code/references/$f" "$references/$extras_subdir/"; done
  log "   ugc-wgw-extras installed to $references/$extras_subdir"
else
  log "   warning: $references not writable; copy code/references/{$extras_files} to $references/$extras_subdir/ yourself"
fi
# WDL read_map() rejects comment lines, so the template's comments are dropped at render time
ref_map=$dest/references/ugc_wgw_ref_map.$ref_build.tsv
sed -e '/^#/d' -e '/^$/d' -e "s|<local_path_prefix>|$references|g" "$dest/code/references/ugc_wgw_ref_map.$ref_build.template.tsv" > "$ref_map"
if [ -f "$tree/manifest.json" ]; then
  # the rendered map must say exactly what the container's manifest says (files under the tree, scalars as strings)
  "$python_exe" - "$ref_map" "$tree" "$ref_data_path" <<'PYEOF' || die "rendered map $ref_map disagrees with $tree/manifest.json"
import json, os, sys
ref_map, tree, data_path = sys.argv[1:4]
have = dict(line.rstrip("\n").split("\t", 1) for line in open(ref_map) if line.strip())
keys = {"ref_name": "name", "ref_fasta": "fasta", "ref_index": "fasta_index"}
bad = []
for e in json.load(open(os.path.join(tree, "manifest.json")))["entries"]:
    key = keys.get(e["name"], e["name"])
    if e["type"] == "File":
        want = os.path.join(tree, os.path.relpath(e["value"][0], os.path.dirname(data_path)))
    else:
        want = json.dumps(e["value"]) if isinstance(e["value"], bool) else str(e["value"])
    if have.get(key) != want:
        bad.append(f"{key}: map has {have.get(key)!r}, manifest says {want!r}")
if bad:
    print("\n".join("[install-bundle]   " + b for b in bad), file=sys.stderr)
    sys.exit(1)
PYEOF
fi
registry_default="ghcr.io/uppsalagenomecenter"  # export-public: registry default
for tpl in "$dest"/code/backends/hpc/*.inputs.template.json; do
  sed -e "s|<local_path_prefix>|$dest/references|g" -e "s|<ugc_wgw_version>|$version|g" \
      -e "s|<ugc_wgw_container_registry>|$registry_default|g" "$tpl" > "$dest/inputs/$(basename "${tpl%.template.json}.json")"
done
missing=0
while read -r f; do
  [ -e "$f" ] || missing=$((missing + 1))
done < <(cut -f2 "$ref_map" | grep '^/')
if [ "$missing" -eq 0 ]; then
  log "   all reference files of $ref_map present"
else
  log "   warning: $missing reference file(s) named by $ref_map not (yet) present under $references"
fi

# ---- 6. activate ---------------------------------------------------------
if [ "$activate" -eq 1 ]; then
  ln -sfn "versions/$version" "$prefix/current.tmp" && mv -T "$prefix/current.tmp" "$prefix/current"
  log "6. $prefix/current -> versions/$version"
else
  log "6. not activated (use --activate to point $prefix/current at this version)"
fi

# ---- 7. report -----------------------------------------------------------
cat <<REPORT
install_dir=$dest
code=$dest/code
miniwdl=$dest/venv/bin/miniwdl
cfg=$dest/miniwdl.cfg
sif_cache=$dest/sif
references=$tree
ref_map=$ref_map
inputs_templates=$dest/inputs
ugc_wgw=$dest/code/bin/ugc-wgw
resources=$TASK_RESOURCES
plugin=ugc-wgw-miniwdl $plugin
# next: $dest/code/bin/ugc-wgw init <project_dir> --install $dest --ref-map $ref_map
REPORT
