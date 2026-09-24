#!/usr/bin/env bash
# Build the offline code bundle for the HPC. See docs/DESIGN.md §14. Dev machine only.
#
# Usage: scripts/make-bundle.sh --version <git ref> --out <dir>
#            [--images all|ugc-wgw|data|none] [--sif-cache <dir>] [--python-version X.Y]
#            [--with-references <dir-with-tars>] [--keep-build]
#
# Steps:
#   1. git archive <ref> -> <out>/build/<name>/code/  (vendor/ included; VERSION inside names the bundle)
#   2. union of code/vendor/hifi-human-wgs-wdl/image_manifest.txt + code/image_manifest.ugc-wgw.txt
#      -> scripts/populate-sif-cache.sh into --sif-cache (incremental), hardlinked into build/sif/;
#      tag-form references are resolved to registry digests with scripts/image-digest.py;
#      --images data takes only upstream's reference data container(s), for tests of the install
#   3. pip download -r code/bundle/requirements.txt (+ pip itself) for --python-version
#      into build/wheels/ (binary wheels only, manylinux x86_64; no build tools needed offline);
#      pip wheel of code/plugins/ugc_wgw_miniwdl (our miniwdl task plugin, pure python) into the same
#      directory, listed in build/wheels/requirements.txt after the PyPI pins
#   4. install-bundle.sh and requirements.txt copied to the top of the build (and next to the tar)
#   5. manifest.json: sha256 + size of every file, versions (ours, git, upstream.lock, references.lock,
#      python, apptainer), images with SIF checksum and digest, build host and date
#   6. tar -> <out>/ugc-pacbio-wgw-<version>.tar with <name>/ as top-level dir; .sha256 alongside
#   7. --with-references: verify the archival reference tars (references.lock archival section, Zenodo)
#      against their md5 and copy them to <out>/ for sites that keep a tar next to the bundle; optional,
#      the installer takes the reference data from the bundled container (step 2, --images all)
# Env: UGC_WGW_PIP (default "python3 -m pip", falling back to the repo .venv), APPTAINER_CACHEDIR (default <out>/apptainer-cache)
set -euo pipefail
cd "$(dirname "$0")/.."
repo=$(pwd)

ref=""; out=""; images="all"; sif_cache=""; pyver=""; with_refs=""; keep_build=0
while [ $# -gt 0 ]; do
  case "$1" in
    --version) ref="$2"; shift 2 ;;
    --out) out="$2"; shift 2 ;;
    --images) images="$2"; shift 2 ;;
    --sif-cache) sif_cache="$2"; shift 2 ;;
    --python-version) pyver="$2"; shift 2 ;;
    --with-references) with_refs="$2"; shift 2 ;;
    --keep-build) keep_build=1; shift ;;
    -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "unknown option $1" >&2; exit 2 ;;
  esac
done
[ -n "$ref" ] && [ -n "$out" ] || { echo "--version <git ref> and --out <dir> are required" >&2; exit 2; }
case "$images" in all|ugc-wgw|data|none) ;; *) echo "--images must be all, ugc-wgw, data or none" >&2; exit 2 ;; esac
[ -n "$pyver" ] || pyver=$(python3 -c 'import sys; print(f"{sys.version_info[0]}.{sys.version_info[1]}")')
out=$(mkdir -p "$out" && cd "$out" && pwd)
sif_cache=${sif_cache:-$out/sif-cache}
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$out/apptainer-cache}"
mkdir -p "$APPTAINER_CACHEDIR"

log() { echo "[make-bundle] $*" >&2; }

# ---- 1. code -------------------------------------------------------------
commit=$(git rev-parse --verify "$ref^{commit}")
version=$(git show "$commit:VERSION" | tr -d '[:space:]')
[ -n "$version" ] || { echo "VERSION file missing at $ref" >&2; exit 1; }
name="ugc-pacbio-wgw-$version"
if ! git describe --tags --exact-match "$commit" 2>/dev/null | grep -qx "v$version"; then
  log "warning: $ref ($commit) is not tagged v$version; building a development bundle of $version"
fi
if [ "$(git rev-parse HEAD)" = "$commit" ] && [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  log "warning: working tree has uncommitted changes; the bundle contains the committed tree only"
fi
build="$out/build/$name"
rm -rf "$out/build"
mkdir -p "$build/code" "$build/sif" "$build/wheels"
log "1. git archive $ref ($commit) -> code/"
git archive --format=tar "$commit" | tar -x -C "$build/code"

# ---- 2. SIF cache --------------------------------------------------------
manifests=()
case "$images" in
  all) manifests=("$build/code/vendor/hifi-human-wgs-wdl/image_manifest.txt" "$build/code/image_manifest.ugc-wgw.txt") ;;
  ugc-wgw) manifests=("$build/code/image_manifest.ugc-wgw.txt") ;;
  data) grep 'workflow-data-container' "$build/code/vendor/hifi-human-wgs-wdl/image_manifest.txt" > "$out/build/data-manifest.txt"
        manifests=("$out/build/data-manifest.txt") ;;
esac
sif_list="$out/build/sif-list.tsv"
: > "$sif_list"
if [ ${#manifests[@]} -gt 0 ]; then
  log "2. SIF cache ($images) in $sif_cache"
  no_https="${UGC_WGW_REGISTRY_INSECURE:-}"  # export-public: insecure registry host (empty = none)
  nh=(); [ -z "$no_https" ] || nh=(--no-https-host "$no_https")
  scripts/populate-sif-cache.sh --dir "$sif_cache" ${nh[@]+"${nh[@]}"} --list "$sif_list" "${manifests[@]}"
  while IFS=$'\t' read -r _ sifname; do
    ln "$sif_cache/$sifname" "$build/sif/$sifname" 2>/dev/null || cp "$sif_cache/$sifname" "$build/sif/$sifname"
  done < "$sif_list"
  log "2b. resolving registry digests"
  digest_args=()
  for m in "${manifests[@]}"; do digest_args+=(--manifest "$m"); done
  scripts/image-digest.py "${digest_args[@]}" > "$out/build/digests.tsv" || log "warning: some digests unresolved"
else
  log "2. no images requested (--images none)"
  : > "$out/build/digests.tsv"
fi

# ---- 3. wheels -----------------------------------------------------------
log "3. wheels for python $pyver"
pip_cmd=${UGC_WGW_PIP:-"python3 -m pip"}
if ! $pip_cmd --version >/dev/null 2>&1; then
  pip_cmd="$repo/.venv/bin/python -m pip"
  $pip_cmd --version >/dev/null 2>&1 || { echo "no pip found; set UGC_WGW_PIP" >&2; exit 1; }
fi
abi="cp${pyver/./}"
pip_platform=(--only-binary=:all: --python-version "$pyver" --implementation cp --abi "$abi" --abi none
              --platform manylinux2014_x86_64 --platform manylinux_2_17_x86_64 --platform manylinux_2_28_x86_64
              --platform manylinux1_x86_64 --platform any)
$pip_cmd download --quiet --dest "$build/wheels" "${pip_platform[@]}" -r "$build/code/bundle/requirements.txt" >&2
$pip_cmd download --quiet --dest "$build/wheels" "${pip_platform[@]}" --no-deps pip >&2
log "3b. ugc-wgw-miniwdl plugin wheel"
$pip_cmd wheel --quiet --no-deps --wheel-dir "$build/wheels" "$build/code/plugins/ugc_wgw_miniwdl" >&2 \
  || { echo "plugin wheel build failed (pip builds it in an isolated env: needs network for setuptools)" >&2; exit 1; }
rm -rf "$build/code/plugins/ugc_wgw_miniwdl/build" "$build/code/plugins/ugc_wgw_miniwdl"/*.egg-info
plugin_whl=$(find "$build/wheels" -maxdepth 1 -name 'ugc_wgw_miniwdl-*.whl' | head -1)
[ -n "$plugin_whl" ] || { echo "plugin wheel missing after pip wheel" >&2; exit 1; }
plugin_version=$(basename "$plugin_whl" | cut -d- -f2)
{ cat "$build/code/bundle/requirements.txt"; echo "# built by make-bundle.sh from code/plugins/ugc_wgw_miniwdl"; echo "ugc-wgw-miniwdl==$plugin_version"; } > "$build/wheels/requirements.txt"

# ---- 4. installer at the top ---------------------------------------------
cp "$build/code/scripts/install-bundle.sh" "$build/install-bundle.sh"
cp "$build/install-bundle.sh" "$out/install-bundle.sh"
chmod +x "$build/install-bundle.sh" "$out/install-bundle.sh"

# ---- 5. manifest.json ----------------------------------------------------
log "5. manifest.json"
python3 - "$build" "$name" "$version" "$ref" "$commit" "$pyver" "$sif_list" "$out/build/digests.tsv" "$images" <<'PYEOF'
import datetime as dt, hashlib, json, os, platform, socket, subprocess, sys
build, name, version, ref, commit, pyver, sif_list, digests_tsv, images_mode = sys.argv[1:10]
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()
files = []
for root, dirs, names in os.walk(build):
    dirs.sort()
    for fn in sorted(names):
        p = os.path.join(root, fn)
        rel = os.path.relpath(p, build)
        if rel == "manifest.json" or os.path.islink(p):
            continue
        files.append({"path": rel, "sha256": sha256(p), "size": os.path.getsize(p)})
digests = {}
if os.path.exists(digests_tsv):
    for line in open(digests_tsv):
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2:
            digests[parts[0]] = parts[1]
images = []
by_path = {f["path"]: f for f in files}
if os.path.exists(sif_list):
    for line in open(sif_list):
        image, sifname = line.rstrip("\n").split("\t")
        entry = by_path.get("sif/" + sifname, {})
        images.append({"image": image, "sif": sifname, "sha256": entry.get("sha256"), "size": entry.get("size"),
                       "digest": digests.get(image, image.split("@", 1)[1] if "@" in image else "unresolved")})
def load(rel, key):
    p = os.path.join(build, "code", rel)
    return json.load(open(p)).get(key, {}) if os.path.exists(p) else {}
def run(cmd):
    try:
        return subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=30).stdout.strip().splitlines()[0]
    except Exception:
        return "unknown"
wheels = sorted(f["path"][len("wheels/"):] for f in files if f["path"].startswith("wheels/") and f["path"].endswith((".whl", ".tar.gz")))
doc = {
    "schema": 1,
    "name": name,
    "ugc_pacbio_wgw": {"version": version, "git_ref": ref, "git_commit": commit},
    "upstream": load("upstream.lock", "sources"),
    "references": load("references.lock", "bundles"),
    "python_version": pyver,
    "wheels": wheels,
    "images_mode": images_mode,
    "images": images,
    "build": {"host": socket.gethostname(), "date": dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
              "apptainer": run(["apptainer", "--version"]), "platform": platform.platform()},
    "files": files,
}
with open(os.path.join(build, "manifest.json"), "w") as fh:
    json.dump(doc, fh, indent=2, sort_keys=True)
    fh.write("\n")
total = sum(f["size"] for f in files)
print(f"[make-bundle]   {len(files)} files, {total/1e9:.2f} GB; {len(images)} images; {len(wheels)} wheels", file=sys.stderr)
PYEOF

# ---- 6. tar --------------------------------------------------------------
tarball="$out/$name.tar"
log "6. $tarball"
mtime=$(git show -s --format=%cI "$commit")
tar --sort=name --owner=0 --group=0 --numeric-owner --mtime="$mtime" -C "$out/build" -cf "$tarball" "$name"
(cd "$out" && sha256sum "$name.tar" > "$name.tar.sha256")
log "   $(du -h "$tarball" | cut -f1)  sha256 $(cut -d' ' -f1 "$tarball.sha256")"

# ---- 7. references -------------------------------------------------------
if [ -n "$with_refs" ]; then
  log "7. reference tars from $with_refs"
  python3 - "$build/code/references.lock" "$with_refs" "$out" <<'PYEOF'
import hashlib, json, os, shutil, sys
lock, src, out = sys.argv[1:4]
doc = json.load(open(lock))
status = 0
for bname, bundle in doc["bundles"].items():
    for f in bundle.get("archival", {}).get("files", []):
        path = os.path.join(src, f["name"])
        if not os.path.exists(path):
            print(f"[make-bundle]   {bname}: {f['name']} not found in {src}", file=sys.stderr); status = 1; continue
        algo = "md5" if "md5" in f else "sha256"
        h = hashlib.new(algo)
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        if h.hexdigest() != f[algo]:
            print(f"[make-bundle]   {bname}: {f['name']} {algo} mismatch", file=sys.stderr); status = 1; continue
        shutil.copy2(path, os.path.join(out, f["name"]))
        print(f"[make-bundle]   {bname}: {f['name']} verified and copied", file=sys.stderr)
sys.exit(status)
PYEOF
fi

[ "$keep_build" -eq 1 ] || rm -rf "$out/build"
log "done: $tarball"
echo "$tarball"
