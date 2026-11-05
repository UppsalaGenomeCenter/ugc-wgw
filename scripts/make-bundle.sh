#!/usr/bin/env bash
# Build the offline code bundle for the HPC. See docs/DESIGN.md §14. Dev machine only.
#
# Usage: scripts/make-bundle.sh --version <git ref> --out <dir>
#            [--images all|ugc-wgw|data|install|none|delta] [--base <manifest.json|bundle.tar|dir>]
#            [--sif-cache <dir>] [--python-version X.Y] [--with-references <dir-with-tars>] [--keep-build]
#            [--without-ancestry]
#
# Steps:
#   1. git archive <ref> -> <out>/build/<name>/code/  (vendor/ included; VERSION inside names the bundle)
#   2. union of code/vendor/hifi-human-wgs-wdl/image_manifest.txt + code/image_manifest.ugc-wgw.txt
#      -> scripts/populate-sif-cache.sh into --sif-cache (incremental); every SIF's sha256 and size
#      (from the cache's .sha256 sidecars, computed once when absent); hardlinked into build/sif/ unless
#      --base names a previous bundle's manifest.json (or its tar, or a directory holding it) that already
#      carries a SIF of the same sha256: then the image is "required but not shipped", and the installer
#      takes it from the site's SIF store. --images delta is the full image set filtered that way (--base
#      required); --base with any other set filters the same way. tag-form references are resolved to
#      registry digests with scripts/image-digest.py;
#      --images data takes only upstream's reference data container(s), for tests of the install;
#      --images install adds the images the installer builds the derived reference files with (references.lock built[])
#   3. pip download -r code/bundle/requirements.txt (+ pip itself) for --python-version
#      into build/wheels/ (binary wheels only, manylinux x86_64; no build tools needed offline);
#      pip wheel of code/plugins/ugc_wgw_miniwdl (our miniwdl task plugin, pure python) into the same
#      directory, listed in build/wheels/requirements.txt after the PyPI pins
#   4. install-bundle.sh and requirements.txt copied to the top of the build (and next to the tar)
#   5. manifest.json (schema 2): sha256 + size of every file, versions (ours, git, upstream.lock,
#      references.lock, python, apptainer), every required image with SIF checksum, size, digest and
#      `shipped`, the base bundle when there is one, build host and date; a delta bundle also carries an
#      empty top-level DELTA-BUNDLE file (an older installer refuses it instead of installing one SIF)
#   6. tar -> <out>/ugc-pacbio-wgw-<version>[-delta].tar with <name>/ as top-level dir; .sha256 alongside;
#      the manifest is copied next to the tar as <tar name>.manifest.json (the --base of the next release)
#   7. --with-references: verify the archival reference tars (references.lock archival section, Zenodo)
#      against their md5 and copy them to <out>/ for sites that keep a tar next to the bundle; optional,
#      the installer takes the reference data from the bundled container (step 2, --images all)
#   7b. fetched references (references.lock <entry>.fetched[]): downloaded into build/references/ and verified
#      by sha256 before the manifest is written, so they travel in the tar (the 1000 Genomes somalier set of
#      ugc-wgw-qc, 83 MB); --without-ancestry leaves out the entries marked optional = ancestry
# Env: UGC_WGW_PIP (default "python3 -m pip", falling back to the repo .venv), APPTAINER_CACHEDIR (default <out>/apptainer-cache)
set -euo pipefail
cd "$(dirname "$0")/.."
repo=$(pwd)

ref=""; out=""; images="all"; base=""; sif_cache=""; pyver=""; with_refs=""; keep_build=0; without_ancestry=0
while [ $# -gt 0 ]; do
  case "$1" in
    --version) ref="$2"; shift 2 ;;
    --out) out="$2"; shift 2 ;;
    --images) images="$2"; shift 2 ;;
    --base) base="$2"; shift 2 ;;
    --sif-cache) sif_cache="$2"; shift 2 ;;
    --python-version) pyver="$2"; shift 2 ;;
    --with-references) with_refs="$2"; shift 2 ;;
    --keep-build) keep_build=1; shift ;;
    --without-ancestry) without_ancestry=1; shift ;;
    -h|--help) sed -n '2,29p' "$0"; exit 0 ;;
    *) echo "unknown option $1" >&2; exit 2 ;;
  esac
done
if [ -z "$ref" ] || [ -z "$out" ]; then echo "--version <git ref> and --out <dir> are required" >&2; exit 2; fi
case "$images" in all|ugc-wgw|data|install|none|delta) ;; *) echo "--images must be all, ugc-wgw, data, install, none or delta" >&2; exit 2 ;; esac
if [ "$images" = delta ] && [ -z "$base" ]; then echo "--images delta needs --base <manifest.json|bundle.tar|dir> of the previous bundle" >&2; exit 2; fi
if [ "$images" = none ] && [ -n "$base" ]; then echo "--base makes no sense with --images none" >&2; exit 2; fi
[ -z "$base" ] || [ -e "$base" ] || { echo "--base not found: $base" >&2; exit 2; }
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
base_manifest=""
if [ -n "$base" ]; then
  base_manifest="$out/build/base-manifest.json"
  python3 - "$base" "$base_manifest" <<'PYEOF' || exit 1
import json, os, sys, tarfile
src, out = sys.argv[1:3]
if os.path.isdir(src):
    src = os.path.join(src, "manifest.json")
if src.endswith(".tar"):
    with tarfile.open(src) as t:   # members are sorted by name: code/ first, then manifest.json, long before sif/
        member = next((m for m in t if m.name.count("/") == 1 and m.name.endswith("/manifest.json")), None)
        if member is None:
            sys.exit(f"[make-bundle] no manifest.json in {src}")
        text = t.extractfile(member).read().decode()
else:
    text = open(src).read()
doc = json.loads(text)
if "images" not in doc or "name" not in doc:
    sys.exit(f"[make-bundle] {src} is not a bundle manifest")
open(out, "w").write(text)
print(f"[make-bundle]   base: {doc['name']} ({len(doc['images'])} images, schema {doc.get('schema', 1)})", file=sys.stderr)
PYEOF
fi
log "1. git archive $ref ($commit) -> code/"
git archive --format=tar "$commit" | tar -x -C "$build/code"

# ---- 2. SIF cache --------------------------------------------------------
manifests=()
case "$images" in
  all|delta) manifests=("$build/code/vendor/hifi-human-wgs-wdl/image_manifest.txt" "$build/code/image_manifest.ugc-wgw.txt") ;;
  ugc-wgw) manifests=("$build/code/image_manifest.ugc-wgw.txt") ;;
  data) grep 'workflow-data-container' "$build/code/vendor/hifi-human-wgs-wdl/image_manifest.txt" > "$out/build/data-manifest.txt"
        manifests=("$out/build/data-manifest.txt") ;;
  install) grep 'workflow-data-container' "$build/code/vendor/hifi-human-wgs-wdl/image_manifest.txt" > "$out/build/data-manifest.txt"
        python3 -c 'import json,sys; print("\n".join(sorted({bt["image"] for b in json.load(open(sys.argv[1]))["bundles"].values() for bt in b.get("built", [])})))' \
          "$build/code/references.lock" >> "$out/build/data-manifest.txt"
        manifests=("$out/build/data-manifest.txt") ;;
esac
sif_list="$out/build/sif-list.tsv"
sif_sha="$out/build/sif-sha.tsv"
: > "$sif_list"; : > "$sif_sha"
if [ ${#manifests[@]} -gt 0 ]; then
  log "2. SIF cache ($images) in $sif_cache"
  no_https="${UGC_WGW_REGISTRY_INSECURE:-}"  # export-public: insecure registry host (empty = none)
  nh=(); [ -z "$no_https" ] || nh=(--no-https-host "$no_https")
  scripts/populate-sif-cache.sh --dir "$sif_cache" ${nh[@]+"${nh[@]}"} --list "$sif_list" "${manifests[@]}"
  # 2a. sha256 and size of every SIF: the cache's sidecar when it is newer than the file, computed once otherwise
  while IFS=$'\t' read -r image sifname; do
    f="$sif_cache/$sifname"
    if [ -f "$f.sha256" ] && [ "$f.sha256" -nt "$f" ]; then
      sha=$(cut -d' ' -f1 "$f.sha256")
    else
      sha=$(sha256sum "$f" | cut -d' ' -f1)
      echo "$sha  $sifname" > "$f.sha256"
    fi
    printf '%s\t%s\t%s\t%s\n' "$image" "$sifname" "$sha" "$(stat -c %s "$f")" >> "$sif_sha"
  done < "$sif_list"
  # 2b. ship every SIF, or only those the base bundle does not already carry (by sha256)
  base_shas=""
  [ -z "$base_manifest" ] || base_shas=$(python3 -c 'import json,sys; print("\n".join(str(i.get("sha256") or "") for i in json.load(open(sys.argv[1])).get("images", [])))' "$base_manifest")
  n_ship=0; n_skip=0; b_ship=0; b_skip=0
  while IFS=$'\t' read -r _ sifname sha size; do
    if [ -n "$base_shas" ] && grep -qx "$sha" <<< "$base_shas"; then
      n_skip=$((n_skip + 1)); b_skip=$((b_skip + size)); continue
    fi
    ln "$sif_cache/$sifname" "$build/sif/$sifname" 2>/dev/null || cp "$sif_cache/$sifname" "$build/sif/$sifname"
    n_ship=$((n_ship + 1)); b_ship=$((b_ship + size))
  done < "$sif_sha"
  if [ -n "$base_manifest" ]; then
    log "   $((n_ship + n_skip)) images required: $n_ship shipped ($(echo "$b_ship" | awk '{printf "%.2f GB", $1/1e9}')), $n_skip already in the base ($(echo "$b_skip" | awk '{printf "%.2f GB", $1/1e9}'))"
    : > "$build/DELTA-BUNDLE"
  else
    log "   $n_ship images shipped ($(echo "$b_ship" | awk '{printf "%.2f GB", $1/1e9}'))"
  fi
  log "2c. resolving registry digests"
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

# ---- 3c. fetched references (references.lock <entry>.fetched[]) ----------
log "3c. fetched references"
python3 - "$build/code/references.lock" "$build/references" "$without_ancestry" <<'PYEOF' || exit 1
import hashlib, json, os, sys, urllib.request
lock, dest, without_ancestry = sys.argv[1], sys.argv[2], sys.argv[3] == "1"
doc = json.load(open(lock))
n = 0
for bname, bundle in doc["bundles"].items():
    for f in bundle.get("fetched", []):
        if without_ancestry and f.get("optional") == "ancestry":
            print(f"[make-bundle]   {bname}: {f['name']} left out (--without-ancestry)", file=sys.stderr); continue
        os.makedirs(dest, exist_ok=True)
        path = os.path.join(dest, f["name"])
        if not os.path.exists(path):
            print(f"[make-bundle]   {bname}: fetching {f['url']}", file=sys.stderr)
            with urllib.request.urlopen(f["url"], timeout=120) as r, open(path + ".partial", "wb") as out:
                for chunk in iter(lambda: r.read(1 << 20), b""):
                    out.write(chunk)
            os.replace(path + ".partial", path)
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        if h.hexdigest() != f["sha256"] or os.path.getsize(path) != f["bytes"]:
            sys.exit(f"[make-bundle]   {bname}: {f['name']} does not match references.lock (sha256/bytes)")
        print(f"[make-bundle]   {bname}: {f['name']} verified ({f['bytes']} bytes)", file=sys.stderr); n += 1
if n == 0:
    print("[make-bundle]   none", file=sys.stderr)
PYEOF

# ---- 4. installer at the top ---------------------------------------------
cp "$build/code/scripts/install-bundle.sh" "$build/install-bundle.sh"
cp "$build/install-bundle.sh" "$out/install-bundle.sh"
chmod +x "$build/install-bundle.sh" "$out/install-bundle.sh"

# ---- 5. manifest.json ----------------------------------------------------
log "5. manifest.json"
python3 - "$build" "$name" "$version" "$ref" "$commit" "$pyver" "$sif_sha" "$out/build/digests.tsv" "$images" "$base_manifest" <<'PYEOF'
import datetime as dt, hashlib, json, os, platform, socket, subprocess, sys
build, name, version, ref, commit, pyver, sif_sha, digests_tsv, images_mode, base_manifest = sys.argv[1:11]
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
        if rel in ("manifest.json", "DELTA-BUNDLE") or os.path.islink(p):
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
if os.path.exists(sif_sha):
    for line in open(sif_sha):
        image, sifname, sha, size = line.rstrip("\n").split("\t")
        entry = by_path.get("sif/" + sifname)
        if entry is not None and (entry["sha256"] != sha or entry["size"] != int(size)):
            sys.exit(f"[make-bundle] {sifname}: the shipped file differs from the cache record")
        images.append({"image": image, "sif": sifname, "sha256": sha, "size": int(size), "shipped": entry is not None,
                       "digest": digests.get(image, image.split("@", 1)[1] if "@" in image else "unresolved")})
base = None
if base_manifest:
    bdoc = json.load(open(base_manifest))
    base = {"name": bdoc.get("name"), "version": bdoc.get("ugc_pacbio_wgw", {}).get("version"),
            "git_commit": bdoc.get("ugc_pacbio_wgw", {}).get("git_commit"), "manifest_sha256": sha256(base_manifest),
            "images": len(bdoc.get("images", []))}
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
    "schema": 2,   # 2: images carry `shipped`; `base`, `images_shipped`, `images_total` (delta bundles)
    "name": name,
    "ugc_pacbio_wgw": {"version": version, "git_ref": ref, "git_commit": commit},
    "upstream": load("upstream.lock", "sources"),
    "references": load("references.lock", "bundles"),
    "python_version": pyver,
    "wheels": wheels,
    "images_mode": images_mode,
    "images": images,
    "images_total": len(images),
    "images_shipped": sum(1 for i in images if i["shipped"]),
    "base": base,
    "build": {"host": socket.gethostname(), "date": dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
              "apptainer": run(["apptainer", "--version"]), "platform": platform.platform()},
    "files": files,
}
with open(os.path.join(build, "manifest.json"), "w") as fh:
    json.dump(doc, fh, indent=2, sort_keys=True)
    fh.write("\n")
total = sum(f["size"] for f in files)
shipped = sum(1 for i in images if i["shipped"])
print(f"[make-bundle]   {len(files)} files, {total/1e9:.2f} GB; {len(images)} images ({shipped} shipped"
      + (f", {len(images) - shipped} in the base {base['name']}" if base else "") + f"); {len(wheels)} wheels", file=sys.stderr)
PYEOF

# ---- 6. tar --------------------------------------------------------------
tarname="$name${base:+-delta}"
tarball="$out/$tarname.tar"
log "6. $tarball"
mtime=$(git show -s --format=%cI "$commit")
tar --sort=name --owner=0 --group=0 --numeric-owner --mtime="$mtime" -C "$out/build" -cf "$tarball" "$name"
(cd "$out" && sha256sum "$tarname.tar" > "$tarname.tar.sha256")
cp "$build/manifest.json" "$out/$tarname.manifest.json"
log "   $(du -h "$tarball" | cut -f1)  sha256 $(cut -d' ' -f1 "$tarball.sha256"); manifest copied to $out/$tarname.manifest.json"

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
