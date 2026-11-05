#!/usr/bin/env bash
# Install a code bundle on the isolated HPC. See docs/DESIGN.md §14.
#
# Usage: install-bundle.sh --bundle <ugc-pacbio-wgw-X.Y.Z[-delta].tar> --prefix <root>
#            [--references <dir>] [--site <site.cfg>] [--python <exe>] [--activate] [--verify-only] [--force]
#            [--sif-source <dir>]...
#        install-bundle.sh --prefix <root> --dedupe-sif
#
# Steps:
#   1. verify the tar's .sha256 (if present), unpack into <root>/.staging-<pid>/, verify every manifest.json entry
#   1b. the SIF store <root>/sif/<sha256>.sif, one copy of each image for every installed version: every image the
#      manifest requires but the bundle does not ship (a delta bundle) must be in the store or adoptable by name and
#      sha256 from <root>/current, the other versions/*/sif or a --sif-source directory (hardlinked into the store,
#      copied across file systems); otherwise the install is refused before anything is placed (the HPC never pulls)
#   2. move to <root>/versions/<version>/{code,sif,wheels}; shipped SIFs are linked into the store, the others
#      linked from it (versions/<v>/sif/ holds hardlinks: miniwdl's image_cache stays per version);
#      python3 -m venv --without-pip venv; bootstrap pip from the bundled wheel;
#      pip install --no-index --find-links wheels -r requirements.txt
#   3. render code/backends/hpc/miniwdl.cfg.template -> <root>/versions/<version>/miniwdl.cfg
#      ({{SIF_CACHE_DIR}}, {{CALL_CACHE_DIR}}, {{TASK_CONCURRENCY}}, {{SLURM_EXTRA_ARGS}},
#      {{TASK_RUNTIME_DEFAULTS}} from TASK_TIME_MINUTES + SLURM_PARTITION + SLURM_ACCOUNT (+ the _GPU pair),
#      {{SINGULARITY_RUN_OPTIONS}} from SINGULARITY_NV,
#      {{TASK_CPU_MAX}}, {{TASK_MEMORY_MAX}}, {{TASK_RESOURCES}}, {{TASK_SCRATCH_ENV}}; site values from --site KEY=VALUE
#      lines, see backends/hpc/site.cfg.example); create <root>/resources.tsv from the example when
#      absent, validate it when present
#   4. references: copy the build's tree and manifest.json out of the bundled reference data
#      container SIF into <references>/<install_subdir>/ (--references, default <root>/references)
#      unless present, verify every md5 in references.lock, install the other lock entries' files (ugc-wgw-extras,
#      ugc-wgw-asm) and build their derived files there once with the bundled images (the reference's pbmm2
#      index, the assembly stages' decoy-free fasta), verified by md5,
#      render references/ugc_wgw_ref_map.<build>.template.tsv into <root>/versions/<version>/references/
#      and check it against the manifest, render the inputs templates into .../inputs/
#   5. --verify-only: step 1, a temporary venv, the ugc_wgw_resources plugin entry point, `miniwdl check
#      --strict` on every entrypoint, whether the reference data container is in the bundle, and for a delta
#      bundle whether every unshipped image is in the store or adoptable (sha256 verified); nothing installed
#   6. --activate: <root>/current -> versions/<version> (never done implicitly)
#   7. print the paths the driver needs and the `ugc-wgw init` command line
#   --dedupe-sif (no --bundle): replace every versions/*/sif/*.sif that is its own copy by a hardlink into the
#      store (adopting it when the store lacks its content); prints the bytes freed
# Needs only bash, coreutils, tar, python3 (matching the bundle's python_version), apptainer for step 4. No network.
set -euo pipefail

bundle=""; prefix=""; references=""; site=""; python_exe=""; activate=0; verify_only=0; force=0; dedupe=0
sif_sources=()
while [ $# -gt 0 ]; do
  case "$1" in
    --bundle) bundle="$2"; shift 2 ;;
    --prefix) prefix="$2"; shift 2 ;;
    --references) references="$2"; shift 2 ;;
    --site) site="$2"; shift 2 ;;
    --python) python_exe="$2"; shift 2 ;;
    --sif-source) sif_sources+=("$2"); shift 2 ;;
    --activate) activate=1; shift ;;
    --verify-only) verify_only=1; shift ;;
    --force) force=1; shift ;;
    --dedupe-sif) dedupe=1; shift ;;
    -h|--help) sed -n '2,32p' "$0"; exit 0 ;;
    *) echo "unknown option $1" >&2; exit 2 ;;
  esac
done
[ -n "$prefix" ] || { echo "--prefix <root> is required" >&2; exit 2; }
if [ "$dedupe" -eq 0 ] && [ -z "$bundle" ]; then echo "--bundle <tar> and --prefix <root> are required (or --prefix with --dedupe-sif)" >&2; exit 2; fi
[ -z "$bundle" ] || [ -f "$bundle" ] || { echo "bundle not found: $bundle" >&2; exit 1; }
python_exe=${python_exe:-python3}
command -v "$python_exe" >/dev/null || { echo "python not found: $python_exe" >&2; exit 1; }

log() { echo "[install-bundle] $*" >&2; }
die() { echo "[install-bundle] error: $*" >&2; exit 1; }

# The SIF store: <prefix>/sif/<sha256>.sif, one file per distinct image content, shared by every installed
# version; versions/<v>/sif/<miniwdl name>.sif are hardlinks into it (a relative symlink where the file system
# refuses hardlinks). Modes: resolve (adopt what the bundle does not ship; refuse on anything missing),
# dry-run (--verify-only: report, write nothing), link (step 2b), dedupe (--dedupe-sif).
read -r -d '' SIF_STORE_PY <<'PYEOF' || true
import hashlib, json, os, shutil, sys, time
mode = sys.argv[1]
def log(msg):
    print("[install-bundle]   " + msg, file=sys.stderr)
def sha256_of(path):
    side = path + ".sha256"
    try:
        if os.path.isfile(side) and os.path.getmtime(side) >= os.path.getmtime(path):
            tok = open(side).read().split()
            if tok and len(tok[0]) == 64:
                return tok[0], "sidecar"
    except OSError:
        pass
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest(), "sha256"
def gb(n):
    return f"{n / 1e9:.2f} GB"
def index_row(store, sha, size, sif, image, note):
    try:
        with open(os.path.join(store, "INDEX.tsv"), "a") as fh:
            if fh.tell() == 0:
                fh.write("# sha256\tsize\tsif\timage\tadded_by\tdate\n")
            fh.write("\t".join([sha, str(size), sif, image, note, time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())]) + "\n")
    except OSError as exc:
        log(f"warning: INDEX.tsv not updated: {exc}")
def place(src, dst):
    """Link src to dst (same file system), else copy through a partial file verified while copying; dst must not exist."""
    try:
        os.link(src, dst)
        return "linked"
    except FileExistsError:
        return "exists"
    except OSError:
        pass
    tmp = f"{dst}.partial-{os.getpid()}"
    h = hashlib.sha256()
    with open(src, "rb") as fi, open(tmp, "wb") as fo:
        for chunk in iter(lambda: fi.read(1 << 22), b""):
            h.update(chunk); fo.write(chunk)
    want = os.path.basename(dst)[:-4] if dst.endswith(".sif") else ""
    if want and len(want) == 64 and h.hexdigest() != want:
        os.unlink(tmp)
        raise SystemExit(f"[install-bundle] error: {src} changed while being copied")
    os.replace(tmp, dst)
    return "copied"
def version_key(name):
    out = []
    for tok in name.replace("-", ".").split("."):
        out.append((0, int(tok)) if tok.isdigit() else (1, tok))
    return out

if mode in ("resolve", "dry-run"):
    manifest, store, prefix, version, plan_out, *sources = sys.argv[2:]
    doc = json.load(open(manifest))
    cand_dirs = []
    cur = os.path.join(prefix, "current")
    if os.path.islink(cur):
        cand_dirs.append(os.path.join(os.path.realpath(cur), "sif"))
    vdir = os.path.join(prefix, "versions")
    if os.path.isdir(vdir):
        for v in sorted(os.listdir(vdir), key=version_key, reverse=True):
            cand_dirs.append(os.path.join(vdir, v, "sif"))
    cand_dirs += sources
    seen = set(); dirs = []
    for d in cand_dirs:
        d = os.path.realpath(d)
        if d not in seen and os.path.isdir(d):
            seen.add(d); dirs.append(d)
    rows = []; n_store = n_adopt = n_ship = n_missing = 0; b_missing = 0; t0 = time.time(); hashed = 0
    for img in doc.get("images", []):
        name, sha, size, shipped = img["sif"], img.get("sha256"), img.get("size"), img.get("shipped", True)
        if shipped:
            rows.append((name, sha, size, "shipped", "bundle")); n_ship += 1; continue
        sp = os.path.join(store, sha + ".sif")
        if os.path.isfile(sp) and os.path.getsize(sp) == size:
            if mode == "dry-run":
                got, how = sha256_of(sp); hashed += size if how == "sha256" else 0
                if got != sha:
                    rows.append((name, sha, size, "corrupt", sp)); n_missing += 1; b_missing += size
                    log(f"{name}: the store's {os.path.basename(sp)} does not hash to its name"); continue
            rows.append((name, sha, size, "store", sp)); n_store += 1; continue
        found = None
        for d in dirs:
            c = os.path.join(d, name)
            if os.path.isfile(c) and os.path.getsize(c) == size:
                got, how = sha256_of(c); hashed += size if how == "sha256" else 0
                if got == sha:
                    found = c; break
        if found is None:
            rows.append((name, sha, size, "missing", "")); n_missing += 1; b_missing += size; continue
        if mode == "dry-run":
            rows.append((name, sha, size, "adoptable", found)); n_adopt += 1; continue
        os.makedirs(store, exist_ok=True)
        how = place(found, sp)
        index_row(store, sha, size, name, img.get("image", ""), f"{version} (from {os.path.dirname(found)})")
        rows.append((name, sha, size, "adopted", found)); n_adopt += 1
    with open(plan_out, "w") as fh:
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    elapsed = time.time() - t0
    where = f" ({gb(hashed)} hashed in {elapsed:.0f} s)" if hashed else ""
    if mode == "dry-run":
        log(f"images: {len(rows)} required; {n_ship} shipped; {n_store} in the SIF store {store}; {n_adopt} adoptable from the installed versions"
            + (" or --sif-source" if sources else "") + f"; {n_missing} missing{where}")
    else:
        log(f"1b. SIF store {store}: {n_ship} shipped by the bundle, {n_store} already in the store, {n_adopt} adopted from the installed versions"
            + (" or --sif-source" if sources else "") + f", {n_missing} missing{where}")
    if n_missing:
        base = doc.get("base") or {}
        log(f"error: {n_missing} image(s) this bundle needs are neither in it nor on this site ({store}, versions/*/sif"
            + (", --sif-source" if sources else "") + f"), {gb(b_missing)}:")
        for name, sha, size, status, _src in rows:
            if status in ("missing", "corrupt"):
                log(f"  {name} ({gb(size)}){' (corrupt in the store)' if status == 'corrupt' else ''}")
        if base:
            log(f"  This is a delta bundle (base {base.get('name')}). Install the full bundle (make-bundle.sh --images all) or a delta built against this site's manifest; the HPC never pulls.")
        sys.exit(1)

elif mode == "link":
    dest, store, plan, version = sys.argv[2:6]
    os.makedirs(store, exist_ok=True)
    added = linked = symlinked = 0; b_added = 0
    for line in open(plan):
        name, sha, size, status, _src = line.rstrip("\n").split("\t")
        size = int(size)
        dsif = os.path.join(dest, "sif", name)
        sp = os.path.join(store, sha + ".sif")
        if status == "shipped":
            how = place(dsif, sp)
            if how == "exists":   # the store already had this content: keep one copy, relink the version's file to it
                tmp = dsif + ".relink"
                os.link(sp, tmp); os.replace(tmp, dsif)
            else:
                added += 1; b_added += size
                index_row(store, sha, size, name, "", f"{version} (shipped)")
        else:
            try:
                os.link(sp, dsif); linked += 1
            except OSError:
                os.symlink(os.path.relpath(sp, os.path.dirname(dsif)), dsif); symlinked += 1
    for line in open(plan):
        name, sha, size, status, _src = line.rstrip("\n").split("\t")
        p = os.path.join(dest, "sif", name)
        if not os.path.isfile(p) or os.path.getsize(p) != int(size):
            sys.exit(f"[install-bundle] error: {p} is missing or has the wrong size after linking")
    n_store = len([f for f in os.listdir(store) if f.endswith(".sif")])
    b_store = sum(os.path.getsize(os.path.join(store, f)) for f in os.listdir(store) if f.endswith(".sif"))
    log(f"2b. SIF store {store}: {added} added ({gb(b_added)}), {linked} linked from the store"
        + (f", {symlinked} symlinked (no hardlinks on this file system: the store must stay)" if symlinked else "")
        + f"; {n_store} images, {gb(b_store)} in the store")

elif mode == "dedupe":
    prefix, store = sys.argv[2:4]
    os.makedirs(store, exist_ok=True)
    vdir = os.path.join(prefix, "versions")
    freed = adopted = relinked = 0
    for v in sorted(os.listdir(vdir)) if os.path.isdir(vdir) else []:
        sdir = os.path.join(vdir, v, "sif")
        if not os.path.isdir(sdir):
            continue
        for name in sorted(os.listdir(sdir)):
            p = os.path.join(sdir, name)
            if not name.endswith(".sif") or os.path.islink(p) or not os.path.isfile(p) or os.stat(p).st_nlink > 1:
                continue
            sha, _how = sha256_of(p)
            sp = os.path.join(store, sha + ".sif")
            if os.path.isfile(sp):
                if os.path.getsize(sp) != os.path.getsize(p):
                    log(f"warning: {sp} has the size of neither {p}; skipped"); continue
                tmp = p + ".relink"
                try:
                    os.link(sp, tmp); os.replace(tmp, p)
                except OSError as exc:
                    log(f"warning: {p}: {exc}"); continue
                freed += os.path.getsize(sp); relinked += 1
            else:
                how = place(p, sp)
                index_row(store, sha, os.path.getsize(p), name, "", f"{v} (dedupe)")
                adopted += 1
    log(f"dedupe: {relinked} file(s) replaced by store hardlinks ({gb(freed)} freed), {adopted} adopted into {store}")
PYEOF
fs_id() { stat -c %d "$1" 2>/dev/null || echo "?"; }
md5_is() { [ "$(md5sum "$1" | cut -d' ' -f1)" = "$2" ]; }
warn_fs() {  # $1 label, $2 path, $3 label, $4 path: warn when the two paths are on different file systems
  local a b
  a=$(fs_id "$2"); b=$(fs_id "$4")
  if [ "$a" != "?" ] && [ "$b" != "?" ] && [ "$a" != "$b" ]; then
    log "   warning: $1 ($2) and $3 ($4) are on different file systems; outputs that miniwdl would hardlink across them become symlinks (guide chapter 04)"
  fi
}

mkdir -p "$prefix"
prefix=$(cd "$prefix" && pwd)
references=${references:-$prefix/references}
store="$prefix/sif"

if [ "$dedupe" -eq 1 ]; then
  [ -z "$bundle" ] || die "--dedupe-sif runs on its own (no --bundle)"
  log "dedupe: hashing every versions/*/sif/*.sif that is its own copy under $prefix"
  "$python_exe" -c "$SIF_STORE_PY" dedupe "$prefix" "$store"
  exit 0
fi

# ---- 1. verify and unpack ----------------------------------------------
if [ -f "$bundle.sha256" ]; then
  log "1. verifying $(basename "$bundle").sha256"
  (cd "$(dirname "$bundle")" && sha256sum --check --status "$(basename "$bundle").sha256") || die "bundle checksum mismatch"
else
  log "1. no .sha256 next to the bundle; relying on manifest.json only"
fi
staging="$prefix/.staging-$$"
rm -rf "$staging"
mkdir -p "$staging"
trap 'rm -rf "$staging"' EXIT
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
if doc.get("schema", 1) > 2:
    bad.append(f"bundle manifest schema {doc['schema']} is newer than this installer understands (2); use the installer shipped with the bundle")
for img in doc.get("images", []):
    shipped = img.get("shipped", True)
    present = any(f["path"] == "sif/" + img["sif"] for f in doc["files"])
    if shipped and not present:
        bad.append(f"image marked shipped but absent from files: {img['sif']}")
    if not shipped and present:
        bad.append(f"image marked not shipped but present in files: {img['sif']}")
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
        if rel not in ("manifest.json", "DELTA-BUNDLE") and rel not in seen:
            bad.append(f"not in manifest: {rel}")
want = doc["python_version"]
have = f"{sys.version_info[0]}.{sys.version_info[1]}"
if want != have:
    bad.append(f"python {have} does not match the bundle's wheels for python {want} (rebuild with make-bundle.sh --python-version {have})")
if bad:
    print("\n".join("[install-bundle]   " + b for b in bad), file=sys.stderr)
    sys.exit(1)
imgs = doc.get("images", [])
shipped = sum(1 for i in imgs if i.get("shipped", True))
base = doc.get("base") or {}
print(f"[install-bundle]   {len(doc['files'])} files OK; {len(imgs)} images ({shipped} shipped"
      + (f", {len(imgs) - shipped} to be taken from the site's SIF store; delta on {base.get('name')}" if base else "")
      + f"); version {doc['ugc_pacbio_wgw']['version']}", file=sys.stderr)
PYEOF
version=$("$python_exe" -c 'import json,sys; print(json.load(open(sys.argv[1]))["ugc_pacbio_wgw"]["version"])' "$src/manifest.json")
sif_plan="$staging/sif-plan.tsv"

make_venv() {  # $1 = venv dir, $2 = wheels dir
  "$python_exe" -m venv --without-pip "$1"
  pipwhl=$(find "$2" -maxdepth 1 -name 'pip-*.whl' | head -1)
  [ -n "$pipwhl" ] || die "no pip wheel in the bundle"
  "$1/bin/python" "$pipwhl/pip" install --quiet --no-index --find-links "$2" pip
  "$1/bin/python" -m pip install --quiet --no-index --find-links "$2" -r "$2/requirements.txt"
}

plugin_version() {  # $1 = venv dir; prints the ugc-wgw-miniwdl version, fails when the task plugin, the
                    # ugc_slurm_singularity backend or its exec script is not in place
  "$1/bin/python" - <<'PYEOF'
import importlib.metadata as m, os, sys
eps = m.entry_points()
def group(name):
    return eps.select(group=name) if hasattr(eps, "select") else eps.get(name, [])
hit = [e for e in group("miniwdl.plugin.task") if e.name == "ugc_wgw_resources"]
if not hit or hit[0].value != "ugc_wgw_miniwdl.resources:task" or not callable(hit[0].load()):
    sys.exit(1)
hit = [e for e in group("miniwdl.plugin.container_backend") if e.name == "ugc_slurm_singularity"]
if not hit or hit[0].value != "ugc_wgw_miniwdl.scratch:SlurmSingularityScratch" or not callable(hit[0].load()):
    sys.exit(1)
import ugc_wgw_miniwdl.scratch as scratch
if not os.path.isfile(scratch.EXEC_SCRIPT) or scratch.SlurmSingularity is None:
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
    log "   plugin: ugc-wgw-miniwdl $pv (task plugin ugc_wgw_resources and container backend ugc_slurm_singularity registered)"
  else
    log "   FAIL ugc-wgw-miniwdl incomplete in the venv (task plugin, container backend or exec script missing)"; status=1
  fi
  if "$python_exe" -c 'import json,sys; sys.exit(0 if all(i.get("shipped", True) for i in json.load(open(sys.argv[1])).get("images", [])) else 1)' "$src/manifest.json"; then
    log "   images: all $("$python_exe" -c 'import json,sys; print(len(json.load(open(sys.argv[1])).get("images", [])))' "$src/manifest.json") shipped"
  elif ! "$python_exe" -c "$SIF_STORE_PY" dry-run "$src/manifest.json" "$store" "$prefix" "$version" "$sif_plan" ${sif_sources[@]+"${sif_sources[@]}"}; then
    log "   FAIL image(s) not shipped and not on this site; a delta bundle needs a site that already holds them (install the full bundle there; on the dev machine verify with --sif-source <sif cache>)"; status=1
  fi
  rm -rf "$staging"
  [ "$status" -eq 0 ] || die "verify-only failed"
  log "verify-only passed for $name"
  exit 0
fi

# ---- site values (read before anything is installed, so a bad site file costs nothing) ------
TASK_CONCURRENCY=50; CALL_CACHE_DIR="$prefix/call_cache"; SLURM_EXTRA_ARGS=""; TASK_TIME_MINUTES=4320
SLURM_PARTITION=""; SLURM_ACCOUNT=""; TASK_CPU_MAX=0; TASK_MEMORY_MAX=0; TASK_RESOURCES="$prefix/resources.tsv"
SLURM_PARTITION_GPU=""; SLURM_ACCOUNT_GPU=""; SINGULARITY_NV=0; TASK_SCRATCH_ENV=""
if [ -n "$site" ]; then
  [ -f "$site" ] || die "site file not found: $site"
  while IFS='=' read -r key value; do
    key=$(echo "$key" | tr -d '[:space:]')
    case "$key" in
      ''|'#'*) ;;
      TASK_CONCURRENCY|CALL_CACHE_DIR|SLURM_EXTRA_ARGS|TASK_TIME_MINUTES|SLURM_PARTITION|SLURM_ACCOUNT|TASK_CPU_MAX|TASK_MEMORY_MAX|TASK_RESOURCES|SLURM_PARTITION_GPU|SLURM_ACCOUNT_GPU|SINGULARITY_NV|TASK_SCRATCH_ENV)
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
if [ -n "$TASK_SCRATCH_ENV" ] && ! printf '%s' "$TASK_SCRATCH_ENV" | grep -qE '^[A-Za-z_][A-Za-z0-9_]*$'; then
  die "TASK_SCRATCH_ENV must be an environment variable name (TMPDIR, SNIC_TMP) or empty, got '$TASK_SCRATCH_ENV'"
fi

# ---- 1b. SIF store: every image the bundle needs but does not ship must be on this site ------
dest="$prefix/versions/$version"
if [ -e "$dest" ] && [ "$force" -ne 1 ]; then die "$dest exists (use --force to replace)"; fi
"$python_exe" -c "$SIF_STORE_PY" resolve "$src/manifest.json" "$store" "$prefix" "$version" "$sif_plan" ${sif_sources[@]+"${sif_sources[@]}"} \
  || die "the bundle cannot be installed here (see above)"

# ---- 2. place and venv ---------------------------------------------------
if [ -e "$dest" ]; then
  rm -rf "$dest"
fi
mkdir -p "$prefix/versions"
mv "$sif_plan" "$staging/sif-plan.keep.tsv"
mv "$src" "$dest"
"$python_exe" -c "$SIF_STORE_PY" link "$dest" "$store" "$staging/sif-plan.keep.tsv" "$version" || die "linking the SIFs failed"
rm -rf "$staging"
log "2. venv from bundled wheels"
make_venv "$dest/venv" "$dest/wheels"
plugin=$(plugin_version "$dest/venv") || die "ugc-wgw-miniwdl is not complete in the venv: task plugin ugc_wgw_resources, container backend ugc_slurm_singularity and its exec script are all required (wheel ugc_wgw_miniwdl missing from the bundle?)"
log "   $("$dest/venv/bin/miniwdl" --version | head -1); miniwdl-slurm $("$dest/venv/bin/python" -c 'import importlib.metadata as m; print(m.version("miniwdl-slurm"))'); ugc-wgw-miniwdl $plugin"
# outputs on another file system (a cached run elsewhere, a passed-through input) are symlinked instead of
# failing with "Invalid cross-device link": the plugin's crossdev module wraps miniwdl's link helper, and this
# .pth imports it at interpreter start so a cached workflow's outputs (linked before any plugin loads) are covered
site_packages=$("$dest/venv/bin/python" -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')
echo "import ugc_wgw_miniwdl.crossdev" > "$site_packages/ugc_wgw_miniwdl.pth"
"$dest/venv/bin/python" -c 'import sys, WDL._util as u; sys.exit(0 if getattr(u.symlink_force, "ugc_wgw_crossdev", False) else 1)' \
  || die "ugc_wgw_miniwdl.crossdev is not active in the venv (ugc_wgw_miniwdl.pth)"
log "   outputs: hardlinked, symlinked across file systems (ugc_wgw_miniwdl.pth)"

# ---- 3. miniwdl.cfg and the resource policy ---------------------------------
log "3. miniwdl.cfg"
mkdir -p "$CALL_CACHE_DIR"
log "   node-local scratch: ${TASK_SCRATCH_ENV:-off}${TASK_SCRATCH_ENV:+ (bound at /mnt/scratch, TMPDIR inside every task)}"
"$python_exe" - "$dest/code/backends/hpc/miniwdl.cfg.template" "$dest/miniwdl.cfg" \
  "SIF_CACHE_DIR=$dest/sif" "CALL_CACHE_DIR=$CALL_CACHE_DIR" "TASK_CONCURRENCY=$TASK_CONCURRENCY" \
  "SLURM_EXTRA_ARGS=$SLURM_EXTRA_ARGS" "TASK_TIME_MINUTES=$TASK_TIME_MINUTES" \
  "SLURM_PARTITION=$SLURM_PARTITION" "SLURM_ACCOUNT=$SLURM_ACCOUNT" "TASK_CPU_MAX=$TASK_CPU_MAX" \
  "TASK_MEMORY_MAX=$TASK_MEMORY_MAX" "TASK_RESOURCES=$TASK_RESOURCES" \
  "SLURM_PARTITION_GPU=$SLURM_PARTITION_GPU" "SLURM_ACCOUNT_GPU=$SLURM_ACCOUNT_GPU" "SINGULARITY_NV=$SINGULARITY_NV" \
  "TASK_SCRATCH_ENV=$TASK_SCRATCH_ENV" <<'PYEOF'
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
# Apptainer exports APPTAINER_BIND and SINGULARITY_BIND into every container with the whole bind list; miniwdl
# binds every input file, and past ~500 files that one string exceeds the kernel's 128 KB limit per environment
# string, so the container's /bin/sh cannot start ("argument list too long"; an 800-sample cohort_call, 2026-10-26).
# Emptying both lets the cohort tasks run; nothing of ours reads them (docs/guide/04-install.md).
run_options = ["--containall", "--no-mount", "hostfs", "--env", "APPTAINER_BIND=", "--env", "SINGULARITY_BIND="]
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
# site profiles for `ugc-wgw init --profile NAME`: the examples are installed once next to the policy and never
# overwritten, so a site's edits survive upgrades
mkdir -p "$prefix/profiles"
for f in "$dest"/code/backends/hpc/profiles/*.json; do
  [ -f "$f" ] || continue
  if [ -e "$prefix/profiles/$(basename "$f")" ]; then
    log "   profile kept: $prefix/profiles/$(basename "$f")"
  else
    cp "$f" "$prefix/profiles/$(basename "$f")"
    log "   profile installed from the example: $prefix/profiles/$(basename "$f")"
  fi
done

# ---- 4. references and inputs templates ---------------------------------
log "4. references under $references"
mkdir -p "$dest/references" "$dest/inputs"
# references.lock names the build, the container paths, the install subdirs and every file's md5;
# manifest.json names the SIF of the reference data container when the bundle carries it
ref_build=""; ref_subdir=""; ref_sif=""; ref_data_path=""; ref_manifest_path=""
# plan of the other lock entries (ugc-wgw-extras, ugc-wgw-asm, ...): the files they ship from code/references/ and the
# files the installer builds there; one record per line, fields separated by the unit separator (commands hold spaces)
plan=$dest/references-plan.txt
eval "$("$python_exe" - "$dest/code/references.lock" "$dest/manifest.json" "$plan" <<'PYEOF'
import json, shlex, sys
b = json.load(open(sys.argv[1]))["bundles"]
r = b["hifi-wdl-resources"]
m = json.load(open(sys.argv[2]))
images = {i.get("image"): i["sif"] for i in m.get("images", [])}
with open(sys.argv[3], "w") as out:
    for key, x in b.items():
        if key == "hifi-wdl-resources":
            continue
        out.write("\x1f".join(["files", key, x["install_subdir"], " ".join(f["name"] for f in x["files"])]) + "\n")
        for ft in x.get("fetched", []):   # downloaded into the bundle's references/ by make-bundle.sh (absent with --without-ancestry)
            out.write("\x1f".join(["fetched", key, x["install_subdir"], ft["name"], ft["sha256"], ft.get("optional", "")]) + "\n")
        for bt in x.get("built", []):   # built by this installer, not shipped
            out.write("\x1f".join(["built", key, x["install_subdir"], bt["name"], bt["md5"], bt["image"], images.get(bt["image"], ""),
                                   bt["source_file"], bt.get("resources", ""), bt["command"]]) + "\n")
for k, v in {"ref_build": r["build"], "ref_subdir": r["install_subdir"], "ref_sif": images.get(r["source"], ""),
             "ref_data_path": r["container_paths"]["data"], "ref_manifest_path": r["container_paths"]["manifest"]}.items():
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
# shipped files: code/references/<name> -> <references>/<install_subdir>/ (ugc-wgw-extras, ugc-wgw-asm, ...)
while IFS=$'\x1f' read -r kind entry subdir names _; do
  [ "$kind" = files ] || continue
  if mkdir -p "$references/$subdir" 2>/dev/null && [ -w "$references/$subdir" ]; then
    for f in $names; do cp "$dest/code/references/$f" "$references/$subdir/"; done
    log "   $entry installed to $references/$subdir ($names)"
  else
    log "   warning: $references not writable; copy code/references/{$names} to $references/$subdir/ yourself"
  fi
done < "$plan"
# fetched files (references.lock <entry>.fetched: the 1000 Genomes somalier set): copied from the bundle's references/ and
# verified by sha256; a bundle built with --without-ancestry lacks them, which the map then names as absent paths
while IFS=$'\x1f' read -r kind entry subdir name sha optional; do
  [ "$kind" = fetched ] || continue
  target=$references/$subdir/$name
  if [ -f "$target" ] && [ "$(sha256sum "$target" | cut -d' ' -f1)" = "$sha" ]; then
    log "   $name present and verified ($entry)"
  elif [ -f "$dest/references/$name" ]; then
    [ "$(sha256sum "$dest/references/$name" | cut -d' ' -f1)" = "$sha" ] || die "$dest/references/$name does not match references.lock (sha256 $sha)"
    if mkdir -p "$references/$subdir" 2>/dev/null && [ -w "$references/$subdir" ]; then
      cp "$dest/references/$name" "$target"
      log "   $name installed to $references/$subdir ($entry)"
    else
      log "   warning: $references/$subdir not writable; copy $dest/references/$name there yourself"
    fi
  else
    log "   warning: $name not in the bundle (make-bundle.sh --without-ancestry${optional:+; optional: $optional}): the map names it, the stage that needs it skips that step"
  fi
done < "$plan"
# built files (references.lock <entry>.built: the reference's pbmm2 index, the assembly stages' decoy-free fasta): built once
# per references directory with the bundled image and the lock's command ({source}, {out} and {dir} substituted, paths as
# seen inside the container), verified against the lock's md5, then renamed from the temporary name. Some tools insist on
# the suffix (pbmm2: .mmi), so the temporary name keeps it: <stem>.building.<suffix>.
while IFS=$'\x1f' read -r kind entry subdir name md5 image sif source resources command; do
  [ "$kind" = built ] || continue
  target=$references/$subdir/$name
  tmp_name=${name%.*}.building.${name##*.}
  cmd=${command//\{source\}//mnt/ugc-references/$source}
  cmd=${cmd//\{out\}//mnt/ugc-references/$subdir/$tmp_name}
  cmd=${cmd//\{dir\}//mnt/ugc-references/$subdir}
  if [ -f "$target" ] && md5_is "$target" "$md5"; then
    log "   $name present and verified ($entry)"
  elif [ -f "$target" ]; then
    die "$target does not match references.lock (md5 $md5); remove it and run the installer again to rebuild it"
  elif [ ! -f "$references/$source" ]; then
    log "   warning: $name not built: its source $references/$source is not present"
  elif [ -z "$sif" ] || [ ! -f "$dest/sif/$sif" ]; then
    log "   warning: $name not built: the bundle carries no $image. Build it with that image (/mnt/ugc-references = $references;"
    log "   the md5 must be $md5), then run the installer again: $cmd"
  elif [ ! -w "$references/$subdir" ]; then
    log "   warning: $name not built: $references/$subdir is not writable"
  else
    command -v apptainer >/dev/null || die "apptainer is needed to build $name with $sif"
    log "   building $name with $sif ($resources)"
    rm -f "$references/$subdir/$tmp_name"
    if ! apptainer exec --bind "$references:/mnt/ugc-references" "$dest/sif/$sif" sh -c "$cmd" > "$dest/build-$name.log" 2>&1; then
      rm -f "$references/$subdir/$tmp_name"
      die "building $name failed (log: $dest/build-$name.log). If the login node killed it for memory, run the same apptainer command in a job (sbatch --mem=16G --cpus-per-task=8), check the md5 ($md5) and rename $tmp_name to $name"
    fi
    md5_is "$references/$subdir/$tmp_name" "$md5" || { rm -f "$references/$subdir/$tmp_name"; die "the built $name does not match references.lock (md5 $md5)"; }
    mv "$references/$subdir/$tmp_name" "$target"
    log "   $name built and verified: $target"
  fi
done < "$plan"
# WDL read_map() rejects comment lines, so the template's comments are dropped at render time
warn_fs "install prefix" "$prefix" "references" "$tree"
warn_fs "install prefix" "$prefix" "call cache" "$CALL_CACHE_DIR"
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
sif_store=$store
references=$tree
ref_map=$ref_map
inputs_templates=$dest/inputs
ugc_wgw=$dest/code/bin/ugc-wgw
resources=$TASK_RESOURCES
profiles=$prefix/profiles
plugin=ugc-wgw-miniwdl $plugin
# next: $dest/code/bin/ugc-wgw init <project_dir> --install $dest --ref-map $ref_map [--profile cpu|parabricks|assembly]
REPORT
