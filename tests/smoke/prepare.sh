#!/usr/bin/env bash
# Prepare the smoke dataset once, on a machine with internet: the reference tree copied out
# of PacBio's reference data container (the SIF the bundle ships) and chr20 slices of the
# GIAB Revio trio (aligned public BAMs, fetched by region; upstream's pbmm2 step realigns
# them). See tests/smoke/PREPARE.md.
#
# Usage: prepare.sh [--references-only | --data-only] [--threads N] [--force]
# Environment: UGC_WGW_SMOKE_DIR (default $HOME/ugc-smoke); UGC_WGW_SIF_CACHE (default
# bundle/out/sif-cache, where the data container SIF is or gets pulled); needs samtools >= 1.10
# with htslib https support, curl, md5sum, python3, tar, apptainer; about 25 GB of disk.
# Writes:
#   $UGC_WGW_SMOKE_DIR/references/hifi-wdl-resources-v4.0.0-GRCh38_GIABv3/   the container's data tree + manifest.json
#   $UGC_WGW_SMOKE_DIR/references/ugc-wgw-extras-0.2.0/                          scatter_regions.GRCh38_GIABv3.tsv
#   $UGC_WGW_SMOKE_DIR/data/<ID>.chr20.hifi_reads.bam           HG002, HG003, HG004, HG002ds
#   $UGC_WGW_SMOKE_DIR/data/samples.tsv, cohort.txt, data.lock  driver inputs and provenance
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
smoke=${UGC_WGW_SMOKE_DIR:-$HOME/ugc-smoke}
sif_cache=${UGC_WGW_SIF_CACHE:-$repo/bundle/out/sif-cache}
threads=4; do_refs=1; do_data=1; force=0
while [ $# -gt 0 ]; do
  case "$1" in
    --references-only) do_data=0; shift ;;
    --data-only) do_refs=0; shift ;;
    --threads) threads="$2"; shift 2 ;;
    --force) force=1; shift ;;
    -h|--help) sed -n '2,14p' "$0"; exit 0 ;;
    *) echo "unknown option $1" >&2; exit 2 ;;
  esac
done
for tool in samtools curl md5sum python3 tar apptainer; do
  command -v "$tool" >/dev/null || { echo "[prepare] $tool not found" >&2; exit 1; }
done
log() { echo "[prepare] $(date -u +%H:%M:%S) $*" >&2; }
mkdir -p "$smoke/references" "$smoke/data" "$smoke/downloads"

# ---- 1. references from the data container (background while the BAM slices download) ----
prepare_references() {
  local image subdir data_path manifest_path extras extras_file
  image=""; subdir=""; data_path=""; manifest_path=""; extras=""; extras_file=""
  eval "$(python3 - "$repo/references.lock" <<'PY'
import json, shlex, sys
b = json.load(open(sys.argv[1]))["bundles"]
r, x = b["hifi-wdl-resources"], b["ugc-wgw-extras"]
for k, v in {"image": r["source"], "subdir": r["install_subdir"], "data_path": r["container_paths"]["data"],
             "manifest_path": r["container_paths"]["manifest"], "extras": x["install_subdir"], "extras_file": x["files"][0]["name"]}.items():
    print(f"{k}={shlex.quote(v)}")
PY
)"
  local tree=$smoke/references/$subdir marker=$smoke/references/.$subdir.ok
  if [ -f "$marker" ] && [ -f "$tree/manifest.json" ] && [ "$force" -eq 0 ]; then
    log "references: $subdir present (marker $marker)"
  else
    # the same SIF the bundle ships (miniwdl's cache name); pulled into the cache when absent
    local sif=$sif_cache/docker___${image//[\/:]/_}.sif
    if [ ! -f "$sif" ]; then
      log "references: pulling the reference data container into $sif_cache"
      printf '%s\n' "$image" > "$smoke/downloads/data-manifest.txt"
      "$repo/scripts/populate-sif-cache.sh" --dir "$sif_cache" "$smoke/downloads/data-manifest.txt"
    fi
    rm -rf "$tree"; mkdir -p "$tree"
    log "references: copying $data_path and $manifest_path out of $(basename "$sif")"
    apptainer exec --bind "$tree:/mnt/ugc-references" "$sif" \
      sh -c "cp -r '$data_path' /mnt/ugc-references/ && cp '$manifest_path' /mnt/ugc-references/"
    log "references: verifying md5s against references.lock"
    python3 - "$repo/references.lock" "$tree" <<'PY'
import hashlib, json, os, sys
r = json.load(open(sys.argv[1]))["bundles"]["hifi-wdl-resources"]
bad = []
for f in r["files"]:
    p = os.path.join(sys.argv[2], f["name"])
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    if os.path.getsize(p) != f["bytes"] or h.hexdigest() != f["md5"]:
        bad.append(f["name"])
if bad:
    print("[prepare] reference files differ from references.lock: " + " ".join(bad), file=sys.stderr)
    sys.exit(1)
print(f"[prepare] {len(r['files'])} reference files verified", file=sys.stderr)
PY
    touch "$marker"
  fi
  mkdir -p "$smoke/references/$extras"
  cp "$repo/references/$extras_file" "$smoke/references/$extras/"
  log "references: ready"
}

# ---- 2. chr20 slices of the GIAB Revio trio ------------------------------------
base=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio
# id  giab_dir  file  sex  father  mother
samples=(
  "HG002 HG002_NA24385_son HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam MALE HG003 HG004"
  "HG003 HG003_NA24149_father HG003_PacBio-HiFi-Revio_20231031_46x_GRCh38-GIABv3.bam MALE - -"
  "HG004 HG004_NA24143_mother HG004_PacBio-HiFi-Revio_20231031_36x_GRCh38-GIABv3.bam FEMALE - -"
)
region=chr20

recorded_md5() {  # recorded_md5 <file basename> -> md5 from data.lock or empty
  [ -f "$smoke/data/data.lock" ] || return 0
  python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); print(d.get("files",{}).get(sys.argv[2],{}).get("md5",""))' \
      "$smoke/data/data.lock" "$1"
}

fetch_slice() {  # fetch_slice <id> <url>
  local id=$1 url=$2 out=$smoke/data/$1.$region.hifi_reads.bam want
  want=$(recorded_md5 "$(basename "$out")")
  if [ -f "$out" ] && [ "$force" -eq 0 ]; then
    if [ -n "$want" ] && [ "$(md5sum "$out" | cut -d' ' -f1)" = "$want" ]; then
      log "data: $id slice present and verified"; return 0
    fi
    [ -z "$want" ] && { log "data: $id slice present (no record to verify against)"; return 0; }
    log "data: $id slice md5 differs from data.lock; refetching"
  fi
  local try
  for try in 1 2 3; do
    log "data: fetching $region of $id (try $try): $url"
    # run from downloads/ so that htslib's copy of the remote .bai lands there, not in the caller's cwd
    if (cd "$smoke/downloads" && samtools view -b -@ "$threads" -o "$out.partial" "$url" "$region"); then
      mv "$out.partial" "$out"; return 0
    fi
    rm -f "$out.partial"; sleep 30
  done
  echo "[prepare] failed to fetch $id" >&2; return 1
}

prepare_data() {
  local line id dir file url
  for line in "${samples[@]}"; do
    read -r id dir file _ _ _ <<<"$line"
    url=$base/$dir/PacBio_HiFi-Revio_20231031/$file
    fetch_slice "$id" "$url"
  done
  local ds=$smoke/data/HG002ds.$region.hifi_reads.bam
  if [ ! -f "$ds" ] || [ "$force" -eq 1 ]; then
    log "data: downsampling HG002 to 50 % as HG002ds"
    samtools view -b -@ "$threads" -s 42.5 -o "$ds.partial" "$smoke/data/HG002.$region.hifi_reads.bam"
    # `samtools view -s` writes a free-text "@CO Sub-sampled fraction=..." header line. HiPhase's header parser
    # (rust-htslib 0.39.5) panics on a header line without TAG:value fields, and pbmm2 in v4.0.0 keeps @CO lines
    # of the input (the v3 path's `samtools reset` dropped them). Strip it; instrument BAMs carry no @CO lines.
    samtools reheader -c 'grep -v "^@CO"' "$ds.partial" > "$ds.partial.reheader"
    mv "$ds.partial.reheader" "$ds"
    rm -f "$ds.partial"
  fi
  # driver inputs
  {
    printf 'sample_id\tsex\thifi_reads\tfail_reads\tfather_id\tmother_id\n'
    for line in "${samples[@]}"; do
      read -r id _ _ sex father mother <<<"$line"
      [ "$father" = - ] && father=""; [ "$mother" = - ] && mother=""
      printf '%s\t%s\t%s\t\t%s\t%s\n' "$id" "$sex" "$smoke/data/$id.$region.hifi_reads.bam" "$father" "$mother"
    done
    printf 'HG002ds\t\t%s\t\t\t\n' "$ds"
  } > "$smoke/data/samples.tsv"
  printf '# smoke cohort: order is frozen with the cohort\nHG002\nHG003\nHG004\nHG002ds\n' > "$smoke/data/cohort.txt"
  # provenance
  log "data: computing md5s for data.lock"
  python3 - "$smoke/data" "$base" "$region" "${samples[@]}" <<'PY'
import datetime, hashlib, json, os, subprocess, sys
data, base, region, *rows = sys.argv[1:]
def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()
files = {}
for row in rows:
    sid, d, f, *_ = row.split()
    name = f"{sid}.{region}.hifi_reads.bam"
    p = os.path.join(data, name)
    files[name] = {"source": f"{base}/{d}/PacBio_HiFi-Revio_20231031/{f}", "region": region,
                   "bytes": os.path.getsize(p), "md5": md5(p)}
name = f"HG002ds.{region}.hifi_reads.bam"
p = os.path.join(data, name)
files[name] = {"source": f"HG002.{region}.hifi_reads.bam", "region": region, "downsample": "samtools view -s 42.5; samtools reheader dropping the @CO line",
               "bytes": os.path.getsize(p), "md5": md5(p)}
samtools = subprocess.run(["samtools", "--version"], capture_output=True, text=True).stdout.splitlines()[0]
doc = {"schema": 1, "prepared_at": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
       "samtools": samtools, "files": files}
with open(os.path.join(data, "data.lock"), "w") as fh:
    json.dump(doc, fh, indent=2, sort_keys=True); fh.write("\n")
for n, f in sorted(files.items()):
    print(f"  {n}: {f['bytes']} bytes md5 {f['md5']}")
PY
  log "data: ready ($smoke/data)"
}

rc=0
if [ "$do_refs" -eq 1 ]; then
  prepare_references > "$smoke/downloads/references.log" 2>&1 &
  refs_pid=$!
fi
if [ "$do_data" -eq 1 ]; then
  prepare_data || rc=1
fi
if [ "$do_refs" -eq 1 ]; then
  wait "$refs_pid" || { rc=1; log "references FAILED; see $smoke/downloads/references.log"; }
  tail -n 3 "$smoke/downloads/references.log" >&2
fi
if [ "$rc" -eq 0 ]; then log "done: $smoke"; else log "finished with errors"; fi
exit "$rc"
