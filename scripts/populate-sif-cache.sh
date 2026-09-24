#!/usr/bin/env bash
# Pull every image of one or more image manifests into a SIF cache directory, named the way
# miniwdl's singularity backend looks them up: docker://<image> with '/' and ':' replaced by '_',
# plus ".sif". Same convention as vendor/hifi-human-wgs-wdl/scripts/populate_miniwdl_singularity_cache.sh
# (verified against miniwdl 1.15.0 WDL/runtime/backend/singularity.py), extended with plain-HTTP
# registries (plain HTTP) and a machine-readable listing. Dev machine only; the HPC never pulls.
#
# Usage: scripts/populate-sif-cache.sh --dir <sif-cache> [--no-https-host <host[:port]>]... [--list <out.tsv>]
#                                      <image_manifest.txt>...
# Manifest lines: "<registry>/<name>@sha256:<digest>" or "<name>:<tag>"; '#' and blank lines ignored.
# Existing SIFs are kept (incremental). --list writes "<image>\t<sif filename>" per image.
set -euo pipefail

sif_dir=""
list_out=""
no_https_hosts=()
manifests=()
while [ $# -gt 0 ]; do
  case "$1" in
    --dir) sif_dir="$2"; shift 2 ;;
    --no-https-host) no_https_hosts+=("$2"); shift 2 ;;
    --list) list_out="$2"; shift 2 ;;
    -h|--help) sed -n '2,13p' "$0"; exit 0 ;;
    -*) echo "unknown option $1" >&2; exit 2 ;;
    *) manifests+=("$1"); shift ;;
  esac
done
[ -n "$sif_dir" ] || { echo "--dir is required" >&2; exit 2; }
[ ${#manifests[@]} -gt 0 ] || { echo "at least one manifest file is required" >&2; exit 2; }
command -v apptainer >/dev/null || { echo "apptainer not found" >&2; exit 1; }
mkdir -p "$sif_dir"
[ -z "$list_out" ] || : > "$list_out"

sif_name() {  # image reference -> SIF file name (miniwdl convention)
  local uri="docker://$1"
  uri="${uri//\//_}"
  uri="${uri//:/_}"
  printf '%s.sif' "$uri"
}

failed=0
while read -r image; do
  [ -n "$image" ] || continue
  name=$(sif_name "$image")
  target="$sif_dir/$name"
  opts=()
  for host in "${no_https_hosts[@]+"${no_https_hosts[@]}"}"; do
    case "$image" in "$host"/*) opts+=(--no-https) ;; esac
  done
  if [ -f "$target" ]; then
    echo "[sif] exists: $name" >&2
  else
    echo "[sif] pulling $image" >&2
    if apptainer pull "${opts[@]+"${opts[@]}"}" "$target.partial" "docker://$image" >&2; then
      mv "$target.partial" "$target"
    else
      echo "[sif] FAILED: $image" >&2
      rm -f "$target.partial"
      failed=$((failed + 1))
      continue
    fi
  fi
  [ -z "$list_out" ] || printf '%s\t%s\n' "$image" "$name" >> "$list_out"
done < <(cat "${manifests[@]}" | sed 's/#.*//' | awk 'NF { print $1 }' | sort -u)

if [ "$failed" -gt 0 ]; then
  echo "[sif] $failed image(s) failed" >&2
  exit 1
fi
echo "[sif] cache complete: $(find "$sif_dir" -maxdepth 1 -name '*.sif' | wc -l) SIF(s), $(du -sh "$sif_dir" | cut -f1)" >&2
