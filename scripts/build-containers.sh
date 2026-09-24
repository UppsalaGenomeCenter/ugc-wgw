#!/usr/bin/env bash
# Build one or all containers under containers/, push them to the registry and
# print the pushed digests. See containers/README.md.
#
# Usage: scripts/build-containers.sh [<tool> ...]   (default: all)
# Env:   UGC_WGW_REGISTRY        registry/namespace (default ghcr.io/uppsalagenomecenter)
#        UGC_WGW_CONTAINER_CLI   docker | podman (default: docker if its daemon answers, else podman)
#        UGC_WGW_PUSH            1 (default) to push and print digests; 0 to only build
#        UGC_WGW_TLS_VERIFY      podman only: false (default) for a plain-HTTP registry
#
# Contract:
#   - reads ARG <TOOL>_VERSION from containers/<tool>/Dockerfile to form the tag
#   - <cli> build --pull; <cli> push; digest from RepoDigests (docker) or --digestfile (podman)
#   - prints "<registry>/<tool>@sha256:<digest>" on stdout, one per tool
#   - does NOT edit image_manifest.ugc-wgw.txt (the human/Claude does, and reviews it)
#   - log in first: `docker login <registry>` or `podman login --tls-verify=false <registry>`
set -euo pipefail
cd "$(dirname "$0")/.."
REGISTRY="${UGC_WGW_REGISTRY:-ghcr.io/uppsalagenomecenter}"  # export-public: registry default
PUSH="${UGC_WGW_PUSH:-1}"
TLS_VERIFY="${UGC_WGW_TLS_VERIFY:-true}"  # export-public: tls default

CLI="${UGC_WGW_CONTAINER_CLI:-}"
if [ -z "$CLI" ]; then
  if docker info >/dev/null 2>&1; then
    CLI=docker
  elif command -v podman >/dev/null 2>&1; then
    CLI=podman
  else
    echo "neither a usable docker daemon nor podman found" >&2
    exit 1
  fi
fi

tools=("$@")
if [ ${#tools[@]} -eq 0 ]; then
  mapfile -t tools < <(find containers -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)
fi

for tool in "${tools[@]}"; do
  df="containers/$tool/Dockerfile"
  [ -f "$df" ] || { echo "no Dockerfile for $tool" >&2; exit 1; }
  version=$(sed -nE 's/^ARG [A-Z_]*VERSION=(.*)$/\1/p' "$df" | head -1)
  [ -n "$version" ] || { echo "$df: missing 'ARG <TOOL>_VERSION=<v>'" >&2; exit 1; }
  image="$REGISTRY/$tool:$version"
  echo "== building $image with $CLI" >&2
  "$CLI" build --pull -t "$image" "containers/$tool" >&2
  if [ "$PUSH" != "1" ]; then
    echo "$image (built, not pushed: UGC_WGW_PUSH=$PUSH)"
    continue
  fi
  echo "== pushing $image" >&2
  if [ "$CLI" = "podman" ]; then
    digestfile=$(mktemp)
    "$CLI" push --tls-verify="$TLS_VERIFY" --digestfile "$digestfile" "$image" >&2
    echo "$REGISTRY/$tool@$(cat "$digestfile")"
    rm -f "$digestfile"
  else
    "$CLI" push "$image" >&2
    "$CLI" inspect --format='{{index .RepoDigests 0}}' "$image"
  fi
done
