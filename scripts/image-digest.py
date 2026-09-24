#!/usr/bin/env python3
"""Resolve container image references to registry manifest digests (stdlib only, dev machine only).

Digest-form references are echoed. Tag-form references are resolved through the registry HTTP
API (anonymous token flow for Docker Hub, quay.io, nvcr.io). Output, one line per input:
    <reference>\t<digest or "unresolved">\t<media type or error>

Usage:
  scripts/image-digest.py google/deepvariant:1.10.0 quay.io/pacbio/pbmm2@sha256:...
  scripts/image-digest.py --manifest image_manifest.txt [--manifest ...]
See docs/DESIGN.md §11.3: tag-form upstream images are pinned by SIF checksum plus this digest.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.error
import urllib.parse
import urllib.request

ACCEPT = ", ".join([
    "application/vnd.oci.image.index.v1+json",
    "application/vnd.docker.distribution.manifest.list.v2+json",
    "application/vnd.oci.image.manifest.v1+json",
    "application/vnd.docker.distribution.manifest.v2+json",
])


def split_ref(ref: str) -> tuple[str, str, str, bool]:
    """-> (host, repository, reference, is_digest)"""
    first, _, rest = ref.partition("/")
    if rest and ("." in first or ":" in first or first == "localhost"):
        host = first
    else:
        host, rest = "registry-1.docker.io", ref
    if "@" in rest:
        repo, tag = rest.split("@", 1)
        is_digest = True
    elif ":" in rest.rsplit("/", 1)[-1]:
        repo, tag = rest.rsplit(":", 1)
        is_digest = False
    else:
        repo, tag = rest, "latest"
        is_digest = False
    if host == "docker.io":
        host = "registry-1.docker.io"
    if host == "registry-1.docker.io" and "/" not in repo:
        repo = "library/" + repo
    return host, repo, tag, is_digest


def _get(url: str, headers: dict[str, str]) -> tuple[dict[str, str], bytes]:
    req = urllib.request.Request(url, headers=headers, method="GET")
    with urllib.request.urlopen(req, timeout=60) as r:
        return {k.lower(): v for k, v in r.headers.items()}, r.read()


def _token(challenge: str, repo: str) -> str:
    fields = dict(re.findall(r'(\w+)="([^"]*)"', challenge))
    query = {"service": fields.get("service", ""), "scope": fields.get("scope") or f"repository:{repo}:pull"}
    _, body = _get(fields["realm"] + "?" + urllib.parse.urlencode(query), {})
    doc = json.loads(body)
    return str(doc.get("token") or doc.get("access_token") or "")


def resolve(ref: str, plain_http: bool = False) -> tuple[str, str]:
    host, repo, tag, is_digest = split_ref(ref)
    if is_digest:
        return tag, "digest-form"
    scheme = "http" if plain_http else "https"
    url = f"{scheme}://{host}/v2/{repo}/manifests/{tag}"
    headers = {"Accept": ACCEPT}
    try:
        h, _ = _get(url, headers)
    except urllib.error.HTTPError as exc:
        if exc.code != 401:
            raise
        headers["Authorization"] = "Bearer " + _token(exc.headers.get("WWW-Authenticate", ""), repo)
        h, _ = _get(url, headers)
    return h.get("docker-content-digest", "unresolved"), h.get("content-type", "")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("refs", nargs="*")
    ap.add_argument("--manifest", action="append", default=[], help="image manifest file (one reference per line)")
    ap.add_argument("--plain-http-host", action="append", default=[], help="registry host to query over http")
    a = ap.parse_args()
    refs = list(a.refs)
    for path in a.manifest:
        with open(path) as fh:
            for line in fh:
                s = line.split("#", 1)[0].strip()
                if s:
                    refs.append(s)
    status = 0
    for ref in refs:
        plain = any(ref.startswith(h + "/") for h in a.plain_http_host)
        try:
            digest, media = resolve(ref, plain)
        except Exception as exc:  # noqa: BLE001 - report and continue
            digest, media = "unresolved", f"error: {exc}"
            status = 1
        print(f"{ref}\t{digest}\t{media}")
    return status


if __name__ == "__main__":
    sys.exit(main())
