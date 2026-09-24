#!/usr/bin/env python3
"""Vendor an upstream repository at a tag into vendor/<name>/ and record it in upstream.lock.

The copy is verbatim (minus .git directories) and includes submodules resolved
at the commits the upstream tag pins. This is the ONLY sanctioned way to write
under vendor/. See docs/DESIGN.md §3, §11.1, §12.

Usage:
  scripts/vendor-upstream.py --name hifi-human-wgs-wdl --tag v4.0.0
  scripts/vendor-upstream.py --list
  scripts/vendor-upstream.py --name hifi-human-wgs-wdl --tag v4.1.0 --dry-run
  scripts/vendor-upstream.py --name <name> --remove        # delete vendor/<name> and its lock entry

A repository that is not in DEFAULT_URLS needs --url the first time. Removal is the
only sanctioned way to delete a vendored tree (it keeps the lockfile in step).

Requires: git, python3 (stdlib only). Network access to the upstream URL.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

DEFAULT_URLS = {
    "hifi-human-wgs-wdl": "https://github.com/PacificBiosciences/HiFi-human-WGS-WDL.git",
}

REPO_ROOT = Path(__file__).resolve().parent.parent
LOCK_DEFAULT = REPO_ROOT / "upstream.lock"
VENDOR_DEFAULT = REPO_ROOT / "vendor"


def run(cmd: list[str], cwd: Path | None = None) -> str:
    res = subprocess.run(cmd, cwd=cwd, check=True, capture_output=True, text=True)
    return res.stdout.strip()


def load_lock(path: Path) -> dict:
    if path.exists():
        with path.open() as fh:
            return json.load(fh)
    return {"schema": 1, "sources": {}}


def save_lock(path: Path, lock: dict) -> None:
    with path.open("w") as fh:
        json.dump(lock, fh, indent=2, sort_keys=True)
        fh.write("\n")


def submodule_commits(clone: Path) -> dict[str, str]:
    """Return {submodule_path: commit} for all (recursively) checked-out submodules."""
    out = run(["git", "submodule", "status", "--recursive"], cwd=clone)
    result: dict[str, str] = {}
    for line in out.splitlines():
        line = line.strip()
        if not line:
            continue
        # format: "[ +-U]<sha> <path> (<describe>)"; strip status char if present
        if line[0] in "+-U ":
            line = line[1:].strip()
        sha, path = line.split()[:2]
        result[path] = sha
    return result


def ignore_git(_dir: str, names: list[str]) -> set[str]:
    return {n for n in names if n == ".git"}


def copy_tree(src: Path, dst: Path) -> None:
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst, ignore=ignore_git, symlinks=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", help="vendor directory name (e.g. hifi-human-wgs-wdl)")
    ap.add_argument("--tag", help="upstream git tag to vendor")
    ap.add_argument("--url", help="override upstream URL (default: lockfile, then built-in map)")
    ap.add_argument("--lock", type=Path, default=LOCK_DEFAULT)
    ap.add_argument("--vendor-dir", type=Path, default=VENDOR_DEFAULT)
    ap.add_argument("--list", action="store_true", help="print the lockfile and exit")
    ap.add_argument("--dry-run", action="store_true", help="clone and report, but do not write")
    ap.add_argument("--remove", action="store_true", help="delete vendor/<name> and its upstream.lock entry")
    args = ap.parse_args()

    lock = load_lock(args.lock)

    if args.list:
        print(json.dumps(lock, indent=2, sort_keys=True))
        return 0

    if args.remove:
        if not args.name:
            ap.error("--remove needs --name")
        dest = args.vendor_dir / args.name
        entry = lock["sources"].pop(args.name, None)
        if entry is None and not dest.exists():
            ap.error(f"{args.name} is neither in {args.lock} nor under {args.vendor_dir}")
        if args.dry_run:
            print(f"[vendor] dry run: would remove {dest} and the lock entry ({(entry or {}).get('tag')})", file=sys.stderr)
            return 0
        if dest.exists():
            shutil.rmtree(dest)
        save_lock(args.lock, lock)
        print(f"[vendor] removed {dest} and its entry ({(entry or {}).get('tag')}, {(entry or {}).get('commit')})\n"
              f"Suggested commit: git add -A {dest.relative_to(REPO_ROOT) if dest.is_relative_to(REPO_ROOT) else dest} "
              f"{args.lock.relative_to(REPO_ROOT) if args.lock.is_relative_to(REPO_ROOT) else args.lock}", file=sys.stderr)
        return 0

    if not args.name or not args.tag:
        ap.error("--name and --tag are required (or use --list or --remove)")

    entry = lock["sources"].get(args.name, {})
    url = args.url or entry.get("url") or DEFAULT_URLS.get(args.name)
    if not url:
        ap.error(f"no URL known for {args.name}; pass --url")

    dest = args.vendor_dir / args.name
    previous = entry.get("tag")

    with tempfile.TemporaryDirectory(prefix="ugc-vendor-") as tmp:
        clone = Path(tmp) / args.name
        print(f"[vendor] cloning {url} @ {args.tag} ...", file=sys.stderr)
        run([
            "git", "clone", "--quiet", "--recursive", "--shallow-submodules",
            "--depth", "1", "--branch", args.tag, url, str(clone),
        ])
        commit = run(["git", "rev-parse", "HEAD"], cwd=clone)
        subs = submodule_commits(clone)

        print(f"[vendor] {args.name}: {args.tag} -> {commit}", file=sys.stderr)
        for path, sha in subs.items():
            print(f"[vendor]   submodule {path} -> {sha}", file=sys.stderr)
        if previous:
            print(f"[vendor] previously vendored: {previous} ({entry.get('commit')})", file=sys.stderr)

        if args.dry_run:
            print("[vendor] dry run; nothing written", file=sys.stderr)
            return 0

        print(f"[vendor] copying into {dest} ...", file=sys.stderr)
        copy_tree(clone, dest)

    new_entry = {
        "url": url,
        "tag": args.tag,
        "commit": commit,
        "submodules": subs,
        "vendored_at": dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "notes": entry.get("notes", ""),
    }
    lock["sources"][args.name] = new_entry
    save_lock(args.lock, lock)

    # Show what changed, if we are inside a git work tree.
    try:
        stat = run(["git", "diff", "--stat", "--", str(dest), str(args.lock)], cwd=REPO_ROOT)
        if stat:
            print("\n" + stat, file=sys.stderr)
    except subprocess.CalledProcessError:
        pass

    def rel(p: Path) -> str:
        try:
            return str(p.resolve().relative_to(REPO_ROOT))
        except ValueError:
            return str(p)

    prev = f" (was {previous})" if previous else ""
    print(
        "\nSuggested commit:\n"
        f"  git add {rel(dest)} {rel(args.lock)}\n"
        f"  git commit -m 'vendor: {args.name} {args.tag}{prev}'\n"
        "\nNext: reconcile every row of DERIVED_FILES.md against the new tree, then tests/check.sh.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
