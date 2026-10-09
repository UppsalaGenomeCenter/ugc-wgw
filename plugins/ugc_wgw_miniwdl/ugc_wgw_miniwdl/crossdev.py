"""Hardlink, else symlink: miniwdl's `output_hardlinks` across file systems (docs/guide/04-install.md).

With `[file_io] output_hardlinks = true` miniwdl builds every run's `out/` with `os.link`, and a target on
another file system (the cached output of an earlier run whose results live elsewhere, or an input file
a workflow passes through as an output) fails with `OSError: [Errno 18] Invalid cross-device link`.
miniwdl's own code marks the spot with "TODO: what if target is an input from a different filesystem?".

`install()` wraps `WDL._util.symlink_force`, the helper behind both hardlinks and symlinks: a hardlink
that fails with EXDEV becomes a symlink, with a WARNING on the `wdl.ugc-wgw` logger (miniwdl's log
handlers sit on the root logger, so it reaches the run's stderr and workflow log). Work directories are
always on the same file system as their run's `out/`, so `delete_work = success` keeps working: only
links to files that outlive the run turn into symlinks, which is what miniwdl itself makes when
hardlinks are off.

Imported, and therefore installed, in two ways: the bundle's venv carries `ugc_wgw_miniwdl.pth` (written
by scripts/install-bundle.sh) that imports this module at interpreter start, which covers a run whose
whole workflow is a cache hit (linked before any plugin loads); and the task plugin module imports it,
which covers any venv where only the plugin is installed.
"""
from __future__ import annotations

import errno
import logging
from typing import Any, Callable, Optional

LOGGER_NAME = "wdl.ugc-wgw"
MESSAGE = "ugc-wgw cross-device output: symlinked instead of hardlinked"
MARK = "ugc_wgw_crossdev"


def _structured(message: str, **kwargs: Any) -> object:
    try:
        from WDL._util import StructuredLogMessage  # type: ignore
    except ImportError:  # pragma: no cover - only outside miniwdl
        return message + " :: " + ", ".join(f"{k}: {v}" for k, v in kwargs.items())
    return StructuredLogMessage(message, **kwargs)


def wrap(original: Callable[..., None], log: logging.Logger) -> Callable[..., None]:
    """`symlink_force(src, dst, hard=False)` that falls back from a cross-device hardlink to a symlink."""

    def symlink_force(src: str, dst: str, hard: bool = False) -> None:
        try:
            original(src, dst, hard=hard)
            return
        except OSError as exc:
            if not hard or exc.errno != errno.EXDEV:
                raise
        log.warning(_structured(MESSAGE, source=src, link=dst))
        original(src, dst, hard=False)

    symlink_force.__wrapped__ = original  # type: ignore[attr-defined]
    setattr(symlink_force, MARK, True)
    return symlink_force


def install(util: Optional[Any] = None) -> bool:
    """Patch `WDL._util.symlink_force` (or `util.symlink_force` when a module is given). Returns True when
    the patch was applied now; False when it was already in place or miniwdl is not importable."""
    if util is None:
        try:
            import WDL._util as util  # type: ignore
        except ImportError:
            return False
    current = getattr(util, "symlink_force", None)
    if current is None or getattr(current, MARK, False):
        return False
    util.symlink_force = wrap(current, logging.getLogger(LOGGER_NAME))
    return True


def installed(util: Optional[Any] = None) -> bool:
    if util is None:
        try:
            import WDL._util as util  # type: ignore
        except ImportError:
            return False
    return bool(getattr(getattr(util, "symlink_force", None), MARK, False))


install()
