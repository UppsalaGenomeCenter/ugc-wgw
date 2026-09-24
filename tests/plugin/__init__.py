"""Unit tests for plugins/ugc_wgw_miniwdl (stdlib unittest, no miniwdl needed). Puts the package on sys.path."""
import pathlib
import sys

_PKG = pathlib.Path(__file__).resolve().parents[2] / "plugins" / "ugc_wgw_miniwdl"
if str(_PKG) not in sys.path:
    sys.path.insert(0, str(_PKG))
