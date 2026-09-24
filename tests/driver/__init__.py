"""Unit tests for the ugc-wgw driver (stdlib unittest). Puts bin/ on sys.path so `ugc_wgw` imports."""
import pathlib
import sys

_BIN = pathlib.Path(__file__).resolve().parents[2] / "bin"
if str(_BIN) not in sys.path:
    sys.path.insert(0, str(_BIN))
