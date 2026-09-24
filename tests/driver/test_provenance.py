"""Tests for the code provenance the driver records: version and commit of a checkout or an installed bundle."""
from __future__ import annotations

import json
import re
import tempfile
import unittest
from pathlib import Path

from ugc_wgw import config, manifest

from .helpers import REPO


class GitCommitTest(unittest.TestCase):
    def test_checkout_uses_git(self):
        commit = config.git_commit(REPO)
        self.assertRegex(commit, r"^[0-9a-f]{40}$")

    def test_installed_bundle_uses_its_manifest(self):
        tmp = Path(tempfile.mkdtemp())
        version_dir = tmp / "versions" / "0.1.0"
        code = version_dir / "code"
        code.mkdir(parents=True)
        (code / "VERSION").write_text("0.1.0\n")
        self.assertEqual(config.git_commit(code), "unknown")
        (version_dir / "manifest.json").write_text(json.dumps(
            {"schema": 1, "ugc_pacbio_wgw": {"version": "0.1.0", "git_ref": "HEAD", "git_commit": "4db1257b" * 5}}))
        self.assertEqual(config.git_commit(code), "4db1257b" * 5)
        info = manifest.load_code_info(code)
        self.assertEqual((info.version, info.git_commit), ("0.1.0", "4db1257b" * 5))
        (version_dir / "manifest.json").write_text("not json")
        self.assertEqual(config.git_commit(code), "unknown")

    def test_code_dir_inside_a_repo_does_not_borrow_its_commit(self):
        tmp = Path(tempfile.mkdtemp(dir=REPO / "tests"))  # under the repo, but no .git of its own
        try:
            (tmp / "VERSION").write_text("9.9.9\n")
            self.assertEqual(config.git_commit(tmp), "unknown")
        finally:
            for f in tmp.iterdir():
                f.unlink()
            tmp.rmdir()


if __name__ == "__main__":
    unittest.main()
