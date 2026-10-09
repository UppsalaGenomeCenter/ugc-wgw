"""The hardlink-else-symlink wrap around miniwdl's symlink_force, driven with a stand-in module (no miniwdl needed)."""
from __future__ import annotations

import errno
import logging
import os
import tempfile
import unittest
import uuid
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from ugc_wgw_miniwdl import crossdev


def symlink_force(src: str, dst: str, hard: bool = False) -> None:
    """miniwdl 1.15.0's WDL._util.symlink_force, verbatim."""
    assert not dst.endswith("/")
    tn = dst + ".tmp." + str(uuid.uuid1())
    if hard:
        os.link(src, tn)
    else:
        os.symlink(src, tn)
    os.rename(tn, dst)


class CrossdevTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.src = self.tmp / "src.txt"
        self.src.write_text("data\n")
        self.util = SimpleNamespace(symlink_force=symlink_force)

    def test_install_is_idempotent_and_marks_the_wrap(self):
        self.assertFalse(crossdev.installed(self.util))
        self.assertTrue(crossdev.install(self.util))
        self.assertTrue(crossdev.installed(self.util))
        self.assertFalse(crossdev.install(self.util))           # already wrapped
        self.assertIs(self.util.symlink_force.__wrapped__, symlink_force)
        self.assertFalse(crossdev.install(SimpleNamespace()))   # nothing to wrap

    def test_hardlink_and_symlink_still_work(self):
        crossdev.install(self.util)
        hard, soft = self.tmp / "hard.txt", self.tmp / "soft.txt"
        self.util.symlink_force(str(self.src), str(hard), hard=True)
        self.util.symlink_force(str(self.src), str(soft))
        self.assertFalse(hard.is_symlink())
        self.assertEqual(os.stat(hard).st_ino, os.stat(self.src).st_ino)
        self.assertTrue(soft.is_symlink())
        self.assertEqual(os.readlink(soft), str(self.src))

    def test_cross_device_hardlink_becomes_a_symlink_with_a_warning(self):
        crossdev.install(self.util)
        dst = self.tmp / "out.txt"
        real_link = os.link

        def exdev(src, dst, *a, **k):
            raise OSError(errno.EXDEV, "Invalid cross-device link", src, None, dst)

        with mock.patch("os.link", side_effect=exdev), self.assertLogs(crossdev.LOGGER_NAME, level="WARNING") as logs:
            self.util.symlink_force(str(self.src), str(dst), hard=True)
        self.assertTrue(dst.is_symlink())
        self.assertEqual(os.readlink(dst), str(self.src))
        self.assertEqual(dst.read_text(), "data\n")
        self.assertEqual(len(logs.records), 1)
        self.assertIn(crossdev.MESSAGE, str(logs.records[0].getMessage()))
        self.assertIn(str(dst), str(logs.records[0].getMessage()))
        self.assertEqual(logs.records[0].name, "wdl.ugc-wgw")
        self.assertIs(os.link, real_link)
        # a stray temp name is never left behind
        self.assertEqual(sorted(p.name for p in self.tmp.iterdir()), ["out.txt", "src.txt"])

    def test_other_errors_propagate(self):
        crossdev.install(self.util)
        with mock.patch("os.link", side_effect=OSError(errno.EACCES, "Permission denied")):
            with self.assertRaises(OSError) as ctx:
                self.util.symlink_force(str(self.src), str(self.tmp / "x.txt"), hard=True)
        self.assertEqual(ctx.exception.errno, errno.EACCES)
        with self.assertRaises(OSError):   # a symlink that fails is not retried
            with mock.patch("os.symlink", side_effect=OSError(errno.EXDEV, "x")):
                self.util.symlink_force(str(self.src), str(self.tmp / "y.txt"))

    def test_real_miniwdl_is_wrapped_when_importable(self):
        try:
            import WDL._util as util  # type: ignore
        except ImportError:
            self.skipTest("miniwdl not importable")
        self.assertTrue(crossdev.installed(util))   # importing crossdev (or resources) installed it
        self.assertTrue(getattr(util.link_force, "__name__", "") == "link_force")

    def test_logger_propagates_to_root_handlers(self):
        """miniwdl attaches its handlers to the root logger; the warning must reach them."""
        crossdev.install(self.util)
        root = logging.getLogger()
        with mock.patch("os.link", side_effect=OSError(errno.EXDEV, "x")), self.assertLogs(root, level="WARNING") as logs:
            self.util.symlink_force(str(self.src), str(self.tmp / "z.txt"), hard=True)
        self.assertEqual([r.name for r in logs.records], ["wdl.ugc-wgw"])


if __name__ == "__main__":
    unittest.main()
