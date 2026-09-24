import io
import logging
import os
import re
import tempfile
import time
import unittest
from pathlib import Path

from ugc_wgw import log as ugclog
from ugc_wgw.util import UgcError

from .helpers import make_project, run_cli, write_ids, write_tsv

ESC = re.compile(r"\x1b\[[0-9;]*m")


def rec(level: int, msg: str) -> logging.LogRecord:
    r = logging.LogRecord("ugc_wgw", level, "", 0, msg, None, None)
    r.created = 1790082199.0
    return r


class _Tty(io.StringIO):
    def isatty(self) -> bool:
        return True


class ColorFormatterTest(unittest.TestCase):
    def setUp(self):
        self.f = ugclog.ColorFormatter()
        self.plain = logging.Formatter(ugclog.LOG_FORMAT, datefmt=ugclog.DATE_FORMAT)
        self.plain.converter = time.gmtime

    def same_when_stripped(self, r: logging.LogRecord) -> str:
        s = self.f.format(r)
        self.assertEqual(ESC.sub("", s), self.plain.format(r))
        return s

    def test_failure_line(self):
        s = self.same_when_stripped(rec(logging.INFO, "run.failed run_id=X-a1 subject=S1 stage=singleton attempt=1 exit_code=101 "
                                                      "error_class=CommandFailed kind=tool stderr=/x/stderr.txt"))
        self.assertIn(f"{ugclog.RED}run.failed{ugclog.RESET}", s)
        self.assertIn(f"subject={ugclog.BOLD}S1{ugclog.RESET}", s)
        self.assertIn(f"stage={ugclog.BOLD}singleton{ugclog.RESET}", s)
        self.assertIn(f"error_class={ugclog.RED}CommandFailed{ugclog.RESET}", s)
        self.assertIn(f"kind={ugclog.RED}tool{ugclog.RESET}", s)
        self.assertIn(f"{ugclog.DIM}run_id=X-a1{ugclog.RESET}", s)
        self.assertIn(f"{ugclog.DIM}stderr=/x/stderr.txt{ugclog.RESET}", s)
        self.assertTrue(s.startswith(f"{ugclog.DIM}2026-09-22T13:03:19Z{ugclog.RESET} {ugclog.DIM}INFO{ugclog.RESET} "))

    def test_success_and_other_events(self):
        s = self.same_when_stripped(rec(logging.INFO, "run.success run_id=X subject=S1 stage=singleton attempt=1 exit_code=0"))
        self.assertIn(f"{ugclog.GREEN}run.success{ugclog.RESET}", s)
        self.assertNotIn(ugclog.RED, s)  # exit_code is not a failure detail on a success line
        s = self.same_when_stripped(rec(logging.INFO, "run.blocked run_id=None subject=SMOKE stage=cohort_merge reason=waiting"))
        self.assertIn(f"{ugclog.DIM}run.blocked{ugclog.RESET}", s)
        s = self.same_when_stripped(rec(logging.INFO, "run.newkind run_id=X subject=S1"))
        self.assertIn(f"{ugclog.CYAN}run.newkind{ugclog.RESET}", s)  # unknown run.* events default to cyan
        s = self.same_when_stripped(rec(logging.INFO, "submit.stop run_id=None launched=4 failures=1"))
        self.assertIn(f"{ugclog.BOLD}submit.stop{ugclog.RESET}", s)
        s = self.same_when_stripped(rec(logging.INFO, "no event token here"))
        self.assertIn(" no event token here", s)

    def test_levels_and_progress_block(self):
        s = self.same_when_stripped(rec(logging.WARNING, "cohort_freq: sample S1 is XY in the sheet but FEMALE was inferred"))
        self.assertIn(f"{ugclog.BOLD}{ugclog.YELLOW}WARNING{ugclog.RESET}", s)
        self.assertIn(f"{ugclog.BOLD}{ugclog.YELLOW}cohort_freq:{ugclog.RESET}", s)
        s = self.same_when_stripped(rec(logging.ERROR, "reconcile: run X has a dead driver pid=7 reason=gone"))
        self.assertIn(f"{ugclog.BOLD}{ugclog.RED}ERROR{ugclog.RESET}", s)
        self.assertIn(f"reason={ugclog.RED}gone{ugclog.RESET}", s)
        s = self.same_when_stripped(rec(logging.INFO, "progress mode=standalone 1/6 (16%) eta≈4h02m\n  singleton     done 1/4  active 2"))
        head, tail = s.split("\n")
        self.assertIn(f"{ugclog.BOLD}{ugclog.CYAN}progress mode=standalone 1/6 (16%) eta≈4h02m{ugclog.RESET}", head)
        self.assertEqual(tail, "  singleton     done 1/4  active 2")  # continuation lines stay plain


class GatingTest(unittest.TestCase):
    def setUp(self):
        self.env = dict(os.environ)
        os.environ.pop("NO_COLOR", None)
        os.environ["TERM"] = "xterm-256color"

    def tearDown(self):
        os.environ.clear()
        os.environ.update(self.env)

    def test_modes(self):
        self.assertTrue(ugclog.color_enabled("always", io.StringIO()))
        self.assertFalse(ugclog.color_enabled("never", _Tty()))
        self.assertFalse(ugclog.color_enabled("auto", io.StringIO()))
        self.assertTrue(ugclog.color_enabled("auto", _Tty()))
        self.assertTrue(ugclog.color_enabled("AUTO", _Tty()))
        os.environ["NO_COLOR"] = "1"
        self.assertFalse(ugclog.color_enabled("auto", _Tty()))
        self.assertTrue(ugclog.color_enabled("always", _Tty()))
        del os.environ["NO_COLOR"]
        os.environ["TERM"] = "dumb"
        self.assertFalse(ugclog.color_enabled("auto", _Tty()))
        with self.assertRaises(UgcError):
            ugclog.color_enabled("sometimes", _Tty())

    def test_cli(self):
        tmp = Path(tempfile.mkdtemp())
        cfg = make_project(tmp)
        proj = str(cfg.project_dir)
        tsv = write_tsv(tmp, [{"sample_id": "S1", "sex": "MALE", "hifi_reads": "s1.bam"}])
        self.assertEqual(run_cli(["--project", proj, "samples", "add", str(tsv)])[0], 0)
        ids = str(write_ids(tmp, ["S1"]))
        code, out, err = run_cli(["--project", proj, "-v", "--color", "always", "cohort", "freeze", "C1", "--samples", ids])
        self.assertEqual(code, 0, err)
        self.assertIn("cohort.frozen", err)
        self.assertIn("\x1b[", err)
        code, out, err = run_cli(["--project", proj, "-v", "--color", "never", "cohort", "freeze", "C2", "--samples", ids])
        self.assertEqual(code, 0, err)
        self.assertIn("cohort.frozen", err)
        self.assertNotIn("\x1b[", err)
        # auto: the captured stderr is not a terminal, so plain; UGC_WGW_COLOR sets the default
        code, out, err = run_cli(["--project", proj, "-v", "cohort", "freeze", "C3", "--samples", ids])
        self.assertEqual(code, 0, err)
        self.assertNotIn("\x1b[", err)
        code, out, err = run_cli(["--project", proj, "-v", "cohort", "freeze", "C4", "--samples", ids], env={"UGC_WGW_COLOR": "always"})
        self.assertEqual(code, 0, err)
        self.assertIn("\x1b[", err)
        code, out, err = run_cli(["--project", proj, "-v", "cohort", "freeze", "C5", "--samples", ids], env={"UGC_WGW_COLOR": "rainbow"})
        self.assertEqual(code, 1)
        self.assertIn("UGC_WGW_COLOR", err)
        # the file log never carries escapes
        self.assertNotIn("\x1b[", (cfg.logs_dir / "ugc-wgw.log").read_text())


if __name__ == "__main__":
    unittest.main()
