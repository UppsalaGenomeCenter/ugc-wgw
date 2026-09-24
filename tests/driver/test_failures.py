import tempfile
import unittest
from pathlib import Path

from ugc_wgw import failures
from ugc_wgw.engine import RunResult

from .helpers import make_project


class ClassifyTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.attempt = self.tmp / "attempt-1"
        self.task = self.attempt / "call-x"
        self.task.mkdir(parents=True)

    def result(self, cls, *, exit_status=None, exit_code=None, message="", task_dir=True, node="call-x"):
        return RunResult("failed", exit_code, cls, message, task_dir=str(self.task) if task_dir else None,
                         node=node, exit_status=exit_status)

    def test_rule_table(self):
        cases = [
            ("version_mismatch", {}, "version"),
            ("driver_lost", {}, "transient"), ("Interrupted", {}, "transient"), ("Terminated", {}, "transient"),
            ("launch_error", {}, "transient"), ("killed", {}, "transient"), ("NoResult", {}, "transient"),
            ("InputError", {}, "input"), ("DownloadFailed", {}, "input"),
            ("CommandFailed", {"exit_status": 137}, "resource"), ("CommandFailed", {"exit_code": 253}, "resource"),
            ("CommandFailed", {"exit_status": 1}, "tool"), ("OutputError", {}, "tool"), ("SomethingElse", {}, "tool"),
            (None, {}, "unknown"), ("Unknown", {}, "unknown"),
        ]
        for cls, kw, kind in cases:
            c = failures.classify(self.result(cls, **kw), self.attempt)
            self.assertEqual(c.kind, kind, (cls, kw, c))

    def test_resource_from_slurm_log_text(self):
        (self.task / "slurm_singularity.log.txt").write_text(
            "12345\nslurmstepd: error: *** JOB 12345 ON n1 CANCELLED AT 2026-09-11T10:00:00 DUE TO TIME LIMIT ***\n")
        c = failures.classify(self.result("CommandFailed", exit_status=143), self.attempt)
        self.assertEqual(c.kind, "resource")
        self.assertTrue(str(c.evidence["source"]).endswith("slurm_singularity.log.txt"))
        self.assertIn("DUE TO TIME LIMIT", str(c.evidence["line"]))
        self.assertEqual(c.evidence["task_dir"], str(self.task))
        self.assertIn("exit_status=143", c.message)
        self.assertIn("node=call-x", c.message)

    def test_sbatch_refusal_is_a_resource_failure(self):
        (self.task / "slurm_singularity.log.txt").write_text(
            "sbatch: error: CPU count per node can not be satisfied\n"
            "sbatch: error: Batch job submission failed: Requested node configuration is not available\n")
        c = failures.classify(self.result("CommandFailed", exit_status=1), self.attempt)
        self.assertEqual(c.kind, "resource")
        self.assertIn("CPU count per node", str(c.evidence["line"]))

    def test_preemption_text_is_transient(self):
        (self.task / "slurm_singularity.log.txt").write_text("777\n*** JOB 777 CANCELLED DUE TO PREEMPTION ***\n")
        c = failures.classify(self.result("CommandFailed", exit_status=1), self.attempt)
        self.assertEqual(c.kind, "transient")

    def test_input_from_stderr_tail(self):
        (self.task / "stderr.txt").write_text("samtools: /data/x.bam: No such file or directory\n")
        c = failures.classify(self.result("CommandFailed", exit_status=1), self.attempt)
        self.assertEqual(c.kind, "input")
        self.assertTrue(str(c.evidence["source"]).endswith("stderr.txt"))

    def test_input_text_only_counts_for_command_failures(self):
        (self.task / "stderr.txt").write_text("No such file or directory\n")
        c = failures.classify(self.result("SomethingElse", exit_status=1), self.attempt)
        self.assertEqual(c.kind, "tool")

    def test_miniwdl_stderr_when_no_task_dir(self):
        (self.attempt / "miniwdl.stderr").write_text("Out of memory\n")
        c = failures.classify(self.result("CommandFailed", task_dir=False, node=None), self.attempt)
        self.assertEqual(c.kind, "resource")
        self.assertTrue(str(c.evidence["source"]).endswith("miniwdl.stderr"))

    def test_message_is_short_and_uses_engine_message(self):
        long = "x" * 500
        c = failures.classify(self.result("InputError", message=long), self.attempt)
        self.assertLessEqual(len(c.message), failures.MESSAGE_CHARS)
        self.assertTrue(c.message.startswith("InputError node=call-x: xxx"))

    def test_backoff_delay(self):
        cfg = make_project(self.tmp)
        self.assertEqual([failures.backoff_delay(cfg, n) for n in (1, 2, 3, 4, 5)], [300.0, 600.0, 1200.0, 2400.0, 3600.0])
        cfg.backoff_seconds = 0
        self.assertEqual(failures.backoff_delay(cfg, 3), 0.0)


if __name__ == "__main__":
    unittest.main()
