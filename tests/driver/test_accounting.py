"""accounting.json: sacct value parsing, folding of parent and step rows, the join with task directories and the
workflow log, and the capture at run finalization."""
import json
import os
import tempfile
import unittest
from pathlib import Path

from ugc_wgw import accounting, config, layout, slurm
from ugc_wgw.db import DB, RunRecord
from ugc_wgw.util import utc_now

from .helpers import VERSION, make_project, run_cli, sacct_env, write_tsv

# pipe-separated, the fields of slurm.SACCT_FIELDS; a plain job with .batch/.extern, a GPU job with a day-long TotalCPU,
# a parent-only job, a cancelled one, and one still running
SACCT_TEXT = """\
1001|call-a|COMPLETED|0:0|2026-09-25T10:00:00|2026-09-25T11:02:10|2026-09-25T11:03:21|71|24|24|1704|20:03.400|131072M||billing=24,cpu=24,mem=131072M,node=1|core|m12|2880||
1001.batch|batch|COMPLETED|0:0|2026-09-25T11:02:10|2026-09-25T11:02:10|2026-09-25T11:03:21|71|24|24|1704|20:03.400||24412.50M|cpu=24,mem=131072M,node=1||m12||812.30M|96.10M
1001.extern|extern|COMPLETED|0:0|2026-09-25T11:02:10|2026-09-25T11:02:10|2026-09-25T11:03:21|71|24|24|1704|00:00.001||1.20M|billing=24,cpu=24,mem=131072M,node=1||m12|||
1002|call-a|FAILED|1:0|2026-09-25T09:00:00|2026-09-25T09:00:05|2026-09-25T09:10:05|600|24|24|14400|02:00:00|131072M||billing=24,cpu=24,mem=131072M,node=1|core|m13|2880||
1002.batch|batch|FAILED|1:0|2026-09-25T09:00:05|2026-09-25T09:00:05|2026-09-25T09:10:05|600|24|24|14400|02:00:00||1024M|cpu=24,mem=131072M,node=1||m13||10M|5M
2001|call-b-0|COMPLETED|0:0|2026-09-25T12:00:00|2026-09-25T12:00:30|2026-09-25T12:12:30|720|48|48|34560|1-02:03:04|192000M||billing=48,cpu=48,gres/gpu=4,gres/gpu:l40s=4,mem=192000M,node=1|gpu|m202|2880||
2001.batch|batch|COMPLETED|0:0|2026-09-25T12:00:30|2026-09-25T12:00:30|2026-09-25T12:12:30|720|48|48|34560|1-02:03:04||137216M|cpu=48,gres/gpu=4,mem=192000M,node=1||m202||4096M|2048M
2002|call-b-0|COMPLETED|0:0|2026-09-25T12:00:00|2026-09-25T12:00:30|2026-09-25T12:12:30|720|8|8|5760|00:10.000|2000Mc||billing=8,cpu=8,mem=16000M,node=1|core|m14|||
3001|call-c|CANCELLED by 1000|0:0|2026-09-25T13:00:00|2026-09-25T13:00:10|2026-09-25T13:00:40|30|4|4|120|00:01.000|8192M||billing=4,cpu=4,mem=8192M,node=1|core|m15|2880||
3002|call-d|RUNNING|0:0|2026-09-25T13:00:00|2026-09-25T13:00:10|Unknown|100|4|4|400||8192M||billing=4,cpu=4,mem=8192M,node=1|core|m15|2880||
"""
NAMES = [f.split("%", 1)[0] for f in slurm.SACCT_FIELDS]


def rows_of(text: str) -> list:
    return [dict(zip(NAMES, line.split("|"))) for line in text.splitlines() if line.strip()]


class ParseTest(unittest.TestCase):
    def test_seconds(self):
        self.assertEqual(accounting.parse_seconds("1-02:03:04"), 93784.0)
        self.assertEqual(accounting.parse_seconds("02:03:04"), 7384.0)
        self.assertAlmostEqual(accounting.parse_seconds("03:04.567"), 184.567)
        self.assertEqual(accounting.parse_seconds("71"), 71.0)
        for bad in ("", "UNLIMITED", "Partition_Limit", "Unknown", "x-1:2:3", "a:b"):
            self.assertIsNone(accounting.parse_seconds(bad), bad)

    def test_memory(self):
        self.assertEqual(accounting.parse_mem("25600M"), (25600.0, False))
        self.assertEqual(accounting.parse_mem("1.5G"), (1536.0, False))
        self.assertEqual(accounting.parse_mem("1024K"), (1.0, False))
        self.assertEqual(accounting.parse_mem("96000Mn"), (96000.0, False))
        self.assertEqual(accounting.parse_mem("2000Mc"), (2000.0, True))
        self.assertEqual(accounting.parse_mem("512"), (512.0, False))
        self.assertEqual(accounting.parse_mem(""), (None, False))
        self.assertIsNone(accounting.parse_mb("lots"))

    def test_tres(self):
        tres = accounting.parse_tres("billing=48,cpu=48,gres/gpu=4,gres/gpu:a100=4,mem=192000M,node=1")
        self.assertEqual(tres["mem"], "192000M")
        self.assertEqual(accounting.gpu_of_tres(tres), (4, "a100"))
        self.assertEqual(accounting.gpu_of_tres(accounting.parse_tres("cpu=8,mem=16000M")), (0, None))
        self.assertEqual(accounting.gpu_of_tres(accounting.parse_tres("gres/gpu:l40s=6")), (6, "l40s"))
        self.assertEqual(accounting.parse_tres(""), {})

    def test_fold(self):
        folded = accounting.fold(rows_of(SACCT_TEXT))
        self.assertEqual(sorted(folded), ["1001", "1002", "2001", "2002", "3001", "3002"])
        j = folded["1001"]
        self.assertEqual(j["parent"]["JobName"], "call-a")
        self.assertEqual(len(j["steps"]), 2)
        self.assertEqual((j["max_rss_mb"], j["disk_read_mb"], j["disk_write_mb"]), (24412.5, 812.3, 96.1))
        self.assertIsNone(folded["3002"]["max_rss_mb"])


class BuildTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.attempt = self.tmp / "attempt-1"
        for rel, jid in (("call-a", "1001"), ("call-a/failed1", "1002"), ("call-sub-0/call-b-0", "2001"),
                         ("call-sub-1/call-b-0", "2002"), ("call-c", "3001"), ("call-d", "3002")):
            d = self.attempt / rel
            d.mkdir(parents=True)
            (d / slurm.LOG_NAME).write_text(f"{jid};cluster\n" if jid == "1001" else f"{jid}\n")
        # the workflow log: `call-b-0` under call-sub-0 is an aliased call of task real_b; call-d has no setup line
        lines = []
        for name, src, rel in (("pbmm2_align_wgs", "wdl.w:ugc_wgw_x.t:call-a", "call-a"),
                               ("real_b", "wdl.w:ugc_wgw_x.w:call-sub-0.t:call-b-0", "call-sub-0/call-b-0"),
                               ("real_b", "wdl.w:ugc_wgw_x.w:call-sub-1.t:call-b-0", "call-sub-1/call-b-0"),
                               ("ctask", "wdl.w:ugc_wgw_x.t:call-c", None)):
            d = {"message": "task setup", "name": name, "source": src, "timestamp": 1.0, "level": "NOTICE"}
            if rel:
                d["dir"] = str(self.attempt / rel)
            lines.append(d)
            lines.append({"message": "done", "source": src, "timestamp": 2.0, "level": "NOTICE"})
        (self.attempt / "workflow.log.json").write_text("".join(json.dumps(l) + "\n" for l in lines))
        self.run = RunRecord(run_id="S1-x-a1", subject_type="sample", subject_id="S1", stage="singleton", mode="standalone",
                             ugc_wgw_version=VERSION, run_dir=str(self.attempt), status="success", attempt=1,
                             meta={"slurm_job_ids": ["1001", "1002", "2001", "2002", "3001", "3002", "4004"]})

    def test_job_refs(self):
        refs = {r.job_id: r for r in slurm.job_refs(self.attempt)}
        self.assertEqual((refs["1001"].rel_dir, refs["1001"].call_path, refs["1001"].call_id, refs["1001"].retry),
                         ("call-a", "call-a", "call-a", 0))
        self.assertEqual((refs["1002"].rel_dir, refs["1002"].call_path, refs["1002"].retry), ("call-a/failed1", "call-a", 1))
        self.assertEqual(refs["2001"].call_path, "call-sub-0/call-b-0")
        self.assertEqual(refs["2002"].call_path, "call-sub-1/call-b-0")
        self.assertEqual(slurm.job_ids(self.attempt), ["1001", "1002", "3001", "3002", "2001", "2002"])  # walk order
        self.assertEqual(slurm.call_base("call-pbmm2_align_wgs-3"), "pbmm2_align_wgs")
        self.assertEqual(slurm.call_base("call-deepvariant-03-chr3"), "deepvariant")
        self.assertEqual(slurm.call_base("call-glnexus-00"), "glnexus")
        self.assertEqual(slurm.call_base("hifiasm"), "hifiasm")

    def test_build_joins_dirs_rows_and_log(self):
        from ugc_wgw.wdllog import parse_workflow_log
        tasks, _ = parse_workflow_log(self.attempt / "workflow.log.json")
        acct = accounting.build(self.run, slurm.job_refs(self.attempt), rows_of(SACCT_TEXT), tasks,
                                self.run.meta["slurm_job_ids"], sacct_format="JobID,JobName")
        by = {j.job_id: j for j in acct.jobs}
        self.assertEqual(sorted(by), ["1001", "1002", "2001", "2002", "3001", "3002"])
        self.assertEqual(acct.missing, ["4004"])
        self.assertEqual(acct.status, "partial")   # a missing id and a running job
        a = by["1001"]
        self.assertEqual((a.task, a.task_source, a.call_path, a.retry, a.state, a.partition, a.nodes),
                         ("pbmm2_align_wgs", "workflow.log", "call-a", 0, "COMPLETED", "core", "m12"))
        self.assertEqual((a.queue_seconds, a.elapsed_seconds, a.alloc_cpus, a.req_cpus), (3730.0, 71.0, 24, 24))
        self.assertEqual((a.cpu_seconds_alloc, a.cpu_seconds_used), (1704.0, 1203.4))
        self.assertEqual((a.req_mem_mb, a.alloc_mem_mb, a.max_rss_mb), (131072.0, 131072.0, 24412.5))
        self.assertEqual((a.disk_read_mb, a.disk_write_mb, a.gpus, a.gpu_type, a.final, a.timelimit_minutes),
                         (812.3, 96.1, 0, None, True, 2880))
        self.assertEqual((by["1002"].task, by["1002"].retry, by["1002"].state), ("pbmm2_align_wgs", 1, "FAILED"))
        # the two call-b-0 are kept apart by their path; both name the aliased task from the log
        self.assertEqual((by["2001"].task, by["2001"].call_path, by["2001"].gpus, by["2001"].gpu_type),
                         ("real_b", "call-sub-0/call-b-0", 4, "l40s"))
        self.assertEqual(by["2001"].cpu_seconds_used, 93784.0)
        self.assertEqual((by["2002"].task, by["2002"].call_path, by["2002"].gpus), ("real_b", "call-sub-1/call-b-0", 0))
        self.assertEqual(by["2002"].req_mem_mb, 16000.0)   # 2000Mc × 8 CPUs
        self.assertEqual(by["2002"].alloc_mem_mb, 16000.0)
        self.assertEqual((by["3001"].task, by["3001"].state, by["3001"].final), ("ctask", "CANCELLED", True))  # matched by call path
        self.assertEqual((by["3002"].task, by["3002"].task_source, by["3002"].final, by["3002"].end), ("d", "call id", False, None))
        self.assertEqual(by["3002"].cpu_seconds_alloc, 400.0)
        self.assertTrue(any("not in a final state" in n for n in acct.notes))
        # round trip
        path = accounting.write(self.attempt, acct)
        self.assertEqual(path, layout.run_files(self.attempt).accounting)
        doc = json.loads(path.read_text())
        self.assertEqual((doc["schema"], doc["subject"], doc["status"]), (1, {"type": "sample", "id": "S1"}, "partial"))
        back = accounting.read(self.attempt)
        self.assertEqual([j.job_id for j in back.jobs], [j.job_id for j in acct.jobs])
        self.assertEqual(back.jobs[0].max_rss_mb, 24412.5)
        self.assertIsNone(accounting.read(self.tmp / "nowhere"))

    def test_build_without_rows_marks_everything_missing(self):
        acct = accounting.build(self.run, slurm.job_refs(self.attempt), [], [], [])
        self.assertEqual(acct.jobs, [])
        self.assertEqual(len(acct.missing), 6)
        self.assertEqual(acct.status, "partial")


class SacctTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.old = os.environ.copy()

    def tearDown(self):
        os.environ.clear()
        os.environ.update(self.old)

    def test_skipped_without_binary(self):
        os.environ["PATH"] = "/nonexistent"
        self.assertEqual(slurm.sacct(["1"])[0], "skipped")
        self.assertEqual(slurm.sacct([])[0], "skipped")

    def test_chunks_and_optional_field_retry(self):
        env, log = sacct_env(self.tmp, SACCT_TEXT)
        os.environ.update(env)
        status, rows, detail = slurm.sacct([str(i) for i in range(900)], chunk=400)
        self.assertEqual(status, "ok", detail)
        calls = log.read_text().strip().splitlines()
        self.assertEqual(len(calls), 3)
        self.assertTrue(all("-P --noheader --units=M --format=JobID,JobName%100," in c for c in calls))
        self.assertEqual(len(rows), 3 * 10)
        self.assertEqual(rows[0]["JobID"], "1001")
        self.assertEqual(rows[0]["MaxDiskWrite"], "")
        # an older sacct rejects TimelimitRaw: the chunk is asked again without the optional fields
        log.write_text("")
        os.environ["UGC_WGW_FAKE_SACCT_REJECT"] = "TimelimitRaw"
        status, rows, detail = slurm.sacct(["1", "2"])
        self.assertEqual(status, "ok", detail)
        calls = log.read_text().strip().splitlines()
        self.assertEqual(len(calls), 2)
        self.assertIn("TimelimitRaw", calls[0])
        self.assertNotIn("TimelimitRaw", calls[1])
        self.assertNotIn("MaxDiskWrite", calls[1])
        self.assertNotIn("MaxDiskWrite", rows[0])
        self.assertIn("fields JobID,JobName,State", detail)


class CaptureTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp)
        self.proj = str(self.cfg.project_dir)
        run_cli(["--project", self.proj, "samples", "add", str(write_tsv(self.tmp, [{"sample_id": "S1", "hifi_reads": "a.bam"}]))])

    def submit(self, extra_env: dict) -> RunRecord:
        env, _ = sacct_env(self.tmp, "31337|call-fake|COMPLETED|0:0|2026-09-25T10:00:00|2026-09-25T10:00:30|2026-09-25T10:05:30|300|8|8|2400|30:00.000|16000M||billing=8,cpu=8,mem=16000M,node=1|core|m1|120||\n"
                                     "31337.batch|batch|COMPLETED|0:0|2026-09-25T10:00:30|2026-09-25T10:00:30|2026-09-25T10:05:30|300|8|8|2400|30:00.000||2048M|cpu=8,mem=16000M,node=1||m1||1M|2M\n")
        env.update({"UGC_WGW_FAKE_SLURM_LOG": "31337", "UGC_WGW_FAKE_TASK_LOG": "1"})
        env.update(extra_env)
        code, _, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--samples", "S1"], env=env)
        self.assertEqual(code, 0, err)
        db = DB(self.cfg.db_path)
        run = db.latest_run("sample", "S1", "singleton", VERSION)
        self.events = db.events_for(run.run_id)
        db.close()
        return run

    def test_finalize_writes_accounting(self):
        run = self.submit({})
        acct = accounting.read(run.run_path)
        self.assertIsNotNone(acct)
        self.assertEqual((acct.status, len(acct.jobs), acct.missing), ("ok", 1, []))
        j = acct.jobs[0]
        self.assertEqual((j.job_id, j.task, j.task_source, j.call_path, j.alloc_cpus, j.elapsed_seconds, j.max_rss_mb),
                         ("31337", "fake_task", "workflow.log", "call-fake", 8, 300.0, 2048.0))
        ev = [e for e in self.events if e["event"] == "run.accounting"]
        self.assertEqual(len(ev), 1)
        self.assertEqual((ev[0]["detail"]["status"], ev[0]["detail"]["jobs"], ev[0]["level"]), ("ok", 1, "info"))
        self.assertEqual(json.loads(layout.run_files(run.run_path).accounting.read_text())["jobs"][0]["raw"]["JobName"], "call-fake")

    def test_disabled_and_no_sacct(self):
        self.cfg.accounting = False
        config.save(self.cfg)
        run = self.submit({})
        self.assertIsNone(accounting.read(run.run_path))
        self.assertEqual([e for e in self.events if e["event"] == "run.accounting"], [])
        # enabled, but no sacct on PATH: the run finishes, nothing is recorded (no file, no event)
        self.cfg.accounting = True
        config.save(self.cfg)
        run_cli(["--project", self.proj, "samples", "remove", "S1", "--force"])
        run_cli(["--project", self.proj, "samples", "add", str(write_tsv(self.tmp, [{"sample_id": "S1", "hifi_reads": "c.bam"}]))])
        code, _, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--samples", "S1"],
                               env={"PATH": "/usr/bin:/bin", "UGC_WGW_FAKE_SLURM_LOG": "31337"})
        self.assertEqual(code, 0, err)
        db = DB(self.cfg.db_path)
        run = db.latest_run("sample", "S1", "singleton", VERSION)
        ev = [e for e in db.events_for(run.run_id) if e["event"] == "run.accounting"]
        db.close()
        self.assertEqual(ev, [])
        self.assertIsNone(accounting.read(run.run_path))
        self.assertEqual(accounting.collect(self.cfg, run)[0], "skipped")


if __name__ == "__main__":
    unittest.main()
