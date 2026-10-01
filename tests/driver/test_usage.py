"""`ugc-wgw usage`: prices, the estimate from the workflow log, aggregation (clean vs as run, cohort shares), sizes,
the CLI tables, --collect/--refresh, and the config keys."""
import json
import os
import tempfile
import unittest
from pathlib import Path

from ugc_wgw import accounting, config, resources, usage
from ugc_wgw.db import DB, RunRecord
from ugc_wgw.util import UgcError
from ugc_wgw.wdllog import TaskStat

from .helpers import VERSION, make_project, run_cli, sacct_env, write_ids, write_tsv

SACCT_ROW = ("31337|call-fake|COMPLETED|0:0|2026-09-25T10:00:00|2026-09-25T10:00:30|2026-09-25T10:05:30|300|8|8|2400|30:00.000|16000M||"
             "billing=8,cpu=8,mem=16000M,node=1|core|m1|120||\n"
             "31337.batch|batch|COMPLETED|0:0|2026-09-25T10:00:30|2026-09-25T10:00:30|2026-09-25T10:05:30|300|8|8|2400|30:00.000||2048M|"
             "cpu=8,mem=16000M,node=1||m1||1M|2M\n")


def run(subject: str, stage: str, attempt: int, status: str, kind: str = "sample") -> RunRecord:
    return RunRecord(run_id=f"{subject}-{stage}-a{attempt}", subject_type=kind, subject_id=subject, stage=stage, mode="standalone",
                     ugc_wgw_version=VERSION, run_dir=f"/nowhere/{subject}/{stage}/attempt-{attempt}", status=status, attempt=attempt)


class PricesTest(unittest.TestCase):
    def setUp(self):
        self.cfg = make_project(Path(tempfile.mkdtemp()))

    def test_parse_and_format(self):
        self.cfg.prices = {"cpu_hour": 0.04, "currency": "EUR"}
        p = usage.parse_prices(self.cfg, ["gpu_hour=2.5", "storage_gb_month=0.02"])
        self.assertEqual((p.cpu_hour, p.gpu_hour, p.mem_gb_hour, p.storage_gb_month, p.currency), (0.04, 2.5, None, 0.02, "EUR"))
        self.assertTrue(p.any())
        self.assertEqual(p.describe(), "cpu_hour 0.04, gpu_hour 2.5, storage_gb_month 0.02 EUR")
        self.assertEqual(p.money(1234.5), "1,234.50 EUR")
        self.assertEqual(p.money(None), "–")
        self.assertFalse(usage.parse_prices(self.cfg, None).any() if not self.cfg.prices else False)
        self.assertEqual(usage.Prices().describe(), "none")
        for bad, msg in (("bogus=1", "unknown price"), ("cpu_hour=lots", "must be a number"), ("gpu_hour=-1", "negative"), ("cpu_hour", "number")):
            with self.assertRaises(UgcError) as cm:
                usage.parse_prices(self.cfg, [bad])
            self.assertIn(msg, str(cm.exception), bad)
        self.cfg.prices = {"nope": 1}
        with self.assertRaises(UgcError) as cm:
            usage.parse_prices(self.cfg)
        self.assertIn("config.json prices", str(cm.exception))

    def test_cost_math_and_basis(self):
        u = usage.Usage(core_h_alloc=10.0, core_h_req=5.0, gpu_h=2.0, mem_gb_h=100.0, out_bytes=10 * 1024 ** 3)
        p = usage.Prices(cpu_hour=0.1, gpu_hour=2.0, mem_gb_hour=0.01, storage_gb_month=0.02)
        self.assertAlmostEqual(u.cost(p, "allocated"), 1.0 + 4.0 + 1.0)
        self.assertAlmostEqual(u.cost(p, "requested"), 0.5 + 4.0 + 1.0)
        self.assertAlmostEqual(u.storage_per_month(p), 0.2)
        self.assertIsNone(u.cost(usage.Prices(), "allocated"))
        self.assertIsNone(usage.Usage().storage_per_month(p))

    def test_config_round_trip(self):
        self.cfg.prices = {"cpu_hour": 0.04}
        self.cfg.accounting = False
        self.cfg.accounting_timeout = 30.0
        config.save(self.cfg)
        back = config.load(self.cfg.project_dir)
        self.assertEqual((back.prices, back.accounting, back.accounting_timeout), ({"cpu_hour": 0.04}, False, 30.0))
        doc = json.loads(self.cfg.config_path.read_text())
        for k in ("prices", "accounting", "accounting_timeout"):
            doc.pop(k)
        self.cfg.config_path.write_text(json.dumps(doc))
        back = config.load(self.cfg.project_dir)
        self.assertEqual((back.prices, back.accounting, back.accounting_timeout), ({}, True, 120.0))


class EstimateTest(unittest.TestCase):
    def test_estimate_from_log_and_inventory(self):
        cfg = make_project(Path(tempfile.mkdtemp()))
        inv, gpu, notes = usage.inventory(cfg)
        self.assertEqual(notes, [])
        self.assertIn("pbmm2_align_wgs", inv)
        mos = inv["mosdepth"]
        tasks = [
            TaskStat("pbmm2_align_wgs", "call-pbmm2-0", 0.0, 600.0, cpu_policy=24, cpu_requested=32, cpu_granted=32),   # policy wins
            TaskStat("pbmm2_align_wgs", "call-pbmm2-1", 0.0, 600.0, cached=True),                                    # cached: skipped
            TaskStat("mosdepth", "call-mosdepth", 0.0, 1800.0),                                                     # inventory cpu
            TaskStat("deepvariant_call_variants_gpu", "call-dv", 0.0, 300.0, gres="gpu:a100:1"),                    # GPU from the log
            TaskStat("run_parabricks_deepvariant", "call-pb", 0.0, 3600.0, failed=True),                           # GPU count from the project
            TaskStat("unknown_task", "call-u", 0.0, 60.0),                                                          # no cpu known
            TaskStat("running", "call-r", 0.0, None),                                                               # unfinished: skipped
        ]
        total, per_task = usage.estimate(tasks, inv, gpu)
        self.assertEqual((total.jobs, total.failed_jobs, total.estimated_attempts), (5, 1, 0))
        self.assertAlmostEqual(per_task["pbmm2_align_wgs"].core_h_alloc, 24 * 600 / 3600)
        self.assertAlmostEqual(per_task["mosdepth"].core_h_alloc, (mos.cpu or 0) * 1800 / 3600)
        self.assertAlmostEqual(per_task["mosdepth"].mem_gb_h, (mos.memory or 0) / 1024 ** 3 * 0.5)
        self.assertAlmostEqual(per_task["deepvariant_call_variants_gpu"].gpu_h, 300 / 3600)
        self.assertAlmostEqual(per_task["run_parabricks_deepvariant"].gpu_h, 4.0)            # parabricks_gpus default 4 × 1 h
        self.assertAlmostEqual(per_task["run_parabricks_deepvariant"].wasted_core_h, per_task["run_parabricks_deepvariant"].core_h_alloc)
        self.assertEqual(per_task["unknown_task"].core_h_alloc, 0.0)
        self.assertNotIn("running", per_task)
        self.assertEqual(total.core_h_alloc, total.core_h_req)
        self.assertIsNone(total.efficiency)
        self.assertEqual(total.queue_s, [])


class AggregateTest(unittest.TestCase):
    def test_clean_vs_as_run_and_shares(self):
        def att(r: RunRecord, core: float, source: str = "sacct", status: str = "ok", **kw) -> usage.AttemptUsage:
            u = usage.Usage(jobs=1, attempts=1, core_h_alloc=core, core_h_req=core, **kw)
            return usage.AttemptUsage(r, source, status, u, {"t": usage.Usage(jobs=1, core_h_alloc=core)})
        attempts = [
            att(run("S1", "singleton", 1, "failed"), 10.0, wasted_core_h=10.0, failed_jobs=1),
            att(run("S1", "singleton", 2, "success"), 20.0, out_bytes=5),
            att(run("S2", "singleton", 1, "success"), 30.0, source="estimate", status="", out_bytes=7),
            att(run("C1", "cohort_merge", 1, "success", kind="cohort"), 6.0, out_bytes=1),
            att(run("S3", "singleton", 1, "cancelled"), 1.0, source="none", status=""),
        ]
        s = usage.aggregate(attempts)
        self.assertEqual((s.totals.attempts, s.totals.jobs, s.totals.failed_jobs), (5, 5, 1))
        self.assertAlmostEqual(s.totals.core_h_alloc, 67.0)
        self.assertAlmostEqual(s.clean.core_h_alloc, 56.0)              # latest successful per (subject, stage): 20 + 30 + 6
        self.assertEqual(s.n_samples, 2)                                 # S3 never succeeded
        self.assertAlmostEqual(s.per_sample_clean.core_h_alloc, 28.0)
        self.assertAlmostEqual(s.per_sample_as_run.core_h_alloc, 33.5)
        self.assertEqual(s.per_sample_clean.out_bytes, 6)                # (5 + 7 + 1) / 2, rounded
        self.assertEqual((s.per_sample_clean.jobs, s.per_sample_clean.attempts), (0, 0))   # counts are not per sample
        self.assertAlmostEqual(s.totals.wasted_core_h, 10.0)
        self.assertEqual(sorted(s.by_stage), ["cohort_merge", "singleton"])
        self.assertAlmostEqual(s.by_stage["singleton"].core_h_alloc, 61.0)
        self.assertAlmostEqual(s.by_task["t"].core_h_alloc, 67.0)
        self.assertEqual(s.by_subject[("cohort", "C1")].attempts, 1)
        self.assertEqual(dict(s.sources), {"sacct": 3, "estimate": 1, "none": 1})
        self.assertTrue(any("estimated" in n for n in s.notes))
        self.assertTrue(any("without task logs" in n for n in s.notes))
        self.assertEqual(usage.aggregate([]).n_samples, 0)
        self.assertIsNone(usage.aggregate([]).per_sample_clean)
        rows = usage.rows_subject(s, usage.Prices(cpu_hour=1.0, currency="X"), "allocated", sizes=True)
        self.assertEqual(rows[-2]["subject"], "mean per sample, clean")
        self.assertEqual((rows[-2]["core-h"], rows[-2]["cost"], rows[-1]["core-h"]), ("28.0", "28.00 X", "33.5"))
        stage_rows = usage.rows_stage(s, usage.Prices(), "allocated", mode="standalone", sizes=False)
        self.assertEqual([r["stage"] for r in stage_rows], ["singleton", "cohort_merge", "total"])
        self.assertNotIn("cost", stage_rows[0])
        self.assertNotIn("out", stage_rows[0])

    def test_from_accounting_and_dir_bytes(self):
        tmp = Path(tempfile.mkdtemp())
        acct = accounting.Accounting(run_id="r", attempt=1, stage="singleton", subject_type="sample", subject_id="S1",
                                     collected_at="now", status="ok", sacct_format="")
        acct.jobs.append(accounting.JobAcct(job_id="1", call_path="call-a", call_id="call-a", retry=0, task="pbmm2_align_wgs",
                                            task_source="workflow.log", job_name="call-a", state="COMPLETED", exit_code="0:0",
                                            partition="core", nodes="m1", submit=None, start=None, end=None, queue_seconds=120.0,
                                            elapsed_seconds=3600.0, timelimit_minutes=None, req_cpus=24, alloc_cpus=48,
                                            cpu_seconds_alloc=48 * 3600.0, cpu_seconds_used=24 * 3600.0, req_mem_mb=131072.0,
                                            alloc_mem_mb=131072.0, max_rss_mb=24000.0, gpus=2, gpu_type="a100", disk_read_mb=10.0,
                                            disk_write_mb=20.0, final=True))
        acct.jobs.append(accounting.JobAcct(job_id="2", call_path="call-b", call_id="call-b", retry=0, task="mosdepth",
                                            task_source="workflow.log", job_name="call-b", state="FAILED", exit_code="1:0",
                                            partition="core", nodes="m1", submit=None, start=None, end=None, queue_seconds=None,
                                            elapsed_seconds=1800.0, timelimit_minutes=None, req_cpus=None, alloc_cpus=4,
                                            cpu_seconds_alloc=None, cpu_seconds_used=None, req_mem_mb=8192.0, alloc_mem_mb=None,
                                            max_rss_mb=None, gpus=0, gpu_type=None, disk_read_mb=None, disk_write_mb=None, final=True))
        total, per_task = usage.from_accounting(acct)
        self.assertEqual((total.jobs, total.failed_jobs), (2, 1))
        self.assertAlmostEqual(total.core_h_alloc, 48.0 + 2.0)      # the second job: AllocCPUS × elapsed
        self.assertAlmostEqual(total.core_h_req, 24.0 + 2.0)        # ReqCPUS, else AllocCPUS
        self.assertAlmostEqual(total.cpu_h_used, 24.0)
        self.assertAlmostEqual(total.core_h_with_used, 48.0)
        self.assertAlmostEqual(total.efficiency, 50.0)
        self.assertAlmostEqual(total.gpu_h, 2.0)
        self.assertAlmostEqual(total.mem_gb_h, 128.0 + 4.0)
        self.assertEqual((total.peak_rss_mb, total.req_mem_mb_max, total.queue_s), (24000.0, 131072.0, [120.0]))
        self.assertAlmostEqual(total.wasted_core_h, 2.0)
        self.assertAlmostEqual(per_task["mosdepth"].mem_use or -1, -1)
        self.assertAlmostEqual(per_task["pbmm2_align_wgs"].mem_use, 100.0 * 24000 / 131072)
        # dir_bytes: a hardlinked file is counted once
        d = tmp / "attempt"
        (d / "out").mkdir(parents=True)
        (d / "work").mkdir()
        big = d / "work" / "big.bam"
        big.write_bytes(b"x" * 4096)
        os.link(big, d / "out" / "big.bam")
        (d / "out" / "small.txt").write_bytes(b"y" * 10)
        os.symlink(big, d / "out" / "link.bam")
        self.assertEqual(usage.dir_bytes(d / "out"), 4096 + 10 + len(str(big)))
        self.assertEqual(usage.dir_bytes(d), 4096 + 10 + len(str(big)))


class CliTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp)
        self.proj = str(self.cfg.project_dir)
        run_cli(["--project", self.proj, "samples", "add", str(write_tsv(self.tmp, [
            {"sample_id": "S1", "hifi_reads": "s1.bam"}, {"sample_id": "S2", "hifi_reads": "s2.bam"}]))])
        run_cli(["--project", self.proj, "cohort", "freeze", "C1", "--samples", str(write_ids(self.tmp, ["S1", "S2"]))])
        self.env, self.log = sacct_env(self.tmp, SACCT_ROW)
        self.env.update({"UGC_WGW_FAKE_SLURM_LOG": "31337", "UGC_WGW_FAKE_TASK_LOG": "1"})
        code, _, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--cohort", "C1"], env=self.env)
        self.assertEqual(code, 0, err)

    def test_tables_json_tsv_and_prices(self):
        code, out, err = run_cli(["--project", self.proj, "usage", "--json"])
        self.assertEqual(code, 0, err)
        doc = json.loads(out)
        self.assertEqual((doc["sources"], doc["n_samples"], doc["basis"]), ({"sacct": 4}, 2, "allocated"))
        self.assertEqual(len(doc["attempts"]), 4)                              # S1, S2 singletons, cohort_merge, cohort_freq
        self.assertAlmostEqual(doc["totals"]["core_h_alloc"], 4 * 8 * 300 / 3600, places=3)
        self.assertAlmostEqual(doc["totals"]["cpu_h_used"], 4 * 1800 / 3600, places=3)   # TotalCPU 30:00 per job
        self.assertEqual(doc["totals"]["efficiency_pct"], 75.0)
        self.assertEqual((doc["totals"]["queue_median_s"], doc["totals"]["peak_rss_mb"]), (30.0, 2048.0))
        self.assertAlmostEqual(doc["per_sample_clean"]["core_h_alloc"], 2 * 8 * 300 / 3600, places=3)
        self.assertIsNone(doc["totals"]["cost"])
        self.assertEqual(doc["by_task"]["fake_task"]["jobs"], 4)
        self.assertEqual(sorted(doc["by_stage"]), ["cohort_freq", "cohort_merge", "singleton"])
        code, out, err = run_cli(["--project", self.proj, "usage", "--price", "cpu_hour=0.5", "--price", "currency=EUR", "--sizes"])
        self.assertEqual(code, 0, err)
        self.assertIn("# 4 attempt(s) (sacct 4), 2 sample(s) with results; core-hours allocated 2.67", out)
        self.assertIn("cost 1.33 EUR at cpu_hour 0.5 EUR (0.67 EUR per sample, clean)", out)
        self.assertIn("queue wait median 30s max 30s", out)
        self.assertIn("results on disk", out)
        self.assertIn("# by stage", out)
        self.assertIn("# by task, top 20 by core-hours", out)
        self.assertIn("mean per sample, clean", out)
        self.assertRegex(out, r"singleton +2 +2 +0 +10m 00s +30s +30s +1\.33 +1\.33 +1\.00 +75% +0\.00")
        self.assertRegex(out, r"\ncost\b|  cost")
        code, out, err = run_cli(["--project", self.proj, "usage", "--by", "task", "--top", "1", "--tsv"])
        self.assertEqual(code, 0, err)
        lines = out.strip().splitlines()
        self.assertEqual(lines[0].split("\t")[:3], ["task", "jobs", "failed"])
        self.assertEqual(lines[1].split("\t")[:2], ["fake_task", "4"])
        self.assertEqual(len(lines), 2)
        code, out, err = run_cli(["--project", self.proj, "usage", "--by", "attempt", "--mode", "standalone", "--cohort", "C1"])
        self.assertEqual(code, 0, err)
        self.assertIn("sacct (ok)", out)
        self.assertEqual(out.count("\nS1 "), 1)
        code, out, err = run_cli(["--project", self.proj, "usage", "--samples", "S2", "--json"])
        self.assertEqual(json.loads(out)["totals"]["attempts"], 1)
        code, _, err = run_cli(["--project", self.proj, "usage", "--price", "bogus=1"])
        self.assertEqual(code, 1)
        self.assertIn("unknown price", err)
        code, _, err = run_cli(["--project", self.proj, "usage", "--cohort", "nope"])
        self.assertEqual(code, 1)
        self.assertIn("unknown cohort", err)

    def test_collect_refresh_and_estimate_fallback(self):
        db = DB(self.cfg.db_path)
        runs = db.all_runs(VERSION)
        db.close()
        paths = [accounting.read(r.run_path) and (r.run_path / "accounting.json") for r in runs]
        self.assertTrue(all(p and p.exists() for p in paths))
        paths[0].unlink()
        # no sacct on PATH: --collect refuses, plain usage falls back to the estimate for that attempt
        code, _, err = run_cli(["--project", self.proj, "usage", "--collect"], env={"PATH": "/usr/bin:/bin"})
        self.assertEqual(code, 1)
        self.assertIn("sacct, which is not on PATH", err)
        code, out, err = run_cli(["--project", self.proj, "usage", "--json"], env={"PATH": "/usr/bin:/bin"})
        self.assertEqual(code, 0, err)
        doc = json.loads(out)
        self.assertEqual(doc["sources"], {"sacct": 3, "estimate": 1})
        self.assertTrue(any("estimated from workflow.log" in n for n in doc["notes"]))
        self.assertIn("note: 1 attempt(s) estimated", err)
        # --collect backfills the missing file; --refresh re-reads all four
        self.log.write_text("")
        code, out, err = run_cli(["--project", self.proj, "usage", "--collect", "--json"], env=self.env)
        self.assertEqual(code, 0, err)
        self.assertIn("accounting: kept 3, ok 1", err)
        self.assertTrue(paths[0].exists())
        self.assertEqual(json.loads(out)["sources"], {"sacct": 4})
        self.assertEqual(len(self.log.read_text().strip().splitlines()), 1)
        code, _, err = run_cli(["--project", self.proj, "usage", "--refresh", "--by", "stage"], env=self.env)
        self.assertEqual(code, 0, err)
        self.assertIn("accounting: ok 4", err)
        self.assertEqual(len(self.log.read_text().strip().splitlines()), 5)
        # a partial record is re-collected without --refresh
        doc = json.loads(paths[1].read_text())
        doc["status"] = "partial"
        paths[1].write_text(json.dumps(doc))
        code, _, err = run_cli(["--project", self.proj, "usage", "--collect", "--by", "stage"], env=self.env)
        self.assertIn("accounting: kept 3, ok 1", err)


if __name__ == "__main__":
    unittest.main()
