"""Resource escalation (escalate.py): what the next attempt carries after a time or memory failure, the plugin's
environment, and the learned requests."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from ugc_wgw import escalate, failures
from ugc_wgw.db import DB, RunRecord
from ugc_wgw.engine import RunResult
from ugc_wgw.log import Events

from .helpers import VERSION, make_project

G = 1024 ** 3


def run_record(tmp: Path, attempt: int, meta: dict | None = None, status: str = "failed") -> RunRecord:
    run_dir = tmp / f"attempt-{attempt}"
    run_dir.mkdir(parents=True, exist_ok=True)
    return RunRecord(run_id=f"S1-assembly_call-a{attempt}-x", subject_type="sample", subject_id="S1", stage="assembly_call",
                     mode="assembly", ugc_wgw_version=VERSION, run_dir=str(run_dir), status=status, attempt=attempt,
                     meta=dict(meta or {}))


def write_log(run_dir: Path, lines: list[dict]) -> None:
    (run_dir / "workflow.log.json").write_text("".join(json.dumps(l) + "\n" for l in lines))


SRC = "wdl.w:ugc_wgw_assembly_call.t:call-ugc_wgw_pav"


class EscalateTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp)

    def classification(self, resource: str | None, node: str = "call-ugc_wgw_pav"):
        return failures.Classification("resource", "CommandFailed", "x", {"node": node}, resource)

    def test_task_name_from_the_log(self):
        run = run_record(self.tmp, 1)
        write_log(run.run_path, [{"message": "task setup", "name": "ugc_wgw_pav", "source": SRC, "timestamp": 1.0, "level": "NOTICE"},
                                 {"message": "task setup", "name": "glnexus", "source": "wdl.w:ugc_wgw_cohort_call.w:call-ugc_wgw_glnexus_scatter.t:call-glnexus-07",
                                  "timestamp": 1.0, "level": "NOTICE"}])
        self.assertEqual(escalate.task_name_for(run.run_path, "call-ugc_wgw_pav"), "ugc_wgw_pav")
        self.assertEqual(escalate.task_name_for(run.run_path, "call-glnexus-07"), "glnexus")
        self.assertIsNone(escalate.task_name_for(run.run_path, "call-nope"))
        self.assertIsNone(escalate.task_name_for(run.run_path, None))

    def test_plan_doubles_then_exhausts(self):
        run = run_record(self.tmp, 1)
        write_log(run.run_path, [{"message": "task setup", "name": "ugc_wgw_pav", "source": SRC, "timestamp": 1.0, "level": "NOTICE"}])
        esc, why = escalate.plan(self.cfg, run, self.classification("memory"))
        self.assertEqual((esc["task"], esc["resource"], esc["factor"], esc["n"], esc["max"]), ("ugc_wgw_pav", "memory", 2.0, 1, 2))
        self.assertEqual(why, "")
        self.assertEqual(escalate.describe(esc), "ugc_wgw_pav memory ×2 (escalation 1 of 2)")
        # the second attempt ran with ×2 and failed again: ×4 is the last one
        run2 = run_record(self.tmp, 2, meta={"escalations": escalate.carry(RunRecord(**{**run.__dict__, "meta": {"escalation": esc}}))})
        self.assertEqual(run2.meta["escalations"], {"ugc_wgw_pav": {"memory_factor": 2.0}})
        write_log(run2.run_path, [{"message": "task setup", "name": "ugc_wgw_pav", "source": SRC, "timestamp": 1.0, "level": "NOTICE"}])
        esc2, _ = escalate.plan(self.cfg, run2, self.classification("memory"))
        self.assertEqual((esc2["factor"], esc2["n"]), (4.0, 2))
        run3 = run_record(self.tmp, 3, meta={"escalations": {"ugc_wgw_pav": {"memory_factor": 4.0}}})
        write_log(run3.run_path, [{"message": "task setup", "name": "ugc_wgw_pav", "source": SRC, "timestamp": 1.0, "level": "NOTICE"}])
        esc3, why3 = escalate.plan(self.cfg, run3, self.classification("memory"))
        self.assertIsNone(esc3)
        self.assertIn("escalation exhausted", why3)
        self.assertIn("×4", why3)
        # a time failure of the same task is its own ladder
        esc4, _ = escalate.plan(self.cfg, run3, self.classification("time"))
        self.assertEqual((esc4["resource"], esc4["factor"], esc4["n"]), ("time", 2.0, 1))

    def test_plan_declines_without_a_resource_a_task_or_a_budget(self):
        run = run_record(self.tmp, 1)
        write_log(run.run_path, [{"message": "task setup", "name": "ugc_wgw_pav", "source": SRC, "timestamp": 1.0, "level": "NOTICE"}])
        esc, why = escalate.plan(self.cfg, run, self.classification(None))
        self.assertIsNone(esc)
        self.assertIn("refused or the resource is unknown", why)
        esc, why = escalate.plan(self.cfg, run, self.classification("memory", node="call-unknown"))
        self.assertIsNone(esc)
        self.assertIn("failed task is unknown", why)
        self.cfg.resource_retry_max = 0
        esc, why = escalate.plan(self.cfg, run, self.classification("memory"))
        self.assertIsNone(esc)
        self.assertIn("resource_retry_max is 0", why)
        self.assertEqual(escalate.plan(self.cfg, run, failures.Classification("tool", "CommandFailed", "x", {})), (None, ""))

    def test_env_merges_learned_minimums_and_run_factors(self):
        run = run_record(self.tmp, 2, meta={"escalations": {"ugc_wgw_pav": {"memory_factor": 2.0}}})
        self.assertEqual(escalate.env_for(self.cfg, run_record(self.tmp, 1)), {})
        env = escalate.env_for(self.cfg, run)
        self.assertEqual(json.loads(env[escalate.ESCALATION_ENV]), {"ugc_wgw_pav": {"memory_factor": 2.0}})
        escalate.save_learned(self.cfg, {"ugc_wgw_dipcall": {"memory": 64 * G, "time_minutes": 480, "run_id": "r"},
                                         "ugc_wgw_pav": {"memory": 100 * G}})
        env = escalate.env_for(self.cfg, run)
        self.assertEqual(json.loads(env[escalate.ESCALATION_ENV]),
                         {"ugc_wgw_dipcall": {"memory_min": 64 * G, "time_min": 480},
                          "ugc_wgw_pav": {"memory_factor": 2.0, "memory_min": 100 * G}})

    def test_learn_from_finished_escalated_tasks_only(self):
        db = DB(self.cfg.db_path)
        events = Events(db, self.cfg.logs_dir / "events.jsonl")
        run = run_record(self.tmp, 2, meta={"escalations": {"ugc_wgw_pav": {"memory_factor": 2.0}, "ugc_wgw_dipcall": {"time_factor": 2.0}}})
        dip = "wdl.w:ugc_wgw_assembly_call.t:call-ugc_wgw_dipcall"
        write_log(run.run_path, [
            {"message": "task setup", "name": "ugc_wgw_pav", "source": SRC, "timestamp": 1.0, "level": "NOTICE"},
            {"message": escalate.ESCALATION_MESSAGE, "task": "ugc_wgw_pav", "summary": "memory 40G→80G",
             "changes": {"memory_reservation": [40 * G, 80 * G]}, "source": SRC, "timestamp": 1.5, "level": "NOTICE"},
            {"message": "done", "source": SRC, "timestamp": 100.0, "level": "NOTICE"},
            {"message": "task setup", "name": "ugc_wgw_dipcall", "source": dip, "timestamp": 1.0, "level": "NOTICE"},
            {"message": escalate.ESCALATION_MESSAGE, "task": "ugc_wgw_dipcall", "summary": "time 04:00:00→08:00:00",
             "changes": {"time_minutes": [240, 480]}, "source": dip, "timestamp": 1.5, "level": "NOTICE"},
            {"message": "task ugc_wgw_dipcall failed", "source": dip, "timestamp": 50.0, "level": "ERROR"},
        ])
        lessons = escalate.learn(self.cfg, events, run)
        self.assertEqual(lessons, ["ugc_wgw_pav memory 85899345920"])   # dipcall failed: nothing learned from it
        learned = escalate.load_learned(self.cfg)
        self.assertEqual(learned["ugc_wgw_pav"]["memory"], 80 * G)
        self.assertEqual(learned["ugc_wgw_pav"]["run_id"], run.run_id)
        self.assertNotIn("ugc_wgw_dipcall", learned)
        self.assertTrue(any(json.loads(l)["event"] == "resources.learned" for l in (self.cfg.logs_dir / "events.jsonl").read_text().splitlines()))
        # a smaller later value never lowers the lesson; a run without escalations teaches nothing
        write_log(run.run_path, [
            {"message": "task setup", "name": "ugc_wgw_pav", "source": SRC, "timestamp": 1.0, "level": "NOTICE"},
            {"message": escalate.ESCALATION_MESSAGE, "task": "ugc_wgw_pav", "summary": "memory 40G→60G",
             "changes": {"memory_reservation": [40 * G, 60 * G]}, "source": SRC, "timestamp": 1.5, "level": "NOTICE"},
            {"message": "done", "source": SRC, "timestamp": 100.0, "level": "NOTICE"}])
        self.assertEqual(escalate.learn(self.cfg, events, run), [])
        self.assertEqual(escalate.learn(self.cfg, events, run_record(self.tmp, 1)), [])
        self.cfg.learn_resources = False
        self.assertEqual(escalate.learn(self.cfg, events, run), [])
        db.close()

    def test_policy_rows(self):
        from ugc_wgw import resources
        policy = resources.policy_module(self.cfg.code_dir)
        self.assertEqual(escalate.policy_rows(self.cfg, policy), [])
        escalate.save_learned(self.cfg, {"ugc_wgw_pav": {"memory": 80 * G}, "ugc_wgw_dipcall": {"time_minutes": 480, "memory": 64 * G}})
        rows = escalate.policy_rows(self.cfg, policy)
        self.assertEqual(rows, ["ugc_wgw_dipcall\t-\t64G\t08:00:00\t-\t-", "ugc_wgw_pav\t-\t80G\t-\t-\t-"])
