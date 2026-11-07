"""The entry kinds (0.16.0): a sample registered from reads, from an aligned BAM, or from an aligned BAM with its
VCFs; the sheet rules, the file checks, the per-member stage sequences, the builders, rebase, the pages and the
manifest (bin/ugc_wgw/samples.py, vcf.py, stages.py, plan.py, inputs.py)."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from ugc_wgw import cohorts, inputs, plan, samples, stages, vcf
from ugc_wgw.db import DB
from ugc_wgw.log import Events
from ugc_wgw.stages import declared_inputs, wdl_path
from ugc_wgw.util import UgcError

from .helpers import BGZF_EOF, REPO, VERSION, bgzf_block, make_project, run_cli, seed_success, write_bam, write_ids, write_tsv


def write_vcf(path: Path, *, samples_: tuple[str, ...] = ("S",), contigs: tuple[str, ...] = ("chr20",), bgzf: bool = True,
              records: int = 2, truncated: bool = False) -> Path:
    """A small BGZF VCF with a header (contigs, a source line), `samples_` columns and a few records; its .tbi."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["##fileformat=VCFv4.2", "##source=DeepVariant_v1.9.0"]
    lines += [f"##contig=<ID={c},length=64444167>" for c in contigs]
    lines.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
    lines.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples_]))
    for i in range(records):
        lines.append("\t".join([contigs[0], str(1000 + i), ".", "A", "G", "30", "PASS", ".", "GT", *(["0/1"] * len(samples_))]))
    data = ("\n".join(lines) + "\n").encode()
    out = b"".join(bgzf_block(data[i:i + 60000]) for i in range(0, len(data), 60000)) + BGZF_EOF if bgzf else __import__("gzip").compress(data)
    if truncated:
        out = out[: len(out) * 2 // 3]
    path.write_bytes(out)
    Path(str(path) + ".tbi").write_bytes(b"TBI\x01" + bytes(28))
    return path


def write_aligned_bam(path: Path, **kw: object) -> Path:
    p = write_bam(path, aligned=True, **kw)
    Path(str(p) + ".bai").write_bytes(b"BAI\x01" + bytes(8))
    return p


class ProjectCase(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp)
        self.proj = str(self.cfg.project_dir)
        self.db = DB(self.cfg.db_path)
        self.events = Events(self.db, self.cfg.logs_dir / "events.jsonl")
        self.data = self.tmp / "data"
        self.data.mkdir(exist_ok=True)

    def tearDown(self):
        self.db.close()

    def files(self, sid: str, *, gvcf: bool = True, vcf_samples: tuple[str, ...] = ("S",)) -> dict[str, str]:
        bam = write_aligned_bam(self.data / f"{sid}.aligned.bam", sample=sid)
        row = {"aligned_bam": str(bam),
               "small_variant_vcf": str(write_vcf(self.data / f"{sid}.small.vcf.gz", samples_=vcf_samples)),
               "sv_vcf": str(write_vcf(self.data / f"{sid}.sv.vcf.gz", samples_=vcf_samples))}
        if gvcf:
            row["gvcf"] = str(write_vcf(self.data / f"{sid}.g.vcf.gz", samples_=vcf_samples))
        return row

    def sheet(self, rows):
        return write_tsv(self.tmp, rows)


class VcfInspectTest(unittest.TestCase):
    def test_inspect(self):
        tmp = Path(tempfile.mkdtemp())
        info = vcf.inspect(write_vcf(tmp / "a.vcf.gz", samples_=("HG002",), contigs=("chr20", "chr21")))
        self.assertTrue(info.ok, info.problems)
        self.assertEqual((info.samples, info.n_contigs, info.contigs, info.fileformat), (["HG002"], 2, ["chr20", "chr21"], "VCFv4.2"))
        self.assertEqual(info.sources, ["source=DeepVariant_v1.9.0"])
        self.assertTrue(info.bgzf and info.eof_ok)
        two = vcf.inspect(write_vcf(tmp / "b.vcf.gz", samples_=("A", "B")))
        self.assertIn("2 sample column(s); exactly one is required", two.problems)
        sites = vcf.inspect(write_vcf(tmp / "c.vcf.gz", samples_=()))
        self.assertIn("a sites-only VCF (no sample column)", sites.problems)
        plain = vcf.inspect(write_vcf(tmp / "d.vcf.gz", bgzf=False))
        self.assertTrue(any("not BGZF-compressed" in p for p in plain.problems), plain.problems)
        cut = vcf.inspect(write_vcf(tmp / "e.vcf.gz", truncated=True))
        self.assertTrue(any("truncated" in p or "unreadable" in p for p in cut.problems), cut.problems)
        (tmp / "f.vcf.gz").write_bytes(b"")
        self.assertIn("empty file", vcf.inspect(tmp / "f.vcf.gz").problems[0])


class SheetTest(ProjectCase):
    def test_kinds_indices_and_listing(self):
        v = self.files("V1")
        b = {"aligned_bam": str(write_aligned_bam(self.data / "B1.bam", sample="B1"))}
        tsv = self.sheet([{"sample_id": "R1", "sex": "MALE", "hifi_reads": "r1.bam"},
                          {"sample_id": "B1", "sex": "FEMALE", **b},
                          {"sample_id": "V1", **v}])
        n = samples.add_samples(self.db, self.events, tsv)
        self.assertEqual(n, 3)
        r1, b1, v1 = (self.db.get_sample(s) for s in ("R1", "B1", "V1"))
        self.assertEqual((r1.entry_kind, b1.entry_kind, v1.entry_kind), ("reads", "bam", "vcf"))
        self.assertEqual(b1.aligned_bam_index, str(self.data / "B1.bam.bai"))          # found by convention
        self.assertEqual(v1.small_variant_vcf_index, str(self.data / "V1.small.vcf.gz.tbi"))
        self.assertTrue(v1.gvcf and v1.gvcf_index)
        self.assertEqual(b1.hifi_reads, [])
        self.assertIn(b1.aligned_bam, b1.input_info)
        self.assertTrue(b1.input_info[b1.aligned_bam]["aligned"])
        self.assertEqual(v1.input_info[v1.small_variant_vcf]["samples"], ["S"])
        self.assertEqual(v1.meta["input_check"]["entry_kind"], "vcf")
        self.assertEqual(sorted(v1.meta["input_check"]["vcfs"]), ["gvcf", "small_variant_vcf", "sv_vcf"])
        warnings = [e for e in self.db.all_events() if e["event"] == "sample.warning"]
        self.assertTrue(any("sample column S is renamed to the sample id" in str(e["detail"]) for e in warnings), warnings)
        rows = samples.list_samples(self.db, "standalone", VERSION)
        self.assertEqual({r["sample_id"]: r["entry_kind"] for r in rows}, {"R1": "reads", "B1": "bam", "V1": "vcf"})
        code, out, err = run_cli(["--project", self.proj, "samples", "check", "V1", "--tsv"])
        self.assertEqual(code, 0, err)
        self.assertIn("small_variant_vcf", out)
        # a sheet with only aligned_bam in the header is a sheet
        tsv2 = self.tmp / "bams.tsv"
        tsv2.write_text("sample_id\taligned_bam\nB2\t" + str(write_aligned_bam(self.data / "B2.bam", sample="B2")) + "\n")
        self.assertEqual(samples.add_samples(self.db, self.events, tsv2), 1)
        self.assertEqual(self.db.get_sample("B2").entry_kind, "bam")

    def test_refusals(self):
        v = self.files("V1")
        cases = [
            ({"sample_id": "X", "hifi_reads": "x.bam", **{"aligned_bam": v["aligned_bam"]}}, "hifi_reads and aligned_bam are both given"),
            ({"sample_id": "X", "small_variant_vcf": v["small_variant_vcf"], "sv_vcf": v["sv_vcf"]}, "VCFs need the aligned_bam"),
            ({"sample_id": "X", "aligned_bam": v["aligned_bam"], "small_variant_vcf": v["small_variant_vcf"]}, "small_variant_vcf without sv_vcf"),
            ({"sample_id": "X", "aligned_bam": v["aligned_bam"], "gvcf": v["gvcf"]}, "gvcf given without small_variant_vcf and sv_vcf"),
            ({"sample_id": "X", "aligned_bam": v["aligned_bam"], "fail_reads": "xf.bam"}, "fail_reads are only taken with hifi_reads"),
            ({"sample_id": "X", "sex": "MALE"}, "no input: give hifi_reads, or aligned_bam"),
            ({"sample_id": "X", "aligned_bam": v["aligned_bam"], "sv_vcf_index": v["sv_vcf"] + ".tbi"}, "sv_vcf_index without sv_vcf"),
        ]
        for row, text in cases:
            with self.assertRaises(UgcError, msg=text) as cm:
                samples.add_samples(self.db, self.events, self.sheet([row]))
            self.assertIn(text, str(cm.exception))
        # an unaligned BAM under aligned_bam, a two-sample VCF, a missing index
        ubam = write_bam(self.data / "u.bam")
        Path(str(ubam) + ".bai").write_bytes(b"BAI\x01")
        with self.assertRaises(UgcError) as cm:
            samples.add_samples(self.db, self.events, self.sheet([{"sample_id": "U", "aligned_bam": str(ubam)}]))
        self.assertIn("not aligned (no @SQ lines in the header)", str(cm.exception))
        two = self.files("T", vcf_samples=("A", "B"))
        with self.assertRaises(UgcError) as cm:
            samples.add_samples(self.db, self.events, self.sheet([{"sample_id": "T", **two}]))
        self.assertIn("2 sample column(s)", str(cm.exception))
        noidx = write_bam(self.data / "n.bam", aligned=True, sample="N")
        with self.assertRaises(UgcError) as cm:
            samples.add_samples(self.db, self.events, self.sheet([{"sample_id": "N", "aligned_bam": str(noidx)}]))
        self.assertIn("no index next to n.bam", str(cm.exception))
        # the read group's sample name must be the sample id (sawfish and TRGT name their outputs after it); no read group at all
        other = write_aligned_bam(self.data / "o.bam", sample="OTHER")
        with self.assertRaises(UgcError) as cm:
            samples.add_samples(self.db, self.events, self.sheet([{"sample_id": "O", "aligned_bam": str(other)}]))
        self.assertIn("@RG SM is OTHER, not the sample id", str(cm.exception))
        self.assertIn("samtools reheader", str(cm.exception))
        norg = write_aligned_bam(self.data / "g.bam", read_group=False)
        with self.assertRaises(UgcError) as cm:
            samples.add_samples(self.db, self.events, self.sheet([{"sample_id": "G", "aligned_bam": str(norg)}]))
        self.assertIn("no @RG SM in the header", str(cm.exception))
        self.assertIsNone(self.db.get_sample("O"))
        # --no-check trusts the conventional index name
        tsv = self.tmp / "nocheck.tsv"
        tsv.write_text("sample_id\taligned_bam\nQ\t/nowhere/q.bam\n")
        samples.add_samples(self.db, self.events, tsv, check_paths=False)
        self.assertEqual(self.db.get_sample("Q").aligned_bam_index, "/nowhere/q.bam.bai")


class SequenceTest(unittest.TestCase):
    def test_member_sequences(self):
        self.assertEqual(stages.sample_stages("standalone", "bam"), ("upstream", "downstream"))
        self.assertEqual(stages.sample_stages("joint", "vcf"), ("intake", "downstream"))
        self.assertEqual(stages.sample_stages("assembly", "vcf"), ())
        self.assertEqual(stages.member_sequence("standalone", "vcf"), ("intake", "downstream", "cohort_qc", "cohort_merge", "cohort_freq", "cohort_methyl"))
        self.assertEqual(stages.member_sequence("joint", "vcf", has_gvcf=False, allow_standalone=True),
                         ("intake", "cohort_qc", "downstream", "cohort_merge", "cohort_freq", "cohort_methyl"))
        self.assertEqual(stages.member_sequence("joint", "vcf", has_gvcf=False), stages.member_sequence("joint", "vcf"))
        with self.assertRaises(UgcError):
            stages.sample_stages("standalone", "fastq")


class BuildersAndRunsTest(ProjectCase):
    """A reads, a bam and a vcf member through the fake in both modes: the plan, the builders, the pages, the manifest."""

    def setUp(self):
        super().setUp()
        v = self.files("V1")
        b = {"aligned_bam": str(write_aligned_bam(self.data / "B1.bam", sample="B1"))}
        samples.add_samples(self.db, self.events, self.sheet([{"sample_id": "R1", "sex": "MALE", "hifi_reads": "r1.bam"},
                                                             {"sample_id": "B1", "sex": "FEMALE", **b},
                                                             {"sample_id": "V1", "sex": "MALE", **v}]))
        self.cohort = cohorts.freeze(self.db, self.events, "C1", write_ids(self.tmp, ["R1", "B1", "V1"]))

    def ctx(self, subject_type, subject_id, mode="standalone", **kw):
        return inputs.BuildContext(self.cfg, self.db, mode, VERSION, subject_type, subject_id, cohort=self.cohort, **kw)

    def check_doc(self, stage, doc):
        declared = declared_inputs(wdl_path(REPO, stage))
        ns = f"ugc_wgw_{stage}."
        plain = {k[len(ns):] for k in doc if "." not in k[len(ns):]}
        self.assertTrue(plain <= set(declared), plain - set(declared))
        self.assertTrue({n for n, r in declared.items() if r} <= plain, {n for n, r in declared.items() if r} - plain)
        return {k[len(ns):]: v for k, v in doc.items()}

    def test_builders(self):
        d = self.check_doc("intake", inputs.generate(self.ctx("sample", "V1"), "intake")[0])
        self.assertTrue(d["aligned_bam"].endswith("V1.aligned.bam") and d["aligned_bam_index"].endswith(".bai"))
        self.assertTrue(d["small_variant_vcf_in"].endswith("V1.small.vcf.gz") and d["sv_vcf_in_index"].endswith(".tbi"))
        self.assertTrue(d["gvcf_in"].endswith("V1.g.vcf.gz"))
        self.assertEqual(d["declared_sex"], "MALE")
        d = self.check_doc("upstream", inputs.generate(self.ctx("sample", "B1"), "upstream")[0])
        self.assertEqual(d["hifi_reads"], [])
        self.assertTrue(d["aligned_bam_in"].endswith("B1.bam") and d["aligned_bam_in_index"].endswith("B1.bam.bai"))
        self.assertNotIn("fail_reads", d)
        d = self.check_doc("upstream", inputs.generate(self.ctx("sample", "R1", "joint"), "upstream")[0])
        self.assertNotIn("aligned_bam_in", d)
        with self.assertRaises(UgcError):
            inputs.generate(self.ctx("sample", "B1"), "singleton")
        with self.assertRaises(UgcError):
            inputs.generate(self.ctx("sample", "V1"), "upstream")
        with self.assertRaises(UgcError):
            inputs.generate(self.ctx("sample", "R1"), "intake")
        # downstream reads intake for the vcf member, upstream for the bam member (standalone)
        seed_success(self.cfg, self.db, "intake", "sample", "V1", {"sample_id": "V1", "gvcf_in": "/g.vcf.gz", "inferred_sex": "MALE"})
        seed_success(self.cfg, self.db, "upstream", "sample", "B1", {"sample_id": "B1", "aligned_bam_in": "/b.bam", "inferred_sex": "FEMALE"})
        d = self.check_doc("downstream", inputs.generate(self.ctx("sample", "V1"), "downstream")[0])
        self.assertIn("/intake/", d["aligned_hifi_reads"])
        self.assertIn("/intake/", d["small_variant_vcf"])
        self.assertEqual(d["sex"], "MALE")
        d = self.check_doc("downstream", inputs.generate(self.ctx("sample", "B1"), "downstream")[0])
        self.assertIn("/upstream/", d["aligned_hifi_reads"])
        # the entry record of the manifest
        rec = inputs.entry_record(self.ctx("sample", "V1"))
        self.assertEqual(rec["kind"], "vcf")
        self.assertEqual([x["role"] for x in rec["inputs"]], ["aligned_bam", "aligned_bam_index", "small_variant_vcf", "small_variant_vcf_index",
                                                             "sv_vcf", "sv_vcf_index", "gvcf", "gvcf_index"])
        self.assertEqual(rec["inputs"][2]["info"]["samples"], ["S"])
        self.assertEqual(inputs.entry_record(self.ctx("sample", "R1"))["inputs"], [])

    def test_standalone_run_and_pages(self):
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--cohort", "C1"])
        self.assertEqual(code, 0, err)
        self.assertIn("mitochondrial merge is left out", err)   # V1 has no mitorsaw call
        latest = {(r.subject_id, r.stage): r.status for r in self.db.all_runs(VERSION)}
        self.assertEqual(latest[("R1", "singleton")], "success")
        self.assertEqual((latest[("B1", "upstream")], latest[("B1", "downstream")]), ("success", "success"))
        self.assertEqual((latest[("V1", "intake")], latest[("V1", "downstream")]), ("success", "success"))
        self.assertNotIn(("B1", "singleton"), latest)
        for stage in ("cohort_qc", "cohort_merge", "cohort_freq", "cohort_methyl"):
            self.assertEqual(latest[("C1", stage)], "success", stage)
        # cohort_freq carries the kinds; cohort_merge's inputs took the bam member's downstream and the vcf member's intake
        freq_run = self.db.latest_success("cohort", "C1", "cohort_freq", VERSION)
        freq_in = json.loads(Path(freq_run.inputs_path).read_text())
        self.assertEqual(freq_in["ugc_wgw_cohort_freq.sample_entry_kinds"], ["reads", "bam", "vcf"])
        merge_in = json.loads(Path(self.db.latest_success("cohort", "C1", "cohort_merge", VERSION).inputs_path).read_text())
        self.assertFalse(merge_in["ugc_wgw_cohort_merge.merge_mitochondria"])
        self.assertTrue(all("/downstream/" in p or "/singleton/" in p for p in merge_in["ugc_wgw_cohort_merge.sv_vcfs"]))
        man = json.loads((Path(self.db.latest_success("sample", "V1", "intake", VERSION).run_dir) / "run_manifest.json").read_text())
        self.assertEqual(man["entry"]["kind"], "vcf")
        self.assertEqual(man["entry"]["inputs"][0]["role"], "aligned_bam")
        man_b = json.loads((Path(self.db.latest_success("sample", "B1", "upstream", VERSION).run_dir) / "run_manifest.json").read_text())
        self.assertEqual(man_b["entry"]["kind"], "bam")
        # the pages: the external_input marker and the overview note
        code, out, err = run_cli(["--project", self.proj, "summary", "--cohort", "C1"])
        self.assertEqual(code, 0, err)
        out_dir = self.cfg.results_dir / "reports" / "summary"
        dv = json.loads((out_dir / "V1.summary.json").read_text())
        self.assertEqual(dv["sample"]["entry_kind"], "vcf")
        self.assertIn("external_input", {f["code"] for f in dv["qc"]["flags"]})
        self.assertIn("input_checks", dv["files"])
        self.assertIn("External input", (out_dir / "V1.summary.html").read_text())
        dr = json.loads((out_dir / "R1.summary.json").read_text())
        self.assertNotIn("external_input", {f["code"] for f in dr["qc"]["flags"]})
        cd = json.loads((out_dir / "C1.summary.json").read_text())
        self.assertEqual({r["id"]: r["present"] for r in cd["members"]}, {"R1": True, "B1": True, "V1": True})
        # samples list shows the kinds and the statuses of each member's own stages
        code, out, err = run_cli(["--project", self.proj, "samples", "list", "--json"])
        self.assertEqual(code, 0, err)
        rows = {r["sample_id"]: r for r in json.loads(out)}
        self.assertEqual(rows["V1"]["entry_kind"], "vcf")

    def test_joint_gvcf_rule(self):
        # V2: a vcf member without a gVCF
        v2 = self.files("V2", gvcf=False)
        samples.add_samples(self.db, self.events, self.sheet([{"sample_id": "V2", **v2}]))
        c2 = cohorts.freeze(self.db, self.events, "C2", write_ids(self.tmp, ["R1", "V2"], "c2.txt"))
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "joint", "--cohort", "C2"])
        self.assertEqual(code, 0, err)
        self.assertIn("joint mode needs a gVCF from every member; vcf member(s) without one: V2", err)
        latest = {(r.subject_id, r.stage): r.status for r in self.db.all_runs(VERSION)}
        self.assertNotIn(("C2", "cohort_call"), latest)
        self.assertEqual((latest[("R1", "upstream")], latest[("V2", "intake")]), ("success", "success"))
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "joint", "--cohort", "C2", "--allow-standalone-members"])
        self.assertEqual(code, 0, err)
        self.assertIn("1 vcf member(s) without a gVCF left out of the joint call", err)
        latest = {(r.subject_id, r.stage): r.status for r in self.db.all_runs(VERSION)}
        self.assertEqual(latest[("C2", "cohort_call")], "success")
        self.assertEqual((latest[("R1", "downstream")], latest[("V2", "downstream")]), ("success", "success"))
        cc_in = json.loads(Path(self.db.latest_success("cohort", "C2", "cohort_call", VERSION).inputs_path).read_text())
        self.assertEqual(cc_in["ugc_wgw_cohort_call.sample_ids"], ["R1"])
        ds_in = json.loads(Path(self.db.latest_success("sample", "V2", "downstream", VERSION).inputs_path).read_text())
        self.assertIn("/intake/", ds_in["ugc_wgw_downstream.small_variant_vcf"])
        man = json.loads((Path(self.db.latest_success("cohort", "C2", "cohort_call", VERSION).run_dir) / "run_manifest.json").read_text())
        self.assertTrue(man["allow_standalone_members"])
        for stage in ("cohort_merge", "cohort_freq", "cohort_methyl"):
            self.assertEqual(latest[("C2", stage)], "success", stage)

    def test_rebase(self):
        seed_success(self.cfg, self.db, "singleton", "sample", "R1", {"sample_id": "R1", "inferred_sex": "MALE"})
        run = self.db.latest_success("sample", "R1", "singleton", VERSION)
        code, out, err = run_cli(["--project", self.proj, "samples", "rebase", "R1", "--from-run", run.run_id])
        self.assertEqual(code, 0, err)
        self.assertIn("entry kind bam", err)
        rec = self.db.get_sample("R1")
        self.assertEqual((rec.entry_kind, rec.hifi_reads), ("bam", []))
        self.assertTrue(rec.aligned_bam.endswith("haplotagged.bam"))
        self.assertEqual(rec.meta["rebase"]["run_id"], run.run_id)
        self.assertEqual(len(rec.meta["rebase"]["previous_hifi_reads"]), 1)
        self.assertEqual(stages.member_sequence("standalone", rec.entry_kind)[:2], ("upstream", "downstream"))
        code, out, err = run_cli(["--project", self.proj, "samples", "rebase", "R1", "--from-run", run.run_id, "--with-vcfs"])
        self.assertEqual(code, 0, err)
        rec = self.db.get_sample("R1")
        self.assertEqual(rec.entry_kind, "vcf")
        self.assertTrue(rec.small_variant_vcf.endswith("small_variants.phased.vcf.gz") and rec.gvcf.endswith(".g.vcf.gz"))
        code, out, err = run_cli(["--project", self.proj, "samples", "rebase", "R1", "--from-run", "nope"])
        self.assertEqual(code, 1)
        self.assertIn("not a run of sample R1", err)


if __name__ == "__main__":
    unittest.main()
