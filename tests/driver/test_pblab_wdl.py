"""The pblab task file: the CPU and GPU variants of the pbdv task are one command, and the image is pinned."""
import re
import unittest

from .helpers import REPO

WDL = (REPO / "workflows" / "ugc_wgw" / "pblab" / "pblab.wdl").read_text()


def task(name: str) -> str:
    start = WDL.index(f"task {name} {{")
    nxt = WDL.find("\ntask ", start + 1)
    return WDL[start:nxt if nxt != -1 else len(WDL)]


def command(name: str) -> str:
    return re.search(r"command <<<\n(.*?)\n  >>>", task(name), re.S).group(1)


class PblabWdlTest(unittest.TestCase):
    def test_pbdv_cpu_and_gpu_are_one_command(self):
        """Two tasks only because the SLURM backend gives any task with gpuCount a GPU: the commands must not drift."""
        gpu, cpu = command("ugc_wgw_pbdv_call_gpu"), command("ugc_wgw_pbdv_call_cpu")
        self.assertIn("--device gpu", gpu)
        self.assertIn("--device cpu", cpu)
        strip = lambda text: [l for l in text.replace("--device gpu", "--device X").replace("--device cpu", "--device X").splitlines()
                              if "nvidia-smi" not in l]
        self.assertEqual(strip(gpu), strip(cpu))
        self.assertIn("gpuCount: gpuCount", task("ugc_wgw_pbdv_call_gpu"))
        self.assertNotIn("gpuCount", task("ugc_wgw_pbdv_call_cpu"))

    def test_outputs_are_upstreams_and_the_shims_are_there(self):
        for name in ("ugc_wgw_pbdv_call_gpu", "ugc_wgw_pbdv_call_cpu"):
            body = task(name)
            self.assertIn('String out_prefix = "~{sample_id}.~{ref_name}.small_variants"', body)
            cmd = command(name)
            self.assertIn("--exclude-uncalled", cmd)                 # as upstream's DeepVariant tasks end
            self.assertEqual(cmd.count("--remove FORMAT/PS"), 2)     # VCF and gVCF: unphased, no phase set
            self.assertEqual(cmd.count("--new-gt u"), 2)
            self.assertIn('--sample-name "~{sample_id}"', cmd)
        align = command("ugc_wgw_pbweave_align")
        for flag in ("--preset HIFI", '--sample "~{sample_id}"', "--sort", "--unmapped", "--min-length 50", "--max-memory"):
            self.assertIn(flag, align)
        self.assertIn("--remove-tag HP,PS,PC", align)                # pbweave has no --strip-tags

    def test_one_image_by_digest_in_the_manifest(self):
        pins = set(re.findall(r'docker: "~\{ugc_wgw_container_registry\}/(pblab@sha256:[0-9a-f]{64})"', WDL))
        self.assertEqual(len(pins), 1, pins)
        self.assertEqual(WDL.count("docker:"), 3)
        self.assertIn("/" + pins.pop(), (REPO / "image_manifest.ugc-wgw.txt").read_text())
