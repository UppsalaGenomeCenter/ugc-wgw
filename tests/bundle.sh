#!/usr/bin/env bash
# Bundle round-trip: make-bundle.sh (--images data when the dev machine's SIF cache already holds the
# reference data container, so the installer's copy-and-verify step runs; --images none otherwise)
# -> install-bundle.sh --verify-only and a real install into a scratch prefix, then `ugc-wgw init --install`;
# then, with fake SIFs and no Apptainer, a full bundle A, a delta bundle B against A's manifest, their installs
# into one prefix (the SIF store: hardlinks, adoption from an older layout, --dedupe-sif, the refusal of a delta
# on a site without its images).
# Needs network for pip download (dev machine or CI). Usage: tests/bundle.sh [<scratch dir>]
# UGC_WGW_BUNDLE_IMAGES=none|data|install forces the images mode (CI has no SIF cache and no Apptainer: `none`).
# A scratch dir it made itself (about 8 GB) is removed after a passing run and kept, with its path printed, after
# a failing one or when UGC_WGW_KEEP_SCRATCH=1; a scratch dir given as the argument is always kept.
set -euo pipefail
cd "$(dirname "$0")/.."
keep=${UGC_WGW_KEEP_SCRATCH:-}
if [ $# -ge 1 ]; then
  scratch=$1; keep=1
else
  scratch=$(mktemp -d "${TMPDIR:-/tmp}/ugc-bundle-test.XXXXXX")
fi
cleanup() {
  rc=$?
  if [ "$rc" -eq 0 ] && [ -z "$keep" ]; then
    rm -rf "$scratch"
  else
    echo "== scratch kept at $scratch" >&2
  fi
}
trap cleanup EXIT
echo "== bundle test in $scratch"
images=none; bundle_args=()
if ls bundle/out/sif-cache/*workflow-data-container-hifi-human-wgs-wdl-grch38_giabv3*.sif >/dev/null 2>&1; then
  images=data; bundle_args=(--sif-cache "$(pwd)/bundle/out/sif-cache")
  # with the pbmm2 image cached too, the installer's index build runs (about a minute, 12 GB of memory)
  if ls bundle/out/sif-cache/*pbmm2@sha256_86e39aa6*.sif >/dev/null 2>&1; then images=install; fi
fi
case "${UGC_WGW_BUNDLE_IMAGES:-}" in
  '') ;;
  none) images=none; bundle_args=() ;;
  data|install) images=$UGC_WGW_BUNDLE_IMAGES ;;
  *) echo "UGC_WGW_BUNDLE_IMAGES must be none, data or install" >&2; exit 2 ;;
esac
echo "== images mode $images"
# --without-ancestry: the 1000 Genomes somalier set (83 MB from Zenodo) stays out of the test bundles
tarball=$(scripts/make-bundle.sh --version HEAD --out "$scratch/out" --images "$images" --without-ancestry "${bundle_args[@]}")
[ -f "$tarball" ] && [ -f "$tarball.sha256" ]
scripts/install-bundle.sh --bundle "$tarball" --prefix "$scratch/prefix" --verify-only 2>&1 | tee "$scratch/verify.log" >&2
grep -q 'plugin: ugc-wgw-miniwdl' "$scratch/verify.log"
grep -q 'reference data container:' "$scratch/verify.log"
# a site file that keeps --partition in SLURM_EXTRA_ARGS next to SLURM_PARTITION must be refused before anything is installed
printf 'SLURM_PARTITION=core\nSLURM_EXTRA_ARGS=--partition core\n' > "$scratch/bad-site.cfg"
if scripts/install-bundle.sh --bundle "$tarball" --prefix "$scratch/prefix-bad" --site "$scratch/bad-site.cfg" >/dev/null 2>&1; then
  echo "installer accepted --partition in SLURM_EXTRA_ARGS next to SLURM_PARTITION" >&2; exit 1
fi
[ ! -d "$scratch/prefix-bad/versions" ]
printf 'SLURM_PARTITION=core\nSLURM_ACCOUNT=acc-1\nSLURM_EXTRA_ARGS=--qos short\nTASK_CPU_MAX=48\nTASK_MEMORY_MAX=360G\nTASK_RESOURCES=<prefix>/resources.tsv\nSLURM_PARTITION_GPU=gpu\nSINGULARITY_NV=1\nTASK_SCRATCH_ENV=TMPDIR\n' > "$scratch/site.cfg"
report=$(scripts/install-bundle.sh --bundle "$tarball" --prefix "$scratch/prefix" --site "$scratch/site.cfg" --activate 2> "$scratch/install.log")
cat "$scratch/install.log" >&2
install_dir=$(echo "$report" | sed -n 's/^install_dir=//p')
[ -x "$install_dir/venv/bin/miniwdl" ]
echo "$report" | grep -q '^profiles='
[ -f "$scratch/prefix/profiles/cpu.json" ] && [ -f "$scratch/prefix/profiles/parabricks.json" ]
[ -f "$install_dir/miniwdl.cfg" ] && ! grep -q '{{' "$install_dir/miniwdl.cfg"
[ "$(readlink "$scratch/prefix/current")" = "versions/$(cat VERSION)" ]
grep -q 'command_shell = /bin/bash' "$install_dir/miniwdl.cfg"
grep -q '^cpu_max = 48$' "$install_dir/miniwdl.cfg"
grep -q '^memory_max = 386547056640$' "$install_dir/miniwdl.cfg"   # 360 GiB in bytes
grep -q '^defaults = {"docker": "ubuntu:20.04", "maxRetries": 1, "time_minutes": 4320, "slurm_partition": "core", "slurm_account": "acc-1", "slurm_partition_gpu": "gpu"}$' "$install_dir/miniwdl.cfg"
grep -q '^run_options = \["--containall", "--no-mount", "hostfs", "--env", "APPTAINER_BIND=", "--env", "SINGULARITY_BIND=", "--nv"\]$' "$install_dir/miniwdl.cfg"
grep -q '^extra_args = "--qos short"$' "$install_dir/miniwdl.cfg"
grep -q "^resources = $scratch/prefix/resources.tsv$" "$install_dir/miniwdl.cfg"
grep -q '^container_backend = ugc_slurm_singularity$' "$install_dir/miniwdl.cfg"
grep -q '^scratch_env = TMPDIR$' "$install_dir/miniwdl.cfg"
grep -q '^scratch_mount = /mnt/scratch$' "$install_dir/miniwdl.cfg"
[ -f "$scratch/prefix/resources.tsv" ] && grep -q '^task	cpu	memory	time	partition	constraint$' "$scratch/prefix/resources.tsv"
echo "$report" | grep -q '^plugin=ugc-wgw-miniwdl '
echo "$report" | grep -q "^sif_store=$scratch/prefix/sif$"
# manifest schema 2: every image carries `shipped` and is in files[] exactly when shipped; no base for a full bundle
python3 - "$install_dir/manifest.json" <<'PYEOF'
import json, sys
m = json.load(open(sys.argv[1]))
assert m["schema"] == 2, m["schema"]
assert m["base"] is None and m["images_total"] == len(m["images"]) and m["images_shipped"] == len(m["images"])
paths = {f["path"] for f in m["files"]}
for i in m["images"]:
    assert isinstance(i["shipped"], bool) and len(i["sha256"]) == 64 and isinstance(i["size"], int), i
    assert i["shipped"] == (("sif/" + i["sif"]) in paths), i["sif"]
PYEOF
# every SIF of the version is a hardlink into the content-addressed store
for f in "$install_dir"/sif/*.sif; do
  [ -f "$f" ] || continue
  sha=$(sha256sum "$f" | cut -d' ' -f1)
  [ "$(stat -c %i "$f")" = "$(stat -c %i "$scratch/prefix/sif/$sha.sif")" ] || { echo "$f is not linked into the store" >&2; exit 1; }
  [ "$(stat -c %h "$f")" -ge 2 ]
done
if ls "$install_dir"/sif/*.sif >/dev/null 2>&1; then
  [ "$(grep -vc '^#' "$scratch/prefix/sif/INDEX.tsv")" -ge 1 ]
fi
# an invalid scratch variable name is refused before anything is installed
# shellcheck disable=SC2016  # a literal dollar on purpose: variable names are not expanded by the installer
echo 'TASK_SCRATCH_ENV=$TMPDIR' > "$scratch/bad-scratch.cfg"
if scripts/install-bundle.sh --bundle "$tarball" --prefix "$scratch/prefix-bad2" --site "$scratch/bad-scratch.cfg" >/dev/null 2>&1; then
  echo "installer accepted TASK_SCRATCH_ENV=\$TMPDIR" >&2; exit 1
fi
[ ! -d "$scratch/prefix-bad2/versions" ]
# the node-local scratch backend: registered, selected by the rendered cfg, its exec script shipped in the wheel and
# working: a fake apptainer shows the bind and TMPDIR inserted after `exec`. The class is resolved through miniwdl's
# plugin registry but instantiated directly: task_container.new() would run the backend's global_init, which probes
# `apptainer --version`, and CI has no Apptainer
"$install_dir/venv/bin/python" - "$install_dir/miniwdl.cfg" "$scratch" <<'PYEOF'
import logging, os, pathlib, stat, subprocess, sys
from WDL.runtime.config import Loader
from WDL.runtime import config as wdlcfg
from ugc_wgw_miniwdl import scratch
cfg = Loader(logging.getLogger("bundle-test"), filenames=[sys.argv[1]])
assert cfg["scheduler"]["container_backend"] == "ugc_slurm_singularity"
backends = dict(wdlcfg.load_plugins(cfg, "container_backend"))
backend_cls = backends[cfg["scheduler"]["container_backend"]]
assert backend_cls is scratch.SlurmSingularityScratch, backend_cls
host = pathlib.Path(sys.argv[2]) / "backend-host"; host.mkdir()
backend = backend_cls(cfg, "call-x-1", str(host))
args = backend._slurm_invocation()
assert args[-4:] == [scratch.EXEC_SCRIPT, "TMPDIR", "/mnt/scratch", "--"], args[-4:]
assert os.path.isfile(scratch.EXEC_SCRIPT) and "site-packages" in scratch.EXEC_SCRIPT, scratch.EXEC_SCRIPT
fake = pathlib.Path(sys.argv[2]) / "fake-apptainer"
fake.write_text('#!/usr/bin/env bash\nfor a in "$@"; do printf "%s\\n" "$a"; done\n'); fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
job_scratch = pathlib.Path(sys.argv[2]) / "job-scratch"; job_scratch.mkdir()
res = subprocess.run(["bash", scratch.EXEC_SCRIPT, "TMPDIR", "/mnt/scratch", "--", str(fake), "exec", "--pwd", "/w", "img", "cmd"],
                     env={"PATH": os.environ["PATH"], "TMPDIR": str(job_scratch)}, capture_output=True, text=True)
assert res.returncode == 0, res.stderr
assert res.stdout.splitlines()[:5] == ["exec", "--bind", f"{job_scratch}:/mnt/scratch", "--env", "TMPDIR=/mnt/scratch"], res.stdout   # the fake prints its arguments
assert "scratch:" in res.stderr, res.stderr
PYEOF
# the hardlink-else-symlink wrap is in place at interpreter start (ugc_wgw_miniwdl.pth), before any plugin loads
ls "$install_dir"/venv/lib/python*/site-packages/ugc_wgw_miniwdl.pth >/dev/null
"$install_dir/venv/bin/python" -c 'import sys, WDL._util as u; sys.exit(0 if getattr(u.symlink_force, "ugc_wgw_crossdev", False) else 1)'
# the rendered cfg parses as miniwdl config and the plugin sees the policy path through it
"$install_dir/venv/bin/python" - "$install_dir/miniwdl.cfg" "$scratch/prefix/resources.tsv" <<'PYEOF'
import logging, sys
from WDL.runtime.config import Loader
from ugc_wgw_miniwdl import resources
cfg = Loader(logging.getLogger("t"), filenames=[sys.argv[1]])
assert resources.policy_path(cfg) == sys.argv[2], resources.policy_path(cfg)
assert cfg["task_runtime"].get_int("cpu_max") == 48
assert cfg["task_runtime"].get_dict("defaults")["slurm_partition"] == "core"
assert cfg["task_runtime"].get_dict("defaults")["slurm_partition_gpu"] == "gpu"
assert cfg["singularity"].get_list("run_options") == ["--containall", "--no-mount", "hostfs", "--env", "APPTAINER_BIND=", "--env", "SINGULARITY_BIND=", "--nv"]
PYEOF
ref_map=$install_dir/references/ugc_wgw_ref_map.GRCh38_GIABv3.tsv
grep -q "^name	GRCh38_GIABv3$" "$ref_map"
grep -q "^scatter_regions	$scratch/prefix/references/ugc-wgw-extras-0.3.0/scatter_regions.GRCh38_GIABv3.tsv$" "$ref_map"
mmi=$scratch/prefix/references/ugc-wgw-extras-0.3.0/human_GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.mmi
grep -q "^pbmm2_index	$mmi$" "$ref_map"
grep -q "^run_starphase	true$" "$ref_map"
[ -f "$scratch/prefix/references/ugc-wgw-extras-0.3.0/scatter_regions.GRCh38_GIABv3.tsv" ]
# the identity references (references.lock ugc-wgw-qc): the two shipped files installed, the fetched set left out and said so
qc=$scratch/prefix/references/ugc-wgw-qc-0.1.0
grep -q "^somalier_sites_vcf	$qc/somalier_sites.hg38.vcf.gz$" "$ref_map"
grep -q "^somalier_1kg_labels	$qc/somalier_ancestry_labels_1kg.tsv$" "$ref_map"
grep -q "^somalier_1kg_tar	$qc/1kg.somalier.tar.gz$" "$ref_map"
[ "$(sha256sum "$qc/somalier_sites.hg38.vcf.gz" | cut -d' ' -f1)" = d1a853b8bb2e5f1a520bc67c8303be699543d225d785bce51e8335a2420489b5 ]
[ -f "$qc/somalier_ancestry_labels_1kg.tsv" ] && [ ! -f "$qc/1kg.somalier.tar.gz" ]
grep -q '1kg.somalier.tar.gz not in the bundle (make-bundle.sh --without-ancestry' "$scratch/install.log"
# the methylation references (references.lock ugc-wgw-methyl): three shipped files, the map's header names its inputs
me=$scratch/prefix/references/ugc-wgw-methyl-0.1.0
grep -q "^methbat_report_regions_tsv	$me/methbat_report_regions.hg38.tsv$" "$ref_map"
grep -q "^methbat_map	$me/methbat_map.GRCh38_GIABv3.tsv$" "$ref_map"
grep -q "^methbat_atlas_tsv	$me/nanomix_39Bisulfite.methbat.tsv$" "$ref_map"
[ "$(sha256sum "$me/methbat_report_regions.hg38.tsv" | cut -d' ' -f1)" = 7d8b0a5e9751e90e877c3bd61c651214363f3ee5920679ef4fcfe47dd21e19f3 ]
grep -q '^##methbat_map_input=gtf	gencode.v50.annotation.gtf.gz	sha256:83fba3e9' "$me/methbat_map.GRCh38_GIABv3.tsv"
[ "$(grep -vc '^#' "$me/methbat_map.GRCh38_GIABv3.tsv")" -eq 32039 ]   # header + 32,038 islands
[ -f "$me/nanomix_39Bisulfite.methbat.tsv" ]
# the tandem-repeat references (references.lock ugc-wgw-tr): two generated tables, both reference-map keys
trr=$scratch/prefix/references/ugc-wgw-tr-0.1.0
grep -q "^tr_disease_loci_tsv	$trr/tr_disease_loci.strchive-2.4.4.tsv$" "$ref_map"
grep -q "^tr_named_overlaps_tsv	$trr/tr_named_overlaps.tsv$" "$ref_map"
grep -q '^##tr_map_input=strchive	STRchive-loci.v2.4.4.json@v2.4.4	sha256:bbfe8823' "$trr/tr_disease_loci.strchive-2.4.4.tsv"
grep -q 'CC BY 4.0' "$trr/tr_disease_loci.strchive-2.4.4.tsv"
[ "$(grep -vc '^#' "$trr/tr_disease_loci.strchive-2.4.4.tsv")" -eq 74 ]   # header + 73 named loci
[ "$(grep -vc '^#' "$trr/tr_named_overlaps.tsv")" -eq 63 ]                # header + 62 twins
# the smoke's truth data (references.lock ugc-wgw-truth): six shipped files, no reference-map key
tr=$scratch/prefix/references/ugc-wgw-truth-0.1.0
[ "$(sha256sum "$tr/HG002_GRCh38_v5.0q_smvar.chr20.vcf.gz" | cut -d' ' -f1)" = 8cb4332ba078d77e4720c22c5de912a3cac933988953611549b30abadb7edcdd ]
[ -s "$tr/HG002_GRCh38_v5.0q_stvar.chr20.vcf.gz.tbi" ] && [ -s "$tr/HG002_GRCh38_v5.0q_stvar.benchmark.chr20.bed" ]
if grep -q 'ugc-wgw-truth' "$ref_map"; then echo "the truth files must not be reference-map keys" >&2; exit 1; fi
# the assembly stages' references (references.lock ugc-wgw-asm): three shipped files, the derived fasta and its index built
asm=$scratch/prefix/references/ugc-wgw-asm-0.1.0
grep -q "^asm_fasta	$asm/human_GRCh38_GIABv3_asm.fasta$" "$ref_map"
grep -q "^asm_fasta_index	$asm/human_GRCh38_GIABv3_asm.fasta.fai$" "$ref_map"
grep -q "^asm_par_bed	$asm/asm_par.GRCh38.bed$" "$ref_map"
grep -q "^asm_graph_chromosomes	$asm/asm_graph_chromosomes.GRCh38.txt$" "$ref_map"
[ "$(wc -l < "$asm/asm_keep_contigs.GRCh38_GIABv3.txt")" -eq 194 ]
[ "$(wc -l < "$asm/asm_graph_chromosomes.GRCh38.txt")" -eq 24 ]
[ "$(md5sum "$asm/asm_par.GRCh38.bed" | cut -d' ' -f1)" = 25276f1418127b3971f3d84495639b29 ]
echo "$report" | grep -q "^ref_map=$ref_map$"
if [ "$images" = data ] || [ "$images" = install ]; then
  # the tree was copied out of the container and every file verified against references.lock
  tree=$scratch/prefix/references/hifi-wdl-resources-v4.0.0-GRCh38_GIABv3
  [ -f "$tree/manifest.json" ] && [ -f "$tree/GRCh38_GIABv3/trgt/adotto_strchive_20250827.hg38.bed.gz" ]
  echo "$report" | grep -q "^references=$tree$"
  if [ "$images" = install ]; then
    # the reference's pbmm2 index and the assembly stages' fasta (+ .fai) were built with the bundled images and match references.lock
    [ -f "$mmi" ]
    [ "$(md5sum "$mmi" | cut -d' ' -f1)" = "$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["bundles"]["ugc-wgw-extras"]["built"][0]["md5"])' references.lock)" ]
    [ "$(md5sum "$asm/human_GRCh38_GIABv3_asm.fasta" | cut -d' ' -f1)" = 1210a749eee62f28de18ba906ae7198f ]
    [ "$(md5sum "$asm/human_GRCh38_GIABv3_asm.fasta.fai" | cut -d' ' -f1)" = 374a75f82ff6f6ae919649aa0322cc20 ]
    [ "$(wc -l < "$asm/human_GRCh38_GIABv3_asm.fasta.fai")" -eq 194 ]
    grep -q 'human_GRCh38_GIABv3_asm.fasta built and verified' "$scratch/install.log"
    # every path the map names exists, except the optional 1000 Genomes set left out by --without-ancestry (asserted absent above)
    while read -r f; do [ -f "$f" ] || { echo "map names a missing file: $f" >&2; exit 1; }; done < <(cut -f2 "$ref_map" | grep '^/' | grep -v '/1kg.somalier.tar.gz$')
  else
    [ ! -f "$mmi" ] && [ ! -f "$asm/human_GRCh38_GIABv3_asm.fasta" ]   # no build images in the bundle: the installer says how to build them and the map still names them
    grep -q 'not built: the bundle carries no quay.io/pacbio/pbmm2' "$scratch/install.log"
  fi
  # a second install must reuse the tree (no copy), verify it, keep the index, and link every SIF from the store (0 added)
  scripts/install-bundle.sh --bundle "$tarball" --prefix "$scratch/prefix" --site "$scratch/site.cfg" --force > "$scratch/reinstall.log" 2>&1
  grep -q 'reference tree present' "$scratch/reinstall.log"
  grep -q 'profile kept' "$scratch/reinstall.log"
  grep -q 'SIF store .*: 0 added' "$scratch/reinstall.log"
  if [ "$images" = install ]; then grep -q '\.mmi present and verified' "$scratch/reinstall.log" && grep -q 'asm.fasta present and verified' "$scratch/reinstall.log"; else grep -q '\.mmi not built' "$scratch/reinstall.log"; fi
fi
if grep -rq 'v3p1p0\|tertiary' "$install_dir/references" "$install_dir/inputs"; then echo "v3 map names survive in the install" >&2; exit 1; fi
# WDL read_map() needs exactly two fields per line: no comments, no blank lines, in any rendered map
for m in "$install_dir"/references/*.tsv; do
  if grep -qE '^#|^$' "$m" || awk -F'\t' 'NF != 2 { bad = 1 } END { exit bad }' "$m"; then :; else echo "rendered map is not a two-column TSV: $m" >&2; exit 1; fi
done
"$install_dir/code/bin/ugc-wgw" init "$scratch/proj" --install "$scratch/prefix/current" --ref-map "$ref_map"
grep -q '"miniwdl_cfg"' "$scratch/proj/.ugc-wgw/config.json"
# the rendered inputs templates point at files the installer really wrote
for f in $(grep -ho '"[^"]*/references/[^"]*"' "$install_dir"/inputs/*.json | tr -d '"' | sort -u); do
  [ -f "$f" ] || { echo "rendered template points at a missing file: $f" >&2; exit 1; }
done

# ---- delta bundles and the SIF store, with fake SIFs (no Apptainer needed: the ugc-wgw image set has no built[]
# image and no data container, every reference is digest-form, and a cache that holds every SIF pulls nothing)
echo "== delta bundles"
fake=$scratch/fake-sif
mkdir -p "$fake"
while read -r image; do
  [ -n "$image" ] || continue
  uri="docker://$image"; uri="${uri//\//_}"; uri="${uri//:/_}"
  printf 'fake sif of %s\n' "$image" > "$fake/$uri.sif"
done < <(sed 's/#.*//' image_manifest.ugc-wgw.txt | awk 'NF { print $1 }' | sort -u)
n_fake=$(find "$fake" -maxdepth 1 -name '*.sif' | wc -l)
[ "$n_fake" -ge 4 ]
tarA=$(scripts/make-bundle.sh --version HEAD --out "$scratch/outA" --images ugc-wgw --without-ancestry --sif-cache "$fake")
first_fake=$(find "$fake" -maxdepth 1 -name '*.sif' | sort | head -1)
[ -f "$first_fake.sha256" ]   # the builder wrote the sidecars
manA="${tarA%.tar}.manifest.json"
[ -f "$manA" ]
scripts/install-bundle.sh --bundle "$tarA" --prefix "$scratch/pfxA" > "$scratch/installA.log" 2>&1
[ "$(find "$scratch/pfxA/sif" -maxdepth 1 -name '*.sif' | wc -l)" -eq "$n_fake" ]
[ "$(grep -vc '^#' "$scratch/pfxA/sif/INDEX.tsv")" -eq "$n_fake" ]
grep -q "SIF store $scratch/pfxA/sif: $n_fake added" "$scratch/installA.log"
# one image changes; the delta against A ships only that one
changed=$(find "$fake" -maxdepth 1 -name '*hifiasm*.sif' | head -1)
printf 'fake sif of hifiasm, rebuilt\n' > "$changed"
tarB=$(scripts/make-bundle.sh --version HEAD --out "$scratch/outB" --images ugc-wgw --base "$manA" --without-ancestry --sif-cache "$fake")
case "$tarB" in *-delta.tar) ;; *) echo "delta tar not named -delta: $tarB" >&2; exit 1 ;; esac
tar -tf "$tarB" > "$scratch/tarB.list"   # (not piped into grep -q: tar would die of SIGPIPE under pipefail)
grep -q '/DELTA-BUNDLE$' "$scratch/tarB.list"
python3 - "${tarB%.tar}.manifest.json" "$manA" "$n_fake" <<'PYEOF'
import json, sys
m = json.load(open(sys.argv[1])); a = json.load(open(sys.argv[2])); n = int(sys.argv[3])
assert m["schema"] == 2 and m["images_total"] == n and m["images_shipped"] == 1, (m["images_total"], m["images_shipped"])
assert m["base"]["name"] == a["name"] and m["base"]["images"] == n, m["base"]
assert sum(1 for f in m["files"] if f["path"].startswith("sif/")) == 1
assert [i["sif"] for i in m["images"] if i["shipped"]][0].find("hifiasm") >= 0
PYEOF
[ "$(stat -c %s "$tarB")" -lt 20000000 ]
# verify-only: fine against the prefix that holds A, refused against a fresh one, fine with --sif-source
scripts/install-bundle.sh --bundle "$tarB" --prefix "$scratch/pfxA" --verify-only > "$scratch/verifyB.log" 2>&1
grep -q "images: $n_fake required; 1 shipped; $((n_fake - 1)) in the SIF store" "$scratch/verifyB.log"
if scripts/install-bundle.sh --bundle "$tarB" --prefix "$scratch/pfxFresh" --verify-only > "$scratch/verifyB-fresh.log" 2>&1; then
  echo "verify-only accepted a delta bundle on a site without its images" >&2; exit 1
fi
grep -q "$((n_fake - 1)) missing" "$scratch/verifyB-fresh.log"
[ ! -d "$scratch/pfxFresh/versions" ] && [ ! -d "$scratch/pfxFresh/sif" ]
scripts/install-bundle.sh --bundle "$tarB" --prefix "$scratch/pfxFresh" --verify-only --sif-source "$fake" > "$scratch/verifyB-src.log" 2>&1
grep -q "$((n_fake - 1)) adoptable" "$scratch/verifyB-src.log"
# install B over A: one SIF added, the rest linked from the store; the old hifiasm content is now an orphan
scripts/install-bundle.sh --bundle "$tarB" --prefix "$scratch/pfxA" --force > "$scratch/installB.log" 2>&1
grep -q "SIF store $scratch/pfxA/sif: 1 added .*, $((n_fake - 1)) linked from the store" "$scratch/installB.log"
[ "$(find "$scratch/pfxA/sif" -maxdepth 1 -name '*.sif' | wc -l)" -eq $((n_fake + 1)) ]
[ "$(find "$scratch/pfxA/sif" -maxdepth 1 -name '*.sif' -links 1 | wc -l)" -eq 1 ]
[ "$(find "$scratch/pfxA/versions/$(cat VERSION)/sif" -maxdepth 1 -name '*.sif' | wc -l)" -eq "$n_fake" ]
cmp -s "$changed" "$scratch/pfxA/versions/$(cat VERSION)/sif/$(basename "$changed")"
# a delta on a site that never installed the full bundle is refused before anything is placed
if scripts/install-bundle.sh --bundle "$tarB" --prefix "$scratch/pfxMissing" > "$scratch/installB-missing.log" 2>&1; then
  echo "installer accepted a delta bundle on a site without its images" >&2; exit 1
fi
grep -q 'Install the full bundle' "$scratch/installB-missing.log"
[ ! -d "$scratch/pfxMissing/versions" ] && [ ! -d "$scratch/pfxMissing/sif" ]
# the old layout (no store, every version its own copy): the first delta install adopts from versions/<v>/sif,
# and --dedupe-sif later replaces copies by hardlinks
scripts/install-bundle.sh --bundle "$tarA" --prefix "$scratch/pfxOld" > /dev/null 2>&1
rm -rf "$scratch/pfxOld/sif"
for f in "$scratch/pfxOld/versions/$(cat VERSION)/sif"/*.sif; do cp -p "$f" "$f.tmp" && mv -f "$f.tmp" "$f"; done
[ "$(stat -c %h "$changed")" -ge 1 ]
scripts/install-bundle.sh --bundle "$tarB" --prefix "$scratch/pfxOld" --force > "$scratch/installB-old.log" 2>&1
grep -q "$((n_fake - 1)) adopted from the installed versions" "$scratch/installB-old.log"
[ "$(find "$scratch/pfxOld/sif" -maxdepth 1 -name '*.sif' | wc -l)" -eq "$n_fake" ]
for f in "$scratch/pfxOld/versions/$(cat VERSION)/sif"/*.sif; do cp -p "$f" "$f.tmp" && mv -f "$f.tmp" "$f"; done
scripts/install-bundle.sh --prefix "$scratch/pfxOld" --dedupe-sif > "$scratch/dedupe.log" 2>&1
grep -q "dedupe: $n_fake file(s) replaced by store hardlinks" "$scratch/dedupe.log"
for f in "$scratch/pfxOld/versions/$(cat VERSION)/sif"/*.sif; do [ "$(stat -c %h "$f")" -ge 2 ]; done
echo "== bundle test passed ($(du -h "$tarball" | cut -f1) bundle; delta $(du -h "$tarB" | cut -f1))"
