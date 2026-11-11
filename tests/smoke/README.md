# Smoke test

Runs every stage on a chr20 dataset through the driver and checks that the
expected outputs exist and are non-empty. Not byte-for-byte: DeepVariant and
hifiasm are not bit-reproducible across versions. Run on demand before
tagging and after an upstream sync; not part of CI.

## Files

| File | Purpose |
|---|---|
| `PREPARE.md`, `prepare.sh` | Prepare the dataset once (reference tree out of the data container + chr20 slices of the GIAB Revio trio). |
| `ship.sh`, `HPC.md` | Stage bundle (references inside), dataset and site file for the share; the checklist for the run on the HPC and what to verify there. |
| `run.sh` | Create a ugc-wgw project on the dataset (`init`), run one stage (`run.sh <stage> [<subject>]`, any of the nine: `assembly_call HG002`, `cohort_pangenome SMOKE`, ...), check outputs (`check`), or run everything (`all`). |
| `miniwdl.local.cfg.template` | miniwdl config for a run on one machine without SLURM; rendered by `run.sh init`. |
| `expected/<stage>.txt` | One glob per line relative to the run's `out/` (`name/<file>`, arrays `name/<i>/<file>`). Presence and non-emptiness only. |

## Quick start

```bash
export UGC_WGW_SMOKE_DIR=/data/ugc-smoke
tests/smoke/prepare.sh --threads 8                        # once; see PREPARE.md
tests/smoke/run.sh init --inflight 2                      # dev machine: repo .venv miniwdl, bundle/out/sif-cache
tests/smoke/run.sh init --install /proj/ugc/current       # HPC: the installed engine, config and maps
tests/smoke/run.sh all                                    # standalone, joint, assembly; checks after each
tests/smoke/run.sh singleton HG002                        # one stage; subject defaults to HG002 / the mode's cohort
tests/smoke/run.sh status
```

`all` runs `ugc-wgw submit --mode standalone --cohort SMOKE` (four `singleton`
runs, then `cohort_merge` and `cohort_freq`), `--mode joint --cohort SMOKEJ`
(`upstream` ×4, `cohort_call`, `downstream` ×4, `cohort_merge`,
`cohort_freq`) and `--mode assembly --cohort SMOKE` (`assembly` ×4, HG002
trio-binned with HG003 and HG004, then `assembly_call` ×4 and
`cohort_pangenome` on SMOKE, a chr20 graph with 17 expected outputs: `init`
sets `stage_inputs.cohort_pangenome.chromosomes` to `["chr20"]`; chr20-only
assemblies, so every sample infers karyotype `unknown`). It is resumable:
the driver skips what already succeeded. Two cohort IDs are used because
results paths do not encode the mode, so the standalone `cohort_merge` would
count as done in joint mode.

Logs of every submit are under `$UGC_WGW_SMOKE_DIR/logs/`; the project is
`$UGC_WGW_SMOKE_DIR/project`, results under `$UGC_WGW_SMOKE_DIR/results`.
`run.sh init` renders the reference map once (`maps/`); after a change to the
map template or to the references (a new `references.lock` entry), run
`prepare.sh --references-only` and `run.sh init --force` (the call cache is
kept), or the stages that read the new keys fail at their first task with
"Map key not found".

## Evaluation profiles

`run.sh eval sawfish23` initialises a second project on the same data,
references, call cache and SIF cache (`project-sawfish23/`,
`results-sawfish23/`; `UGC_WGW_SMOKE_VARIANT` names it) with the
`eval-sawfish23` profile, runs the standalone phase there (the cache
replays everything up to sawfish), and writes to `logs/eval-sawfish23-<stamp>/`
a per-member truvari comparison of the phased SV VCFs against the
production project's (`sv.tsv`: TP, FP, FN, precision, recall, F1,
genotype concordance) and the sawfish tasks' resources of both projects
(`sawfish.tsv`). Paste both tables into the run entry below; `eval
sawfish23 compare` repeats the comparison alone on the finished runs.
`init --profile NAME` applies any shipped profile to the smoke project.

## Validation against GIAB

`run.sh validate [<member>]` (default HG002) scores the member's newest
successful per-sample run of each route (`singleton`, and `downstream` of
the joint mode) against GIAB's HG002 v5.0q benchmark on chr20: the truth
slices are the `ugc-wgw-truth` entry of `references.lock`
(`prepare.sh --references-only` installs them), aardvark and truvari come
from the SIF cache (both are in `image_manifest.ugc-wgw.txt`; neither is
used by a workflow). It writes `logs/validate-<stamp>/` (`small.tsv`,
`sv.tsv`, the tools' directories with the TP/FP/FN VCFs) and prints the
rows for `VALIDATION.md`, which is the record: add them there with the
release. About a minute. `--min-f1 snv=0.995` (also `indel`, `sv`) makes a
figure a gate; `--truth-dir` and `--region` point it at GIAB's full files
for a whole-genome HG002.

## Entry kinds

`run.sh entry` (0.16.0) needs the main project's joint phase: it writes a
sheet under `logs/entry-<stamp>/` registering `HG002b` as the `bam` kind
(HG002's aligned BAM from its `upstream` run) and `HG003v` as the `vcf`
kind (HG003's aligned BAM, DeepVariant VCF, sawfish VCF and gVCF, with
the `_index` columns because the run keeps indices in their own output
directories). The two BAMs are reheadered copies (`samtools reheader` from
the pb_wdl_base SIF, `SM:HG002` → `SM:HG002b`, `SM:HG003` → `SM:HG003v`,
then indexed): the read group's sample name must be the sample id. It
initialises `project-entry/` with the sheet (`UGC_WGW_SMOKE_SHEET`
and `UGC_WGW_SMOKE_COHORT` point `init` at another sheet and cohort file),
runs the standalone phase (HG002b: `upstream` + `downstream`; HG003v:
`intake` + `downstream`; the four cohort stages), compares HG003v's phased
VCFs with HG003's singleton ones (`compare.tsv`: record counts and whether
the bodies are identical), then the negative cases: `HG003r` (HG003v's BAM
under another id) must be refused by `samples add`, and `HG003x` (its own
reheadered BAM, HG003's small-variant VCF with the contigs renamed `chr20`
→ `20`) must fail `intake` as an input error. The checks of `all standalone` are kind-aware: a member is
checked on the stages it ran (`sample_stages_of`); `check cohort_merge` and
`check cohort_freq` do not expect the mitochondrial outputs when the run's
inputs show the branch off (a `vcf`-kind member has no mitorsaw call), and
`summary` expects the pages of the first member with a run and of the modes
whose phase ran in the results.

