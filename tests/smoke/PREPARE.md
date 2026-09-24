# Preparing the smoke dataset

Done once, on a machine with internet (the dev machine), by
`tests/smoke/prepare.sh`. The result is copied to the HPC share like any
other input data. About 25 GB of disk and, depending on the connection, one
to two hours.

## What the dataset is

Four samples, chromosome 20 only, HiFi Revio reads with 5mC tags:

| Sample | Source | Role |
|---|---|---|
| `HG002` | GIAB Ashkenazi son, 48x | Every stage; child of the trio for assembly |
| `HG003` | GIAB father, 46x | Cohort member; father for trio binning |
| `HG004` | GIAB mother, 36x | Cohort member; mother for trio binning |
| `HG002ds` | `HG002` downsampled to 50 % (`samtools view -s 42.5`, then reheadered without the `@CO` line samtools adds: HiPhase panics on it, guide chapter 09) | Fourth cohort member |

The reads come from NIST's GIAB release `PacBio_HiFi-Revio_20231031`
(`https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/`):
the BAMs aligned to GRCh38-GIABv3 are indexed, so `samtools view <url> chr20`
fetches only the chr20 blocks over HTTPS (about 1.5, 1.4 and 1.0 GB) instead
of 70 GB per sample. The slices are still aligned BAMs. That is fine: the
pipeline's pbmm2 task (v4.0.0, pbmm2 26.2) detects alignment, strips the
alignments and the `HP`/`PS`/`PC` haplotype tags (`--strip-tags`), disables
chunking for that input and realigns, writing `Input … is already aligned`
to the messages file. `MM`/`ML` survive, so methylation is intact and the
reads go through the same path as production uBAMs.

The references are the `GRCh38_GIABv3` tree of PacBio's reference data
container (`quay.io/pacbio/workflow-data-container-hifi-human-wgs-wdl-grch38_giabv3`,
digest and md5s in `references.lock`; the SIF the bundle ships), copied out
with `apptainer exec` exactly as `install-bundle.sh` does on the HPC
(`docs/guide/04-installing.md`), plus this repository's
`scatter_regions.GRCh38_GIABv3.tsv` as `ugc-wgw-extras-0.2.0/`. Nothing is
downloaded from Zenodo.

## Running it

```bash
export UGC_WGW_SMOKE_DIR=/data/ugc-smoke          # default: $HOME/ugc-smoke
tests/smoke/prepare.sh --threads 8            # references and data
tests/smoke/prepare.sh --data-only            # or one half
```

`prepare.sh` is idempotent: the reference tree is skipped when its marker
and `manifest.json` exist, each slice is skipped when its md5 matches
`data/data.lock`, `--force` redoes everything. The data container SIF is
taken from `$UGC_WGW_SIF_CACHE` (default `bundle/out/sif-cache`) and pulled
there when absent. It leaves:

```
$UGC_WGW_SMOKE_DIR/
├── downloads/                 references.log
├── references/
│   ├── hifi-wdl-resources-v4.0.0-GRCh38_GIABv3/   GRCh38_GIABv3/ (fasta, trgt/, sawfish/, methbat/) + manifest.json
│   └── ugc-wgw-extras-0.2.0/scatter_regions.GRCh38_GIABv3.tsv
└── data/
    ├── HG002.chr20.hifi_reads.bam   HG003..., HG004..., HG002ds...
    ├── samples.tsv                  driver sample sheet with the pedigree
    ├── cohort.txt                   HG002 HG003 HG004 HG002ds
    └── data.lock                    per file: source URL, region, bytes, md5, samtools version
```

Commit `data.lock` values into `tests/smoke/README.md` when the dataset is
regenerated, so a run on another machine can verify it got the same files.

## Copying to the HPC

The HPC has no network. `tests/smoke/ship.sh --out <dir>` stages the bundle
(which carries the reference data container), `data/` and a site file with
checksums in one directory for the share; `tests/smoke/HPC.md` is the
checklist from there on.

## What a local run does and does not verify

A run on the dev machine (`run.sh init` without `--install`) uses miniwdl's
plain Singularity backend: every task runs on the machine, cpu requests
above the core count are rounded down with a warning, memory requests are
reservations only. It verifies the workflow wiring, the driver, the images,
the reference maps, the expected outputs and the manifests. It does not
verify anything SLURM-specific: `time_minutes` reaching `sbatch`, memory
kills, `slurm_job_ids` in manifests, orphan cancellation. Those need the
HPC run through `--install`.

The chr20 subset also means: pbmm2 runs as one chunk per sample; DeepVariant
still schedules its genome-wide shards, most of which finish at once;
Paraphase, mitorsaw, kivvi and StarPhase see no reads in their regions and
must produce empty results rather than fail; hifiasm assembles chr20 only.
