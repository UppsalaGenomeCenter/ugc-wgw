# ugc-wgw

<p align="center">
    <img src="docs/images/ugc_wgw_driver.jpeg" alt="ugc-wgw driver" width="480">
</p>

[![CI](https://github.com/UppsalaGenomeCenter/ugc-wgw/actions/workflows/ci.yml/badge.svg)](https://github.com/UppsalaGenomeCenter/ugc-wgw/actions/workflows/ci.yml)

PacBio HiFi whole-genome analysis for cohorts of hundreds to thousands of
samples on an HPC without network access. ugc-wgw composes PacBio's
[HiFi-human-WGS-WDL](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL)
into stages that a small driver sequences across the cohort, and packages
everything the cluster needs, workflows, containers, engine and reference
data, into one verifiable offline bundle.

## What it adds to the upstream workflow

- **Stages instead of one run per family.** Per-sample alignment, variant
  calling, phasing and repeat genotyping run as independent, restartable
  miniwdl runs; cohort stages merge structural variants (svx), tandem
  repeats (`trgt merge`, trgt-lps) and small variants (GLnexus), and compute
  allele frequencies. Two modes: standalone (per-sample DeepVariant, cohort
  merges afterwards) and joint (GLnexus genotypes before phasing, as
  upstream's family workflow does). De novo assembly (hifiasm, trio-binned
  when parents are present) is a third mode.
- **A driver, `ugc-wgw`,** that registers samples, freezes cohorts, generates
  every input file, submits runs at a chosen concurrency, retries transient
  failures, records provenance in a manifest per run, and reports
  progress, failures and resource usage. Stdlib Python, SQLite state, one
  SLURM job per task through miniwdl-slurm.
- **Offline bundles.** `scripts/make-bundle.sh` builds a tar with the code,
  every container as a SIF, the engine's wheels and PacBio's reference data
  container, each pinned by digest or checksum; `install-bundle.sh` verifies
  and installs it on the cluster, renders the site's configuration and
  never touches the network.
- **Site resource policy.** Per-task cores, memory, wall time and partition
  are rewritten by a miniwdl plugin from one table, so no WDL is edited to
  fit a cluster. GPU DeepVariant and NVIDIA Parabricks are one project
  switch.
- **Reports.** A self-contained HTML run report (timings, cache hits,
  failures, provenance) and per-sample and per-cohort analysis summaries
  (coverage, variants, phasing, repeats, methylation, QC flags), drawn with
  inline SVG, no viewer needed.

## Quick start

On the HPC, with a bundle on the share (`docs/guide/00-quick-start.md` has
the full walkthrough):

```bash
./install-bundle.sh --bundle ugc-pacbio-wgw-X.Y.Z.tar --prefix /proj/ugc --verify-only
./install-bundle.sh --bundle ugc-pacbio-wgw-X.Y.Z.tar --prefix /proj/ugc \
    --references /proj/ugc/references --site site.cfg --activate
export PATH=/proj/ugc/current/code/bin:$PATH

ugc-wgw init /proj/ugc/projects/cohort2026 --install /proj/ugc/current \
    --ref-map /proj/ugc/current/references/ugc_wgw_ref_map.GRCh38_GIABv3.tsv
cd /proj/ugc/projects/cohort2026
ugc-wgw samples add samples.tsv
ugc-wgw cohort freeze C1 --samples cohort.txt
ugc-wgw -v submit --mode standalone --cohort C1 --max-inflight 40
ugc-wgw summary --cohort C1
```

Building the bundle happens on a machine with network access, Docker or
Podman and Apptainer:

```bash
python3 -m venv .venv && .venv/bin/pip install -r bundle/requirements.txt
.venv/bin/pip install -e plugins/ugc_wgw_miniwdl
scripts/make-bundle.sh --version vX.Y.Z --out bundle/out --python-version 3.12
```

The three images this repository builds itself (svx, trgt-lps, hifiasm) are
published on GHCR under `ghcr.io/uppsalagenomecenter/`; `image_manifest.ugc-wgw.txt` pins
them by digest, and `containers/` holds their Dockerfiles.

## Documentation

The operator guide is [`docs/guide/`](docs/guide/README.md): a quick start,
the concepts, choosing a mode, installing, project setup, running, results,
the stages, troubleshooting, upgrading, a command reference and the task
resource model. `bin/README.md` is the driver's reference page. The design
notes and the per-workflow derivation notes of the internal repository are
not part of this copy; ask if you need them.

## Layout

```
workflows/        the seven entrypoints (ugc_wgw_*.wdl), our tasks (ugc_wgw/), one override
vendor/           PacBio's HiFi-human-WGS-WDL at the pinned tag (upstream.lock)
bin/              the ugc-wgw driver (stdlib Python)
plugins/          the miniwdl task plugin that applies the site resource policy
scripts/          bundle build and install, SIF cache, image digests, manifest conversion
backends/hpc/     miniwdl.cfg and inputs templates rendered by the installer
references/       reference map template and scatter regions (references.lock)
containers/       Dockerfiles of the images we build
tests/            unit tests, static checks, bundle round trip, smoke harness
docs/guide/       the operator guide
```

## Development

`tests/check.sh` runs `miniwdl check` on every entrypoint and the drift
guards of the generated docs; `python3 -m unittest discover -s tests/driver
-t .` (and `tests/plugin`, `tests/docs`) are the unit suites;
`tests/bundle.sh` round-trips a bundle; `tests/smoke/` runs every stage on a
chr20 GIAB trio. Lockfiles (`upstream.lock`, `references.lock`,
`image_manifest.ugc-wgw.txt`, `DERIVED_FILES.md`) are the truth about what
is pinned; images are referenced by digest only.

## Licence and acknowledgements

ugc-wgw is released under the MIT licence (`LICENSE`). The vendored
workflows are PacBio's HiFi-human-WGS-WDL, BSD-3-Clause-Clear
(`vendor/hifi-human-wgs-wdl/LICENSE`); the tools it runs (pbmm2,
DeepVariant, HiPhase, sawfish, TRGT, Paraphase, mitorsaw, MethBat, StarPhase,
kivvi, pbjam, GLnexus, hifiasm, svx, trgt-lps and others) are cited in the
analysis summaries and in `vendor/hifi-human-wgs-wdl/docs/`.

Logo fixed by **[Gemini](https://gemini.google.com/)** <br>
Code fixed by **[Claude](https://claude.ai/)** <br>
Ideas and guidance by **[Iggy](https://github.com/iggyB)** <br>
