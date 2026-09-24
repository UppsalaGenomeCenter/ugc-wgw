# Changelog

Release notes of the public copy. Each version corresponds to a verified
offline bundle; the upstream tag it is built on is named first.

## 0.4.0

First public release. Built on PacBio HiFi-human-WGS-WDL v4.0.0
(`15e82cb9`) and its reference data container `GRCh38_GIABv3`.

- Seven entrypoints: `singleton`, `upstream`, `cohort_call`, `downstream`,
  `cohort_merge`, `cohort_freq`, `assembly`; standalone, joint and assembly
  modes; cohort stages scattered over 26 region shards at every cohort size.
- The `ugc-wgw` driver: sample sheet and frozen cohorts, plan and submit with
  a concurrency cap, automatic retries of transient failures, orphan
  cancellation, SQLite state, per-run provenance manifests, `status`,
  `progress`, `logs`, `report` and `summary`.
- Offline bundle build and install with digest-pinned containers, engine
  wheels and the reference data container; site caps and a per-task resource
  policy applied by a miniwdl plugin; GPU DeepVariant and Parabricks as a
  project switch.
- Per-sample and per-cohort HTML analysis summaries with QC flags; a
  self-contained run report.
- `scripts/manifest-to-samples.py` turns a per-file delivery manifest into
  the sample sheet.
