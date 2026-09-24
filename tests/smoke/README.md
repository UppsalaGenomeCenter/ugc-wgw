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
| `run.sh` | Create a ugc-wgw project on the dataset (`init`), run one stage (`run.sh <stage> [<subject>]`), check outputs (`check`), or run everything (`all`). |
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
`cohort_freq`) and `--mode assembly
--samples HG002` (trio-binned with HG003 and HG004). It is resumable: the
driver skips what already succeeded. Two cohort IDs are used because results
paths do not encode the mode, so the standalone `cohort_merge` would count as
done in joint mode.

Logs of every submit are under `$UGC_WGW_SMOKE_DIR/logs/`; the project is
`$UGC_WGW_SMOKE_DIR/project`, results under `$UGC_WGW_SMOKE_DIR/results`.

