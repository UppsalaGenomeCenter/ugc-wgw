#!/usr/bin/env bash
# Smoke test runner: drives every stage through bin/ugc-wgw on the chr20 dataset
# prepared by tests/smoke/prepare.sh, then checks tests/smoke/expected/<stage>.txt.
# See tests/smoke/README.md and docs/DESIGN.md §15.
#
# Usage:
#   run.sh init [--install DIR | --code DIR --miniwdl EXE [--cfg FILE]] [--sif-cache DIR]
#               [--inflight N] [--deepvariant cpu|gpu|parabricks] [--gpu-type TYPE]
#               [--parabricks-gpus N] [--force]  create the ugc-wgw project on the smoke data; a GPU
#                                                flavour adds --nv to the local miniwdl.cfg and, without
#                                                SLURM, runs one sample at a time (--inflight overrides)
#   run.sh <stage> [<subject>] [--mode M] [--retry [--max-attempts N]]   submit (or retry) one stage, then check
#   run.sh check <stage> [<subject>]           check expected outputs of the current attempt
#   run.sh all [standalone] [joint] [assembly]  the phases in that order (default: all three):
#                                              standalone (SMOKE) -> joint (SMOKEJ) -> assembly HG002
#   run.sh status                              ugc-wgw status for the three modes
#   run.sh report                              ugc-wgw report (HTML) and the analysis summaries into $UGC_WGW_SMOKE_DIR/logs/
#   run.sh summary                             ugc-wgw summary only (per-sample and cohort analysis pages)
#   run.sh sacct                               SLURM accounting of every job the runs submitted (HPC only)
#
# Environment:
#   UGC_WGW_SMOKE_DIR   root with data/ and references/ from prepare.sh (default: $HOME/ugc-smoke);
#                   references/ is the data container's tree and manifest.json plus ugc-wgw-extras;
#                   run.sh adds project/, results/, cache/, maps/, logs/ and miniwdl.cfg
# Stages default to subject HG002 (sample stages) or the cohort of their mode (cohort stages);
# singleton/cohort_merge run in standalone, upstream/cohort_call/downstream in joint,
# assembly in assembly. Joint mode uses cohort SMOKEJ (same members) because results paths
# do not encode the mode and cohort_merge would otherwise count as done after standalone.
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
smoke=${UGC_WGW_SMOKE_DIR:-$HOME/ugc-smoke}
project=$smoke/project
results=$smoke/results
ugc=()

log() { echo "[smoke] $*" >&2; }
die() { echo "[smoke] error: $*" >&2; exit 1; }
usage() { sed -n '2,22p' "$0"; exit "${1:-0}"; }

mode_of_stage() {
  case "$1" in
    singleton|cohort_merge|cohort_freq) echo standalone ;;
    upstream|cohort_call|downstream) echo joint ;;
    assembly) echo assembly ;;
    *) die "unknown stage $1" ;;
  esac
}
cohort_of_mode() { case "$1" in joint) echo SMOKEJ ;; *) echo SMOKE ;; esac; }
is_cohort_stage() { case "$1" in cohort_call|cohort_merge|cohort_freq) return 0 ;; *) return 1 ;; esac; }

members() { sed 's/#.*//' "$smoke/data/cohort.txt" | awk 'NF { print $1 }'; }

set_ugc() {
  [ -f "$project/.ugc-wgw/config.json" ] || die "no smoke project at $project; run: $0 init"
  ugc=(python3 "$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["code_dir"])' "$project/.ugc-wgw/config.json")/bin/ugc-wgw" --project "$project")
}

# ---- init --------------------------------------------------------------------
cmd_init() {
  local install="" code="$repo" miniwdl="$repo/.venv/bin/miniwdl" cfg="" sif="$repo/bundle/out/sif-cache" inflight="" force=0
  local deepvariant=cpu gpu_type="" parabricks_gpus=""
  while [ $# -gt 0 ]; do
    case "$1" in
      --install) install="$2"; shift 2 ;;
      --code) code="$2"; shift 2 ;;
      --miniwdl) miniwdl="$2"; shift 2 ;;
      --cfg) cfg="$2"; shift 2 ;;
      --sif-cache) sif="$2"; shift 2 ;;
      --inflight) inflight="$2"; shift 2 ;;
      --deepvariant) deepvariant="$2"; shift 2 ;;
      --gpu-type) gpu_type="$2"; shift 2 ;;
      --parabricks-gpus) parabricks_gpus="$2"; shift 2 ;;
      --force) force=1; shift ;;
      -h|--help) usage ;;
      *) die "init: unknown option $1" ;;
    esac
  done
  [ -f "$smoke/data/samples.tsv" ] && [ -f "$smoke/data/cohort.txt" ] || die "no smoke data under $smoke/data; run tests/smoke/prepare.sh"
  # default runs in flight: two, except a GPU flavour without SLURM, where every task sees every GPU and two
  # call_variants processes on one card end in CUDA_ERROR_OUT_OF_MEMORY (SLURM confines devices per job)
  if [ -z "$inflight" ]; then
    if [ -z "$install" ] && [ "$deepvariant" != cpu ]; then inflight=1; else inflight=2; fi
  fi
  local ref_map build subdir
  if [ -n "$install" ]; then
    install=$(cd "$install" && pwd)
    code=$install/code; miniwdl=$install/venv/bin/miniwdl; cfg=${cfg:-$install/miniwdl.cfg}
    build=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["bundles"]["hifi-wdl-resources"]["build"])' "$code/references.lock")
    ref_map=$install/references/ugc_wgw_ref_map.$build.tsv
    [ -f "$ref_map" ] || die "no rendered reference map in the install: $ref_map (install-bundle.sh step 4)"
  else
    read -r build subdir < <(python3 -c 'import json,sys; b=json.load(open(sys.argv[1]))["bundles"]["hifi-wdl-resources"]; print(b["build"], b["install_subdir"])' "$code/references.lock")
    [ -f "$smoke/references/$subdir/manifest.json" ] || die "reference tree missing under $smoke/references/$subdir; run tests/smoke/prepare.sh"
    [ -d "$sif" ] || die "SIF cache not found: $sif (--sif-cache)"
    mkdir -p "$smoke/maps" "$smoke/cache"
    # read_map() rejects comment lines: drop them like install-bundle.sh does
    sed -e '/^#/d' -e '/^$/d' -e "s|<local_path_prefix>|$smoke/references|g" "$code/references/ugc_wgw_ref_map.$build.template.tsv" > "$smoke/maps/ugc_wgw_ref_map.$build.tsv"
    ref_map=$smoke/maps/ugc_wgw_ref_map.$build.tsv
    if [ -z "$cfg" ]; then
      cfg=$smoke/miniwdl.cfg
      [ -e "$smoke/resources.tsv" ] || cp "$code/backends/hpc/resources.tsv.example" "$smoke/resources.tsv"
      sed -e "s|{{TASK_CONCURRENCY}}|8|" -e "s|{{CALL_CACHE_DIR}}|$smoke/cache|" -e "s|{{SIF_CACHE_DIR}}|$sif|" \
          -e "s|{{RESOURCES_FILE}}|$smoke/resources.tsv|" "$here/miniwdl.local.cfg.template" > "$cfg"
      grep -q '{{' "$cfg" && die "unrendered placeholder in $cfg"
      # a GPU flavour needs the host driver in the container (the installer does this from SINGULARITY_NV=1)
      [ "$deepvariant" = cpu ] || sed -i 's/"hostfs"\]/"hostfs", "--nv"]/' "$cfg"
    fi
    local missing=0 p
    while IFS=$'\t' read -r _ p; do
      case "$p" in /*) [ -e "$p" ] || { log "missing reference: $p"; missing=$((missing + 1)); } ;; esac
    done < "$ref_map"
    [ "$missing" -eq 0 ] || die "$missing reference path(s) missing; check the reference tree under $smoke/references"
  fi
  [ -x "$miniwdl" ] || die "miniwdl not executable: $miniwdl"
  "$(dirname "$miniwdl")/python" -c 'import ugc_wgw_miniwdl.resources' 2>/dev/null \
    || log "warning: the ugc_wgw_resources task plugin is not installed in $(dirname "$miniwdl"); the resource policy will not apply (pip install -e plugins/ugc_wgw_miniwdl)"
  if [ -d "$project/.ugc-wgw" ]; then
    [ "$force" -eq 1 ] || die "project exists: $project (use --force to recreate .ugc-wgw; results are kept)"
    rm -rf "$project/.ugc-wgw"
  fi
  mkdir -p "$project" "$results" "$smoke/logs"
  local -a gpu_args=(--deepvariant "$deepvariant")
  [ -z "$gpu_type" ] || gpu_args+=(--gpu-type "$gpu_type")
  [ -z "$parabricks_gpus" ] || gpu_args+=(--parabricks-gpus "$parabricks_gpus")
  if [ "$deepvariant" != cpu ] && ! grep -q -- '--nv' "$cfg"; then
    log "warning: $cfg has no --nv in [singularity] run_options; the GPU tasks will not see a device (site.cfg SINGULARITY_NV=1)"
  fi
  "$code/bin/ugc-wgw" init "$project" --code "$code" --miniwdl "$miniwdl" --cfg "$cfg" \
      --ref-map "$ref_map" --results "$results" \
      --max-inflight "$inflight" --poll-interval 15 "${gpu_args[@]}"
  set_ugc
  # the sheet from prepare.sh names the paths of the machine that prepared the data;
  # point every read at this machine's data/ (columns 3 and 4: hifi_reads, fail_reads)
  awk -F'\t' -v OFS='\t' -v d="$smoke/data" 'NR == 1 { print; next }
    { for (i = 3; i <= 4; i++) if ($i != "") { n = split($i, parts, ","); out = ""
        for (j = 1; j <= n; j++) { sub(/.*\//, "", parts[j]); out = out (j > 1 ? "," : "") d "/" parts[j] }
        $i = out }
      print }' "$smoke/data/samples.tsv" > "$project/samples.tsv"
  "${ugc[@]}" samples add "$project/samples.tsv"
  local ids=$smoke/data/cohort.txt
  "${ugc[@]}" cohort freeze SMOKE --samples "$ids"
  "${ugc[@]}" cohort freeze SMOKEJ --samples "$ids"
  log "project ready: $project (code $code, cfg $cfg, results $results)"
}

# ---- check -------------------------------------------------------------------
cmd_check() {
  local stage=$1 subject=${2:-}
  [ -n "$stage" ] || usage 2
  set_ugc
  local mode; mode=$(mode_of_stage "$stage")
  if [ -z "$subject" ]; then
    if is_cohort_stage "$stage"; then subject=$(cohort_of_mode "$mode"); else subject=HG002; fi
  fi
  local kind=samples; is_cohort_stage "$stage" && kind=cohorts
  local version; version=$(cat "$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["code_dir"])' "$project/.ugc-wgw/config.json")/VERSION")
  local cur=$results/$kind/$subject/$version/$stage/current
  [ -d "$cur" ] || { log "check $stage $subject: no finished attempt at $cur"; return 1; }
  local status
  status=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get("status"))' "$cur/run_manifest.json" 2>/dev/null || echo none)
  # the mode the run was made in decides the mode-specific expectations (cohort_merge differs per mode)
  mode=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get("mode") or sys.argv[2])' "$cur/run_manifest.json" "$mode" 2>/dev/null || echo "$mode")
  local -a expected=("$here/expected/$stage.txt")
  [ -f "$here/expected/$stage.$mode.txt" ] && expected+=("$here/expected/$stage.$mode.txt")
  local bad=0 n=0 pattern
  while IFS= read -r pattern; do
    case "$pattern" in ''|'#'*) continue ;; esac
    n=$((n + 1))
    local hit=0 f
    for f in "$cur"/out/$pattern; do
      if [ -s "$f" ]; then hit=1; break; fi
    done
    if [ "$hit" -eq 0 ]; then log "  MISSING $pattern"; bad=$((bad + 1)); fi
  done < <(cat "${expected[@]}")
  if [ "$status" = success ] && [ "$bad" -eq 0 ]; then
    log "check $stage $subject: ok ($n patterns, manifest success, $cur)"
    return 0
  fi
  log "check $stage $subject: FAILED (manifest status $status, $bad of $n patterns missing)"
  return 1
}

# ---- one stage ---------------------------------------------------------------
cmd_stage() {
  local stage=$1; shift
  local subject="" mode="" verb=submit
  local -a extra=()
  while [ $# -gt 0 ]; do
    case "$1" in
      --mode) mode="$2"; shift 2 ;;
      --retry) verb=retry; shift ;;
      --max-attempts) extra+=(--max-attempts "$2"); shift 2 ;;
      -*) die "unknown option $1" ;;
      *) subject="$1"; shift ;;
    esac
  done
  set_ugc
  mode=${mode:-$(mode_of_stage "$stage")}
  local cohort; cohort=$(cohort_of_mode "$mode")
  local -a sel
  if is_cohort_stage "$stage"; then
    subject=${subject:-$cohort}; sel=(--cohort "$subject")
  else
    subject=${subject:-HG002}; sel=(--samples "$subject")
    [ "$mode" = joint ] && sel+=(--cohort "$cohort")
  fi
  mkdir -p "$smoke/logs"
  local logf stamp
  stamp=$(date -u +%Y%m%dT%H%M%SZ)
  logf=$smoke/logs/$stage-$subject-$stamp.log
  log "$verb --mode $mode --stage $stage ${sel[*]} (log $logf)"
  local rc=0
  "${ugc[@]}" -v "$verb" --mode "$mode" --stage "$stage" "${sel[@]}" "${extra[@]}" 2>&1 | tee "$logf" || rc=${PIPESTATUS[0]}
  [ "$rc" -eq 0 ] || log "ugc-wgw $verb exited $rc"
  cmd_check "$stage" "$subject"
}

# ---- all ---------------------------------------------------------------------
phase() {  # phase <name> <ugc-wgw submit args...>
  local name=$1; shift
  local t0; t0=$(date +%s)
  mkdir -p "$smoke/logs"
  local logf stamp
  stamp=$(date -u +%Y%m%dT%H%M%SZ)
  logf=$smoke/logs/all-$name-$stamp.log
  log "phase $name: ugc-wgw submit $* (log $logf)"
  local rc=0
  "${ugc[@]}" -v submit "$@" 2>&1 | tee "$logf" || rc=${PIPESTATUS[0]}
  log "phase $name finished in $(( $(date +%s) - t0 )) s, exit $rc"
}

cmd_all() {  # cmd_all [standalone] [joint] [assembly]: the named phases in the fixed order (default all)
  set_ugc
  local t0 fails=0 m p want=" $* "
  [ $# -gt 0 ] || want=" standalone joint assembly "
  for p in "$@"; do case "$p" in standalone|joint|assembly) ;; *) die "all: unknown phase $p (standalone, joint, assembly)" ;; esac; done
  t0=$(date +%s)
  if [[ "$want" == *" standalone "* ]]; then
    phase standalone --mode standalone --cohort SMOKE
    for m in $(members); do cmd_check singleton "$m" || fails=$((fails + 1)); done
    cmd_check cohort_merge SMOKE || fails=$((fails + 1))
    cmd_check cohort_freq SMOKE || fails=$((fails + 1))
  fi
  if [[ "$want" == *" joint "* ]]; then
    phase joint --mode joint --cohort SMOKEJ
    for m in $(members); do cmd_check upstream "$m" || fails=$((fails + 1)); done
    cmd_check cohort_call SMOKEJ || fails=$((fails + 1))
    for m in $(members); do cmd_check downstream "$m" || fails=$((fails + 1)); done
    cmd_check cohort_merge SMOKEJ || fails=$((fails + 1))
    cmd_check cohort_freq SMOKEJ || fails=$((fails + 1))
  fi
  if [[ "$want" == *" assembly "* ]]; then
    phase assembly --mode assembly --samples HG002
    cmd_check assembly HG002 || fails=$((fails + 1))
  fi
  log "phases [${want# }] done in $(( $(date +%s) - t0 )) s; $fails check(s) failed"
  cmd_report || true
  [ "$fails" -eq 0 ]
}

cmd_status() {
  set_ugc
  "${ugc[@]}" status --mode standalone --cohort SMOKE
  "${ugc[@]}" status --mode joint --cohort SMOKEJ
  "${ugc[@]}" status --mode assembly
}

cmd_sacct() {  # SLURM accounting for every job id the manifests recorded; requested vs used memory per task
  set_ugc
  command -v sacct >/dev/null || die "sacct not found: this needs SLURM"
  mkdir -p "$smoke/logs"
  local out
  out=$smoke/logs/sacct-$(date -u +%Y%m%dT%H%M%SZ).tsv
  local ids
  ids=$(find "$results" -name run_manifest.json -exec python3 -c '
import json, sys
for p in sys.argv[1:]:
    try: print(*json.load(open(p)).get("slurm_job_ids", []))
    except (OSError, ValueError): pass' {} + | tr " " "\n" | grep -E "^[0-9]+$" | sort -un | paste -sd, -)
  [ -n "$ids" ] || die "no slurm_job_ids in any manifest under $results"
  sacct -j "$ids" -P --units=M --format=JobID,JobName%60,State,ExitCode,Elapsed,Timelimit,AllocCPUS,ReqMem,MaxRSS,NodeList,AllocTRES > "$out"
  log "wrote $out"
  python3 - "$out" <<'PY'
import collections, re, sys
rows = [l.rstrip("\n").split("|") for l in open(sys.argv[1])]
head, rows = rows[0], rows[1:]
ix = {h: i for i, h in enumerate(head)}
def mb(v):
    m = re.match(r"([0-9.]+)([KMGT]?)", v or "")
    if not m:
        return None
    n, u = float(m.group(1)), m.group(2)
    return n * {"K": 1 / 1024, "": 1, "M": 1, "G": 1024, "T": 1024 * 1024}[u]
jobs = {}
for r in rows:
    jid = r[ix["JobID"]]
    base = jid.split(".")[0]
    j = jobs.setdefault(base, {"name": "", "state": "", "elapsed": "", "limit": "", "cpus": "", "req": None, "max": 0.0, "gres": "-"})
    if "." not in jid:
        tres = r[ix["AllocTRES"]] if "AllocTRES" in ix else ""
        gres = ",".join(t.split("gres/", 1)[1] for t in tres.split(",") if t.startswith("gres/gpu"))
        j.update(name=r[ix["JobName"]], state=r[ix["State"]], elapsed=r[ix["Elapsed"]], limit=r[ix["Timelimit"]],
                 cpus=r[ix["AllocCPUS"]], req=mb(r[ix["ReqMem"]]), gres=gres or "-")
    m = mb(r[ix["MaxRSS"]])
    if m:
        j["max"] = max(j["max"], m)
by = collections.defaultdict(list)
for j in jobs.values():
    by[re.sub(r"-\d+(-[^-]+)?$", "", j["name"])].append(j)  # call-x-03-chr3 -> call-x
print(f"{'task':40} {'jobs':>4} {'state':10} {'limit':>11} {'cpus':>4} {'gres':>12} {'max elapsed':>11} {'req MB':>8} {'max RSS MB':>10} {'use':>5}")
def secs(t):
    d, _, hms = t.rpartition("-")
    parts = [int(x) for x in hms.split(":")] if hms.replace(":", "").isdigit() else [0]
    while len(parts) < 3:
        parts.insert(0, 0)
    return (int(d) if d else 0) * 86400 + parts[0] * 3600 + parts[1] * 60 + parts[2]
for name, js in sorted(by.items(), key=lambda kv: -max(j["max"] for j in kv[1])):
    req = max((j["req"] or 0) for j in js)
    mx = max(j["max"] for j in js)
    states = ",".join(sorted({j["state"].split()[0] for j in js}))
    longest = max(js, key=lambda j: secs(j["elapsed"]))["elapsed"]
    gres = ",".join(sorted({j["gres"] for j in js}))
    print(f"{name[:40]:40} {len(js):4d} {states[:10]:10} {js[0]['limit']:>11} {js[0]['cpus']:>4} {gres[:12]:>12} {longest:>11} {req:8.0f} {mx:10.0f} {(100 * mx / req if req else 0):4.0f}%")
PY
}

cmd_report() {  # self-contained HTML summary of every run of the smoke project, then the analysis summaries
  set_ugc
  mkdir -p "$smoke/logs"
  local out
  out=$smoke/logs/report-$(date -u +%Y%m%dT%H%M%SZ).html
  "${ugc[@]}" report --any-version --sizes --out "$out"
  cmd_summary
}

cmd_summary() {  # ugc-wgw summary for both cohorts: per-sample pages plus the cohort pages under logs/summary-<stamp>/
  set_ugc
  local dir fails=0 f
  dir=$smoke/logs/summary-$(date -u +%Y%m%dT%H%M%SZ)
  "${ugc[@]}" summary --cohort SMOKE --mode standalone --any-version --out-dir "$dir/standalone" || fails=$((fails + 1))
  "${ugc[@]}" summary --cohort SMOKEJ --mode joint --any-version --out-dir "$dir/joint" || fails=$((fails + 1))
  for f in standalone/SMOKE.summary.html standalone/SMOKE.summary.json standalone/HG002.summary.html standalone/HG002.summary.json \
           joint/SMOKEJ.summary.html joint/HG002.summary.html; do
    [ -s "$dir/$f" ] || { log "summary: missing or empty $dir/$f"; fails=$((fails + 1)); }
  done
  if [ "$fails" -eq 0 ]; then log "summary: ok ($dir)"; else log "summary: FAILED ($fails problem(s), $dir)"; fi
  echo "$dir"
  [ "$fails" -eq 0 ]
}

[ $# -ge 1 ] || usage 2
case "$1" in
  init) shift; cmd_init "$@" ;;
  check) shift; cmd_check "$@" ;;
  all) shift; cmd_all "$@" ;;
  status) cmd_status ;;
  report) cmd_report ;;
  summary) cmd_summary ;;
  sacct) cmd_sacct ;;
  -h|--help) usage ;;
  singleton|upstream|cohort_call|downstream|cohort_merge|cohort_freq|assembly) cmd_stage "$@" ;;
  *) die "unknown command $1" ;;
esac
