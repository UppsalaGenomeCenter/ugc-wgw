#!/usr/bin/env bash
# Smoke test runner: drives every stage through bin/ugc-wgw on the chr20 dataset
# prepared by tests/smoke/prepare.sh, then checks tests/smoke/expected/<stage>.txt.
# See tests/smoke/README.md and docs/DESIGN.md §15.
#
# Usage:
#   run.sh init [--install DIR | --code DIR --miniwdl EXE [--cfg FILE]] [--sif-cache DIR]
#               [--inflight N] [--deepvariant cpu|gpu|parabricks] [--gpu-type TYPE]
#               [--parabricks-gpus N] [--hifiasm-threads N] [--force]
#                                                create the ugc-wgw project on the smoke data; a GPU
#                                                flavour adds --nv to the local miniwdl.cfg and, without
#                                                SLURM, runs one sample at a time (--inflight overrides);
#                                                --hifiasm-threads sets stage_inputs.assembly.hifiasm_threads
#                                                (the fat partition allocates whole 96-core nodes)
#   run.sh <stage> [<subject>] [--mode M] [--retry [--max-attempts N]] [--any-version]   submit (or retry) one stage, then check
#   run.sh check <stage> [<subject>]           check expected outputs of the current attempt
#   run.sh all [standalone] [joint] [assembly]  the phases in that order (default: all three):
#                                              standalone (SMOKE) -> joint (SMOKEJ) -> assembly (all four
#                                              samples, then the SMOKE cohort's chr20 pangenome)
#   run.sh status                              ugc-wgw status for the three modes
#   run.sh report                              ugc-wgw report (HTML) and the analysis summaries into $UGC_WGW_SMOKE_DIR/logs/
#   run.sh summary                             ugc-wgw summary only (per-sample and cohort analysis pages)
#   run.sh sacct                               SLURM accounting of every job the runs submitted (HPC only)
#   run.sh entry                               the entry kinds (0.16.0): a project-entry/ whose sheet registers HG002b (bam kind:
#                                              HG002's aligned BAM from the main project's upstream run) and HG003v (vcf kind:
#                                              HG003's aligned BAM, DeepVariant VCF, sawfish VCF and gVCF), runs its standalone
#                                              phase, compares HG003v's phased VCFs with HG003's singleton ones, and registers
#                                              HG003x with contigs renamed to see intake fail as an input error; needs the joint phase
#   run.sh validate [<member>] [--truth-dir DIR] [--region chr20] [--min-f1 snv|indel|sv=VALUE ...]
#                                              score the member's (default HG002) newest per-sample run against GIAB's HG002
#                                              v5.0q benchmark on chr20: small variants with aardvark, SVs with truvari bench
#                                              and refine; tables under logs/validate-<stamp>/ and the row for VALIDATION.md;
#                                              reports only, unless --min-f1 sets a floor
#   run.sh eval <name> [compare]               a second project (project-<name>/, results-<name>/) with other tools on the
#                                              same data: sawfish23 (the eval-sawfish23 profile), pblab (pbweave + pbdv on
#                                              a GPU), pblab-cpu, pbdv, pbdv-cpu, pbweave. Its standalone phase and
#                                              validate, then per member against the production project: truvari bench of
#                                              the phased SV VCFs, the PASS small variants and TRGT allele lengths shared
#                                              and not, the per-sample digests side by side (compare_runs.py) and the
#                                              changed tasks' resources (ugc-wgw usage); "compare" repeats the comparison
#                                              alone on the finished runs; init accepts --profile NAME and --aligner too
#
# Environment:
#   UGC_WGW_SMOKE_DIR   root with data/ and references/ from prepare.sh (default: $HOME/ugc-smoke);
#                   references/ is the data container's tree and manifest.json plus ugc-wgw-extras;
#                   run.sh adds project/, results/, cache/, maps/, logs/ and miniwdl.cfg
# Stages default to subject HG002 (sample stages) or the cohort of their mode (cohort stages);
# singleton/cohort_merge run in standalone, upstream/cohort_call/downstream in joint,
# assembly/assembly_call/cohort_pangenome in assembly (cohort SMOKE, chromosomes chr20 only: init
# sets stage_inputs.cohort_pangenome.chromosomes, and cohort_tr.outlier_min_called = 4 for the four members). Joint mode uses cohort SMOKEJ (same members) because results paths
# do not encode the mode and cohort_merge would otherwise count as done after standalone.
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
smoke=${UGC_WGW_SMOKE_DIR:-$HOME/ugc-smoke}
variant=${UGC_WGW_SMOKE_VARIANT:-}   # a second project on the same data, references, cache and SIFs: project-<v>/, results-<v>/ (run.sh eval)
project=$smoke/project${variant:+-$variant}
results=$smoke/results${variant:+-$variant}
ugc=()

log() { echo "[smoke] $*" >&2; }
# the teed driver logs are files first: never coloured, whatever UGC_WGW_COLOR or the terminal say
die() { echo "[smoke] error: $*" >&2; exit 1; }
usage() { sed -n '2,/^# Environment:/p' "$0" | sed '$d'; exit "${1:-0}"; }   # the header's description and Usage

mode_of_stage() {
  case "$1" in
    singleton|cohort_qc|cohort_merge|cohort_freq|cohort_tr|cohort_methyl|intake) echo standalone ;;   # the cohort stages exist in joint too: --mode joint; intake: run.sh entry
    upstream|cohort_call|downstream) echo joint ;;
    assembly|assembly_call|cohort_pangenome) echo assembly ;;
    *) die "unknown stage $1" ;;
  esac
}
cohort_of_mode() { case "$1" in joint) echo SMOKEJ ;; *) echo SMOKE ;; esac; }
is_cohort_stage() { case "$1" in cohort_qc|cohort_call|cohort_merge|cohort_freq|cohort_tr|cohort_methyl|cohort_pangenome) return 0 ;; *) return 1 ;; esac; }

sheet=${UGC_WGW_SMOKE_SHEET:-$smoke/data/samples.tsv}     # run.sh entry registers another sheet in its variant project
cohort_file=${UGC_WGW_SMOKE_COHORT:-$smoke/data/cohort.txt}
members() { sed 's/#.*//' "$cohort_file" | awk 'NF { print $1 }'; }
# the per-sample stages a member ran in standalone mode: singleton (reads), or upstream/intake + downstream (bam/vcf kinds)
sample_stages_of() {  # <member>
  local m=$1 version; version=$(cat "$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["code_dir"])' "$project/.ugc-wgw/config.json")/VERSION")
  if [ -d "$results/samples/$m/$version/singleton/current" ]; then echo singleton; return; fi
  local st; for st in upstream intake downstream; do [ -d "$results/samples/$m/$version/$st/current" ] && echo "$st"; done
}

set_ugc() {
  [ -f "$project/.ugc-wgw/config.json" ] || die "no smoke project at $project; run: $0 init"
  ugc=(python3 "$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["code_dir"])' "$project/.ugc-wgw/config.json")/bin/ugc-wgw" --project "$project")
}
# the SIF of an image in the project's image cache: <manifest file> <image name as the manifest writes it, without the digest>.
# The cache is the one the project's miniwdl.cfg names (config.json's miniwdl_cfg: the smoke's own cfg, or the install's).
sif_of() {
  local manifest=$1 image=$2 cfg dir ref name
  cfg=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["miniwdl_cfg"])' "$project/.ugc-wgw/config.json")
  dir=$(python3 -c 'import configparser,sys; c=configparser.ConfigParser(); c.read(sys.argv[1]); print(c.get("singularity", "image_cache"))' "$cfg")
  ref=$(grep -o "$image@sha256:[0-9a-f]*" "$manifest" | head -1 || true)
  [ -n "$ref" ] || die "no image $image in $manifest"
  name="docker___${ref//\//_}"; name="${name//:/_}.sif"
  [ -f "$dir/$name" ] || die "the SIF of $ref is not in $dir ($name); populate the SIF cache from the image manifests"
  echo "$dir/$name"
}

# the deepvariant flavours whose calling task asks for a GPU (config.GPU_MODES)
gpu_flavour() { case "$1" in gpu|parabricks|pbdv) return 0 ;; *) return 1 ;; esac; }

# ---- init --------------------------------------------------------------------
cmd_init() {
  local install="" code="$repo" miniwdl="$repo/.venv/bin/miniwdl" cfg="" sif="$repo/bundle/out/sif-cache" inflight="" force=0
  local deepvariant=cpu gpu_type="" parabricks_gpus="" hifiasm_threads="" profile="" aligner=""
  while [ $# -gt 0 ]; do
    case "$1" in
      --install) install="$2"; shift 2 ;;
      --code) code="$2"; shift 2 ;;
      --miniwdl) miniwdl="$2"; shift 2 ;;
      --cfg) cfg="$2"; shift 2 ;;
      --sif-cache) sif="$2"; shift 2 ;;
      --inflight) inflight="$2"; shift 2 ;;
      --deepvariant) deepvariant="$2"; shift 2 ;;
      --aligner) aligner="$2"; shift 2 ;;
      --gpu-type) gpu_type="$2"; shift 2 ;;
      --parabricks-gpus) parabricks_gpus="$2"; shift 2 ;;
      --hifiasm-threads) hifiasm_threads="$2"; shift 2 ;;
      --profile) profile="$2"; shift 2 ;;
      --force) force=1; shift ;;
      -h|--help) usage ;;
      *) die "init: unknown option $1" ;;
    esac
  done
  if [ ! -f "$sheet" ] || [ ! -f "$cohort_file" ]; then die "no smoke data under $smoke/data ($sheet, $cohort_file); run tests/smoke/prepare.sh"; fi
  # default runs in flight: two, except a GPU flavour without SLURM, where every task sees every GPU and two
  # call_variants processes on one card end in CUDA_ERROR_OUT_OF_MEMORY (SLURM confines devices per job)
  if [ -z "$inflight" ]; then
    if [ -z "$install" ] && gpu_flavour "$deepvariant"; then inflight=1; else inflight=2; fi
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
      if gpu_flavour "$deepvariant"; then
        sed -i 's/^run_options = \[/run_options = ["--nv", /' "$cfg"
        grep -q '^run_options = \["--nv", ' "$cfg" || die "could not add --nv to [singularity] run_options in $cfg"
      fi
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
  [ -z "$profile" ] || gpu_args+=(--profile "$profile")
  [ -z "$aligner" ] || gpu_args+=(--aligner "$aligner")
  if gpu_flavour "$deepvariant" && ! grep -q -- '--nv' "$cfg"; then
    log "warning: $cfg has no --nv in [singularity] run_options; the GPU tasks will not see a device (site.cfg SINGULARITY_NV=1)"
  fi
  "$code/bin/ugc-wgw" init "$project" --code "$code" --miniwdl "$miniwdl" --cfg "$cfg" \
      --ref-map "$ref_map" --results "$results" \
      --max-inflight "$inflight" --poll-interval 15 "${gpu_args[@]}"
  set_ugc
  # the smoke data is chr20 only: the pangenome of the SMOKE cohort builds that chromosome alone
  python3 - "$project/.ugc-wgw/config.json" <<'PY'
import json, sys
path = sys.argv[1]
doc = json.load(open(path))
doc.setdefault("stage_inputs", {}).setdefault("cohort_pangenome", {})["chromosomes"] = ["chr20"]
# four members: cohort_tr's outlier scan needs its floor of called members lowered to run at all
doc["stage_inputs"].setdefault("cohort_tr", {})["outlier_min_called"] = 4
json.dump(doc, open(path, "w"), indent=2, sort_keys=True)
PY
  log "stage_inputs.cohort_pangenome.chromosomes = [chr20]; stage_inputs.cohort_tr.outlier_min_called = 4"
  if [ -n "$hifiasm_threads" ]; then
    # the request and the command's -t both come from this workflow input; a policy row would change only the request
    python3 - "$project/.ugc-wgw/config.json" "$hifiasm_threads" <<'PY'
import json, sys
path, n = sys.argv[1], int(sys.argv[2])
doc = json.load(open(path))
doc.setdefault("stage_inputs", {}).setdefault("assembly", {})["hifiasm_threads"] = n
json.dump(doc, open(path, "w"), indent=2, sort_keys=True)
PY
    log "stage_inputs.assembly.hifiasm_threads = $hifiasm_threads"
  fi
  # the sheet from prepare.sh names the paths of the machine that prepared the data;
  # point every read at this machine's data/ (columns 3 and 4: hifi_reads, fail_reads)
  if [ "$sheet" = "$smoke/data/samples.tsv" ]; then
    awk -F'\t' -v OFS='\t' -v d="$smoke/data" 'NR == 1 { print; next }
      { for (i = 3; i <= 4; i++) if ($i != "") { n = split($i, parts, ","); out = ""
          for (j = 1; j <= n; j++) { sub(/.*\//, "", parts[j]); out = out (j > 1 ? "," : "") d "/" parts[j] }
          $i = out }
        print }' "$sheet" > "$project/samples.tsv"
  else
    cp "$sheet" "$project/samples.tsv"   # another sheet (run.sh entry): absolute paths already
  fi
  "${ugc[@]}" samples add "$project/samples.tsv"
  local ids=$cohort_file
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
  # 0.15.0 invariant: the cohort small-variant VCF carries no chrM record (mitorsaw is the chrM caller; GLnexus skips the shard)
  local jv
  for jv in "$cur"/out/cohort_small_variant_vcf/*.vcf.gz; do
    [ -s "$jv" ] || continue
    if zcat "$jv" | awk -F'\t' '!/^#/ && $1 == "chrM" { found = 1; exit } END { exit !found }'; then
      log "check $stage $subject: chrM records in $(basename "$jv") (GLnexus must skip chrM)"; return 1
    fi
  done
  # the mitochondrial branch is off when the cohort holds a vcf-kind member (no mitorsaw call): its outputs are not expected then
  local mt_on=1 skipped=0
  case "$stage" in
    cohort_merge|cohort_freq)
      mt_on=$(python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); off=d.get("ugc_wgw_cohort_merge.merge_mitochondria") is False or (sys.argv[2] == "cohort_freq" and not d.get("ugc_wgw_cohort_freq.mt_small_variant_vcf")); print(0 if off else 1)' "$cur/inputs.json" "$stage" 2>/dev/null || echo 1) ;;
  esac
  local bad=0 n=0 pattern
  while IFS= read -r pattern; do
    case "$pattern" in ''|'#'*) continue ;; esac
    if [ "$mt_on" -eq 0 ] && [[ "$pattern" == *mt_* ]]; then skipped=$((skipped + 1)); continue; fi
    n=$((n + 1))
    local hit=0 f
    for f in "$cur"/out/$pattern; do
      if [ -s "$f" ]; then hit=1; break; fi
    done
    if [ "$hit" -eq 0 ]; then log "  MISSING $pattern"; bad=$((bad + 1)); fi
  done < <(cat "${expected[@]}")
  if [ "$status" = success ] && [ "$bad" -eq 0 ]; then
    log "check $stage $subject: ok ($n patterns, manifest success, $cur)${skipped:+$([ "$skipped" -gt 0 ] && echo "; $skipped mitochondrial pattern(s) not expected: the branch is off (a vcf-kind member)")}"
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
      --any-version) extra+=(--any-version); shift ;;   # a cohort stage over per-sample runs of an earlier version (the repo's VERSION moved on)
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
  "${ugc[@]}" --color never -v "$verb" --mode "$mode" --stage "$stage" "${sel[@]}" "${extra[@]}" 2>&1 | tee "$logf" || rc=${PIPESTATUS[0]}
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
  "${ugc[@]}" --color never -v submit "$@" 2>&1 | tee "$logf" || rc=${PIPESTATUS[0]}
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
    local st
    for m in $(members); do
      for st in $(sample_stages_of "$m"); do cmd_check "$st" "$m" || fails=$((fails + 1)); done
    done
    cmd_check cohort_qc SMOKE || fails=$((fails + 1))
    cmd_check cohort_merge SMOKE || fails=$((fails + 1))
    cmd_check cohort_freq SMOKE || fails=$((fails + 1))
    cmd_check cohort_tr SMOKE || fails=$((fails + 1))
    cmd_check cohort_methyl SMOKE || fails=$((fails + 1))
  fi
  if [[ "$want" == *" joint "* ]]; then
    phase joint --mode joint --cohort SMOKEJ
    for m in $(members); do cmd_check upstream "$m" || fails=$((fails + 1)); done
    cmd_check cohort_qc SMOKEJ || fails=$((fails + 1))
    cmd_check cohort_call SMOKEJ || fails=$((fails + 1))
    for m in $(members); do cmd_check downstream "$m" || fails=$((fails + 1)); done
    cmd_check cohort_merge SMOKEJ || fails=$((fails + 1))
    cmd_check cohort_freq SMOKEJ || fails=$((fails + 1))
    cmd_check cohort_tr SMOKEJ || fails=$((fails + 1))
    cmd_check cohort_methyl SMOKEJ || fails=$((fails + 1))
  fi
  if [[ "$want" == *" assembly "* ]]; then
    # every smoke sample is assembled (HG002 trio-binned, the others not): chr20 only, enough to wire
    # assembly -> assembly_call -> cohort_pangenome; every sample infers karyotype unknown from a chr20-only
    # assembly, so the SMOKE cohort's chr20 graph has no sex strata; its SVs are compared with the
    # standalone cohort_merge of SMOKE when that phase ran
    phase assembly --mode assembly --cohort SMOKE
    for m in $(members); do cmd_check assembly "$m" || fails=$((fails + 1)); done
    for m in $(members); do cmd_check assembly_call "$m" || fails=$((fails + 1)); done
    cmd_check cohort_pangenome SMOKE || fails=$((fails + 1))
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
  # the driver's own view: accounting.json per attempt (read now for runs that finished without it), per task
  "${ugc[@]}" --color never usage --collect --by task --any-version
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
  "${ugc[@]}" summary --cohort SMOKE --mode standalone --any-version --out-dir "$dir/standalone" \
      --pass-list "$dir/standalone/identity-pass.txt" --pass-on identity || fails=$((fails + 1))   # HG002ds duplicates HG002: both out
  "${ugc[@]}" summary --cohort SMOKEJ --mode joint --any-version --out-dir "$dir/joint" || fails=$((fails + 1))
  "${ugc[@]}" summary --cohort SMOKE --mode assembly --any-version --out-dir "$dir/assembly" || fails=$((fails + 1))
  # the pages of the first member and of the cohorts; a mode's pages are expected only when that phase ran in these results
  local first="" version d
  version=$(cat "$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["code_dir"])' "$project/.ugc-wgw/config.json")/VERSION")
  for d in "$results"/samples/*/; do   # the first member with a finished run (the variant projects register other members)
    if [ -n "$(ls -d "$d$version"/*/current 2>/dev/null)" ]; then first=$(basename "$d"); break; fi
  done
  [ -n "$first" ] || first=$(members | head -1)
  local -a pages=(standalone/SMOKE.summary.html standalone/SMOKE.summary.json "standalone/$first.summary.html" "standalone/$first.summary.json" standalone/identity-pass.txt)
  if [ -d "$results/cohorts/SMOKEJ/$version/cohort_call/current" ]; then pages+=(joint/SMOKEJ.summary.html "joint/$first.summary.html"); else log "summary: no joint phase in $results; its pages are not expected"; fi
  if [ -d "$results/samples/$first/$version/assembly/current" ]; then pages+=(assembly/SMOKE.summary.html "assembly/$first.summary.html"); else log "summary: no assembly phase in $results; its pages are not expected"; fi
  for f in "${pages[@]}"; do
    [ -s "$dir/$f" ] || { log "summary: missing or empty $dir/$f"; fails=$((fails + 1)); }
  done
  if [ "$fails" -eq 0 ]; then log "summary: ok ($dir)"; else log "summary: FAILED ($fails problem(s), $dir)"; fi
  echo "$dir"
  [ "$fails" -eq 0 ]
}

# ---- validate ---------------------------------------------------------------
cmd_validate() {  # [<member>] [--truth-dir DIR] [--region R] [--min-f1 KIND=VALUE ...]: the member's calls against GIAB's HG002 v5.0q benchmark
  local member=HG002 truth="" region=chr20
  local -a floors=()
  while [ $# -gt 0 ]; do
    case "$1" in
      --truth-dir) truth="$2"; shift 2 ;;
      --region) region="$2"; shift 2 ;;
      --min-f1) floors+=("$2"); shift 2 ;;
      -*) die "validate: unknown option $1" ;;
      *) member="$1"; shift ;;
    esac
  done
  set_ugc
  local code fasta refs_root run
  code=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["code_dir"])' "$project/.ugc-wgw/config.json")
  fasta=$(python3 -c 'import json,sys
m = dict(l.rstrip("\n").split("\t", 1) for l in open(json.load(open(sys.argv[1]))["ref_map_file"]) if "\t" in l)
print(m["fasta"])' "$project/.ugc-wgw/config.json")
  refs_root=${fasta%%/hifi-wdl-resources-*}
  if [ -z "$truth" ]; then
    truth=$refs_root/$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["bundles"]["ugc-wgw-truth"]["install_subdir"])' "$code/references.lock")
  fi
  local smvar=$truth/HG002_GRCh38_v5.0q_smvar.$region.vcf.gz smbed=$truth/HG002_GRCh38_v5.0q_smvar.benchmark.$region.bed
  local stvar=$truth/HG002_GRCh38_v5.0q_stvar.$region.vcf.gz stbed=$truth/HG002_GRCh38_v5.0q_stvar.benchmark.$region.bed
  local f
  for f in "$smvar" "$smvar.tbi" "$smbed" "$stvar" "$stvar.tbi" "$stbed"; do
    [ -s "$f" ] || die "validate: no $f (the ugc-wgw-truth entry of references.lock: tests/smoke/prepare.sh --references-only installs it; --truth-dir names another directory)"
  done
  # the member's newest successful run of each per-sample stage with phased calls, at any version: the standalone route
  # (singleton: DeepVariant's own call set) and the joint one (downstream: the member's slice of the GLnexus VCF) score differently
  local runs
  runs=$(python3 -c 'import glob, json, os, sys
for stage in ("singleton", "downstream"):
    best = None
    for man in glob.glob(os.path.join(sys.argv[1], "samples", sys.argv[2], "*", stage, "current", "run_manifest.json")):
        try:
            ok = json.load(open(man)).get("status") == "success"
        except (OSError, ValueError):
            ok = False
        if ok and (best is None or os.path.getmtime(man) > best[0]):
            best = (os.path.getmtime(man), os.path.dirname(os.path.realpath(man)))
    if best:
        print(stage + "\t" + best[1])' "$results" "$member")
  [ -n "$runs" ] || die "validate: $member has no successful singleton or downstream run under $results"
  local base_sif aardvark_sif truvari_sif
  base_sif=$(sif_of "$code/vendor/hifi-human-wgs-wdl/image_manifest.txt" quay.io/pacbio/pb_wdl_base)
  aardvark_sif=$(sif_of "$code/image_manifest.ugc-wgw.txt" quay.io/biocontainers/aardvark)
  truvari_sif=$(sif_of "$code/image_manifest.ugc-wgw.txt" quay.io/biocontainers/truvari)
  local stamp out threads stage run small sv tag
  stamp=$(date -u +%Y%m%dT%H%M%SZ); out=$smoke/logs/validate-$stamp; mkdir -p "$out"
  threads=$(nproc); [ "$threads" -le 8 ] || threads=8
  local -a bind=(--bind "$smoke,$code,$truth,$refs_root")
  while IFS=$'\t' read -r stage run; do
    small=$(find "$run/out/phased_small_variant_vcf" -name '*.vcf.gz' 2>/dev/null | head -1)
    sv=$(find "$run/out/phased_sv_vcf" -name '*.vcf.gz' 2>/dev/null | head -1)
    if [ -z "$small" ] || [ -z "$sv" ]; then log "validate: $run has no phased small-variant or SV VCF; skipped"; continue; fi
    tag=$member-$stage
    mkdir -p "$out/in-$tag"
    log "validate $member ($stage): $run against GIAB HG002 v5.0q on $region ($truth)"
    # the indices live in their own output directories: both next to each VCF here
    ln -sfn "$small" "$out/in-$tag/small.vcf.gz"; ln -sfn "${small%/phased_small_variant_vcf/*}/phased_small_variant_vcf_index/$(basename "$small").tbi" "$out/in-$tag/small.vcf.gz.tbi"
    ln -sfn "$sv" "$out/in-$tag/sv.vcf.gz"; ln -sfn "${sv%/phased_sv_vcf/*}/phased_sv_vcf_index/$(basename "$sv").tbi" "$out/in-$tag/sv.vcf.gz.tbi"
    # the query: PASS records of the region where the member carries an alternate allele, genotypes only (the joint slice
    # keeps hom-ref sites and GLnexus's FORMAT fields, which aardvark does not read)
    apptainer exec "${bind[@]}" "$base_sif" bcftools view --no-version --apply-filters PASS,. --regions "$region" --include 'GT="alt"' \
        --output-type u "$out/in-$tag/small.vcf.gz" \
      | apptainer exec "${bind[@]}" "$base_sif" bcftools annotate --no-version --remove 'INFO,^FORMAT/GT' --output-type z \
        --output "$out/in-$tag/small.query.vcf.gz" -
    apptainer exec "${bind[@]}" "$base_sif" tabix --force --preset vcf "$out/in-$tag/small.query.vcf.gz"
    apptainer exec "${bind[@]}" "$aardvark_sif" aardvark compare --threads "$threads" --reference "$fasta" --truth-vcf "$smvar" \
        --query-vcf "$out/in-$tag/small.query.vcf.gz" --regions "$smbed" --output-dir "$out/aardvark-$tag" > "$out/aardvark-$tag.log" 2>&1 \
      || die "validate: aardvark compare failed (see $out/aardvark-$tag.log)"
    apptainer exec "${bind[@]}" "$truvari_sif" truvari bench -b "$stvar" -c "$out/in-$tag/sv.vcf.gz" --includebed "$stbed" --passonly \
        --sizemin 50 --sizefilt 50 --pick ac -o "$out/truvari-$tag" > "$out/truvari-$tag.log" 2>&1 \
      || die "validate: truvari bench failed (see $out/truvari-$tag.log)"
    # GIAB's README for v5.0q recommends the refine step: the benchmark keeps complex SVs, whose representations differ
    apptainer exec "${bind[@]}" "$truvari_sif" truvari refine --use-original-vcfs --reference "$fasta" --threads "$threads" \
        "$out/truvari-$tag" >> "$out/truvari-$tag.log" 2>&1 \
      || die "validate: truvari refine failed (see $out/truvari-$tag.log)"
    printf '%s\t%s\n' "$stage" "$run" >> "$out/runs.tsv"
  done <<< "$runs"
  [ -s "$out/runs.tsv" ] || die "validate: nothing scored"
  python3 - "$out" "$member" "$code" "$region" "${floors[@]+"${floors[@]}"}" <<'PY'
import csv, datetime, gzip, json, sys
out, member, code, region = sys.argv[1:5]
floors = dict(f.split("=", 1) for f in sys.argv[5:])
f4 = lambda v: "%.4f" % float(v)

def header_value(path, prefix):
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("##"):
                return ""
            if line.startswith(prefix):
                return line[len(prefix):].strip().strip('"')
    return ""

upstream = json.load(open(f"{code}/upstream.lock")).get("sources") or {}
tag = next((v.get("tag", "") for v in upstream.values() if isinstance(v, dict)), "")
small_rows, sv_rows, md, bad = [], [], [], []
unknown = [k for k in floors if k not in ("snv", "indel", "sv")]
if unknown:
    sys.exit("validate: --min-f1 kinds are snv, indel, sv; not " + ", ".join(unknown))
for line in open(f"{out}/runs.tsv"):
    stage, run = line.rstrip("\n").split("\t")
    t = f"{member}-{stage}"
    by = {}
    for r in csv.DictReader(open(f"{out}/aardvark-{t}/summary.tsv"), delimiter="\t"):
        if r["comparison"] == "GT" and r["region_label"] == "ALL" and r["variant_type"] in ("Snv", "Insertion", "Deletion", "JointIndel"):
            kind = {"Snv": "snv", "JointIndel": "indel"}.get(r["variant_type"], r["variant_type"].lower())
            by[kind] = [stage, kind, r["truth_total"], r["query_total"], r["truth_tp"], r["truth_fn"], r["query_fp"], f4(r["metric_recall"]),
                        f4(r["metric_precision"]), f4(r["metric_f1"])]
    small_rows += [by[k] for k in ("snv", "indel", "insertion", "deletion") if k in by]
    sv = []
    for step, path in (("bench", "summary.json"), ("refine", "refine.variant_summary.json")):
        d = json.load(open(f"{out}/truvari-{t}/{path}"))
        gt = d.get("gt_concordance")
        sv.append([stage, step, str(d["base cnt"]), str(d["comp cnt"]), str(d["TP-base"]), str(d["FN"]), str(d["FP"]), f4(d["recall"]), f4(d["precision"]),
                   f4(d["f1"]), f4(gt) if isinstance(gt, (int, float)) else "NA"])
    sv_rows += sv
    man = json.load(open(f"{run}/run_manifest.json"))
    dv = header_value(f"{out}/in-{t}/small.vcf.gz", "##DeepVariant_version=")
    glnexus = header_value(f"{out}/in-{t}/small.vcf.gz", "##GLnexusVersion=")
    caller = ("DeepVariant " + dv if dv else "") + (" GLnexus " + glnexus if glnexus else "")
    # pbdv writes DeepVariant's header verbatim: the run manifest's tools block (0.19.0) says which of them called, and
    # which aligner made the BAM (named in the run column when it is not pbmm2)
    tools = man.get("tools") or {}
    if tools.get("small_variant_caller") == "pbdv":
        caller = f"pbdv ({tools.get('small_variant_device', '?')}; DeepVariant {dv} model)"
    elif tools.get("small_variant_caller") == "parabricks_deepvariant":
        caller = f"Parabricks {caller}".strip()
    route = f", {tools['aligner']}" if tools.get("aligner") not in (None, "pbmm2") else ""
    md.append("| " + " | ".join([datetime.date.today().isoformat(), str((man.get("ugc_pacbio_wgw") or {}).get("version", "")), tag, f"{member} {stage} ({region}{route})",
                                  caller.strip() or "-", header_value(f"{out}/in-{t}/sv.vcf.gz", "##source=") or "-",
                                  " / ".join(by["snv"][7:10]), " / ".join(by["indel"][7:10]), " / ".join(sv[0][7:10]), " / ".join(sv[1][7:10]), sv[0][10]]) + " |")
    f1 = {"snv": float(by["snv"][9]), "indel": float(by["indel"][9]), "sv": float(sv[1][9])}
    bad += [f"{stage}: {k} F1 {f1[k]:.4f} below the floor {v}" for k, v in floors.items() if f1[k] < float(v)]
with open(f"{out}/small.tsv", "w") as fh:
    fh.write("stage\tkind\ttruth\tquery\tTP\tFN\tFP\trecall\tprecision\tF1\n" + "".join("\t".join(r) + "\n" for r in small_rows))
with open(f"{out}/sv.tsv", "w") as fh:
    fh.write("stage\tstep\ttruth\tquery\tTP\tFN\tFP\trecall\tprecision\tF1\tgt_concordance\n" + "".join("\t".join(r) + "\n" for r in sv_rows))
with open(f"{out}/validation_rows.md", "w") as fh:
    fh.write("\n".join(md) + "\n")
if bad:
    sys.exit("validate: " + "; ".join(bad))
PY
  log "validate: small variants ($member against the v5.0q small-variant benchmark on $region, genotype match, aardvark):"
  column -t -s $'\t' "$out/small.tsv" >&2
  log "validate: structural variants >= 50 bp (truvari bench, then refine):"
  column -t -s $'\t' "$out/sv.tsv" >&2
  log "validate: the rows for tests/smoke/VALIDATION.md:"
  cat "$out/validation_rows.md" >&2
  log "validate: written to $out"
}

# ---- eval -------------------------------------------------------------------
# the evaluations: <name> -> the init options of its project, and the tasks whose resources the comparison prints
eval_init_args() {
  case "$1" in
    sawfish23) echo "--profile eval-sawfish23" ;;
    pblab) echo "--aligner pbweave --deepvariant pbdv" ;;        # what the pblab profile sets, as flags (a flag wins over a profile)
    pblab-cpu) echo "--aligner pbweave --deepvariant pbdv-cpu" ;;
    pbdv) echo "--deepvariant pbdv" ;;                           # the caller alone: alignments replay from the cache
    pbdv-cpu) echo "--deepvariant pbdv-cpu" ;;
    pbweave) echo "--aligner pbweave" ;;                         # the aligner alone, DeepVariant on its BAM
    *) return 1 ;;
  esac
}
eval_tasks() { case "$1" in sawfish23) echo "sawfish" ;; *) echo "pbmm2|pbweave|deepvariant|pbdv|pbsamoa_merge" ;; esac; }

cmd_eval() {  # cmd_eval <name> [compare]: a second project with the evaluation's tools, its standalone phase, validate, then the comparison; "compare" = the comparison alone, on finished runs
  local what=${1:-} step=${2:-} init_args
  init_args=$(eval_init_args "$what") || die "eval: the known evaluations: sawfish23 pblab pblab-cpu pbdv pbdv-cpu pbweave"
  [ -z "$step" ] || [ "$step" = compare ] || die "eval: the second word is 'compare' or nothing"
  [ -z "$variant" ] || die "eval: run it without UGC_WGW_SMOKE_VARIANT (it sets the variant itself)"
  set_ugc   # the production project is the base of the comparison
  if [ -z "$step" ]; then
    log "eval $what: project-$what (init $init_args; the same data, references, cache and SIFs), then the comparison with $project"
    # shellcheck disable=SC2086  # the options are words
    UGC_WGW_SMOKE_VARIANT=$what "$0" init --force $init_args
    UGC_WGW_SMOKE_VARIANT=$what "$0" all standalone || log "eval $what: the standalone phase reported failures (see above); comparing what finished"
    UGC_WGW_SMOKE_VARIANT=$what "$0" validate || log "eval $what: validate did not finish (see above)"
  else
    [ -d "$smoke/project-$what/.ugc-wgw" ] || die "eval: no project-$what yet; run 'eval $what' first"
    log "eval $what: the comparison alone (project-$what against $project)"
  fi
  eval_compare "$what" "$results" "$smoke/results-$what" "$smoke/project-$what"
}

cmd_entry() {  # the entry kinds: a variant project whose members enter as an aligned BAM (HG002b) and as BAM + VCFs (HG003v); then a bad VCF
  [ -z "$variant" ] || die "entry: run it without UGC_WGW_SMOKE_VARIANT (it sets the variant itself)"
  set_ugc
  local version; version=$(cat "$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["code_dir"])' "$project/.ugc-wgw/config.json")/VERSION")
  local up2=$results/samples/HG002/$version/upstream/current/out up3=$results/samples/HG003/$version/upstream/current/out
  if [ ! -d "$up2" ] || [ ! -d "$up3" ]; then die "entry: the main project needs successful upstream runs of HG002 and HG003 (run.sh all joint)"; fi
  local stamp out; stamp=$(date -u +%Y%m%dT%H%M%SZ); out=$smoke/logs/entry-$stamp; mkdir -p "$out"
  local base_sif; base_sif=$(sif_of "$repo/vendor/hifi-human-wgs-wdl/image_manifest.txt" quay.io/pacbio/pb_wdl_base)
  # the read-group sample name must be the sample id (the callers name their outputs after it): the main project's BAMs
  # carry SM:HG002 and SM:HG003, so the variant members get reheadered copies (samtools reheader rewrites the header only)
  local b2 b2i b3 b3i sv svi smv smvi gv gvi
  b2=$(ls "$up2"/aligned_hifi_reads/*.bam); b3=$(ls "$up3"/aligned_hifi_reads/*.bam)
  reheader_bam() {  # <in.bam> <old SM> <new SM> <out.bam>
    apptainer exec --bind "$smoke","$repo" "$base_sif" samtools reheader -c "sed -E 's/(\tSM:)$2(\t|\$)/\1$3\2/'" "$1" > "$4"
    apptainer exec --bind "$smoke","$repo" "$base_sif" samtools index "$4"
  }
  reheader_bam "$b2" HG002 HG002b "$out/HG002b.aligned.bam"; b2=$out/HG002b.aligned.bam; b2i=$b2.bai
  reheader_bam "$b3" HG003 HG003v "$out/HG003v.aligned.bam"; b3=$out/HG003v.aligned.bam; b3i=$b3.bai
  smv=$(ls "$up3"/small_variant_vcf/*.vcf.gz); smvi=$(ls "$up3"/small_variant_vcf_index/*.tbi)
  sv=$(ls "$up3"/sv_vcf/*.vcf.gz); svi=$(ls "$up3"/sv_vcf_index/*.tbi)
  gv=$(ls "$up3"/small_variant_gvcf/*.g.vcf.gz); gvi=$(ls "$up3"/small_variant_gvcf_index/*.tbi)
  printf 'sample_id\tsex\taligned_bam\taligned_bam_index\tsmall_variant_vcf\tsmall_variant_vcf_index\tsv_vcf\tsv_vcf_index\tgvcf\tgvcf_index\n' > "$out/samples.tsv"
  printf 'HG002b\tMALE\t%s\t%s\t\t\t\t\t\t\n' "$b2" "$b2i" >> "$out/samples.tsv"
  printf 'HG003v\tMALE\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$b3" "$b3i" "$smv" "$smvi" "$sv" "$svi" "$gv" "$gvi" >> "$out/samples.tsv"
  printf 'HG002b\nHG003v\n' > "$out/cohort.txt"
  log "entry: project-entry with HG002b (bam kind) and HG003v (vcf kind) from the main project's upstream outputs, BAMs reheadered to their ids (sheet $out/samples.tsv)"
  UGC_WGW_SMOKE_VARIANT=entry UGC_WGW_SMOKE_SHEET=$out/samples.tsv UGC_WGW_SMOKE_COHORT=$out/cohort.txt "$0" init --force
  local rc=0
  UGC_WGW_SMOKE_VARIANT=entry UGC_WGW_SMOKE_SHEET=$out/samples.tsv UGC_WGW_SMOKE_COHORT=$out/cohort.txt "$0" all standalone || rc=1
  # HG003v's phased VCFs against HG003's singleton ones: the same BAM and the same unphased calls went into HiPhase
  local eres=$smoke/results-entry f3 fv n3 nv m3 mv
  log "entry: HG003v (vcf kind, downstream) against HG003 (singleton), bodies without header lines:"
  printf 'file\trecords_singleton\trecords_entry\tidentical\n' > "$out/compare.tsv"
  for name in phased_small_variant_vcf phased_sv_vcf; do
    local -a c3=("$results"/samples/HG003/"$version"/singleton/current/out/"$name"/*.vcf.gz) cv=("$eres"/samples/HG003v/"$version"/downstream/current/out/"$name"/*.vcf.gz)
    f3=${c3[0]}; fv=${cv[0]}
    if [ ! -e "$f3" ] || [ ! -e "$fv" ]; then printf '%s\t-\t-\tmissing\n' "$name" >> "$out/compare.tsv"; continue; fi
    n3=$(zcat "$f3" | grep -vc '^#'); nv=$(zcat "$fv" | grep -vc '^#')
    m3=$(zcat "$f3" | grep -v '^#' | sed 's/HG003v/HG003/g' | md5sum | cut -c1-32); mv=$(zcat "$fv" | grep -v '^#' | sed 's/HG003v/HG003/g' | md5sum | cut -c1-32)
    printf '%s\t%s\t%s\t%s\n' "$name" "$n3" "$nv" "$([ "$m3" = "$mv" ] && echo yes || echo no)" >> "$out/compare.tsv"
  done
  column -t -s $'\t' "$out/compare.tsv" >&2
  # the negative cases: a BAM whose read group is another sample is refused at registration; HG003's small-variant VCF
  # with its contigs renamed (chr20 -> 20) must fail intake as an input error
  local -a eugc=(python3 "$repo/bin/ugc-wgw" --project "$smoke/project-entry")
  printf 'sample_id\tsex\taligned_bam\taligned_bam_index\n' > "$out/badrg.tsv"
  printf 'HG003r\tMALE\t%s\t%s\n' "$b3" "$b3i" >> "$out/badrg.tsv"
  local rrc=0
  "${eugc[@]}" --color never samples add "$out/badrg.tsv" > "$out/HG003r.add.log" 2>&1 || rrc=$?
  if [ "$rrc" -ne 0 ] && grep -q '@RG SM is HG003v, not the sample id' "$out/HG003r.add.log"; then
    log "entry: HG003r (HG003v's BAM under another id) refused at registration, as expected"
  else
    log "entry: HG003r (HG003v's BAM under another id) was NOT refused (exit $rrc; log $out/HG003r.add.log)"; rc=1
  fi
  reheader_bam "$b3" HG003v HG003x "$out/HG003x.aligned.bam"
  zcat "$smv" | sed -e 's/^chr20\t/20\t/' -e 's/ID=chr20,/ID=20,/' | apptainer exec --bind "$smoke","$repo" "$base_sif" bgzip -c > "$out/HG003x.small_variants.vcf.gz"
  apptainer exec --bind "$smoke","$repo" "$base_sif" tabix -p vcf "$out/HG003x.small_variants.vcf.gz"
  printf 'sample_id\tsex\taligned_bam\taligned_bam_index\tsmall_variant_vcf\tsmall_variant_vcf_index\tsv_vcf\tsv_vcf_index\n' > "$out/bad.tsv"
  printf 'HG003x\tMALE\t%s\t%s\t%s\t%s\t%s\t%s\n' "$out/HG003x.aligned.bam" "$out/HG003x.aligned.bam.bai" "$out/HG003x.small_variants.vcf.gz" "$out/HG003x.small_variants.vcf.gz.tbi" "$sv" "$svi" >> "$out/bad.tsv"
  "${eugc[@]}" samples add "$out/bad.tsv"
  local xrc=0
  "${eugc[@]}" --color never submit --mode standalone --stage intake --samples HG003x > "$out/HG003x.intake.log" 2>&1 || xrc=$?
  local kind; kind=$(python3 -c 'import json,sys,glob; p=sorted(glob.glob(sys.argv[1]))[-1]; d=json.load(open(p)); print(d.get("error", {}).get("kind"), d.get("status"))' "$eres/samples/HG003x/$version/intake/attempt-*/run_manifest.json" 2>/dev/null || echo none)
  log "entry: HG003x (contigs renamed) intake -> $kind (expected 'input failed'; the driver exited $xrc; log $out/HG003x.intake.log)"
  case "$kind" in "input failed") ;; *) rc=1 ;; esac
  log "entry: done, exit $rc ($out)"
  return $rc
}

eval_compare() {  # <what> <base results> <comp results> <comp project>: per member truvari bench of the phased SV VCFs, the small-variant and TRGT call sets and the digests side by side; the changed tasks' resources
  local what=$1 base=$2 comp=$3 comp_project=$4
  local code version stamp out truvari_sif base_sif tasks
  code=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["code_dir"])' "$project/.ugc-wgw/config.json")
  version=$(cat "$code/VERSION")
  stamp=$(date -u +%Y%m%dT%H%M%SZ)
  out=$smoke/logs/eval-$what-$stamp
  mkdir -p "$out"
  truvari_sif=$(sif_of "$code/image_manifest.ugc-wgw.txt" quay.io/biocontainers/truvari)
  base_sif=$(sif_of "$code/vendor/hifi-human-wgs-wdl/image_manifest.txt" quay.io/pacbio/pb_wdl_base)
  tasks=$(eval_tasks "$what")
  # the benchmark regions of the truth entry (the smoke's chr20): the call sets are compared inside and outside them,
  # because outside (centromere, satellites, segmental duplications) two aligners place reads differently and neither is scored
  local region=chr20 fasta truth smbed="" stbed=""
  fasta=$(python3 -c 'import json,sys
m = dict(l.rstrip("\n").split("\t", 1) for l in open(json.load(open(sys.argv[1]))["ref_map_file"]) if "\t" in l)
print(m["fasta"])' "$project/.ugc-wgw/config.json")
  truth=${fasta%%/hifi-wdl-resources-*}/$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["bundles"]["ugc-wgw-truth"]["install_subdir"])' "$code/references.lock")
  [ ! -s "$truth/HG002_GRCh38_v5.0q_smvar.benchmark.$region.bed" ] || smbed=$truth/HG002_GRCh38_v5.0q_smvar.benchmark.$region.bed
  [ ! -s "$truth/HG002_GRCh38_v5.0q_stvar.benchmark.$region.bed" ] || stbed=$truth/HG002_GRCh38_v5.0q_stvar.benchmark.$region.bed
  local m b c
  printf 'member\tscope\tbase_calls\tcomp_calls\tTP_base\tTP_comp\tFP\tFN\tprecision\trecall\tF1\tgt_concordance\n' > "$out/sv.tsv"
  for m in $(members); do
    b=$(find "$base/samples/$m/$version/singleton/current/out/phased_sv_vcf" -name '*.vcf.gz' 2>/dev/null | head -1)
    c=$(find "$comp/samples/$m/$version/singleton/current/out/phased_sv_vcf" -name '*.vcf.gz' 2>/dev/null | head -1)
    if [ -z "$b" ] || [ -z "$c" ]; then log "eval: $m has no phased SV VCF in both projects; skipped"; continue; fi
    # truvari wants the index next to the VCF; the run keeps it in its own output directory (phased_sv_vcf_index)
    mkdir -p "$out/in-$m"
    ln -sfn "$b" "$out/in-$m/base.vcf.gz"; ln -sfn "${b%/phased_sv_vcf/*}/phased_sv_vcf_index/$(basename "$b").tbi" "$out/in-$m/base.vcf.gz.tbi"
    ln -sfn "$c" "$out/in-$m/comp.vcf.gz"; ln -sfn "${c%/phased_sv_vcf/*}/phased_sv_vcf_index/$(basename "$c").tbi" "$out/in-$m/comp.vcf.gz.tbi"
    local scope
    for scope in all benchmark; do
      local -a bed=()
      if [ "$scope" = benchmark ]; then
        [ -n "$stbed" ] || continue
        bed=(--includebed "$stbed")
      fi
      apptainer exec --bind "$smoke","$code","$truth" "$truvari_sif" truvari bench -b "$out/in-$m/base.vcf.gz" -c "$out/in-$m/comp.vcf.gz" -o "$out/truvari-$m-$scope" \
          --sizemin 50 --sizefilt 50 --pctseq 0 "${bed[@]}" > "$out/truvari-$m-$scope.log" 2>&1 || { log "eval: truvari bench failed for $m (see $out/truvari-$m-$scope.log)"; continue; }
      python3 - "$m" "$scope" "$out/truvari-$m-$scope/summary.json" >> "$out/sv.tsv" <<'PY'
import json, sys
m, scope, d = sys.argv[1], sys.argv[2], json.load(open(sys.argv[3]))
f = lambda k, n=4: ("%.*f" % (n, d[k])) if isinstance(d.get(k), (int, float)) and d.get(k) is not None else "NA"
print("\t".join([m, scope, str(d.get("base cnt", "NA")), str(d.get("comp cnt", "NA")), str(d.get("TP-base", "NA")), str(d.get("TP-comp", "NA")), str(d.get("FP", "NA")),
                 str(d.get("FN", "NA")), f("precision"), f("recall"), f("f1"), f("gt_concordance")]))
PY
    done
  done
  log "eval: SV sets, production (base) vs $what (comp), truvari bench --sizemin 50 --pctseq 0, everywhere and inside the SV benchmark regions of $region:"
  column -t -s $'\t' "$out/sv.tsv" >&2
  # the call sets themselves: PASS small variants (site, alleles, genotype without phase) and TRGT allele lengths per locus
  printf 'member\twhat\tshared\tonly_a\tonly_b\n' > "$out/pairs.tsv"
  local kind q side run f
  for m in $(members); do
    for side in base comp; do
      if [ "$side" = base ]; then run=$base/samples/$m/$version/singleton/current/out; else run=$comp/samples/$m/$version/singleton/current/out; fi
      f=$(find "$run/phased_small_variant_vcf" -name '*.vcf.gz' 2>/dev/null | head -1)
      [ -z "$f" ] || apptainer exec --bind "$smoke","$code" "$base_sif" bcftools query --include '(FILTER="PASS" || FILTER=".") && GT="alt"' \
          --format '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$f" 2>/dev/null | tr '|' '/' | sed 's#\t1/0$#\t0/1#' | LC_ALL=C sort > "$out/in-$m/small.$side.txt"
      f=$(find "$run/phased_trgt_vcf" -name '*.vcf.gz' 2>/dev/null | head -1)
      [ -z "$f" ] || apptainer exec --bind "$smoke","$code" "$base_sif" bcftools query --include 'GT!="mis"' \
          --format '%INFO/TRID[\t%AL]\n' "$f" 2>/dev/null | LC_ALL=C sort > "$out/in-$m/trgt.$side.txt"
    done
    # the small variants in three scopes: inside the benchmark regions, the rest of the region, the rest of the genome
    if [ -s "$out/in-$m/small.base.txt" ] && [ -s "$out/in-$m/small.comp.txt" ]; then
      for side in base comp; do
        if [ -n "$smbed" ]; then
          awk -F'\t' 'NR == FNR { n[$1]++; s[$1, n[$1]] = $2; e[$1, n[$1]] = $3; next }
                      { hit = 0; for (i = 1; i <= n[$1]; i++) if ($2 > s[$1, i] && $2 <= e[$1, i]) { hit = 1; break } print > (hit ? inside : outside) }' \
              inside="$out/in-$m/small-bench.$side.txt" outside="$out/in-$m/small-rest.$side.txt" "$smbed" "$out/in-$m/small.$side.txt"
          touch "$out/in-$m/small-bench.$side.txt" "$out/in-$m/small-rest.$side.txt"
        else
          cp "$out/in-$m/small.$side.txt" "$out/in-$m/small-rest.$side.txt"
        fi
        awk -F'\t' -v r="$region" '$1 == r' "$out/in-$m/small-rest.$side.txt" > "$out/in-$m/small-region.$side.txt"
        awk -F'\t' -v r="$region" '$1 != r' "$out/in-$m/small-rest.$side.txt" > "$out/in-$m/small-off.$side.txt"
      done
    fi
    for kind in small-bench small-region small-off trgt; do
      if [ ! -s "$out/in-$m/$kind.base.txt" ] || [ ! -s "$out/in-$m/$kind.comp.txt" ]; then continue; fi
      case "$kind" in
        small-bench) q="PASS small variants inside the benchmark regions of $region" ;;
        small-region) q="PASS small variants on $region outside them" ;;
        small-off) q="PASS small variants off $region" ;;
        *) q="TRGT loci (allele lengths)" ;;
      esac
      printf '%s\t%s\t%s\t%s\t%s\n' "$m" "$q" "$(LC_ALL=C comm -12 "$out/in-$m/$kind.base.txt" "$out/in-$m/$kind.comp.txt" | wc -l)" \
          "$(LC_ALL=C comm -23 "$out/in-$m/$kind.base.txt" "$out/in-$m/$kind.comp.txt" | wc -l)" \
          "$(LC_ALL=C comm -13 "$out/in-$m/$kind.base.txt" "$out/in-$m/$kind.comp.txt" | wc -l)" >> "$out/pairs.tsv"
    done
  done
  # the digests of both projects, made here so they are of the same code and thresholds
  local -a ids; mapfile -t ids < <(members)
  "$code/bin/ugc-wgw" --project "$project" summary --mode standalone --any-version --samples "${ids[@]}" --out-dir "$out/digests-base" > "$out/summary-base.log" 2>&1 \
    || log "eval: summary of the production project failed (see $out/summary-base.log)"
  "$code/bin/ugc-wgw" --project "$comp_project" summary --mode standalone --any-version --samples "${ids[@]}" --out-dir "$out/digests-$what" > "$out/summary-$what.log" 2>&1 \
    || log "eval: summary of project-$what failed (see $out/summary-$what.log)"
  log "eval: production against $what, member by member (the digests and the call sets; $out/compare.tsv):"
  python3 "$here/compare_runs.py" --a-label production --b-label "$what" --digests-a "$out/digests-base" --digests-b "$out/digests-$what" \
      --pairs "$out/pairs.tsv" --out "$out/compare.tsv" "${ids[@]}" >&2 || log "eval: the digest comparison did not run"
  # the changed tasks' resources in both projects (the call cache replayed everything else)
  "$code/bin/ugc-wgw" --project "$project" usage --by task --any-version --json > "$out/usage-base.json" 2>/dev/null || true
  "$code/bin/ugc-wgw" --project "$comp_project" usage --by task --any-version --json > "$out/usage-$what.json" 2>/dev/null || true
  python3 - "$out/usage-base.json" "$out/usage-$what.json" "$what" "$tasks" > "$out/tasks.tsv" <<'PY'
import json, re, sys
want = re.compile(sys.argv[4])
rows = ["project\ttask\tjobs\telapsed_s\tcore_h_alloc\tpeak_rss_mb"]
for label, path in (("production", sys.argv[1]), (sys.argv[3], sys.argv[2])):
    try:
        by_task = json.load(open(path)).get("by_task", {})
    except (OSError, ValueError):
        continue
    for task, u in sorted(by_task.items()):
        if want.search(task):
            rows.append("\t".join([label, task, str(u.get("jobs", "")), "%.0f" % float(u.get("elapsed_s") or 0), "%.3f" % float(u.get("core_h_alloc") or 0), "%.0f" % float(u.get("peak_rss_mb") or 0)]))
print("\n".join(rows))
PY
  log "eval: the changed tasks' resources (ugc-wgw usage --by task, every run of each project; without sacct the memory is an estimate):"
  column -t -s $'\t' "$out/tasks.tsv" >&2
  log "eval: written to $out (sv.tsv, compare.tsv, pairs.tsv, tasks.tsv, truvari-<member>/); paste the tables into tests/smoke/README.md"
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
  eval) shift; cmd_eval "$@" ;;
  entry) cmd_entry ;;
  validate) shift; cmd_validate "$@" ;;
  -h|--help) usage ;;
  singleton|upstream|intake|cohort_qc|cohort_call|downstream|cohort_merge|cohort_freq|cohort_tr|cohort_methyl|assembly|assembly_call|cohort_pangenome) cmd_stage "$@" ;;
  *) die "unknown command $1" ;;
esac
