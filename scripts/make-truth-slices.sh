#!/usr/bin/env bash
# Make the chr20 slices of GIAB's HG002 v5.0q benchmark that tests/smoke/run.sh validate scores against
# (references.lock entry ugc-wgw-truth; docs/DESIGN.md §18 2026-11-11). Run once per benchmark release on the
# dev machine; the outputs are committed under references/.
#
#   scripts/make-truth-slices.sh --src DIR --out DIR --sif PB_WDL_BASE.sif [--region chr20]
#
#   --src   a directory holding GIAB's files as downloaded from
#           https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/
#           (HG002_GRCh38_v5.0q_{smvar,stvar}.vcf.gz with .tbi, the two .benchmark.bed files, checksum.md5)
#   --sif   the pb_wdl_base SIF of the vendored manifest (bcftools, bgzip, tabix)
#
# Every source file is checked against GIAB's checksum.md5 first. The slices keep the source header (so the
# contig list is the benchmark's) and add one header line naming the source and its md5; the BEDs are the
# region's rows, unchanged. The script prints the md5 and sha256 of what it wrote, for references.lock.
set -euo pipefail

src="" out="" sif="" region="chr20"
while [ $# -gt 0 ]; do
  case "$1" in
    --src) src="$2"; shift 2 ;;
    --out) out="$2"; shift 2 ;;
    --sif) sif="$2"; shift 2 ;;
    --region) region="$2"; shift 2 ;;
    -h|--help) sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown argument $1" >&2; exit 2 ;;
  esac
done
if [ -z "$src" ] || [ -z "$out" ] || [ -z "$sif" ]; then echo "--src, --out and --sif are required" >&2; exit 2; fi
[ -f "$sif" ] || { echo "no SIF at $sif" >&2; exit 1; }
src=$(cd "$src" && pwd -P)
mkdir -p "$out"
out=$(cd "$out" && pwd -P)
release="v5.0q"
tool() { apptainer exec --bind "$src","$out" "$sif" "$@"; }

for kind in smvar stvar; do
  for f in "HG002_GRCh38_${release}_$kind.vcf.gz" "HG002_GRCh38_${release}_$kind.vcf.gz.tbi" "HG002_GRCh38_${release}_$kind.benchmark.bed"; do
    [ -s "$src/$f" ] || { echo "missing $src/$f" >&2; exit 1; }
    want=$(awk -v f="$f" '$2 == f { print $1 }' "$src/checksum.md5")
    have=$(md5sum "$src/$f" | cut -d' ' -f1)
    if [ -z "$want" ] || [ "$want" != "$have" ]; then echo "md5 of $f is $have, GIAB's checksum.md5 says ${want:-nothing}" >&2; exit 1; fi
  done
  vcf="$src/HG002_GRCh38_${release}_$kind.vcf.gz"
  md5=$(md5sum "$vcf" | cut -d' ' -f1)
  dest="$out/HG002_GRCh38_${release}_$kind.$region.vcf.gz"
  echo "##ugc_wgw_truth_slice=region $region of HG002_GRCh38_${release}_$kind.vcf.gz (md5 $md5), GIAB HG002 $release, NIST; made by scripts/make-truth-slices.sh" > "$out/.slice-header.$kind.txt"
  tool bcftools view --no-version --regions "$region" "$vcf" \
    | tool bcftools annotate --no-version --header-lines "$out/.slice-header.$kind.txt" --output-type z --output "$dest" -
  tool tabix --force --preset vcf "$dest"
  rm -f "$out/.slice-header.$kind.txt"
  awk -v r="$region" -F'\t' '$1 == r' "$src/HG002_GRCh38_${release}_$kind.benchmark.bed" > "$out/HG002_GRCh38_${release}_$kind.benchmark.$region.bed"
  [ -s "$out/HG002_GRCh38_${release}_$kind.benchmark.$region.bed" ] || { echo "no $region rows in the $kind benchmark BED" >&2; exit 1; }
  n=$(tool bcftools view --no-header "$dest" | wc -l)
  echo "$kind: $n records, $(wc -l < "$out/HG002_GRCh38_${release}_$kind.benchmark.$region.bed") benchmark regions on $region" >&2
done
( cd "$out" && for f in HG002_GRCh38_"${release}"_*; do printf '%s\t%s\t%s\t%s\n' "$f" "$(stat -c %s "$f")" "$(md5sum "$f" | cut -d' ' -f1)" "$(sha256sum "$f" | cut -d' ' -f1)"; done )
