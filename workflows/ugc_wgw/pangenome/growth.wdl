version 1.0
# ugc-pacbio-wgw — how one chromosome's graph grows with its members: panacus coverage histograms and the growth in cohort order
# origin: none (panacus 0.5.2 from the Cactus image; commands verified on the dev machine 2026-11-13)
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8.2
import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_pangenome_growth {
  meta {
    description: "panacus on one chromosome's GFA, the members' paths grouped per sample (PanSN sample name) and the reference's path left out of the groups: `panacus hist` gives, for base pairs and for nodes, how much of the graph is carried by exactly c members, once for everything the members carry and once without what the reference path touches (--exclude: the sequence the reference lacks); `panacus ordered-histgrowth` gives the cumulative amount in the cohort's order, which is the ranking's order when the cohort was frozen from one. The histograms are all a growth curve needs: the expected amount in k of N members drawn at random follows from them (the summary computes it with N = the cohort's size, so the chromosomes add up even when a member has no sequence in one, as an XX member on chrY). One long table per chromosome."
    outputs: {
      growth: {
        description: "<cohort_id>.pangenome.<ref_name>.<chromosome>.growth.tsv: table (info, hist, ordered), chromosome, set (all, nonref), count (bp, node, members), key (a coverage, a member's sample ID, or cohort/present), value"
      }
    }
  }

  parameter_meta {
    cohort_id: {
      description: "Cohort ID; output prefix"
    }
    ref_name: {
      description: "Reference name; part of the output name"
    }
    chromosome: {
      description: "The chromosome of this graph"
    }
    gfa: {
      description: "The chromosome's graph (GFA 1.1, gzipped): W lines per path fragment, sample#haplotype#contig"
    }
    name_map: {
      description: "sample_id, safe_name, trio per member in cohort order (the names task): the order of the ordered growth, and the way back from the graph's names to sample IDs"
    }
    reference: {
      description: "The reference genome's name in the graph"
    }
    threads: {
      description: "panacus threads"
    }
    mem_gb: {
      description: "Memory (GB); default 8 plus ten times the gzipped GFA"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String cohort_id
    String ref_name
    String chromosome
    File gfa
    File name_map
    String reference = "GRCh38"
    Int threads = 4
    Int? mem_gb
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{cohort_id}.pangenome.~{ref_name}.~{chromosome}"
  Int mem = select_first([mem_gb, ceil(size(gfa, "GB") * 10 + 8)])
  Int disk_size = ceil(size(gfa, "GB") * 2 + 20)

  command <<<
    set -euo pipefail
    panacus --version | head -1 >&2
    out="~{out_prefix}.growth.tsv"
    printf 'table\tchromosome\tset\tcount\tkey\tvalue\n' > "$out"

    # the graph's paths: the members' (the groups) and the reference's (what --exclude takes away)
    : > ref.paths; : > member.paths; : > present.raw
    zcat -f "~{gfa}" | awk -F'\t' -v ref="~{reference}" '
      $1 == "W" { p = $2 "#" $3 "#" $4; s = $2 }
      $1 == "P" { p = $2; split($2, a, "#"); s = a[1] }
      $1 == "W" || $1 == "P" { if (s == ref) print p > "ref.paths"; else { print p > "member.paths"; print s > "present.raw" } }'
    sort -u present.raw > present.txt
    n_cohort=$(awk 'NR > 1' "~{name_map}" | wc -l)
    printf 'info\t~{chromosome}\tall\tmembers\tcohort\t%s\n' "$n_cohort" >> "$out"
    printf 'info\t~{chromosome}\tall\tmembers\tpresent\t%s\n' "$(wc -l < present.txt)" >> "$out"
    if [ ! -s member.paths ]; then
      echo "no member path in ~{chromosome}: the table holds the two info rows only" >&2
      exit 0
    fi

    # panacus prints '#' comment lines, four header lines, then one row per coverage (hist) or per group (ordered)
    hist() {  # hist <set> <count> [panacus options...]
      local set=$1 count=$2; shift 2
      panacus hist -S -c "$count" -s member.paths "$@" -t ~{threads} "~{gfa}" > hist.raw
      grep -v '^#' hist.raw | tail -n +5 | awk -F'\t' -v OFS='\t' -v s="$set" -v c="$count" \
        '$1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ { print "hist", "~{chromosome}", s, c, $1, $2 }' >> "$out"
    }
    hist all bp
    hist all node
    hist nonref bp -e ref.paths
    hist nonref node -e ref.paths

    # the cohort's order, members with a path in this chromosome only; panacus names a group by the graph's sample name
    awk -F'\t' 'NR > 1 { print $2 }' "~{name_map}" | { grep -F -x -f present.txt || true; } > order.txt
    ordered() {  # ordered <set> [panacus options...]
      local set=$1; shift
      panacus ordered-histgrowth -S -O order.txt -c bp -s member.paths "$@" -t ~{threads} "~{gfa}" > ordered.raw
      grep -v '^#' ordered.raw | tail -n +5 | awk -F'\t' -v OFS='\t' -v s="$set" '
        NR == FNR { if (FNR > 1) id[$2] = $1; next }
        $1 != "No_Path" && $2 ~ /^[0-9]+$/ { print "ordered", "~{chromosome}", s, "bp", (($1 in id) ? id[$1] : $1), $2 }' "~{name_map}" - >> "$out"
    }
    if [ -s order.txt ]; then
      ordered all
      ordered nonref -e ref.paths
    fi
    grep -c . "$out" >&2
  >>>

  output {
    File growth = "~{out_prefix}.growth.tsv"
  }

  runtime {
    docker: "quay.io/comparative-genomics-toolkit/cactus@sha256:fb4f1bab95d9c3b003d1036989f977e8cde2ffcf90a27691b992bd7d2d1a253d"  # v3.3.0 (panacus 0.5.2)
    cpu: threads
    memory: "~{mem} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}
