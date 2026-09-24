version 1.0

# ugc-pacbio-wgw — turn references/scatter_regions.<ref>.tsv into per-shard BEDs, contig lists and svx target positions
# origin: none
# see docs/DESIGN.md §7.1

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_regions_prepare {
  meta {
    description: "Group the rows of a scatter_regions TSV (#shard chrom start end) by shard and emit, in file order, one BED per shard plus the contig list and svx --target-positions string for each shard"
    outputs: {
      shard_names: {
        description: "Shard names, in TSV order"
      },
      shard_beds: {
        description: "One 3-column BED per shard, same order"
      },
      shard_contigs: {
        description: "Comma-separated contig names per shard (trgt merge --contig), same order"
      },
      shard_targets: {
        description: "svx --target-positions value per shard, 1-based inclusive, same order"
      }
    }
  }

  parameter_meta {
    scatter_regions: {
      description: "TSV with columns shard, chrom, start (0-based), end (exclusive); '#' lines ignored"
    }
    require_whole_contigs: {
      description: "Fail if a contig is split across shards (needed by tools that only scope by contig, e.g. trgt merge)"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File scatter_regions
    Boolean require_whole_contigs = false
    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  Int mem_gb = 1
  Int disk_size = 1

  command <<<
    set -euo pipefail

    python3 --version >&2

    mkdir -p beds

    python3 - "~{scatter_regions}" "~{require_whole_contigs}" << 'EOF'
    import sys
    tsv, whole = sys.argv[1], sys.argv[2] == "true"
    shards, order = {}, []
    for line in open(tsv):
        if not line.strip() or line.startswith("#"):
            continue
        shard, chrom, start, end = line.rstrip("\n").split("\t")[:4]
        if shard not in shards:
            shards[shard] = []
            order.append(shard)
        shards[shard].append((chrom, int(start), int(end)))
    if not order:
        sys.exit(f"{tsv}: no shards")
    seen = {}
    for shard in order:
        for chrom, _, _ in shards[shard]:
            if whole and seen.get(chrom, shard) != shard:
                sys.exit(f"contig {chrom} appears in shards {seen[chrom]} and {shard}; whole-contig shards required")
            seen[chrom] = shard
    with open("shard_names.txt", "w") as names, open("shard_beds.txt", "w") as beds, \
         open("shard_contigs.txt", "w") as contigs, open("shard_targets.txt", "w") as targets:
        for shard in order:
            rows = shards[shard]
            with open(f"beds/{shard}.bed", "w") as bed:
                for chrom, s, e in rows:
                    bed.write(f"{chrom}\t{s}\t{e}\n")
            names.write(shard + "\n")
            beds.write(f"beds/{shard}.bed\n")
            contigs.write(",".join(dict.fromkeys(c for c, _, _ in rows)) + "\n")
            targets.write(",".join(f"({c}:{s + 1}-{e})" for c, s, e in rows) + "\n")
    EOF
  >>>

  output {
    Array[String] shard_names = read_lines("shard_names.txt")
    Array[File] shard_beds = read_lines("shard_beds.txt")
    Array[String] shard_contigs = read_lines("shard_contigs.txt")
    Array[String] shard_targets = read_lines("shard_targets.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}
