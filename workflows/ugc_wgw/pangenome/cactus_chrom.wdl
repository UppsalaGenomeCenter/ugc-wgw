version 1.0
# ugc-pacbio-wgw — Minigraph-Cactus steps 4-7 for one chromosome: construct, map, align, join (VCF, vcfwave VCF, GBZ, GFA)
# origin: none (cactus v3.3.0's own batch commands, doc/pangenome.md "--mgSplit"; verified on the dev machine 2026-10-20)
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8.2
import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_cactus_chrom {
  meta {
    description: "One chromosome's graph from the split step's tar: the absolute paths of its chromfile line and seqfile rewritten to this task, then cactus-minigraph --batch, cactus-graphmap --batch, cactus-align --pangenome --outVG --batch and cactus-graphmap-join --vcf --vcfwave --gbz --gfa with the minigraph GFA, every call single-machine with --maxCores/--maxMemory at the allocation and HOME, TMPDIR and --workDir in the task directory. No --inputContigSizes (the exclusion report needs cactus-pangenome's own table); the chromosome's graph-stats.tsv is kept."
    outputs: {
      vcf: {
        description: "<cohort_id>.pangenome.<ref_name>.<chromosome>.vcf.gz: the chromosome's cohort VCF (sample columns = Cactus-safe names, haplotypes phased in .1|.2 order)"
      },
      vcf_index: {
        description: "Tabix index"
      },
      wave_vcf: {
        description: "The same after vcfwave (alleles realigned, complex records decomposed)"
      },
      wave_vcf_index: {
        description: "Tabix index"
      },
      gbz: {
        description: "The chromosome's GBZ graph"
      },
      gfa: {
        description: "The chromosome's GFA (gzipped)"
      },
      graph_stats: {
        description: "graph-stats.tsv: nodes, edges, length"
      },
      path_stats: {
        description: "<cohort_id>.pangenome.<ref_name>.<chromosome>.path-stats.tsv.gz: cactus-graphmap-join's path-stats (ref_chrom, path, length), one row per path fragment of every haplotype in the graph"
      },
      refgaps: {
        description: "<cohort_id>.pangenome.<ref_name>.<chromosome>.refgaps.bed.gz: the join's reference intervals the graph does not cover"
      },
      genome_bp: {
        description: "<cohort_id>.pangenome.<ref_name>.<chromosome>.genome-bp.tsv: paths and bp per genome (sample#haplotype) in the chromosome's graph"
      },
      degree_stats: {
        description: "<cohort_id>.pangenome.<ref_name>.<chromosome>.degree.tsv: vg stats -D of the GBZ (node degree histogram)"
      },
      logs: {
        description: "The four Cactus logs, tarred"
      }
    }
  }

  parameter_meta {
    cohort_id: {
      description: "Cohort ID; output prefix"
    }
    ref_name: {
      description: "Reference name; part of the output names"
    }
    chromosome: {
      description: "The chromosome to build (its tar is picked by name)"
    }
    chrom_tars: {
      description: "The split step's tars"
    }
    chromfile: {
      description: "The split step's chromfile"
    }
    split_root: {
      description: "The split step's split/ path, rewritten to this task's"
    }
    reference: {
      description: "The reference genome's name in the graph"
    }
    cpu: {
      description: "Cores (--maxCores, --mgCores, --mapCores, --consCores, --indexCores)"
    }
    mem_gb: {
      description: "Memory (GB, --maxMemory)"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String cohort_id
    String ref_name
    String chromosome
    Array[File] chrom_tars
    File chromfile
    String split_root
    String reference = "GRCh38"
    Int cpu = 16
    Int mem_gb = 64
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{cohort_id}.pangenome.~{ref_name}.~{chromosome}"
  Int disk_size = ceil(size(chrom_tars, "GB") * 2 + 100)

  command <<<
    set -euo pipefail
    export HOME="$PWD" TMPDIR="$PWD/tmp"
    mkdir -p tmp work split
    cactus-pangenome --version 2>&1 | tail -1 >&2
    tar_path=""
    for t in ~{sep=" " chrom_tars}; do
      if [ "$(basename "$t")" = "~{chromosome}.tar" ]; then tar_path=$t; fi
    done
    [ -n "$tar_path" ] || { echo "no tar for ~{chromosome} among the split step's outputs" >&2; exit 1; }
    tar -C split -xf "$tar_path"
    grep -P "^~{chromosome}\t" "~{chromfile}" | sed "s#~{split_root}#$PWD/split#g" > chromfile.txt
    [ -s chromfile.txt ] || { echo "~{chromosome} is not in the chromfile" >&2; exit 1; }
    sed -i "s#~{split_root}#$PWD/split#g" "split/seqfiles/~{chromosome}.seqfile"
    cat chromfile.txt "split/seqfiles/~{chromosome}.seqfile" >&2
    common=(--reference ~{reference} --maxCores ~{cpu} --maxMemory ~{mem_gb}G --workDir "$PWD/work")

    cactus-minigraph "$PWD/js4" chromfile.txt construct --batch --mgCores ~{cpu} "${common[@]}" --logFile step4.minigraph.log
    cactus-graphmap "$PWD/js5" construct/chromfile.mg.txt map --batch --mapCores ~{cpu} "${common[@]}" --logFile step5.graphmap.log
    cactus-align "$PWD/js6" map/chromfile.gm.txt align --batch --pangenome --outVG --consCores ~{cpu} "${common[@]}" --logFile step6.align.log
    cactus-graphmap-join "$PWD/js7" --vg align/*.vg --hal align/*.hal --sv-gfa construct/*.gfa.gz \
      --outDir join --outName "~{chromosome}" --vcf --vcfwave --gbz --gfa --indexCores ~{cpu} "${common[@]}" --logFile step7.join.log

    mv "join/~{chromosome}.vcf.gz" "~{out_prefix}.vcf.gz"; mv "join/~{chromosome}.vcf.gz.tbi" "~{out_prefix}.vcf.gz.tbi"
    mv "join/~{chromosome}.wave.vcf.gz" "~{out_prefix}.wave.vcf.gz"; mv "join/~{chromosome}.wave.vcf.gz.tbi" "~{out_prefix}.wave.vcf.gz.tbi"
    mv "join/~{chromosome}.gbz" "~{out_prefix}.gbz"
    mv "join/~{chromosome}.gfa.gz" "~{out_prefix}.gfa.gz"
    cp "join/~{chromosome}.stats/graph-stats.tsv" "~{out_prefix}.graph-stats.tsv"
    # the join's other statistics (deleted with join/ otherwise): every path fragment's length and the reference gaps
    cp "join/~{chromosome}.stats/path-stats.tsv.gz" "~{out_prefix}.path-stats.tsv.gz"
    cp "join/~{chromosome}.stats/refgaps.bed.gz" "~{out_prefix}.refgaps.bed.gz"
    # bp per genome (PanSN sample#haplotype) over its path fragments, and the node degree histogram of the GBZ
    zcat "~{out_prefix}.path-stats.tsv.gz" | awk -F'\t' 'NR > 1 && NF >= 3 { n = split($2, a, "#"); g = (n >= 2) ? a[1] "#" a[2] : $2; paths[g]++; bp[g] += $3 }
      END { print "genome\tpaths\tbp"; for (g in paths) print g "\t" paths[g] "\t" bp[g] }' | { read -r h; echo "$h"; sort; } > "~{out_prefix}.genome-bp.tsv"
    vg stats -D "~{out_prefix}.gbz" > "~{out_prefix}.degree.tsv"
    tar -czf "~{out_prefix}.cactus-logs.tar.gz" step4.minigraph.log step5.graphmap.log step6.align.log step7.join.log
    cat "~{out_prefix}.graph-stats.tsv" "~{out_prefix}.genome-bp.tsv" >&2
    rm -rf work tmp split construct map align join js4 js5 js6 js7
  >>>

  output {
    File vcf = "~{out_prefix}.vcf.gz"
    File vcf_index = "~{out_prefix}.vcf.gz.tbi"
    File wave_vcf = "~{out_prefix}.wave.vcf.gz"
    File wave_vcf_index = "~{out_prefix}.wave.vcf.gz.tbi"
    File gbz = "~{out_prefix}.gbz"
    File gfa = "~{out_prefix}.gfa.gz"
    File graph_stats = "~{out_prefix}.graph-stats.tsv"
    File path_stats = "~{out_prefix}.path-stats.tsv.gz"
    File refgaps = "~{out_prefix}.refgaps.bed.gz"
    File genome_bp = "~{out_prefix}.genome-bp.tsv"
    File degree_stats = "~{out_prefix}.degree.tsv"
    File logs = "~{out_prefix}.cactus-logs.tar.gz"
  }

  runtime {
    docker: "quay.io/comparative-genomics-toolkit/cactus@sha256:fb4f1bab95d9c3b003d1036989f977e8cde2ffcf90a27691b992bd7d2d1a253d"  # v3.3.0
    cpu: cpu
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
