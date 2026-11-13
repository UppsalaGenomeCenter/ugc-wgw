version 1.0
# ugc-pacbio-wgw — Minigraph-Cactus steps 1 and 2 in shards: the reference-only graph once, then one member's two haplotypes mapped to it per task
# origin: none (cactus v3.3.0's own commands, doc/pangenome.md "--mgSplit"; the per-member invocation verified on the dev machine 2026-11-13)
# see docs/DESIGN.md §5, §13, §17; docs/ENTRYPOINTS.md §8.2
import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_cactus_refgraph {
  meta {
    description: "cactus-minigraph --refOnly: the minigraph graph of the reference alone, the target every haplotype is mapped to before the contigs are split by chromosome. It depends on the reference only, so the call cache serves every cohort of a project from one run. Cactus refuses a seqfile with fewer than two leaf genomes; the second line names the reference again under another name, and --refOnly never reads it into the graph."
    outputs: {
      graph: {
        description: "ref.sv.gfa.gz: the reference-only minigraph graph"
      },
      log: {
        description: "The step's Toil log"
      }
    }
  }

  parameter_meta {
    ref_fasta: {
      description: "The assembly reference FASTA (the graph's reference)"
    }
    reference: {
      description: "The reference genome's name in the graph"
    }
    cpu: {
      description: "Cores (--maxCores, --mgCores)"
    }
    mem_gb: {
      description: "Memory (GB, --maxMemory); 3.3 GB measured for GRCh38"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File ref_fasta
    String reference = "GRCh38"
    Int cpu = 8
    Int mem_gb = 16
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(ref_fasta, "GB") * 4 + 20)

  command <<<
    set -euo pipefail
    # Toil's directories on the node-local scratch the backend names through TMPDIR, else in the task directory
    # (ugc_wgw_cactus_split has the reasons)
    if [ -n "${TMPDIR:-}" ] && [ -d "$TMPDIR" ] && [ -w "$TMPDIR" ]; then toil="$TMPDIR/cactus"; else toil="$PWD/toil"; fi
    mkdir -p "$toil/work" "$toil/tmp" "$toil/coord"
    export HOME="$PWD" TMPDIR="$toil/tmp" TOIL_COORDINATION_DIR="$toil/coord"
    cactus-pangenome --version 2>&1 | tail -1 >&2
    printf '~{reference}\t%s\nREFONLY\t%s\n' "~{ref_fasta}" "~{ref_fasta}" > seqfile.txt
    cactus-minigraph "$toil/js1" seqfile.txt ref.sv.gfa.gz --refOnly --reference ~{reference} --mgCores ~{cpu} \
      --maxCores ~{cpu} --maxMemory ~{mem_gb}G --workDir "$toil/work" --logFile step1.minigraph.log
    rm -rf "$toil"
  >>>

  output {
    File graph = "ref.sv.gfa.gz"
    File log = "step1.minigraph.log"
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

task ugc_wgw_cactus_map {
  meta {
    description: "cactus-graphmap --mgSplit of one member's two haplotypes against the reference-only graph: the seqfile names the reference and the two haplotypes, the mappings run one after the other (--mapCores at the task's cores, so the memory is one mapping's). A haplotype's alignment depends on that haplotype and the reference graph alone, so the task's inputs are the member's and the call cache serves every cohort the member is in: a larger cohort nested around a smaller one maps only its new members. The PAF also holds the reference's own alignment to the graph; the split keeps one copy."
    outputs: {
      paf: {
        description: "<member>.graphmap.paf: both haplotypes' contigs (id=<name>|<contig>) against the reference graph, and the reference's own line"
      },
      graph_fasta: {
        description: "ref.sv.gfa.fa.gz: the graph's node sequences as cactus-graphmap --outputFasta writes them (the same file from every member; the split takes one)"
      },
      log: {
        description: "<member>.graphmap.log: the step's Toil log"
      }
    }
  }

  parameter_meta {
    member: {
      description: "The member's Cactus-safe name (file names)"
    }
    genome_names: {
      description: "<safe>.1 and <safe>.2, the two haplotypes' genome names"
    }
    hap_fastas: {
      description: "The member's two haplotype FASTAs, in the order of genome_names"
    }
    ref_fasta: {
      description: "The assembly reference FASTA"
    }
    ref_graph: {
      description: "The reference-only graph (ugc_wgw_cactus_refgraph)"
    }
    reference: {
      description: "The reference genome's name in the graph"
    }
    cpu: {
      description: "Cores (--maxCores, --mapCores): one mapping at a time; minigraph uses few cores of one mapping (15-22 % CPU efficiency measured at 48), so 8 is plenty"
    }
    mem_gb: {
      description: "Memory (GB, --maxMemory); default 60 plus 15 times the two gzipped FASTAs together: 55 GB measured for a chr20-only haplotype (the reference graph's index), 65-80 GB for a whole human haplotype"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String member
    Array[String] genome_names
    Array[File] hap_fastas
    File ref_fasta
    File ref_graph
    String reference = "GRCh38"
    Int cpu = 8
    Int? mem_gb
    RuntimeAttributes runtime_attributes
  }

  Int mem = select_first([mem_gb, ceil(size(hap_fastas, "GB") * 15 + 60)])
  Int disk_size = ceil(size(hap_fastas, "GB") * 4 + size(ref_fasta, "GB") * 2 + 30)

  command <<<
    set -euo pipefail
    if [ -n "${TMPDIR:-}" ] && [ -d "$TMPDIR" ] && [ -w "$TMPDIR" ]; then toil="$TMPDIR/cactus"; else toil="$PWD/toil"; fi
    mkdir -p "$toil/work" "$toil/tmp" "$toil/coord"
    export HOME="$PWD" TMPDIR="$toil/tmp" TOIL_COORDINATION_DIR="$toil/coord"
    cactus-pangenome --version 2>&1 | tail -1 >&2
    { printf '~{reference}\t%s\n' "~{ref_fasta}"; paste '~{write_lines(genome_names)}' '~{write_lines(hap_fastas)}'; } > seqfile.txt
    cat seqfile.txt >&2
    # cactus-graphmap edits the seqfile it is given (the tree and the graph's own sequence): it gets a copy
    cp seqfile.txt seqfile.gm.txt
    ln -s "~{ref_graph}" ref.sv.gfa.gz
    cactus-graphmap "$toil/js2" seqfile.gm.txt ref.sv.gfa.gz "~{member}.graphmap.paf" --mgSplit --outputFasta ref.sv.gfa.fa.gz \
      --reference ~{reference} --mapCores ~{cpu} --maxCores ~{cpu} --maxMemory ~{mem}G --workDir "$toil/work" \
      --logFile "~{member}.graphmap.log"
    wc -l "~{member}.graphmap.paf" >&2
    rm -rf "$toil"
  >>>

  output {
    File paf = "~{member}.graphmap.paf"
    File graph_fasta = "ref.sv.gfa.fa.gz"
    File log = "~{member}.graphmap.log"
  }

  runtime {
    docker: "quay.io/comparative-genomics-toolkit/cactus@sha256:fb4f1bab95d9c3b003d1036989f977e8cde2ffcf90a27691b992bd7d2d1a253d"  # v3.3.0
    cpu: cpu
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
