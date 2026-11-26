version 1.0
# ugc-pacbio-wgw — Minigraph-Cactus steps 4-7 for one chromosome, one task each: construct, map, align, join (VCF, vcfwave VCF, GBZ, GFA)
# origin: none (cactus v3.3.0's own batch commands, doc/pangenome.md "--mgSplit"; verified on the dev machine 2026-10-20; one task per command since 0.26.0)
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8.2
import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

# The four tasks are the four commands of Cactus's per-chromosome recipe, cut where one command's output directory is the
# next one's input, so that each has its own request and the call cache is a checkpoint between them (a cohort of 200
# chr20 haplotypes: 10 h of construct on one or two cores, 5 h of mappings, 2 h of align at the memory peak, 2 h of join).
# Cactus writes absolute paths into its chromfiles and seqfiles: a task that hands one on replaces its own working
# directory by the token @UGC_WGW_WORK@, and the task that reads it puts its own directory back.

task ugc_wgw_cactus_chrom_construct {
  meta {
    description: "Step 4 for one chromosome: cactus-minigraph construct --batch on the split step's tar, the absolute paths of its seqfile rewritten to this task. minigraph adds the haplotypes to the SV graph one after another (the reference first, the rest by mash distance), on one or two cores whatever is offered, so the task's time is the number of haplotypes in the construct times a few minutes. With construct_genomes naming fewer genomes than the chromosome has, only those (and the reference) build the SV graph; the chromfile this task hands on names every genome again, so that all of them are mapped, aligned and joined. When none of the named genomes has sequence on the chromosome (chrY and a construct of XX members), every genome of the chromosome builds it."
    outputs: {
      sv_gfa: {
        description: "<cohort_id>.pangenome.<ref_name>.<chromosome>.sv.gfa.gz: minigraph's SV graph of the chromosome"
      },
      chromfile_mg: {
        description: "The chromfile for the map step (chromosome, seqfile of every genome, SV graph), this task's directory replaced by @UGC_WGW_WORK@"
      },
      construct_list: {
        description: "<cohort_id>.pangenome.<ref_name>.<chromosome>.construct_genomes.txt: the genomes (haplotypes) whose sequence built the SV graph, one per line"
      },
      log: {
        description: "step4.minigraph.log"
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
    construct_genomes: {
      description: "Genome names (<safe>.1, <safe>.2) whose haplotypes build the SV graph; empty: every genome of the chromosome"
    }
    reference: {
      description: "The reference genome's name in the graph"
    }
    cpu: {
      description: "Cores (--maxCores, --mgCores); minigraph's construction keeps one or two busy"
    }
    mem_gb: {
      description: "Memory (GB, --maxMemory); 3.8 GB measured for 200 chr20 haplotypes"
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
    Array[String] construct_genomes = []
    String reference = "GRCh38"
    Int cpu = 4
    Int mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{cohort_id}.pangenome.~{ref_name}.~{chromosome}"
  Int disk_size = ceil(size(chrom_tars, "GB") * 2 + 50)

  command <<<
    set -euo pipefail
    # Toil's directories as in ugc_wgw_cactus_split: on the node-local scratch the backend names through
    # TMPDIR (lock files there are what Toil refuses on a network file system), else in the task directory
    if [ -n "${TMPDIR:-}" ] && [ -d "$TMPDIR" ] && [ -w "$TMPDIR" ]; then toil="$TMPDIR/cactus"; else toil="$PWD/toil"; fi
    mkdir -p "$toil/work" "$toil/tmp" "$toil/coord" split
    export HOME="$PWD" TMPDIR="$toil/tmp" TOIL_COORDINATION_DIR="$toil/coord"
    echo "toil directories under $toil" >&2
    cactus-pangenome --version 2>&1 | tail -1 >&2
    tar_path=""
    for t in ~{sep=" " chrom_tars}; do
      if [ "$(basename "$t")" = "~{chromosome}.tar" ]; then tar_path=$t; fi
    done
    [ -n "$tar_path" ] || { echo "no tar for ~{chromosome} among the split step's outputs" >&2; exit 1; }
    tar -C split -xf "$tar_path"
    grep -P "^~{chromosome}\t" "~{chromfile}" | sed "s#~{split_root}#$PWD/split#g" > chromfile.txt
    [ -s chromfile.txt ] || { echo "~{chromosome} is not in the chromfile" >&2; exit 1; }
    seqfile="split/seqfiles/~{chromosome}.seqfile"
    sed -i "s#~{split_root}#$PWD/split#g" "$seqfile"
    cat chromfile.txt "$seqfile" >&2

    # the genomes of the construct: every haplotype line of the chromosome's seqfile, or the named ones that are in it
    genome_lines() { awk -F'\t' -v ref="~{reference}" 'NF >= 2 && $1 != ref && $1 != "_MINIGRAPH_" && $0 !~ /^\(/' "$1"; }
    names='~{write_lines(construct_genomes)}'
    construct_chromfile=chromfile.txt
    genome_lines "$seqfile" | cut -f1 > "~{out_prefix}.construct_genomes.txt"
    n_all=$(wc -l < "~{out_prefix}.construct_genomes.txt")
    if [ -s "$names" ]; then
      subset="split/seqfiles/~{chromosome}.construct.seqfile"
      awk -F'\t' -v ref="~{reference}" -v names="$names" 'BEGIN { while ((getline l < names) > 0) keep[l] }
        $0 ~ /^\(/ { next }
        NF < 2 || $1 == ref || $1 == "_MINIGRAPH_" || ($1 in keep)' "$seqfile" > "$subset"
      n_in=$(genome_lines "$subset" | wc -l)
      if [ "$n_in" -eq 0 ]; then
        echo "construct: none of the named genomes has sequence on ~{chromosome}; all $n_all genomes of the chromosome build its SV graph" >&2
      elif [ "$n_in" -lt "$n_all" ]; then
        awk -F'\t' -v OFS='\t' -v s="$subset" '{ $2 = s; print }' chromfile.txt > chromfile.construct.txt
        construct_chromfile=chromfile.construct.txt
        genome_lines "$subset" | cut -f1 > "~{out_prefix}.construct_genomes.txt"
        echo "construct: $n_in of the $n_all genomes of ~{chromosome} build its SV graph" >&2
        cat "$subset" >&2
      fi
    fi
    common=(--reference ~{reference} --maxCores ~{cpu} --maxMemory ~{mem_gb}G --workDir "$toil/work")

    cactus-minigraph "$toil/js4" "$construct_chromfile" construct --batch --mgCores ~{cpu} "${common[@]}" --logFile step4.minigraph.log

    # every genome is mapped to the SV graph, in or out of its construction: the chromfile names the whole seqfile again
    awk -F'\t' -v OFS='\t' -v s="$(cut -f2 chromfile.txt)" '{ $2 = s; print }' construct/chromfile.mg.txt \
      | sed "s#$PWD#@UGC_WGW_WORK@#g" > "~{out_prefix}.chromfile.mg.txt"
    mv "construct/~{chromosome}.sv.gfa.gz" "~{out_prefix}.sv.gfa.gz"
    cat "~{out_prefix}.chromfile.mg.txt" >&2
    rm -rf "$toil" split construct
  >>>

  output {
    File sv_gfa = "~{out_prefix}.sv.gfa.gz"
    File chromfile_mg = "~{out_prefix}.chromfile.mg.txt"
    File construct_list = "~{out_prefix}.construct_genomes.txt"
    File log = "step4.minigraph.log"
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

task ugc_wgw_cactus_chrom_map {
  meta {
    description: "Step 5 for one chromosome: cactus-graphmap map --batch, every genome of the chromosome's seqfile mapped to the SV graph with minigraph, cpu / map_cores mappings at once, then Cactus's filters and the merged PAF. The mappings are independent of one another: more cores are more mappings at a time (200 chr20 haplotypes: 287 s each at 4 threads, 4 h 45 min four at a time)."
    outputs: {
      map_tar: {
        description: "<cohort_id>.pangenome.<ref_name>.<chromosome>.map.tar: the step's map/ directory (PAF, GAF, the graph's FASTA, the seqfile and chromfile for the align step with this task's directory replaced by @UGC_WGW_WORK@)"
      },
      log: {
        description: "step5.graphmap.log"
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
      description: "The chromosome (its tar is picked by name)"
    }
    chrom_tars: {
      description: "The split step's tars (the genomes' FASTAs of the chromosome)"
    }
    split_root: {
      description: "The split step's split/ path, rewritten to this task's"
    }
    sv_gfa: {
      description: "The construct step's SV graph"
    }
    chromfile_mg: {
      description: "The construct step's chromfile"
    }
    reference: {
      description: "The reference genome's name in the graph"
    }
    cpu: {
      description: "Cores (--maxCores)"
    }
    map_cores: {
      description: "Cores per minigraph mapping of the chromosome's haplotypes (--mapCores): cpu / map_cores mappings run at once. Default 4: a chromosome's mappings are small, and minigraph uses few cores on one of them (chr21 of five samples: 280 minutes of mapping one after another, 73 for the slowest haplotype, 15 % CPU efficiency at 16 cores each)"
    }
    mem_gb: {
      description: "Memory (GB, --maxMemory); default 16 plus 10 per mapping at once (3 GB measured per mapping, 9 GiB asked by Toil)"
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
    String split_root
    File sv_gfa
    File chromfile_mg
    String reference = "GRCh38"
    Int cpu = 16
    Int map_cores = 4
    Int? mem_gb
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{cohort_id}.pangenome.~{ref_name}.~{chromosome}"
  Int mem = select_first([mem_gb, 10 * (cpu / map_cores) + 16])
  Int disk_size = ceil(size(chrom_tars, "GB") * 2 + size(sv_gfa, "GB") * 4 + 50)

  command <<<
    set -euo pipefail
    if [ -n "${TMPDIR:-}" ] && [ -d "$TMPDIR" ] && [ -w "$TMPDIR" ]; then toil="$TMPDIR/cactus"; else toil="$PWD/toil"; fi
    mkdir -p "$toil/work" "$toil/tmp" "$toil/coord" split construct
    export HOME="$PWD" TMPDIR="$toil/tmp" TOIL_COORDINATION_DIR="$toil/coord"
    echo "toil directories under $toil" >&2
    cactus-pangenome --version 2>&1 | tail -1 >&2
    tar_path=""
    for t in ~{sep=" " chrom_tars}; do
      if [ "$(basename "$t")" = "~{chromosome}.tar" ]; then tar_path=$t; fi
    done
    [ -n "$tar_path" ] || { echo "no tar for ~{chromosome} among the split step's outputs" >&2; exit 1; }
    tar -C split -xf "$tar_path"
    sed -i "s#~{split_root}#$PWD/split#g" "split/seqfiles/~{chromosome}.seqfile"
    cp "~{sv_gfa}" "construct/~{chromosome}.sv.gfa.gz"
    sed "s#@UGC_WGW_WORK@#$PWD#g" "~{chromfile_mg}" > construct/chromfile.mg.txt
    cat construct/chromfile.mg.txt >&2
    common=(--reference ~{reference} --maxCores ~{cpu} --maxMemory ~{mem}G --workDir "$toil/work")

    cactus-graphmap "$toil/js5" construct/chromfile.mg.txt map --batch --mapCores ~{map_cores} "${common[@]}" --logFile step5.graphmap.log

    sed -i "s#$PWD#@UGC_WGW_WORK@#g" map/chromfile.gm.txt map/*.seqfile
    cat map/chromfile.gm.txt >&2
    tar -cf "~{out_prefix}.map.tar" map
    rm -rf "$toil" split construct map
  >>>

  output {
    File map_tar = "~{out_prefix}.map.tar"
    File log = "step5.graphmap.log"
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

task ugc_wgw_cactus_chrom_align {
  meta {
    description: "Step 6 for one chromosome: cactus-align --pangenome --outVG --batch on the map step's alignments. cactus_consolidated is the chromosome graph's memory peak (33 GiB for 40 chr20 haplotypes, 55 GiB for 200: about 28 GiB plus 2.3 GiB per Gb of haplotype sequence); Cactus's own estimate is several times that and is cut to --maxMemory."
    outputs: {
      hal: {
        description: "<cohort_id>.pangenome.<ref_name>.<chromosome>.hal: the chromosome's alignment"
      },
      vg: {
        description: "<cohort_id>.pangenome.<ref_name>.<chromosome>.raw.vg: the chromosome's graph before the join's clipping and indexing"
      },
      log: {
        description: "step6.align.log"
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
      description: "The chromosome (its tar is picked by name)"
    }
    chrom_tars: {
      description: "The split step's tars (the genomes' FASTAs of the chromosome)"
    }
    split_root: {
      description: "The split step's split/ path, rewritten to this task's"
    }
    map_tar: {
      description: "The map step's map/ directory"
    }
    reference: {
      description: "The reference genome's name in the graph"
    }
    cpu: {
      description: "Cores (--maxCores, --consCores)"
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
    String split_root
    File map_tar
    String reference = "GRCh38"
    Int cpu = 16
    Int mem_gb = 64
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{cohort_id}.pangenome.~{ref_name}.~{chromosome}"
  Int disk_size = ceil(size(chrom_tars, "GB") * 6 + size(map_tar, "GB") * 4 + 100)

  command <<<
    set -euo pipefail
    if [ -n "${TMPDIR:-}" ] && [ -d "$TMPDIR" ] && [ -w "$TMPDIR" ]; then toil="$TMPDIR/cactus"; else toil="$PWD/toil"; fi
    mkdir -p "$toil/work" "$toil/tmp" "$toil/coord" split
    export HOME="$PWD" TMPDIR="$toil/tmp" TOIL_COORDINATION_DIR="$toil/coord"
    echo "toil directories under $toil" >&2
    cactus-pangenome --version 2>&1 | tail -1 >&2
    tar_path=""
    for t in ~{sep=" " chrom_tars}; do
      if [ "$(basename "$t")" = "~{chromosome}.tar" ]; then tar_path=$t; fi
    done
    [ -n "$tar_path" ] || { echo "no tar for ~{chromosome} among the split step's outputs" >&2; exit 1; }
    tar -C split -xf "$tar_path"
    sed -i "s#~{split_root}#$PWD/split#g" "split/seqfiles/~{chromosome}.seqfile"
    tar -xf "~{map_tar}"
    sed -i "s#@UGC_WGW_WORK@#$PWD#g" map/chromfile.gm.txt map/*.seqfile
    cat map/chromfile.gm.txt >&2
    common=(--reference ~{reference} --maxCores ~{cpu} --maxMemory ~{mem_gb}G --workDir "$toil/work")

    cactus-align "$toil/js6" map/chromfile.gm.txt align --batch --pangenome --outVG --consCores ~{cpu} "${common[@]}" --logFile step6.align.log

    mv align/*.hal "~{out_prefix}.hal"
    mv align/*.vg "~{out_prefix}.raw.vg"
    rm -rf "$toil" split map align
  >>>

  output {
    File hal = "~{out_prefix}.hal"
    File vg = "~{out_prefix}.raw.vg"
    File log = "step6.align.log"
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

task ugc_wgw_cactus_chrom_join {
  meta {
    description: "Step 7 for one chromosome: cactus-graphmap-join --vcf --vcfwave --gbz --gfa on the align step's graph and alignment with the construct step's SV graph, then the statistics the summaries read. vcfwave runs in chunks, and one chunk over a complex region can be most of the step (67 of 117 minutes for 200 chr20 haplotypes)."
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
        description: "The four Cactus logs of the chromosome (construct, map, align, join), tarred"
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
      description: "The chromosome; the join's output name"
    }
    vg: {
      description: "The align step's graph"
    }
    hal: {
      description: "The align step's alignment"
    }
    sv_gfa: {
      description: "The construct step's SV graph"
    }
    step_logs: {
      description: "The logs of the construct, map and align steps, for the chromosome's log tar"
    }
    reference: {
      description: "The reference genome's name in the graph"
    }
    cpu: {
      description: "Cores (--maxCores, --indexCores)"
    }
    mem_gb: {
      description: "Memory (GB, --maxMemory); 17 GiB measured for 200 chr20 haplotypes (vg convert)"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String cohort_id
    String ref_name
    String chromosome
    File vg
    File hal
    File sv_gfa
    Array[File] step_logs
    String reference = "GRCh38"
    Int cpu = 16
    Int mem_gb = 64
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{cohort_id}.pangenome.~{ref_name}.~{chromosome}"
  Int disk_size = ceil((size(vg, "GB") + size(hal, "GB")) * 8 + 100)

  command <<<
    set -euo pipefail
    if [ -n "${TMPDIR:-}" ] && [ -d "$TMPDIR" ] && [ -w "$TMPDIR" ]; then toil="$TMPDIR/cactus"; else toil="$PWD/toil"; fi
    mkdir -p "$toil/work" "$toil/tmp" "$toil/coord" align construct
    export HOME="$PWD" TMPDIR="$toil/tmp" TOIL_COORDINATION_DIR="$toil/coord"
    echo "toil directories under $toil" >&2
    cactus-pangenome --version 2>&1 | tail -1 >&2
    # under the names the one-task recipe gave them: the join takes the chromosome's name from them
    ln -s "~{vg}" "align/~{chromosome}.raw.vg"
    ln -s "~{hal}" "align/~{chromosome}.hal"
    cp "~{sv_gfa}" "construct/~{chromosome}.sv.gfa.gz"
    while read -r l; do cp "$l" .; done < '~{write_lines(step_logs)}'
    common=(--reference ~{reference} --maxCores ~{cpu} --maxMemory ~{mem_gb}G --workDir "$toil/work")

    cactus-graphmap-join "$toil/js7" --vg align/*.vg --hal align/*.hal --sv-gfa construct/*.gfa.gz \
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
    rm -rf "$toil" construct align join
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
