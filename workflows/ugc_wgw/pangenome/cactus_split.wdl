version 1.0
# ugc-pacbio-wgw — Minigraph-Cactus step 3 genome-wide: the members' alignments to the reference-only graph gathered, contigs split by chromosome
# origin: none (cactus v3.3.0's own commands, doc/pangenome.md "--mgSplit"; verified on the dev machine 2026-10-20; the mapping moved to cactus_map.wdl on 2026-11-13)
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8.2
import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_cactus_split {
  meta {
    description: "cactus-graphmap-split --mgSplit --refContigsFile --otherContig on the alignments the per-member mapping tasks made (ugc_wgw_cactus_map): the members' PAFs are concatenated in the seqfile's order, the reference's own alignment to the graph (one copy in every member's PAF) kept once, and the seqfile is written as one cactus-graphmap over every haplotype would have left it (the star tree over the reference, the haplotypes and _MINIGRAPH_, and the graph's node sequences). The result equals the one-task route's: contig assignments, per-chromosome FASTAs and PAF (checked on the dev machine). Single-machine inside the full image; Toil's work directory, job store, temporaries and coordination (lock) directory on the node-local scratch the backend names through TMPDIR, else in the task directory. Every chromosome directory and its seqfile go into one tar; the chromfile, the split log, contig_sizes.tsv and the per-member assigned/ambiguous contig counts come out beside them. Ambiguous contigs (the config's coverage fractions by length, best/second ratio, MAPQ) are left out of every graph."
    outputs: {
      chrom_tars: {
        description: "<chromosome>.tar per reference contig with data (and the other contig): split/<chr>/ and split/seqfiles/<chr>.seqfile"
      },
      chromfile: {
        description: "chromfile.txt: chromosome, seqfile path, PAF path (absolute paths inside this task's container)"
      },
      split_root: {
        description: "The absolute path of split/ in this task's container, the prefix the chromosome tasks rewrite"
      },
      chromosomes_found: {
        description: "The chromosomes the split produced, from the chromfile"
      },
      split_log: {
        description: "minigraph.split.log"
      },
      cactus_logs: {
        description: "The Toil logs, tarred: step1.minigraph (the reference graph), every member's graphmap log, step3.split"
      },
      contig_sizes: {
        description: "contig_sizes.tsv: bases per genome and chromosome"
      },
      ambiguous_contigs: {
        description: "genome, assigned, ambiguous: contig counts per genome from the split log"
      }
    }
  }

  parameter_meta {
    ref_name: {
      description: "The reference map's name (output file names)"
    }
    cohort_id: {
      description: "Cohort ID"
    }
    ref_fasta: {
      description: "The assembly reference FASTA (the graph's reference, named GRCh38 in the seqfile)"
    }
    genome_names: {
      description: "<safe>.1, <safe>.2 per member, in the order of hap_fastas"
    }
    hap_fastas: {
      description: "The haplotype FASTAs, hap1 then hap2 of each member, cohort order"
    }
    ref_graph: {
      description: "The reference-only graph (ugc_wgw_cactus_refgraph)"
    }
    graph_fasta: {
      description: "The graph's node sequences (any member's ugc_wgw_cactus_map.graph_fasta)"
    }
    pafs: {
      description: "Every member's alignments (ugc_wgw_cactus_map.paf), cohort order"
    }
    step_logs: {
      description: "The reference graph's and the members' Toil logs, for the log tar"
    }
    ref_contigs: {
      description: "The reference contigs to split by, one per line (the map's asm_graph_chromosomes)"
    }
    other_contig: {
      description: "Name of the contig that collects the rest of the reference (chrM, unplaced)"
    }
    reference: {
      description: "The reference genome's name in the graph"
    }
    cpu: {
      description: "Cores (--maxCores)"
    }
    mem_gb: {
      description: "Memory (GB, --maxMemory)"
    }
    disk_gb: {
      description: "Disk (GB); default from the input sizes: the tars hold every member's FASTA again"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String cohort_id
    String ref_name
    File ref_fasta
    Array[String] genome_names
    Array[File] hap_fastas
    File ref_contigs
    String other_contig = "chrOther"
    File ref_graph
    File graph_fasta
    Array[File] pafs
    Array[File] step_logs
    String reference = "GRCh38"
    Int cpu = 16
    Int mem_gb = 64
    Int? disk_gb
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{cohort_id}.pangenome.~{ref_name}"
  Int disk_size = select_first([disk_gb, ceil(size(hap_fastas, "GB") * 3 + size(ref_fasta, "GB") * 2 + 50)])

  command <<<
    set -euo pipefail
    # Toil, Cactus's engine, keeps lock files in a coordination directory and refuses a network file
    # system for it ("file locks are known to break ceph filesystems", seen on the cluster 2026-10-22);
    # its job stores and work directory are thousands of small files. When the backend binds node-local
    # scratch it sets TMPDIR to it (TASK_SCRATCH_ENV, docs/guide/04-install.md): everything of Toil's
    # goes there and is removed with the job. Without scratch (the dev machine) it stays in the task
    # directory, which must then be a local file system.
    if [ -n "${TMPDIR:-}" ] && [ -d "$TMPDIR" ] && [ -w "$TMPDIR" ]; then toil="$TMPDIR/cactus"; else toil="$PWD/toil"; fi
    mkdir -p "$toil/work" "$toil/tmp" "$toil/coord" tars
    export HOME="$PWD" TMPDIR="$toil/tmp" TOIL_COORDINATION_DIR="$toil/coord"
    echo "toil directories under $toil" >&2
    cactus-pangenome --version 2>&1 | tail -1 >&2
    { printf '~{reference}\t%s\n' "~{ref_fasta}"; paste '~{write_lines(genome_names)}' '~{write_lines(hap_fastas)}'; } > seqfile.txt
    cat seqfile.txt >&2
    common=(--reference ~{reference} --maxCores ~{cpu} --maxMemory ~{mem_gb}G --workDir "$toil/work")

    # the seqfile as one cactus-graphmap over every haplotype leaves it: the star tree and the graph's own sequences
    ln -s "~{ref_graph}" ref.sv.gfa.gz
    ln -s "~{graph_fasta}" ref.sv.gfa.fa.gz
    tree="(~{reference}:1.0,$(sed 's/$/:1.0/' '~{write_lines(genome_names)}' | paste -sd, -),_MINIGRAPH_:1.0);"
    { cat seqfile.txt; echo "$tree"; printf '_MINIGRAPH_\tfile://%s/ref.sv.gfa.fa.gz\n' "$PWD"; } > seqfile.gm.txt
    # every member's alignments in the seqfile's order; each PAF also holds the reference's own alignment to the
    # graph, which one run over all haplotypes writes once
    n=0
    while read -r paf; do
      n=$((n + 1))
      if [ "$n" -eq 1 ]; then cat "$paf"; else awk -F'\t' -v p="id=~{reference}|" 'index($1, p) != 1' "$paf"; fi
    done < '~{write_lines(pafs)}' > all.paf
    echo "gathered $n PAFs, $(wc -l < all.paf) alignments" >&2
    cactus-graphmap-split "$toil/js3" seqfile.gm.txt ref.sv.gfa.gz all.paf --outDir split --mgSplit \
      --refContigsFile "~{ref_contigs}" --otherContig ~{other_contig} "${common[@]}" --logFile step3.split.log

    echo "$PWD/split" > split_root.txt
    cut -f1 split/chromfile.txt > chromosomes_found.txt
    while read -r chr; do
      [ -d "split/$chr" ] || continue
      tar -C split -cf "tars/$chr.tar" "$chr" "seqfiles/$chr.seqfile"
    done < chromosomes_found.txt
    # per genome: contigs (and bp) assigned to a chromosome, and contigs left out as ambiguous
    awk -F'\t' '{ n = split($0, w, " "); id = ""; len = 0
        for (i = 1; i <= n; i++) { if (w[i] ~ /^id=/) { id = substr(w[i], 4); sub(/\|.*/, "", id) } if (w[i] ~ /^len=/) len = substr(w[i], 5) + 0 }
        if (id == "" || $0 ~ /ref-contig/) next; seen[id]
        if ($0 ~ /^Assigned contig to /) { as[id]++; asbp[id] += len } else if ($0 ~ /ambiguous/) { am[id]++; ambp[id] += len } }
      END { print "genome\tassigned\tassigned_bp\tambiguous\tambiguous_bp"
        for (g in seen) print g "\t" (as[g]+0) "\t" (asbp[g]+0) "\t" (am[g]+0) "\t" (ambp[g]+0) }' \
      split/minigraph.split.log | { read -r h; echo "$h"; sort; } > "~{out_prefix}.ambiguous_contigs.tsv"
    cp split/minigraph.split.log "~{out_prefix}.split.log"
    mkdir -p logs
    while read -r l; do cp "$l" logs/; done < '~{write_lines(step_logs)}'
    cp step3.split.log logs/
    tar -czf "~{out_prefix}.split.cactus-logs.tar.gz" -C logs .
    cp split/contig_sizes.tsv "~{out_prefix}.contig_sizes.tsv"
    cat "~{out_prefix}.ambiguous_contigs.tsv" >&2
    rm -rf "$toil"
  >>>

  output {
    Array[File] chrom_tars = glob("tars/*.tar")
    File chromfile = "split/chromfile.txt"
    String split_root = read_string("split_root.txt")
    Array[String] chromosomes_found = read_lines("chromosomes_found.txt")
    File split_log = "~{out_prefix}.split.log"
    File cactus_logs = "~{out_prefix}.split.cactus-logs.tar.gz"
    File contig_sizes = "~{out_prefix}.contig_sizes.tsv"
    File ambiguous_contigs = "~{out_prefix}.ambiguous_contigs.tsv"
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
