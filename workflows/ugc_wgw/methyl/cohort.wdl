version 1.0

# ugc-pacbio-wgw — the cohort's methylation layer: a MethBat background of the members, every member scored against it, group comparisons, joint segmentation, deconvolution, and the tables the driver reads
# origin: none (MethBat 1.1.0 from upstream's methbat image; the gene annotation after nallo's annotate_methylation, MIT, genomic-medicine-sweden/nallo)
# see docs/DESIGN.md §5, §8.1, §18 (2026-11-06); docs/ENTRYPOINTS.md §10

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_methbat_build {
  meta {
    description: "methbat build: the cohort's background profile from the members' MethBat profiles, with the labels the driver gives each member (ALL is implicit). The collection file is written here, where the profiles are mounted. With exclude_index >= 0 one member is left out (the leave-one-out backgrounds)."
    outputs: {
      background: {
        description: "<out_prefix>.background.tsv: per region and label (data_category) the dataset counts per summary label and the means and standard deviations of the combined methylation and of the haplotype delta"
      },
      collection: {
        description: "<out_prefix>.collection.tsv: identifier, filename, labels per member (the filename is the profile's name in this task)"
      }
    }
  }

  parameter_meta {
    sample_ids: {
      description: "Members, cohort order"
    }
    profiles: {
      description: "Their MethBat profiles (methbat_profile), aligned with sample_ids"
    }
    labels: {
      description: "Per member, aligned: labels ';'-joined (may be empty); ALL is added by MethBat"
    }
    exclude_index: {
      description: "Index of a member to leave out (-1: none)"
    }
    out_prefix: {
      description: "Output prefix"
    }
    mem_gb: {
      description: "Memory (GB)"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[String] sample_ids
    Array[File] profiles
    Array[String] labels
    Int exclude_index = -1
    String out_prefix
    Int mem_gb = 8
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(profiles, "GB") * 2 + 5)

  command <<<
    set -euo pipefail

    methbat --version

    mkdir -p in
    python3 - "~{write_lines(sample_ids)}" "~{write_lines(profiles)}" "~{write_lines(labels)}" ~{exclude_index} > "~{out_prefix}.collection.tsv" <<'EOF'
    import os, sys
    ids, files, labels = ([l.rstrip("\n") for l in open(p)] for p in sys.argv[1:4])
    exclude = int(sys.argv[4])
    print("identifier\tfilename\tlabels")
    for i, (sid, path, lab) in enumerate(zip(ids, files, labels)):
        if i == exclude:
            continue
        link = os.path.join("in", f"{i}_{os.path.basename(path)}")
        os.symlink(path, link)
        print(f"{sid}\t{os.path.abspath(link)}\t{lab}")
    EOF

    methbat build \
      --input-collection "~{out_prefix}.collection.tsv" \
      --output-profile "~{out_prefix}.background.tsv"
  >>>

  output {
    File background = "~{out_prefix}.background.tsv"
    File collection = "~{out_prefix}.collection.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/methbat@sha256:281569947c0a6154f9a6bdaa1bb5f1e67dd160ccb49028d3218a49a37a0488fb"  # 1.1.0_build1
    cpu: 1
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

task ugc_wgw_methbat_profile_bg {
  meta {
    description: "methbat profile of one member's 5mC pileup against the cohort background, compared under every label the member carries (ALL and its own), then annotated with the CpG-island gene map: the scored profile as MethBat writes it, the same with hgnc_id, gene_symbol and region_type appended, the member's events (rows whose compare_label is HypoMethylated, HyperMethylated, HypoASM or HyperASM) and its counts per label."
    outputs: {
      scored: {
        description: "<out_prefix>.methbat.scored.tsv: MethBat's single-dataset profile with the background columns (one row per region per label)"
      },
      scored_annotated: {
        description: "<out_prefix>.methbat.scored.annotated.tsv: the same rows with hgnc_id, gene_symbol, region_type appended"
      },
      events: {
        description: "<out_prefix>.methbat.events.tsv: sample_id, the scored row, the three gene columns; only the Hypo/Hyper rows"
      },
      counts: {
        description: "<out_prefix>.methbat.counts.tsv: sample_id, label, metric, value (regions, each compare_label, and for ALL each summary_label as summary_<label>)"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Member"
    }
    cpg_pileup_bed: {
      description: "Its 5mCpG pileup BED (methbat pileup)"
    }
    background: {
      description: "The cohort background (ugc_wgw_methbat_build; the leave-one-out one when that is on)"
    }
    labels: {
      description: "The member's labels ';'-joined (may be empty); ALL is always compared"
    }
    methbat_map: {
      description: "The CpG-island gene map (reference map key methbat_map): chrom, start, end, region_label, hgnc_id, gene_symbol, region_type"
    }
    zscore_min: {
      description: "methbat profile --min-abs-zscore"
    }
    methyl_delta_min: {
      description: "methbat profile --min-abs-delta (percentage points)"
    }
    background_min_datasets: {
      description: "methbat profile --min-datasets"
    }
    out_prefix: {
      description: "Output prefix (<sample_id>.<ref_name>)"
    }
    mem_gb: {
      description: "Memory (GB)"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    File cpg_pileup_bed
    File background
    String labels = ""
    File methbat_map
    Float zscore_min = 3.0
    Float methyl_delta_min = 20
    Int background_min_datasets = 10
    String out_prefix
    Int mem_gb = 4
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(cpg_pileup_bed, "GB") * 2 + size(background, "GB") * 3 + 20)

  command <<<
    set -euo pipefail

    methbat --version

    ln --symbolic --verbose "~{cpg_pileup_bed}" .

    label_args=(--profile-label ALL)
    IFS=';' read -r -a own <<< "~{labels}"
    for l in "${own[@]}"; do
      [ -n "$l" ] && [ "$l" != ALL ] && label_args+=(--profile-label "$l")
    done

    methbat profile \
      --input-pileup "~{basename(cpg_pileup_bed)}" \
      --input-regions "~{background}" \
      --output-region-profile "~{out_prefix}.methbat.scored.tsv" \
      --min-abs-zscore ~{zscore_min} \
      --min-abs-delta ~{methyl_delta_min} \
      --min-datasets ~{background_min_datasets} \
      "${label_args[@]}"

    python3 - "~{out_prefix}.methbat.scored.tsv" "~{methbat_map}" "~{sample_id}" "~{out_prefix}" <<'EOF'
    # ---- methbat_annotate.py (begin) ----
    import sys
    from collections import Counter
    scored, map_path, sample_id, prefix = sys.argv[1:5]
    EVENTS = ("HypoMethylated", "HyperMethylated", "HypoASM", "HyperASM")
    genes = {}
    with open(map_path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            c = line.rstrip("\n").split("\t")
            if c[0] == "chrom":
                continue
            genes[(c[0], c[1], c[2])] = (c[4], c[5], c[6])
    counts = {}
    header = None
    with open(scored) as fh, open(prefix + ".methbat.scored.annotated.tsv", "w") as ann, open(prefix + ".methbat.events.tsv", "w") as ev:
        for line in fh:
            if line.startswith("##"):
                ann.write(line)
                continue
            c = line.rstrip("\n").split("\t")
            if header is None:
                header = c
                idx = {k: i for i, k in enumerate(c)}
                ann.write("\t".join(c + ["hgnc_id", "gene_symbol", "region_type"]) + "\n")
                ev.write("\t".join(["sample_id"] + c + ["hgnc_id", "gene_symbol", "region_type"]) + "\n")
                continue
            g = genes.get((c[0], c[1], c[2]), ("", "", "intergenic" if genes else ""))
            ann.write("\t".join(c + list(g)) + "\n")
            label, compare, summary = c[idx["background_category"]], c[idx["compare_label"]], c[idx["summary_label"]]
            cnt = counts.setdefault(label, Counter())
            cnt["regions"] += 1
            cnt[compare] += 1
            if label == "ALL":
                cnt["summary_" + summary] += 1
            if compare in EVENTS:
                ev.write("\t".join([sample_id] + c + list(g)) + "\n")
    with open(prefix + ".methbat.counts.tsv", "w") as out:
        out.write("sample_id\tlabel\tmetric\tvalue\n")
        for label in sorted(counts, key=lambda l: (l != "ALL", l)):
            cnt = counts[label]
            for metric in ["regions", *EVENTS, "Uncategorized", "InsufficientData"]:
                out.write(f"{sample_id}\t{label}\t{metric}\t{cnt.get(metric, 0)}\n")
            for metric in sorted(k for k in cnt if k.startswith("summary_")):
                out.write(f"{sample_id}\t{label}\t{metric}\t{cnt[metric]}\n")
    print(f"{sample_id}: " + "; ".join(f"{l}: {sum(counts[l].get(e, 0) for e in EVENTS)} events" for l in counts), file=sys.stderr)
    # ---- methbat_annotate.py (end) ----
    EOF
  >>>

  output {
    File scored = "~{out_prefix}.methbat.scored.tsv"
    File scored_annotated = "~{out_prefix}.methbat.scored.annotated.tsv"
    File events = "~{out_prefix}.methbat.events.tsv"
    File counts = "~{out_prefix}.methbat.counts.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/methbat@sha256:281569947c0a6154f9a6bdaa1bb5f1e67dd160ccb49028d3218a49a37a0488fb"  # 1.1.0_build1
    cpu: 1
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

task ugc_wgw_methbat_compare {
  meta {
    description: "methbat compare of two label groups of the background (baseline:compare, e.g. XX:XY). A group absent from the background gives a header-only file and a note; a group below min_samples is compared anyway (every row InsufficientData, as MethBat labels it) and noted."
    outputs: {
      comparison: {
        description: "<out_prefix>.compare.<baseline>_vs_<compare>.tsv: MethBat's comparison (summary_comparison, z-scores and deltas, member counts per side)"
      },
      note: {
        description: "One line when a side is missing or below min_samples, else empty"
      }
    }
  }

  parameter_meta {
    background: {
      description: "The cohort background"
    }
    pair: {
      description: "baseline:compare labels"
    }
    zscore_min: {
      description: "methbat compare --min-zscore"
    }
    methyl_delta_min: {
      description: "methbat compare --min-delta"
    }
    min_samples: {
      description: "methbat compare --min-samples"
    }
    out_prefix: {
      description: "Output prefix"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File background
    String pair
    Float zscore_min = 3.0
    Float methyl_delta_min = 20
    Int min_samples = 10
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(background, "GB") * 3 + 5)

  command <<<
    set -euo pipefail

    methbat --version

    PAIR="~{pair}"
    baseline=${PAIR%%:*}
    compare=${PAIR#*:}
    if [ -z "$baseline" ] || [ -z "$compare" ] || [ "$baseline" = "$PAIR" ]; then
      echo "compare_pairs entry '~{pair}' is not baseline:compare" >&2
      exit 1
    fi
    out="~{out_prefix}.compare.${baseline}_vs_${compare}.tsv"
    echo "$out" > out_name.txt
    : > note.txt

    # members per label: the dataset counts of the first region's rows
    python3 - "~{background}" "$baseline" "$compare" ~{min_samples} <<'EOF' > status.txt
    import sys
    path, baseline, compare, min_samples = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    counts, first = {}, None
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("##"):
                continue
            c = line.rstrip("\n").split("\t")
            if header is None:
                header = c
                idx = {k: i for i, k in enumerate(c)}
                continue
            key = (c[0], c[1], c[2])
            if first is None:
                first = key
            if key != first:
                break
            counts[c[idx["data_category"]]] = sum(int(c[idx[k]]) for k in ("NoData", "Uncategorized", "Methylated", "Unmethylated", "AlleleSpecificMethylation"))
    missing = [l for l in (baseline, compare) if l not in counts]
    if missing:
        print("skip\t" + f"compare {baseline}:{compare}: no member carries the label {', '.join(missing)}; comparison not made")
    elif min(counts[baseline], counts[compare]) < min_samples:
        print("run\t" + f"compare {baseline}:{compare}: {counts[baseline]} vs {counts[compare]} members, below min_samples {min_samples}: every row is InsufficientData")
    else:
        print("run\t")
    EOF
    status=$(cut -f1 status.txt)
    cut -f2 status.txt | sed '/^$/d' > note.txt

    if [ "$status" = run ]; then
      methbat compare \
        --input-profile "~{background}" \
        --output-comparison "$out" \
        --baseline-category "$baseline" \
        --compare-category "$compare" \
        --min-zscore ~{zscore_min} \
        --min-delta ~{methyl_delta_min} \
        --min-samples ~{min_samples}
    else
      printf 'chrom\tstart\tend\tbaseline_category\tcompare_category\tsummary_comparison\tzscore_avg_abs_meth_deltas\tdelta_avg_abs_meth_deltas\tbaseline_num_phased\tcompare_num_phased\tzscore_avg_combined_methyls\tdelta_avg_combined_methyls\tbaseline_num_samples\tcompare_num_samples\n' > "$out"
    fi
  >>>

  output {
    File comparison = read_string("out_name.txt")
    String note = read_string("note.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/methbat@sha256:281569947c0a6154f9a6bdaa1bb5f1e67dd160ccb49028d3218a49a37a0488fb"  # 1.1.0_build1
    cpu: 1
    memory: "2 GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}

task ugc_wgw_methbat_joint_segment {
  meta {
    description: "methbat joint-segment over every member's 5mC pileup: the cohort's methylome segmented into Methylated, Unmethylated and AlleleSpecificMethylation regions without predefined regions, with the cohort-averaged combined methylation and ASM spread per segment. Whole genome in one task; its memory grows with the members (off by default until measured)."
    outputs: {
      regions_bed: {
        description: "<out_prefix>.joint.meth_regions.bed: merged segments, the fourth column the category"
      },
      combined_methyl_bedgraph: {
        description: "<out_prefix>.joint.combined_methyl.bedgraph: cohort-averaged combined methylation per segment"
      },
      asm_bedgraph: {
        description: "<out_prefix>.joint.asm.bedgraph: the spread between haplotypes per segment"
      }
    }
  }

  parameter_meta {
    sample_ids: {
      description: "Members, cohort order"
    }
    cpg_pileup_beds: {
      description: "Their 5mC pileups, aligned"
    }
    out_prefix: {
      description: "Output prefix"
    }
    threads: {
      description: "Threads"
    }
    mem_gb: {
      description: "Memory (GB)"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[String] sample_ids
    Array[File] cpg_pileup_beds
    String out_prefix
    Int threads = 8
    Int mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(cpg_pileup_beds, "GB") * 1.5 + 20)

  command <<<
    set -euo pipefail

    methbat --version

    mkdir -p in
    python3 - "~{write_lines(sample_ids)}" "~{write_lines(cpg_pileup_beds)}" > collection.tsv <<'EOF'
    import os, sys
    ids, files = ([l.rstrip("\n") for l in open(p)] for p in sys.argv[1:3])
    print("identifier\tfilename\tlabels")
    for i, (sid, path) in enumerate(zip(ids, files)):
        link = os.path.join("in", f"{i}_{os.path.basename(path)}")
        os.symlink(path, link)
        print(f"{sid}\t{os.path.abspath(link)}\t")
    EOF

    methbat joint-segment \
      --threads ~{threads} \
      --input-collection collection.tsv \
      --output-prefix "~{out_prefix}.joint"
  >>>

  output {
    File regions_bed = "~{out_prefix}.joint.meth_regions.bed"
    File combined_methyl_bedgraph = "~{out_prefix}.joint.combined_methyl.bedgraph"
    File asm_bedgraph = "~{out_prefix}.joint.asm.bedgraph"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/methbat@sha256:281569947c0a6154f9a6bdaa1bb5f1e67dd160ccb49028d3218a49a37a0488fb"  # 1.1.0_build1
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

task ugc_wgw_methbat_deconvolve {
  meta {
    description: "methbat deconvolve: cell-type proportions of one member's 5mC pileup against a reference atlas (non-negative least squares). The shipped atlas is nanomix's 39Bisulfite set, an example not validated for our tissue."
    outputs: {
      estimate: {
        description: "<out_prefix>.methbat.deconvolution.json: MethBat's estimate (metrics.normalized_fraction per cell type, the atlas md5)"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Member"
    }
    cpg_pileup_bed: {
      description: "Its 5mCpG pileup BED"
    }
    atlas_tsv: {
      description: "The cell-type atlas in MethBat format (reference map key methbat_atlas_tsv)"
    }
    out_prefix: {
      description: "Output prefix (<sample_id>.<ref_name>)"
    }
    mem_gb: {
      description: "Memory (GB)"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    File cpg_pileup_bed
    File atlas_tsv
    String out_prefix
    Int mem_gb = 4
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(cpg_pileup_bed, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    methbat --version
    echo "~{sample_id}" >&2

    ln --symbolic --verbose "~{cpg_pileup_bed}" .

    methbat deconvolve \
      --input-pileup "~{basename(cpg_pileup_bed)}" \
      --atlas-regions "~{atlas_tsv}" \
      --output-estimate "~{out_prefix}.methbat.deconvolution.json"
  >>>

  output {
    File estimate = "~{out_prefix}.methbat.deconvolution.json"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/methbat@sha256:281569947c0a6154f9a6bdaa1bb5f1e67dd160ccb49028d3218a49a37a0488fb"  # 1.1.0_build1
    cpu: 1
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

task ugc_wgw_methyl_summary {
  meta {
    description: "The methylation tables the driver reads: a long-format summary (members, labels, events per member, the summary-label distribution, imprinting counts, the genes with the most members carrying an event, the events per region type, the comparisons' label counts, deconvolution's top cell type, notes), the per-region table (every region and label of the background with the dataset counts, the event counts across members and the gene columns), the cohort events file (every member's Hypo/Hyper rows, annotated, sorted by gene then member), the imprinting table (member x region from methbat report, with where the report was made) and one row per member."
    outputs: {
      summary: {
        description: "<out_prefix>.summary.tsv: resource, metric, key, value"
      },
      regions: {
        description: "<out_prefix>.regions.tsv: chrom, start, end, region_label, label, n_datasets, n_nodata, n_uncategorized, n_methylated, n_unmethylated, n_asm, avg_combined_methyl, stdev_combined_methyl, avg_abs_meth_delta, stdev_abs_meth_delta, n_hypo, n_hyper, n_hypo_asm, n_hyper_asm, hgnc_id, gene_symbol, region_type"
      },
      cohort_events: {
        description: "<out_prefix>.events.tsv: sample_id, the scored row, hgnc_id, gene_symbol, region_type; Hypo/Hyper rows only, sorted by gene, position, member"
      },
      imprinting: {
        description: "<out_prefix>.imprinting.tsv: sample_id, the methbat report row, report_source (run or stage)"
      },
      samples: {
        description: "<out_prefix>.samples.tsv: one row per member: labels, event counts, summary-label counts, imprinting counts, report_source, deconvolution's top cell type"
      }
    }
  }

  parameter_meta {
    sample_ids: {
      description: "Members, cohort order"
    }
    labels: {
      description: "Their labels, ';'-joined, aligned"
    }
    events: {
      description: "Their events files (ugc_wgw_methbat_profile_bg), aligned"
    }
    counts: {
      description: "Their counts files, aligned"
    }
    reports: {
      description: "Their methbat reports, aligned"
    }
    report_sources: {
      description: "run or stage per member, aligned"
    }
    comparisons: {
      description: "The comparison files (ugc_wgw_methbat_compare), one per pair"
    }
    compare_notes: {
      description: "Their notes, aligned (empty strings for none)"
    }
    background: {
      description: "The cohort background"
    }
    methbat_map: {
      description: "The CpG-island gene map"
    }
    deconvolution: {
      description: "The members' deconvolution estimates, aligned with sample_ids, when run"
    }
    background_min_datasets: {
      description: "The minimum members per label for a comparison, for the notes"
    }
    out_prefix: {
      description: "Output prefix"
    }
    mem_gb: {
      description: "Memory (GB)"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[String] sample_ids
    Array[String] labels
    Array[File] events
    Array[File] counts
    Array[File] reports
    Array[String] report_sources
    Array[File] comparisons = []
    Array[String] compare_notes = []
    File background
    File methbat_map
    Array[File] deconvolution = []
    Int background_min_datasets = 10
    String out_prefix
    Int mem_gb = 4
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(background, "GB") * 2 + size(events, "GB") * 3 + 10)

  command <<<
    set -euo pipefail

    python3 - "~{out_prefix}" "~{write_lines(sample_ids)}" "~{write_lines(labels)}" "~{write_lines(events)}" "~{write_lines(counts)}" \
      "~{write_lines(reports)}" "~{write_lines(report_sources)}" "~{write_lines(comparisons)}" "~{write_lines(compare_notes)}" \
      "~{background}" "~{methbat_map}" "~{write_lines(deconvolution)}" ~{background_min_datasets} <<'EOF'
    # ---- methyl_summary.py (begin) ----
    import json, sys
    from collections import Counter, defaultdict
    (prefix, ids_f, labels_f, events_f, counts_f, reports_f, sources_f, compares_f, notes_f,
     background, map_path, deconv_f, min_datasets) = sys.argv[1:14]
    min_datasets = int(min_datasets)
    read = lambda p: [l.rstrip("\n") for l in open(p)]
    ids, labels, events, counts, reports, sources = read(ids_f), read(labels_f), read(events_f), read(counts_f), read(reports_f), read(sources_f)
    compares, notes, deconvs = [l for l in read(compares_f) if l], read(notes_f), [l for l in read(deconv_f) if l]
    EVENTS = ("HypoMethylated", "HyperMethylated", "HypoASM", "HyperASM")
    EVENT_COL = {"HypoMethylated": "n_hypo", "HyperMethylated": "n_hyper", "HypoASM": "n_hypo_asm", "HyperASM": "n_hyper_asm"}

    def table(path):
        """(header, rows) of a TSV with ## lines."""
        header, rows = None, []
        with open(path) as fh:
            for line in fh:
                if line.startswith("##") or not line.strip():
                    continue
                c = line.rstrip("\n").split("\t")
                if header is None:
                    header = c
                else:
                    rows.append(c)
        return header or [], rows

    genes = {}
    for c in table(map_path)[1]:
        genes[(c[0], c[1], c[2])] = (c[4], c[5], c[6])

    out_rows = []
    def add(resource, metric, key, value):
        out_rows.append((resource, metric, str(key), str(value)))

    # ---- events: concatenated, sorted by gene, position, member; counted per region x label, per gene, per region type
    ev_header, ev_rows = None, []
    per_region_events = defaultdict(Counter)        # (chrom, start, end, label) -> compare_label -> n
    gene_members = defaultdict(set)                 # gene symbol -> members with an event
    type_events = Counter()
    for sid, path in zip(ids, events):
        h, rows = table(path)
        if ev_header is None and h:
            ev_header = h
            idx = {k: i for i, k in enumerate(h)}
        for c in rows:
            ev_rows.append(c)
            per_region_events[(c[idx["chrom"]], c[idx["start"]], c[idx["end"]], c[idx["background_category"]])][c[idx["compare_label"]]] += 1
            type_events[c[idx["region_type"]] or "unannotated"] += 1
            for g in (c[idx["gene_symbol"]] or "").split(";"):
                if g:
                    gene_members[g].add(c[idx["sample_id"]])
    if ev_header is None:
        ev_header = ["sample_id"]
    else:
        def ev_key(c):
            first_gene = (c[idx["gene_symbol"]] or "").split(";")[0]
            return (first_gene == "", first_gene, c[idx["chrom"]], int(c[idx["start"]]), c[idx["sample_id"]], c[idx["background_category"]])
        ev_rows.sort(key=ev_key)
    with open(prefix + ".events.tsv", "w") as out:
        out.write("\t".join(ev_header) + "\n")
        for c in ev_rows:
            out.write("\t".join(c) + "\n")

    # ---- per-member counts
    member_counts = {sid: defaultdict(Counter) for sid in ids}   # sid -> label -> metric -> value
    for sid, path in zip(ids, counts):
        for c in table(path)[1]:
            member_counts[sid][c[1]][c[2]] += int(c[3])
    label_members = Counter()
    for sid, lab in zip(ids, labels):
        label_members["ALL"] += 1
        for l in lab.split(";"):
            if l and l != "ALL":
                label_members[l] += 1

    # ---- the background: per region x label
    bg_header, bg_rows = table(background)
    bi = {k: i for i, k in enumerate(bg_header)}
    with open(prefix + ".regions.tsv", "w") as out:
        out.write("\t".join(["chrom", "start", "end", "region_label", "label", "n_datasets", "n_nodata", "n_uncategorized", "n_methylated",
                             "n_unmethylated", "n_asm", "avg_combined_methyl", "stdev_combined_methyl", "avg_abs_meth_delta", "stdev_abs_meth_delta",
                             "n_hypo", "n_hyper", "n_hypo_asm", "n_hyper_asm", "hgnc_id", "gene_symbol", "region_type"]) + "\n")
        n_regions = set()
        for c in bg_rows:
            key = (c[bi["chrom"]], c[bi["start"]], c[bi["end"]])
            n_regions.add(key)
            label = c[bi["data_category"]]
            cats = [int(c[bi[k]]) for k in ("NoData", "Uncategorized", "Methylated", "Unmethylated", "AlleleSpecificMethylation")]
            ev = per_region_events.get(key + (label,), Counter())
            g = genes.get(key, ("", "", ""))
            out.write("\t".join([*key, c[bi["region_label"]], label, str(sum(cats)), *map(str, cats),
                                 c[bi["avg_combined_methyls"]], c[bi["stdev_combined_methyls"]], c[bi["avg_abs_meth_deltas"]], c[bi["stdev_abs_meth_deltas"]],
                                 *[str(ev.get(e, 0)) for e in EVENTS], *g]) + "\n")

    # ---- imprinting
    imp_header, imp_rows = None, []
    imp_counts = {sid: Counter() for sid in ids}
    for sid, path, source in zip(ids, reports, sources):
        h, rows = table(path)
        if imp_header is None and h:
            imp_header = h
            ii = {k: i for i, k in enumerate(h)}
        for c in rows:
            imp_rows.append([sid] + c + [source])
            imp_counts[sid][c[ii["report_summary"]]] += 1
    with open(prefix + ".imprinting.tsv", "w") as out:
        out.write("\t".join(["sample_id"] + (imp_header or []) + ["report_source"]) + "\n")
        for c in imp_rows:
            out.write("\t".join(c) + "\n")
    region_anomalous = Counter(c[1 + ii["region_label"]] for c in imp_rows if c[1 + ii["report_summary"]] in ("Anomalous", "AnomalousQcWarning")) if imp_rows else Counter()

    # ---- deconvolution
    top_cell = {}
    for sid, path in zip(ids, deconvs):
        try:
            frac = json.load(open(path)).get("metrics", {}).get("normalized_fraction", {})
        except (OSError, ValueError):
            frac = {}
        if frac:
            cell, value = max(frac.items(), key=lambda kv: kv[1])
            top_cell[sid] = (cell, value)
            add("deconvolution", "top_cell_type", sid, cell)
            add("deconvolution", "top_fraction", sid, f"{value:.3f}")

    # ---- samples table
    with open(prefix + ".samples.tsv", "w") as out:
        out.write("\t".join(["sample_id", "labels", "n_hypo", "n_hyper", "n_hypo_asm", "n_hyper_asm", "n_events", "n_insufficient", "n_regions_scored",
                             "n_methylated", "n_unmethylated", "n_asm", "n_uncategorized", "n_nodata",
                             "imprinting_pass", "imprinting_inconclusive", "imprinting_anomalous", "imprinting_anomalous_qc", "report_source",
                             "deconvolution_top_cell_type", "deconvolution_top_fraction"]) + "\n")
        for sid, lab, source in zip(ids, labels, sources):
            cnt = member_counts[sid]
            all_c = cnt.get("ALL", Counter())
            per_event = {e: sum(cnt[l].get(e, 0) for l in cnt) for e in EVENTS}
            n_events = sum(per_event.values())
            scored = all_c.get("regions", 0) - all_c.get("summary_NoData", 0)
            top = top_cell.get(sid, ("", ""))
            out.write("\t".join(map(str, [sid, lab, per_event["HypoMethylated"], per_event["HyperMethylated"], per_event["HypoASM"], per_event["HyperASM"], n_events,
                                          sum(cnt[l].get("InsufficientData", 0) for l in cnt), scored,
                                          all_c.get("summary_Methylated", 0), all_c.get("summary_Unmethylated", 0), all_c.get("summary_AlleleSpecificMethylation", 0),
                                          all_c.get("summary_Uncategorized", 0), all_c.get("summary_NoData", 0),
                                          imp_counts[sid].get("PASS", 0), imp_counts[sid].get("Inconclusive", 0), imp_counts[sid].get("Anomalous", 0),
                                          imp_counts[sid].get("AnomalousQcWarning", 0), source, top[0], f"{top[1]:.3f}" if top[0] else ""])) + "\n")
            add("members", "n_events", sid, n_events)
            for e in EVENTS:
                add("members", EVENT_COL[e], sid, per_event[e])
            add("members", "n_regions_scored", sid, scored)
            for k in ("Methylated", "Unmethylated", "AlleleSpecificMethylation", "Uncategorized", "NoData"):
                add("members", "summary_" + k, sid, all_c.get("summary_" + k, 0))
            for k in ("PASS", "Inconclusive", "Anomalous", "AnomalousQcWarning"):
                add("imprinting", k, sid, imp_counts[sid].get(k, 0))

    # ---- the summary
    add("cohort", "n_members", "", len(ids))
    add("cohort", "n_regions", "", len(n_regions))
    add("cohort", "n_events", "", len(ev_rows))
    add("cohort", "background_min_datasets", "", min_datasets)
    add("cohort", "reports_made_in_stage", "", sum(1 for s in sources if s == "stage"))
    for label, n in sorted(label_members.items(), key=lambda kv: (kv[0] != "ALL", kv[0])):
        add("labels", "n_members", label, n)
        if n < min_datasets:
            add("notes", "label_below_min", label, f"{n} member(s): comparisons against {label} are InsufficientData")
    for region, n in sorted(region_anomalous.items(), key=lambda kv: (-kv[1], kv[0])):
        add("imprinting_regions", "n_anomalous", region, n)
    for g, members in sorted(gene_members.items(), key=lambda kv: (-len(kv[1]), kv[0]))[:25]:
        add("genes", "members_with_event", g, len(members))
    for t, n in sorted(type_events.items(), key=lambda kv: (-kv[1], kv[0])):
        add("region_type", "n_events", t, n)
    for path, note in zip(compares, notes + [""] * (len(compares) - len(notes))):
        h, rows = table(path)
        if not h:
            continue
        ci = {k: i for i, k in enumerate(h)}
        pair = (rows[0][ci["baseline_category"]] + ":" + rows[0][ci["compare_category"]]) if rows else path.rsplit(".compare.", 1)[-1].rsplit(".tsv", 1)[0].replace("_vs_", ":")
        for label, n in sorted(Counter(c[ci["summary_comparison"]] or "none" for c in rows).items()):
            add("compare", pair, label, n)
        if note:
            add("notes", "compare", pair, note)
    with open(prefix + ".summary.tsv", "w") as out:
        out.write("resource\tmetric\tkey\tvalue\n")
        for r in out_rows:
            out.write("\t".join(r) + "\n")
    print(f"{len(ids)} members, {len(n_regions)} regions, {len(ev_rows)} events, {len(imp_rows)} imprinting rows", file=sys.stderr)
    # ---- methyl_summary.py (end) ----
    EOF
  >>>

  output {
    File summary = "~{out_prefix}.summary.tsv"
    File regions = "~{out_prefix}.regions.tsv"
    File cohort_events = "~{out_prefix}.events.tsv"
    File imprinting = "~{out_prefix}.imprinting.tsv"
    File samples = "~{out_prefix}.samples.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
    cpu: 1
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
