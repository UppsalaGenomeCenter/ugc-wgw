version 1.0

# ugc-pacbio-wgw — somalier: per-sample extract and sex resolution (upstream), relatedness, ancestry and the QC summary (cohort_qc)
# origin: none
# see docs/DESIGN.md §5, §8.1; docs/ENTRYPOINTS.md §9.1

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_somalier_extract {
  meta {
    description: "Genotype-like information of one sample at somalier's ~17 000 informative sites, from the aligned BAM: one small .somalier file the cohort stage reads instead of the BAM. The sample name inside the file is sample_id (SOMALIER_SAMPLE_NAME), not the BAM's read group."
    outputs: {
      somalier: {
        description: "<sample_id>.somalier"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID (the name stored in the extract)"
    }
    aligned_bam: {
      description: "Aligned BAM"
    }
    aligned_bam_index: {
      description: "Index for aligned_bam"
    }
    sites_vcf: {
      description: "somalier sites VCF (reference map somalier_sites_vcf)"
    }
    ref_fasta: {
      description: "Reference genome FASTA"
    }
    ref_index: {
      description: "Reference genome FASTA index"
    }
    mem_gb: {
      description: "Memory (GB)"
    }
    ugc_wgw_container_registry: {
      description: "Registry holding the ugc-built images"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    File aligned_bam
    File aligned_bam_index
    File sites_vcf
    File ref_fasta
    File ref_index
    Int mem_gb = 2
    String ugc_wgw_container_registry
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(ref_fasta, "GB") + 5)

  command <<<
    set -euo pipefail

    somalier 2>&1 | head -1 >&2

    # somalier opens the index next to the BAM and the .fai next to the FASTA
    ln --symbolic --verbose "~{aligned_bam}" .
    ln --symbolic --verbose "~{aligned_bam_index}" .
    ln --symbolic --verbose "~{ref_fasta}" .
    ln --symbolic --verbose "~{ref_index}" .

    SOMALIER_SAMPLE_NAME="~{sample_id}" somalier extract \
      --sites "~{sites_vcf}" \
      --fasta "~{basename(ref_fasta)}" \
      --out-dir . \
      "~{basename(aligned_bam)}"
    test -s "~{sample_id}.somalier"
  >>>

  output {
    File somalier = "~{sample_id}.somalier"
  }

  runtime {
    docker: "~{ugc_wgw_container_registry}/somalier@sha256:8b4c7d98c83fdd1e459baad80876f6eeb9a607a488ef7ce2b3db0a4ae707fe64"  # somalier 0.3.7
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

task ugc_wgw_somalier_relate_one {
  meta {
    description: "somalier relate on one extract: the per-sample table (depth, het and hom-alt counts, chrX and chrY statistics) the sex resolution reads; no pairs"
    outputs: {
      samples_tsv: {
        description: "<sample_id>.samples.tsv: somalier's per-sample statistics (one row)"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID"
    }
    somalier: {
      description: "The sample's .somalier extract"
    }
    sites_vcf: {
      description: "somalier sites VCF (the one the extract was made with)"
    }
    min_depth: {
      description: "somalier relate --min-depth: sites below this depth are not counted"
    }
    min_ab: {
      description: "somalier relate --min-ab: heterozygote allele-balance threshold"
    }
    ugc_wgw_container_registry: {
      description: "Registry holding the ugc-built images"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    File somalier
    File sites_vcf
    Int min_depth = 7
    Float min_ab = 0.3
    String ugc_wgw_container_registry
    RuntimeAttributes runtime_attributes
  }

  command <<<
    set -euo pipefail

    somalier 2>&1 | head -1 >&2

    somalier relate \
      --sites "~{sites_vcf}" \
      --min-depth ~{min_depth} \
      --min-ab ~{min_ab} \
      --output-prefix "~{sample_id}." \
      "~{somalier}"
    test -s "~{sample_id}.samples.tsv"
  >>>

  output {
    File samples_tsv = "~{sample_id}.samples.tsv"
  }

  runtime {
    docker: "~{ugc_wgw_container_registry}/somalier@sha256:8b4c7d98c83fdd1e459baad80876f6eeb9a607a488ef7ce2b3db0a4ae707fe64"  # somalier 0.3.7
    cpu: 1
    memory: "1 GiB"
    disk: "5 GB"
    disks: "local-disk 5 HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}

task ugc_wgw_resolve_sex {
  meta {
    description: "One sex for the sex-dependent callers, from three sources in fixed precedence: declared (the sample sheet), coverage (mosdepth's chrY inference), somalier (chrX heterozygosity from the one-sample relate table, somalier's own rule: X_het/X_hom_alt below 0.05 male, above 0.4 female, given enough X sites and depth). Every disagreement between available sources is a message; no source at all leaves the sex empty, which the callers treat as upstream does."
    outputs: {
      resolved_sex: {
        description: "MALE, FEMALE or empty"
      },
      resolved_sex_source: {
        description: "declared, coverage, somalier or none"
      },
      somalier_sex: {
        description: "somalier's call: MALE, FEMALE or empty (too few X sites, too little depth, or an ambiguous ratio)"
      },
      sex_tsv: {
        description: "<sample_id>.sex.tsv: one row with every source, the resolved value and its reason, and the X/Y statistics"
      },
      msg: {
        description: "Messages: one per disagreement, one when no source could sex the sample"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID"
    }
    declared_sex: {
      description: "The sample sheet's sex (MALE or FEMALE), when it has one"
    }
    coverage_sex: {
      description: "mosdepth's inferred sex (MALE, FEMALE or empty)"
    }
    somalier_samples_tsv: {
      description: "The one-sample somalier relate table"
    }
    min_x_sites: {
      description: "Minimum chrX sites with a genotype before somalier's call is used"
    }
    min_x_depth: {
      description: "Minimum mean chrX depth before somalier's call is used"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    String? declared_sex
    String coverage_sex
    File somalier_samples_tsv
    Int min_x_sites = 50
    Float min_x_depth = 5.0
    RuntimeAttributes runtime_attributes
  }

  command <<<
    set -euo pipefail

    python3 --version >&2

    python3 - "~{sample_id}" "~{default='' declared_sex}" "~{coverage_sex}" "~{somalier_samples_tsv}" ~{min_x_sites} ~{min_x_depth} <<'EOF'
    import csv, sys
    sid, declared, coverage, table, min_sites, min_depth = sys.argv[1], sys.argv[2].upper(), sys.argv[3].upper(), sys.argv[4], int(sys.argv[5]), float(sys.argv[6])
    norm = {"MALE": "MALE", "FEMALE": "FEMALE", "M": "MALE", "F": "FEMALE", "XY": "MALE", "XX": "FEMALE", "1": "MALE", "2": "FEMALE"}
    declared, coverage = norm.get(declared, ""), norm.get(coverage, "")
    with open(table) as fh:
        rows = list(csv.DictReader((line.lstrip("#") for line in fh), delimiter="\t"))
    row = rows[0] if rows else {}
    def num(key):
        try:
            return float(row.get(key, ""))
        except (TypeError, ValueError):
            return 0.0
    x_n, x_het, x_hom_alt, x_depth = int(num("X_n")), int(num("X_het")), int(num("X_hom_alt")), num("X_depth_mean")
    gt_depth, y_depth = num("gt_depth_mean"), num("Y_depth_mean")
    somalier, reason = "", ""
    ratio = x_het / x_hom_alt if x_hom_alt else (float("inf") if x_het else 0.0)
    if x_n < min_sites:
        reason = f"somalier: {x_n} chrX sites genotyped, below {min_sites}"
    elif x_depth < min_depth:
        reason = f"somalier: chrX depth {x_depth:.1f} below {min_depth}"
    elif ratio < 0.05:
        somalier, reason = "MALE", f"somalier: X_het/X_hom_alt {ratio:.3f} below 0.05"
    elif ratio > 0.4:
        somalier, reason = "FEMALE", f"somalier: X_het/X_hom_alt {ratio:.3f} above 0.4"
    else:
        reason = f"somalier: X_het/X_hom_alt {ratio:.3f} is ambiguous (0.05 to 0.4)"
    if declared:
        resolved, source = declared, "declared"
    elif coverage:
        resolved, source = coverage, "coverage"
    elif somalier:
        resolved, source = somalier, "somalier"
    else:
        resolved, source = "", "none"
    messages = []
    have = [(n, v) for n, v in (("declared", declared), ("coverage", coverage), ("somalier", somalier)) if v]
    for i, (n1, v1) in enumerate(have):
        for n2, v2 in have[i + 1:]:
            if v1 != v2:
                messages.append(f"{sid}: sex sources disagree: {n1} {v1}, {n2} {v2} (using {source} {resolved})")
    if not resolved:
        tried = "; ".join(f"{n} {'empty' if not v else v}" for n, v in (("declared", declared), ("coverage", coverage), ("somalier", somalier)))
        messages.append(f"{sid}: no source could sex the sample ({tried}; {reason}); the sex-dependent callers assume XX")
    y_norm = 2 * y_depth / gt_depth if gt_depth else 0.0
    with open(f"{sid}.sex.tsv", "w") as out:
        out.write("sample_id\tdeclared\tcoverage\tsomalier\tresolved\tsource\tx_n\tx_het\tx_hom_alt\tx_het_ratio\tx_depth_mean\ty_depth_norm\treason\n")
        out.write("\t".join([sid, declared, coverage, somalier, resolved, source, str(x_n), str(x_het), str(x_hom_alt),
                             f"{ratio:.4f}" if ratio != float("inf") else "inf", f"{x_depth:.2f}", f"{y_norm:.3f}", reason]) + "\n")
    with open("resolved_sex.txt", "w") as out:
        out.write(resolved + "\n")
    with open("resolved_sex_source.txt", "w") as out:
        out.write(source + "\n")
    with open("somalier_sex.txt", "w") as out:
        out.write(somalier + "\n")
    with open("messages.txt", "w") as out:
        for m in messages:
            out.write(m + "\n")
    EOF
  >>>

  output {
    String resolved_sex = read_string("resolved_sex.txt")
    String resolved_sex_source = read_string("resolved_sex_source.txt")
    String somalier_sex = read_string("somalier_sex.txt")
    File sex_tsv = "~{sample_id}.sex.tsv"
    Array[String] msg = read_lines("messages.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
    cpu: 1
    memory: "1 GiB"
    disk: "5 GB"
    disks: "local-disk 5 HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}

task ugc_wgw_somalier_ped {
  meta {
    description: "A PED file of the cohort from the sample sheet's sex and parents: family = cohort, parents outside the cohort kept as founders of their own, sex 1/2/-9, phenotype -9"
    outputs: {
      ped: {
        description: "<cohort_id>.qc.ped"
      }
    }
  }

  parameter_meta {
    cohort_id: {
      description: "Cohort ID (the PED family ID)"
    }
    sample_ids: {
      description: "Cohort members, cohort order"
    }
    sample_sexes: {
      description: "XX, XY or empty per member, aligned with sample_ids"
    }
    father_ids: {
      description: "Father per member (empty where the sheet has none), aligned with sample_ids"
    }
    mother_ids: {
      description: "Mother per member (empty where none), aligned with sample_ids"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String cohort_id
    Array[String] sample_ids
    Array[String] sample_sexes
    Array[String] father_ids
    Array[String] mother_ids
    RuntimeAttributes runtime_attributes
  }

  command <<<
    set -euo pipefail

    python3 --version >&2

    python3 - "~{cohort_id}" "~{write_lines(sample_ids)}" "~{write_lines(sample_sexes)}" "~{write_lines(father_ids)}" "~{write_lines(mother_ids)}" <<'EOF'
    import sys
    cohort, ids_f, sexes_f, fathers_f, mothers_f = sys.argv[1:6]
    def lines(p):
        return [l.rstrip("\n") for l in open(p)]
    ids, sexes, fathers, mothers = lines(ids_f), lines(sexes_f), lines(fathers_f), lines(mothers_f)
    assert len(ids) == len(sexes) == len(fathers) == len(mothers), "arrays are not aligned"
    code = {"XY": "1", "XX": "2"}
    members = set(ids)
    rows, founders = [], []
    for sid, sex, fa, mo in zip(ids, sexes, fathers, mothers):
        rows.append([cohort, sid, fa or "-9", mo or "-9", code.get(sex.upper(), "-9"), "-9"])
        for parent, psex in ((fa, "1"), (mo, "2")):
            if parent and parent not in members and parent not in founders:
                founders.append(parent)
                rows.append([cohort, parent, "-9", "-9", psex, "-9"])
    with open(f"{cohort}.qc.ped", "w") as out:
        for r in rows:
            out.write("\t".join(r) + "\n")
    print(f"{len(ids)} members, {len(founders)} parent(s) outside the cohort", file=sys.stderr)
    EOF
  >>>

  output {
    File ped = "~{cohort_id}.qc.ped"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
    cpu: 1
    memory: "1 GiB"
    disk: "5 GB"
    disks: "local-disk 5 HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}

task ugc_wgw_somalier_relate {
  meta {
    description: "somalier relate over the cohort's extracts with the sheet's pedigree: pairwise relatedness, IBS0/IBS2 and expected relatedness, per-sample depth, heterozygosity and chrX/chrY statistics, sample groups, somalier's HTML. Every pair is reported (SOMALIER_REPORT_ALL_PAIRS), with a fixed seed."
    outputs: {
      pairs: {
        description: "<out_prefix>.pairs.tsv: one row per pair"
      },
      samples: {
        description: "<out_prefix>.samples.tsv: one row per sample"
      },
      groups: {
        description: "<out_prefix>.groups.tsv: samples that look identical"
      },
      html: {
        description: "<out_prefix>.html: somalier's interactive report"
      }
    }
  }

  parameter_meta {
    out_prefix: {
      description: "Output prefix"
    }
    somalier_files: {
      description: "One .somalier extract per member"
    }
    ped: {
      description: "PED file of the cohort (ugc_wgw_somalier_ped)"
    }
    sites_vcf: {
      description: "somalier sites VCF (the one the extracts were made with)"
    }
    min_depth: {
      description: "somalier relate --min-depth"
    }
    min_ab: {
      description: "somalier relate --min-ab"
    }
    report_all_pairs: {
      description: "Report every pair (somalier subsamples unrelated pairs to 200 000 otherwise)"
    }
    threads: {
      description: "Threads"
    }
    mem_gb: {
      description: "Memory (GB)"
    }
    ugc_wgw_container_registry: {
      description: "Registry holding the ugc-built images"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String out_prefix
    Array[File] somalier_files
    File ped
    File sites_vcf
    Int min_depth = 7
    Float min_ab = 0.3
    Boolean report_all_pairs = true
    Int threads = 2
    Int mem_gb = 4
    String ugc_wgw_container_registry
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(somalier_files, "GB") * 2 + 5)

  command <<<
    set -euo pipefail

    somalier 2>&1 | head -1 >&2

    cp "~{write_lines(somalier_files)}" extracts.list
    export SOMALIER_SEED=1
    ~{if report_all_pairs then "export SOMALIER_REPORT_ALL_PAIRS=1" else ""}
    somalier relate \
      --sites "~{sites_vcf}" \
      --ped "~{ped}" \
      --min-depth ~{min_depth} \
      --min-ab ~{min_ab} \
      --output-prefix "~{out_prefix}." \
      extracts.list
    test -s "~{out_prefix}.samples.tsv"
  >>>

  output {
    File pairs = "~{out_prefix}.pairs.tsv"
    File samples = "~{out_prefix}.samples.tsv"
    File groups = "~{out_prefix}.groups.tsv"
    File html = "~{out_prefix}.html"
  }

  runtime {
    docker: "~{ugc_wgw_container_registry}/somalier@sha256:8b4c7d98c83fdd1e459baad80876f6eeb9a607a488ef7ce2b3db0a4ae707fe64"  # somalier 0.3.7
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

task ugc_wgw_somalier_ancestry {
  meta {
    description: "somalier ancestry: the cohort's extracts projected onto a PCA of the labelled 1000 Genomes extracts (superpopulations AFR, AMR, EAS, EUR, SAS); a description of structure, never a QC flag"
    outputs: {
      ancestry_tsv: {
        description: "<out_prefix>.somalier-ancestry.tsv: labelled and query samples with the predicted superpopulation, its probabilities and the principal components"
      },
      ancestry_html: {
        description: "<out_prefix>.somalier-ancestry.html"
      }
    }
  }

  parameter_meta {
    out_prefix: {
      description: "Output prefix"
    }
    somalier_files: {
      description: "One .somalier extract per member"
    }
    labelled_tar: {
      description: "The 1000 Genomes extracts as one tar (reference map somalier_1kg_tar)"
    }
    labels: {
      description: "The labels of the 1000 Genomes samples (reference map somalier_1kg_labels)"
    }
    n_pcs: {
      description: "Principal components kept"
    }
    mem_gb: {
      description: "Memory (GB)"
    }
    ugc_wgw_container_registry: {
      description: "Registry holding the ugc-built images"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String out_prefix
    Array[File] somalier_files
    File labelled_tar
    File labels
    Int n_pcs = 5
    Int mem_gb = 8
    String ugc_wgw_container_registry
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(labelled_tar, "GB") * 10 + size(somalier_files, "GB") + 5)

  command <<<
    set -euo pipefail

    somalier 2>&1 | head -1 >&2

    mkdir labelled
    tar -xzf "~{labelled_tar}" -C labelled
    n_labelled=$(find labelled -name '*.somalier' | wc -l)
    echo "labelled extracts: $n_labelled" >&2
    find labelled -name '*.somalier' | sort > labelled.list
    cp "~{write_lines(somalier_files)}" query.list
    # somalier takes the labelled set, then `++`, then the query set
    # shellcheck disable=SC2046
    somalier ancestry \
      --labels "~{labels}" \
      --n-pcs ~{n_pcs} \
      --output-prefix "~{out_prefix}." \
      $(cat labelled.list) ++ $(cat query.list)
    test -s "~{out_prefix}.somalier-ancestry.tsv"
  >>>

  output {
    File ancestry_tsv = "~{out_prefix}.somalier-ancestry.tsv"
    File ancestry_html = "~{out_prefix}.somalier-ancestry.html"
  }

  runtime {
    docker: "~{ugc_wgw_container_registry}/somalier@sha256:8b4c7d98c83fdd1e459baad80876f6eeb9a607a488ef7ce2b3db0a4ae707fe64"  # somalier 0.3.7
    cpu: 2
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

task ugc_wgw_qc_summary {
  meta {
    description: "The identity tables the driver reads: a long-format summary (the relatedness histogram of every pair in 0.05 bins, pair and member counts), one row per member (closest relative with its relatedness, IBS0 and expected relatedness, heterozygosity ratio, somalier's chrX call, the ancestry label and its probability, whether the extract was made in this stage), and every pair with relatedness at or above 0.1 or an expectation from the pedigree"
    outputs: {
      summary: {
        description: "<out_prefix>.summary.tsv: section, metric, key, value"
      },
      members: {
        description: "<out_prefix>.members.tsv: one row per member"
      },
      related_pairs: {
        description: "<out_prefix>.related_pairs.tsv: the pairs worth a look, pairs.tsv columns plus expected"
      }
    }
  }

  parameter_meta {
    out_prefix: {
      description: "Output prefix"
    }
    pairs: {
      description: "somalier relate pairs table"
    }
    samples: {
      description: "somalier relate samples table"
    }
    ancestry_tsv: {
      description: "somalier ancestry table, when the step ran"
    }
    sample_ids: {
      description: "Cohort members, cohort order"
    }
    sample_sexes: {
      description: "XX, XY or empty per member"
    }
    sample_sex_sources: {
      description: "declared, coverage, somalier or empty per member"
    }
    legacy_sample_ids: {
      description: "Per member, its ID when the extract was made in this stage (the per-sample run predates it), else empty"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String out_prefix
    File pairs
    File samples
    File? ancestry_tsv
    Array[String] sample_ids
    Array[String] sample_sexes
    Array[String] sample_sex_sources
    Array[String] legacy_sample_ids
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(pairs, "GB") * 2 + 5)

  command <<<
    set -euo pipefail

    python3 --version >&2

    python3 - "~{out_prefix}" "~{pairs}" "~{samples}" "~{default='' ancestry_tsv}" "~{write_lines(sample_ids)}" "~{write_lines(sample_sexes)}" "~{write_lines(sample_sex_sources)}" "~{write_lines(legacy_sample_ids)}" <<'EOF'
    import csv, sys
    prefix, pairs_f, samples_f, ancestry_f, ids_f, sexes_f, sources_f, legacy_f = sys.argv[1:9]
    def lines(p):
        return [l.rstrip("\n") for l in open(p)]
    def table(path):
        with open(path) as fh:
            return list(csv.DictReader((line.lstrip("#") for line in fh), delimiter="\t"))
    ids, sexes, sources, legacy = lines(ids_f), lines(sexes_f), lines(sources_f), {x for x in lines(legacy_f) if x}
    sources = sources + [""] * (len(ids) - len(sources))
    def num(v, default=0.0):
        try:
            return float(v)
        except (TypeError, ValueError):
            return default
    # pairs: one pass for the histogram, the closest relative per member and the pairs worth keeping
    hist = {}
    closest = {}
    keep = []
    n_pairs = 0
    with open(pairs_f) as fh:
        reader = csv.DictReader((line.lstrip("#") for line in fh), delimiter="\t")
        cols = reader.fieldnames or []
        for r in reader:
            n_pairs += 1
            rel = num(r.get("relatedness"))
            expected = num(r.get("expected_relatedness"), -1.0)
            b = f"{int(rel * 20) / 20:.2f}" if rel >= 0 else "<0"
            hist[b] = hist.get(b, 0) + 1
            for a, c in ((r["sample_a"], r["sample_b"]), (r["sample_b"], r["sample_a"])):
                if a not in closest or rel > closest[a][1]:
                    closest[a] = (c, rel, num(r.get("ibs0")), expected)
            if rel >= 0.1 or expected > 0:
                keep.append(r)
    with open(f"{prefix}.related_pairs.tsv", "w") as out:
        out.write("\t".join(cols) + "\texpected\n")
        for r in sorted(keep, key=lambda r: -num(r.get("relatedness"))):
            out.write("\t".join(r.get(c, "") for c in cols) + "\t" + ("" if num(r.get("expected_relatedness"), -1.0) < 0 else r.get("expected_relatedness", "")) + "\n")
    # samples: heterozygosity and somalier's X call
    srows = {r["sample_id"]: r for r in table(samples_f)}
    anc = {}
    if ancestry_f:
        for r in table(ancestry_f):
            probs = {k[:-5]: num(v) for k, v in r.items() if k.endswith("_prob")}
            anc[r.get("sample_id", "")] = (r.get("predicted_ancestry", ""), max(probs.values()) if probs else 0.0)
    with open(f"{prefix}.members.tsv", "w") as out:
        out.write("sample_id\tsex_group\tsex_source\tsomalier_sex\tx_n\tx_het_ratio\thet_ratio\tdepth_mean\tclosest_relative\tclosest_relatedness\tclosest_ibs0\tclosest_expected\tancestry\tancestry_prob\textracted_here\n")
        for sid, sex, src in zip(ids, sexes, sources):
            r = srows.get(sid, {})
            x_n, x_het, x_hom = int(num(r.get("X_n"))), int(num(r.get("X_het"))), int(num(r.get("X_hom_alt")))
            ratio = x_het / x_hom if x_hom else (float("inf") if x_het else 0.0)
            som = "" if x_n < 50 else ("MALE" if ratio < 0.05 else "FEMALE" if ratio > 0.4 else "")
            n_het, n_hom = num(r.get("n_het")), num(r.get("n_hom_alt"))
            het_ratio = n_het / (n_het + n_hom) if n_het + n_hom else 0.0
            c = closest.get(sid, ("", 0.0, 0.0, -1.0))
            a = anc.get(sid, ("", 0.0))
            out.write("\t".join([sid, sex, src, som, str(x_n), f"{ratio:.4f}" if ratio != float("inf") else "inf", f"{het_ratio:.4f}",
                                 f"{num(r.get('gt_depth_mean')):.1f}", c[0], f"{c[1]:.3f}", f"{c[2]:.0f}", "" if c[3] < 0 else f"{c[3]:.2f}",
                                 a[0], f"{a[1]:.4f}" if a[0] else "", "yes" if sid in legacy else "no"]) + "\n")
    with open(f"{prefix}.summary.tsv", "w") as out:
        out.write("section\tmetric\tkey\tvalue\n")
        out.write(f"pairs\tn\ttotal\t{n_pairs}\n")
        out.write(f"pairs\tn\tkept\t{len(keep)}\n")
        for b in sorted(hist, key=lambda x: (x == '<0', x)):
            out.write(f"pairs\thist\t{b}\t{hist[b]}\n")
        out.write(f"members\tn\ttotal\t{len(ids)}\n")
        out.write(f"members\tn\textracted_here\t{len(legacy)}\n")
        out.write(f"members\tn\twith_ancestry\t{sum(1 for sid in ids if sid in anc)}\n")
    print(f"{n_pairs} pairs, {len(keep)} kept, {len(ids)} members", file=sys.stderr)
    EOF
  >>>

  output {
    File summary = "~{out_prefix}.summary.tsv"
    File members = "~{out_prefix}.members.tsv"
    File related_pairs = "~{out_prefix}.related_pairs.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
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
