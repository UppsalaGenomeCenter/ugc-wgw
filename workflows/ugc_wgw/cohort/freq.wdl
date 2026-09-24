version 1.0

# ugc-pacbio-wgw — cohort allele frequencies: sites-only VCF shards with bcftools +fill-tags (sex strata) and a summary
# origin: none
# see docs/DESIGN.md §7.7

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_bcftools_freq {
  meta {
    description: "Allele counts of one region shard of a cohort VCF: AC, AN, AF, NS, AC_Hom, AC_Het, AC_Hemi overall and per sex group (XX, XY), genotypes dropped. Optionally counts fully missing genotypes as hom-ref first (svx non-carriers) and adds custom tags."
    outputs: {
      shard_vcf: {
        description: "Sites-only shard VCF with the frequency tags"
      },
      shard_vcf_index: {
        description: "Index for the shard VCF"
      },
      vcf_samples: {
        description: "Sample names of the input VCF, one per line, sorted"
      }
    }
  }

  parameter_meta {
    vcf: {
      description: "Cohort VCF with per-sample genotypes (indexed)"
    }
    vcf_index: {
      description: "Index for vcf"
    }
    regions_bed: {
      description: "BED of the shard; a record belongs to the shard holding its POS (--regions-overlap 0)"
    }
    sample_ids: {
      description: "Cohort members; every one must be a sample of the VCF"
    }
    sample_sexes: {
      description: "XX, XY or empty per member, aligned with sample_ids; only XX and XY rows form the strata"
    }
    header_lines: {
      description: "##key=value lines added to the header (provenance of the frequency file)"
    }
    missing_to_homref: {
      description: "Count fully missing genotypes (./.) as hom-ref before tagging; for the svx merge, where a sample without the variant in its own call set is ./."
    }
    remove_tags: {
      description: "bcftools annotate --remove list applied first, e.g. INFO/AF,INFO/AC,INFO/AN to drop GLnexus's own values and header lines"
    }
    extra_tags: {
      description: "Extra +fill-tags expressions run after the strata, e.g. CF:1=INFO/SUPP/N_SAMPLES (kept separate: custom tags inside the stratified call get unstratified _XX/_XY copies)"
    }
    out_prefix: {
      description: "Output prefix"
    }
    threads: {
      description: "Threads"
    }
    mem_gb: {
      description: "Memory (GB); the pipeline streams one record at a time"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File vcf
    File vcf_index
    File regions_bed
    Array[String] sample_ids
    Array[String] sample_sexes
    Array[String] header_lines
    Boolean missing_to_homref = false
    String remove_tags = ""
    String extra_tags = ""
    String out_prefix
    Int threads = 2
    Int mem_gb = 4
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(vcf, "GB") + 20)

  command <<<
    set -euo pipefail
    export LC_ALL=C

    bcftools --version >&2

    mkdir -p in
    ln --symbolic "~{vcf}" "in/~{basename(vcf)}"
    ln --symbolic "~{vcf_index}" "in/~{basename(vcf_index)}"

    # every cohort member must be a sample of the VCF, by name: GLnexus orders samples lexicographically,
    # svx in cohort order, and +fill-tags only warns about a name it cannot find
    bcftools query --list-samples "in/~{basename(vcf)}" | sort > vcf_samples.txt
    sort "~{write_lines(sample_ids)}" | comm -23 - vcf_samples.txt > missing.txt
    if [ -s missing.txt ]; then
      echo "cohort members absent from ~{basename(vcf)}: $(tr '\n' ' ' < missing.txt)" >&2
      exit 1
    fi

    # strata: sample <TAB> XX|XY; members of unknown sex are counted in the unstratified tags only
    paste "~{write_lines(sample_ids)}" "~{write_lines(sample_sexes)}" \
      | awk -F'\t' -v OFS='\t' '$2 == "XX" || $2 == "XY" { print $1, $2 }' > groups.tsv
    # bcftools +fill-tags --samples-file (1.20 to 1.23) rejects sample names shorter than three characters
    # ("Could not parse the file", a bounds check in its parser); fail here with the names instead
    awk -F'\t' 'length($1) < 3 { short = short " " $1 }
      END { if (short != "") { print "sample names too short for bcftools +fill-tags --samples-file (three characters or more):" short > "/dev/stderr"; exit 1 } }' groups.tsv

    remove_args=()
    if [ -n "~{remove_tags}" ]; then remove_args=(--remove "~{remove_tags}"); fi
    group_args=()
    if [ -s groups.tsv ]; then group_args=(--samples-file groups.tsv); fi

    homref() {
      if [ "~{missing_to_homref}" = "true" ]; then
        bcftools +setGT --output-type u - -- --target-gt ./. --new-gt 0
      else
        cat
      fi
    }
    extra() {
      if [ -n "~{extra_tags}" ]; then
        bcftools +fill-tags --output-type u - -- --tags '~{extra_tags}'
      else
        cat
      fi
    }

    bcftools view \
      --threads ~{threads} \
      --regions-file "~{regions_bed}" \
      --regions-overlap 0 \
      --output-type u \
      "in/~{basename(vcf)}" \
    | bcftools annotate \
      ${remove_args[@]+"${remove_args[@]}"} \
      --header-lines "~{write_lines(header_lines)}" \
      --output-type u \
      - \
    | homref \
    | bcftools +fill-tags --output-type u - -- \
      ${group_args[@]+"${group_args[@]}"} \
      --tags AC,AN,AF,NS,AC_Hom,AC_Het,AC_Hemi \
    | extra \
    | bcftools view \
      --drop-genotypes \
      --output-type z \
      --output "~{out_prefix}.vcf.gz" \
      -
    bcftools index \
      --threads ~{threads} \
      --tbi \
      "~{out_prefix}.vcf.gz"
  >>>

  output {
    File shard_vcf = "~{out_prefix}.vcf.gz"
    File shard_vcf_index = "~{out_prefix}.vcf.gz.tbi"
    File vcf_samples = "vcf_samples.txt"
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

task ugc_wgw_freq_summary {
  meta {
    description: "Summary of the frequency files (records by type, contig class, filter, AF bin; singletons; call rate; strata) and the table of members with the sex used"
    outputs: {
      summary: {
        description: "Long-format TSV: resource, metric, key, value"
      },
      samples: {
        description: "Per member: group (XX, XY, unknown), sex source, presence in each VCF"
      }
    }
  }

  parameter_meta {
    small_variant_vcf: {
      description: "Sites-only small-variant frequency VCF (absent when the cohort has no joint VCF)"
    }
    small_variant_vcf_index: {
      description: "Index"
    }
    sv_vcf: {
      description: "Sites-only SV frequency VCF"
    }
    sv_vcf_index: {
      description: "Index"
    }
    small_variant_samples: {
      description: "Sample list of the small-variant input VCF (from the shard task)"
    }
    sv_samples: {
      description: "Sample list of the SV input VCF"
    }
    sample_ids: {
      description: "Cohort members"
    }
    sample_sexes: {
      description: "XX, XY or empty per member"
    }
    sample_sex_sources: {
      description: "sheet, inferred or empty per member"
    }
    out_prefix: {
      description: "Output prefix"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File? small_variant_vcf
    File? small_variant_vcf_index
    File sv_vcf
    File sv_vcf_index
    File? small_variant_samples
    File sv_samples
    Array[String] sample_ids
    Array[String] sample_sexes
    Array[String] sample_sex_sources
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  Boolean have_small_variants = defined(small_variant_vcf)
  Int threads = 1
  Int mem_gb = 2
  Int disk_size = ceil(size(sv_vcf, "GB") + size(small_variant_vcf, "GB") + 10)

  command <<<
    set -euo pipefail
    export LC_ALL=C

    bcftools --version >&2

    # the strata tags exist only when at least one member has a known sex
    has_tag() { bcftools view --header-only "$1" | grep -q "^##INFO=<ID=$2,"; }
    fmt_for() {
      local f='%CHROM\t%TYPE\t%FILTER\t%INFO/AC\t%INFO/AN\t%INFO/AF\t%INFO/NS' t
      for t in AC_XX AN_XX AC_XY AN_XY; do
        if has_tag "$1" "$t"; then f="$f\t%INFO/$t"; else f="$f\t."; fi
      done
      echo "$f"
    }

    if [ "~{have_small_variants}" = "true" ]; then
      bcftools query --format "$(fmt_for "~{default="" small_variant_vcf}")\n" "~{default="" small_variant_vcf}" > small.tsv
      cp "~{default="" small_variant_samples}" small_samples.txt
    else
      : > small.tsv
      : > small_samples.txt
    fi
    bcftools query --format "$(fmt_for "~{sv_vcf}")\t%INFO/SVTYPE\n" "~{sv_vcf}" > sv.tsv

    python3 - "~{out_prefix}" "~{have_small_variants}" small.tsv sv.tsv small_samples.txt "~{sv_samples}" \
      "~{write_lines(sample_ids)}" "~{write_lines(sample_sexes)}" "~{write_lines(sample_sex_sources)}" <<'PYEOF'
    import collections
    import sys

    (out_prefix, have_small, small_tsv, sv_tsv, small_samples, sv_samples,
     ids_file, sexes_file, sources_file) = sys.argv[1:10]

    def lines(path):
        try:
            with open(path) as fh:
                return [l.rstrip("\n") for l in fh]
        except OSError:
            return []

    ids = lines(ids_file)
    sexes = lines(sexes_file)
    sources = lines(sources_file)
    sexes += [""] * (len(ids) - len(sexes))
    sources += [""] * (len(ids) - len(sources))
    in_small = set(lines(small_samples)) if have_small == "true" else set()
    in_sv = set(lines(sv_samples))
    n = len(ids)
    rows = []

    def add(resource, metric, key, value):
        rows.append((resource, metric, key, value))

    add("samples", "n", "total", n)
    for g in ("XX", "XY"):
        add("samples", "n", g, sum(1 for s in sexes if s == g))
    add("samples", "n", "unknown", sum(1 for s in sexes if s not in ("XX", "XY")))
    for src in ("sheet", "inferred"):
        add("samples", "n_sex_source", src, sum(1 for s in sources if s == src))
    add("samples", "n_sex_source", "none", sum(1 for s in sources if s not in ("sheet", "inferred")))
    add("samples", "n_in_vcf", "small_variants", sum(1 for s in ids if s in in_small))
    add("samples", "n_in_vcf", "structural_variants", sum(1 for s in ids if s in in_sv))

    with open(out_prefix + ".samples.tsv", "w") as fh:
        fh.write("sample_id\tgroup\tsex_source\tin_small_variant_vcf\tin_sv_vcf\n")
        for sid, sex, src in zip(ids, sexes, sources):
            fh.write("\t".join([sid, sex if sex in ("XX", "XY") else "unknown", src or "none",
                                "yes" if sid in in_small else "no", "yes" if sid in in_sv else "no"]) + "\n")

    def contig_class(chrom):
        c = chrom[3:] if chrom.lower().startswith("chr") else chrom
        if c == "X":
            return "chrX"
        if c == "Y":
            return "chrY"
        if c in ("M", "MT"):
            return "chrM"
        return "autosome" if c.isdigit() else "other"

    def nums(field):
        out = []
        for x in field.split(","):
            try:
                out.append(float(x))
            except ValueError:
                pass
        return out

    def summarise(resource, path, sv):
        n_rec = 0
        by_type = collections.Counter()
        by_class = collections.Counter()
        by_filter = collections.Counter()
        bins = collections.Counter()
        strata = {"XX": 0, "XY": 0}
        ac_gt0 = singletons = low_ns = ns_n = 0
        ns_sum = 0.0
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 11:
                    continue
                chrom, vtype, filt, ac, an, af, ns, ac_xx, an_xx, ac_xy, an_xy = f[:11]
                if sv:
                    vtype = f[11] if len(f) > 11 else "."
                n_rec += 1
                by_type["MIXED" if "," in vtype else vtype] += 1
                by_class[contig_class(chrom)] += 1
                by_filter[filt] += 1
                total_ac = sum(nums(ac))
                afs = nums(af)
                if total_ac > 0:
                    ac_gt0 += 1
                if total_ac == 1:
                    singletons += 1
                if afs and total_ac > 0:
                    m = max(afs)
                    bins["<0.001" if m < 0.001 else "<0.01" if m < 0.01 else "<0.05" if m < 0.05 else ">=0.05"] += 1
                nsv = nums(ns)
                if nsv:
                    ns_sum += nsv[0]
                    ns_n += 1
                    if n and nsv[0] < 0.9 * n:
                        low_ns += 1
                if sum(nums(ac_xx)) > 0:
                    strata["XX"] += 1
                if sum(nums(ac_xy)) > 0:
                    strata["XY"] += 1
        add(resource, "records", "total", n_rec)
        for k, v in sorted(by_type.items()):
            add(resource, "records_by_type", k, v)
        for k, v in sorted(by_class.items()):
            add(resource, "records_by_contig_class", k, v)
        for k, v in sorted(by_filter.items()):
            add(resource, "records_by_filter", k, v)
        add(resource, "records", "ac_gt0", ac_gt0)
        add(resource, "records", "singletons", singletons)
        for k in ("<0.001", "<0.01", "<0.05", ">=0.05"):
            add(resource, "records_by_af_bin", k, bins[k])
        for g in ("XX", "XY"):
            add(resource, "records_ac_gt0_in", g, strata[g])
        add(resource, "mean_call_rate", "all", "%.4f" % (ns_sum / ns_n / n) if ns_n and n else "NA")
        add(resource, "records", "call_rate_below_0.9", low_ns)

    if have_small == "true":
        summarise("small_variants", small_tsv, False)
    summarise("structural_variants", sv_tsv, True)

    with open(out_prefix + ".summary.tsv", "w") as fh:
        fh.write("# ugc-pacbio-wgw cohort_freq: AF bins use the largest allele frequency of a record with AC > 0; "
                 "call rate = NS / n_samples; strata count records with AC > 0 in that group\n")
        fh.write("resource\tmetric\tkey\tvalue\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    PYEOF
  >>>

  output {
    File summary = "~{out_prefix}.summary.tsv"
    File samples = "~{out_prefix}.samples.tsv"
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
