version 1.0

# ugc-pacbio-wgw — mitorsaw's per-haplotype chrM genotypes as a standard diploid call set, and a haplogroup per haplotype
# origin: none (the collapse rule after nallo's collapse_mitochondrial_gt.py 1.0.0, MIT, genomic-medicine-sweden/nallo)
# see docs/DESIGN.md §5, §8.1, §18 (2026-11-05); docs/ENTRYPOINTS.md §4, §7

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_mito_collapse {
  meta {
    description: "mitorsaw writes one GT entry per detected mitochondrial haplotype (0|1|1|1), which no VCF consumer can merge or count. This task collapses GT to diploid by nallo's rule (every haplotype reference -> 0/0, every haplotype alternate -> a/a, mixed -> 0/a: heteroplasmic; the read fraction of the allele is mitorsaw's own FORMAT/VAF), splits the records into small variants and SVs at mito_sv_min_size, keeps DP, AD and VAF, records the rule and the haplotypes' abundances in the header, and writes a third VCF with one haploid column per haplotype for the haplogroup classification. A VCF without records (no chrM reads) gives empty, valid files and a message."
    outputs: {
      mt_small_variant_vcf: {
        description: "<out_prefix>.mt.small_variants.vcf.gz: the collapsed records whose alleles differ by less than mito_sv_min_size"
      }
      mt_small_variant_vcf_index: {
        description: "Index"
      }
      mt_sv_vcf: {
        description: "<out_prefix>.mt.svs.vcf.gz: the collapsed records with an allele-length difference at or above mito_sv_min_size, or a symbolic ALT"
      }
      mt_sv_vcf_index: {
        description: "Index"
      }
      mt_haplotypes_vcf: {
        description: "<out_prefix>.mt.haplotypes.vcf.gz: every record, one haploid GT column per mitorsaw haplotype named <vcf sample>.<haplotype label>"
      }
      mt_haplotypes_vcf_index: {
        description: "Index"
      }
      mt_haplotypes_tsv: {
        description: "<out_prefix>.mt.haplotypes.tsv: haplotype, estimated abundance, variants carried"
      }
      mt_summary: {
        description: "<out_prefix>.mt.summary.tsv: one row: haplotypes, records, homoplasmic, heteroplasmic, SVs, mean and max VAF of the heteroplasmic sites, the major haplotype and its abundance"
      }
      msg: {
        description: "Messages: no chrM variant, a haplotype count that differs between the JSON and the GT, records with more than two alleles"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID (written to the summary; the VCF's own sample name is kept in the VCFs)"
    }
    mitorsaw_vcf: {
      description: "mitorsaw haplotype's VCF (GT per haplotype, FORMAT GT:DP:AD:VAF)"
    }
    mitorsaw_vcf_index: {
      description: "Index for mitorsaw_vcf"
    }
    mitorsaw_hap_stats: {
      description: "mitorsaw haplotype's JSON (haplotype labels and estimated abundances, in GT order)"
    }
    mito_sv_min_size: {
      description: "Allele-length difference (bp) from which a record is an SV (nallo's 50)"
    }
    out_prefix: {
      description: "Output prefix"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    File mitorsaw_vcf
    File mitorsaw_vcf_index
    File mitorsaw_hap_stats
    Int mito_sv_min_size = 50
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  command <<<
    set -euo pipefail

    python3 --version >&2
    bcftools --version | head -1 >&2

    python3 - "~{mitorsaw_vcf}" "~{mitorsaw_hap_stats}" "~{sample_id}" "~{out_prefix}" ~{mito_sv_min_size} <<'EOF'
    # ---- mito_collapse.py (begin) ----
    import gzip, json, sys
    vcf_path, json_path, sample_id, prefix, sv_min = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
    haps = json.load(open(json_path)).get("haplotypes", [])
    labels = [str(h.get("label", f"hap_{i}")) for i, h in enumerate(haps)]
    abund = [float(h.get("estimated_abundance") or 0.0) for h in haps]
    messages, header, records = [], [], []
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line:
                (header if line.startswith("#") else records).append(line)
    cols = header[-1].split("\t") if header and header[-1].startswith("#CHROM") else []
    if len(cols) < 10:
        sys.exit(f"{vcf_path}: no sample column")
    vcf_sample = cols[9]
    mito_version = next((l.split("=", 1)[1].strip('"') for l in header if l.startswith("##mitorsaw_version=")), "")

    def entries_of(record):
        return record.split("\t")[9].split(":")[0].replace("|", "/").split("/")

    width = None
    for r in records:
        n = len(entries_of(r))
        width = n if width is None else width
        if n != width:
            messages.append(f"mitorsaw: GT width varies between records ({width} and {n}); the haplotype columns follow the first record")
            break
    if width is not None and width != len(labels):
        messages.append(f"mitorsaw: {len(labels)} haplotypes in the JSON but {width} GT entries per record; the haplotype columns follow the GT")
    n_hap = width if width is not None else max(len(labels), 1)
    labels = (labels + [f"hap_{i}" for i in range(len(labels), n_hap)])[:n_hap]
    abund = (abund + [0.0] * n_hap)[:n_hap]

    def collapse(gt):
        """nallo's rule extended to several ALT alleles: the alleles present across the haplotypes, 0 first, then by
        the abundance of the haplotypes carrying them; one allele -> a/a, two -> a/b, more -> the two most abundant."""
        entries = gt.replace("|", "/").split("/")
        known = [(e, i) for i, e in enumerate(entries) if e != "."]
        if not known:
            return "./.", entries, False
        weight = {}
        for e, i in known:
            weight[e] = weight.get(e, 0.0) + (abund[i] if i < len(abund) else 0.0)
        alleles = sorted(weight, key=lambda a: (a != "0", -weight[a], int(a) if a.isdigit() else 0))
        if len(alleles) == 1:
            return f"{alleles[0]}/{alleles[0]}", entries, False
        return f"{alleles[0]}/{alleles[1]}", entries, len(alleles) > 2

    rule = ("every haplotype reference -> 0/0; every haplotype alternate -> a/a; mixed -> 0/a (heteroplasmic), after nallo "
            "collapse_mitochondrial_gt.py 1.0.0; FORMAT/VAF is the read fraction of each ALT allele")
    extra = [f"##ugc_wgw_mito_collapse_rule={rule}", f"##ugc_wgw_mito_haplotypes={n_hap}",
             "##ugc_wgw_mito_abundances=" + ",".join(f"{l}:{a:.4f}" for l, a in zip(labels, abund)),
             f"##ugc_wgw_mito_sv_min_size={sv_min}", f"##ugc_wgw_mito_source=mitorsaw {mito_version}".rstrip()]
    head_common = header[:-1] + extra
    small, svs, hap_rows = [], [], []
    n_hom = n_het = n_sv = n_more = 0
    het_vafs = []
    carried = [0] * n_hap
    for r in records:
        f = r.split("\t")
        ref, alts = f[3], f[4].split(",")
        keys, vals = f[8].split(":"), f[9].split(":")
        if "GT" not in keys:
            sys.exit(f"{vcf_path}: a record without GT at {f[0]}:{f[1]}")
        gt_i = keys.index("GT")
        new_gt, entries, more = collapse(vals[gt_i])
        n_more += 1 if more else 0
        vals[gt_i] = new_gt
        f[9] = ":".join(vals)
        is_sv = any(a.startswith("<") or abs(len(a) - len(ref)) >= sv_min for a in alts)
        (svs if is_sv else small).append("\t".join(f))
        a, b = new_gt.split("/")
        if is_sv:
            n_sv += 1
        elif a == b and a not in (".", "0"):
            n_hom += 1
        elif a != b and "." not in (a, b):
            n_het += 1
            k = b if a == "0" else a
            if "VAF" in keys and k.isdigit() and int(k) >= 1:
                parts = vals[keys.index("VAF")].split(",")
                try:
                    het_vafs.append(float(parts[int(k) - 1]))
                except (ValueError, IndexError):
                    pass
        padded = (entries + ["."] * n_hap)[:n_hap]
        for i, e in enumerate(padded):
            if e not in ("0", "."):
                carried[i] += 1
        hap_rows.append("\t".join(f[:8] + ["GT"] + padded))
    if n_more:
        messages.append(f"mitorsaw: {n_more} record(s) with more than two alleles across the haplotypes; the two most abundant were kept")
    if not records:
        messages.append("mitorsaw: no chrM variant called (no chrM reads?)")
    with open(f"{prefix}.mt.small_variants.vcf", "w") as out:
        out.write("\n".join(head_common + [header[-1]] + small) + "\n")
    with open(f"{prefix}.mt.svs.vcf", "w") as out:
        out.write("\n".join(head_common + [header[-1]] + svs) + "\n")
    hap_header = ([l for l in head_common if not l.startswith("##FORMAT=")]
                  + ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                     "\t".join(cols[:9] + [f"{vcf_sample}.{l}" for l in labels])])
    with open(f"{prefix}.mt.haplotypes.vcf", "w") as out:
        out.write("\n".join(hap_header + hap_rows) + "\n")
    with open(f"{prefix}.mt.haplotypes.tsv", "w") as out:
        out.write("haplotype\tabundance\tn_variants\n")
        for l, a, c in zip(labels, abund, carried):
            out.write(f"{l}\t{a:.4f}\t{c}\n")
    major = max(range(n_hap), key=lambda i: abund[i])
    with open(f"{prefix}.mt.summary.tsv", "w") as out:
        out.write("sample_id\tvcf_sample\thaplotypes\trecords\thomoplasmic\theteroplasmic\tsvs\thet_vaf_mean\thet_vaf_max\tmajor_haplotype\tmajor_abundance\n")
        out.write("\t".join([sample_id, vcf_sample, str(n_hap), str(len(records)), str(n_hom), str(n_het), str(n_sv),
                             f"{sum(het_vafs) / len(het_vafs):.4f}" if het_vafs else "", f"{max(het_vafs):.4f}" if het_vafs else "",
                             labels[major], f"{abund[major]:.4f}"]) + "\n")
    with open("messages.txt", "w") as out:
        out.write("".join(m + "\n" for m in messages))
    print(f"{len(records)} records: {n_hom} homoplasmic, {n_het} heteroplasmic, {n_sv} SVs; {n_hap} haplotypes", file=sys.stderr)
    # ---- mito_collapse.py (end) ----
    EOF

    for name in mt.small_variants mt.svs mt.haplotypes; do
      bcftools view --output-type z --output "~{out_prefix}.${name}.vcf.gz" "~{out_prefix}.${name}.vcf"
      bcftools index --tbi "~{out_prefix}.${name}.vcf.gz"
      rm "~{out_prefix}.${name}.vcf"
    done
  >>>

  output {
    File mt_small_variant_vcf = "~{out_prefix}.mt.small_variants.vcf.gz"
    File mt_small_variant_vcf_index = "~{out_prefix}.mt.small_variants.vcf.gz.tbi"
    File mt_sv_vcf = "~{out_prefix}.mt.svs.vcf.gz"
    File mt_sv_vcf_index = "~{out_prefix}.mt.svs.vcf.gz.tbi"
    File mt_haplotypes_vcf = "~{out_prefix}.mt.haplotypes.vcf.gz"
    File mt_haplotypes_vcf_index = "~{out_prefix}.mt.haplotypes.vcf.gz.tbi"
    File mt_haplotypes_tsv = "~{out_prefix}.mt.haplotypes.tsv"
    File mt_summary = "~{out_prefix}.mt.summary.tsv"
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

task ugc_wgw_haplogrep {
  meta {
    description: "A mitochondrial haplogroup for each mitorsaw haplotype: haplogrep3 classifies every column of the per-haplotype VCF against a phylogenetic tree installed in the image, and the result is joined with the haplotypes' abundances. The sample's haplogroup is its most abundant haplotype's. A VCF without records gives an empty table and no haplogroup."
    outputs: {
      mt_haplogroups: {
        description: "<out_prefix>.mt.haplogroups.tsv: haplotype, abundance, haplogroup, quality, range, found_polys, remaining_polys (haplogrep3's columns of the same names)"
      }
      mt_haplogroup: {
        description: "The haplogroup of the most abundant haplotype; empty when nothing was classified"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID"
    }
    mt_haplotypes_vcf: {
      description: "The per-haplotype VCF of ugc_wgw_mito_collapse"
    }
    mt_haplotypes_tsv: {
      description: "The haplotypes' abundances of ugc_wgw_mito_collapse"
    }
    tree: {
      description: "haplogrep3 tree id; the image has phylotree-fu-rcrs@1.3 (PhyloTree 17 forensic update, the default) and phylotree-rcrs@17.3"
    }
    out_prefix: {
      description: "Output prefix"
    }
    ugc_wgw_container_registry: {
      description: "Registry holding the ugc-built images (haplogrep3)"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    File mt_haplotypes_vcf
    File mt_haplotypes_tsv
    String tree = "phylotree-fu-rcrs@1.3"
    String out_prefix
    String ugc_wgw_container_registry
    RuntimeAttributes runtime_attributes
  }

  command <<<
    set -euo pipefail

    # haplogrep3 prints its banner and exits 2 on --version: the version line must not fail the pipeline
    haplogrep3 --version 2>&1 | grep -i 'haplogrep 3' >&2 || true
    gzip -dc "~{mt_haplotypes_vcf}" > haplotypes.vcf

    if grep -qv '^#' haplotypes.vcf; then
      haplogrep3 classify --in haplotypes.vcf --tree "~{tree}" --out classified.txt --hits 1 --extend-report >&2
    else
      printf '"SampleID"\t"Haplogroup"\t"Rank"\t"Quality"\t"Range"\t"Not_Found_Polys"\t"Found_Polys"\t"Remaining_Polys"\t"AAC_In_Remainings"\t"Input_Sample"\n' > classified.txt
      echo "no chrM variant in ~{sample_id}: nothing to classify" >&2
    fi

    # join on the haplotype label (the column name's part after the last dot) with the abundances; quotes stripped
    awk -F'\t' -v OFS='\t' '
      NR == FNR { if (FNR > 1) { ab[$1] = $2; order[++n] = $1 } next }
      FNR == 1 { print "haplotype", "abundance", "haplogroup", "quality", "range", "found_polys", "remaining_polys"; next }
      { gsub(/"/, ""); m = split($1, p, "."); h = p[m]; print h, ab[h], $2, $4, $5, $7, $8 }
    ' "~{mt_haplotypes_tsv}" classified.txt > "~{out_prefix}.mt.haplogroups.tsv"

    # the sample's haplogroup: the row of the most abundant haplotype
    sort -t $'\t' -k2,2gr "~{mt_haplotypes_tsv}" | awk -F'\t' 'NR == 1 && $1 != "haplotype" { print $1 }' > major.txt
    awk -F'\t' -v major="$(cat major.txt)" 'NR > 1 && $1 == major { print $3; exit }' "~{out_prefix}.mt.haplogroups.tsv" > haplogroup.txt
    [ -s haplogroup.txt ] || : > haplogroup.txt
  >>>

  output {
    File mt_haplogroups = "~{out_prefix}.mt.haplogroups.tsv"
    String mt_haplogroup = read_string("haplogroup.txt")
  }

  runtime {
    docker: "~{ugc_wgw_container_registry}/haplogrep3@sha256:12de1423817083013a974595cabb4ca8aadbd97a3109fafe4e1ddf10a973c09e"  # haplogrep3 3.3.2
    cpu: 1
    memory: "2 GiB"
    disk: "5 GB"
    disks: "local-disk 5 HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}
