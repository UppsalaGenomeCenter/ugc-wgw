version 1.0
# ugc-pacbio-wgw — cohort_graph_gt's tasks: the graph's SV alleles as a site list, a pileup index and the genotypes per member (kanpig), the members' columns joined per region shard, the gather
# origin: none
# see docs/DESIGN.md §5, §18 (2026-11-19); docs/ENTRYPOINTS.md §9.6; docs/expansion/pangenome-genotyping.design.md
import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_graph_gt_sites {
  meta {
    description: "The site list: the graph's vcfwave VCF split to biallelic records and left-aligned (bcftools norm), the alleles that change the length by size_min bp or more kept, each with a stable ID, its type and length, and the allele count of the graph's own members; no genotype columns. Alleles longer than size_max stay in the list marked SKIP=too_long (the genotyper leaves them uncalled). The graph members' genotypes at the same sites go to a second, indexed table for the comparison with the reads."
    outputs: {
      sites: {
        description: "<out_prefix>.sites.vcf.gz: one record per SV allele, INFO SVTYPE, SVLEN, GRAPH_AC, GRAPH_AN, GRAPH_NS, IN_TR, IN_EXCL, SKIP"
      },
      sites_index: {
        description: "Its index"
      },
      graph_calls: {
        description: "<out_prefix>.graph_calls.tsv.gz: chrom, pos, id, then per graph member alt:called (alternate alleles : called alleles)"
      },
      graph_calls_index: {
        description: "Its tabix index"
      },
      skipped: {
        description: "<out_prefix>.skipped.tsv: id, chrom, pos, svtype, length, reason, for every allele in the list that is not genotyped or was left out of it"
      },
      counts: {
        description: "<out_prefix>.sites_counts.tsv: metric, key, value: records read, alleles kept by type and size class, left out by reason, the graph's members"
      }
    }
  }

  parameter_meta {
    graph_vcf: {
      description: "cohort_pangenome's pangenome_wave_vcf"
    }
    graph_vcf_index: {
      description: "Its index"
    }
    reference: {
      description: "Reference FASTA (bcftools norm left-aligns against it)"
    }
    reference_index: {
      description: "Its index"
    }
    tandem_repeat_bed: {
      description: "Tandem-repeat regions (the TRGT catalogue): an allele that starts in one is marked IN_TR, where a repeat allele from an assembly and the same allele in a read alignment are often written differently"
    }
    excluded_regions_bed: {
      description: "Regions where SV calls from read alignments are not trusted (sawfish's exclusion list: gaps, centromeres, common copy-number regions): an allele that starts in one is marked IN_EXCL"
    }
    size_min: {
      description: "Smallest length change (bp) of an allele in the list"
    }
    size_max: {
      description: "Alleles longer than this are listed and not genotyped"
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
    File graph_vcf
    File graph_vcf_index
    File reference
    File reference_index
    File? tandem_repeat_bed
    File? excluded_regions_bed
    Int size_min
    Int size_max
    String out_prefix
    Int mem_gb = 4
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(graph_vcf, "GB") * 6 + size(reference, "GB") + 10)

  command <<<
    set -euo pipefail

    bcftools --version | head -1 >&2
    ln -s "~{graph_vcf}" graph.vcf.gz
    ln -s "~{graph_vcf_index}" graph.vcf.gz.tbi

    # the split records are read as a stream: the graph VCF of a large cohort is too big to write out uncompressed
    mkfifo norm.fifo
    bcftools norm --multiallelics -any --fasta-ref "~{reference}" --output-type v graph.vcf.gz > norm.fifo 2> norm.log &
    norm_pid=$!

    python3 - norm.fifo "~{size_min}" "~{size_max}" "~{out_prefix}" "~{select_first([tandem_repeat_bed, ""])}" "~{select_first([excluded_regions_bed, ""])}" <<'EOF'
    # ---- graph_gt_sites.py (begin) ----
    import bisect
    import gzip
    import hashlib
    import re
    import sys

    source, size_min, size_max, prefix = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    tr_bed = sys.argv[5] if len(sys.argv) > 5 else ""
    excluded_bed = sys.argv[6] if len(sys.argv) > 6 else ""
    SPLIT = re.compile(r"[|/]")


    def load_regions(path):
        """A BED as chrom -> (sorted starts, the furthest end among the regions up to each start); empty without a path."""
        if not path:
            return {}
        regions = {}
        with (gzip.open if path.endswith(".gz") else open)(path, "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 3 and f[1].isdigit() and f[2].isdigit():
                    regions.setdefault(f[0], []).append((int(f[1]), int(f[2])))
        out = {}
        for chrom, ivs in regions.items():
            ivs.sort()
            ends, far = [], 0
            for _, e in ivs:
                far = max(far, e)
                ends.append(far)
            out[chrom] = ([s for s, _ in ivs], ends)
        return out


    def inside(regions, chrom, pos):
        """A 1-based position inside a region of the BED (0-based, half-open)."""
        if chrom not in regions:
            return False
        starts, ends = regions[chrom]
        i = bisect.bisect_right(starts, pos - 1) - 1
        return i >= 0 and ends[i] > pos - 1


    tandem_repeats, excluded = load_regions(tr_bed), load_regions(excluded_bed)


    def size_class(n):
        return "50-999" if n < 1000 else ("1k-9.9k" if n < 10000 else "10k+")


    def allele_counts(field):
        """(alternate alleles, called alleles) of one genotype of a biallelic record."""
        parts = SPLIT.split(field.split(":", 1)[0])
        return sum(1 for a in parts if a == "1"), sum(1 for a in parts if a != ".")


    counts = {}


    def count(metric, key, n=1):
        counts[(metric, key)] = counts.get((metric, key), 0) + n


    sites = open(prefix + ".sites.vcf", "w")
    calls = open(prefix + ".graph_calls.tsv", "w")
    skipped = open(prefix + ".skipped.tsv", "w")
    skipped.write("id\tchrom\tpos\tsvtype\tlength\treason\n")
    samples = []
    pending = {}      # (ref, alt) -> [chrom, pos, per-sample [alt, called]] for the position being read
    pending_at = None


    def flush():
        for (ref, alt), (chrom, pos, per_sample) in sorted(pending.items()):
            d = len(alt) - len(ref)
            svtype, length = ("INS" if d > 0 else "DEL"), abs(d)
            sid = "%s-%s-%s-%d-%s" % (chrom, pos, svtype, length, hashlib.sha1((ref + ">" + alt).encode()).hexdigest()[:8])
            ac = sum(a for a, _ in per_sample)
            an = sum(c for _, c in per_sample)
            ns = sum(1 for _, c in per_sample if c)
            info = "SVTYPE=%s;SVLEN=%d;GRAPH_AC=%d;GRAPH_AN=%d;GRAPH_NS=%d" % (svtype, d, ac, an, ns)
            cls = size_class(length)
            if inside(excluded, chrom, int(pos)):
                info += ";IN_EXCL"
                count("sites_in_excluded_regions", svtype)
            if inside(tandem_repeats, chrom, int(pos)):
                info += ";IN_TR"
                count("sites_in_tandem_repeats", svtype)
            if length > size_max:
                info += ";SKIP=too_long"
                skipped.write("\t".join((sid, chrom, pos, svtype, str(length), "too_long")) + "\n")
                count("not_genotyped", "too_long")
            count("sites", svtype + ":" + cls)
            sites.write("\t".join((chrom, pos, sid, ref, alt, ".", "PASS", info)) + "\n")
            calls.write("\t".join([chrom, pos, sid] + ["%d:%d" % (a, c) for a, c in per_sample]) + "\n")
        pending.clear()


    for line in open(source):
        if line.startswith("##"):
            if line.startswith(("##fileformat", "##contig")):
                sites.write(line)
            continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#"):
            samples = f[9:]
            for s in samples:
                count("graph_sample", s)
            for tag, number, typ, text in (
                    ("SVTYPE", "1", "String", "INS or DEL: the allele is longer or shorter than the reference"),
                    ("SVLEN", "1", "Integer", "Length of the alternate allele minus the reference allele"),
                    ("GRAPH_AC", "1", "Integer", "Alternate alleles among the graph's own members (from their assemblies)"),
                    ("GRAPH_AN", "1", "Integer", "Called alleles among the graph's own members"),
                    ("GRAPH_NS", "1", "Integer", "Graph members with a called allele"),
                    ("IN_TR", "0", "Flag", "The allele starts inside a tandem repeat of the TRGT catalogue"),
                    ("IN_EXCL", "0", "Flag", "The allele starts inside a region the SV caller's exclusion list names (gaps, centromeres, common copy-number regions)"),
                    ("SKIP", "1", "String", "Not genotyped from the reads: too_long (longer than size_max)")):
                sites.write('##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % (tag, number, typ, text))
            sites.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            calls.write("#chrom\tpos\tid\t" + "\t".join(samples) + "\n")
            continue
        count("records", "read")
        chrom, pos, ref, alt = f[0], f[1], f[3], f[4]
        if (chrom, pos) != pending_at:
            flush()
            pending_at = (chrom, pos)
        if alt.startswith("<") or alt in ("*", ".") or "[" in alt or "]" in alt:
            if alt.startswith("<"):
                skipped.write("\t".join((".", chrom, pos, alt.strip("<>"), ".", "symbolic")) + "\n")
                count("left_out", "symbolic")
            continue
        d = abs(len(alt) - len(ref))
        if d < size_min:
            if min(len(ref), len(alt)) >= size_min:   # an inversion or a replacement: long on both sides, no length change
                skipped.write("\t".join((".", chrom, pos, "COMPLEX", str(max(len(ref), len(alt))), "no_length_change")) + "\n")
                count("left_out", "no_length_change")
            else:
                count("left_out", "below_size_min")
            continue
        per_sample = [list(allele_counts(g)) for g in f[9:]]
        if (ref, alt) in pending:   # the same allele from two records of the graph VCF: one site, the genotypes added
            count("left_out", "duplicate_merged")
            for mine, theirs in zip(pending[(ref, alt)][2], per_sample):
                mine[1] = max(mine[1], theirs[1])
                mine[0] = min(mine[1], mine[0] + theirs[0])
        else:
            pending[(ref, alt)] = [chrom, pos, per_sample]
    flush()
    for f in (sites, calls, skipped):
        f.close()
    with open(prefix + ".sites_counts.tsv", "w") as out:
        out.write("metric\tkey\tvalue\n")
        for (metric, key), n in sorted(counts.items()):
            out.write("%s\t%s\t%d\n" % (metric, key, n))
    # ---- graph_gt_sites.py (end) ----
    EOF
    wait "$norm_pid" || { tail -n 20 norm.log >&2; exit 1; }
    tail -n 3 norm.log >&2
    # bcftools norm keeps its output sorted; the site VCF goes through sort once more because left-alignment may
    # have moved a record across a neighbour that this script wrote in one position group
    bcftools sort --output-type z --output "~{out_prefix}.sites.vcf.gz" "~{out_prefix}.sites.vcf" 2> /dev/null
    tabix --preset vcf "~{out_prefix}.sites.vcf.gz"
    (head -n 1 "~{out_prefix}.graph_calls.tsv"; tail -n +2 "~{out_prefix}.graph_calls.tsv" | LC_ALL=C sort -k1,1 -k2,2n) | bgzip > "~{out_prefix}.graph_calls.tsv.gz"
    tabix --sequence 1 --begin 2 --end 2 --comment '#' "~{out_prefix}.graph_calls.tsv.gz"
    rm "~{out_prefix}.sites.vcf" "~{out_prefix}.graph_calls.tsv"
    cat "~{out_prefix}.sites_counts.tsv" >&2
  >>>

  output {
    File sites = "~{out_prefix}.sites.vcf.gz"
    File sites_index = "~{out_prefix}.sites.vcf.gz.tbi"
    File graph_calls = "~{out_prefix}.graph_calls.tsv.gz"
    File graph_calls_index = "~{out_prefix}.graph_calls.tsv.gz.tbi"
    File skipped = "~{out_prefix}.skipped.tsv"
    File counts = "~{out_prefix}.sites_counts.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
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

task ugc_wgw_kanpig_plup {
  meta {
    description: "kanpig plup: a member's alignment as a pileup index, every insertion and deletion of size_min to size_max bp its reads carry, sorted, compressed and indexed. It depends on the member's BAM alone, so the call cache keeps it: a second graph genotypes the member from this file without reading the BAM again."
    outputs: {
      plup: {
        description: "<sample_id>.<ref_name>.plup.gz"
      },
      plup_index: {
        description: "Its tabix index"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID"
    }
    bam: {
      description: "The member's alignment to the reference (the haplotagged BAM of singleton or downstream)"
    }
    bam_index: {
      description: "Its index"
    }
    ref_name: {
      description: "Reference name, for the file name"
    }
    size_min: {
      description: "Smallest variant indexed"
    }
    size_max: {
      description: "Largest variant indexed"
    }
    threads: {
      description: "Threads (kanpig says it is limited by reading: little gain past 4 to 8)"
    }
    mem_gb: {
      description: "Memory (GB); kanpig's own figure is about 2 GB per thread"
    }
    ugc_wgw_container_registry: {
      description: "Registry holding the ugc-pacbio-wgw images"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    File bam
    File bam_index
    String ref_name
    Int size_min
    Int size_max
    Int threads = 8
    Int mem_gb = 16
    String ugc_wgw_container_registry
    RuntimeAttributes runtime_attributes
  }

  String out = "~{sample_id}.~{ref_name}.plup.gz"
  Int disk_size = ceil(size(bam, "GB") + 20)

  command <<<
    set -euo pipefail

    kanpig --version >&2
    ln -s "~{bam}" reads.bam
    ln -s "~{bam_index}" reads.bam.bai

    kanpig plup --bam reads.bam --threads ~{threads} --sizemin ~{size_min} --sizemax ~{size_max} --output plup.raw 2> plup.log \
      || { tail -n 20 plup.log >&2; exit 1; }
    grep -E 'processed|regions' plup.log >&2 || true

    # kanpig writes the index unsorted and uncompressed: header lines first, then by position
    mkdir -p sorttmp
    (awk '/^#/' plup.raw; awk '!/^#/' plup.raw | LC_ALL=C sort -k1,1 -k2,2n -S 2G -T sorttmp) | bgzip --threads ~{threads} > "~{out}"
    tabix --preset bed "~{out}"
    rm -r plup.raw sorttmp
  >>>

  output {
    File plup = "~{out}"
    File plup_index = "~{out}.tbi"
  }

  runtime {
    docker: "~{ugc_wgw_container_registry}/kanpig@sha256:59504c7d557c97bba4e17fcd4cdfbd2b98445c1e6d6c286e81ab122ca4fad92d"  # kanpig 2.0.2 on pb_wdl_base
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

task ugc_wgw_kanpig_gt {
  meta {
    description: "kanpig gt: one member's genotype at every site of the list, from its pileup index. Ploidy follows the member's sex through the reference bundle's expected-copy-number BED (chrX outside the pseudoautosomal regions is haploid in an XY member, chrY is absent in an XX member); a member of unknown sex is genotyped as diploid everywhere. Sites outside the size bounds are written uncalled. The VCF keeps every field kanpig writes."
    outputs: {
      vcf: {
        description: "<out_prefix>.vcf.gz: the member's genotypes at the sites, sorted"
      },
      vcf_index: {
        description: "Its index"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID; the VCF's sample column"
    }
    sex: {
      description: "XX, XY or empty"
    }
    plup: {
      description: "The member's pileup index"
    }
    plup_index: {
      description: "Its tabix index"
    }
    sites: {
      description: "The site list"
    }
    sites_index: {
      description: "Its index"
    }
    reference: {
      description: "Reference FASTA"
    }
    reference_index: {
      description: "Its index"
    }
    expected_cn_xx: {
      description: "Expected copy number per region for an XX sample (chrom, start, end, name, copies)"
    }
    expected_cn_xy: {
      description: "The same for an XY sample"
    }
    size_min: {
      description: "Smallest variant genotyped"
    }
    size_max: {
      description: "Largest variant genotyped"
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
    ugc_wgw_container_registry: {
      description: "Registry holding the ugc-pacbio-wgw images"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    String sex
    File plup
    File plup_index
    File sites
    File sites_index
    File reference
    File reference_index
    File expected_cn_xx
    File expected_cn_xy
    Int size_min
    Int size_max
    String out_prefix
    Int threads = 4
    Int mem_gb = 8
    String ugc_wgw_container_registry
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(plup, "GB") * 2 + size(sites, "GB") * 20 + size(reference, "GB") + 10)

  command <<<
    set -euo pipefail

    kanpig --version >&2
    ln -s "~{plup}" reads.plup.gz
    ln -s "~{plup_index}" reads.plup.gz.tbi
    ln -s "~{sites}" sites.vcf.gz
    ln -s "~{sites_index}" sites.vcf.gz.tbi
    ln -s "~{reference}" ref.fa
    ln -s "~{reference_index}" ref.fa.fai

    # kanpig's ploidy BED lists the regions that are not diploid: chrom, start, end, copies
    ploidy=()
    sex="~{sex}"
    case "$sex" in
      XX) awk -F'\t' 'NF >= 5 && $5 != 2 {print $1 "\t" $2 "\t" $3 "\t" $5}' "~{expected_cn_xx}" > ploidy.bed ;;
      XY) awk -F'\t' 'NF >= 5 && $5 != 2 {print $1 "\t" $2 "\t" $3 "\t" $5}' "~{expected_cn_xy}" > ploidy.bed ;;
      *) echo "sex unknown: genotyped as diploid everywhere" >&2 ;;
    esac
    if [ -s ploidy.bed ]; then
      ploidy=(--ploidy-bed ploidy.bed)
      cat ploidy.bed >&2
    fi

    kanpig gt --input sites.vcf.gz --reads reads.plup.gz --reference ref.fa --sample "~{sample_id}" \
      --threads ~{threads} --sizemin ~{size_min} --sizemax ~{size_max} "${ploidy[@]}" --out gt.raw.vcf 2> gt.log \
      || { tail -n 20 gt.log >&2; exit 1; }
    grep -E 'variants|genotype counts|Non:|Ref:|Hom:|Het:' gt.log >&2 || true

    mkdir -p sorttmp
    bcftools sort --temp-dir sorttmp --output-type z --output "~{out_prefix}.vcf.gz" gt.raw.vcf 2> /dev/null
    tabix --preset vcf "~{out_prefix}.vcf.gz"
    rm -r gt.raw.vcf sorttmp
  >>>

  output {
    File vcf = "~{out_prefix}.vcf.gz"
    File vcf_index = "~{out_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "~{ugc_wgw_container_registry}/kanpig@sha256:59504c7d557c97bba4e17fcd4cdfbd2b98445c1e6d6c286e81ab122ca4fad92d"  # kanpig 2.0.2 on pb_wdl_base
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

task ugc_wgw_graph_gt_join {
  meta {
    description: "One region shard of the cohort VCF: the sites that start in the shard, a column per cohort member in cohort order read from the members' VCFs one at a time (genotype unphased, GQ, DP, AD; a member without reads is uncalled everywhere), allele counts and frequencies over the cohort and per sex group beside the graph's own; per member the counts of its genotypes; and for members that are also in the graph, the genotype from the reads against the genotype from the assembly."
    outputs: {
      vcf: {
        description: "<out_prefix>.vcf.gz: the shard's sites with every member's genotype"
      },
      members: {
        description: "<out_prefix>.members.tsv: per member and SV type the genotype counts of the shard"
      },
      concordance: {
        description: "<out_prefix>.concordance.tsv: per graph member, SV type, size class and context (excluded_region, tandem_repeat or other) the comparison counts of the shard"
      },
      site_counts: {
        description: "<out_prefix>.site_counts.tsv: per SV type and size class the shard's sites, genotyped, variable in the cohort, carried in the graph only"
      }
    }
  }

  parameter_meta {
    members_bundle: {
      description: "ugc_wgw_files_bundle tar of the members' VCFs followed by their indexes, both in the order of read_sample_ids"
    }
    read_sample_ids: {
      description: "The members that were genotyped (those with an alignment), in the bundle's order"
    }
    sample_ids: {
      description: "The cohort's members, in order: the VCF's columns"
    }
    sample_sexes: {
      description: "XX, XY or empty per member"
    }
    sites: {
      description: "The site list"
    }
    sites_index: {
      description: "Its index"
    }
    graph_calls: {
      description: "The graph members' genotypes at the sites"
    }
    graph_calls_index: {
      description: "Its tabix index"
    }
    regions_bed: {
      description: "The shard's regions"
    }
    out_prefix: {
      description: "Output prefix"
    }
    mem_gb: {
      description: "Memory (GB): the shard's sites times the members, a few bytes each"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File members_bundle
    Array[String] read_sample_ids
    Array[String] sample_ids
    Array[String] sample_sexes
    File sites
    File sites_index
    File graph_calls
    File graph_calls_index
    File regions_bed
    String out_prefix
    Int mem_gb = 8
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(members_bundle, "GB") * 2.5 + 10)

  command <<<
    set -euo pipefail

    bcftools --version | head -1 >&2
    ln -s "~{sites}" sites.vcf.gz
    ln -s "~{sites_index}" sites.vcf.gz.tbi
    ln -s "~{graph_calls}" graph_calls.tsv.gz
    ln -s "~{graph_calls_index}" graph_calls.tsv.gz.tbi
    mkdir members
    tar -xf "~{members_bundle}" -C members

    python3 - sites.vcf.gz graph_calls.tsv.gz "~{regions_bed}" members "~{write_lines(sample_ids)}" "~{write_lines(sample_sexes)}" "~{write_lines(read_sample_ids)}" "~{out_prefix}" <<'EOF' | bgzip > "~{out_prefix}.vcf.gz"
    # ---- graph_gt_join.py (begin) ----
    import gzip
    import os
    import re
    import subprocess
    import sys

    sites_vcf, graph_calls, regions_bed, members_dir, ids_file, sexes_file, read_ids_file, prefix = sys.argv[1:9]
    SPLIT = re.compile(r"[|/]")


    def lines_of(path):
        return [l.rstrip("\n") for l in open(path)]


    def size_class(n):
        return "50-999" if n < 1000 else ("1k-9.9k" if n < 10000 else "10k+")


    def run(cmd):
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            sys.exit("%s failed: %s" % (" ".join(cmd[:3]), res.stderr.strip()[-500:]))
        return res.stdout.splitlines()


    sample_ids = [s for s in lines_of(ids_file) if s]
    sexes = lines_of(sexes_file)
    sexes = (sexes + [""] * len(sample_ids))[:len(sample_ids)]
    read_ids = [s for s in lines_of(read_ids_file) if s]
    regions = []
    for line in open(regions_bed):
        f = line.split()
        if len(f) >= 3 and not line.startswith("#"):
            regions.append((f[0], int(f[1]), int(f[2])))


    def in_shard(chrom, pos):
        return any(c == chrom and start < pos <= end for c, start, end in regions)


    # the bundle: <index>__<file name>, the VCFs first, then their indexes, both in read_ids order
    entries = sorted(os.listdir(members_dir))
    if len(entries) != 2 * len(read_ids):
        sys.exit("the bundle has %d files for %d genotyped members" % (len(entries), len(read_ids)))
    member_vcf = {}
    for i, name in enumerate(entries):
        plain = name.split("__", 1)[1]
        os.rename(os.path.join(members_dir, name), os.path.join(members_dir, plain))
        if i < len(read_ids):
            member_vcf[read_ids[i]] = os.path.join(members_dir, plain)

    # ---- the shard's sites
    header = [l for l in run(["bcftools", "view", "-h", sites_vcf]) if l.startswith("##")]
    site_index, site_rows = {}, []
    if regions:
        for line in run(["bcftools", "view", "-H", "-R", regions_bed, sites_vcf]):
            f = line.split("\t")
            if f[2] in site_index or not in_shard(f[0], int(f[1])):
                continue
            info = dict(kv.split("=", 1) if "=" in kv else (kv, "") for kv in f[7].split(";"))
            site_index[f[2]] = len(site_rows)
            site_rows.append({"fixed": f[:7], "info": f[7], "svtype": info.get("SVTYPE", "."), "len": abs(int(info.get("SVLEN", "0"))),
                              "skip": info.get("SKIP", ""), "graph_ac": int(info.get("GRAPH_AC", "0")), "cols": [],
                              "context": "excluded_region" if "IN_EXCL" in info else ("tandem_repeat" if "IN_TR" in info else "other"),
                              "ac": 0, "an": 0, "ns": 0, "sex": {"XX": [0, 0], "XY": [0, 0]}, "calls": {}})

    # ---- the members' columns, one member's VCF at a time
    member_counts = {}   # (sample, svtype) -> [called, missing, hom_ref, het, hom_alt, dp_sum, dp_n]
    for sid, sex in zip(sample_ids, sexes):
        got = {}
        if sid in member_vcf and site_rows:
            for line in run(["bcftools", "query", "-R", regions_bed, "-f", "%ID[\t%GT\t%GQ\t%DP\t%AD]\n", member_vcf[sid]]):
                f = line.split("\t")
                if f[0] in site_index and len(f) >= 5:
                    got[f[0]] = f[1:5]
        for sid_site, i in site_index.items():
            row = site_rows[i]
            gt, gq, dp, ad = got.get(sid_site, (".", ".", ".", "."))
            parts = SPLIT.split(gt)
            alt = sum(1 for a in parts if a == "1")
            called = sum(1 for a in parts if a != ".")
            if len(parts) == 1:
                text = parts[0] if called or sid_site in got else "./."
            else:
                text = "/".join(sorted(parts, key=lambda a: (a != ".", a)))
            row["cols"].append(text if (gq, dp, ad) == (".", ".", ".") or not called else "%s:%s:%s:%s" % (text, gq, dp, ad))
            if row["skip"]:
                continue
            row["calls"][sid] = (alt, called)
            key = (sid, row["svtype"])
            c = member_counts.setdefault(key, [0, 0, 0, 0, 0, 0, 0])
            if not called:
                c[1] += 1
                continue
            c[0] += 1
            c[2 if alt == 0 else (4 if alt == called else 3)] += 1
            if dp.isdigit():
                c[5] += int(dp)
                c[6] += 1
            row["ac"] += alt
            row["an"] += called
            row["ns"] += 1
            if sex in row["sex"]:
                row["sex"][sex][0] += alt
                row["sex"][sex][1] += called

    # ---- reads against assemblies, for the members that are in the graph
    conc = {}   # (sample, svtype, size class, context) -> [compared, concordant, carrier_missed, carrier_added, zygosity, read_missing, not_comparable]
    with gzip.open(graph_calls, "rt") as fh:
        graph_samples = fh.readline().rstrip("\n").split("\t")[3:]
    both = [(j, s) for j, s in enumerate(graph_samples) if s in member_vcf and s in sample_ids]
    if both and site_rows:
        for line in run(["tabix", "-R", regions_bed, graph_calls]):
            f = line.split("\t")
            i = site_index.get(f[2])
            if i is None or site_rows[i]["skip"]:
                continue
            row = site_rows[i]
            if row.get("compared"):
                continue
            row["compared"] = True
            for j, s in both:
                g_alt, g_called = (int(x) for x in f[3 + j].split(":"))
                r_alt, r_called = row["calls"].get(s, (0, 0))
                c = conc.setdefault((s, row["svtype"], size_class(row["len"]), row["context"]), [0, 0, 0, 0, 0, 0, 0])
                if not g_called:
                    continue
                if not r_called:
                    c[5] += 1
                elif r_called != g_called:   # a half-missing assembly genotype, or a ploidy the two do not share
                    c[6] += 1
                else:
                    c[0] += 1
                    if r_alt == g_alt:
                        c[1] += 1
                    elif g_alt and not r_alt:
                        c[2] += 1
                    elif r_alt and not g_alt:
                        c[3] += 1
                    else:
                        c[4] += 1

    # ---- the shard's VCF
    out = sys.stdout
    for line in header:
        out.write(line + "\n")
    for tag, number, typ, text in (
            ("AC", "A", "Integer", "Alternate alleles in the cohort's called genotypes"),
            ("AN", "1", "Integer", "Called alleles in the cohort"),
            ("AF", "A", "Float", "AC / AN"),
            ("NS", "1", "Integer", "Members with a called genotype"),
            ("AC_XX", "A", "Integer", "AC among the XX members"), ("AN_XX", "1", "Integer", "AN among the XX members"),
            ("AF_XX", "A", "Float", "AC_XX / AN_XX"),
            ("AC_XY", "A", "Integer", "AC among the XY members"), ("AN_XY", "1", "Integer", "AN among the XY members"),
            ("AF_XY", "A", "Float", "AC_XY / AN_XY")):
        out.write('##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % (tag, number, typ, text))
    for tag, number, typ, text in (
            ("GT", "1", "String", "Genotype from the member's reads (kanpig), unphased"),
            ("GQ", "1", "Integer", "kanpig's genotype quality"),
            ("DP", "1", "Integer", "Reads over the site"),
            ("AD", "R", "Integer", "Reads supporting the reference and the alternate allele")):
        out.write('##FORMAT=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % (tag, number, typ, text))
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_ids) + "\n")


    def af(ac, an):
        return ("%.6g" % (ac / an)) if an else "."


    site_counts = {}   # (svtype, size class) -> [sites, genotyped, variable, graph_only, no_call, in a tandem repeat, in an excluded region]
    for row in site_rows:
        c = site_counts.setdefault((row["svtype"], size_class(row["len"])), [0, 0, 0, 0, 0, 0, 0])
        c[0] += 1
        c[5] += 1 if row["context"] == "tandem_repeat" else 0
        c[6] += 1 if row["context"] == "excluded_region" else 0
        info = row["info"]
        if not row["skip"]:
            c[1] += 1
            c[2] += 1 if row["ac"] else 0
            c[3] += 1 if row["an"] and not row["ac"] and row["graph_ac"] else 0
            c[4] += 0 if row["an"] else 1
            xx, xy = row["sex"]["XX"], row["sex"]["XY"]
            info += ";AC=%d;AN=%d;AF=%s;NS=%d;AC_XX=%d;AN_XX=%d;AF_XX=%s;AC_XY=%d;AN_XY=%d;AF_XY=%s" % (
                row["ac"], row["an"], af(row["ac"], row["an"]), row["ns"], xx[0], xx[1], af(*xx), xy[0], xy[1], af(*xy))
        out.write("\t".join(row["fixed"] + [info, "GT:GQ:DP:AD"] + row["cols"]) + "\n")

    with open(prefix + ".members.tsv", "w") as fh:
        fh.write("sample_id\tsvtype\tcalled\tmissing\thom_ref\thet\thom_alt\tdp_sum\tdp_n\n")
        for (sid, svtype), c in sorted(member_counts.items()):
            fh.write("\t".join([sid, svtype] + [str(x) for x in c]) + "\n")
    with open(prefix + ".concordance.tsv", "w") as fh:
        fh.write("sample_id\tsvtype\tsize_class\tcontext\tcompared\tconcordant\tcarrier_missed\tcarrier_added\tzygosity\tread_missing\tnot_comparable\n")
        for (sid, svtype, cls, context), c in sorted(conc.items()):
            fh.write("\t".join([sid, svtype, cls, context] + [str(x) for x in c]) + "\n")
    with open(prefix + ".site_counts.tsv", "w") as fh:
        fh.write("svtype\tsize_class\tsites\tgenotyped\tvariable\tgraph_only\tno_call\tin_tandem_repeat\tin_excluded_region\n")
        for (svtype, cls), c in sorted(site_counts.items()):
            fh.write("\t".join([svtype, cls] + [str(x) for x in c]) + "\n")
    print("%d sites, %d members with reads of %d, %d compared with their assembly" % (
        len(site_rows), len(member_vcf), len(sample_ids), len(both)), file=sys.stderr)
    # ---- graph_gt_join.py (end) ----
    EOF
  >>>

  output {
    File vcf = "~{out_prefix}.vcf.gz"
    File members = "~{out_prefix}.members.tsv"
    File concordance = "~{out_prefix}.concordance.tsv"
    File site_counts = "~{out_prefix}.site_counts.tsv"
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

task ugc_wgw_graph_gt_gather {
  meta {
    description: "The shards joined: the cohort VCF (concatenated without recompression) and its sites-only copy, the member table, the reads-against-assembly table with its totals, and the summary."
    outputs: {
      vcf: {
        description: "<out_prefix>.vcf.gz: every site with every member's genotype"
      },
      vcf_index: {
        description: "Its index"
      },
      sites_vcf: {
        description: "<out_prefix>.sites.vcf.gz: the same without the member columns (the frequency resource)"
      },
      sites_vcf_index: {
        description: "Its index"
      },
      members: {
        description: "<out_prefix>.members.tsv: per member its sex, whether it has reads and is in the graph, genotyped sites, call rate, genotype counts, mean depth"
      },
      concordance: {
        description: "<out_prefix>.concordance.tsv: per graph member, SV type, size class and context (excluded_region, tandem_repeat or other; the totals as ALL) the sites compared and how they agree"
      },
      summary: {
        description: "<out_prefix>.summary.tsv: resource, metric, key, value"
      }
    }
  }

  parameter_meta {
    vcfs: {
      description: "The shards' VCFs, in shard order"
    }
    member_counts: {
      description: "The shards' member counts"
    }
    concordance_counts: {
      description: "The shards' comparison counts"
    }
    site_counts: {
      description: "The shards' site counts"
    }
    sites_counts: {
      description: "The site step's counts"
    }
    sample_ids: {
      description: "The cohort's members, in order"
    }
    sample_sexes: {
      description: "XX, XY or empty per member"
    }
    read_sample_ids: {
      description: "The members that were genotyped"
    }
    graph_cohort_id: {
      description: "The cohort whose graph gave the sites"
    }
    parameters: {
      description: "The run's parameters, key=value separated by semicolons, written into the summary"
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
    Array[File] vcfs
    Array[File] member_counts
    Array[File] concordance_counts
    Array[File] site_counts
    File sites_counts
    Array[String] sample_ids
    Array[String] sample_sexes
    Array[String] read_sample_ids
    String graph_cohort_id
    String parameters
    String out_prefix
    Int mem_gb = 4
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(vcfs, "GB") * 3 + 10)

  command <<<
    set -euo pipefail

    bcftools --version | head -1 >&2
    bcftools concat --naive --file-list "~{write_lines(vcfs)}" --output "~{out_prefix}.vcf.gz"
    tabix --preset vcf "~{out_prefix}.vcf.gz"
    bcftools view --drop-genotypes --output-type z --output "~{out_prefix}.sites.vcf.gz" "~{out_prefix}.vcf.gz"
    tabix --preset vcf "~{out_prefix}.sites.vcf.gz"

    python3 - "~{write_lines(member_counts)}" "~{write_lines(concordance_counts)}" "~{write_lines(site_counts)}" "~{sites_counts}" "~{write_lines(sample_ids)}" "~{write_lines(sample_sexes)}" "~{write_lines(read_sample_ids)}" "~{graph_cohort_id}" "~{parameters}" "~{out_prefix}" <<'EOF'
    # ---- graph_gt_gather.py (begin) ----
    import sys

    members_list, conc_list, site_list, sites_counts, ids_file, sexes_file, read_ids_file, graph_cohort, parameters, prefix = sys.argv[1:11]


    def lines_of(path):
        return [l.rstrip("\n") for l in open(path)]


    def rows(list_file):
        for path in lines_of(list_file):
            if not path:
                continue
            with open(path) as fh:
                cols = fh.readline().rstrip("\n").split("\t")
                for line in fh:
                    yield dict(zip(cols, line.rstrip("\n").split("\t")))


    def share(a, b):
        return ("%.4f" % (a / b)) if b else ""


    sample_ids = [s for s in lines_of(ids_file) if s]
    sexes = (lines_of(sexes_file) + [""] * len(sample_ids))[:len(sample_ids)]
    read_ids = set(s for s in lines_of(read_ids_file) if s)
    graph_samples, site_step = [], []
    with open(sites_counts) as fh:
        fh.readline()
        for line in fh:
            metric, key, value = line.rstrip("\n").split("\t")
            if metric == "graph_sample":
                graph_samples.append(key)
            else:
                site_step.append((metric, key, value))

    FIELDS = ("called", "missing", "hom_ref", "het", "hom_alt", "dp_sum", "dp_n")
    per_member = {}   # sample -> svtype -> totals
    for r in rows(members_list):
        t = per_member.setdefault(r["sample_id"], {}).setdefault(r["svtype"], dict.fromkeys(FIELDS, 0))
        for k in FIELDS:
            t[k] += int(r[k])

    tot_called = tot_sites = 0
    with open(prefix + ".members.tsv", "w") as out:
        out.write("sample_id\tsex\treads\tin_graph\tsites\tcalled\tmissing\tcall_rate\thom_ref\thet\thom_alt\tins_carried\tdel_carried\tmean_dp\n")
        for sid, sex in zip(sample_ids, sexes):
            by_type = per_member.get(sid, {})
            t = {k: sum(v[k] for v in by_type.values()) for k in FIELDS}
            n = t["called"] + t["missing"]
            carried = {st: by_type.get(st, {}).get("het", 0) + by_type.get(st, {}).get("hom_alt", 0) for st in ("INS", "DEL")}
            if sid in read_ids:
                tot_called += t["called"]
                tot_sites += n
            out.write("\t".join([sid, sex, "yes" if sid in read_ids else "no", "yes" if sid in graph_samples else "no", str(n),
                                 str(t["called"]), str(t["missing"]), share(t["called"], n) if sid in read_ids else "",
                                 str(t["hom_ref"]), str(t["het"]), str(t["hom_alt"]), str(carried["INS"]), str(carried["DEL"]),
                                 ("%.1f" % (t["dp_sum"] / t["dp_n"])) if t["dp_n"] else ""]) + "\n")

    CFIELDS = ("compared", "concordant", "carrier_missed", "carrier_added", "zygosity", "read_missing", "not_comparable")
    conc = {}
    for r in rows(conc_list):
        s, t_, c_, x = r["sample_id"], r["svtype"], r["size_class"], r["context"]
        for key in ((s, t_, c_, x), (s, "ALL", "ALL", x), (s, "ALL", "ALL", "ALL"), ("ALL", t_, c_, "ALL"), ("ALL", t_, "ALL", "ALL"),
                    ("ALL", "ALL", c_, "ALL"), ("ALL", "ALL", "ALL", x), ("ALL", "ALL", "ALL", "ALL")):
            t = conc.setdefault(key, dict.fromkeys(CFIELDS, 0))
            for k in CFIELDS:
                t[k] += int(r[k])
    with open(prefix + ".concordance.tsv", "w") as out:
        out.write("sample_id\tsvtype\tsize_class\tcontext\t" + "\t".join(CFIELDS) + "\tconcordance\n")
        for key in sorted(conc, key=lambda k: (k[0] != "ALL", k[0], k[1] != "ALL", k[1], k[2] != "ALL", k[2], k[3] != "ALL", k[3])):
            t = conc[key]
            out.write("\t".join(list(key) + [str(t[k]) for k in CFIELDS] + [share(t["concordant"], t["compared"])]) + "\n")

    SFIELDS = ("sites", "genotyped", "variable", "graph_only", "no_call", "in_tandem_repeat", "in_excluded_region")
    by_class = {}
    for r in rows(site_list):
        for key in ((r["svtype"], r["size_class"]), (r["svtype"], "ALL"), ("ALL", "ALL")):
            t = by_class.setdefault(key, dict.fromkeys(SFIELDS, 0))
            for k in SFIELDS:
                t[k] += int(r[k])

    total = by_class.get(("ALL", "ALL"), dict.fromkeys(SFIELDS, 0))
    listed = sum(int(v) for m, _, v in site_step if m == "sites")
    all_conc = conc.get(("ALL", "ALL", "ALL", "ALL"), dict.fromkeys(CFIELDS, 0))
    in_both = [s for s in sample_ids if s in graph_samples and s in read_ids]
    with open(prefix + ".summary.tsv", "w") as out:
        out.write("resource\tmetric\tkey\tvalue\n")

        def w(metric, key, value):
            out.write("graph_gt\t%s\t%s\t%s\n" % (metric, key, value))

        w("graph", "cohort", graph_cohort)
        w("graph", "members", len(graph_samples))
        w("members", "cohort", len(sample_ids))
        w("members", "with_reads", sum(1 for s in sample_ids if s in read_ids))
        w("members", "without_reads", sum(1 for s in sample_ids if s not in read_ids))
        w("members", "in_graph_with_reads", len(in_both))
        for metric, key, value in site_step:
            w("site_list_" + metric, key, value)
        w("sites", "listed", listed)
        w("sites", "in_region_shards", total["sites"])
        w("sites", "outside_region_shards", listed - total["sites"])
        for (svtype, cls), t in sorted(by_class.items(), key=lambda kv: (kv[0][0] != "ALL", kv[0][0], kv[0][1] != "ALL", kv[0][1])):
            for k in SFIELDS:
                w("sites_" + k, svtype + ":" + cls, t[k])
        w("genotypes", "called", tot_called)
        w("genotypes", "possible", tot_sites)
        w("genotypes", "call_rate", share(tot_called, tot_sites))
        for k in CFIELDS:
            w("concordance", k, all_conc[k])
        w("concordance", "share", share(all_conc["concordant"], all_conc["compared"]))
        for context in ("other", "tandem_repeat", "excluded_region"):
            t = conc.get(("ALL", "ALL", "ALL", context))
            if t:
                w("concordance_" + context, "compared", t["compared"])
                w("concordance_" + context, "share", share(t["concordant"], t["compared"]))
        for kv in parameters.split(";"):
            if "=" in kv:
                key, value = kv.split("=", 1)
                w("parameter", key, value)
    # ---- graph_gt_gather.py (end) ----
    EOF
    cat "~{out_prefix}.summary.tsv" >&2
  >>>

  output {
    File vcf = "~{out_prefix}.vcf.gz"
    File vcf_index = "~{out_prefix}.vcf.gz.tbi"
    File sites_vcf = "~{out_prefix}.sites.vcf.gz"
    File sites_vcf_index = "~{out_prefix}.sites.vcf.gz.tbi"
    File members = "~{out_prefix}.members.tsv"
    File concordance = "~{out_prefix}.concordance.tsv"
    File summary = "~{out_prefix}.summary.tsv"
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
