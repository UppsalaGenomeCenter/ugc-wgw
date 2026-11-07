version 1.0

# ugc-pacbio-wgw — checks of externally produced inputs (an aligned BAM, pre-called VCFs) against the installed reference, and the VCF sample column renamed to the sample id
# origin: none
# see docs/DESIGN.md §5, §8.1, §18 (2026-11-06); docs/ENTRYPOINTS.md §3, §11

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_input_checks {
  meta {
    description: "The input checks of the bam and vcf entry kinds, before any job is spent on them: the BAM is coordinate-sorted, indexed, aligned to the installed reference (every @SQ contig in the FASTA index with the same length; contigs the BAM lacks are noted) and carries MM/ML methylation tags in its first reads (a failure unless require_methylation_tags is off, then a note); every VCF given is BGZF-compressed, indexed, has exactly one sample and only contigs of the reference; the BAM's read-group sample name is the sample id (the callers name their outputs after it). Writes the checks table and the provenance table (the BAM header's @PG and @RG lines, the VCF headers' source lines) first, then exits non-zero with INPUT CHECK FAILED lines on stderr when a check fails: the driver classifies that as an input error (blocks the sample, no retry)."
    outputs: {
      checks: {
        description: "<sample_id>.input_checks.tsv: check, target, status (pass, warn, fail, skip), detail"
      },
      provenance: {
        description: "<sample_id>.input_provenance.tsv: source (bam, small_variant_vcf, sv_vcf, gvcf), key, value: the BAM's @PG (ID, PN, VN, CL) and @RG (ID, SM, PU) fields, the VCFs' ##fileformat, ##source and tool version lines"
      },
      vcf_sample: {
        description: "The sample column of small_variant_vcf (empty without one)"
      },
      msg: {
        description: "One line per warning"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample"
    }
    aligned_bam: {
      description: "The aligned BAM given in the sample sheet"
    }
    aligned_bam_index: {
      description: "Its index"
    }
    small_variant_vcf: {
      description: "The small-variant VCF given (vcf kind), if any"
    }
    small_variant_vcf_index: {
      description: "Its index (checked present and non-empty)"
    }
    sv_vcf: {
      description: "The SV VCF given, if any"
    }
    sv_vcf_index: {
      description: "Its index"
    }
    gvcf: {
      description: "The gVCF given, if any"
    }
    gvcf_index: {
      description: "Its index"
    }
    ref_index: {
      description: "The reference FASTA index (names and lengths)"
    }
    require_methylation_tags: {
      description: "Fail when the first reads carry no MM/ML tags (true); false notes it and the methylation steps downstream are skipped by the driver"
    }
    sample_reads: {
      description: "How many reads to scan for MM/ML tags"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    File aligned_bam
    File aligned_bam_index
    File? small_variant_vcf
    File? small_variant_vcf_index
    File? sv_vcf
    File? sv_vcf_index
    File? gvcf
    File? gvcf_index
    File ref_index
    Boolean require_methylation_tags = true
    Int sample_reads = 2000
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = 10

  command <<<
    set -euo pipefail

    samtools --version | head -1 >&2
    bcftools --version | head -1 >&2

    python3 - "~{sample_id}" "~{aligned_bam}" "~{aligned_bam_index}" "~{ref_index}" "~{require_methylation_tags}" ~{sample_reads} \
      "~{default="" small_variant_vcf}" "~{default="" sv_vcf}" "~{default="" gvcf}" \
      "~{default="" small_variant_vcf_index}" "~{default="" sv_vcf_index}" "~{default="" gvcf_index}" <<'EOF'
    # ---- input_checks.py (begin) ----
    import os, subprocess, sys
    sample_id, bam, bai, fai, require_mm, n_reads, small_vcf, sv_vcf, gvcf, small_idx, sv_idx, gvcf_idx = sys.argv[1:13]
    require_mm, n_reads = require_mm == "true", int(n_reads)
    checks, prov, messages = [], [], []

    def check(name, target, status, detail=""):
        checks.append((name, target, status, detail))
        if status == "warn":
            messages.append(f"{name} ({target}): {detail}")

    def run(cmd):
        return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout

    ref = {}
    with open(fai) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 2:
                ref[c[0]] = int(c[1])

    # ---- the BAM
    check("bam_index", "bam", "pass" if os.path.getsize(bai) > 0 else "fail", os.path.basename(bai) if os.path.getsize(bai) > 0 else "the index is empty")
    header = run(["samtools", "view", "-H", bam]).splitlines()
    sq = {}
    sort_order = ""
    for line in header:
        f = line.split("\t")
        if f[0] == "@SQ":
            d = dict(x.split(":", 1) for x in f[1:] if ":" in x)
            sq[d.get("SN", "")] = int(d.get("LN", "0"))
        elif f[0] == "@HD":
            d = dict(x.split(":", 1) for x in f[1:] if ":" in x)
            sort_order = d.get("SO", "")
        elif f[0] == "@PG":
            d = dict(x.split(":", 1) for x in f[1:] if ":" in x)
            for k in ("ID", "PN", "VN", "CL"):
                if k in d:
                    prov.append(("bam", f"PG.{d.get('ID', '')}.{k}", d[k][:500]))
        elif f[0] == "@RG":
            d = dict(x.split(":", 1) for x in f[1:] if ":" in x)
            for k in ("ID", "SM", "PU", "PL"):
                if k in d:
                    prov.append(("bam", f"RG.{d.get('ID', '')}.{k}", d[k]))
    check("bam_sorted", "bam", "pass" if sort_order == "coordinate" else "fail", f"@HD SO:{sort_order or 'absent'}")
    if not sq:
        check("bam_reference", "bam", "fail", "no @SQ lines: an unaligned BAM (register it under hifi_reads)")
    else:
        bad = [f"{n} (BAM {l}, reference {ref.get(n, 'absent')})" for n, l in sq.items() if ref.get(n) != l]
        if bad:
            check("bam_reference", "bam", "fail", f"{len(bad)} of {len(sq)} contigs differ from the reference: " + "; ".join(bad[:5]))
        else:
            missing = [n for n in ref if n not in sq]
            check("bam_reference", "bam", "pass" if not missing else "warn",
                  f"{len(sq)} contigs match the reference" + (f"; {len(missing)} reference contig(s) absent from the BAM" if missing else ""))
    # the callers name their outputs after the read group (sawfish its sample directory and VCF column, TRGT its VCF column;
    # HiPhase is told the sample id): the read-group sample name must be the sample id
    samples = sorted({v for s, k, v in prov if s == "bam" and k.endswith(".SM")})
    if samples == [sample_id]:
        check("bam_sample", "bam", "pass", f"@RG SM {sample_id}")
    elif not samples:
        check("bam_sample", "bam", "fail", "no @RG SM in the header; the callers name their outputs after it: add a read group (samtools addreplacerg)")
    else:
        check("bam_sample", "bam", "fail", f"@RG SM {', '.join(samples)} is not the sample id {sample_id}: register the sample under that name, "
              "or rewrite the header (samtools reheader) before registering")
    # the first reads only: stream them and stop (a whole-genome BAM is tens of GB of SAM text)
    proc = subprocess.Popen(["samtools", "view", bam], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    lines = []
    assert proc.stdout is not None
    for line in proc.stdout:
        lines.append(line)
        if len(lines) >= n_reads:
            break
    proc.kill()
    proc.wait()
    with_mm = sum(1 for l in lines if "\tMM:Z:" in l or "\tMm:Z:" in l)
    with_ml = sum(1 for l in lines if "\tML:B:" in l or "\tMl:B:" in l)
    if not lines:
        check("bam_methylation_tags", "bam", "fail", "no reads")
    elif with_mm and with_ml:
        check("bam_methylation_tags", "bam", "pass", f"{with_mm} of {len(lines)} first reads carry MM/ML")
    else:
        check("bam_methylation_tags", "bam", "fail" if require_mm else "warn",
              f"none of the first {len(lines)} reads carries MM/ML (5mC) tags" + ("" if require_mm else "; the methylation steps are skipped"))

    # ---- the VCFs
    for role, path, idx in (("small_variant_vcf", small_vcf, small_idx), ("sv_vcf", sv_vcf, sv_idx), ("gvcf", gvcf, gvcf_idx)):
        if not path:
            check("vcf", role, "skip", "not given")
            continue
        kind = run(["htsfile", path]).strip()
        if "BGZF-compressed" not in kind:
            check("vcf_bgzf", role, "fail", kind.split("\t")[-1] if "\t" in kind else kind)
            continue
        check("vcf_bgzf", role, "pass", kind.split("\t")[-1])
        # the index is a task input of its own (miniwdl mounts every input under its own directory, so "next to the VCF" means nothing here)
        if idx and os.path.getsize(idx) > 0:
            check("vcf_index", role, "pass", os.path.basename(idx))
        else:
            check("vcf_index", role, "fail", "no index given" if not idx else "the index is empty")
        vcf_samples = run(["bcftools", "query", "-l", path]).split()
        check("vcf_samples", role, "pass" if len(vcf_samples) == 1 else "fail",
              f"sample column {vcf_samples[0]}" if len(vcf_samples) == 1 else f"{len(vcf_samples)} sample columns; exactly one is required")
        contigs = {}
        for line in run(["bcftools", "view", "-h", path]).splitlines():
            if line.startswith("##contig=<"):
                body = line[len("##contig=<"):].rstrip(">")
                d = dict(x.split("=", 1) for x in body.split(",") if "=" in x)
                contigs[d.get("ID", "")] = int(d["length"]) if d.get("length", "").isdigit() else None
            elif line.startswith(("##fileformat", "##source", "##DeepVariant", "##sawfish", "##GLnexus", "##bcftools_", "##reference")):
                k, _, v = line[2:].partition("=")
                prov.append((role, k, v[:500]))
        unknown = [n for n in contigs if n not in ref]
        wrong = [n for n, l in contigs.items() if n in ref and l is not None and l != ref[n]]
        if unknown or wrong:
            check("vcf_contigs", role, "fail", (f"{len(unknown)} contig(s) not in the reference: {', '.join(unknown[:5])}" if unknown else "")
                  + ("; " if unknown and wrong else "") + (f"{len(wrong)} contig(s) with another length: {', '.join(wrong[:5])}" if wrong else ""))
        else:
            check("vcf_contigs", role, "pass", f"{len(contigs)} header contigs, all in the reference")
        if role == "small_variant_vcf":
            with open("vcf_sample.txt", "w") as fh:
                fh.write((vcf_samples[0] if len(vcf_samples) == 1 else "") + "\n")
    if not os.path.exists("vcf_sample.txt"):
        with open("vcf_sample.txt", "w") as fh:
            fh.write("\n")

    with open(sample_id + ".input_checks.tsv", "w") as fh:
        fh.write("check\ttarget\tstatus\tdetail\n")
        for c in checks:
            fh.write("\t".join(c) + "\n")
    with open(sample_id + ".input_provenance.tsv", "w") as fh:
        fh.write("source\tkey\tvalue\n")
        for p in prov:
            fh.write("\t".join(p) + "\n")
    with open("messages.txt", "w") as fh:
        for m in messages:
            fh.write(f"{sample_id}: input check {m}\n")
    failed = [c for c in checks if c[2] == "fail"]
    for c in failed:
        print(f"INPUT CHECK FAILED: {c[0]} ({c[1]}): {c[3]}", file=sys.stderr)
    print(f"{sample_id}: {len(checks)} checks, {len(failed)} failed, {len(messages)} warning(s)", file=sys.stderr)
    sys.exit(1 if failed else 0)
    # ---- input_checks.py (end) ----
    EOF
  >>>

  output {
    File checks = "~{sample_id}.input_checks.tsv"
    File provenance = "~{sample_id}.input_provenance.tsv"
    String vcf_sample = read_string("vcf_sample.txt")
    Array[String] msg = read_lines("messages.txt")
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

task ugc_wgw_vcf_sample_rename {
  meta {
    description: "A given VCF with its single sample column named after the sample id (bcftools reheader when it differs, a copy otherwise), bgzipped and indexed (.tbi), under the name the per-sample stages give that file."
    outputs: {
      vcf: {
        description: "<out_name>"
      },
      vcf_index: {
        description: "<out_name>.tbi"
      }
    }
  }

  parameter_meta {
    vcf: {
      description: "The given VCF (BGZF)"
    }
    vcf_index: {
      description: "Its index"
    }
    sample_id: {
      description: "The sample id the column must carry"
    }
    out_name: {
      description: "Output file name (<sample_id>.<ref_name>.<suffix>.vcf.gz)"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File vcf
    File vcf_index
    String sample_id
    String out_name
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(vcf, "GB") * 3 + 5)

  command <<<
    set -euo pipefail

    bcftools --version | head -1 >&2
    echo "~{vcf_index}" >&2

    current=$(bcftools query -l "~{vcf}" | head -1)
    if [ "$current" = "~{sample_id}" ]; then
      cp "~{vcf}" "~{out_name}"
    else
      echo "renaming the sample column $current to ~{sample_id}" >&2
      echo "~{sample_id}" > names.txt
      bcftools reheader --samples names.txt --output "~{out_name}" "~{vcf}"
    fi
    bcftools index --tbi --force "~{out_name}"
  >>>

  output {
    File vcf_out = "~{out_name}"
    File vcf_out_index = "~{out_name}.tbi"
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
