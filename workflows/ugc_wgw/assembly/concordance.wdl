version 1.0
# ugc-pacbio-wgw — the assembly's SNV genotypes against the sample's read-based calls inside the reliable region: genotype concordance, het loss, assembly-only hets
# origin: none
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8.1
import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_asm_concordance {
  meta {
    description: "PASS SNVs of the assembly VCF and of the read-based VCF within the region BED, merged site by site (bcftools merge --force-samples) and compared as unphased genotypes: sites called by both and their concordance; het loss, the read-based hets the assembly calls homozygous or not at all (collapsed haplotypes, the dominant failure without trio or Hi-C data); and the assembly's hets the reads do not call. One TSV of metric, value."
    outputs: {
      concordance_tsv: {
        description: "<sample_id>.asm.concordance.tsv"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID; output prefix"
    }
    asm_vcf: {
      description: "The assembly caller's finished VCF (dipcall)"
    }
    asm_vcf_index: {
      description: "Its index"
    }
    read_vcf: {
      description: "The sample's read-based small-variant VCF"
    }
    read_vcf_index: {
      description: "Its index"
    }
    region_bed: {
      description: "Where to compare: reliable.bed, or dip.bed when there is no Flagger mask"
    }
    region_name: {
      description: "What the region is, recorded in the TSV"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    File asm_vcf
    File asm_vcf_index
    File read_vcf
    File read_vcf_index
    File region_bed
    String region_name
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil((size(asm_vcf, "GB") + size(read_vcf, "GB")) * 4 + 10)

  command <<<
    set -euo pipefail
    bcftools --version | head -1 >&2
    ln -s "~{asm_vcf}" asm.vcf.gz; ln -s "~{asm_vcf_index}" "asm.vcf.gz.$(basename "~{asm_vcf_index}" | sed 's/.*\.//')"
    ln -s "~{read_vcf}" reads.vcf.gz; ln -s "~{read_vcf_index}" "reads.vcf.gz.$(basename "~{read_vcf_index}" | sed 's/.*\.//')"
    sort -k1,1 -k2,2n "~{region_bed}" > region.bed
    for v in asm reads; do
      bcftools view --apply-filters PASS,. --types snps --regions-file region.bed -Oz -o "$v.pass.vcf.gz" "$v.vcf.gz"
      bcftools index "$v.pass.vcf.gz"
    done
    bcftools merge --force-samples -Ou asm.pass.vcf.gz reads.pass.vcf.gz \
      | bcftools query -f '%CHROM\t%POS[\t%GT]\n' \
      > merged.tsv
    python3 - merged.tsv "~{sample_id}.asm.concordance.tsv" "~{region_name}" <<'PY'
import sys
def alleles(gt):
    if gt in (".", "./.", ".|."):
        return None
    parts = gt.replace("|", "/").split("/")
    if any(p == "." for p in parts):
        return None
    return tuple(sorted(parts))
n = {"sites_merged": 0, "both_called": 0, "concordant": 0, "discordant": 0, "read_het": 0, "het_loss": 0,
     "read_only": 0, "asm_only": 0, "asm_only_het": 0}
with open(sys.argv[1]) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        a, r = alleles(f[2]), alleles(f[3])
        n["sites_merged"] += 1
        if r is not None and len(set(r)) == 2:
            n["read_het"] += 1
            if a is None or len(set(a)) == 1:
                n["het_loss"] += 1
        if a is not None and r is not None:
            n["both_called"] += 1
            n["concordant" if a == r else "discordant"] += 1
        elif a is None:
            n["read_only"] += 1
        else:
            n["asm_only"] += 1
            if len(set(a)) == 2:
                n["asm_only_het"] += 1
with open(sys.argv[2], "w") as out:
    out.write(f"region\t{sys.argv[3]}\n")
    for k, v in n.items():
        out.write(f"{k}\t{v}\n")
    out.write("genotype_concordance\t%s\n" % (f"{n['concordant'] / n['both_called']:.4f}" if n["both_called"] else "NA"))
    out.write("het_loss_rate\t%s\n" % (f"{n['het_loss'] / n['read_het']:.4f}" if n["read_het"] else "NA"))
PY
    cat "~{sample_id}.asm.concordance.tsv" >&2
  >>>

  output {
    File concordance_tsv = "~{sample_id}.asm.concordance.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
    cpu: 2
    memory: "8 GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}
