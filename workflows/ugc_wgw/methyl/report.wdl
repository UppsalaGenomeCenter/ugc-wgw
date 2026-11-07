version 1.0

# ugc-pacbio-wgw — methbat report: a sample's imprinting regions labelled against their expected methylation state
# origin: none (the task shape follows upstream's methbat_profile, vendor/…/tasks/methbat.wdl @ v4.0.0; MethBat 1.1.0)
# see docs/DESIGN.md §5, §13, §18 (2026-11-06); docs/ENTRYPOINTS.md §3, §5

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_methbat_report {
  meta {
    description: "methbat report on the sample's 5mC pileup: every region of the report regions file (MethBat's imprinting targets, expected AlleleSpecificMethylation) labelled PASS, Inconclusive, Anomalous or AnomalousQcWarning against its expected category, with the QC warnings (LowPhasedSites, LowHaplotypeCoverage, WeakASM, ...) that qualify the call. Per sample; cohort_methyl tabulates the reports across members and makes one here for a member whose run predates the task."
    outputs: {
      report: {
        description: "<out_prefix>.methbat.report.tsv: chrom, start, end, region_label, report_summary, qc_warnings, expected_category, summary_label, the region's methylation means and site counts"
      }
    }
  }

  parameter_meta {
    cpg_pileup_bed: {
      description: "The sample's 5mCpG pileup BED (methbat pileup, bgzipped)"
    }
    report_regions_tsv: {
      description: "MethBat report regions: chrom, start, end, cpg_label, expected_category, anomalous_categories (the reference map's methbat_report_regions_tsv)"
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
    File cpg_pileup_bed
    File report_regions_tsv
    String out_prefix
    Int mem_gb = 4
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(cpg_pileup_bed, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    ln --symbolic --verbose "~{cpg_pileup_bed}" .

    methbat --version

    methbat report \
      --input-pileup "~{basename(cpg_pileup_bed)}" \
      --input-regions "~{report_regions_tsv}" \
      --output-report "~{out_prefix}.methbat.report.tsv"
  >>>

  output {
    File report = "~{out_prefix}.methbat.report.tsv"
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
