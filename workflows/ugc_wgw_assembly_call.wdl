version 1.0
# ugc-pacbio-wgw — stage `assembly_call`: one sample's haplotype assemblies aligned to the assembly reference, called with dipcall and PAV 3, masked with Flagger's regions, compared with the read-based calls
# origin: none
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8.1
import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "ugc_wgw/assembly/align_asm.wdl" as Align
import "ugc_wgw/assembly/chrom_coverage.wdl" as ChromCoverage
import "ugc_wgw/assembly/concordance.wdl" as Concordance
import "ugc_wgw/assembly/dipcall.wdl" as Dipcall
import "ugc_wgw/assembly/masking.wdl" as Masking
import "ugc_wgw/assembly/pav.wdl" as Pav
import "ugc_wgw/assembly/vcf_finish.wdl" as VcfFinish
import "ugc_wgw/assembly/call_spectra.wdl" as CallSpectra
import "ugc_wgw/provenance.wdl" as Provenance

workflow ugc_wgw_assembly_call {
  meta {
    description: "ugc-pacbio-wgw stage assembly_call: minimap2 2.31 aligns each haplotype of the assembly stage's contract to the decoy-free GRCh38 (asm_fasta); from the PAFs the stage reports what was assembled per chromosome and infers the karyotype (XX, XY or unknown); dipcall v0.3 (male mode with the PAR BED when XY; chrX/chrY contigs regrouped first unless the haplotypes are trio-binned) and PAV 3 call variants; the primary VCFs carry the sample's name and unphased genotypes unless trio. With flagger_bed, Flagger's non-Hap regions are lifted to the reference through each haplotype's PAF and reliable.bed = dip.bed minus them; records outside it get FILTER ASM_UNRELIABLE in both primary VCFs (annotated, never dropped). With read_small_variant_vcf, dipcall's PASS SNVs are compared with the read-based ones inside reliable.bed (dip.bed without Flagger): genotype concordance, het loss, assembly-only hets."
    outputs: {
      haplotype_labels: {
        description: "The assembly's haplotype labels, the order of every per-haplotype array"
      },
      asm_bams: {
        description: "Contigs aligned to the assembly reference, per haplotype"
      },
      asm_bam_indices: {
        description: "Indices for asm_bams"
      },
      asm_pafs: {
        description: "The same alignments as gzipped PAF with cs tags, per haplotype"
      },
      chrom_coverage: {
        description: "Aligned bases and contigs per haplotype and chromosome"
      },
      karyotype_json: {
        description: "The inferred karyotype with its evidence and thresholds"
      },
      inferred_karyotype: {
        description: "XX, XY or unknown"
      },
      regroup_log: {
        description: "Contigs moved between the dipcall haplotype copies (non-trio XY samples)"
      },
      dipcall_vcf: {
        description: "dipcall's diploid calls: <sample_id>.dipcall.<ref_name>.vcf.gz, sample renamed, unphased unless trio"
      },
      dipcall_vcf_index: {
        description: "Tabix index"
      },
      dipcall_vcf_stats: {
        description: "bcftools stats of dipcall_vcf"
      },
      dipcall_raw_vcf: {
        description: "dipcall's VCF as written (sample syndip, phased)"
      },
      dipcall_bed: {
        description: "dipcall's confident region: covered by exactly one contig of each haplotype"
      },
      dipcall_hap_beds: {
        description: "dipcall's per-haplotype covered regions"
      },
      pav_vcf: {
        description: "PAV's calls: <sample_id>.pav.<ref_name>.vcf.gz, unphased unless trio"
      },
      pav_vcf_index: {
        description: "Tabix index"
      },
      pav_vcf_stats: {
        description: "bcftools stats of pav_vcf"
      },
      call_spectra: {
        description: "<sample_id>.asm.call_spectra.tsv: both callers' records by FILTER class (PASS, ASM_UNRELIABLE, OTHER, BOTH), type, size bin and contig class"
      },
      pav_raw_vcf: {
        description: "PAV's VCF as written (phased in haplotype order)"
      },
      pav_raw_vcf_index: {
        description: "Its .csi"
      },
      pav_callable_beds: {
        description: "Reference regions PAV could call, per haplotype"
      },
      unreliable_beds: {
        description: "Flagger's non-Hap regions lifted to the reference, per haplotype (with flagger_bed)"
      },
      reliable_bed: {
        description: "dip.bed minus both haplotypes' unreliable regions (with flagger_bed)"
      },
      reliable_summary: {
        description: "Bases of dip.bed, unreliable, reliable, and the reliable fraction of the reference (with flagger_bed)"
      },
      concordance: {
        description: "Genotype concordance, het loss and assembly-only hets of dipcall's SNVs against the read-based calls (with read_small_variant_vcf)"
      },
      ugc_wgw_workflow_version: {
        description: "ugc_wgw_version echoed"
      },
      ugc_wgw_manifest: {
        description: "Workflow-side provenance record"
      }
    }
  }

  parameter_meta {
    ugc_wgw_version: {
      description: "ugc-pacbio-wgw version the run was submitted with"
    }
    sample_id: {
      description: "Sample ID; output prefix"
    }
    hap_fastas: {
      description: "The assembly stage's hap_fastas: two bgzipped haplotype FASTAs"
    }
    haplotype_labels: {
      description: "The assembly stage's haplotype_labels (hap1, hap2; pat, mat when trio-binned)"
    }
    trio: {
      description: "The assembly stage's trio: hap1 is paternal, genotypes keep their phase and dipcall takes the haplotypes as they are"
    }
    read_small_variant_vcf: {
      description: "The sample's read-based small-variant VCF (singleton, or downstream in joint mode), when one exists: the concordance check runs"
    }
    read_small_variant_vcf_index: {
      description: "Its index"
    }
    flagger_bed: {
      description: "HMM-Flagger's prediction BED from the assembly stage: the masking runs"
    }
    align_threads: {
      description: "minimap2 task threads per haplotype"
    }
    align_mem_gb: {
      description: "minimap2 task memory (GB) per haplotype"
    }
    dipcall_threads: {
      description: "dipcall's -t (the task takes twice this: both haplotypes align at once)"
    }
    dipcall_mem_gb: {
      description: "dipcall task memory (GB)"
    }
    pav_threads: {
      description: "PAV's --cores"
    }
    pav_align_threads: {
      description: "Threads of each PAV align_map job"
    }
    pav_mem_gb: {
      description: "PAV task memory (GB)"
    }
    chrx_min_fraction: {
      description: "Below this fraction of chrX assembled (non-PAR, both haplotypes) the karyotype is unknown"
    }
    chry_min_fraction: {
      description: "At or above this fraction of chrY assembled (non-PAR) the karyotype is XY"
    }
    ref_map_file: {
      description: "ugc-wgw reference map (docs/DESIGN.md §10.4); keys used here: name, asm_fasta, asm_fasta_index, asm_par_bed"
    }
    ugc_wgw_container_registry: {
      description: "Registry holding the ugc-built images (minimap2)"
    }
    backend: {
      description: "Backend where the workflow will be executed",
      choices: [
        "GCP",
        "Azure",
        "AWS-HealthOmics",
        "HPC"
      ]
    }
    zones: {
      description: "Zones where compute will take place; required if backend is set to 'GCP'"
    }
    cpuPlatform: {
      description: "Optional minimum CPU platform to use for tasks on GCP"
    }
    gpuType: {
      description: "Optional type of GPU/Accelerator to use (no GPU task in this stage)"
    }
    container_registry: {
      description: "Container registry for upstream images; default quay.io/pacbio"
    }
    preemptible: {
      description: "Where possible, run tasks preemptibly (no effect on HPC)"
    }
  }

  input {
    String ugc_wgw_version
    String sample_id
    Array[File] hap_fastas
    Array[String] haplotype_labels
    Boolean trio
    File? read_small_variant_vcf
    File? read_small_variant_vcf_index
    File? flagger_bed
    Int align_threads = 16
    Int align_mem_gb = 32
    Int dipcall_threads = 8
    Int dipcall_mem_gb = 48
    Int pav_threads = 8
    Int pav_align_threads = 4
    Int pav_mem_gb = 40
    Float chrx_min_fraction = 0.5
    Float chry_min_fraction = 0.1
    File ref_map_file
    String ugc_wgw_container_registry

    # Backend configuration
    String backend = "HPC"
    String? zones
    String? cpuPlatform
    String? gpuType
    String? container_registry
    Boolean preemptible = true
  }

  call BackendConfiguration.backend_configuration { input:
    backend = backend,
    zones = zones,
    cpuPlatform = cpuPlatform,
    gpuType = gpuType,
    container_registry = container_registry
  }

  RuntimeAttributes default_runtime_attributes = if preemptible
    then backend_configuration.spot_runtime_attributes
    else backend_configuration.on_demand_runtime_attributes

  #@ except: DeclarationName
  Map[String, String] ref_map = read_map(ref_map_file)
  String ref_name = ref_map["name"]

  scatter (i in range(length(hap_fastas))) {
    call Align.ugc_wgw_minimap2_align_asm { input:
      sample_id = sample_id,
      haplotype = haplotype_labels[i],
      query_fasta = hap_fastas[i],
      ref_fasta = ref_map["asm_fasta"],  # !FileCoercion
      ref_name = ref_name,
      threads = align_threads,
      mem_gb = align_mem_gb,
      ugc_wgw_container_registry = ugc_wgw_container_registry,
      runtime_attributes = default_runtime_attributes
    }
  }

  call ChromCoverage.ugc_wgw_chrom_coverage { input:
    sample_id = sample_id,
    pafs = ugc_wgw_minimap2_align_asm.paf,
    hap_fastas = hap_fastas,
    haplotype_labels = haplotype_labels,
    ref_index = ref_map["asm_fasta_index"],  # !FileCoercion
    par_bed = ref_map["asm_par_bed"],  # !FileCoercion
    trio = trio,
    chrx_min_fraction = chrx_min_fraction,
    chry_min_fraction = chry_min_fraction,
    runtime_attributes = default_runtime_attributes
  }

  Boolean xy = ugc_wgw_chrom_coverage.inferred_karyotype == "XY"
  # dipcall's male mode wants chrY in hap1 and non-PAR chrX in hap2: true of trio-binned haplotypes, arranged for the others
  Array[File] dipcall_fastas = if length(ugc_wgw_chrom_coverage.regrouped_fastas) == 2
    then ugc_wgw_chrom_coverage.regrouped_fastas
    else hap_fastas

  call Dipcall.ugc_wgw_dipcall { input:
    sample_id = sample_id,
    ref_fasta = ref_map["asm_fasta"],  # !FileCoercion
    ref_index = ref_map["asm_fasta_index"],  # !FileCoercion
    hap1_fasta = dipcall_fastas[0],
    hap2_fasta = dipcall_fastas[1],
    xy = xy,
    par_bed = ref_map["asm_par_bed"],  # !FileCoercion
    threads = dipcall_threads,
    mem_gb = dipcall_mem_gb,
    runtime_attributes = default_runtime_attributes
  }

  call Pav.ugc_wgw_pav { input:
    sample_id = sample_id,
    ref_fasta = ref_map["asm_fasta"],  # !FileCoercion
    hap1_fasta = hap_fastas[0],
    hap2_fasta = hap_fastas[1],
    trio = trio,
    threads = pav_threads,
    align_threads = pav_align_threads,
    mem_gb = pav_mem_gb,
    runtime_attributes = default_runtime_attributes
  }

  # masking: Flagger's unreliable regions in reference coordinates, and dip.bed without them
  if (defined(flagger_bed)) {
    scatter (i in range(length(hap_fastas))) {
      call Masking.ugc_wgw_liftover_unreliable { input:
        sample_id = sample_id,
        haplotype = haplotype_labels[i],
        flagger_bed = select_first([flagger_bed]),
        paf = ugc_wgw_minimap2_align_asm.paf[i],
        ugc_wgw_container_registry = ugc_wgw_container_registry,
        runtime_attributes = default_runtime_attributes
      }
    }
    call Masking.ugc_wgw_reliable_regions { input:
      sample_id = sample_id,
      dip_bed = ugc_wgw_dipcall.bed,
      unreliable_beds = ugc_wgw_liftover_unreliable.unreliable_bed,
      ref_index = ref_map["asm_fasta_index"],  # !FileCoercion
      runtime_attributes = default_runtime_attributes
    }
  }

  call VcfFinish.ugc_wgw_vcf_finish as finish_dipcall { input:
    vcf = ugc_wgw_dipcall.vcf,
    sample_id = sample_id,
    caller = "dipcall",
    ref_name = ref_name,
    ref_fasta = ref_map["asm_fasta"],  # !FileCoercion
    ref_index = ref_map["asm_fasta_index"],  # !FileCoercion
    unphase = !trio,
    rename_sample_from = "syndip",
    reliable_bed = ugc_wgw_reliable_regions.reliable_bed,
    runtime_attributes = default_runtime_attributes
  }

  call VcfFinish.ugc_wgw_vcf_finish as finish_pav { input:
    vcf = ugc_wgw_pav.vcf,
    sample_id = sample_id,
    caller = "pav",
    ref_name = ref_name,
    ref_fasta = ref_map["asm_fasta"],  # !FileCoercion
    ref_index = ref_map["asm_fasta_index"],  # !FileCoercion
    unphase = !trio,
    reliable_bed = ugc_wgw_reliable_regions.reliable_bed,
    runtime_attributes = default_runtime_attributes
  }

  call CallSpectra.ugc_wgw_asm_call_spectra { input:
    sample_id = sample_id,
    dipcall_vcf = finish_dipcall.vcf_out,
    dipcall_vcf_index = finish_dipcall.vcf_index,
    pav_vcf = finish_pav.vcf_out,
    pav_vcf_index = finish_pav.vcf_index,
    runtime_attributes = default_runtime_attributes
  }

  # concordance with the read-based calls, inside the reliable region (dip.bed without a Flagger mask)
  if (defined(read_small_variant_vcf)) {
    call Concordance.ugc_wgw_asm_concordance { input:
      sample_id = sample_id,
      asm_vcf = finish_dipcall.vcf_out,
      asm_vcf_index = finish_dipcall.vcf_index,
      read_vcf = select_first([read_small_variant_vcf]),
      read_vcf_index = select_first([read_small_variant_vcf_index]),
      region_bed = select_first([ugc_wgw_reliable_regions.reliable_bed, ugc_wgw_dipcall.bed]),
      region_name = if defined(flagger_bed) then "reliable.bed" else "dip.bed",
      runtime_attributes = default_runtime_attributes
    }
  }

  call Provenance.ugc_wgw_manifest_write { input:
    ugc_wgw_version = ugc_wgw_version,
    stage = "assembly_call",
    subject_type = "sample",
    subject_id = sample_id,
    upstream_workflow_name = "dipcall v0.3 and PAV 3.0.0.dev26 on the assembly stage's contract",
    runtime_attributes = default_runtime_attributes
  }

  output {
    Array[String] haplotype_labels_out = haplotype_labels
    Array[File] asm_bams = ugc_wgw_minimap2_align_asm.bam
    Array[File] asm_bam_indices = ugc_wgw_minimap2_align_asm.bam_index
    Array[File] asm_pafs = ugc_wgw_minimap2_align_asm.paf
    File chrom_coverage = ugc_wgw_chrom_coverage.coverage_tsv
    File karyotype_json = ugc_wgw_chrom_coverage.karyotype_json
    String inferred_karyotype = ugc_wgw_chrom_coverage.inferred_karyotype
    File regroup_log = ugc_wgw_chrom_coverage.regroup_log
    File dipcall_vcf = finish_dipcall.vcf_out
    File dipcall_vcf_index = finish_dipcall.vcf_index
    File dipcall_vcf_stats = finish_dipcall.stats
    File dipcall_raw_vcf = ugc_wgw_dipcall.vcf
    File dipcall_bed = ugc_wgw_dipcall.bed
    Array[File] dipcall_hap_beds = ugc_wgw_dipcall.hap_beds
    File pav_vcf = finish_pav.vcf_out
    File pav_vcf_index = finish_pav.vcf_index
    File pav_vcf_stats = finish_pav.stats
    File call_spectra = ugc_wgw_asm_call_spectra.call_spectra
    File pav_raw_vcf = ugc_wgw_pav.vcf
    File pav_raw_vcf_index = ugc_wgw_pav.vcf_index
    Array[File] pav_callable_beds = ugc_wgw_pav.callable_beds
    Array[File]? unreliable_beds = ugc_wgw_liftover_unreliable.unreliable_bed
    File? reliable_bed = ugc_wgw_reliable_regions.reliable_bed
    File? reliable_summary = ugc_wgw_reliable_regions.summary_tsv
    File? concordance = ugc_wgw_asm_concordance.concordance_tsv
    String ugc_wgw_workflow_version = ugc_wgw_version
    File ugc_wgw_manifest = ugc_wgw_manifest_write.ugc_wgw_manifest
  }
}
