# `assembly`

A phased de novo assembly with hifiasm and its quality control.

| | |
|---|---|
| Subject | per sample |
| Modes | assembly |
| All stages | [08 Stages](../08-stages.md) |

hifiasm per sample, trio-binned with the parents' reads when the pedigree
allows, then each haplotype converted to FASTA with statistics. Ported from
PacBio's HiFi-human-assembly-WDL v1.0.2. The stage ends at an **output
contract** that every later stage reads: `assembler_used`,
`assembler_version`, `haplotype_labels`, `hap_fastas`, statistics and graphs.
The assembler behind it is named by the `assembler` input; only `hifiasm`
exists, and any other value fails at the start of the run.

```mermaid
flowchart TD
  subgraph child["per child uBAM"]
    c[/"uBAM"/]:::data --> cf["samtools fasta"]
  end
  subgraph trio["if both parents' reads are given"]
    subgraph pf["per parent uBAM"]
      p[/"parent uBAM"/]:::data --> pfa["samtools fasta"]
    end
    pfa --> depth{"each parent FASTA under trio_low_depth_gb, 90 GB?"}:::decision
    depth -->|"yes"| y1["yak count -b0, 70 GB, hifiasm -c1 -d1"]
    depth -->|"no"| y2["yak count -b37, 50 GB, hifiasm -c2 -d5"]
  end
  cf --> asm["hifiasm 0.25.0 --telo-m CCCTAA: hap1 and hap2 contig graphs"]
  y1 --> asm
  y2 --> asm
  subgraph hap["per haplotype; hap1 is the father's, hap2 the mother's in trio mode"]
    fa["gfatools gfa2fa, bgzip, calN50"]
  end
  asm --> fa
  fa --> out[/"hap_fastas, haplotype_labels, assembler_version"/]:::result
  fa --> fl["HMM-Flagger: reads to hap1+hap2 (minimap2 map-hifi), correct_bam, bam2cov, hmm_flagger"]
  cf --> fl
  fl --> fb[/"flagger_bed, flagger_summary"/]:::data
  fb --> cq["contig QC: seqtk comp, gap, telo per haplotype + Flagger bp per contig"]
  fa --> cq
  cq --> cqo[/"contig_qc"/]:::result
  y1 --> te["yak trioeval per haplotype (trio only)"]
  y2 --> te
  fa --> te
  classDef data fill:#e9edf1,stroke:#8a8f95,color:#22272e
  classDef result fill:#d5ead4,stroke:#3f7d3c,color:#17361a
  classDef decision fill:#e3dcf5,stroke:#6b56b0,color:#241a4d
  classDef box fill:none,stroke:#8a8f95
  class child,trio,pf,hap box
```

Text version: `samtools fasta` per child uBAM (and per parent uBAM in trio
mode, followed by `yak count` per parent with depth-dependent parameters);
hifiasm with `--telo-m CCCTAA`, trio-binned when the yak tables exist; then
`gfatools gfa2fa` and calN50 per haplotype. The outputs are the contract.
Two QC steps read it: HMM-Flagger aligns the sample's own reads back to the
diploid assembly (hap1 and hap2 together, minimap2 `-x map-hifi -I8g`),
cleans the alignment with secphase's `correct_bam` (primary alignments,
reads and alignments of 5 kb or more, divergence at most 0.1), and models
the per-base coverage into `Err`, `Dup`, `Hap` and `Col` regions
(`flagger_bed`, with the shares in `flagger_summary`); trio-binned samples
also get `yak trioeval` per haplotype against the parental k-mer tables,
with the thresholds the binning used (switch and Hamming error rates).
Since 0.12.0 a small task (seqtk, in the Flagger image) writes one row per
contig with its length, GC, N bases and runs, telomeric bases at each end
(the assembly's motif) and Flagger's bases per label (`contig_qc`).

| Driver fills | You may override |
|---|---|
| `sample_id`, `hifi_reads`; in trio mode `father_id`, `mother_id`, `father_hifi_reads`, `mother_hifi_reads` | `assembler` (`hifiasm`), `hifiasm_telomere_motif` (`CCCTAA`; empty = off), `hifiasm_extra_params`, `hifiasm_threads` (48), `hifiasm_mem_gb` (288), `yak_threads` (24), `yak_mem_gb`, `yak_params`, `hifiasm_trio_params`, `trio_low_depth_gb` (90), `keep_read_alignment` (false: the reads-to-assembly BAM is tens of GB per sample), `flagger_align_threads` (32), `flagger_align_mem_gb` (64), `flagger_threads` (16), `flagger_mem_gb` (32), `trioeval_c`, `trioeval_d` (the binning's) |

Resources: hifiasm 48 CPUs and 288 GB for up to a day at 30x (about 135 GB
observed); yak 24 CPUs and 70 or 50 GB per parent; the Flagger alignment
32 CPUs and 64 GB (the whole diploid assembly indexed at once), bam2cov and
hmm_flagger 16 CPUs and 32 GB; trioeval 8 CPUs and 48 GB. The site's default
wall time of three days covers hifiasm. The yak tables are not exported; the
call cache re-uses them for siblings. Flagger's read alignment is an
intermediate unless `keep_read_alignment`.

**The hifiasm command is frozen.** miniwdl keys its call cache on a task's
command and inputs, so anything that changes the hifiasm command line (the
motif does; the thread count does not) costs every sample's assembly again,
up to a day each. 0.11.0 made the one change (`--telo-m`, the version file);
later releases change only the stages after the contract, and
`assembly_call` with `--any-version` builds on the assemblies of the previous
version (chapter 10).

Outputs: per-haplotype FASTAs and statistics, graphs, Flagger's BED and
summary, the per-contig QC table, trioeval per haplotype (trio), `trio`,
`hap1_parent`, `hap2_parent` (chapter 07). `haplotype_labels` gives the
order of every per-haplotype array: `hap1`, `hap2`, or `pat`, `mat` when
trio-binned. The QC flags the summary derives from these (chapter 11,
`asm_*` thresholds) feed
`ugc-wgw summary --pass-list`, the sample list a pangenome cohort is frozen
from (chapter 06).
