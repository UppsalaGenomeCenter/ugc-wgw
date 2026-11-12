# 13 QC flags

Every flag the pages can show, in one place: what raises it, how serious
it is, which threshold moves it, and what to do. Flags are advisory checks
made by the driver when a page is written, not by the workflows: a flagged
sample has all its results, and nothing downstream is blocked by a flag
(the one gate is the pass list you choose to freeze a cohort from, chapter
06).

## How flags work

- **Severity.** `fail`: the result should not be used as it is. `warn`:
  look before using. `info`: a note. A sample's QC status is its worst
  flag (`pass` without any).
- **Where a flag is decided.** Per sample, from the sample's own results,
  or per cohort, because it needs the other members (identity, the
  cohort-relative outliers). Cohort flags join the member's row on the
  cohort page and in `ugc-wgw export`; they are not written into the
  sample's own digest.
- **Thresholds** come from `summary_thresholds` in `config.json`, a site
  profile, or `--threshold KEY=VALUE` on `summary` and `export`; unknown
  keys are refused. Flags are evaluated again every time a page is
  written, so a changed threshold needs no recomputation.
- **Finding flagged members.** The cohort page's member finder (the box
  above the list on the left) takes a flag name or a status: `fail`,
  `related_hub`, `tr_outlier_count_high`.

## Per-sample flags

| Flag | Severity | Raised when | Threshold (default) | What to do |
|---|---|---|---|---|
| `low_depth` | fail | Mean depth is below the floor. | `depth_mean_min` (20) | Top up the sample or accept it knowingly; cohort statistics count it like any member. |
| `low_mapped` | fail | The share of mapped reads is below the floor. | `mapped_read_percent_min` (95) | Suspect the wrong species or a contaminated library; check the reads section. |
| `low_read_quality` | fail | Median read quality is below the floor. | `read_quality_median_min` (25) | A sequencing-run problem; ask the facility. |
| `sex_not_inferred` | fail, or warn | No source (sheet, coverage, somalier) could sex the sample, so sawfish and TRGT ran with the XX default. A warning when only the coverage inference is missing. | `sex_unknown_fails` (1; 0 makes it a warning) | Add the sex to the sample sheet and rerun the sample's per-sample stage. |
| `sex_mismatch` | warn | The sex sources disagree; the message names each source and what the run used. | none | A sample swap or a sheet error until shown otherwise: check the identity section of the cohort page. |
| `tool_status` | warn | A targeted caller reported a problem: kivvi without reads, Paraphase regions failed for coverage, mitorsaw without haplotypes. | none | Usually low coverage at that region; read the message, the rest of the sample is unaffected. |
| `mito_unexplained_high` | warn | mitorsaw's haplotypes leave too much of the fingerprint reads unexplained: contamination or a nuclear copy of mitochondrial DNA. | `mito_unexplained_max` (0.02) | Check `het_outlier` and `related_hub` for the same member. |
| `imprinting_anomalous` | warn | An imprinting region lost its expected allele-specific methylation without a QC warning of MethBat's own; the message names the regions. | none | A finding to follow up, not a sample problem: look at the region in the member's methylation section. |
| `external_input` | info | The sample entered as an aligned BAM, or as a BAM with its calls: that part was made outside the project. | none | Nothing; the run manifest names the files. |
| `messages` | info | The workflow wrote messages (for example "fail reads included"). | none | Read them at the bottom of the sample page. |
| `parser_warnings` | info | A result file could not be summarised. | none | Report it: the file's format changed or the run is incomplete. |

## Assembly flags (samples with an assembly)

| Flag | Severity | Raised when | Threshold (default) | What to do |
|---|---|---|---|---|
| `asm_flagger_err` | fail | Flagger marks too much of the assembly as erroneous, duplicated or collapsed. | `asm_flagger_err_pct_max` (10) | Do not put the sample into a pangenome cohort. |
| `asm_hap_length` | fail | A haplotype's total length is outside the expected range. | `asm_hap_length_min`, `asm_hap_length_max` (2.5 and 3.5 Gb) | The same. |
| `asm_switch_error` | fail | A trio-binned haplotype's switch-error rate is too high. | `asm_switch_error_max` (0.05) | Check that the parents are the right samples. |
| `asm_reliable_fraction` | fail | Too little of the reference is reliably assembled. | `asm_reliable_fraction_min` (0.5) | Coverage is the usual cause. |
| `asm_het_loss` | fail | Too many of the sample's read-based heterozygous sites are homozygous or missing in the assembly: collapsed haplotypes. | `asm_het_loss_max` (0.2) | Do not use the assembly's genotypes. |
| `asm_concordance_low` | fail | The assembly disagrees with the sample's own read-based SNVs; a swap scores 0.5 to 0.7. | `asm_genotype_concordance_min` (0.90) | Treat as a sample swap between the assembly and the reads. |
| `asm_karyotype_mismatch` | warn | The sheet's sex and the assembly's karyotype disagree. | none | As `sex_mismatch`. |
| `asm_karyotype_unknown` | warn | Too little non-PAR chrX assembled to tell: the sample joins no sex stratum of the pangenome frequencies. | none | Nothing, unless strata matter. |

The assembly `fail` flags are what `summary --pass-list --pass-on assembly`
gates on.

## Cohort flags: identity (`cohort_qc`)

| Flag | Severity | Raised when | Threshold (default) | What to do |
|---|---|---|---|---|
| `duplicate` | fail (both members) | A pair is the same individual by relatedness: a re-delivered sample or monozygotic twins. | `duplicate_relatedness_min` (0.9) | Keep one; `--pass-on identity` drops the later one in the sheet. |
| `pedigree_mismatch` | fail (child and parent) | A parent-child pair the sheet declares is not related like one. | `pedigree_relatedness_min`, `pedigree_relatedness_max` (0.35, 0.65) | A sheet error as often as a swap: check the pedigree first. |
| `related_hub` | fail | A member is related, without a pedigree link, to many members: contamination makes a sample look related to everyone. Its partners are not flagged on its account. | `related_hub_min` (10 members) | Exclude the member; expect `het_outlier` on it too. |
| `related` | warn | An undeclared pair is related (neither member a hub). | `related_relatedness_max` (0.2) | Declare the relation in the sheet, or keep one member for unrelated-cohort statistics. |
| `het_outlier` | warn | The member's heterozygosity ratio is far from the cohort's. | `het_ratio_mad_max` (3 scaled MADs) | Contamination or a sample mix when high; consanguinity or ancestry when low. |

`sex_mismatch` is also raised here for a member whose own run predates the
somalier extract, when the stage's chrX call disagrees with its sheet or
coverage.

## Cohort flags: relative to the other members

These three compare a member with its cohort, so they describe the member
*in this cohort*: the same sample can be flagged in one cohort and not in
another.

| Flag | Severity | Raised when | Threshold (default) | What it usually means |
|---|---|---|---|---|
| `methyl_outlier` | warn | The member has many more methylation events against the cohort background than the others (`cohort_methyl`). | `methyl_events_mad_max` (3 MADs above the median) | Check `cell_composition_outlier` first: another mixture of cells differs at many regions. |
| `cell_composition_outlier` | warn | The member's cell-type composition is far from the cohort's median composition (`cohort_methyl`'s deconvolution). | `cell_composition_mad_max` (4 scaled MADs) and `cell_composition_distance_min` (0.25) | Another tissue or sample type, a cell line, contamination, or an unusual blood picture. |
| `tr_outlier_count_high` | warn | The member has many more expansion outliers at ordinary tandem-repeat loci than the others (`cohort_tr`). | `tr_outlier_count_mad_max` (3 scaled MADs above the median) | An ancestry that is rare in the cohort, or a lower depth than the others; only then a sample problem. |

No flag is raised for a class of the disease-locus screen of `cohort_tr`:
it is a research screen, and a QC status must not carry a finding.

## Cohort warnings (the whole cohort, not a member)

| Warning | Where | Raised when | What to do |
|---|---|---|---|
| `mixed_tools` | Cohort page, Overview | The members were made with more than one aligner or small-variant caller (chapter 05, "Choosing pblab's tools"). | Keep a cohort on one set unless you are comparing: joint calls and frequencies can show the difference as a batch effect. |
| members below `pg_member_call_rate_min` (0.90) | Cohort page, pangenome section | A member has little assembled sequence in the graph. | A signal for the next cohort freeze, not a sample flag. |
