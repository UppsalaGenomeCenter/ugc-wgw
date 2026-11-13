# Changelog

Release notes of the public copy. Each version corresponds to a verified
offline bundle. An entry opens with the upstream tag it is built on and a
short summary, with anything an existing project has to do; the rest is
under Details.

## 0.21.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin (0.6.0), the
images and the reference entries unchanged.

The `assembly` stage keeps more of what it computes, and a new command
chooses the members of a pangenome graph. The stage's two expensive tasks
change, so **a sample assembled with an earlier version is assembled
again at 0.21.0**; this is the one such change, made before the stage is
run on whole cohorts. The read-based stages are untouched.

<details>
<summary>Details</summary>

- **`assembly` keeps the primary contigs.** hifiasm writes them beside
  the two haplotypes in the same run; they are now outputs
  (`primary_fasta`, `primary_stats`): one more contiguous assembly per
  person that switches between the haplotypes. HiFi-only samples; hifiasm
  writes none when it trio-bins.
- **The reads on the sample's own assembly, kept, with methylation.** The
  stage aligns the reads back to the two haplotypes for HMM-Flagger. They
  are now taken from the read BAMs as FASTQ with their `MM`/`ML` tags, so
  the alignment has base qualities and 5mC; it is an output by default
  (`read_alignment_bam`, about the size of the read BAMs;
  `keep_read_alignment = false` turns that off), without minimap2's
  secondary alignments. The placements, and Flagger's result, are the
  same as before. `methbat pileup` writes `asm_cpg_pileup_bed` from it:
  5mC per CpG in assembly coordinates, including sequence the reference
  lacks.
- **`ugc-wgw cohort rank <cohort_id>`.** Orders a cohort's members for a
  pangenome graph by the structural variants each adds to the ones before
  it, from the read-based SV merge the cohort already has: common alleles
  first, relatives of an earlier pick held back, chosen members (trio
  children, a control) placed first on request. It writes the ranking,
  saturation curves (the ranking against members drawn at random, per
  frequency class), how much of a member outside the first 20 to 400 is
  already covered, and a page. `--emit K --pass-list FILE` writes the
  first K that passed assembly QC, so graph cohorts of different sizes
  are nested.
- **Measured and not adopted:** pbweave for the alignment to the own
  assembly. It is 35 times faster there, but where the two haplotypes are
  identical it sends every read to the first one, which distorts both
  Flagger and a methylation pileup (guide chapter 15).

</details>

## 0.20.1

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin (0.6.0), the
images and the reference entries unchanged.

A documentation release: the frames around grouped steps in the guide's
diagrams lose Mermaid's default pale-yellow fill. Workflows, driver and
results are those of 0.20.0; a site that runs 0.20.0 has no reason to
install this one, and a site still on 0.19.0 can install it in place of
0.20.0 with the same commands.

<details>
<summary>Details</summary>

- **Diagram frames.** A frame (a scatter, a conditional section, a mode,
  a machine) now has a grey border and no fill, in the hand-written
  diagrams and in the generated call graphs, on light and dark pages.
  `scripts/wdl-graph.py` writes the class for every section it draws, and
  its test checks that none is left out.

</details>

## 0.20.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin (0.6.0), the
images and the reference entries unchanged.

Cell-type composition becomes a result of `cohort_methyl`: a table, a
section on the cohort and sample pages, and a flag for a member whose
composition is unlike the cohort's. The guide is reorganised (a page per
stage, a QC flag reference, a chapter on reading the pages) and its
diagrams are coloured. An existing project runs `ugc-wgw upgrade
--install`, then `submit --stage cohort_methyl --any-version` once per
cohort to get the composition; commands on a project whose samples ran at
an earlier version need `--any-version`, and the driver now says so.

<details>
<summary>Details</summary>

- **Cell-type composition.** `cohort_methyl` runs `methbat deconvolve` for
  every member by default (`run_deconvolve`; a minute and 2.5 GB per
  member, measured on 50 genomes) against the shipped atlas of 39 cell
  types, and writes `<cohort>.methyl.<ref>.deconvolution.tsv`: a row per
  member, a column per cell type, with MethBat's residual and raw
  fraction sum. The cohort page shows a ribbon with one column per member,
  the cohort's composition per cell type, the members' spread and the
  members far from the cohort; a member's page shows its fractions beside
  the cohort's medians; `ugc-wgw export --domains methyl` copies the
  table. New advisory flag `cell_composition_outlier` (thresholds
  `cell_composition_mad_max` 4 and `cell_composition_distance_min` 0.25):
  another tissue or sample type, a cell line, contamination or an unusual
  blood picture; relative to the cohort, so a saliva cohort is judged
  against saliva. Checked on both sides: whole-blood samples come out as
  granulocytes, T cells and monocytes, and the Genome in a Bottle
  lymphoblastoid lines as 59 to 78% B cells. An estimate for research,
  not a blood count: fractions of a few percent in unrelated tissues are
  the method's noise. The per-sample digests are untouched.
- **pblab on production genomes** (measured with 0.19.0 on 50 genomes at
  about 20x): the whole joint run cost 20.5 core-hours per sample against
  about 83 with pbmm2 and Parabricks; pbweave 4 min 15 s and 1.7
  core-hours per sample, pbdv 5 min 52 s on 16 cores and one GPU. A
  fourth comparison on chr20 (pbweave with DeepVariant) scores and
  differs exactly as pbweave with pbdv: what the route changes is the
  aligner's, and pbdv reproduces DeepVariant. The defaults stay pbmm2 and
  DeepVariant.
- **Driver fixes.** A stage blocked on runs that exist only at another
  version says so and names `--any-version`. `submit --rerun` started a
  stage twice when the stage had no run at the running version; it starts
  once. A profile's `input_thresholds` reach the project (they were
  accepted and dropped).
- **Guide.** One page per stage under `docs/guide/stages/` with chapter 08
  as the index; chapter 13, every QC flag in tables (what raises it, the
  threshold, what to do); chapter 14, reading the cohort and sample pages
  section by section; chapter 15, measurements and evaluations, taken out
  of the resources chapter; an operator and an analyst path in the index.
  Diagrams carry colour as well as shape (per-sample stage blue, cohort
  stage amber, optional dashed, data grey, result green, decision violet;
  in the generated call graphs this pipeline's tasks teal, upstream's
  grey), readable on light and dark pages.

</details>

## 0.19.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin (0.6.0)
unchanged; one new image of ours, `pblab` (PacBio's pblab release 20261021
on upstream's base image); the reference entry `ugc-wgw-tr` moves to 0.2.0.

PacBio's preview aligner and small-variant caller, pbweave and pbdv, can be
chosen per project in the place of pbmm2 and DeepVariant; the defaults do
not change. The tandem-repeat stage is corrected after its first cohort of
864 genomes, and a finished stage can be run again with `submit --rerun`.
An existing project runs `ugc-wgw upgrade --install` (the reference map
names a new table) and, where `cohort_tr` has run, `submit --stage
cohort_tr --any-version` once.

<details>
<summary>Details</summary>

- **pblab's tools as a project choice.** PacBio's pblab release 20261021
  (every tool marked Preview) ships pbweave, an aligner, and pbdv, a
  reimplementation of DeepVariant 1.10.0's PacBio model. New project key
  `aligner: pbmm2|pbweave` (`init --aligner`) and two more values of
  `deepvariant`, `pbdv` (one GPU) and `pbdv-cpu`; either works alone, and
  the profile `pblab` sets pbweave and pbdv. They run inside `singleton`
  and `upstream` with the same outputs and names, so every later tool
  reads their BAM and VCF unchanged; a default project's inputs and call
  cache do not change. Each run's manifest records the aligner and caller
  (`tools`), the sample page names them, and a cohort that mixes sets
  shows the warning `mixed_tools`. Measured on chr20 of the GIAB trio:
  against the HG002 v5.0q truth, SNV F1 0.9985 (DeepVariant), 0.9985
  (pbdv on the same alignments) and 0.9984 (pbweave + pbdv), indel 0.950
  for all three, SV 0.928; pbdv's calls differ from DeepVariant's in a
  handful of records per 94,000. With pbweave, calls agree inside the
  benchmark regions (99% of small variants, SV F1 0.98) and differ in the
  centromeric and duplicated sequence no benchmark scores (54% of small
  variants shared there), so one cohort should stay on one aligner. Speed,
  from PacBio's figures for a 60 Gb genome: alignment 30 s instead of 15
  minutes on 64 cores, calling under 3 minutes on one GPU. The defaults
  stay pbmm2 and DeepVariant until whole genomes have been compared.
  `tests/smoke/run.sh eval pblab` and `tests/smoke/compare_runs.py` make
  that comparison.
- **`cohort_tr` after its first production cohort.** The cohort page
  listed every screen row with any class, note or borderline value (30,482
  rows for 864 members); it now lists the pathogenic and intermediate rows
  and the calls in doubt, and counts the rest per locus. A locus can be
  `not_assessable`: ADTKD_MUC1's thresholds count a cytosine tract inside
  one repeat unit, not units, so its alleles are reported and nothing is
  classed (reference table columns `assessable`, `assessable_note`). The
  per-locus table separates `below_pathogenic` (loci without a benign
  maximum) from "above benign" and marks a locus
  `threshold_or_frame_suspect` when more than a fifth of the alleles fall
  above or below the benign range. `hom_call` is written only on a
  homozygous pathogenic or intermediate call. `tr_outlier_count_high` uses
  scaled MADs (62 flagged members of 864 become 29); the count follows
  ancestry and depth before it says anything about a sample. Measured
  cost at 864 members: half a core-hour.
- **`ugc-wgw submit --stage S --rerun`** starts a new attempt of one stage
  for the selected subjects whose latest run succeeded (a changed
  `stage_inputs` value, reference or threshold). A per-sample stage needs
  `--samples`; the call cache still replays unchanged tasks; later stages
  that already succeeded are named in order and left alone.
- **The cohort page's side list holds every member** (it stopped at 200),
  each with a dot for its QC status, under a box that narrows the list by
  id, status or flag name.
- Release notes: each version is now a "Built on" paragraph, a short
  summary with anything an existing project has to do, and the details
  folded.

</details>

## 0.18.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin (0.6.0)
unchanged; two new reference entries (`ugc-wgw-tr`, `ugc-wgw-truth`); one
new image, used by the smoke test only.

A new optional stage, `cohort_tr`, screens the 73 disease loci of the TRGT
catalog per member against STRchive's thresholds, gives every locus its
allele-length distribution and lists expansion outliers across the cohort;
it is a research screen, not a diagnostic result. `tests/smoke/run.sh
validate` scores the smoke's HG002 against the Genome in a Bottle v5.0q
benchmark on chr20 (SNV F1 0.9985, indel 0.950, SV 0.928). An existing
project runs `ugc-wgw upgrade --install` first, since the stage reads two
new reference-map keys; no per-sample run or digest is rebuilt.

<details>
<summary>Details</summary>

- **A new stage, `cohort_tr`**, from the merged TRGT VCF and trgt-lps table
  of `cohort_merge` alone (no per-sample task, optional: nothing waits for
  it). The 73 named disease loci of the TRGT catalog are screened per
  member against STRchive's thresholds (STRchive v2.4.4, CC BY 4.0): copies
  of the pathogenic motif per allele, a class (`normal`, `intermediate`,
  `pathogenic`, the gaps by name), the fit to the locus's inheritance,
  notes such as an expansion of another motif. Every locus of the catalog
  gets its allele-length distribution and spectrum, an in-house resource,
  and the members whose longer allele is far above the cohort's and rare
  in it are listed as expansion outliers (at least 50 bp and 5 scaled MADs
  above the median, never less than five repeat units, a common long
  allele listing nobody). The cohort page gains a Tandem repeats section
  and the advisory flag `tr_outlier_count_high`; `ugc-wgw export --domains
  tr` copies the tables. A research screen, not a diagnostic result: no
  flag is raised for a class. The per-sample digests are untouched, so
  nothing is rebuilt. trgt-denovo is not part of the stage.
- **`tests/smoke/run.sh validate`** scores the smoke's HG002 against the
  Genome in a Bottle HG002 v5.0q benchmark on chr20, small variants with
  aardvark and structural variants with truvari (bench, then refine), for
  the standalone and the joint route, and `tests/smoke/VALIDATION.md`
  keeps the rows per release. First scores (DeepVariant 1.10.0, sawfish
  2.3.0, 48.7x): SNV F1 0.9985 for DeepVariant's own call set and 0.9977
  for the member's slice of the GLnexus joint VCF, indel F1 0.950 and
  0.939, SV F1 0.928. The chr20 truth slices ship as a reference entry.

</details>

## 0.17.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin (0.6.0) and the
reference entries unchanged; no new image.

Fixes from the first production cohort of 864 members, and sawfish 2.3.0 as
the default structural-variant caller (calls identical to 2.2.1's on three
genomes, at a fraction of the time and memory). Identity QC gains
`related_hub`, which fails a contaminated sample instead of excluding
everyone it resembles; haplogroups are counted for runs that predate 0.15.0;
`upgrade` moves the reference map too, and `submit` refuses a map that lacks
a key. Existing results are untouched: run `ugc-wgw upgrade --install`, then
`summary` and `export` again.

<details>
<summary>Details</summary>

- **sawfish 2.3.0 is the default** for structural-variant discovery and
  calling (upstream v4.0.0 pins 2.2.1). On three whole genomes the call
  sets were identical to 2.2.1's, record for record, with `discover` at
  6 minutes and 3.6 GB instead of 28 minutes and 68 GB. Every per-sample
  run started after the upgrade uses it; existing results are untouched.
- **Identity: `related_hub`.** A member related, without a pedigree link,
  to ten or more others (`related_hub_min`) is a contaminated or mixed
  sample, not kin: it fails, and its pairs no longer count against its
  partners. On an 864-member cohort two such samples had excluded about
  600 members from the identity pass list.
- **Mitochondria.** For members whose runs predate 0.15.0 the most
  abundant mitorsaw haplotype is the major one, so the haplogroups that
  `cohort_merge` classifies for them are counted on the cohort page and in
  `mito.freq.tsv`.
- **Methylation.** `ugc-wgw export --domains methyl` copies the comparison
  files, and the cohort page says where each comparison's differing regions
  lie (chrX, chrY, autosomes): the default `XX:XY` pair is a check as much
  as a result.
- **Moving a project to a new version.** `ugc-wgw upgrade --install` now
  moves the rendered reference map with the engine paths; every command
  refuses a project whose paths or map belong to another installed
  version, and `submit` refuses a reference map that lacks a key the
  stage's workflow reads, naming the keys, instead of failing inside the
  engine.
- **Summaries.** `summary` and `export` write each per-sample digest as it
  is built and show progress, so an interrupted rebuild of a large cohort
  resumes where it stopped (`--jobs N` builds in parallel).

</details>

## 0.16.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin (0.6.0)
unchanged; no new image; one new reference entry, `ugc-wgw-methyl` 0.1.0
(MethBat's imprinting regions, a CpG-island gene map generated from GENCODE
50, the nanomix cell-type atlas).

A new stage, `cohort_methyl`, builds the cohort's methylation background
with MethBat and scores every member against it: outlier regions by gene,
group comparisons, an imprinting screen. Samples can now enter as an aligned
BAM or as pre-called VCFs, each kind with its own stage sequence. `ugc-wgw
upgrade --install` moves a project to another installed version. On first
use the state database moves to schema 4, and every per-sample digest
rebuilds once.

<details>
<summary>Details</summary>

- **A new stage, `cohort_methyl`.** MethBat 1.1.0's cohort layer, from
  upstream's image: a background profile of the cohort from the members'
  MethBat profiles, labelled by the resolved sex group and any sample-sheet
  columns the project names (`methyl_label_columns`); every member's 5mC
  pileup scored against it (Hypo- and Hyper-methylated regions, allele-
  specific methylation gained or lost, with z-scores); group comparisons
  (`XX:XY` by default); the imprinting regions of every member labelled by
  `methbat report` (made per sample since this release, inside the stage
  for older runs: no per-sample re-run); every region annotated with its
  genes and type from the shipped map; optional joint segmentation and
  cell-type deconvolution. The analyst's file is
  `<cohort>.methyl.<ref>.events.tsv`, every member's events sorted by gene;
  the cohort page gets a Methylation section (labels, events per member,
  the genes with the most members carrying an event, recurrent regions, the
  imprinting screen, comparisons), the sample page its imprinting table;
  `ugc-wgw export --domains methyl` copies the stage's tables. New advisory
  flags `imprinting_anomalous` and `methyl_outlier`. Digest schema 5: every
  per-sample digest rebuilds once.
- **Three entry kinds, derived from the sample sheet.** `hifi_reads` as
  before (`reads`); `aligned_bam` alone (`bam`): alignment is skipped, the
  BAM is checked against the installed reference (contigs and lengths,
  index, sort order, a read group named after the sample id, methylation
  tags) and used by every caller;
  `aligned_bam` with `small_variant_vcf` and `sv_vcf` (and `gvcf`; `vcf`):
  calling is skipped too, the new stage `intake` makes what phasing needs
  (checks, depth, the sex resolution, the VCFs' sample column renamed) and
  the sample enters at `downstream`. The sheet rules are explicit refusals;
  each member runs its own stage sequence in both cohort modes; joint mode
  needs a gVCF from every member unless `--allow-standalone-members`.
  `samples rebase <id> --from-run RUN [--with-vcfs]` registers a member's
  own earlier aligned BAM (and VCFs) as its inputs, so a tool change after
  alignment costs one `downstream` per member. Run manifests name the
  external files and the producing programs; the pages mark
  `external_input`. The database schema moves to 4 on first open.
- **`ugc-wgw upgrade --install <root>/current`** moves a project to another
  installed version (the engine paths and the reference map; the previous
  `config.json` is kept), and every command first checks that the project's
  paths are one installed version, the running driver's, stopping with the
  remedy instead of failing late in a stage.

</details>

## 0.15.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin (0.6.0)
unchanged; one new image, haplogrep3 3.3.2, and one image for evaluation
only, PacBio's sawfish 2.3.0_build1.

Each sample gets a diploid mitochondrial call set and a haplogroup
(mitorsaw's haplotypes collapsed, haplogrep3 for the classification), and
the cohort stages merge and count them. A result change: chrM leaves the
cohort small-variant VCF and its frequencies, and the mitochondrial sets
take its place. A tool version ahead of upstream can be tried on a copy of
the data first (`run.sh eval`, here sawfish 2.3.0). An existing cohort needs
no per-sample re-run; every per-sample digest rebuilds once.

<details>
<summary>Details</summary>

- **A diploid mitochondrial call set per sample.** mitorsaw writes one
  genotype entry per haplotype it finds (`0|1|1|1`), which no VCF consumer
  merges or counts. After mitorsaw, `singleton` and `upstream` collapse
  those entries to a diploid genotype (nallo's rule: every haplotype
  alternate `1/1`, mixed `0/1` heteroplasmic, `VAF` the read fraction),
  split the records at 50 bp into small variants and SVs, and write one
  haploid column per haplotype for haplogrep3, which assigns each a
  haplogroup against PhyloTree 17's forensic update (`phylotree-fu-rcrs@1.3`;
  PhyloTree 17 is also installed). The sample's haplogroup is its most
  abundant haplotype's; the sample page shows every haplotype with its
  haplogroup and quality. A re-run of an older sample replays everything
  but the two new tasks from the call cache.
- **chrM leaves the cohort small-variant VCF.** GLnexus no longer merges
  DeepVariant's diploid chrM records (they counted heteroplasmies as
  nuclear alleles): a result change for anyone who read chrM from the
  cohort set or from `cohort_freq`. Instead `cohort_merge` merges the
  members' collapsed mitochondrial sets (`cohort_mt_small_variant_vcf`,
  `cohort_mt_sv_vcf`; a member without the variant is `0/0`, stated in the
  header) and gathers every member's haplogroups into one table; a member
  whose run predates this release is collapsed and classified inside the
  stage from its mitorsaw outputs, so an existing cohort needs no
  per-sample re-run. `cohort_freq` counts the mitochondrial sets without
  sex strata (`AC_Het` heteroplasmic members, `AC_Hom` homoplasmic) and
  the frequency summary gains the `mt_small_variants` and `mt_svs`
  resources; the cohort page gets a Mitochondria block with the haplogroup
  counts and the heteroplasmic sites per member; `ugc-wgw export` fills
  `mito.tsv`'s haplogroup columns and writes `mito.freq.tsv` (haplogroups
  and macro-haplogroups per stratum). Digest schema 4: every per-sample
  digest rebuilds once.
- **Evaluating a tool version ahead of upstream.** `workflows/overrides/sawfish.wdl`
  is upstream's sawfish task file with the image as a task input, and the
  profile `eval-sawfish23` sets PacBio's sawfish 2.3.0_build1 for
  `singleton`, `upstream` and `cohort_call`; `tests/smoke/run.sh eval
  sawfish23` runs a second project on the smoke data with it and compares
  the SV sets (truvari) and the sawfish tasks' resources with the
  production project's. Production behaviour is unchanged; the override
  changes the sawfish task digest once, so a re-run of an older sample
  recomputes sawfish onward (guide 12).
- The identity pass list (`summary --pass-on identity`) keeps a pair the
  sheet declares (parent and child); only undeclared relatives and
  duplicates drop their later member.

</details>

## 0.14.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin (0.6.0)
unchanged; one new image, somalier 0.3.7, and one new reference bundle,
`ugc-wgw-qc` 0.1.0.

A new stage, `cohort_qc`, relates every pair of members with somalier,
checks the pedigree and the sex and projects ancestry; it never blocks the
stages after it. The sex a run uses is resolved from the sheet, the coverage
and chrX before sawfish and TRGT see it. `ugc-wgw export` writes the
cohort's numbers as TSV tables, one per domain. Releases can ship as delta
bundles that carry only the images the previous one lacked.

<details>
<summary>Details</summary>

- **A new stage, `cohort_qc`.** somalier relates every pair of members from
  one 200 KB genotype sketch per member (made by `singleton`/`upstream`
  since this release, or by the stage itself from the aligned BAM for
  members whose run predates it), checks the sheet's pedigree, calls sex
  from chrX heterozygosity and projects the members onto the 1000 Genomes
  superpopulations when that set is installed. Minutes at any cohort size,
  first of the cohort stages in both read modes, and the one stage that
  never blocks the stages after it. The cohort page gets an Identity
  section and the flags `duplicate` and `pedigree_mismatch` (fail),
  `related` and `het_outlier` (warn); `summary --pass-list --pass-on
  identity` writes the members to freeze from; `cohort freeze` warns about
  known duplicates.
- **The sex a run uses is resolved before sawfish.** The sample sheet's sex
  (now also `M`/`F`, `XX`/`XY`, `1`/`2`) is taken when given, else
  mosdepth's coverage inference, else somalier's chrX call; sawfish and
  TRGT take that value, so an XY sample whose coverage inference came back
  empty is no longer called as XX silently. Every source is kept: the
  manifest's `sex` block, the digest's identity block, the pages' sex
  columns; `sex_mismatch` names the disagreeing sources and
  `sex_not_inferred` fails when no source could sex the sample. A re-run of
  an older sample replays every other task from the call cache.
- **`ugc-wgw export --cohort ID`** writes the cohort's numbers as TSV files,
  one per domain (members, coverage, variants, phasing, Paraphase copy
  number, StarPhase diplotypes, HLA alleles, mitochondria, assembly, and
  copies of `cohort_qc`'s tables), with frequencies per stratum (`ALL`, sex,
  the ancestry projection) and an `export_manifest.json`. The cohort page's
  new "Targeted callers" section draws from the same computation. Digest
  schema 3 carries the per-gene rows; `mito_unexplained_high` warns on
  mitorsaw's unexplained fraction; StarPhase's `./.`, `NO_MATCH` and
  `Multiple` count as uncalled.
- **Delta bundles and one copy of each SIF per prefix.** `make-bundle.sh
  --images delta --base <previous manifest>` ships only the images the
  previous bundle did not carry (about 10 MB plus the changed SIFs instead
  of 27 GB); the installer keeps a content-addressed SIF store under the
  prefix and hardlinks every version's images into it (`--dedupe-sif`
  converts older versions; `--verify-only` reports what is shipped, in the
  store or missing). A delta needs its base installed; the full tar stays
  on the share for new sites.

</details>

## 0.13.1

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin (0.6.0) and the
image manifest unchanged.

Fixes from the first full-size joint run (864 samples, 1,772 runs). A task
killed at its time limit is now classified as a resource failure, so the
0.13.0 escalation fires. trgt-lps runs per merge shard instead of once over
the cohort VCF (13h52m before), the per-sample split uses 16 threads, and
the report pages scale to hundreds of members and thousands of runs.

<details>
<summary>Details</summary>

- **Time-limit kills are now classified as resource failures.** SLURM
  writes its "DUE TO TIME LIMIT" and "oom_kill" lines to the job's own
  `slurm-<jobid>.out` in the task directory, which the failure classifier
  did not read; a task killed at its wall time (exit 143) was reported as a
  tool failure and the 0.13.0 escalation did not fire. The classifier reads
  that file first (newest job first), then the sbatch log and every try's
  stderr and stdout. The accounting also takes job ids from those files, so
  a retried task's earlier jobs are counted.
- **trgt-lps runs per TRGT merge shard.** One trgt-lps over the 864-sample
  cohort VCF took 13h52m; `ugc_wgw_cohort_merge` now runs it on every
  shard's VCF (`trgt_lps_threads`, default 4) and a gather task joins the
  tables. Each table is put into its VCF's record order first, so the result
  equals the one-task table byte for byte.
- **The per-sample split uses 16 threads** by default (`split_threads`; 2h14m
  at 8 threads for 864 samples) and `split_mem_gb` scales with them.
- **The reports scale to hundreds of members and thousands of runs.** The
  SUPP histogram bins by share of the cohort with a cumulative line and
  cards for singletons and variants in every member; the distribution and
  sex strips become bands (whisker, interquartile box, median, density, and
  only outliers as points) from 150 members; tables of 100 rows or more fold
  into a collapsed block with a type-to-filter box, failures first; the run
  report's timeline above 200 runs is one lane per stage, with runs in
  flight as an area and failures as ticks. The 864-member cohort page went
  from 2.1 to 0.7 MB; the pages stay self-contained HTML without a charting
  library.
- Docs: the 864-sample run's task sizes and durations in chapter 12; chapter
  09 names `slurm-<jobid>.out` among the evidence and miniwdl's side-by-side
  retry layout (`stderr2.txt`, `work2/`).

</details>

## 0.13.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.6.0 and
the image manifest unchanged.

A task that runs out of wall time or memory is re-attempted at once with
that limit doubled, up to twice per task and resource, while every other
task of the run replays from the call cache. The request that worked is
learned per project, so later runs start there; `ugc-wgw resources
--learned` prints the learned values as policy rows. The miniwdl plugin
moves to 0.6.0.

<details>
<summary>Details</summary>

- **A task that ran out of wall time or memory is re-attempted with its
  limit raised, and the working request is learned.** The failure
  classifier now says which resource ran out; the driver re-attempts the
  run at once with that task's limit multiplied by `resource_retry_factor`
  (2), up to `resource_retry_max` (2) escalations per task and resource.
  The multiplier reaches the `ugc_wgw_resources` plugin through the
  environment and is applied after the policy rows, so the task's
  call-cache key is unchanged and every other task of the run replays from
  the cache; memory stops at the site cap, and a refused `sbatch` request is
  never escalated. What an escalated task finished with goes into the
  project's `.ugc-wgw/learned_resources.json` and every later run starts
  the task there, so one sample's lesson covers the rest of the cohort;
  `ugc-wgw resources` marks learned values and `--learned` prints them as
  rows of the site policy for promotion (chapters 06, 09, 11 and 12). The
  manifest records the escalations, the run report counts them and lists
  the applied changes, and `ugc-wgw retry` carries the last escalation on.
- Plugin ugc-wgw-miniwdl 0.6.0: `[ugc_wgw] escalations`, the driver's
  per-run factors and minimums, applied after the policy rows.

</details>

## 0.12.6

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.5.0 and
the image manifest unchanged.

The concatenation of the cohort VCF and its per-sample split no longer parse
every genotype. On an 800-sample `cohort_call` each took about four hours;
expect minutes and tens of minutes. The records are byte-identical to the
previous commands'.

<details>
<summary>Details</summary>

- **The cohort VCF's concatenation and per-sample split stop parsing every
  genotype.** On an 800-sample `cohort_call` both took about four hours on
  four threads, one after the other, between the last GLnexus shard and
  the first `downstream` run: `bcftools concat` parses every record unless
  told not to, and `bcftools +split` parses every record in one thread and
  writes one single-threaded bgzip stream per sample. The concat task gains
  `naive` (`bcftools concat --naive`: compressed blocks joined, sample
  columns checked), on for the GLnexus shards, the frequency shards and
  the phased small-variant shards. The split task runs `+split` once per
  contig that holds records, `split_threads` (default 8) at a time with
  its pieces on the node-local scratch, then one concatenation per sample
  in contig order; `split_mem_gb` defaults to 8 + 0.04 GB per sample.
  Records byte-identical to the one-pass commands' in both filter modes
  (chapter 12). Expect minutes and tens of minutes instead of four hours
  each.

</details>

## 0.12.5

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.5.0 and
the image manifest unchanged.

The pangenome's Cactus tasks map several haplotypes at once instead of one
at a time (`split_map_cores`, `chrom_map_cores`). On the chr20 smoke the
split went from 75 to 55 minutes and the chromosome task from 51 to 39, with
identical contig assignment. `ugc_wgw_chrom_coverage` no longer needs
`bgzip`, whose absence had turned the public CI of 0.12.4 red.

<details>
<summary>Details</summary>

- **The pangenome's Cactus tasks map several haplotypes at once.** Both
  tasks passed all their cores to minigraph's mapping, so Toil mapped one
  haplotype at a time while minigraph used a few cores of it (15 and 22 %
  CPU efficiency on the first whole-genome run; chr21 spent 4h43m of its
  5h51m mapping ten haplotypes' satellite-rich contigs one after another).
  New inputs `split_map_cores` (default: the task's cores, unchanged
  behaviour) and `chrom_map_cores` (default 4: four at once on 16 cores);
  the `assembly` profile gives the split a whole standard node with three
  mappings at once, each holding the reference graph's index (65–80 GB).
  On the chr20 smoke the split went from 75 to 55 minutes and the
  chromosome task from 51 to 39, with identical contig assignment; a
  chromosome graph's node count varies by about 0.1 % between any two
  builds, concurrency or not (chapter 12). The split keeps its three Toil
  step logs as a new output.
- `ugc_wgw_chrom_coverage` writes the regrouped haplotype copies with
  Python's gzip at level 4: 0.12.4's `bgzip` is not on every machine the
  tests run on (the public CI of v0.12.4 is red for that), and level 4 is
  25 times faster than 0.12.3's level 9 for a fifth more bytes.

</details>

## 0.12.4

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.5.0 and
the image manifest unchanged.

Two fixes from the first whole-genome run of the assembly track (five 20x
samples, every chromosome). The haplotype regrouping of a non-trio XY sample
takes minutes instead of an hour and a half, and `sex_not_inferred` no
longer warns on assembly-only cohorts. Chapter 12 records the run's task
sizes.

<details>
<summary>Details</summary>

- `ugc_wgw_chrom_coverage` writes the regrouped haplotype copies of a
  non-trio XY sample through `bgzip` instead of Python's gzip at level 9:
  an hour and a half per sample before, minutes now, with byte-identical
  output; the task gains a `threads` input (default 2).
- `sex_not_inferred` needs read-based data: every page of an assembly-only
  cohort carried the warning.
- Chapter 12 holds the run's sizes: the per-sample tasks on 20x genomes
  (hifiasm 67 GB on 48 cores, PAV 42 GB at 8 cores, dipcall 41.5 of 48 GB),
  the genome-wide split at 82 GB and 23 minutes per haplotype, the
  chromosome task at a quarter of its formula's request. The split's time
  per haplotype means the mapping has to be sharded before cohorts of
  hundreds of samples; a long-limit partition carries the smaller ones.

</details>

## 0.12.3

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.5.0 and
the image manifest unchanged.

A fix release for cohorts of hundreds of samples and for the pangenome's sex
chromosomes (0.12.2 was cut and superseded the same day, unpublished). Every
container now runs with Apptainer's exported bind list emptied: an
800-sample `cohort_call` had failed with "argument list too long". The
pangenome's cohort VCF carries every member on every chromosome, chrY
included. An existing install needs a four-word patch of its `miniwdl.cfg`
(chapter 09); a new site profile, `assembly`, carries the assembly track's
requests.

<details>
<summary>Details</summary>

- **The pangenome's cohort VCF has every member on every chromosome.**
  Cactus writes a sample column only for the members with sequence in a
  chromosome's graph, so the chrY VCF of the first real-sample cohort
  carried the XY members alone and the gather step's `bcftools concat`
  refused the mixed sample sets. `ugc_wgw_pangenome_vcf_post` now brings
  every chromosome VCF to the cohort's columns first: an absent member
  gets missing genotypes, every column in cohort order. A missing
  genotype keeps its one meaning, no assembled sequence there; the
  chromosome graphs replay from the cache on a retry.
- **Every container runs with Apptainer's exported bind list emptied.** The
  first 800-sample `cohort_call` failed in every `glnexus` shard before any
  tool ran, with "exec /bin/sh failed: argument list too long". Apptainer
  exports `APPTAINER_BIND` and `SINGULARITY_BIND` into the container with
  the whole bind list, miniwdl binds every input file (1,600 per shard),
  and past about 500 files that one environment string exceeds the kernel's
  128 KB limit. The rendered `miniwdl.cfg` now adds `--env APPTAINER_BIND=
  --env SINGULARITY_BIND=` to `[singularity] run_options`; nothing reads
  the variables. `ugc-wgw resources` warns when a SLURM install lacks the
  override, and chapter 09 says how to patch an existing install (four
  words in the installed `miniwdl.cfg`, then `ugc-wgw retry`). The next
  ceiling is the total argument size, about 4,000 samples.
- A third site profile, `assembly`, carries the assembly track's pilot
  requests from the first cluster run (PAV on 16 cores and 120 GiB, dipcall
  64 GiB, the pangenome's genome-wide split 192 GiB, ten runs in flight).
  Chapter 12 explains why the split's wall time, not its memory, grows with
  the cohort, the two routes for large cohorts (a long-limit partition, or
  sharding the mapping), and how memory above the standard node reaches a
  task (a policy row; the site cap rounds the workflow's own request down).

</details>

## 0.12.1

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.5.0 and
the image manifest unchanged.

A fix release for `cohort_pangenome` on clusters whose shared file system is
Ceph. The Cactus tasks keep Toil's directories and lock files on the
node-local scratch, since Toil refuses to keep them on Ceph. Without scratch
configured they stay in the task directory, which then has to be a local
file system.

<details>
<summary>Details</summary>

- **The Cactus tasks keep Toil's directories on the node-local scratch.**
  The first cluster run of `cohort_pangenome` failed in the genome-wide
  split: Toil, the engine inside Cactus, refuses to keep its lock files on
  Ceph ("file locks are known to break ceph filesystems"), and both Cactus
  tasks had put every Toil directory in the task directory on the shared
  file system. They now honour the `TMPDIR` the scratch backend sets
  (`TASK_SCRATCH_ENV`, chapter 04): with scratch, Toil's work directory, job
  stores, temporaries and lock files go to the node's disk and vanish with
  the job; without it they stay in the task directory, which then has to be
  a local file system. Both tasks' cache keys change; a graph is rebuilt per
  frozen cohort anyway.
- Docs: the assembly track's first cluster run is recorded
  (`tests/smoke/README.md`): `assembly` and `assembly_call` at the first
  attempt, the chr20 graph after the fix. Chapter 12 holds the first
  measured sizes of the two Cactus tasks (the split at 82 of 128 GB with
  eight haplotypes, the chromosome task at 4.5 of 64 GB), and the worked
  policy gains time rows for the assembly-track tasks: the catch-all's two
  hours would kill the split.

</details>

## 0.12.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.5.0,
unchanged.

The assembly track's reporting. The sample and cohort pages read every
assembly-track output: contiguity, Flagger's composition, both callers, the
pangenome's structure and the members' genotypes in it. New QC flags catch
an assembly that disagrees with its read-based SNVs (`asm_concordance_low`).
No expensive task is touched; every per-sample digest rebuilds once.

<details>
<summary>Details</summary>

- **The analysis summaries read every assembly-track output.** The sample
  page's Assembly section shows contiguity with the NGx curve and hifiasm's
  low-quality regions, Flagger's composition overall and per haplotype, a
  sequence QC table (GC, N bases, telomeric contig ends) with the contigs
  carrying the most non-Hap bases, trio QC with the contigs carrying the
  most switches, the chromosomes (assembled fraction per haplotype and a
  track of `reliable.bed` against `dip.bed`), the karyotype's evidence, both
  callers with Ts/Tv, het/hom, spectra, the PASS-versus-masked split and
  PAV's structural variants by size, the reliable region and the concordance
  in full. The cohort page draws the assembly metrics as distributions with
  outliers and threshold failures, lists one row per member, counts
  karyotypes, and for the pangenome shows the reference length and
  haplotype sequence per chromosome, the graph's structure (reference gaps,
  node degrees, nesting levels, vcfwave types and allele lengths), truvari's
  genotype concordance and the members' genotypes in the graph. Digest
  schema 2: every per-sample digest rebuilds once.
- New QC: `asm_concordance_low` fails a sample whose assembly disagrees with
  its read-based SNVs (`asm_genotype_concordance_min`, 0.90: a sample swap
  scores 0.5–0.7), `asm_karyotype_unknown` warns, and
  `pg_member_call_rate_min` (0.90) lists members with little assembled
  sequence in the graph on the cohort page.
- New small tasks on existing outputs, none of the expensive tasks touched:
  `assembly` writes a per-contig QC table (`contig_qc`: GC, N bases and
  runs, telomeric bases at each end, Flagger's bases per label; seqtk in the
  Flagger image), `assembly_call` the callers' spectra (`call_spectra`),
  `cohort_pangenome` the members' genotype statistics (`member_stats`) and,
  per chromosome, the join's path statistics and reference gaps, the bases
  per genome and the node degree histogram; its frequency summary counts
  records by nesting level, vcfwave type and allele length.
- The image manifest drops the retired `align_hifiasm` image (32 images in
  the bundle); the smoke harness drives `assembly_call` and `cohort_pangenome`
  as single stages; the masking task's pipeline passes shellcheck 0.9.0 (the
  0.11.0 tag's CI failure).

</details>

## 0.11.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.5.0,
unchanged. There was no 0.10.0.

The assembly track: the `assembly` mode becomes `assembly → assembly_call →
cohort_pangenome`, with assembly-based calling (dipcall, PAV) and a
Minigraph-Cactus pangenome per cohort. Breaking: the `assembly` stage ends
at a new output contract and several old outputs are gone. Expensive:
hifiasm's command changes, so every sample's assembly is recomputed once.
The tree at tag v0.11.0 fails its own CI on a shellcheck 0.9.0 finding in a
task command; public commit 50a73a0 on main carries the one-line fix.

<details>
<summary>Details</summary>

- **Breaking, `assembly`.** The stage ends at an output contract
  (`assembler_used`, `assembler_version`, `haplotype_labels`, `hap_fastas`,
  statistics, graphs); `haplotypes`, `zipped_assembly_fastas`, `asm_bams`,
  `merged_asm_bam` and the `paftools_*` outputs are gone; the paftools call is
  gone with them, while the `align_hifiasm` image, now unused, still ships in
  this release's manifest.
- **Expensive.** hifiasm runs with `--telo-m CCCTAA` (`hifiasm_telomere_motif`,
  empty = off) and writes its version file, so every sample's assembly is
  recomputed once on this version; no later release touches the hifiasm
  command. HMM-Flagger v1.2.0 (own tasks around the published images) and,
  for trio-binned samples, `yak trioeval` give per-assembly QC.
- New stage `assembly_call` (per sample): minimap2 2.31 aligns each
  haplotype to a decoy-free GRCh38 (BAM and PAF), the stage reports assembled
  bases per chromosome, infers the karyotype, and calls with dipcall v0.3 and
  PAV 3. Primary VCFs carry the sample's name, unphased genotypes unless
  trio-binned, and `FILTER=ASM_UNRELIABLE` outside the Flagger-reliable
  region (added, never dropped); with a read-based VCF a concordance table.
- New stage `cohort_pangenome` (per cohort): a Minigraph-Cactus v3.3.0 graph
  of the members' assemblies, built per chromosome behind one genome-wide
  split (`chromosomes` under `stage_inputs` limits the build); the cohort
  VCFs (raw and vcfwave) under the sample IDs with non-trio members
  unphased, allele frequencies with karyotype strata, GBZ and GFA per
  chromosome, optional `truvari bench` against the cohort's read-based SV
  merge (`pangenome_sv_cohort`, a profile key). The graph's paths of members
  that were not trio-binned are not haplotypes: no LD, imputation or
  PanGenie-style use without statistical phasing (guide chapter 02).
- References: the `ugc-wgw-asm` 0.1.0 lock entry (the decoy-free contig list,
  the chromosome list, dipcall's PAR BED; the installer derives the FASTA and
  index from the data container); the map gains `asm_fasta`,
  `asm_fasta_index`, `asm_par_bed`, `asm_graph_chromosomes`.
- Images: a new `minimap2` image of ours (minimap2 2.31, samtools 1.24,
  paftools.js) and the published dipcall, HMM-Flagger, secphase, PAV,
  Minigraph-Cactus and Truvari images pinned by digest.
- Summary: assembly QC thresholds and flags, an "Assembly QC" block per
  sample, QC columns in the cohort table, `summary --pass-list` for freezing
  a QC-gated pangenome cohort, and `summary --mode assembly` with the cohort
  page's pangenome section.

</details>

## 0.9.4

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.5.0,
unchanged.

The `parabricks` profile asks for 4 GPUs of any type and gives Parabricks
DeepVariant 256 GiB, after the GPU shapes were measured on five samples.
Site copies of the profile under `<prefix>/profiles/` are never overwritten
by the installer and need the same edit. The only workflow change is a
parameter description.

<details>
<summary>Details</summary>

- The `parabricks` profile asks for 4 GPUs of any type (`gpu_type` empty)
  and gives `run_parabricks_deepvariant` 256 GiB through the `mem_gb` stage
  input of `singleton` and `upstream`, after measuring the GPU shapes on
  five samples: four GPUs per job with an empty type use every GPU node,
  six GPUs take as long as four, three 2-GPU jobs give a 6-GPU node its
  best throughput, and the calls do not depend on the GPU count or model.
  Site copies of the profile under `<prefix>/profiles/` are never
  overwritten by the installer and need the same edit. Chapter 12 and the
  HPC checklist carry the numbers.
- The `gpuType` descriptions of `ugc_wgw_singleton` and `ugc_wgw_upstream`
  say the plugin applies the type at submit; they still described the
  behaviour before 0.9.3.

</details>

## 0.9.3

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.5.0. No
workflow changed since 0.9.1.

The GPU type is applied by the plugin at submit instead of being a workflow
input, so it is no longer part of every call-cache key and one project can
use nodes of different GPU types. Projects created before this release align
once more. `parabricks_gpus` must be even and is refused otherwise before a
job is submitted.

<details>
<summary>Details</summary>

- The GPU type is applied by the plugin at submit instead of being passed
  as the workflow's `gpuType` input. That input reached every task through
  `RuntimeAttributes`, so the type was part of every call-cache key and
  moving a project between GPU node types re-aligned every sample. The
  driver now hands the project's `gpu_type` to the plugin per run
  (`MINIWDL__UGC_WGW__GPU_TYPE`, miniwdl's override of `[ugc_wgw]
  gpu_type`); an empty type requests any GPU of the partition, so one
  project can use nodes of different GPU types; alignments are shared across
  flavours. The run manifest records `gpu_type`; `--dry-run` prints the
  environment line. Projects created before this release align once more.
- `parabricks_gpus` must be even: upstream's task runs `pbrun deepvariant
  --run-partition`, which refuses an odd GPU count; `init`, `config.json`
  and profiles refuse it before a job is submitted. Chapter 12 explains the
  shapes a 6-GPU node allows and how a project shapes its Parabricks jobs
  through stage inputs.

</details>

## 0.9.2

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.4.0,
unchanged since 0.9.0. No workflow changed since 0.9.1.

The run report's Provenance names the container backend, the node-local
scratch and the resource policy file. `tests/bundle.sh` runs without
Apptainer again, which is what had failed the CI of 0.9.0 and 0.9.1.

<details>
<summary>Details</summary>

- The run report's Provenance names the container backend, the node-local
  scratch (`scratch_env` and its mount, or `off`) and the resource policy
  file of the install's `miniwdl.cfg`, so a report says whether tasks had
  the node's local disk.
- `tests/bundle.sh` runs without Apptainer again: its check of the
  `ugc_slurm_singularity` backend resolves the class through miniwdl's
  plugin registry instead of constructing it, whose initialisation probes
  `apptainer --version`. This is what failed the CI of 0.9.0 and 0.9.1
  after every other step had passed.

</details>

## 0.9.1

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.4.0,
unchanged since 0.9.0.

The pbmm2 reference index is built once at install instead of once per input
BAM, where a campaign's first wave started hundreds of identical 5.8 GB
builds. The reference map gains `pbmm2_index`. Alignment inputs change, so a
sample replayed into a project on this version aligns once more.

<details>
<summary>Details</summary>

- The pbmm2 reference index is built once at install instead of once per
  input BAM. Upstream's alignment subworkflow and its fail-reads path build
  the same 5.8 GB index for every BAM; the call cache catches it only after
  the first build, so a campaign's first wave starts hundreds of identical
  builds. `install-bundle.sh` builds it once with the bundled pbmm2 image
  (about a minute, 12 GB of memory, verified against the md5 in
  `references.lock`; it says how to build it in a job when the login node
  cannot), the reference map gains `pbmm2_index`, and two derived files,
  `workflows/overrides/pbmm2.wdl` and `workflows/overrides/upstream.wdl`,
  pass it to every alignment in the `singleton` and `upstream` stages;
  without an index they behave as upstream. Alignment inputs change, so a
  sample replayed into a project on this version aligns once more.

</details>

## 0.9.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.4.0.

Node-local scratch for tasks: when the site key `TASK_SCRATCH_ENV` names the
node's per-job scratch variable, every task's container gets that directory
as `TMPDIR`, so sort spill files and GLnexus's database leave the shared
file system. A job without usable scratch runs as before and says so in its
SLURM log. The miniwdl plugin moves to 0.4.0.

<details>
<summary>Details</summary>

- Node-local scratch for tasks. The plugin registers a container backend,
  `ugc_slurm_singularity`, which the rendered `miniwdl.cfg` selects: it is
  miniwdl-slurm's `slurm_singularity`, and when the new site key
  `TASK_SCRATCH_ENV` names the node's per-job scratch variable (`TMPDIR`,
  `SNIC_TMP`) it binds that directory at `/mnt/scratch` inside every task's
  container and sets `TMPDIR` to it, so pbmm2's sort spill files, GLnexus's
  database and other TMPDIR-aware temporaries leave the shared file system.
  A job whose variable is unset, unwritable or tmpfs runs as before, and
  one `[ugc-wgw] scratch:` line in the job's `slurm-<jobid>.out` says which.
  The installer validates the key and prints the setting; `ugc-wgw
  resources` shows `scratch=` in its summary line. Until now `TMPDIR` was
  never set inside a container, so GLnexus's database had always been on
  the shared work directory; the guide is corrected.

</details>

## 0.8.2

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.3.0,
unchanged since 0.8.0.

A CUDA host-memory failure of Parabricks is classified as a resource
failure: it logs to stdout, exits 255 and leaves no SLURM out-of-memory
state. The worked resource policy gives `run_parabricks_deepvariant` 256G,
since 192 GiB failed on deeper samples.

<details>
<summary>Details</summary>

- Failure classification reads the failed task's `stdout.txt` as well
  (after the SLURM log and `stderr.txt`) and counts a CUDA host-memory
  failure (`cudaErrorMemoryAllocation`, `CUDA_ERROR_OUT_OF_MEMORY`) as kind
  `resource`: Parabricks logs to stdout, leaves stderr empty and exits 255
  without a SLURM out-of-memory state when a pinned allocation does not
  fit.
- The worked resource policy gives `run_parabricks_deepvariant` 256G: the
  declared 192 GiB held at 19–25x and failed on deeper samples. Chapters 09
  and 12 explain where Parabricks' log is and why SLURM shows nothing.

</details>

## 0.8.1

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.3.0,
unchanged since 0.8.0.

`samples add` and `samples check` show the progress of the input check. The
guide's examples and quick start are brought up to date.

<details>
<summary>Details</summary>

- `samples add` and `samples check` show the progress of the input check
  (`inspecting read files N/M (p%), elapsed, about T left`, rewritten in
  place on a terminal, one line every ten seconds otherwise) and end with
  `inspected N read file(s) in T`; `-v` logs each sample's yield.
- Guide: the examples show the current version and the installer's current
  report; the quick start describes the input check, `--profile cpu` and
  `ugc-wgw usage`; the command reference lists exit code 2.

</details>

## 0.8.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is 0.3.0.

`ugc-wgw samples add` checks every input BAM and refuses an empty, truncated
or read-less one; `samples check` lists size, reads and bases per file.
`submit` checks the raw read files again just before each run and fails the
attempt without a SLURM job when they do not hold. Outputs that cannot be
hardlinked across file systems become symlinks with a warning. The miniwdl
plugin moves to 0.3.0.

<details>
<summary>Details</summary>

- Input BAM check: `ugc-wgw samples add` opens every read file with a
  stdlib BGZF reader. An empty, truncated or read-less BAM is refused
  (`--drop-empty` registers the sample without its empty or read-less
  files; `--no-inspect` keeps only the existence check); aligned inputs and
  missing read groups warn; the read count and bases are estimated from
  the first records (exact from a `.pbi`) and recorded with the sample;
  `input_thresholds` (`file_reads_min`, `sample_gbases_min`; a config and
  profile key) warn about thin files and samples. New `ugc-wgw samples
  check [ID...] [--stored]` lists size, reads, bases, movie and problems
  per file.
- `submit` checks the raw read files just before each run that names them
  (present, same size as registered, BGZF end marker) and fails the attempt
  as an input error without a SLURM job when they do not hold; `--dry-run`
  prints the result.
- Outputs across file systems: a hardlink miniwdl cannot make across file
  systems (a cached output of a run whose results live elsewhere, an input
  passed through) becomes a symlink with a warning instead of
  `Invalid cross-device link`; the installer writes `ugc_wgw_miniwdl.pth`
  into the venv for it and warns when the prefix, references and call cache
  are on different file systems, as `ugc-wgw init` does for the results.

</details>

## 0.7.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is unchanged.

Site profiles: `ugc-wgw init --profile` applies a JSON profile of the
project keys that depend on the cluster and the sample set. Two profiles
ship, `cpu` and `parabricks`. The installer fills `<prefix>/profiles/` once
and never overwrites it.

<details>
<summary>Details</summary>

- Site profiles: `ugc-wgw init --profile NAME|FILE` applies a JSON profile
  of the project keys that depend on the cluster and the sample set
  (`stage_inputs`, the DeepVariant flavour and GPU choice, `max_inflight`,
  prices, QC thresholds) from `<prefix>/profiles/`, which the installer
  fills from `backends/hpc/profiles/` once and never overwrites, or from
  the shipped examples; flags win over the profile; `config.json` records
  the profile's name, checksum and the keys it set, `ugc-wgw stage-inputs`
  marks its overrides and the report names it.
- Shipped profiles `cpu` and `parabricks`: no alignment chunking, pbmm2 on
  24 threads so two alignments share a 48-core node, 20 runs in flight.

</details>

## 0.6.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is unchanged.

`ugc-wgw usage` and a Resource usage section in the run report: core-hours,
CPU efficiency, GPU-hours, memory, queue wait, disk and a cost estimate from
prices you supply. The driver reads `sacct` for every run's jobs when the
run finishes; `usage --collect` does it later for older runs.

<details>
<summary>Details</summary>

- `ugc-wgw usage` and a "Resource usage" section in the run report:
  core-hours (allocated and requested), CPU-hours used with the efficiency,
  GPU-hours, memory GB-hours, queue wait, jobs, disk usage (`--sizes`) and a
  cost estimate from prices you supply (`--price KEY=VALUE`, `config.json`
  `prices`; `--basis allocated|requested`), per stage, task and subject,
  with the mean per sample (clean and as run).
- The driver reads `sacct` for every run's jobs when the run finishes and
  keeps it as `accounting.json` next to the manifest (compute time without
  the queue wait, the wait itself, CPUs, `TotalCPU`, peak memory, GPUs,
  disk); `usage --collect` reads it later for older runs; attempts without
  accounting are estimated from the workflow log and labelled so.
- The smoke harness's `sacct` step also prints the driver's per-task table.

</details>

## 0.5.0

Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`) and its reference
data container `GRCh38_GIABv3`, as 0.4.0; the miniwdl plugin is unchanged.

`ugc-wgw stage-inputs` lists every input of every stage with its type,
default, who fills it and the project's override; `--nested` adds the inputs
of the tasks inside each entrypoint. Upstream's task thread counts are
settable as such nested keys. The guide gains a worked site resource policy.

<details>
<summary>Details</summary>

- `ugc-wgw stage-inputs`: every input of every stage with its type, default,
  description, who fills it (the driver, a `config.json` key, or nobody, so
  it is yours under `stage_inputs`) and the project's current override; it
  warns about overrides `submit` would refuse. `--nested` adds the
  call-qualified inputs of the tasks inside each entrypoint (thread counts,
  memory, tool options), listed by miniwdl itself, and checks dotted
  overrides against them.
- Upstream's task thread counts are settable as such nested keys
  (`upstream.pbmm2.pbmm2_align_wgs.threads`), which size both the SLURM
  request and the command. `ugc-wgw resources` says so when a policy row
  sets `cpu` for a task that interpolates its thread count, and names the
  entrypoint input where one exists (`hifiasm_threads`); the task inventory
  records those inputs.
- Guide: `examples/resources.tsv`, a worked site policy for 48-core, 384 GB
  nodes with its reasoning in chapter 12; every entrypoint input now carries
  a `parameter_meta` description.
- Smoke harness: `run.sh init --hifiasm-threads N`; the drivers it tees run
  with `--color never`.
- A non-object `stage_inputs.<stage>` in `config.json` is a clean error.

</details>

## 0.4.1

Patch release: the shell scripts pass shellcheck 0.9.0, which the CI runs.
No workflow or driver change; a 0.4.0 bundle stays valid.

## 0.4.0

First public release. Built on PacBio HiFi-human-WGS-WDL v4.0.0 (`15e82cb9`)
and its reference data container `GRCh38_GIABv3`.

Seven workflow entrypoints over PacBio's HiFi-human-WGS-WDL in standalone,
joint and assembly modes, with the cohort stages scattered by region at
every cohort size. The `ugc-wgw` driver plans, submits and retries runs and
writes provenance, reports and analysis summaries. Everything installs from
an offline bundle with digest-pinned containers.

<details>
<summary>Details</summary>

- Seven entrypoints: `singleton`, `upstream`, `cohort_call`, `downstream`,
  `cohort_merge`, `cohort_freq`, `assembly`; standalone, joint and assembly
  modes; cohort stages scattered over 26 region shards at every cohort size.
- The `ugc-wgw` driver: sample sheet and frozen cohorts, plan and submit with
  a concurrency cap, automatic retries of transient failures, orphan
  cancellation, SQLite state, per-run provenance manifests, `status`,
  `progress`, `logs`, `report` and `summary`.
- Offline bundle build and install with digest-pinned containers, engine
  wheels and the reference data container; site caps and a per-task resource
  policy applied by a miniwdl plugin; GPU DeepVariant and Parabricks as a
  project switch.
- Per-sample and per-cohort HTML analysis summaries with QC flags; a
  self-contained run report.
- `scripts/manifest-to-samples.py` turns a per-file delivery manifest into
  the sample sheet.

</details>
