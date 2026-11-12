# 14 Reading the pages

`ugc-wgw summary` writes one page per sample and one per cohort (chapter
07 has the command and where the files go). This chapter is for the person
who opens them: what each section shows, how to read it, and where the
explanation of an analysis lives. The pages are self-contained HTML:
inline styles and charts, nothing fetched from the network, so they can be
copied and mailed.

## A reading order for a new cohort

1. **Overview of the cohort page.** The cards (members, with results, QC
   pass, warn, fail), the stages that ran, and any warning for the cohort
   as a whole.
2. **The failures.** Type `fail` in the member finder on the left; chapter
   13 says what each flag means and what to do.
3. **Identity.** Duplicates, pedigree errors and contaminated samples
   before anything else is interpreted: they change who is in the cohort.
4. **Distributions.** One strip per metric across the members; the
   outliers are the members to open.
5. **The analyses**, in any order: frequencies, tandem repeats,
   methylation and cell-type composition, the targeted callers, and for
   an assembly cohort the assemblies and the pangenome.

## The cohort page

### Finding a member

The menu on the left holds the section links and, below them, every
member with a dot for its QC status (red fail, orange warn). The box above
the list narrows it as you type an id, a status or a flag name (`fail`,
`related_hub`, `gos_08`, and the name of a tool for members that were not
made with the default ones); Enter opens the first match.

### Overview and members

The cards, the cohort runs summarised (stage, run, version, attempt), the
members without a result, and the cohort's warnings: `mixed_tools` when
the members were made with more than one aligner or small-variant caller
(nothing is refused for it; joint calls, frequencies and outlier
statistics can show the difference as a batch effect). The members table
has one row per member with its QC status and flags and the main numbers;
click a header to sort, type in its box to filter, each id links to the
sample page.

### Distributions across members

Every metric on its own scale: the median, the interquartile box, the
members beyond the median ± 3 scaled MADs in red, QC thresholds as dashed
lines. Above 150 members the points become a band with only the outliers
drawn. The assembly metrics are part of it when the cohort has assemblies.

### Sex consistency

The sheet, the coverage inference, somalier's chrX call and the sex each
run used, side by side, with the disagreements listed. The assemblies'
karyotypes are counted next to them.

### Identity (`cohort_qc`)

The relatedness histogram of every pair with the thresholds drawn, the
pairs worth a look with what the sheet expects of them, the identity
flags, the ancestry projection (superpopulation counts and a PC1/PC2
scatter over the 1000 Genomes background) and the members whose extract
was made in the stage. Ancestry here is a description of structure, not
an origin. [The stage's page](stages/cohort_qc.md) explains the checks;
chapter 13 the flags.

### Cohort call sets and allele frequencies

The joint small-variant VCF, the merged SVs with the support histogram,
the tandem-repeat merge, then `cohort_freq`'s summary: allele-frequency
bins, singletons, call rate, sex strata.
[`cohort_freq`'s page](stages/cohort_freq.md) says how to read the counts
and use the files.

### Tandem repeats (`cohort_tr`)

The disease-locus screen's pathogenic and intermediate rows and any call
in doubt, each with the member's own dropout class at the locus; the
per-locus table; the expansion outliers per member; the loci with the
most listed members. Read the per-locus table first: it says where the
classes can be trusted
([`cohort_tr`'s page](stages/cohort_tr.md), "Reading the screen on a real
cohort"). A research screen, not a diagnostic result.

### Methylation and cell-type composition (`cohort_methyl`)

The stage's labels and their member counts, the events per member, the
genes with the most members carrying an event, the recurrent regions, the
imprinting screen and the comparisons; then the cell-type composition: a
ribbon with one column per member, the cohort's composition per cell
type, the members' spread, and the members far from the cohort. A member
with many events and an unusual composition is first of all a different
mixture of cells. [`cohort_methyl`'s page](stages/cohort_methyl.md) has
"Reading the scores" and "Cell-type composition".

### Targeted callers

Across the cohort: Paraphase copy-number histograms of the regions with
the most variation and a table of every region, StarPhase diplotype
frequencies per gene with the called rate, HLA allele frequencies,
mitorsaw's unexplained fraction with the members above its threshold, and
the members' haplogroups and macro-haplogroups. These are the tables
`ugc-wgw export` writes.

### Assemblies and the pangenome (assembly cohorts)

One row per member with both haplotypes side by side, the trio QC,
Flagger, the reliable fraction, the concordance and both callers' counts.
For `cohort_pangenome`: the chromosomes built with their nodes, edges and
length, the contigs the split assigned to a chromosome or left ambiguous
per haplotype, the allele frequencies of the vcfwave VCF, truvari's
recall, precision, F1 and genotype concordance against the read-based SV
merge, the graph's structure (haplotype sequence and reference gaps per
chromosome, the node degree histogram, records by nesting level, vcfwave
type and allele length) and the members' genotypes in the graph VCF (call
rate, singletons, het and hom-alt counts; members below
`pg_member_call_rate_min`, 0.90, are listed as a signal for the next
cohort freeze, not as a sample flag).

## The sample page

A fixed menu on the left and these sections:

| Section | What it shows |
|---|---|
| Overview | The QC flags, the sex from every source (sheet, coverage, somalier) beside the one the run used, the runs summarised, and the aligner and small-variant caller behind the sample (preview tools are named as such). |
| Reads | Counts, lengths, quality, mapped share. |
| Coverage | Mean depth per chromosome, the chrX and chrY ratios behind the sex inference, depth along the genome in 1 Mb bins with sawfish's copy-number segments. |
| Small variants | Counts, Ts/Tv, het/hom, substitution spectrum, indel lengths, site quality and depth histograms, runs of homozygosity. |
| Structural variants | By type and size. |
| Phasing | Blocks, phased share, block N50. |
| Tandem repeats | Coverage dropouts and the catalog's known disease loci with their genotypes; with `--cohort`, the member's rows of the cohort's screen and its strongest expansion outliers. |
| Methylation | CpG counts, MethBat's region classes, the imprinting regions; with `--cohort`, the member's cell-type composition beside the cohort's medians. |
| Targeted callers | Paraphase, mitorsaw with the collapsed chrM counts and a haplogroup per haplotype, StarPhase, kivvi. |
| Assembly | When the sample has one: below. |
| Identity | With `--cohort`: closest relative, relatedness and expectation, heterozygosity ratio, somalier's chrX call, ancestry, its identity flags. |
| Messages, files, software | The workflow messages, the files read, tool versions as pinned, engine versions, the reference build and its checksums, citations. |

The blocks marked "with `--cohort`" come from the cohort's runs and are
joined when the page is made, so they appear on pages written by `summary
--cohort` and say which cohort they are from.

### The Assembly section

Every file the assembly stages write: the haplotypes' contiguity (length,
contigs, NG50, LG50 and auNG against a 3.1 Gb genome, hifiasm's
low-quality regions, the NGx curve); Flagger's composition of the diploid
assembly (Err, Dup, Hap, Col, Unk as a stacked bar, with the bases and the
flagged blocks), also per haplotype, with GC, N bases, telomeric contig
ends and the contigs with the most non-Hap bases; for trio-binned samples
the switch and Hamming errors with the contigs that carry the most
switches; the chromosomes (assembled fraction per chromosome and
haplotype, a track of `reliable.bed` and `dip.bed` coverage in 1 Mb bins,
a per-chromosome table); the karyotype with its evidence; both callers'
statistics (records, SNVs, indels, Ts/Tv, het and hom-alt genotypes,
missing genotypes, the indel-length and substitution spectra, records by
FILTER class and type, PAV's PASS structural variants by size); the
reliable region in full; and the concordance with the sample's read-based
SNVs. A sample with only an assembly gets these sections, provenance and
the files table like any other.
