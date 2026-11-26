# `cohort_pangenome`

A Minigraph-Cactus pangenome graph of the members' assemblies and its
variants.

| | |
|---|---|
| Subject | cohort |
| Modes | assembly |
| All stages | [08 Stages](../08-stages.md) |

The last stage of the assembly mode builds a Minigraph-Cactus graph of the
cohort's haplotype assemblies on the decoy-free GRCh38 and turns it into
cohort VCFs with allele frequencies. The graph is built one chromosome at a
time (Cactus's `--mgSplit` route): every member's two haplotypes are mapped
to a reference-only graph in a task of their own, one step bins the contigs
by chromosome, then every requested chromosome is constructed, mapped,
aligned and joined, one task for each of the four steps. The tools are
cactus v3.3.0 (minigraph
0.21, vg 1.76.1, vcfwave, panacus 0.5.2) from the
`comparative-genomics-toolkit/cactus` image, bcftools 1.23 from
`pb_wdl_base` and truvari 5.5.0.

```mermaid
flowchart TD
  asm[/"hap_fastas and trio of every member (assembly)"/]:::data --> names["Cactus-safe names SAMPLE.1 and SAMPLE.2;<br/>per-chromosome BEDs and memory from the reference index"]
  kar[/"inferred_karyotype of every member (assembly_call)"/]:::data --> freq
  ref[/"asm_fasta"/]:::data --> rg["cactus-minigraph --refOnly:<br/>the reference-only graph, once per reference"]
  names --> map
  rg --> map
  subgraph mem["per member"]
    map["cactus-graphmap --mgSplit:<br/>its two haplotypes against that graph"]
  end
  map --> split["alignments gathered; cactus-graphmap-split:<br/>one tar per chromosome (genome-wide)"]
  split --> cons
  names --> cons
  subgraph per["per requested chromosome"]
    cons["cactus-minigraph construct: the SV graph,<br/>from every member or the first construct_members"]
    cons --> cmap["cactus-graphmap: every haplotype mapped to it"]
    cmap --> align["cactus-align --pangenome"]
    align --> chrom["cactus-graphmap-join --vcf --vcfwave --gbz --gfa"]
    chrom --> grow["panacus: coverage histograms,<br/>growth in the cohort's order"]
  end
  grow --> gr[/"chrom_growth"/]:::result
  chrom --> post["concat; sample IDs restored; non-trio members unphased;<br/>header lines (the raw and the vcfwave VCF)"]
  post --> freq["ugc_wgw_bcftools_freq per chromosome on the vcfwave VCF,<br/>missing genotypes kept; concat; summary"]
  post --> ms["bcftools stats per chromosome: each member's genotypes, gathered"]
  chrom --> gs[/"path-stats, refgaps, genome-bp, degree per chromosome"/]:::result
  sv[/"cohort_merge SV VCF of pangenome_sv_cohort, else of this cohort, if any"/]:::data --> tv["truvari bench, 50 bp and longer, sequence comparison off"]
  post --> tv
  post --> out[/"pangenome_vcf, pangenome_wave_vcf, chrom_gbzs, chrom_gfas,<br/>freq_vcf, freq_summary, member_stats, truvari_summary"/]:::result
  freq --> out
  ms --> out
  tv --> out
  classDef default fill:none,stroke:#8a8f95
  classDef data fill:#5d6875,stroke:#5d6875,color:#ffffff
  classDef result fill:#2f7d4a,stroke:#2f7d4a,color:#ffffff
  classDef box fill:#1f7a70,fill-opacity:0.10,stroke:#1f7a70,stroke-opacity:0.45
  class mem,per box
```

Text version: `ugc_wgw_pangenome_names` maps each sample ID to a Cactus-safe
genome name (dots and dashes to `_`; a collision fails the run), writes the
rename table and, from the reference index, a BED and a memory size per
requested chromosome; `ugc_wgw_cactus_refgraph` runs `cactus-minigraph
--refOnly` on the reference; per member `ugc_wgw_cactus_map` runs
`cactus-graphmap --mgSplit` of its two haplotypes (`<name>.1`, `<name>.2`)
against that graph; `ugc_wgw_cactus_split` writes the seqfile (`GRCh38`,
then every haplotype in cohort order), concatenates the members'
alignments, runs `cactus-graphmap-split` over the whole genome, and tars
each chromosome's directory with its seqfile; per requested chromosome
four tasks follow one another, each untarring what it needs and rewriting
the paths to itself: `ugc_wgw_cactus_chrom_construct` runs
`cactus-minigraph` (the SV graph, from every haplotype of the chromosome
or from those of the first `construct_members` members),
`ugc_wgw_cactus_chrom_map` runs `cactus-graphmap` (every haplotype mapped
to that graph, several at a time), `ugc_wgw_cactus_chrom_align` runs
`cactus-align --pangenome` and `ugc_wgw_cactus_chrom_join` runs
`cactus-graphmap-join` with `--vcf --vcfwave --gbz --gfa`; then
`ugc_wgw_pangenome_growth` runs panacus on the chromosome's GFA; two
`ugc_wgw_pangenome_vcf_post` calls bring every chromosome VCF to the
cohort's sample columns (a member without sequence in a chromosome's
graph, an XX member on chrY, gets missing genotypes there), concatenate
the raw and the vcfwave VCFs, restore the sample IDs (`bcftools
reheader`), unphase the members that were not trio-binned and add the
header lines (cohort, chromosomes, the ambiguity rule, the phase policy);
the `cohort_freq` tasks run per
chromosome on the vcfwave VCF with missing genotypes kept, their shards are
concatenated and `ugc_wgw_pangenome_freq_summary` writes the summary (since
0.12.0 also records by nesting level, vcfwave type and allele length);
`ugc_wgw_pangenome_member_stats` runs `bcftools stats --samples -` per
chromosome on the same VCF and the gather task sums each member's
genotypes, singletons, missing calls and call rate (`member_stats`); the
join task keeps the join's path statistics and reference gaps and adds
the bases per genome and the node degree histogram; with a read-based SV
merge, a biallelic copy of the vcfwave VCF (`bcftools norm -m -any`; truvari
refuses multiallelic records) goes through `ugc_wgw_truvari_bench`; the
manifest.

Contigs enter a chromosome's graph only when enough of their length aligns
to one reference contig (75, 50 or 25 % below 100 kb, below 1 Mb and from 1
Mb, with a best-to-second ratio of 2 at MAPQ 5 or more, Cactus's defaults);
the rest is ambiguous and in no graph. `ambiguous_contigs` counts them per
haplotype, and the summary's cohort page shows the totals. The other
contigs (chrM and the unplaced ones, lumped together as `chrOther`) are
built only with `build_other_contigs = true`.

| Driver fills | You may override |
|---|---|
| `cohort_id`, `sample_ids`, `hap1_fastas`, `hap2_fastas` and `trios` (from every member's `assembly`), `sample_sexes` and `sample_sex_sources` (`XX`, `XY` or empty from every member's `assembly_call` karyotype; the sheet's sex is not used), `sv_vcf` and its index (the `cohort_merge` SV VCF of `pangenome_sv_cohort`, else of this cohort, any version and either read mode, when one exists) | `chromosomes` (the map's chr1–22, X, Y), `build_other_contigs` (false), `split_map_cores` (8: cores of each member's mapping task), `map_mem_gb` (60 GB plus 15 times the member's gzipped FASTAs), `split_cpu` (16), `split_mem_gb` (64), `split_disk_gb` (from the inputs), `construct_members` (unset: every member builds the SV graph), `construct_cpu` (4), `construct_mem_gb` (32), `chrom_cpu` (16: the map, align and join tasks), `chrom_map_cores` (4), `chrom_mem_gb_floor` (64), `chrom_mem_gb_per_100mb` (32), `chrom_mem_gb_per_gb` (4; the largest of the three rules wins), `freq_threads` (2), `freq_mem_gb` (4) |

`chromosomes = ["chr20"]` under `stage_inputs.cohort_pangenome` is how the
ladder runs (and the smoke): every haplotype is still mapped genome-wide
and the split covers every chromosome, only the per-chromosome tasks are
limited. Every
member needs a successful `assembly` and `assembly_call` at the version
(or `--any-version`); freeze the pangenome cohort from `ugc-wgw summary
--pass-list` (chapter 06) so that no assembly with a QC failure enters the
graph. `ugc-wgw cohort rank` (chapter 06) chooses which members and how
many: it orders them by the structural variants each adds and shows where
the cohort's common variation is covered.

**The ranking is recommended, not required.** The stage builds a graph of
any frozen cohort whose members are assembled:

```bash
ugc-wgw summary --mode assembly --pass-list asm-pass.txt     # the members whose assemblies pass QC
ugc-wgw cohort freeze G1 --samples asm-pass.txt              # or any list of yours
ugc-wgw submit --mode assembly --stage cohort_pangenome --cohort G1
```

The ranking needs the read-based track (the cohort's `cohort_merge` and
`cohort_qc`), and that is the route to prefer: it is also where a
contaminated sample and a pair of relatives are found, which nothing in
the assembly track notices. Without it the choice of members is yours,
and the only thing this stage reads from the read-based results, the SV
merge it compares the graph with, is skipped with a note when the cohort
has none. Chapter 01, "How the tracks relate", has the whole picture.

A graph of a subset has genotypes for its own members only.
[`cohort_graph_gt`](cohort_graph_gt.md) takes the graph's structural
variants and genotypes them in every member of a read-based cohort.

## One mapping per member, and the ladder

A haplotype's alignment to the reference-only graph depends on that
haplotype and the reference, not on the other members. So the mapping is
one task per member, with the member's files as its only inputs, and
miniwdl's call cache keeps it: a second cohort that holds the same member
does not map it again. That is what makes a ladder of nested cohorts
affordable. Freeze them as prefixes of one ranking (`ugc-wgw cohort rank
--emit`, chapter 06), set `chromosomes` to `["chr20"]`, and run them
smallest first:

```bash
for k in 20 100 200 400; do
  ugc-wgw cohort rank C1 --emit $k --pass-list asm-pass.txt
  ugc-wgw cohort freeze P$k --samples <results>/reports/rank/C1.top$k.txt
  ugc-wgw submit --mode assembly --stage cohort_pangenome --cohort P$k
done
```

P100 maps 80 members, P200 another hundred, and so on; every step runs the
split and the chr20 tasks again, which is what is being measured. After
each step `ugc-wgw usage --cohort P$k --by task` has the time and memory of
`ugc_wgw_cactus_split` and the four `ugc_wgw_cactus_chrom_*` tasks against
the number of haplotypes, and the cohort page's growth curves say what
the members added. The same cohorts built on every chromosome are then graphs of
those sizes; nothing is mapped again for them.

The split's result is the same as when one task mapped everything (checked
on the dev machine: contig assignments, per-chromosome FASTAs and
alignments). The two routes are not a choice: there is one.

## The SV graph from the first members

Of a chromosome's four steps one cannot be divided. minigraph builds the
SV graph, the skeleton of large variation that everything else is aligned
onto, by adding one haplotype after another, each mapped to the graph the
earlier ones left: about three minutes per haplotype of chr20, on one or
two cores whatever the task is given. With 200 haplotypes that is ten
hours of the chromosome's nineteen; with 1,700 it is most of a week.

`construct_members` limits that step to the cohort's first members:

```json
"stage_inputs": {"cohort_pangenome": {"chromosomes": ["chr20"], "construct_members": 100}}
```

The haplotypes of the first hundred members then build the SV graph, and
every member of the cohort is mapped to it, aligned and joined as before:
the graph, its VCFs and its growth curves hold everyone. "First" is the
cohort's own order, which is the ranking's when the cohort was frozen from
`ugc-wgw cohort rank --emit` (chapter 06): the members that carry the most
structural variation not yet covered come first, which is the order a
skeleton wants. Unset, every member builds it, as before 0.26.0. On a
chromosome where none of the first members has sequence (chrY and a
construct of XX members) every haplotype of that chromosome builds it.

What it costs, as far as it has been measured (the dev machine's four
chr20 assemblies, two of them in the construct, against the same graph
built with all four):

- **A member left out keeps its variants.** Its sequence is still aligned
  base by base between the anchors the SV graph gives: 99.1 % of its
  alleles of 50 bp or more are found again away from the two stretches
  below (the members inside the construct: 99.2 and 99.6 %), 97.5 % of
  those that no member of the construct carries, 28 of 29 of those of 1 kb
  or more, and 99.92 % of its small variants.
- **Satellite-rich stretches come out differently for everyone.** At the
  centromere's flank and the chromosome's end, 86 to 90 % of the alleles
  are found again, for members in and out of the construct alike, and
  some haplotypes have less of their sequence in the graph there (350 kb
  of 61 Mb). A smaller SV graph gives the aligner fewer anchors where
  alignment is hard anyway.
- **Not measured:** insertions above 10 kb that only a member outside the
  construct carries. Cactus clips sequence that stays unaligned over more
  than 10 kb, so that class is the one to look at first on real data.

That is one family and one chromosome, so the cohort page states what was
done ("SV graph from the first N of M members"), the VCF header carries
`##ugc_wgw_cohort_pangenome_construct`, and `chrom_construct_genomes`
lists the haplotypes of every chromosome's construct. To judge it on your
own data, build one rung twice: with and without `construct_members`, the
second run reusing the mappings and the split.

## Growth curves

`chrom_growth` holds, per chromosome, panacus's coverage histograms: how
many base pairs and nodes of the graph are carried by exactly one member,
by two, and so on, once for everything the members carry and once for the
sequence that is not on the reference's path. From them the cohort page
draws what k members drawn at random carry between them (the growth
curve), what all of them share (the core), and, beside the expectation,
the cumulative amount in the cohort's own order. The chromosomes add up
because each is counted against the whole cohort: an XX member carries
nothing on chrY. Read the non-reference curve for where a further member
stops bringing new sequence. It never goes quite flat: satellite arrays
and other sequence that aligns poorly look new in every genome.

## Resources

Dev-machine measurements on chr20-only assemblies: one member's mapping
task about an hour for its two haplotypes at 55 GB (the reference graph's
index; a whole-genome haplotype took 23 minutes and 65 to 80 GB on the
cluster), the reference-only graph 4 minutes, the split 17 minutes, one
chromosome task about 50 minutes with `cactus-align` peaking at 11 GB; the
growth task, the VCF post-processing, the frequencies and truvari are
minutes.

Each of a chromosome's four tasks has its own request, and the call cache
is a checkpoint between them: a task that fails or runs out of time is
run again alone, with the earlier ones replayed.

| task | cores | memory | what its time follows |
|---|---|---|---|
| `ugc_wgw_cactus_chrom_construct` | `construct_cpu` (4) | `construct_mem_gb` (32) | haplotypes in the construct, one after another |
| `ugc_wgw_cactus_chrom_map` | `chrom_cpu` (16) | 16 GB plus 10 per mapping at once | haplotypes, `chrom_cpu / chrom_map_cores` at a time |
| `ugc_wgw_cactus_chrom_align` | `chrom_cpu` | the rule below | haplotypes, faster than in proportion |
| `ugc_wgw_cactus_chrom_join` | `chrom_cpu` | the rule below | the graph's most complex region (vcfwave) |

The memory of the align and join tasks is the largest of three rules:
`chrom_mem_gb_floor` (64), `chrom_mem_gb_per_100mb` (32) times the
chromosome's length in 100 Mb, and 32 GB plus `chrom_mem_gb_per_gb` (4)
per Gb of haplotype sequence in the chromosome's graph (the haplotypes
and the reference, times its length). The third knows the cohort: chr20
gets 64 GB up to 40 members, 84 at 100 and 136 at 200; the align step
peaked at 33 and 55 GiB at 20 and 100 members (chapter 15). None of the
four has a wall-time default: give them rows in `resources.tsv` (chapter
12; one row `ugc_wgw_cactus_chrom_*` covers all four).

Outputs: the graph's cohort VCF and its vcfwave form, the GBZ and GFA per
chromosome with their statistics (graph size, path fragments, reference
gaps, bases per genome, node degrees, growth) and logs, the name map and the
split's tables, the frequency VCF with its summary and samples tables, the
members' genotype statistics, truvari's figures (chapter 07).
