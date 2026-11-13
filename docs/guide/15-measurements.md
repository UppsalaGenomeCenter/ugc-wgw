# 15 Measurements and evaluations

What the pipeline's tasks used when they ran, and the comparisons that
decided which tools it runs. Chapter 12 explains how requests are set and
capped; this chapter is the evidence to set them from. Every figure names
where it was measured: the chr20 smoke data on the development machine,
or whole genomes on the cluster.

## The pblab tools on 50 production genomes

Measured on the cluster (2026-11-12; 50 genomes at about 20x, joint mode,
A100 and L40S nodes): `ugc_wgw_pbweave_align` 4 min 15 s on its 24 cores
with a peak of 58.8 GB (96 asked), 1.7 core-hours per sample where pbmm2
took 46; `ugc_wgw_pbdv_call_gpu` 5 min 52 s on 16 cores and one GPU with
31.8 GB (64 asked), 1.6 core-hours and 0.1 GPU-hours per sample where
Parabricks took 11 and 0.9, and a median queue wait of 12 minutes. The
whole joint run cost 20.5 core-hours per sample against about 83 before;
HiPhase (8.4 core-hours per sample) is now the largest task. The defaults
fit; `mem_gb` of the pbdv task could go to 48 to pack six jobs on a
64-core, six-GPU node (with `threads` 10).

## The cohort VCF's concatenation and split

On the 800-sample `cohort_call` of 2026-10-26, `ugc_wgw_bcftools_concat`
(the 134 GLnexus shards into the cohort VCF) and `ugc_wgw_bcftools_split`
(the cohort VCF into 800 single-sample VCFs) took about four hours each on
four threads, one after the other, between the last GLnexus shard and the
first `downstream` run. Neither tool used the threads for the work that
cost: `bcftools concat` parsed every record, about a microsecond per
genotype, and only recompressed in parallel; `bcftools +split` parsed every
record in one thread and wrote one single-threaded bgzip stream per sample,
having no thread option at all. Since 0.12.6 the concat task joins
same-source shards with `bcftools concat --naive` (compressed blocks copied,
sample columns checked), minutes instead of hours, and the split task runs
`+split` once per contig that holds records, `split_threads` (default 8) at
a time, with its pieces on the node-local scratch, then one concatenation
per sample in contig order. On the 864-sample cohort (2026-10-30) the naive
concat took 22 minutes and the split 2h14m at 8 threads, so `split_threads`
now defaults to 16 and `split_mem_gb` to 8 + 0.005 GB per sample per thread
(every `+split` process holds one output stream per sample): about an hour
for 864 samples, 77 GB. The records are byte-identical to the one-pass
commands' (verified on the joint smoke's GLnexus VCF in both filter modes).

The same run's `cohort_merge` spent 13h52m in one `trgt-lps` task over the
whole cohort TRGT VCF (and its first attempt died at the catch-all two
hours). Since 0.13.1 trgt-lps runs once per TRGT merge shard, `trgt_lps_threads`
(4) each, and a gather task joins the tables in the cohort VCF's record
order; with 26 shards the step is about half an hour.

The cohort tasks of that run, for the policy's rows: GLnexus 54 minutes
median and 89.7 of 119 GB at 864 samples (the formula 32 + 0.1 GB per
sample held); the gVCF slices 10 minutes each (864 small jobs, queue-bound);
svx 34 minutes per shard at 15.3 of 52 GB; the TRGT merge 20 minutes per
shard. The campaign's clock was the GPU queue: Parabricks ran 14 minutes
per sample and waited 48h39m in the median.

## The tandem-repeat stage

`cohort_tr` scatters one `ugc_wgw_tr_locus_stats` per shard of
`scatter_regions` (2 cpu, `stats_mem_gb` 2 GB): `bcftools query` streams
the shard's records of the merged TRGT VCF, and the script holds one line
at a time (a line is about 200 bytes per member, 10 kB and more at
hundreds of members), so the task's time is bcftools's reading time and
its memory flat. `ugc_wgw_tr_screen` (1 cpu, `screen_mem_gb` 2 GB) reads
135 records and one pass of the trgt-lps table; `ugc_wgw_tr_gather` (2 cpu,
`gather_mem_gb` 4 GB) concatenates and bgzips, holding the outlier pairs in
memory. On the chr20 smoke cohort the whole stage is under two minutes.
Measured on 864 whole genomes (2026-11-11): 26 `ugc_wgw_tr_locus_stats`
jobs of 32 s on average with a peak of 444 MB, 0.46 core-hours together;
the screen 26 s; the gather 32 s with 57 MB and 590 MB written. The
defaults are generous by a factor of four and can stay.

## The methylation layer and the intake stage

`cohort_methyl`'s tasks are small: `ugc_wgw_methbat_build` reads every
member's profile (about 3 MB each; `build_mem_gb`, 8 GB, covers thousands),
each `ugc_wgw_methbat_profile_bg` reads one 5mC pileup like upstream's
`methbat_profile` (1 cpu, `profile_mem_gb` 4 GB), `ugc_wgw_methbat_compare`
and the summary are seconds. The one task whose size grows with the cohort
is `ugc_wgw_methbat_joint_segment` (off by default): it loads every
member's pileup sites into memory (343 MB for four chr20 pileups;
`joint_segment_threads` 8, `joint_segment_mem_gb` 32 are first guesses
until a whole-genome cohort has measured it). The per-sample
`ugc_wgw_methbat_report` is 1 cpu, 4 GB, seconds. `intake` runs mosdepth
(upstream's request), the somalier extract and the checks (1 cpu, 2 GB,
minutes on a whole-genome BAM: the checks read the header and the first
2,000 reads, the VCF headers and sample lists).

## The pangenome's sizes

`cohort_pangenome` has two tasks whose rows the chr20 scaling pilot sets
(design brief, WP4): `ugc_wgw_cactus_split` maps every haplotype of the
cohort genome-wide in one job (`split_cpu`, `split_mem_gb`, `split_disk_gb`
as workflow inputs), and `ugc_wgw_cactus_chrom` builds one chromosome per
job (`chrom_cpu`; memory the larger of `chrom_mem_gb_floor` and
`chrom_mem_gb_per_100mb` per 100 Mb of the chromosome). The measured sizes
per cohort size go here as the pilots run; the first row is the cluster
smoke of 2026-10-22 (four samples, eight chr20-only haplotypes mapped
against the whole decoy-free GRCh38):

| task | haplotypes | cores | request | peak RSS | wall time |
|---|---|---|---|---|---|
| `ugc_wgw_cactus_split` | 8 (chr20 only) | 32 | 128 GB | 82.2 GB (64 %) | 57 min |
| `ugc_wgw_cactus_chrom` chr20 | 8 (chr20 only) | 16 | 64 GB (the floor) | 4.5 GB (7 %) | 58 min |
| `ugc_wgw_cactus_split` | 10 (whole genomes, 2026-10-26) | 32 | 128 GB | 81.9 GB (64 %) | 3h49m, 23 min per haplotype |
| `ugc_wgw_cactus_chrom`, 24 chromosomes | 10 (whole genomes) | 16 | 64–80 GB (the formula) | 19.9 GB (25 % of chr1's) | 56 min median; chr21 5h51m, chr1 1h46m, chr14 33 min |

The split's peak is the same with eight chr20 haplotypes and ten whole
genomes: its memory is minigraph's index of the reference graph plus Toil,
not the haplotype sequence, and 128 GB holds one mapping. Its wall time was
the haplotype count times 23 minutes at 32 cores, because both Cactus
tasks mapped one haplotype at a time (`--mapCores` equal to the cores) and
minigraph uses only a few cores of one haplotype's mapping: the split ran
at 22 % CPU efficiency, the chromosome task at 15 %. The chromosome task's
time does not follow the chromosome's length either: chr21, 47 Mb, took
5h51m where chr1 took 1h46m and chr14 33 minutes, and its logs say why:
4h43m of it was minigraph mapping the ten haplotypes' chr21 contigs one
after another, 73 and 58 minutes for the two slowest (satellite and rDNA
arrays of the acrocentric short arm; the current Minigraph-Cactus
documentation ignores softmasking and advises against hardmasking, so
masking is no remedy), half a minute for the reference.

Since 0.12.5 the two tasks take the cores per mapping as inputs
(`split_map_cores`, default the task's cores; `chrom_map_cores`, default 4)
and Toil runs `cores / map_cores` mappings at once, within the memory: a
chromosome's mappings are small, so the default already maps four haplotypes
at once on 16 cores (chr21 of five samples: an expected 2h20m instead of
5h51m, the construct step then the larger part); the split's mappings each
hold the reference graph's index, 65–80 GB, so the `assembly` profile gives
it a whole standard node (48 cores, 360 GiB) with `split_map_cores` 16,
three at once. The pilot builds every chromosome at its first point rather
than extrapolating from chr20, and measures how both tasks grow with the
haplotype count. Give both tasks a day of wall time (the rows of the worked
policy): the catch-all two hours would kill the split.

The per-sample tasks on five 20x genomes (2026-10-26, the `cpu` profile):

| task | cores | request | peak RSS | wall time |
|---|---|---|---|---|
| `ugc_wgw_hifiasm_assemble` | 48 | 288 GiB (a whole node) | 67 GB | 1h39m median |
| `ugc_wgw_flagger_align` | 32 | 64 GiB | 47 GB (73 %) | 1h37m |
| `ugc_wgw_minimap2_align_asm` | 16 | 32 GiB | 22 GB | 17 min per haplotype |
| `ugc_wgw_dipcall` | 8 | 48 GiB | 41.5 GB (87 %) | 47 min |
| `ugc_wgw_pav` | 8 | 40 GiB, then 80 | 41.7 GB | 3h19m; 40 GiB failed |
| `ugc_wgw_chrom_coverage` | 2 | 8 GiB | 93 MB | 1h31m before 0.12.4 |

The `assembly` profile's PAV (16 cores, 120 GiB) and dipcall (64 GiB)
requests rest on these peaks.

The `assembly` profile (chapter 05) carries the pilots' requests: PAV on 16
cores and 120 GiB (its two alignment jobs hold the whole reference's index,
about 16 GiB each, whatever the assembly's size, and the calling steps grow
with the variant count), dipcall 64 GiB, the split 192 GiB. Both the split
and the chromosome task read their sizes from workflow inputs
(`split_cpu`, `split_mem_gb`, `chrom_cpu`, `chrom_mem_gb_floor`,
`chrom_mem_gb_per_100mb`), not from policy rows.

What grows with the cohort is the split's wall time, not its memory: the
mappings run one after another, each with every core of the task, and on
the smoke each chr20-only haplotype took about six minutes at 32 threads,
most of it indexing the 3.1 Gb reference graph anew. A whole-genome
haplotype took 23 minutes on the first whole-genome cohort (five samples,
2026-10-26, one mapping at a time): 200 haplotypes were about three days
in one job and 1,600 about 26 days. Three mappings at once on a whole
standard node (the `assembly` profile) divide that by about three, and a
fat node with 2 TB holds twenty mappings at once, which brings 1,600
haplotypes to about two days. Beyond that, two routes exist. A
partition with a long limit and large nodes (on the smoke cluster `fat`:
96 cores, 2 or 4 TB, 14 days) takes the serial job as it is, through a row
`ugc_wgw_cactus_split` with that partition and limit plus `split_cpu` 96
in the profile; the whole mapping restarts if the node fails or the limit
hits, since miniwdl reruns the task; it carries cohorts of up to a few
hundred samples. Sharding the mapping into one job per haplotype on
standard nodes turns weeks of serial work into a day of parallel work and
is what the full cohort needs; it is a workflow change, planned.

The same nodes are the answer for a chromosome whose graph outgrows a
standard node: with 2 and 4 TB the design's gate (chr1 at 1,600 haplotypes
must fit the largest node) is about time, not memory. Two mechanics
matter there. `TASK_MEMORY_MAX` rounds a workflow's own request down to
the standard node, so a `chrom_mem_gb_per_100mb` that asks above it is
silently capped; a policy row is not capped, so memory above 375 GiB goes
through a row `ugc_wgw_cactus_chrom` that also names the partition. And a
row matches the task's name, not the chromosome, so it moves every
chromosome's job to that partition, whole nodes each; a smaller cohort
graph stays on standard nodes with the formula.

## The reads on their own assembly: minimap2 and pbweave

The `assembly` stage aligns a sample's reads back to its own two
haplotypes with minimap2 (an hour and a half per 20x genome on 32
cores). pbweave was tried in its place before 0.21.0, on one smoke
member's chr20 reads against its 137 Mb diploid assembly:

| | minimap2 2.31 | pbweave (pblab 20261021) |
|---|---|---|
| Alignment | 5 min 48 s | 10 s (index 1 s) |
| BAM, primary and supplementary | 1.32 GB | 1.42 GB |
| Reads within 1 kb of the other's position | 95.1% | |
| Reads with mapping quality 0, on hap1 / hap2 | 8,922 / 8,704 | 17,361 / 204 |
| HMM-Flagger: Hap / Err / Dup / Col | 98.63 / 0.01 / 1.21 / 0.15% | 96.54 / 0.94 / 2.46 / 0.06% |

Where the two haplotypes are identical a read fits both equally.
minimap2 picks one at random, so each keeps half the coverage; pbweave
picks the first sequence of the reference, always. Flagger then sees one
haplotype without reads and the other with too many, and a methylation
pileup would give one haplotype every read of such a stretch. That is
why this step has no aligner setting. (pbweave also writes none of the
`MD`, `cs` and `de` tags secphase needs.) Against GRCh38, where there is
one copy of everything, the same aligner agrees with pbmm2 (the section
on the pblab tools above).

The same test fixed three other things of the stage. The reads' `MM` and
`ML` tags reach the alignment through `samtools fastq -T` and
`minimap2 -y` with every placement unchanged. minimap2's secondary
alignments, one per read on the other haplotype, doubled the file and
made a pileup count each read twice; they are not written. And the
pileup's mapping-quality floor is 0: with 1, 4.4% of the CpGs of the
assembly had no row.

## Evaluation profiles

An evaluation profile runs the production entrypoints with one tool
swapped, so that a second project over the same members can be compared
with a production one before anything changes for everyone. The first is
`eval-sawfish23` (`backends/hpc/profiles/eval-sawfish23.json`): it set the
`sawfish_image` input of the sawfish tasks (`workflows/overrides/sawfish.wdl`,
upstream's task file with the image as an input) to PacBio's published
2.3.0_build1 for `singleton`, `upstream` and `cohort_call` while the
default was upstream's 2.2.1_build2; every other task is unchanged, the
call cache replays everything up to sawfish, and sawfish onward
recomputes. The bundle ships both images. Since 0.17.0 the default is
2.3.0_build1 and the profile changes nothing (kept as the record of how
the evaluation was run). On the dev machine `tests/smoke/run.sh eval
sawfish23` initialises `project-sawfish23/` with the profile, runs the
standalone phase, then compares each member's phased SV VCF with the
production project's with truvari (`--sizemin 50 --pctseq 0`: recall,
precision, F1, genotype concordance) and prints the sawfish tasks' wall
time, core-hours and peak memory from `ugc-wgw usage` for both. On the
cluster, `ugc-wgw init --profile eval-sawfish23` on a project of its own
does the same at whole-genome scale over any two or three members that
already have production results (no trio is needed: the measure is each
member's agreement with its own production call set, plus sawfish's wall
time and peak memory); in joint mode run `--stage upstream` only and
compare `upstream`'s `sv_vcf` files, since the eval project needs no joint
call. The numbers go into `tests/smoke/README.md` and the decision into
the design log. sawfish 2.3.0 promises `discover` at
16.5 minutes and 4 GiB instead of 48.9 minutes and 72 GiB per sample with
"usually identical" output: the override's default moves to 2.3 when the
SV sets match within run-to-run noise and the memory drop holds on our
data; otherwise it waits for the upstream sync. On the chr20 smoke data
(2026-11-06) the SV sets were identical for all four members (truvari F1
and genotype concordance 1.0) and `discover` ran five to seven times
faster; on three production genomes (2026-11-09) the call sets were again
identical, record for record, with `discover` at 6 min 25 s and 3.6 GB
against 28 min and 68 GB for 2.2.1 over the 864, so 2.3.0 became the
default in 0.17.0. A `resources.tsv` row of 16 GB for `sawfish_discover`
is safe with it (the task's own request is 128 GB). An evaluation profile
is never applied to a production project.

The same harness command compares the routes a project can choose:
`run.sh eval pblab` (pbweave and pbdv), `eval pbdv` (the caller alone, on
the production alignments), `eval pbweave`, and the `-cpu` variants of the
first two. Each makes its own project, runs the standalone phase and
`validate`, and then prints, per member against the production project,
the SV sets (truvari) and the PASS small variants shared and not, inside
and outside the benchmark regions, TRGT's allele lengths, the digests side
by side (`tests/smoke/compare_runs.py`, which also takes two digest
directories brought back from a cluster) and the changed tasks' resources.
The first results are in chapter 12 ("The pblab tools"), in the section
on 50 production genomes above, and in the smoke README.

A comparison between two versions of a tool says whether something
changed, not whether the calls are right. For that the smoke carries a
truth-set score: `tests/smoke/run.sh validate` compares the smoke's HG002
calls on chr20 with GIAB's HG002 v5.0q benchmark, small variants with
aardvark and structural variants with truvari (bench, then refine), for
the standalone route and the joint one, and prints the rows that
`tests/smoke/VALIDATION.md` keeps per release. On the first run (0.17.0,
DeepVariant 1.10.0, sawfish 2.3.0) the SNV F1 was 0.9985 for DeepVariant's
own set and 0.9977 for the GLnexus slice, the indel F1 0.950 and 0.939, the
SV F1 0.928 after refinement. It reports and records; `--min-f1` turns a
figure into a gate once a few releases have shown its spread.
