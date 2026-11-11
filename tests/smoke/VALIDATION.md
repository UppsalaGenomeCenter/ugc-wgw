# Truth-set scores of the smoke's HG002

One row per route and scoring run of `tests/smoke/run.sh validate`: the
smoke's HG002 (chr20, 48.7x) against the Genome in a Bottle HG002 v5.0q
benchmark for GRCh38, restricted to chr20 (`references.lock` entry
`ugc-wgw-truth`). Small variants: PASS records where the member carries an
alternate allele, genotype match inside the small-variant benchmark regions
with aardvark 1.1.0 (recall / precision / F1). Structural variants of 50 bp
and more: `truvari bench --passonly --pick ac` inside the SV benchmark
regions, then `truvari refine`, which GIAB recommends for this benchmark
because it keeps complex SVs whose representations differ between callers.
`singleton` is DeepVariant's own call set; `downstream` is the member's
slice of the GLnexus joint VCF.

The numbers are a record, not a gate: no floor is set until a few releases
have shown the run-to-run spread (`--min-f1` exists for then). They are for
chr20 of one genome at 48.7x; v5.0q is a harder benchmark than v4.2.1 (it
includes tandem repeats and complex regions), so the indel and SV figures
are lower than the ones quoted against older truth sets.

Run it before each tag and after each upstream sync, and add the rows.

| date | ugc-wgw | upstream | run | small-variant caller | SV caller | SNV R / P / F1 | indel R / P / F1 | SV bench R / P / F1 | SV refine R / P / F1 | SV GT concordance (bench) |
|---|---|---|---|---|---|---|---|---|---|---|
| 2026-11-11 | 0.16.0 | v4.0.0 | HG002 singleton (chr20) | DeepVariant 1.10.0 | sawfish 2.3.0 | 0.9983 / 0.9988 / 0.9985 | 0.9503 / 0.9495 / 0.9499 | 0.8092 / 0.8671 / 0.8372 | 0.9252 / 0.9314 / 0.9283 | 0.9077 |
| 2026-11-11 | 0.16.0 | v4.0.0 | HG002 downstream (chr20) | GLnexus v1.4.3-0-gcecf42e | sawfish 2.3.0 | 0.9964 / 0.9990 / 0.9977 | 0.9288 / 0.9497 / 0.9391 | 0.8092 / 0.8671 / 0.8372 | 0.9252 / 0.9314 / 0.9283 | 0.9077 |
| 2026-11-11 | 0.17.0 | v4.0.0 | HG002 singleton (chr20) | DeepVariant 1.10.0 | sawfish 2.3.0 | 0.9983 / 0.9988 / 0.9985 | 0.9503 / 0.9495 / 0.9499 | 0.8092 / 0.8671 / 0.8372 | 0.9252 / 0.9314 / 0.9283 | 0.9077 |
| 2026-11-11 | 0.17.0 | v4.0.0 | HG002 downstream (chr20) | GLnexus v1.4.3-0-gcecf42e | sawfish 2.3.0 | 0.9964 / 0.9990 / 0.9977 | 0.9288 / 0.9497 / 0.9391 | 0.8092 / 0.8671 / 0.8372 | 0.9252 / 0.9314 / 0.9283 | 0.9077 |

The smoke runs before a release's version bump, so a row carries the
version its runs were made at: the 0.16.0 rows are the smoke of the 0.17.0
release (sawfish 2.3.0 as the default), the 0.17.0 rows the smoke of
0.18.0. Identical scores across the two are expected: no caller changed
between them.
