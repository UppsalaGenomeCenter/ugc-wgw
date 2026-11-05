somalier — genotypes a BAM at ~17 000 informative sites (`extract`), then relatedness, sex from chrX
heterozygosity and heterozygosity per sample (`relate`) and a projection onto 1000 Genomes
superpopulations (`ancestry`). https://github.com/brentp/somalier (MIT).
Used by `ugc_wgw_upstream`/`ugc_wgw_singleton` (the extract and a one-sample relate for the sex
resolution, `workflows/ugc_wgw/qc/somalier.wdl`) and by `ugc_wgw_cohort_qc`. Pinned release binary
and sha256 in the Dockerfile; the sites file `sites.hg38.vcf.gz`, the 1000 Genomes set and its labels
ship as references (`references.lock`, `ugc-wgw-qc`).
