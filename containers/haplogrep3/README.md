haplogrep3 — assigns a mitochondrial haplogroup to a variant profile (a VCF column) against a
phylogenetic tree. https://github.com/genepi/haplogrep3 (MIT). Used by `ugc_wgw_haplogrep`
(`workflows/ugc_wgw/mito/collapse.wdl`): one column per mitorsaw haplotype of a sample, classified
in `ugc_wgw_upstream`/`ugc_wgw_singleton` after mitorsaw, or inside `ugc_wgw_cohort_merge` for a
member whose run predates the task. A Java tool, so the image is eclipse-temurin's JRE 21 by digest
with the release tarball pinned by sha256. Two trees are installed at build time, because haplogrep3
fetches a tree from GitHub the first time it is named and the HPC has no network:
`phylotree-fu-rcrs@1.3` (PhyloTree 17 forensic update, the tasks' default) and `phylotree-rcrs@17.3`
(PhyloTree 17); `haplogrep_tree` is a workflow input. `/usr/local/bin/haplogrep3` is our wrapper
(quoted arguments, `HAPLOGREP3_XMX` heap, default 1500m); `haplogrep3 trees` lists what is installed.
