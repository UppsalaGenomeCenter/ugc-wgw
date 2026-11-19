# kanpig — long-read genotyper of sequence-resolved structural variants (ACEnglish/kanpig, MIT): `kanpig plup`
# turns an alignment into a pileup index, `kanpig gt` genotypes a VCF's SVs from the alignment or from that index.
# One static (musl) binary per release. The base is upstream's own pb_wdl_base, by digest, for the bgzip and tabix
# the pileup index needs (kanpig writes it unsorted and uncompressed) and the bcftools the genotype VCF is sorted
# with; the published quay.io/biocontainers/kanpig image holds the binary alone. See containers/README.md and
# docs/DESIGN.md §18 (2026-11-19), docs/expansion/pangenome-genotyping.design.md.
FROM debian:bookworm-slim@sha256:88200866dfff7ea7f5cbcb6ec7c8a701889efe6fe859fe64d6990e4b07ea4171 AS fetch
ARG KANPIG_VERSION=2.0.2
ARG KANPIG_SHA256=48705b7c61c6cc13043081c3c0350e63cc0a62cd8126bd40de4348a43615528d
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl \
 && rm -rf /var/lib/apt/lists/*
RUN set -eux; cd /tmp; \
    curl -fsSLo kanpig.tar.gz "https://github.com/ACEnglish/kanpig/releases/download/v${KANPIG_VERSION}/kanpig-v${KANPIG_VERSION}-x86_64-unknown-linux-musl.tar.gz"; \
    echo "${KANPIG_SHA256}  kanpig.tar.gz" | sha256sum -c -; \
    tar -xzf kanpig.tar.gz --strip-components=1; \
    install -m 0755 kanpig /usr/local/bin/kanpig

FROM quay.io/pacbio/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae
ARG KANPIG_VERSION=2.0.2
COPY --from=fetch /usr/local/bin/kanpig /usr/local/bin/kanpig
RUN set -eux; kanpig --version | grep -qx "kanpig ${KANPIG_VERSION}"; \
    bgzip --version | head -1; tabix --version | head -1; bcftools --version | head -1
