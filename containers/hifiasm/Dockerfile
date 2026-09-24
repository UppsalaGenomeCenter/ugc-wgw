# hifiasm — haplotype-resolved de novo assembler for PacBio HiFi reads (chhylp123/hifiasm)
# Built from the release tarball so the assembly image is ours to rebuild (docs/DESIGN.md §13);
# replaces quay.io/pacbio/hifiasm 0.25.0_build1. See containers/README.md.
# 0.25.0 (2025-03-18) is the newest release and the tip of master at build time; the tag has no "v" prefix.
FROM debian:bookworm-slim@sha256:88200866dfff7ea7f5cbcb6ec7c8a701889efe6fe859fe64d6990e4b07ea4171 AS build
ARG HIFIASM_VERSION=0.25.0
ARG HIFIASM_SHA256=51633138865207a9d41630da9377d46e4921ad4fc5facaa1740ceccae8611f1f
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl g++ make zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*
RUN set -eux; cd /tmp; \
    curl -fsSL -o "hifiasm-${HIFIASM_VERSION}.tar.gz" "https://github.com/chhylp123/hifiasm/archive/refs/tags/${HIFIASM_VERSION}.tar.gz"; \
    echo "${HIFIASM_SHA256}  hifiasm-${HIFIASM_VERSION}.tar.gz" | sha256sum -c -; \
    tar -xzf "hifiasm-${HIFIASM_VERSION}.tar.gz"; \
    make -C "hifiasm-${HIFIASM_VERSION}" -j"$(nproc)"; \
    install -m 0755 "hifiasm-${HIFIASM_VERSION}/hifiasm" /usr/local/bin/hifiasm; \
    rm -rf /tmp/*

# Runtime: the binary links only libz, libstdc++, libgcc, libm and libpthread, all in the slim base
FROM debian:bookworm-slim@sha256:88200866dfff7ea7f5cbcb6ec7c8a701889efe6fe859fe64d6990e4b07ea4171
COPY --from=build /usr/local/bin/hifiasm /usr/local/bin/hifiasm
RUN hifiasm --version
RUN useradd -m -u 1000 user
USER user
