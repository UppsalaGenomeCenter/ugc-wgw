# minimap2 — long-read and assembly aligner, with samtools and k8 + paftools.js (lh3/minimap2, samtools/samtools,
# attractivechaos/k8). Assembly v2 (docs/DESIGN.md §13): the read-to-assembly alignment HMM-Flagger needs, the
# haplotype-to-reference alignment of assembly_call, and paftools.js liftover of Flagger's masks. Replaces
# quay.io/pacbio/align_hifiasm (minimap2 2.17, deleted upstream in 2024). See containers/README.md.
# Release assets pinned by sha256; samtools' tarball carries its own htslib. paftools.js is the misc/ script of the
# minimap2 release, run by the k8 JavaScript shell (static x86_64 binary from the k8 release). tar --no-same-owner:
# rootless podman cannot map the tarballs' owner ids.
FROM debian:bookworm-slim@sha256:88200866dfff7ea7f5cbcb6ec7c8a701889efe6fe859fe64d6990e4b07ea4171 AS build
ARG MINIMAP2_VERSION=2.31
ARG MINIMAP2_SHA256=c1351de6319c123369c2f4f37ba0ccf18c7ace47e2b1c0a35e30056b4a3bd9c9
ARG SAMTOOLS_VERSION=1.24
ARG SAMTOOLS_SHA256=89b2a440123eeaa400392ce1736e7d60ce9041843027d76819753c5a8246bfdd
ARG K8_VERSION=1.2
ARG K8_SHA256=a86b160a82f3233a21235d21170f3719a600dbd96bf1ec705a6eb57d770953c9
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl bzip2 gcc g++ make \
      zlib1g-dev libbz2-dev liblzma-dev libdeflate-dev \
 && rm -rf /var/lib/apt/lists/*
RUN set -eux; cd /tmp; \
    curl -fsSL -o minimap2.tar.bz2 "https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}.tar.bz2"; \
    echo "${MINIMAP2_SHA256}  minimap2.tar.bz2" | sha256sum -c -; \
    tar --no-same-owner -xjf minimap2.tar.bz2; \
    make -C "minimap2-${MINIMAP2_VERSION}" -j"$(nproc)"; \
    install -m 0755 "minimap2-${MINIMAP2_VERSION}/minimap2" /usr/local/bin/minimap2; \
    install -m 0755 "minimap2-${MINIMAP2_VERSION}/misc/paftools.js" /usr/local/bin/paftools.js; \
    rm -rf /tmp/*
RUN set -eux; cd /tmp; \
    curl -fsSL -o samtools.tar.bz2 "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2"; \
    echo "${SAMTOOLS_SHA256}  samtools.tar.bz2" | sha256sum -c -; \
    tar --no-same-owner -xjf samtools.tar.bz2; \
    cd "samtools-${SAMTOOLS_VERSION}"; \
    ./configure --prefix=/usr/local --without-curses --disable-libcurl; \
    make -j"$(nproc)"; make install; \
    rm -rf /tmp/*
RUN set -eux; cd /tmp; \
    curl -fsSL -o k8.tar.bz2 "https://github.com/attractivechaos/k8/releases/download/v${K8_VERSION}/k8-${K8_VERSION}.tar.bz2"; \
    echo "${K8_SHA256}  k8.tar.bz2" | sha256sum -c -; \
    tar --no-same-owner -xjf k8.tar.bz2; \
    install -m 0755 "k8-${K8_VERSION}/k8-x86_64-Linux" /usr/local/bin/k8; \
    rm -rf /tmp/*

# Runtime: minimap2 links libz and libm; samtools links libz, libbz2, liblzma and libdeflate; k8 is static
FROM debian:bookworm-slim@sha256:88200866dfff7ea7f5cbcb6ec7c8a701889efe6fe859fe64d6990e4b07ea4171
RUN apt-get update && apt-get install -y --no-install-recommends zlib1g libbz2-1.0 liblzma5 libdeflate0 \
 && rm -rf /var/lib/apt/lists/*
COPY --from=build /usr/local/bin/minimap2 /usr/local/bin/paftools.js /usr/local/bin/k8 /usr/local/bin/samtools /usr/local/bin/
RUN minimap2 --version && samtools --version | head -1 && paftools.js version
RUN useradd -m -u 1000 user
USER user
