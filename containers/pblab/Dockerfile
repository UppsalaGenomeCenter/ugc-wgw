# pblab — PacBio's preview HiFi tools (PacificBiosciences/pblab, BSD-3-Clause-Clear): pbweave (aligner), pbdv (a
# reimplementation of DeepVariant 1.10.0's PACBIO model, CPU or NVIDIA GPU), pbsamoa, pbhappy, pbtruvari.
# One release tarball of binaries that need only glibc 2.17+ (pbdv carries its CUDA kernels and both models; on a
# GPU it needs the host's NVIDIA driver, which Apptainer's --nv provides). The base is upstream's own pb_wdl_base,
# by digest, for the samtools and bcftools (with plugins) the tasks use beside the tools. See containers/README.md
# and docs/DESIGN.md §17 (pblab 20261021), §18 (2026-11-11).
FROM debian:bookworm-slim@sha256:88200866dfff7ea7f5cbcb6ec7c8a701889efe6fe859fe64d6990e4b07ea4171 AS fetch
ARG PBLAB_VERSION=20261021
ARG PBLAB_SHA256=46cd78f3759a35be9c05d00035244cccda5300ea9df9f5f5ca749d179e9f0737
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl xz-utils \
 && rm -rf /var/lib/apt/lists/*
RUN set -eux; cd /tmp; \
    curl -fsSLo pblab.tar.xz "https://github.com/PacificBiosciences/pblab/releases/download/${PBLAB_VERSION}/pblab-${PBLAB_VERSION}-linux-x86_64.tar.xz"; \
    echo "${PBLAB_SHA256}  pblab.tar.xz" | sha256sum -c -; \
    mkdir /opt/pblab; tar -xJf pblab.tar.xz -C /opt/pblab --strip-components=1; \
    cd /opt/pblab/bin; sha256sum -c ../SHA256SUMS; \
    rm -rf /tmp/*

FROM quay.io/pacbio/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae
COPY --from=fetch /opt/pblab/bin/ /usr/local/bin/
COPY --from=fetch /opt/pblab/LICENSE /usr/share/licenses/pblab/LICENSE
RUN set -eux; pbweave --version | grep -q 'pbweave 26.3.99'; pbdv --version | grep -q 'pbdv 26.3.99'; \
    samtools --version | head -1; bcftools plugin -l | grep -qx setGT
