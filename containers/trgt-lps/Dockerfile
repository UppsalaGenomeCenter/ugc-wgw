# trgt-lps — longest pure segment detection on merged TRGT VCFs (PacificBiosciences/trgt-lps)
# Runs after `trgt merge` in the cohort_merge stage. See containers/README.md.
# The release asset is a gzipped, statically linked x86-64 binary; the tag has no "v" prefix.
FROM debian:bookworm-slim@sha256:88200866dfff7ea7f5cbcb6ec7c8a701889efe6fe859fe64d6990e4b07ea4171
ARG TRGT_LPS_VERSION=0.10.0
ARG TRGT_LPS_ASSET=trgt-lps-v${TRGT_LPS_VERSION}.gz
ARG TRGT_LPS_SHA256=204e3a34e10c7bfe2e851ff3f02154a6c11f4dcc501626192f8583ab1279c941
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl \
 && rm -rf /var/lib/apt/lists/*
RUN set -eux; cd /tmp; \
    curl -fsSLO "https://github.com/PacificBiosciences/trgt-lps/releases/download/${TRGT_LPS_VERSION}/${TRGT_LPS_ASSET}"; \
    echo "${TRGT_LPS_SHA256}  ${TRGT_LPS_ASSET}" | sha256sum -c -; \
    gunzip -c "${TRGT_LPS_ASSET}" > /usr/local/bin/trgt-lps; \
    chmod 0755 /usr/local/bin/trgt-lps; \
    rm -rf /tmp/*; \
    trgt-lps --version
RUN useradd -m -u 1000 user
USER user
