# somalier — relatedness, sex and ancestry from ~17 000 informative sites (brentp/somalier, MIT)
# One static binary per release (no shared libraries); the sites file and the 1000 Genomes labelled set are
# reference files (references.lock, ugc-wgw-qc), not part of the image. See containers/README.md.
# `somalier --version` prints the version and exits 1 ("unknown program"), hence the `|| true`.
FROM debian:bookworm-slim@sha256:88200866dfff7ea7f5cbcb6ec7c8a701889efe6fe859fe64d6990e4b07ea4171
ARG SOMALIER_VERSION=0.3.7
ARG SOMALIER_SHA256=c4495dbd6f549d3923bed4eda578967cd19043241d115c64bfbda13275dd4345
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl \
 && rm -rf /var/lib/apt/lists/*
RUN set -eux; cd /tmp; \
    curl -fsSLo somalier "https://github.com/brentp/somalier/releases/download/v${SOMALIER_VERSION}/somalier"; \
    echo "${SOMALIER_SHA256}  somalier" | sha256sum -c -; \
    install -m 0755 somalier /usr/local/bin/somalier; \
    rm -rf /tmp/*; \
    somalier --version 2>&1 | grep -q "somalier version: ${SOMALIER_VERSION}"
RUN useradd -m -u 1000 user
USER user
