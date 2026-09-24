# svx — structural variant merging for sawfish output (PacificBiosciences/svx)
# Early-stage tool: pin hard, expect to bump often. See containers/README.md.
# svx 0.6.0 sorts and tabix-indexes its own compressed output (--sort-output), so no htslib is needed here.
# Release tag has no "v" prefix; the asset name does.
FROM debian:bookworm-slim@sha256:88200866dfff7ea7f5cbcb6ec7c8a701889efe6fe859fe64d6990e4b07ea4171
ARG SVX_VERSION=0.6.0
ARG SVX_ASSET=svx-v${SVX_VERSION}-x86_64-unknown-linux-gnu.tar.gz
ARG SVX_SHA256=86746d20012aec97ab16312175608bf63ee7f5b41b8cf2fc72b76a8c82ce513e
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl \
 && rm -rf /var/lib/apt/lists/*
RUN set -eux; cd /tmp; \
    curl -fsSLO "https://github.com/PacificBiosciences/svx/releases/download/${SVX_VERSION}/${SVX_ASSET}"; \
    echo "${SVX_SHA256}  ${SVX_ASSET}" | sha256sum -c -; \
    tar -xzf "${SVX_ASSET}"; \
    install -m 0755 "svx-v${SVX_VERSION}-x86_64-unknown-linux-gnu/svx" /usr/local/bin/svx; \
    rm -rf /tmp/*; \
    svx --version
RUN useradd -m -u 1000 user
USER user
