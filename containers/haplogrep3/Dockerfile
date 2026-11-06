# haplogrep3 — mitochondrial haplogroup classification of mitorsaw's haplotypes (genepi/haplogrep3, MIT)
# A Java tool: eclipse-temurin's JRE 21 (Ubuntu 24.04) by digest, the release tarball by sha256, and the two
# phylogenetic trees the tasks use installed at build time: haplogrep3 downloads a tree from GitHub the first
# time it is named, and the HPC has no network. The wrapper in /usr/local/bin quotes its arguments and sizes the
# heap (the shipped launcher word-splits them and asks for 10 GB). See containers/README.md.
FROM eclipse-temurin:21-jre-noble@sha256:0eda5830ed7b669abc246a8dc520ea8a23cab06f9a8d4b6964109e2f2394e46f
ARG HAPLOGREP3_VERSION=3.3.2
ARG HAPLOGREP3_SHA256=957165143ddd7d6ad9f34cc3af56b784fc236081e031c99bea4a574d6ad70959
RUN set -eux; mkdir -p /opt/haplogrep3; cd /opt/haplogrep3; \
    curl -fsSLo haplogrep3.tar.gz "https://github.com/genepi/haplogrep3/releases/download/v${HAPLOGREP3_VERSION}/haplogrep3-${HAPLOGREP3_VERSION}-linux.tar.gz"; \
    echo "${HAPLOGREP3_SHA256}  haplogrep3.tar.gz" | sha256sum -c -; \
    tar --no-same-owner -xzf haplogrep3.tar.gz; \
    rm haplogrep3.tar.gz; \
    java -jar haplogrep3.jar install-tree phylotree-fu-rcrs@1.3 phylotree-rcrs@17.3; \
    chmod -R a+rX /opt/haplogrep3
COPY haplogrep3 /usr/local/bin/haplogrep3
RUN set -eux; chmod 0755 /usr/local/bin/haplogrep3; \
    haplogrep3 --version 2>&1 | grep -q "Haplogrep 3 ${HAPLOGREP3_VERSION}"; \
    haplogrep3 trees 2>&1 | grep -q 'phylotree-fu-rcrs@1.3'; \
    haplogrep3 trees 2>&1 | grep -q 'phylotree-rcrs@17.3'
# the base image's uid 1000 (ubuntu); the tasks run as the submitting user under Apptainer anyway
USER 1000
