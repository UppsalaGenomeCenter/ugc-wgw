# containers/

One directory per tool we build ourselves; upstream images are never rebuilt.
Each `Dockerfile` pins its base image by digest and the tool's release asset
by sha256, so a rebuild is reproducible. `scripts/build-containers.sh [tool]`
builds, pushes and prints the pushed `<registry>/<tool>@sha256:<digest>`;
the digest is what the WDL tasks and `image_manifest.ugc-wgw.txt` reference,
never a tag. The registry is `UGC_WGW_REGISTRY` (the script's default names
the registry the images were last pushed to); on the HPC it is irrelevant,
because miniwdl reads the SIF cache the bundle installed. Adding a tool: a
directory here, the build, the manifest line, the WDL task with the digest,
and the changelog.
