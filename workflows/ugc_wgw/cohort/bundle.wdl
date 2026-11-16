version 1.0

# ugc-pacbio-wgw — many small per-member files as one tar, so that the task reading them has one input
# origin: none
# see docs/DESIGN.md §7, §18 (2026-11-16)

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_files_bundle {
  meta {
    description: "One array of files as one uncompressed tar, the entries named <index>__<file name> with the index zero-padded so that their sorted order is the array's. miniwdl binds every input file of a task into its container, and the whole bind list is on the command line Slurm is given: past about a megabyte (some 3,400 files at this pipeline's path lengths) sbatch refuses the job. A task that reads several per-member arrays takes them as bundles instead; each bundle task has one array, so the ceiling is one file per member."
    outputs: {
      bundle: {
        description: "<name>.tar"
      }
    }
  }

  parameter_meta {
    files: {
      description: "The files, in the order to keep (it may be empty: the tar then has no entry)"
    }
    name: {
      description: "The tar's name without the extension"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[File] files
    String name
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(files, "GB") * 2 + 10)

  command <<<
    set -euo pipefail

    python3 - "~{write_lines(files)}" "~{name}.tar" <<'EOF'
    import os, sys, tarfile
    paths = [l.rstrip("\n") for l in open(sys.argv[1]) if l.strip()]
    with tarfile.open(sys.argv[2], "w", dereference=True) as tar:
        for i, p in enumerate(paths):
            tar.add(p, arcname=f"{i:06d}__{os.path.basename(p)}", recursive=False)
    print(f"{len(paths)} files", file=sys.stderr)
    EOF
  >>>

  output {
    File bundle = "~{name}.tar"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
    cpu: 1
    memory: "2 GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}
