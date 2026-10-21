version 1.0
# ugc-pacbio-wgw — a task that fails with a message: WDL 1.0 has no assert, so an entrypoint calls it under a condition it must refuse
# origin: none
# see docs/DESIGN.md §13; docs/ENTRYPOINTS.md §8
import "../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_fail {
  meta {
    description: "Exit 1 with the message on stderr. Called under a condition the workflow refuses (an unknown assembler); every other call that needs its result then never runs."
    outputs: {
      never: {
        description: "Never produced: the task always fails"
      }
    }
  }

  parameter_meta {
    message: {
      description: "What to say before failing"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String message
    RuntimeAttributes runtime_attributes
  }

  command <<<
    set -euo pipefail
    echo "ugc-wgw: ~{message}" >&2
    exit 1
  >>>

  output {
    String never = read_string(stdout())
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
    cpu: 1
    memory: "1 GiB"
    disk: "1 GB"
    disks: "local-disk 1 HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: 0
    awsBatchRetryAttempts: 0  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}
