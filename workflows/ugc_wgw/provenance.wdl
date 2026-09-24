version 1.0

# ugc-pacbio-wgw — write <subject>.<stage>.ugc_wgw_manifest.json, the workflow-side half of run provenance
# origin: none
# see docs/DESIGN.md §9.2

import "../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_manifest_write {
  meta {
    description: "Write the workflow-side provenance record for one run of a ugc-wgw stage; the driver embeds it in run_manifest.json"
    outputs: {
      ugc_wgw_manifest: {
        description: "JSON: schema, ugc_wgw_version, stage, subject, cohort members, upstream workflow name and version, host, timestamp"
      }
    }
  }

  parameter_meta {
    ugc_wgw_version: {
      description: "ugc-pacbio-wgw version the run was submitted with; the driver checks it"
    }
    stage: {
      description: "Stage name (singleton, upstream, cohort_call, downstream, cohort_merge, cohort_freq, assembly)"
    }
    subject_type: {
      description: "sample or cohort"
    }
    subject_id: {
      description: "Sample ID or cohort ID"
    }
    member_ids: {
      description: "Sample IDs of the cohort members (cohort stages only)"
    }
    upstream_workflow_name: {
      description: "Upstream workflow name, when the stage wraps one (e.g. humanwgs_singleton)"
    }
    upstream_workflow_version: {
      description: "Upstream workflow version string, when the stage wraps one (e.g. v3.3.1)"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String ugc_wgw_version
    String stage
    String subject_type
    String subject_id
    Array[String] member_ids = []
    String? upstream_workflow_name
    String? upstream_workflow_version
    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  Int mem_gb = 1
  Int disk_size = 1

  command <<<
    set -euo pipefail

    python3 --version >&2

    python3 - "~{write_lines(member_ids)}" > "~{subject_id}.~{stage}.ugc_wgw_manifest.json" << 'EOF'
    import datetime, json, socket, sys
    members = [line.strip() for line in open(sys.argv[1]) if line.strip()]
    doc = {
        "schema": 1,
        "ugc_wgw_version": "~{ugc_wgw_version}",
        "stage": "~{stage}",
        "subject": {"type": "~{subject_type}", "id": "~{subject_id}"},
        "cohort_members": members,
        "upstream": {
            "workflow_name": "~{default="" upstream_workflow_name}",
            "workflow_version": "~{default="" upstream_workflow_version}",
        },
        "written_at": datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "host": socket.gethostname(),
    }
    json.dump(doc, sys.stdout, indent=2, sort_keys=True)
    sys.stdout.write("\n")
    EOF
  >>>

  output {
    File ugc_wgw_manifest = "~{subject_id}.~{stage}.ugc_wgw_manifest.json"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}
