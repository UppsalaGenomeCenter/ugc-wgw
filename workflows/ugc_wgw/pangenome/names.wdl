version 1.0
# ugc-pacbio-wgw — cohort pangenome bookkeeping: Cactus-safe genome names, the sample rename and unphase lists, per-chromosome BEDs and memory
# origin: none
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8.2
import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_pangenome_names {
  meta {
    description: "Cactus reserves the dot for the haplotype suffix (SAMPLE.1, SAMPLE.2) while ugc-wgw sample IDs may contain dots: every member gets a safe name (every character outside [A-Za-z0-9_] becomes an underscore), the task fails on a collision, and writes the genome-name list of the seqfile (<safe>.1, <safe>.2 per member, in cohort order), the map back to the sample IDs (for bcftools reheader), the non-trio members to unphase, one BED per requested chromosome for the frequency shards, and the memory each chromosome task gets (floor, or per-100-Mb rate times the chromosome length, whichever is larger)."
    outputs: {
      genome_names: {
        description: "<safe>.1 and <safe>.2 per member, cohort order: the seqfile's genome names"
      },
      name_map: {
        description: "sample_id, safe_name, trio per member"
      },
      rename_samples: {
        description: "safe_name<TAB>sample_id per member, for bcftools reheader --samples"
      },
      unphase_samples: {
        description: "Sample IDs of the members that are not trio-binned (their genotypes are unphased)"
      },
      chrom_beds: {
        description: "One BED per requested chromosome (whole chromosome), in the chromosomes' order"
      },
      chrom_mem_gb: {
        description: "JSON object chromosome -> memory (GB) for its Cactus task"
      }
    }
  }

  parameter_meta {
    cohort_id: {
      description: "Cohort ID"
    }
    sample_ids: {
      description: "Members in cohort order"
    }
    trios: {
      description: "Whether each member's assembly is trio-binned, same order"
    }
    ref_index: {
      description: "The assembly reference's .fai (chromosome lengths)"
    }
    chromosomes: {
      description: "The chromosomes to build"
    }
    chrom_mem_gb_floor: {
      description: "Memory (GB) of a chromosome task at least"
    }
    chrom_mem_gb_per_100mb: {
      description: "Memory (GB) per 100 Mb of chromosome length"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String cohort_id
    Array[String] sample_ids
    Array[Boolean] trios
    File ref_index
    Array[String] chromosomes
    Int chrom_mem_gb_floor = 64
    Int chrom_mem_gb_per_100mb = 32
    RuntimeAttributes runtime_attributes
  }

  command <<<
    set -euo pipefail
    mkdir -p beds
    python3 - "~{write_lines(sample_ids)}" "~{write_lines(trios)}" "~{ref_index}" "~{write_lines(chromosomes)}" \
      ~{chrom_mem_gb_floor} ~{chrom_mem_gb_per_100mb} "~{cohort_id}" <<'PY'
import json, re, sys
ids = [l.strip() for l in open(sys.argv[1]) if l.strip()]
trios = [l.strip() == "true" for l in open(sys.argv[2]) if l.strip()]
lengths = {}
for line in open(sys.argv[3]):
    f = line.split("\t")
    lengths[f[0]] = int(f[1])
chroms = [l.strip() for l in open(sys.argv[4]) if l.strip()]
floor, per100 = int(sys.argv[5]), int(sys.argv[6])
missing = [c for c in chroms if c not in lengths]
if missing:
    sys.exit(f"chromosomes absent from the reference index: {' '.join(missing)}")
safe = [re.sub(r"[^A-Za-z0-9_]", "_", s) for s in ids]   # Cactus-safe: alphanumerics and underscore
seen = {}
for sid, name in zip(ids, safe):
    if name in seen:
        sys.exit(f"Cactus name collision: {sid} and {seen[name]} both become {name}")
    if name == "GRCh38" or name.startswith("_"):
        sys.exit(f"sample {sid}: the name {name} is reserved")
    seen[name] = sid
with open("genome_names.txt", "w") as fh:
    for name in safe:
        fh.write(f"{name}.1\n{name}.2\n")
with open("~{cohort_id}.pangenome.name_map.tsv", "w") as fh:
    fh.write("sample_id\tsafe_name\ttrio\n")
    for sid, name, trio in zip(ids, safe, trios):
        fh.write(f"{sid}\t{name}\t{'true' if trio else 'false'}\n")
with open("rename_samples.txt", "w") as fh:
    for sid, name in zip(ids, safe):
        fh.write(f"{name}\t{sid}\n")
with open("unphase_samples.txt", "w") as fh:
    for sid, trio in zip(ids, trios):
        if not trio:
            fh.write(sid + "\n")
mem = {}
for i, c in enumerate(chroms):
    with open(f"beds/{i:03d}.{c}.bed", "w") as fh:
        fh.write(f"{c}\t0\t{lengths[c]}\n")
    mem[c] = max(floor, -(-lengths[c] * per100 // 100_000_000))
with open("chrom_mem_gb.json", "w") as fh:
    json.dump(mem, fh, indent=1, sort_keys=True)
print(f"{sys.argv[7]}: {len(ids)} members, {sum(trios)} trio-binned, {len(chroms)} chromosomes; memory {mem}", file=sys.stderr)
PY
  >>>

  output {
    Array[String] genome_names = read_lines("genome_names.txt")
    File name_map = "~{cohort_id}.pangenome.name_map.tsv"
    File rename_samples = "rename_samples.txt"
    Array[String] unphase_samples = read_lines("unphase_samples.txt")
    Array[File] chrom_beds = glob("beds/*.bed")
    Map[String, Int] chrom_mem_gb = read_json("chrom_mem_gb.json")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
    cpu: 1
    memory: "2 GiB"
    disk: "2 GB"
    disks: "local-disk 2 HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}
