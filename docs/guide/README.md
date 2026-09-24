# Operator guide

This guide is for the people who run ugc-pacbio-wgw on the HPC and for the
analysts who read its results. It explains what the pipeline produces, how the
driver sequences the work, and what to do when something fails. The driver's
reference page is [`../../bin/README.md`](../../bin/README.md).

## Chapters

| Chapter | Read it when |
|---|---|
| [00 Quick start](00-quick-start.md) | You have a bundle and want the first cohort's results: the eight steps with what to expect at each. |
| [01 Overview](01-overview.md) | You want the one-page picture: what runs where, and the five things an operator does. |
| [02 Concepts](02-concepts.md) | You need the vocabulary: stage, subject, mode, version, attempt, manifest, frozen cohort. |
| [03 Choosing a mode](03-choosing-a-mode.md) | You are deciding between standalone and joint calling, or adding assembly. |
| [04 Install](04-install.md) | A bundle has arrived and must be installed on the HPC. |
| [05 Project setup](05-project-setup.md) | You are creating a project, registering samples and freezing a cohort. |
| [06 Running](06-running.md) | You are submitting, watching, cancelling or retrying runs. |
| [07 Results](07-results.md) | You are looking for a file, or want to know what a manifest says. |
| [08 Stages](08-stages.md) | You want to know what each stage runs, which switches it has and what it costs. |
| [09 Troubleshooting](09-troubleshooting.md) | Something failed. |
| [10 Upgrading](10-upgrading.md) | A new bundle version must coexist with, or replace, the current one. |
| [11 Reference](11-reference.md) | Every command, flag, config key and file name in one place. |
| [12 Task resources](12-resources.md) | A job was refused or killed for cores, memory or time, or the requests must fit the cluster. |
| [Workflow call graphs](workflow-graphs.md) | Generated exact call graphs of the seven entrypoints. |
| [Task resources, inventory](task-resources.md) | Generated table of what every task asks SLURM for. |

## Reading the diagrams

Diagrams are Mermaid blocks. GitHub, Forgejo and most editors render them in
place; on the HPC there is no browser, so the core diagrams are
followed by a plain-text version introduced by "Text version:". The guide is
shipped inside every bundle: on the HPC it is at
`<prefix>/versions/<v>/code/docs/guide/`.

## Conventions

| Placeholder | Meaning |
|---|---|
| `<prefix>` | The install root on the HPC, e.g. `/proj/ugc`. Holds `versions/`, `current`, `call_cache/`, `references/`. |
| `<v>` | A ugc-pacbio-wgw version, e.g. `0.2.0`. Appears in install paths and in every results path. |
| `<project>` | A directory holding `.ugc-wgw/` (config, state database, logs). One project per campaign. |
| `<results>` | The results root, `<project>` unless `ugc-wgw init --results` said otherwise. |
| `ugc-wgw` | `<prefix>/current/code/bin/ugc-wgw`. Put it on `PATH` or define an alias. |
| `<sample_id>`, `<cohort_id>` | Identifiers matching `[A-Za-z0-9._-]+`; both become directory names. |

Commands shown in `bash` blocks are meant to be copied. Blocks without a
language tag are file layouts, synopses or program output.

## The quick path

Chapter 00 is the walkthrough: install, project, samples, cohort, a dry
run, `submit` under `tmux`, watching and retrying, the result pages. In
three lines:

```bash
./install-bundle.sh --bundle ugc-pacbio-wgw-0.3.0.tar --prefix /proj/ugc --references /proj/ugc/references --site site.cfg --activate
ugc-wgw init <project> --install /proj/ugc/current --ref-map <ref_map from the installer's report>
ugc-wgw samples add samples.tsv && ugc-wgw cohort freeze C1 --samples cohort.txt && ugc-wgw -v submit --mode standalone --cohort C1
```
