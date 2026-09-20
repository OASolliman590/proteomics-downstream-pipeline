# Implementation entry point

**Objective:** mature the recovered baseline into a streamlined, scientifically defensible, full downstream protein-abundance analysis pipeline. Deliver development through a corrected and frozen Spec Kit, then execute one packet at a time under maintainer review with independent audits at the specification-freeze and risky-diff gates.

This package is a specification and implementation mandate. It does not claim the future functionality already exists. The recovered Python intake and ten baseline tests exist; the new CLI, R package, contracts and validation described here are implementation targets.

## Read in this order

1. [Constitution](../../.specify/memory/constitution.md) and [parent specification](spec.md).
2. [Roadmap](roadmap.md), [architecture plan](plan.md), [scientific contract](contracts/scientific-methods.md), [data model](data-model.md) and [validation strategy](validation-strategy.md).
3. [CLI/file contract](contracts/cli-and-artifacts.md), [configuration schema](contracts/analysis.schema.json), and [research decisions](research.md).
4. [Corrected packet index](packet-index.md), then the next numbered slice's spec.md, plan.md and tasks.md.
5. [Progress record](progress.md) and [traceability register](traceability.json). These must reflect verified work throughout implementation.

## Completion boundary

All twelve slices are mandatory software deliverables, including tested method-specific adapters, blocked-design inference, valid pathway dispatch, descriptive treatment-response analysis, independent-score validation, reproducible reporting, CI/calibration and archived-study regression. Some methods or analyses are optional *per study* because their scientific prerequisites are not always satisfied. Their code and applicability tests are still required. Missing permissions or study metadata must be reported accurately and must not cause unrelated development to stop.

Explicitly excluded: raw instrument conversion/search, peptide identification, de novo protein inference, phosphosite localization/PTM workflows, single-cell proteomics, untargeted cross-omics integration, biomarker classifiers, causal drug-mechanism claims, network-service UI and public deployment. Peptide/PSM **counts already supplied with quantified proteins** are supported for DEqMS; peptide-level abundance modeling is a future separately specified extension. Generic protein-level matrix intake supports other species and assays only through validated contracts.

## Execution location

Use a source-only checkout rooted at `<project-root>` for implementation. The package at `<private-archive-root>` remains the historical snapshot. Work in an isolated checkout and write any reviewed, versioned successor under `<release-root>/releases/`. Do not mutate the original ZIP or `<private-evidence-root>/original/`.

External storage may require a repository-specific Git `safe.directory` setting when its filesystem lacks ownership metadata. Never add a global wildcard trust rule. Only the maintainer accesses private study data for final regression.

## Spec Kit format and provenance

The structure follows the official [Spec Kit workflow](https://github.github.com/spec-kit/) and [spec-of-specs roadmap convention](https://github.github.com/spec-kit/concepts/spec-of-specs.html), with a project constitution and numbered features containing spec/plan/tasks. Artifacts were authored for this project; the Spec Kit CLI/integration has not been installed or executed. No generated upstream template text or slash-command registration is claimed. Use ordinary file-based execution if no Spec Kit integration is available; never run `specify init --force` over this package.
