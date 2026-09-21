# Parent implementation plan — v1.2.0-frozen

This is a future implementation plan, not implementation evidence. [Scientific methods](contracts/scientific-methods.md), [schema](contracts/analysis.schema.json), [semantic validation](contracts/semantic-validation.md), [data model](data-model.md) and [CLI](contracts/cli-and-artifacts.md) freeze the scientific/API contracts. [PHASES.md](PHASES.md) controls release scope.

Python owns safe configuration, validation orchestration, immutable filesystem state, provenance and offline report assembly. R owns the actual statistical computations through one installed proteomicsCore package. Numeric files and typed stage JSON are the interface; no service, database, cluster workflow engine or user plugin framework is required.

R01 creates only a genuine package/CLI/I/O/runtime foundation. R02–R04 supply canonical inputs, preprocessing/QC and validated plans. R05 supplies primary limma/TREAT and exact contrasts. R10a provides minimal offline HTML/status integration. Qualified alternative engines, resources, pathways and response belong to Phase 2. Full reporting, locked reproduction/calibration and legacy reconciliation belong to Phase 3. The exact files, dependencies and only legal parallel group are in [packet-index.md](packet-index.md) and [ownership](packet-ownership.json); do not reopen shared files.

An AnalysisPlan is frozen after deterministic intake/preprocessing/design calculations and before the first fit. Every fit consumes its validated plan/input hashes. Required errors stop dependants; optional scientific inapplicability and operational failure stay distinct. The report consumes artifacts and does not recalculate statistical claims in templates. A zero-discovery analysis can be completed; a missing/failed analysis cannot be filled with healthy empty tables.

Phase 1 qualification uses small direct-reference/scientific-negative tests and actual solved environment records. Full lock restoration, 1,000-dataset calibration, cross-platform evidence and performance gates are R11. No package version is called tested until it was actually installed/executed; no gene-set snapshot is called pinned without local bytes and a verified digest. The [validation strategy](validation-strategy.md) fixes thresholds before execution.

A frozen packet implementer writes tests from its explicit oracles, implements only its allowlist and returns evidence. The Maintainer independently verifies and records completion. Contract contradictions pause the dependent packet for an ADR/correction, not an invented convenient default. A missing external prerequisite is reported; the rest of the packet can proceed, but no missing gate is marked passed.

Current next step is independent audit and Maintainer freeze, then [R01 only](IMPLEMENTATION_BRIEF_R01.md). No maintained engine code, pyproject.toml, R package internals or renv.lock is created in this hardening pass.
