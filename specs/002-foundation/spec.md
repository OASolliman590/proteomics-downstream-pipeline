# Feature Specification: Runtime, package and contract foundation

**Phase:** 1. **Packet:** R01. **Status:** 1.2.0-frozen; implementation pending explicit GO.

## Scope

Deliver only FR-001–FR-010 for US1. [Packet index](../001-downstream-proteomics/packet-index.md), [ownership](../001-downstream-proteomics/packet-ownership.json), [data model](../001-downstream-proteomics/data-model.md) and [scientific contract](../001-downstream-proteomics/contracts/scientific-methods.md) define the shared interfaces. No implementation dispatch is authorized yet.

## Requirements

- **FR-001 — Installable maintained CLI:** The system MUST expose proteomics and python -m proteomics_pipeline; --help and --version exit 0 with identical version and list the five Phase 1 commands without reading data or invoking R.
- **FR-002 — Strict versioned configuration:** The system MUST accept the three resolved examples; reject unknown fields, old versions, malformed contrasts and conflicting primary fields with a JSON pointer and typed error; insert only the four documented defaults.
- **FR-003 — Typed run and stage states:** The system MUST keep failed child as FAILED, missing implementation as NOT_RUN and scientific inapplicability as INAPPLICABLE; no such stage becomes COMPLETED or PASS.
- **FR-004 — Safe Python-to-R bridge:** The system MUST pass literal argv with shell=false; roundtrip exact UTF-8 paths and preserve stdout/stderr, nonzero child exit and authentic session metadata.
- **FR-005 — Atomic artifact publication:** The system MUST permit at most one exclusive writer; an interrupted temporary never appears as a completed artifact or reusable cache; promotion preserves the verified hash.
- **FR-006 — Environment doctor:** The system MUST report actual paths/versions and per-capability availability; required missing R/package/resource exits 3; --help remains independent of doctor.
- **FR-007 — Maintained R package skeleton:** The system MUST install the real package where prerequisites permit; dispatch io_roundtrip and reproduce IDs, values, NA and metadata exactly; export no fabricated analysis function.
- **FR-008 — Independent test harness:** The system MUST retain all ten historical tests unchanged; the new harness discovers its actual tests, returns nonzero for a failing assertion and records missing R as NOT_RUN.
- **FR-009 — Baseline preservation registry:** The system MUST verify all pipeline/legacy/audit bytes against baseline; mutation of any protected file is detected without changing the original registry.
- **FR-010 — One-command developer bootstrap:** The system MUST bootstrap only project-local dependencies, record solved versions and invoke the real baseline/foundation gates; analysis commands never bootstrap automatically.

## Acceptance Scenarios

<a id="V001"></a>

### V001: Installable maintained CLI

**Fixture:** Fresh temporary Python environment, installed project, empty working directory and PATH without R.

**Oracle:** Subprocess return code/stdout from both installed entry points.

**Exact assertion:** Expose proteomics and python -m proteomics_pipeline; --help and --version exit 0 with identical version and list the five Phase 1 commands without reading data or invoking R.

**Negative case:** An unknown command exits 2; a missing future run capability exits 3/NOT_RUN rather than 0.

**Contract:** SM25; **owner:** R01; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V002"></a>

### V002: Strict versioned configuration

**Fixture:** The three contract examples; copies with typo primary_engnie, old schema_version, duplicate IDs and mismatched primary model.

**Oracle:** Draft202012Validator plus explicit reference/primary-role checks in semantic-validation.

**Exact assertion:** Accept the three resolved examples; reject unknown fields, old versions, malformed contrasts and conflicting primary fields with a JSON pointer and typed error; insert only the four documented defaults.

**Negative case:** Remove score_test: it resolves to off; explicitly set invalid score_test: reject rather than replacing it with off.

**Contract:** SM07; **owner:** R01; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V003"></a>

### V003: Typed run and stage states

**Fixture:** A foundation stage with a real child exit 7, an unavailable capability and a scientifically inapplicable optional request.

**Oracle:** The CLI/state transition table and actual subprocess exit.

**Exact assertion:** Keep failed child as FAILED, missing implementation as NOT_RUN and scientific inapplicability as INAPPLICABLE; no such stage becomes COMPLETED or PASS.

**Negative case:** A forged COMPLETED record with exit_code=7 fails stage schema/semantic verification.

**Contract:** SM25; **owner:** R01; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V004"></a>

### V004: Safe Python-to-R bridge

**Fixture:** Real R I/O stage in a directory named space Ω ; literal.txt, followed by an R stop() fixture.

**Oracle:** R-reported normalized argv and real stderr/sessionInfo from the same process.

**Exact assertion:** Pass literal argv with shell=false; roundtrip exact UTF-8 paths and preserve stdout/stderr, nonzero child exit and authentic session metadata.

**Negative case:** A metacharacter path cannot create a sentinel file outside the temp directory; missing R is NOT_RUN, not a mocked pass.

**Contract:** SM25; **owner:** R01; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V005"></a>

### V005: Atomic artifact publication

**Fixture:** Two actual subprocess writers target the same temporary run; terminate the first before promotion.

**Oracle:** Filesystem inspection plus independently computed SHA-256.

**Exact assertion:** Permit at most one exclusive writer; an interrupted temporary never appears as a completed artifact or reusable cache; promotion preserves the verified hash.

**Negative case:** Mutate a promoted artifact: verify exits 5 and downstream consumption is refused.

**Contract:** SM25; **owner:** R01; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V006"></a>

### V006: Environment doctor

**Fixture:** Actual foundation Python/R installation; repeat with R removed from PATH and with a configured missing resource.

**Oracle:** Independent which/version calls and resource SHA-256 calculation.

**Exact assertion:** Report actual paths/versions and per-capability availability; required missing R/package/resource exits 3; --help remains independent of doctor.

**Negative case:** A package name in configuration is not evidence of installation and must not be reported AVAILABLE without inspection.

**Contract:** SM10; **owner:** R01; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V007"></a>

### V007: Maintained R package skeleton

**Fixture:** Minimal real proteomicsCore package, an 8×12 TSV containing NA and a UTF-8 metadata field.

**Oracle:** R CMD INSTALL exit and independent read/write comparison, not a Python mock.

**Exact assertion:** Install the real package where prerequisites permit; dispatch io_roundtrip and reproduce IDs, values, NA and metadata exactly; export no fabricated analysis function.

**Negative case:** Invoke limma capability before R05: typed unavailable failure. Missing runtime or installer-required unresolved metadata remains NOT_RUN, never invented license content.

**Contract:** SM25; **owner:** R01; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V008"></a>

### V008: Independent test harness

**Fixture:** Existing tests/test_runner.py and a temporary intentionally failing assertion in an isolated harness fixture.

**Oracle:** unittest, pytest and real R testthat exit codes and collected test IDs.

**Exact assertion:** Retain all ten historical tests unchanged; the new harness discovers its actual tests, returns nonzero for a failing assertion and records missing R as NOT_RUN.

**Negative case:** An empty test selection is not a passed scientific suite; absent testthat cannot be converted to skipped-as-PASS.

**Contract:** SM25; **owner:** R01; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V009"></a>

### V009: Baseline preservation registry

**Fixture:** Fresh baseline files plus a scratch-only copy with one byte changed under each protected tree.

**Oracle:** SHA-256 of actual bytes compared with source_provenance and an independently captured protected-tree inventory.

**Exact assertion:** Verify all pipeline/legacy/audit bytes against baseline; mutation of any protected file is detected without changing the original registry.

**Negative case:** A copied registry regenerated from the mutated tree is rejected as an untrusted baseline replacement.

**Contract:** SM25; **owner:** R01; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V010"></a>

### V010: One-command developer bootstrap

**Fixture:** Clean project-local environment with an available documented R binary and network-enabled installation explicitly requested by the developer.

**Oracle:** Actual installation commands, solved-version inventory and foundation/baseline tests.

**Exact assertion:** Bootstrap only project-local dependencies, record solved versions and invoke the real baseline/foundation gates; analysis commands never bootstrap automatically.

**Negative case:** Missing R/permission/license-required metadata exits nonzero with unresolved/NOT_RUN and leaves global configuration unchanged.

**Contract:** SM25; **owner:** R01; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

## Boundary

No recovered source/audit changes, private-data transfer, invented license, fake backend or unsupported completion claim. A missing prerequisite is NOT_RUN, not a successful acceptance case. Only the Maintainer updates traceability after independent gate execution.

The first real files and commands are fixed in [IMPLEMENTATION_BRIEF_R01.md](../001-downstream-proteomics/IMPLEMENTATION_BRIEF_R01.md). Do not begin a later packet to make R01 look complete.
