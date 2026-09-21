# Feature Specification: Reproducibility, calibration and continuous validation

**Phase:** 3. **Packet:** R11. **Status:** 1.2.0-frozen; implementation pending prerequisites and separate authorization.

## Scope

Deliver only FR-101–FR-110 for US6. [Packet index](../001-downstream-proteomics/packet-index.md), [ownership](../001-downstream-proteomics/packet-ownership.json), [data model](../001-downstream-proteomics/data-model.md) and [scientific contract](../001-downstream-proteomics/contracts/scientific-methods.md) define the shared interfaces. No implementation dispatch is authorized yet.

## Requirements

- **FR-101 — Pinned compatible environments:** The system MUST restore compatible Python/R/Bioconductor environments and exercise the declared APIs; environment.yml alone is not a lock.
- **FR-102 — Immutable resource-offline reproduction:** The system MUST reproduce semantic tables/masks/universes using unchanged resource/code/environment hashes; separate timestamps/cosmetic image metadata from scientific values.
- **FR-103 — Resume and cache invalidation:** The system MUST reuse only compatible verified stage outputs; invalidate each affected stage and dependant; never reuse an incomplete temporary.
- **FR-104 — Scientific golden reference matrix:** The system MUST freeze versions/tolerances before candidate comparison and record actual fieldwise differences and unsupported combinations.
- **FR-105 — Null and mixture calibration:** The system MUST meet prespecified core error upper-bound ≤.075 and CI-coverage target .93–.97 with adequate independent-dataset precision, reporting power/bias/failures as well; retain failed qualifications.
- **FR-106 — Correlated pathway calibration:** The system MUST assess CAMERA/ROAST under their actual hypotheses and record Monte Carlo precision; retain fgsea exploratory gene-set-null limits.
- **FR-107 — Cross-platform CI:** The system MUST run real core Python/R examples on both systems and distinguish expensive release calibration from smoke CI; absent jobs remain NOT_RUN.
- **FR-108 — Performance and resource benchmark:** The system MUST assess ≤30 minutes and ≤8 GiB with hardware/settings recorded, timing specialized engines/calibration separately; report any target failure.
- **FR-109 — Report numerical/visual validation:** The system MUST verify numeric labels/coordinates/source links and record actual visual inspection, readability and failed-layout findings.
- **FR-110 — Machine-readable validation ledger:** The system MUST record requirement/task/acceptance IDs, expected oracle, actual evidence and reviewed tree; preserve distinct PASS/FAIL/NOT_RUN.

## Acceptance Scenarios

<a id="V101"></a>

### V101: Pinned compatible environments

**Fixture:** Fresh local Python environment and R library restored from actual requirements.lock/renv.lock.

**Oracle:** Installed versions and hashes compared with lock contents and real adapter API probes.

**Exact assertion:** Restore compatible Python/R/Bioconductor environments and exercise the declared APIs; environment.yml alone is not a lock.

**Negative case:** A floating dependency, incompatible R/Bioconductor pair or unsupported output column fails qualification rather than being silently accepted.

**Contract:** SM25; **owner:** R11; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V102"></a>

### V102: Immutable resource-offline reproduction

**Fixture:** Complete verified run, local resource snapshots and a second fresh run with network disabled.

**Oracle:** Sorted scientific keys and independently computed fieldwise semantic comparison under frozen tolerances.

**Exact assertion:** Reproduce semantic tables/masks/universes using unchanged resource/code/environment hashes; separate timestamps/cosmetic image metadata from scientific values.

**Negative case:** A changed resource byte or attempted network fetch fails; changed sample n is not ignored as a cosmetic difference.

**Contract:** SM14; **owner:** R11; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V103"></a>

### V103: Resume and cache invalidation

**Fixture:** Completed run plus copies with changed config, input bytes, one code file, resource hash and interrupted temporary stage.

**Oracle:** Dependency graph reachability and immutable input/output hashes.

**Exact assertion:** Reuse only compatible verified stage outputs; invalidate each affected stage and dependant; never reuse an incomplete temporary.

**Negative case:** A resume mixing stale and new scientific lineage or treating directory existence as a cache hit fails.

**Contract:** SM25; **owner:** R11; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V104"></a>

### V104: Scientific golden reference matrix

**Fixture:** Every eligible backend/hypothesis/design case in the declared reference matrix, each using an independent reference script.

**Oracle:** Direct official packages or mathematical calculations authored separately from production adapters.

**Exact assertion:** Freeze versions/tolerances before candidate comparison and record actual fieldwise differences and unsupported combinations.

**Negative case:** An oracle importing the production adapter under test or widened post-failure tolerance without an ADR fails.

**Contract:** SM08; **owner:** R11; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V105"></a>

### V105: Null and mixture calibration

**Fixture:** Future R11 release simulator with ≥1000 independent fast-core datasets per declared null/mixture scenario and fixed seed bank. No large dataset is created in this spec pass.

**Oracle:** Dataset-level FDP, all-null any-rejection probability and independent binomial/bootstrap Monte Carlo intervals defined in validation-strategy.

**Exact assertion:** Meet prespecified core error upper-bound ≤.075 and CI-coverage target .93–.97 with adequate independent-dataset precision, reporting power/bias/failures as well; retain failed qualifications.

**Negative case:** Selecting favorable seeds or replacing dataset-level FDR with per-protein false-positive rate fails. Phase 1 is not required to run this release calibration.

**Contract:** SM12; **owner:** R11; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V106"></a>

### V106: Correlated pathway calibration

**Fixture:** Future R11 prespecified correlated gene-set null/spike-in simulator with competitive versus self-contained hypotheses explicitly distinct.

**Oracle:** Known generating memberships/effects and method-null-specific dataset-level error/power calculations.

**Exact assertion:** Assess CAMERA/ROAST under their actual hypotheses and record Monte Carlo precision; retain fgsea exploratory gene-set-null limits.

**Negative case:** Demanding or claiming a universal sample-permutation guarantee for fgsea, or relabeling a failed core test stress-only afterward, fails.

**Contract:** SM16; **owner:** R11; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V107"></a>

### V107: Cross-platform CI

**Fixture:** Fresh Linux and Windows jobs with actual Python/R installations and public synthetic examples.

**Oracle:** Actual job logs, command exits and artifact manifests.

**Exact assertion:** Run real core Python/R examples on both systems and distinguish expensive release calibration from smoke CI; absent jobs remain NOT_RUN.

**Negative case:** Python-only CI or skipped R tests cannot support a cross-platform statistical PASS claim.

**Contract:** SM25; **owner:** R11; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V108"></a>

### V108: Performance and resource benchmark

**Fixture:** Future R11 generated 20000-feature×100-observation, eight-contrast, 2000-set benchmark on documented four-core hardware.

**Oracle:** Actual wall clock and peak process-tree memory measurement for QC+limma+CAMERA+report.

**Exact assertion:** Assess ≤30 minutes and ≤8 GiB with hardware/settings recorded, timing specialized engines/calibration separately; report any target failure.

**Negative case:** Downsizing the benchmark to make the target pass or extrapolating timing without an actual run fails.

**Contract:** SM25; **owner:** R11; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V109"></a>

### V109: Report numerical/visual validation

**Fixture:** Rendered full/null/partial reports and requested publication figures from real accepted runs.

**Oracle:** Exact table-to-figure data comparison plus human/independent reviewer inspection of rendered output.

**Exact assertion:** Verify numeric labels/coordinates/source links and record actual visual inspection, readability and failed-layout findings.

**Negative case:** A generated filename alone is not rendered visual evidence; uninspected reports remain NOT_RUN for the visual gate.

**Contract:** SM25; **owner:** R11; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V110"></a>

### V110: Machine-readable validation ledger

**Fixture:** One real passing, one failing and one deliberately unavailable acceptance execution record.

**Oracle:** Observed command/exit/version/hash and reviewer conclusion for each ID.

**Exact assertion:** Record requirement/task/acceptance IDs, expected oracle, actual evidence and reviewed tree; preserve distinct PASS/FAIL/NOT_RUN.

**Negative case:** A proposed path, worker assertion or copied expected output cannot populate a PASS evidence record.

**Contract:** SM25; **owner:** R11; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

## Boundary

No recovered source/audit changes, private-data transfer, invented license, fake backend or unsupported completion claim. A missing prerequisite is NOT_RUN, not a successful acceptance case. Only the Maintainer updates traceability after independent gate execution.
