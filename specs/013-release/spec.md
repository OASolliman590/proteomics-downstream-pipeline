# Feature Specification: Archived-study regression, documentation and versioned successor

**Phase:** 3. **Packet:** R12. **Status:** 1.2.0-frozen; implementation pending prerequisites and separate authorization.

## Scope

Deliver only FR-111–FR-120 for US6. [Packet index](../001-downstream-proteomics/packet-index.md), [ownership](../001-downstream-proteomics/packet-ownership.json), [data model](../001-downstream-proteomics/data-model.md) and [scientific contract](../001-downstream-proteomics/contracts/scientific-methods.md) define the shared interfaces. No implementation dispatch is authorized yet.

## Requirements

- **FR-111 — Private local legacy adapter regression:** The system MUST verify private reconstruction locally without exposing private inputs/logs to implementers/reviewers/docs workers and without modifying originals.
- **FR-112 — Legacy core numerical comparison:** The system MUST compare declared coefficients/BH/settings, preserve differences and record NOT_RUN/NOT_REPRODUCED when historical runtime/resources are unavailable.
- **FR-113 — Scientific-change reconciliation:** The system MUST explain the Negr1 auxiliary exception, 2/70 versus 3/71 and 6/70 versus 7/71 arithmetic, and intentional model/pathway/response differences without production constants.
- **FR-114 — Full maintained study run:** The system MUST execute and verify the local study as new artifacts, retaining inherited unknown metadata and separation from the archived analysis.
- **FR-115 — Comprehensive usage and methods docs:** The system MUST document supported assays/designs/engines, rejection messages and exact tested commands; separate phase milestones from completed v1.0.
- **FR-116 — Spec-task-evidence reconciliation:** The system MUST retain exactly 120 unique requirement/task/acceptance identities, reviewed evidence or explicit unresolved status, and no silently retired ID.
- **FR-117 — Source release hygiene:** The system MUST exclude private data, credentials, environments and private execution logs; retain only permitted resource content/terms and explicitly unresolved ownership/license.
- **FR-118 — Versioned successor export:** The system MUST export only a new versioned successor directory after gates and permission; never overwrite the recovered archive package.
- **FR-119 — Final integrated acceptance run:** The system MUST re-execute release gates on that exact tree, including R10a compatibility, and report external limits separately without calling missing core work complete.
- **FR-120 — Review-and-handoff completion:** The system MUST close only verified requested scope, record real delivery paths and remaining blockers, and perform no unauthorized commit/push/license/private-data action.

## Acceptance Scenarios

<a id="V111"></a>

### V111: Private local legacy adapter regression

**Fixture:** Maintainer-only external original archive/manifest; absent from worker inputs and source checkout.

**Oracle:** Maintainer local original ZIP/165-file hashes and actual recovered input values/masks/metadata.

**Exact assertion:** Verify private reconstruction locally without exposing private inputs/logs to implementers/reviewers/docs workers and without modifying originals.

**Negative case:** Missing archive or permission is NOT_RUN, not an invented reconstructed matrix or requested private-data upload.

**Contract:** SM01; **owner:** R12; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V112"></a>

### V112: Legacy core numerical comparison

**Fixture:** Maintainer-local matched-universe legacy settings and an actually provisioned compatible historical R/resource environment.

**Oracle:** Independent historical table arithmetic and real reference fits when reproducible.

**Exact assertion:** Compare declared coefficients/BH/settings, preserve differences and record NOT_RUN/NOT_REPRODUCED when historical runtime/resources are unavailable.

**Negative case:** Old stored values are not evidence that R was rerun; copying them into a new result fails.

**Contract:** SM12; **owner:** R12; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V113"></a>

### V113: Scientific-change reconciliation

**Fixture:** Public audit findings plus Maintainer-approved local comparison summaries, never raw private logs.

**Oracle:** Explicit change ledger linking old versus corrected hypotheses, family/filter/mapping and response definitions.

**Exact assertion:** Explain the Negr1 auxiliary exception, 2/70 versus 3/71 and 6/70 versus 7/71 arithmetic, and intentional model/pathway/response differences without production constants.

**Negative case:** Calling the Negr1 treated-versus-control significance rescue, or principal null results positive discoveries, fails.

**Contract:** SM19; **owner:** R12; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V114"></a>

### V114: Full maintained study run

**Fixture:** Maintainer-only qualified local study inputs in a fresh output directory after all required engines are implemented.

**Oracle:** Actual maintained run manifest, independent fit checks and unchanged original hashes.

**Exact assertion:** Execute and verify the local study as new artifacts, retaining inherited unknown metadata and separation from the archived analysis.

**Negative case:** Unavailable private study or R produces NOT_RUN; synthetic data cannot substitute for this acceptance record.

**Contract:** SM25; **owner:** R12; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V115"></a>

### V115: Comprehensive usage and methods docs

**Fixture:** Final implemented capability inventory and real examples/method reports.

**Oracle:** Documentation claims compared with actual accepted handlers and tests.

**Exact assertion:** Document supported assays/designs/engines, rejection messages and exact tested commands; separate phase milestones from completed v1.0.

**Negative case:** Any unsupported method/backend or unselected license described as shipped/selected fails.

**Contract:** SM25; **owner:** R12; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V116"></a>

### V116: Spec-task-evidence reconciliation

**Fixture:** All FR-001–120/T001–120/V001–120 rows, including R10a/b ownership and unresolved evidence.

**Oracle:** Independent set/cardinality/key checks against specs/tasks and reviewed evidence.

**Exact assertion:** Retain exactly 120 unique requirement/task/acceptance identities, reviewed evidence or explicit unresolved status, and no silently retired ID.

**Negative case:** A rename/merge/split without redirect and ADR, or a pending path represented as evidence, fails.

**Contract:** SM25; **owner:** R12; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V117"></a>

### V117: Source release hygiene

**Fixture:** Source export whitelist plus deliberately seeded dummy secret/private-path files in a scratch test directory.

**Oracle:** Independent file inventory and exact whitelist/protected-source comparison.

**Exact assertion:** Exclude private data, credentials, environments and private execution logs; retain only permitted resource content/terms and explicitly unresolved ownership/license.

**Negative case:** A fake open-source license, private workbook or unapproved artifact in the package blocks public release.

**Contract:** SM14; **owner:** R12; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V118"></a>

### V118: Versioned successor export

**Fixture:** Verified source/spec/test/docs tree and an explicitly selected fresh release-root destination.

**Oracle:** Independent full file manifest/checksums and unchanged historical/private original hashes.

**Exact assertion:** Export only a new versioned successor directory after gates and permission; never overwrite the recovered archive package.

**Negative case:** Existing destination collision, absent release authorization or unresolved required license/disclosure gate blocks public export/publication.

**Contract:** SM25; **owner:** R12; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V119"></a>

### V119: Final integrated acceptance run

**Fixture:** Final reviewed tree in a clean environment with all eligible real core/golden/integration/calibration gates.

**Oracle:** Independent rerun of documented commands and final acceptance ledger.

**Exact assertion:** Re-execute release gates on that exact tree, including R10a compatibility, and report external limits separately without calling missing core work complete.

**Negative case:** A prior-commit green result or a Phase 1 suite alone cannot support v1.0 completion.

**Contract:** SM25; **owner:** R12; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V120"></a>

### V120: Review-and-handoff completion

**Fixture:** Actual final staged/unstaged/untracked diff, verified tree/commit identities and delivery manifest.

**Oracle:** Maintainer review of all changes, evidence and authorization records.

**Exact assertion:** Close only verified requested scope, record real delivery paths and remaining blockers, and perform no unauthorized commit/push/license/private-data action.

**Negative case:** A worker completion receipt or an uncreated delivery path cannot close the roadmap.

**Contract:** SM25; **owner:** R12; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

## Boundary

No recovered source/audit changes, private-data transfer, invented license, fake backend or unsupported completion claim. A missing prerequisite is NOT_RUN, not a successful acceptance case. Only the Maintainer updates traceability after independent gate execution.

**Private regression:** V111 and V114 are Maintainer-only. The original archive and private logs never enter a worker context. Any public summary requires explicit disclosure review; missing private access remains NOT_RUN.
