# Implementation plan — R13

**Phase:** 2. **Status:** operator-authorized amendment (2026-10-01); see [ADR 0006](../../docs/adr/0006-permanova-scope-amendment.md).

Read [spec.md](spec.md), [tasks.md](tasks.md), [contracts/permanova.md](contracts/permanova.md), SM26–SM30 in the [scientific contract](../001-downstream-proteomics/contracts/scientific-methods.md#multivariate-permanova) and the [ownership record](../001-downstream-proteomics/packet-ownership.json).

Dispatch: after accepted R05 and R10a, serially, before the R06 ∥ R07 ∥ R09 group (operator order). R13 adds the `permanova` capability through the fixed maps (Maintainer amendment A-2026-10-01-01 in `runtime.py` and `dispatch.R`), its own Python service, R handler, report-section provider, fixtures, tests and evidence. It consumes plan artifacts read-only and verifies their hashes; DEP-derived sets read the completed R05 table and verify its hash.

Use each V case's independent oracle (direct vegan/base-R calls written in the test) before trusting the adapter. Record command, versions, seeds, hashes and exit codes in packet-local evidence. Missing vegan/permute is NOT_RUN.
