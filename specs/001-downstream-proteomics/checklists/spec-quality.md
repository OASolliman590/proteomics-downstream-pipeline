# Spec quality checklist — v1.2.0 frozen

Evidence for this hardening pass is in the external bundle verification record and [hardening report](../../../docs/SPEC_HARDENING_REPORT.md), not in a fabricated implementation-evidence path. A document/schema check cannot close a software V case.

- [x] Check phase language consistently distinguishes baseline, v0.1/R10a, v0.2 and full v1.0/R10b.
- [x] Check all packet write allowlists are disjoint and only R06/R07/R09 parallelism is permitted after accepted R05 and R10a; R08 explicitly waits for all three group members.
- [x] Parse every changed JSON; validate Draft 2020-12 schemas and all three examples, with negative cases.
- [x] Retain exactly FR-001–120/T001–120/V001–120; each case has fixture, independent oracle, exact assertion and negative case; evidence stays empty/pending.
- [x] Resolve local links in START_HERE, roadmap, packet-index and README, identifying unchanged baseline-only targets separately from local overlay targets.
- [x] Verify every historical source/audit byte and registered source SHA-256 on the actual checkout; do not regenerate provenance from changed files.
- [x] Run the original repository scripts/check_spec_kit.py on the complete checkout/overlay before acceptance.
- [x] Record the initial independent Astra audit against candidate manifest `e65317b7873d63361399fe8b55f7b90c893f0ae92b5193ff119f6c642a9ee90d` and disposition its three P0 and two P1 findings.
- [x] Complete the fresh independent re-audit of the corrected shared contracts and record the Maintainer freeze; R01 still requires explicit GO.
- [ ] Resolve license/public-description/disclosure gates before any public successor release; do not invent a license.

No maintained engine was implemented. Current first action is independent audit and Maintainer freeze, then R01 only. The hardening report is deliberately outside docs/audit/ to honor the immutable-tree constraint.
