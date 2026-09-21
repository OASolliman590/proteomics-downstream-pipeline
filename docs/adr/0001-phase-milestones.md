# ADR 0001 — Phase milestones and immutable evidence

Status: proposed for freeze; date 2026-09-20.

## Context

The draft treated twelve large packets as one release and contained conflicting approval language. R10 shared files across many packets.

## Decision

Use the four phases in PHASES.md. Split dispatch into R10a (FR/T/V091–094) and R10b (095–100), retaining one R10 slice and all 120 identities; no redirect is required because no requirement was merged/split/retired. All packet write ownership is disjoint via fixed file interfaces. Review/freeze remain outstanding.

## Consequences

Phase 1 is a named milestone, never v1.0 completion. roadmap.json describes slice completion; packet ownership describes subpacket dispatch. The requested hardening report path under immutable docs/audit conflicted with the no-audit-diff constraint: use docs/SPEC_HARDENING_REPORT.md, leaving the audit tree untouched. The new draft config field layout is explicit in semantic-validation; no maintained old-schema consumer exists.

## Verification

V001–V050/V091–V094 plus phase/ownership/identity checks; final full release still requires all 120.
