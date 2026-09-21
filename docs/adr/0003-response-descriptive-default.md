# ADR 0003 — Descriptive treatment response by default

Status: proposed for freeze; date 2026-09-20.

## Context

In-sample selected disease axes, shared untreated controls, unlimited reversal labels and missing equivalence could support invalid rescue claims.

## Decision

Default response_mode=descriptive_only and score_test=off. Freeze vocabulary/RI bounds in SM19–SM24. Separate descriptive output with no P/q fields from eligible TOST, independently directed conjunction and independently frozen scores. Phase 1 has no response-score module.

## Consequences

RI>1 crosses control; RI>1.2 is the overshoot class. Nonsignificance is not equivalence. Exact enumeration uses k/N; Monte Carlo uses its separately defined plus-one protocol. Historical fractions/Negr1 stay fixtures.

## Verification

V035, V050, V081–V090, V092.
