# ADR 0002 — Limma as default primary engine

Status: proposed for freeze; date 2026-09-20.

## Context

Several scientifically different methods were named without a single operational default. A convenient method switch could become hit-count selection.

## Decision

Default primary_engine=limma. Require exactly one matching primary model/hypothesis. DEqMS/proDA remain mandatory Phase 2 adapters, scientifically eligibility-gated and separate in result type/family; do not rank them by discoveries.

## Consequences

Zero-null and TREAT remain different endpoints. General contrasts use the exact SM08 route, not an approximate shortcut assumed safe. No selected package version or environment is called qualified before actual execution.

## Verification

V038, V041–V050, V051–V060.
