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

## Amendment A-2026-10-01-44 (ADR 0011, packet R15a; 2026-10-13)

Recorded by the implementing route under the Maintainer scope decision of 2026-10-08 (ADR 0011, D-66). Pending Maintainer confirmation.

- `welch` is added to the `primary_engine` enumeration of `analysis.schema.json`. It is accepted only when `methods_extension.welch.enabled` is true and the single primary model has `engine: welch`; otherwise `E_WELCH_NOT_DECLARED` is raised.
- The default remains `limma`. The single-primary rule is unchanged: exactly one model has `role: primary`, and a second primary is refused by the schema's single-item rule before any engine is considered.
- A Welch result carries its own family and is never presented as a limma result. A failed or missing Welch result is never replaced by limma under the Welch name (ADR 0005).
- The Welch engine is Phase 5 only (`runtime.phase` 5). Phase 1 to 4 configurations cannot enable `methods_extension`.
