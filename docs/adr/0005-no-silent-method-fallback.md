# ADR 0005 — No silent method fallback

Status: proposed for freeze; date 2026-09-20.

## Context

An optional engine's absence or failure could be concealed by substituting limma while retaining another method label or by calling installation failure inapplicability.

## Decision

Keep requested engine and result identity fixed. Scientific ineligibility, missing installation, unimplemented capability and numerical failure get distinct typed reasons and actual states. A separately prespecified model may run under its own name, never as a substitute.

## Consequences

Required failure stops dependent work; optional operational failure yields visible partial execution. Reports derive states from artifacts, not hardcoded PASS/biology. No new placeholder backend or automatic install is authorized.

## Verification

V003, V006, V038, V059, V080, V091–V099.
