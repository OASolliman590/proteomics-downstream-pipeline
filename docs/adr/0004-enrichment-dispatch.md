# ADR 0004 — Design-compatible enrichment dispatch

Status: proposed for freeze; date 2026-09-20.

## Context

A single pathway label concealed competitive, self-contained and preranked nulls; paired CAMERA could ignore covariance and ortholog projection could be mislabeled.

## Decision

Use SM16 dispatch: independent finite linear→CAMERA; paired/repeated→ROAST with correct fixed or random block; fgsea exploratory gene-set null; ORA measured mapped universe with all eligible sets adjusted before foreground-zero display removal. Use SM14–SM15 hashed resources and label-independent mapping.

## Consequences

No paired fixed-subject CAMERA exception in this version. Alternative-engine linear pathways are separate sensitivities, not equivalent tests of the native model. Human→rat collections retain projection evidence.

## Verification

V061–V080 and V096.
