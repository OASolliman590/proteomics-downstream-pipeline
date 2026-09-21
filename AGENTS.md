# Proteomics downstream orchestration

Repository: OASolliman590/proteomics-downstream-pipeline
EPIC: specs/001-downstream-proteomics

## Roles and routing

Operator-selected routes:
- Parent: GPT-5.6 Sol high.
- Independent audit: GPT-6 Astra low.
- Packet implementation: GPT-5.6 Luna high.

Sol high owns the kit, freeze, dispatch, independent verification and closure.
Astra low audits the candidate or a bounded risky diff, in fresh context, then stops.
Luna high implements exactly one frozen packet in a fresh worker thread.

Use the model routes selected by the operator. Do not substitute another route
silently or claim a prompt changed the model. No max/ultra/fast escalation.
Threads: `sol / proteomics-kit`, `astra / spec-freeze`,
`luna / <packet-id>`, `astra / <packet-id>-diff`.
No inherited worker conversation, thread forks, recursive delegation or automatic
spawning. Sol writes starters; the operator opens the intended threads.

## Source of truth

Read the live constitution, EPIC/PHASES.md, packet-index.md,
packet-ownership.json and relevant shared contracts/slice files.
Scientific rules live in contracts/scientific-methods.md, not this routing file.
Any contradiction returns to Sol before dependent work continues.
An unresolved Sol/Astra disagreement on a requirement goes to the user;
neither reviewer silently overrides a product decision.
A downloaded overlay, successful checker or worker receipt is not a freeze.

## Scope and ownership

1. No Luna implementation before independent audit, recorded Maintainer freeze
   and explicit GO for the named packet/group. R01 is always first.
2. Preserve FR-001–FR-120, T001–T120 and V001–V120. No scope expansion,
   weakened eligibility, silent renumbering or placeholder scientific backend.
3. Never edit pipeline/**, legacy/**, docs/audit/**, recovered provenance/hashes
   or original private evidence. Never request/upload private study inputs.
4. Luna writes only its exact packet-ownership allowlist. Everything else,
   including AGENTS.md, .specify/** and specs/**, is read-only for Luna.
   Shared-interface repairs return to Sol; they are not allowlist exceptions.
5. Sol coordinates specifications and verification, not unassigned engine code.
   High risk requires review; it does not prohibit a bounded Luna implementation.
6. No commits, pushes, publication, license selection, destructive cleanup,
   global installations or remote-setting changes without explicit authorization.
   Protect existing edits; never reset/stash them merely to satisfy a preflight.

## Scheduling

Phase 1: R01 -> R02 -> R03 -> R04 -> R05 -> R10a.
Phase 2: R06 || R07 || R09, then a whole-group barrier, then R08.
Phase 3: R10b -> R11 -> R12.

The parallel group starts only after accepted R05 AND the Phase 1/R10a gate.
Only that group may overlap. R08 waits for the entire group, even though its
scientific dependencies are R06/R07. No concurrent reporting/documentation job.
Sol prepares isolated workers from the same verified prerequisite snapshot;
separate chats alone do not establish filesystem isolation. Worker outputs,
project-local environments and build directories must not collide.

## Evidence and closure

A receipt means ready for review, never accepted. Sol inspects actual diffs,
reruns relevant gates and records evidence before advancing dependencies.
Keep logs in packet-local evidence files; no transcript dumps into Sol.
Astra reads specifications, source/test diffs and oracle definitions, not test
logs, private evidence or private-derived contents. Astra never implements,
marks completion, authorizes dispatch or closes the loop.

PASS, FAIL, NOT_RUN, SKIPPED and INAPPLICABLE remain distinct. An unavailable
runtime is not PASS; absent software is not scientific inapplicability.
Do not invent a verified_commit for uncommitted work: identify the actual
working-source manifest/hash and leave verified_commit null.

Phase 1 is not full R10 or v1.0. R10a covers V091–V094; R10b covers V095–V100.
R11 owns release locks/calibration. R12 private gates are Maintainer-only.
Stop at the authorized packet/phase boundary. Sol alone closes that scope.
