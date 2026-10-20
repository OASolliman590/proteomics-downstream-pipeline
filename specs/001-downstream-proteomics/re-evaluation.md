# Local specification reconciliation

Date: 2026-10-18. Base: `claude/methods-016` at `4b3fb9f4c0c00b34290e2bbee55e40a23cabd931`, with pre-existing uncommitted R15 evidence and ledger edits. This review does not certify the dirty working tree as that commit. No scientific requirements, frozen receipts, task checkboxes or acceptance records were changed.

## Evidence and scope

The file-based Spec Kit comprises `.specify/memory/constitution.md`, slices 001–017, shared contracts, plans, tasks, research, quickstart, ADRs, ownership and traceability. No installed Spec Kit command definitions or templates were found in the inspected `.specify`, `.claude` and `.codex` locations. No installer or generator was run.

Executed from the repository root with bytecode disabled:

- `PYTHONDONTWRITEBYTECODE=1 .venv/bin/python scripts/check_spec_kit.py`: exit 0, PASS; 222 FR/T/V identities, 16 slices, 73 Markdown files at initial check. This checks structure, links and protected-source hashes, not scientific correctness.
- `PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest tests/regression/test_traceability_reconciliation.py -q -p no:cacheprovider`: exit 0, 11 passed. This checks identity, status, checkbox and evidence reconciliation and negative mutations.
- An independent read-only SHA-256 comparison checked all 1,042 artifact references in the acceptance ledger: zero missing or mismatched hashes. This does not independently reproduce their numerical results.

The ledger and traceability agree: 207 PASS, V180/V181 FAIL, 13 NOT_RUN. R15c's receipt explicitly says planted-signal gate confirmation is pending; retained PASS records must not be described as independent approval of that policy. R15d's framework has recorded synthetic tests, but real-study slots remain NOT_RUN. Phase 6 has recorded PASS evidence, contradicting its former unimplemented status prose. No release is authorized by those records.

## Corrections made

Updated only status prose in PHASES, packet-index, progress, quickstart and slice 016/017 spec introductions. Preserved the frozen scientific contract, FR/T/V identifiers, thresholds, historical freeze receipt, task checkboxes, source, tests and pre-existing dirty files. Historical progress entries remain inspectable; the dated status block takes precedence only for current status.

## Decisions and defects still open

1. FR-122/SM27 explicitly require `< floor` for a Monte Carlo p-value at its floor. The current code follows this rule; reporting a numeric estimate of 0.0001 as `<0.0001` needs a proposed scientific/reporting amendment, not an undocumented code fix. Preserve exact-enumeration semantics. Decision and synchronized test/report changes remain pending.
2. FR-126/FR-129 and the R13 output contract use centroid terminology, while `betadisper` defaults to spatial medians. Decide precise labels without silently changing the dispersion algorithm or renaming artifact columns.
3. V180/V181 calibration failures remain FAIL. Diagnose the stored outputs, pre-registration and implementation before any rerun; do not relax gates to manufacture acceptance. Resolve proposed planted-signal gate approval separately.
4. `workflow.py` registers every file under the figure catalogue through `rglob`. Metadata sidecars can therefore be treated as figure artifacts on filesystems that create them. Prepare a synthetic sidecar regression and bounded fix after implementation authorization; artifact integrity alone does not establish artifact type correctness.
5. Constitution frozen/not-yet-frozen wording and obsolete scope restrictions coexist with later accepted ADR amendments. The governance owner must reconcile precedence explicitly; this review did not amend the constitution or freeze manifest.

## Orchestration and ownership

This review is the sole writer of the seven documentation paths listed above plus this record (six existing documents and one new record). Pre-existing modified traceability, tasks, acceptance, reason-code and calibration files and untracked R15 receipts/evidence remain owned by their earlier coordinator; their live ownership was not independently established. No worktree was created or pruned.

AGENTS names Sol/Astra/Luna, while accepted ADR 0007 records the operator-selected Claude Opus route for the covered packets; later R15 receipts record Haiku builder/Opus review but are not themselves authorization for future substitution. Keep those histories distinct. Do not infer a current model change from this review. Native Codex agents and external CLI agents are distinct harnesses.

Local readiness on 2026-10-18: Codex CLI 0.163.0-alpha.10.2 is authenticated with ChatGPT; Claude Code 2.1.299 reports `loggedIn: false`. The global fleet has feature/review-debate `gpt-6.1-sol`, docs `gpt-5.6-terra`, review-main Grok `grok-4.6`, complex Claude and UI agy. Global lanes do not override project policy. Grok's exhausted balance is dispatch-supplied evidence and was not retested. No paid/model-generation readiness probe ran; authentication and cached model listings do not prove available quota. No reviewer or implementer was dispatched, and no cross-model review is claimed.

## Prioritized handoff and reviewer brief

Next authorized work is documentation/evidence review only. Dependency order: confirm scope and ownership → independent documentation review → Maintainer decisions on reporting/provenance/governance → separately authorized bounded repairs → relevant verification → final review. Calibration and real-study execution require a separate decision; neither is part of this handoff's authorization.

Reviewer: read AGENTS, ADR 0007 and the applicable later authorization, this record, the exact diff of the six changed existing documents and this file, original requirements and check outputs. Review only specification/source/test/oracle material; do not receive private study data, logs or derived results. Return reproducible findings with file/line, trigger, impact, confidence and blocking/nonblocking disposition. Do not edit, close gates, freeze, merge or expand scope. Check that status changes do not imply independent scientific acceptance. Review remains NOT_RUN because the applicable role/route needs confirmation and the recorded Claude route is unauthenticated. A Codex substitution requires authorization; no silent fallback.

Future repair briefs must pin the then-current base and dirty-file manifest; this base alone is insufficient. Require a short understanding check before coding, a single owned allowlist or isolated checkout, explicit dependencies, acceptance criteria, non-goals and exact commands/results in the return. After two unsuccessful repair rounds, stop for re-scoping. Suggested bounded briefs are: (a) R13 display semantics after amendment, owning only its formatter/report tests; (b) figure-catalogue sidecar exclusion, owning enumeration and a synthetic boundary test; (c) calibration diagnosis, read-only existing synthetic evidence until a new run is authorized. No overlapping writers, merges or automatic dispatch.
