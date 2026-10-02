# R10a Maintainer-route receipt — Thin offline report and Phase 1 workflow (R10a)

Packet R10a (FR-091–FR-094 / T091–T094 / V091–V094) was implemented and verified on 2026-10-01 (all gates re-executed 2026-10-02 after the Phase 2/3 amendments) by the **Claude (Opus) route**, which the operator authorized to act for the Sol/Astra/Luna roles ([ADR 0007](../../../adr/0007-claude-route-and-integration-amendments.md)). This is not work by the routes named in AGENTS.md.

## Identities

- Contract versions: kit 1.2.0 with the 2026-10-01 amendments (ADR 0006/0007); scientific methods contract v1.2.0.
- Working-source manifest (packet allowlist, 8 files): `82a0918bbe29740dd0b3ebf824f812b510f81f9814fcc2bc08b2e0f5e35c6196` — lines in `evidence/working-source-manifest.txt` (`sha256  path`, sorted, LF).
- `verified_commit`: **null** — nothing is committed (branch `claude/full-pipeline`, uncommitted working tree).
- Environment: macOS-26.7-x86_64-i386-64bit-Mach-O; Python 3.13.15; R 4.6.1 with limma 3.68.5, statmod 1.5.2, impute 1.86.0, vegan 2.7-6, permute 0.9-10, DEqMS, proDA and fgsea at the versions pinned in `renv.lock`, in a project-local library (`.r-lib/`); `LANG=en_US.UTF-8`.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/integration/test_report_stub.py -q -rs -p no:cacheprovider` | 0 | 9 passed in 39.58s | `python-tests-1.log` |

Phase-level regression (R01 suite, historical unittest, all R tests, baseline registry, spec checker, whole Python suite) is in [`docs/validation/phase-gates`](../../phase-gates/README.md).

## Acceptance cases

| Case | Status | Tests | Positive oracle | Negative case |
|---|---|---|---|---|
| V091 | PASS | test_report_stub.py::test_v091_* | only verified values; zero discoveries stay zero; missing DEA stays unknown; failures stay failed | tampered artifact shown UNVERIFIED with no values; verify exits 5 |
| V092 | PASS | test_report_stub.py::test_v092_* | one CLI command runs intake→preprocessing→design→plan→limma→report; plan precedes the fit | genuine R failure gives FAILED/exit 4 with a partial report; Phase 1 score request E_PHASE_CAPABILITY |
| V093 | PASS | test_report_stub.py::test_v093_* | embedded CSS, captions, th scope, navigation; local links resolve | no script/remote/stylesheet links; no absolute local paths |
| V094 | PASS | test_report_stub.py::test_v094_* | biological vs technical n, exclusion reason, unknown provenance, exact source links | every canonical observation is analysed or explicitly excluded |


## Notes

- Decisions taken on ambiguous points are D-01–D-30 in [decisions.md](../../../../specs/001-downstream-proteomics/decisions.md).
- Independent audit (Astra-equivalent): **NOT_RUN**. Private-study regression: **NOT_RUN** (Maintainer-only, out of scope). No private data or private-derived numbers were used.

## Acceptance

All listed gates exited 0 and every acceptance case above, including its negative case, passed on re-execution. R10a is recorded as accepted by the operator-authorized Claude route (self-verified, uncommitted); independent review remains open.
