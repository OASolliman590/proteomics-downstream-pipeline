# R10b Maintainer-route receipt — Full offline report and user workflow (R10b)

Packet R10b (FR-095–FR-100 / T095–T100 / V095–V100) was implemented and verified on 2026-10-01 (all gates re-executed 2026-10-02 after the Phase 2/3 amendments) by the **Claude (Opus) route**, which the operator authorized to act for the Sol/Astra/Luna roles ([ADR 0007](../../../adr/0007-claude-route-and-integration-amendments.md)). This is not work by the routes named in AGENTS.md.

## Identities

- Contract versions: kit 1.2.0 with the 2026-10-01 amendments (ADR 0006/0007); scientific methods contract v1.2.0.
- Working-source manifest (packet allowlist, 11 files): `61f2a927313893d4965470c6bcacb932b401171940bec047c97b399072e18334` — lines in `evidence/working-source-manifest.txt` (`sha256  path`, sorted, LF).
- `verified_commit`: **null** — nothing is committed (branch `claude/full-pipeline`, uncommitted working tree).
- Environment: macOS-26.7-x86_64-i386-64bit-Mach-O; Python 3.13.15; R 4.6.1 with limma 3.68.5, statmod 1.5.2, impute 1.86.0, vegan 2.7-6, permute 0.9-10, DEqMS, proDA and fgsea at the versions pinned in `renv.lock`, in a project-local library (`.r-lib/`); `LANG=en_US.UTF-8`.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/integration/test_reports.py tests/integration/test_report_stub.py -q -rs -p no:cacheprovider` | 0 | 16 passed in 110.60s (0:01:50) | `python-tests-1.log` |

Phase-level regression (R01 suite, historical unittest, all R tests, baseline registry, spec checker, whole Python suite) is in [`docs/validation/phase-gates`](../../phase-gates/README.md).

## Acceptance cases

| Case | Status | Tests | Positive oracle | Negative case |
|---|---|---|---|---|
| V095 | PASS | test_reports.py::test_v095_* | complete zero_null table rendered with all 24 planned rows incl. secondary and excluded rows; family totals from families.tsv | no display-filtered counts |
| V096 | PASS | test_reports.py::test_v096_* | null-type labels, exploratory fgsea/ORA, descriptive response without P columns, PERMANOVA view | no confirmation/rescue/mechanism wording |
| V097 | PASS | test_reports.py::test_v097_* | SVG points equal source-table rows; coordinates equal DEA effect and -log10 P; valid PDF/PNG/SVG headers | no stars, no extra points |
| V098 | PASS | test_reports.py::test_v098_* | methods differ by executed hypothesis; unknown tissue explicit | unexecuted engines never mentioned |
| V099 | PASS | test_reports.py::test_v099_* | QC-only, cancelled and crashed runs render truthfully | no COMPLETED/zero substitution |
| V100 | PASS | test_reports.py::test_v100_* | the seven documented commands run verbatim (validate, plan, run, verify, report, run Phase 2 demo, compare) | only implemented commands documented |


## Notes

- Decisions taken on ambiguous points are D-01–D-30 in [decisions.md](../../../../specs/001-downstream-proteomics/decisions.md).
- Independent audit (Astra-equivalent): **NOT_RUN**. Private-study regression: **NOT_RUN** (Maintainer-only, out of scope). No private data or private-derived numbers were used.

## Acceptance

All listed gates exited 0 and every acceptance case above, including its negative case, passed on re-execution. R10b is recorded as accepted by the operator-authorized Claude route (self-verified, uncommitted); independent review remains open.
