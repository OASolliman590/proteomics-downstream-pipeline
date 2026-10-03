# Independent audit of fd8dfa9: fixes and fail-before/pass-after evidence (2026-10-02)

Verdict of the audit: **ACCEPT WITH FIXES**. Decisions are in [ADR 0008](../../adr/0008-independent-audit-fixes.md) and D-31 to D-38.

How the evidence was produced. `git archive fd8dfa9` was extracted to a scratch directory, and its own `proteomicsCore` was built into a separate library placed first on `R_LIBS_USER`. The new and amended test files were copied in, and the tests were run with `PYTHONPATH` pointing at that tree. The result is [`evidence/before.log`](evidence/before.log). The same tests were then run on the fixed working tree: [`evidence/after.log`](evidence/after.log). Nothing was stashed, reset or committed.

| Finding | New or amended tests | fd8dfa9 | fixed |
|---|---|---|---|
| MAJOR 1: primary model on an imputed matrix | `test_limma.py::test_audit_m1_*` (5) | FAIL (5) | PASS |
| MAJOR 2: degenerate subject-blocked PERMANOVA | `test-permanova.R` "audit M2" (3); `test_permanova.py::test_audit_m2_subject_blocked_*[between]`, `test_audit_m2_mixed_*` | FAIL | PASS |
| MAJOR 2: within-subject case | `test_permanova.py::test_audit_m2_subject_blocked_*[within]` | PASS (the old within-subject scheme was already valid for this design) | PASS |
| m1: calibration fit path | `test_calibration_evidence.py::test_v105_calibration_uses_the_production_fit_and_eligibility_rule` | FAIL | PASS |
| m2: weak oracles | `test-permanova.R` "audit m2: covariate …", V124 (amended) | FAIL | PASS |
| m2: exact labels | `test-permanova.R` "audit m2: exact interpretation labels", `test_permanova.py::test_audit_m2_exact_location_label_*` | PASS (the old labelling was correct; these tests strengthen the oracle) | PASS |
| m3: silent pair/group drops | `test-permanova.R` "audit m3" | FAIL | PASS |
| m4: interaction refusal | `test-permanova.R` "audit m4", V125 (amended); `test_permanova.py::test_audit_m4_*`, `test_v125_*` | FAIL | PASS |
| m5: contract and R hash | `test-permanova.R` "audit m5"; `test_permanova.py::test_audit_m5_*` | FAIL | PASS |
| m6: covariate names | `test-permanova.R` "audit m6"; `test_permanova.py::test_audit_m6_*` | FAIL | PASS |
| R03/R04 vectorization (performance, behaviour-preserving) | `test-preprocessing.R`/`test-design.R` "audit: vectorized …" (the old loop kept as the oracle) | PASS (regression guard) | PASS |

The re-run calibration (production path, 1000 datasets per scenario) and the re-run 20000×100 benchmark are recorded in the [R11 receipt](../012-validation/receipt.md).
