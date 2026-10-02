# Decisions and open questions — 2026-10-01 continuation (Claude route)

Recorded by the Claude (Opus) route acting for the Maintainer role under the operator's authorization ([ADR 0007](../../docs/adr/0007-claude-route-and-integration-amendments.md)). Where a requirement was ambiguous the conservative reading was chosen; nothing here weakens an eligibility rule.

## Decisions

| ID | Topic | Decision (conservative reading) |
|---|---|---|
| D-01 | Route | Claude (Opus) filled Sol/Astra/Luna roles; no separate independent audit was run (`independent_audit: NOT_RUN` in every receipt). |
| D-02 | Run layout | Each stage is promoted into its semantic directory (`inputs/`, `preprocessing/`, `designs/`, `dea/`, `permanova/`, `report/`) with its `stage-result.json` inside, because R01 `verify` requires outputs to sit beside their result. No separate `stages/<id>/` tree. Failed stages keep their partial outputs under `logs/<stage>-failed/`. |
| D-03 | Planning stage hash | intake/preprocessing/design run before the plan exists and carry `plan_hash = null`; the plan then freezes their artifact hashes (A-2026-10-01-02). |
| D-04 | Observed-data fit | Primary limma fits only genuinely observed cells; prior-imputed values are set to NA for the primary model ("observed-data limma adapter", FR-041). Sensitivity models use their declared sensitivity matrices. |
| D-05 | Interaction grammar | Factors entering an interaction use reference-omitted columns (first declared level / declared reference), giving full-rank conventional designs; aliasing is still checked by QR. |
| D-06 | Contrast direction | For a two-group contrast, `required_groups[0]` is the numerator; weights whose implied group-mean weights contradict that order fail `E_CONTRAST_DIRECTION` (V035 negative). |
| D-07 | Omnibus | Moderated F uses limma's convention df2 = df.prior + df.residual (matches `eBayes`/`classifyTestsF`); reduced-design omnibus tests are refused (`E_OMNIBUS_REDUCED_UNSUPPORTED`). |
| D-08 | Numerical failure | Features whose fit has non-finite residual variance or coefficients are removed from the moderation universe and exported as `numerical_failure` (family incomplete), because limma would otherwise report t = 0, P = 1. All-feature failure is `E_ENGINE_FAILED`. |
| D-09 | Detection on paired designs | Requested detection on a blocked design is documented as INAPPLICABLE (`E_DETECTION_DESIGN_UNSUPPORTED`) inside preprocessing with no P column, rather than rejecting the whole plan (detection is exploratory). |
| D-10 | Vendor profiles | The four vendor mappings are header-set profiles qualified only on synthetic exports; they do not claim compatibility with any vendor release. `input.features = "from_mapping"` marks mapped/legacy inputs. |
| D-11 | Legacy importer | Only `load_workbook`, `sample_columns`, `canonical_sample`, `SHEET_CONTRASTS` and `GROUP_MAP` of the recovered parser are reused; its `validate_and_reconstruct` hardcodes the historical study size and is not called. Method B names follow docs/PIPELINE.md. |
| D-12 | Serialization | Numbers are written with exact round-trip precision; infinite diagnostic values (for example an infinite empirical-Bayes prior df) are written literally as `Inf`; NaN is never written. |
| D-13 | Report requiredness | The R10a report stage is optional: its failure makes a run PARTIAL but never upgrades or hides an analysis state. |
| D-14 | PERMANOVA feature universe | Genuinely observed, complete and non-constant features within each test's observations; z-scoring within each test (as in the reference workflow); no imputation. |
| D-15 | PERMANOVA phase | R13 is Phase 2 (v0.2) and dispatched after R05+R10a; Phase 1 configurations cannot enable it. |
| D-16 | PERMANOVA production resolution | Production profile requires ≥ 999 permutations (default 9,999); smoke-test profiles may use fewer. The actual number of permutations used (which permute may reduce under complete enumeration) defines the floor. |
| D-17 | Per-feature contribution | Only univariate R² and group means are exported; no per-feature P values (protein-wise inference belongs to R05). |
| D-18 | Optional PERMANOVA ineligibility | An empty interaction cell or covariates in a subject-blocked design make an optional PERMANOVA INAPPLICABLE with the typed code, and a required one a plan rejection (exit 2). |
| D-19 | Figures | Base R graphics; PNG always, PDF when requested, SVG only when svglite is installed (cairo SVG is unavailable on this Mac); otherwise the SVG is recorded NOT_RUN. |
| D-20 | Environment | R tests need a UTF-8 locale (`LANG=en_US.UTF-8`) on macOS; the R01 suite fails under the C locale with non-ASCII paths. Verified on macOS x86_64, R 4.6.1, Python 3.13.15. |
| D-21 | proDA reproducibility | proDA uses random numbers while fitting; the adapter calls `set.seed(seed, "L'Ecuyer-CMRG")` with the plan seed before `proDA()` and records seed and RNG kind. |
| D-22 | Engine eligibility timing | Adapter eligibility (for example proDA on TMT or prior-imputed data) is checked declaratively at plan time, before any data stage runs, so an ineligible optional engine is INAPPLICABLE and never causes a preprocessing error. |
| D-23 | Pathway gene model when the primary is not limma | Gene-level pathway tests always use a limma model; if the primary engine is DEqMS/proDA, the named linear sensitivity model is used and every pathway row is labelled sensitivity (V078). |
| D-24 | Score randomization phase | The R05 guard against score testing applies only to Phase 1 plans; Phase 2 plans may run the R09 score randomization. |
| D-25 | Report figures | The full report draws SVG in Python directly from the saved figure-source tables (exact coordinates) and PNG/PDF through R base graphics from the same tables; a format whose renderer is unavailable is listed NOT_RUN on the figure. Supersedes D-19 for R10b. |
| D-26 | Cache scope | A stage is reused only when its fingerprint (normalized request, full code manifest, environment) matches and its outputs verify; any code change therefore invalidates every stage (conservative). |
| D-27 | Validation ledger tree | Uncommitted Claude-route records use the packet working-source manifest SHA-256 as `reviewed_tree`; `verified_commit` stays null. |
| D-28 | Report visual inspection | Automated V109 checks re-derive every plotted coordinate; the implementer's own visual inspection is recorded with its findings (overplotted identical points are not marked; SVG points at the range extremes touch the axes; the SVG has no zero reference line while the R PNG does). The independent reviewer's visual inspection is NOT_RUN. |
| D-29 | CI | `.github/workflows/maintained.yml` defines Linux and Windows smoke jobs plus a dispatch-only release-calibration job. It has never run (no push authorized), so V107 is NOT_RUN. |
| D-30 | TREAT calibration | Calibration covers the zero-null core and CAMERA/ROAST; TREAT error control and the response/PERMANOVA tests are not separately calibrated in R11 (their acceptance relies on oracle tests). |

## Open questions for Omar

1. **p display at the floor.** Per your instruction a P at 1/(nperm+1) is shown as "< 1/(nperm+1)". Strictly the estimate equals the floor ("≤"). Keep "<", or switch to "≤"?
2. **Independent audit.** Do you want an independent reviewer (another model or a person) to audit R02–R11 and R13 (including the V109 visual inspection) before you treat them as accepted?
3. **Commit plan.** Nothing is committed. Do you want one commit per packet (as the repository history suggests) when you authorize commits?
4. **`.r-lib/` is not git-ignored.** `scripts/maintained/bootstrap.sh` uses it, but `.gitignore` (R12-owned) only ignores `/.Rlib/`. Should R12 add `/.r-lib/`?
5. **PERMANOVA milestone.** R13 is placed in Phase 2 (v0.2). Should it instead be released as a separately labelled extension?
6. **Sex/sample key.** The private reanalysis that motivated PERMANOVA still has an unconfirmed sex key; nothing from it is in this repository, but any future private regression of PERMANOVA should wait for that confirmation.
7. **R12 tooling.** Work stopped at the R11 boundary. R12 (legacy regression tooling, release builder, release docs) is not implemented and all of V111–V120 are NOT_RUN. Its criteria reference private-study details, so the private comparison must be run by you locally. Should R12 tooling be the next packet?
8. **CI (V107).** `.github/workflows/maintained.yml` has never run. V107 can only pass after you push a branch and the Linux and Windows jobs succeed.
9. **Performance headroom.** The 20000×100 benchmark passed in 884 s, and 382 s of that was the design stage (a per-feature, per-contrast estimability loop with a large `rbind`). Do you want this optimised (an R04 amendment) before release?
10. **Report layout findings (V109).** The SVG volcano does not mark overplotted identical points, has no padding at the range extremes and has no zero reference line. Should these be fixed in R10b before an independent visual review?
