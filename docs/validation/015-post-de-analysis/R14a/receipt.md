# R14a Maintainer-route receipt — Post-DE result structure and concordance (operator-authorized amendment, ADR 0009)

Packet R14a (FR-131–FR-136 / T131–T136 / V131–V136) was implemented and verified on 2026-10-03 by the **Claude (Opus) route**, which the operator authorized to act for the Sol/Astra/Luna roles ([ADR 0007](../../../adr/0007-claude-route-and-integration-amendments.md)). This is not work by the routes named in AGENTS.md. Independent audit: **NOT_RUN**.

## Identities

- Contract versions: kit 1.3.0 with the ADR 0009 amendment; scientific methods contract v1.3.0 (SM31, SM41); `specs/015-post-de-analysis/contracts/post-de.md`.
- Working-source manifest (packet allowlist, 8 files): `a9dc49634dc1a87a2ce0ae42f5fdc35ec4b839292152bfcf0669b05706603662` — `evidence/working-source-manifest.txt` (`sha256  path`, sorted, LF).
- Shared-interface amendment made first: A-2026-10-01-14 (capability maps, optional `post_de` schema block and Phase 4, design-service re-plan API, workflow integration, report providers, R dependencies).
- `verified_commit`: **null** (local commit on `claude/post-de`, not pushed).
- Environment: macOS x86_64; Python 3.13.15; R 4.6.1 with limma 3.68.5 in the project-local `.r-lib/`; `LANG=en_US.UTF-8`.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/scientific/test_post_de_sets.py -q -rs -p no:cacheprovider` | 0 | 11 passed | `evidence/python-tests-1.log` |
| `Rscript --vanilla -e "testthat::test_file('r/proteomicsCore/tests/testthat/test-post-de-sets.R', ...)"` | 0 | 146 expectations, 0 failed, 0 skipped | `evidence/r-test-post-de-sets.log` |

Before the commit the whole Python suite (287 passed, 0 skipped; `evidence/full-python-suite.log`), all R testthat files (502 expectations, 0 failed, 0 skipped; `evidence/r-all-testthat.log`) and the spec kit (PASS; `evidence/spec-kit-check.log`) were re-run on this tree.

## Acceptance cases

| Case | Status | Positive oracle | Negative case |
|---|---|---|---|
| V131 | PASS | membership of 10 declared sets (leaves, union, intersect, difference, complement within tested, direction, raw P, TREAT) equals Python set operations on `dea/zero_null.tsv`/`treat.tsv` with q < 0.05; one membership row per tested feature; raw-P set labelled `exploratory_raw_p`; R-level oracle with base R `union`/`intersect`/`setdiff` | a free-form string, an unknown operator and an unknown node key fail `E_SETRULE_GRAMMAR` before any stage runs |
| V132 | PASS | every set row carries rule (JSON), criteria, thresholds, family IDs, endpoints, null type and input matrix | mixing zero-null and TREAT `E_SETRULE_MIXED_NULL`; an undefined contrast `E_SETRULE_UNKNOWN_CONTRAST`; primary- and sensitivity-matrix families `E_SETRULE_MIXED_MATRIX` |
| V133 | PASS | regions equal brute-force enumeration of the 2^k − 1 patterns for k = 5 (pipeline) and k = 3, 5 (R); counts sum to the union; protein-group and unique-gene counts differ where two protein groups share a gene; Venn and UpSet with source tables | a Venn for k = 4 is refused `E_VENN_K` (planner sub-analysis record and `refusals.tsv`) while the UpSet is produced |
| V134 | PASS | concordant-up, concordant-down, discordant and one-only memberships equal the hand-labelled planted features | discordant features are never counted as concordant overlap |
| V135 | PASS | Pearson, Spearman and origin slope equal R `cor()`/`lm(y ~ 0 + x)` within 1e-10; sign agreement equals the oracle; label `descriptive` | no P column is emitted |
| V136 | PASS | the unit-disjoint pair (B-A, D-C) reports the `phyper` upper tail with the universe tested in both (BH within the descriptive-overlap family) | the shared-control pair (B-A, C-A) has no P row and reason `shared_biological_units` |

## Notes

- Sets use the strict criterion q < threshold (and p < threshold for `exploratory_raw_p`), as V131 states; family rejection counts elsewhere use q ≤ cutoff (D-45).
- Overlap and direction-aware overlap are computed for pairs of single-endpoint sets; a composite set has no single effect direction (D-45).
- No private data or private-derived numbers were used; all fixtures are synthetic (`tests/fixtures/post_de/`).

## Acceptance

Both gates exited 0 and every acceptance case above, including its negative case, passed. R14a is recorded as accepted by the operator-authorized Claude route (self-verified; independent audit NOT_RUN).
