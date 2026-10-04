# R14e Maintainer-route receipt — Post-DE co-abundance modules and interaction networks (operator-authorized amendment, ADR 0009)

Packet R14e (FR-161–FR-165 / T161–T165 / V161–V165) was implemented and verified on 2026-10-04 by the **Claude (Opus) route**, operator-authorized for the Sol/Astra/Luna roles ([ADR 0007](../../../adr/0007-claude-route-and-integration-amendments.md)). Independent audit: **NOT_RUN**.

## Identities

- Contracts: kit 1.3.0 (ADR 0009); scientific methods v1.3.0 (SM39, SM40 with SM14); `specs/015-post-de-analysis/contracts/post-de.md`.
- Working-source manifest (packet allowlist with the shared design factory, 7 files): `b3b093c05ce92455f6cf82130e8fef3c2fe52a9a8580ee9adb54248a2053b35a` — `evidence/working-source-manifest.txt`.
- Shared changes: A-2026-10-01-18 (R07 snapshot kind `ppi`, `E_RESOURCE_DOWNLOAD`, resource inputs in the shared post-DE input check, `cut_height`). Decisions D-44 (no WGCNA package dependency) and D-49.
- `verified_commit`: **null** (local commit on `claude/post-de`, not pushed). Environment: macOS x86_64, Python 3.13.15, R 4.6.1, dynamicTreeCut 1.63.1, limma 3.68.5, `.r-lib/`, `LANG=en_US.UTF-8`.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/scientific/test_post_de_networks.py -q -rs -p no:cacheprovider` | 0 | 11 passed | `evidence/python-tests-1.log` |
| `Rscript --vanilla -e "testthat::test_file('r/proteomicsCore/tests/testthat/test-post-de-networks.R', ...)"` | 0 | 28 expectations, 0 failed, 0 skipped | `evidence/r-test-post-de-networks.log` |
| `.venv/bin/python -m pytest tests/scientific/test_limma.py -q -rs -p no:cacheprovider` (shared limma change A-2026-10-01-19) | 0 | 14 passed | `evidence/python-limma-1.log` |
| `.venv/bin/python -m pytest tests -q -p no:cacheprovider -rs` (at commit e1345cc, before the follow-up below; re-run with R14f) | 0 | 339 passed, 0 skipped | `evidence/full-python-suite.log` |
| all R testthat files (`testthat::test_dir(..., load_package = "installed")`) | 0 | 677 expectations, 0 failed, 0 skipped (after the follow-up) | `evidence/r-all-testthat.log` |
| `.venv/bin/python scripts/check_spec_kit.py` | 0 | PASS | `evidence/spec-kit-check.log` |

## Acceptance cases

| Case | Status | Positive oracle | Negative case |
|---|---|---|---|
| V161 | PASS | n = 12 units is refused (INAPPLICABLE, `E_COABUNDANCE_SMALL_N`) against the declared minimum of 20 and writes no module output; n = 60 runs; a paired design is refused with its design reason (`E_COABUNDANCE_DESIGN_UNSUPPORTED`) before the unit count is considered | an enabled networks block that declares no analysis is refused (`E_NETWORKS_NOT_DECLARED`, exit 2) and a PPI-only declaration builds no modules: co-abundance never runs by default |
| V162 | PASS | three planted modules recovered with adjusted Rand index ≥ 0.9 (computed in the test; noise features unassigned); eigengenes equal `prcomp()` PC1 of each module's scaled data up to sign; stability reports its seed and draw count; R level: signed adjacency, TOM, ARI and the scale-free fit index equal hand formulas | the unit-bootstrap Jaccard matrix equals a hand loop that resamples units with the same L'Ecuyer stream, every draw is a set of units, and a feature-level bootstrap on the same stream gives a different result (so feature resampling would fail the test) |
| V163 | PASS | module-trait statistics equal a direct `lmFit`/`eBayes(trend = FALSE, robust = TRUE)` on the eigengenes for the group trait and for a numeric phenotype (design `~ group + centred phenotype`); families `module_trait__group` and `module_trait__score`; label `module_level`; the planted group module and the planted phenotype module have q < 0.05; a single module (one-row fit) is tested and equals the direct fit | the protein DEA families contain only `protein-primary`: module rows never enter protein families |
| V164 | PASS | a synthetic PPI snapshot prepared by `resources prepare` (kind `ppi`) loads with its source, release, species, score type and manifest hash recorded; the run executes inside `reproduction.offline()` (sockets refused) | a hash mismatch fails `E_RESOURCE_HASH`; a resource path naming a remote location fails `E_RESOURCE_DOWNLOAD` before analysis |
| V165 | PASS | the planted dense cluster has P at the floor 1/200 and the null-contrast set a non-small P; observed edges, k, null mean, set size and universe size equal a degree-preserving sampling loop over the measured mapped universe written in the test (same seed); the universe is smaller than the snapshot; hub/degree results are `descriptive`; undirected edges are counted once (R-level hand case) | `null_universe: snapshot` fails `E_NETWORK_UNIVERSE`; a null asked to place a gene outside the measured universe errors |

## Notes on the takeover review

This packet was completed after the first builder session was interrupted. The takeover review changed five things. The not-declared case is now a typed declaration refusal; it was previously reported as INAPPLICABLE with the unrelated code `E_COABUNDANCE_SMALL_N`. The scale-free fit index now follows the published pickSoftThreshold definition. Snapshot edges are canonicalised as undirected pairs before counting, and the snapshot source is recorded. The bootstrap became a tested function (`pd_module_stability`) with a hand-loop oracle. The phenotype trait and the repeated-design refusal gained tests.

During a test run, runtime socket attempts are refused by the R11 offline guard with `E_NETWORK_FORBIDDEN` (exit 5). `E_RESOURCE_DOWNLOAD` is the typed refusal of a remote resource location. The analysis code has no download path.

## Follow-up after the first commit (A-2026-10-01-19)

The R14f generality matrix found a crash: with exactly one co-abundance module, the module-trait fit failed (`E_STAGE_EXCEPTION`), because limma drops the row names of a one-row fit and the shared `moderation_diagnostics` then built a frame with no feature IDs. The shared fix restores the row names from the tested universe (`model_limma.R`). The R test "V163 a single module ..." fails before the fix (`evidence/single-module-fail-before.log`) and passes after. In the same change the module-trait moderation was set to trend = FALSE: eigengenes are centred, so their means are identically 0 and a mean-variance trend is undefined. The robust setting follows the primary model. The setting is recorded in `module_trait.tsv` (`moderation`) and in D-49.

## Calibration

The calibration of the connectivity null (dense cluster versus random sets) is recorded with R14f (`docs/validation/015-post-de-analysis/R14f`).

## Acceptance

All gates exited 0 and every acceptance case above, including its negative case, passed. R14e is accepted by the operator-authorized Claude route (self-verified; independent audit NOT_RUN). Synthetic fixtures and synthetic snapshots only; network results make no mechanism or causal claim.
