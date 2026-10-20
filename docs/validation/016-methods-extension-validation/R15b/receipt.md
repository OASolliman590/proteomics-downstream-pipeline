# R15b receipt — Offline GMT pathway libraries (ADR 0011, FR-174 to FR-178)

Packet R15b (FR-174 to FR-178 / T174 to T178 / V174 to V178) was implemented on 2026-10-13 by the **Claude route (Haiku builder, Opus review)**, after R15a. The packet is ready for review; it is not accepted. Independent audit: **NOT_RUN**.

## Identities

- Contracts: scientific methods v1.4.0 (SM44 with SM14 and SM17); `specs/016-methods-extension-validation/contracts/methods-extension.md` (library and licence rules, `E_LIBRARY_REDISTRIBUTION`, `E_ORA_UNIVERSE`, `E_RANK_*`, `E_ENRICHMENT_FAMILY_POOLED`).
- Working-source manifest (9 files: the R15b allowlist and the two shared registrations): aggregate SHA-256 `702e7b148260ca57da4fcd5cf90c95b4576c5d6d8b644169290f8d6ff16273a7` — `evidence/working-source-manifest.txt`.
- `verified_commit`: **null** (local commit on `claude/methods-016`, not pushed).
- Environment: as R15a (macOS, Python 3.13.15, R 4.6.1, limma 3.68.5, testthat 3.3.2, `.r-lib/`, `LANG=en_US.UTF-8`).

## What was built

- `src/proteomics_pipeline/gene_set_libraries.py` (capability `offline_gene_set_libraries`): hashes every declared local library before any parse (`E_RESOURCE_HASH` on a mismatch), records its source, version, species, identifier type, licence identifier, redistributable flag and a synthetic marker, refuses network locations (`E_RESOURCE_DOWNLOAD`), refuses a restricted library declared redistributable, and checks any artifact member against the licence rule before anything is written (`E_LIBRARY_REDISTRIBUTION`). The R stage runs under a socket guard that makes any socket creation fail with `E_RESOURCE_DOWNLOAD`.
- `r/proteomicsCore/R/gene_set_libraries.R`: the GMT reader, the measured-background check (`E_ORA_UNIVERSE`), ORA through the existing exact adapter `ora_test` (k = 0 sets are kept), the licence rule, family discipline (`enrichment_family_adjust`, refusing a family that spans libraries), and the stage handler `gene_set_libraries_stage`. The outputs are `pathways/libraries.json` (manifest records, no library content, no local path), `ora_<library>.tsv` and `prerank_<library>.tsv` (exploratory).
- `r/proteomicsCore/R/pathway_prerank_welch.R`: ranks from the Welch t of R15a (`E_RANK_DUPLICATE`; ties refused as `E_RANK_TIE` unless the SM15 stable-ID rule breaks them; `E_RANK_INVALID` for non-finite values); the running-sum score (hit weights |t|, miss step 1 / (N − k)); the exact permutation P by enumerating every subset of the set size, which is the named alternative used when `choose(N, k) <= 200000`; and the existing exploratory fgsea adapter (`run_fgsea`) for larger rankings. The estimator used is recorded per row.
- Shared, recorded under A-2026-10-01-45: `dispatch.R` (handler), `runtime.py` (capability map), `docs/methods/README.md` (capability row) and the regenerated reason-code catalogue.
- Fixtures: `tests/fixtures/methods_extension/pathway_libraries/` (a generator and one committed synthetic GMT; every row carries the marker `SYNTHETIC_FIXTURE`).

## Decisions to review

1. **Named alternative for pre-rank P.** FR-177 allows the existing fgsea adapter or a documented named alternative. The exact subset enumeration is used when feasible (it is the only way to meet V177's exact-P assertion); fgsea is used otherwise, labelled `prerank_welch_fgsea_exploratory`.
2. **Licence rule.** A restricted licence (KEGG or MSigDB prefix) cannot be declared redistributable. A library that is not synthetic cannot be part of an artifact (only marked synthetic fixtures may be redistributed), even under an open licence.
3. **Family key.** A family is `library|method|contrast` (`gene_set_family_key`); methods are `ora_up`, `ora_down` and `prerank_welch`.
4. **Foreground.** ORA foreground is the Welch family q within the threshold (default 0.05) with the stated direction, on the measured universe of tested, mapped genes.
5. **Wired into the run (F2).** The `methods_extension.pathways` block (both schema copies; Phase 5 only) declares the libraries. The plan verifies every library by SHA-256 and licence before it is frozen, and the capability is planned with `required` from `execution_requirement`. The stage runs after the primary differential stage on `dea/zero_null.tsv`, with one representative gene per feature and contrast (the smallest stable feature ID wins a tie, and features without a mapped gene are left out). Its measured universe is the mapped genes of the run (one universe for every contrast, SM17). Outputs: `methods_extension/pathways/libraries.json` (manifest records: hash, source, version, species, identifier type, licence, redistributable and synthetic flags; no path), `ora_<library>.tsv` and `prerank_<library>.tsv`. The run catalogue draws `path_enrichment_dotplot` for each library and direction with at least one foreground hit.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `python -m pytest tests/scientific/test_gene_set_libraries.py -q -rs -p no:cacheprovider` | 0 | 14 passed, 0 skipped | `evidence/python-gene-set-libraries.log` |
| testthat `test-gene-set-libraries.R` (`stop_on_failure = TRUE`) | 0 | 40 expectations, 0 failed, 0 skipped | `evidence/r-gene-set-libraries-testthat.log` |
| `python -m pytest tests/regression tests/unit tests/scientific/reference/test_text_encoding.py -q -p no:cacheprovider` | 0 | 91 passed (shared with R15a) | `evidence/gates-regression-unit-encoding.log` |
| `python scripts/check_spec_kit.py` | 0 | PASS | `evidence/spec-kit-check.log` |

Full-suite results are recorded in the section "Full suites" below.

## Acceptance cases

| Case | Status | Positive oracle | Negative case |
|---|---|---|---|
| V174 | PASS | `hashlib` SHA-256 over the file bytes equals the manifest hash; the run's `libraries.json` records the same hash with source, version, species, identifier type, licence and redistributable flag; no local path is recorded. | A one-byte edit after declaration fails `E_RESOURCE_HASH`. A network location fails `E_RESOURCE_DOWNLOAD`; a socket created under the run guard fails `E_RESOURCE_DOWNLOAD` and the guard is removed afterwards. |
| V175 | PASS | Every GMT file that git would publish (tracked or not ignored) carries the synthetic marker; a valid artifact made of the synthetic fixture is accepted. | The KEGG-licensed, the MSigDB-licensed and the non-synthetic redistributable variants each fail `E_LIBRARY_REDISTRIBUTION`, and no artifact directory is created. |
| V176 | PASS | Measured universe N = 500, set A K = 30 with k = 8 and n = 40: p equals `stats::phyper(k − 1, K, N − K, n, lower.tail = FALSE)` within 1e-12 relative; the k = 0 set has p = 1 and is kept; BH equals `p.adjust` within 1e-12 inside the family. | A whole-genome background (N = 20,000) fails `E_ORA_UNIVERSE`, as does a declared size that differs from the measured universe. |
| V177 | PASS | Hand-computed running sum: the six-gene ranking {5, 4, 3, 2, 1, 0.5} with set {1, 3} gives ES = 0.75; the stage's ES and exact P equal a Python brute force over all `combn(10, 3)` subsets within 1e-12; the rows carry `exploratory = true`, `claim_label = exploratory` and `estimator = prerank_welch_exact_enumeration`. | Duplicated gene IDs fail `E_RANK_DUPLICATE`; unbroken ties fail `E_RANK_TIE`, and the SM15 stable-ID rule resolves them (smallest ID first). |
| V178 | PASS | Two libraries × three methods × two contrasts give 12 distinct families, each of the form `library|method|contrast`; every family's q equals `p.adjust` on that family's P values within 1e-12; no family spans two libraries. | A BH adjustment whose family key spans libraries fails `E_ENRICHMENT_FAMILY_POOLED`. |

## Full suites

Final run on the reviewed tree (after the follow-up group F1–F4):

| Suite | Command | Exit | Result | Log |
|---|---|---|---|---|
| Full Python suite | `python -m pytest tests -q -rs -p no:cacheprovider` | 0 | 596 passed, 0 skipped (36 min) | `evidence/full-python-suite.log` |
| Full R suite | `testthat::test_dir("r/proteomicsCore/tests/testthat", load_package = "installed", stop_on_failure = TRUE)` | 0 | 2249 expectations, 0 failed, 0 skipped, 20 files | `evidence/r-all-testthat.log` |

The first full Python run of this stage (before F1 to F4) was stopped at about 24 percent: the post-DE generality test (`test_v167_*`) runs one full pipeline per cell and takes minutes per cell, so the run had not hung; it was stopped to start the final run. No result from that run is recorded.

## Gaps and doubts (for the reviewer)

- **Representative gene.** The stage keeps the smallest stable feature ID per gene and contrast. The full SM15 rule (greatest coverage, then median abundance, then the stable ID) needs the coverage and median columns of the R07 mapping, which the run does not pass to this stage; the simpler rule is a documented gap.
- **Universe.** The measured universe counts only tested features with a mapped gene (the eligible, mapped and measured set of SM17).
- **New reason codes beyond the spec's table:** `E_LIBRARY_FORMAT` (malformed GMT) and `E_RANK_INVALID` (non-finite ranks). Both are documented in the methods-extension contract and in the regenerated catalogue. The Maintainer should confirm them.
- **Licence identifiers.** The fixtures use the synthetic identifiers `SYNTHETIC-CC0` and `KEGG-restricted`/`MSigDB-restricted`. Real identifiers must be confirmed by the Maintainer before any library is used (D-69).
- **Socket guard scope.** The guard covers socket creation in the Python process that runs the stage. Rscript runs as a subprocess and is not sandboxed; no network use is present in the stage code.
- **Prerank scale.** The exact enumeration is limited to 200,000 subsets; larger rankings use the exploratory fgsea adapter whose P values are randomised and are not exact.
- **Maintainer ratification.** Amendments A-2026-10-01-45 and A-2026-10-01-46 are drafted by the implementing route and are pending Maintainer confirmation.

## Acceptance

All packet gates exited 0 and every acceptance case above, including its negative case, passed. R15b is ready for review; it is not accepted until the reviewer and the Maintainer confirm the amendments and the licence identifiers.
