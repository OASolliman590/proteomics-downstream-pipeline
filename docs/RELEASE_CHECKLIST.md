# Release preparation

## Completed in this recovery

- [x] Identify both Method A and Method B and reconstruct stage/contrast mappings.
- [x] Preserve the original ZIP, extracted files and source hashes outside the proposed repository.
- [x] Independently compare workbook input, matrix exports, fold changes and BH adjustments.
- [x] Execute Method A intake and compare reproduced tables.
- [x] Restore Method B TSV inputs with original group/rat metadata and missingness.
- [x] Provide scientific audit, technical audit, appraisal and evidence records.
- [x] Provide a portable execution wrapper, input validation and ten passing Python tests.
- [x] Add code/data exclusions and a CI workflow limited to Python checks.

## Before publishing any repository

- [ ] Confirm code ownership, authorship, citation text and license. None was supplied or inferred.
- [ ] Decide public versus private repository and whether audit-level study results may be disclosed. Original workbook, abundance, individual-rat tables and generated results remain outside the Git working tree.
- [ ] Review the initial file list and first diff, then choose repository name/remote. No push is part of this recovery.
- [ ] Review MSigDB/KEGG and other external-resource attribution and redistribution terms before adding snapshots or result/data releases.

An archival code release may retain known limitations prominently. It must not claim a validated complete pipeline.

## Before claiming reproducible or scientifically validated analysis

- [ ] Provision R and all dependencies on the target platform, execute every stage from a clean run directory, compare outputs to archived references, and record expected numerical differences.
- [ ] Capture exact versions from the actual analysis processes; create and test a lockfile or immutable environment.
- [ ] Retrieve the historical gene-set release if possible; otherwise define a new release, save mapping/snapshot provenance and document that enrichment is a reanalysis.
- [ ] Recover tissue, rat independence confirmation, batch/randomization, sex/age, dose/duration and original MS processing/identification metadata.
- [ ] Prespecify primary model, filters, estimands, effect threshold and multiple-testing families.
- [ ] Correct the broad no-DEP summary, selected-score inference, exact permutation arithmetic and overshoot terminology in a separately versioned successor.
- [ ] Add residual/design diagnostics, finite/estimable contrast checks, gene-mapping sensitivity and correlation-aware pathway sensitivity.
- [ ] Add scientifically meaningful R tests: independent null simulations with feature selection, exact enumeration ties, overshoot/equivalence semantics, all-missing groups and empty-result reporting.
- [ ] Handle no significant pathways / no candidates as valid outcomes and preserve warnings and numerical diagnostics.

The Python CI workflow alone cannot satisfy these scientific or R reproducibility requirements.
