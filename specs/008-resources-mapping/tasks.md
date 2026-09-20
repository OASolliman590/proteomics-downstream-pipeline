# Tasks: Versioned annotation, protein groups and gene sets

**Input:** [spec.md](spec.md), [plan.md](plan.md), parent shared contracts.
**Roadmap:** R07. **Dependencies:** R05.

## Phase A: Preconditions

Read the live constitution, frozen packet contract, prerequisite evidence/commits and existing interfaces. Maintainer issues one packet starter only after audit, freeze and explicit authorization; parallel work is limited to a frozen disjoint-file group.

## Phase B: Independently testable implementation

- [ ] T061 [US3] Implement and independently test **Local resource snapshot registry** in `src/proteomics_pipeline/resources.py; schemas/resource.schema.json`. Acceptance V061: Snapshots identify source/release/species/terms/hash; changed content or missing required files fails before production analysis. Evidence: command/output/version/hash in `docs/validation/008-resources-mapping/V061.json`. Covers FR-061.
- [ ] T062 [US3] Implement and independently test **Explicit resource preparation command** in `src/proteomics_pipeline/resources.py; configs/resources/`. Acceptance V062: Fetch/build occurs only in a separate explicit preparation action, records exact versions and never auto-fetches during a frozen run. Evidence: command/output/version/hash in `docs/validation/008-resources-mapping/V062.json`. Covers FR-062.
- [ ] T063 [US3] Implement and independently test **Species-aware identifier mapping** in `r/proteomicsCore/R/mapping.R`. Acceptance V063: Stable IDs/species keys align; symbol collisions, retired IDs and missing mappings have explicit records. Evidence: command/output/version/hash in `docs/validation/008-resources-mapping/V063.json`. Covers FR-063.
- [ ] T064 [US3] Implement and independently test **Orthology evidence preservation** in `r/proteomicsCore/R/mapping.R`. Acceptance V064: Source/target species and mapping evidence/version are retained; projected sets are labeled projected and ambiguous mappings are reported. Evidence: command/output/version/hash in `docs/validation/008-resources-mapping/V064.json`. Covers FR-064.
- [ ] T065 [US3] Implement and independently test **Label-independent gene representatives** in `r/proteomicsCore/R/gene_matrix.R`. Acceptance V065: Coverage/median/ID representative choice is invariant to sample-label permutation and never chooses minimum P or maximum statistic. Evidence: command/output/version/hash in `docs/validation/008-resources-mapping/V065.json`. Covers FR-065.
- [ ] T066 [US3] Implement and independently test **Ambiguous protein-group policy** in `r/proteomicsCore/R/gene_matrix.R`. Acceptance V066: Multi-gene groups are excluded from default mapping with reasons; declared aggregation sensitivity retains its own universe and estimand. Evidence: command/output/version/hash in `docs/validation/008-resources-mapping/V066.json`. Covers FR-066.
- [ ] T067 [US3] Implement and independently test **Finite pathway matrix and model** in `r/proteomicsCore/R/gene_matrix.R`. Acceptance V067: Gene matrix is finite for required design, models are refit consistently and missing-value loss is quantified rather than silently imputed. Evidence: command/output/version/hash in `docs/validation/008-resources-mapping/V067.json`. Covers FR-067.
- [ ] T068 [US3] Implement and independently test **Gene-set overlap eligibility** in `r/proteomicsCore/R/genesets.R`. Acceptance V068: Declared overlap-size filters are applied before tests, full/eligible/mapped membership counts are retained and tiny universe is explicit. Evidence: command/output/version/hash in `docs/validation/008-resources-mapping/V068.json`. Covers FR-068.
- [ ] T069 [US3] Implement and independently test **Measured gene background for ORA** in `r/proteomicsCore/R/genesets.R`. Acceptance V069: Universe equals eligible tested mapped genes and changes predictably under filtering; whole-genome background is never default. Evidence: command/output/version/hash in `docs/validation/008-resources-mapping/V069.json`. Covers FR-069.
- [ ] T070 [US3] Implement and independently test **Resource/mapping regression fixtures** in `tests/contract/test_resources.py; r/proteomicsCore/tests/testthat/test-mapping.R`. Acceptance V070: Tests cover one-to-many mapping, duplicated symbols, zero overlap, ortholog projection, offline checksum errors and stable representative selection. Evidence: command/output/version/hash in `docs/validation/008-resources-mapping/V070.json`. Covers FR-070.

## Phase C: Review checkpoint

All checkboxes start unchecked. Complete a checkbox only when Maintainer has read the diff, rerun the relevant tests and recorded acceptance evidence plus commit. A task can span multiple Packet Implementer correction tasks. If a runtime/resource is unavailable, keep its verification pending and continue independent work; do not write a successful fixture-shaped substitute.

## Dependencies & Execution Order

The order above follows shared contract → implementation → integration → adversarial verification. Preserve data-model and method boundaries. Resolve packet prerequisites before dependent tasks and obey the exact allowlist in the frozen packet index. A thin vertical demonstration is a checkpoint, not the end of the parent roadmap.
