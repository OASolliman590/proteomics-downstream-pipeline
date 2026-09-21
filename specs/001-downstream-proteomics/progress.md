# Progress — evidence, not anticipated completion

Kit: **1.2.0-frozen** on **2026-09-20**. Baseline: recovered study code at `3aefdd95c46a1c56dae89e79ddf441b6cc978148`. Maintained successor: **NOT_STARTED**.

| Gate / scope | State |
|---|---|
| Specification hardening artifact | FROZEN v1.2.0 as uncommitted working-tree files; actual checkout checker, schema cases, protected trees and registered source hashes pass. |
| Independent candidate audit | Initial Astra audit found three P0 and two P1 items. Sol corrected dispatch barriers, explicit engine requiredness/status rules, executable R03 V027/V029 seams and the frozen R function map. Fresh bounded Astra re-audit found no P0; its final P1 registration-wording ambiguity was resolved without changing interface meaning. |
| Maintainer freeze / implementation authorization | Freeze recorded 2026-09-20. The user authorized R01 only; no commit or push is authorized. |
| R01–R05 / R10a (Phase 1) | R01 accepted on 2026-09-21 with V001–V010 PASS. R02 is the next dependency-eligible packet but remains unstarted pending explicit authorization; R03–R05/R10a remain blocked on serial prerequisites. |
| R06–R09 (Phase 2) | Not started; requires accepted Phase 1/R10a and separate authorization. |
| R10b / R11 / R12 (Phase 3) | Not started; requires exact prerequisites and separate authorization. |
| R execution, maintained numerical qualification, locks/calibration | R01 foundation R execution PASS; later scientific qualification and R11 locks/calibration remain NOT_RUN. |
| Private legacy regression | NOT_RUN; Maintainer-only |
| License / ownership / public-release authorization | Unresolved; no license selected |

[Traceability](traceability.json) contains no implementation evidence and every `verified_commit` remains null. A schema/document check is not a package or scientific PASS. The next action is an explicit user decision to authorize R01; later packets remain separate. The GitHub repository has not been modified by this uncommitted working tree.

<!-- FREEZE_RECEIPT_START -->
## Freeze receipt

Freeze version/date: `1.2.0`, `2026-09-20`.

Baseline HEAD: `3aefdd95c46a1c56dae89e79ddf441b6cc978148`.

Working-source manifest SHA-256: `84149fd7526e07346ba914d349939f9945cd3f0f0e6653fc1591f8b807101120` over the 80 overlay paths plus `AGENTS.md`, sorted by repository-relative path. Raw bytes are hashed except that this self-referential receipt block is removed from `progress.md` before hashing.

Frozen contract-manifest SHA-256: `0619003741ee23f93ed08c332c1bfce3e21b446f277073be7010fa68286a6ffa` over the 80 normative paths listed below. `checklists/structure-validation.json` is the only whole-file exclusion because it is solely mutable verification/status evidence. For the contract digest only: this receipt block is removed; progress-table state cells become `<STATUS>`; completion checkbox states in task/checklist files become unchecked while their text remains; and only `status`, `evidence` and `verified_commit` values in `traceability.json` are normalized (`<STATUS>`, `[]`, `null`) before canonical sorted compact JSON serialization. No schema, scientific contract, ownership rule, slice requirement or acceptance definition is excluded. Each line is `normalized_sha256  path`; the aggregate is SHA-256 of the lexicographically sorted UTF-8 lines including LF.

Independent audit references:

- Initial specification audit: Codex task `astra / spec-freeze`, thread `01a0c015-859a-76f1-ae8a-47fddc179414`, reviewed working manifest `e65317b7873d63361399fe8b55f7b90c893f0ae92b5193ff119f6c642a9ee90d`.
- Corrected-contract re-audit: fresh Codex task `astra / spec-freeze`, thread `01a0c020-d2a2-7f40-b2af-947fa29ff169`, reviewed working manifest `be73544c0999d85adb46c55d5d9e56219b5c90c2e50ec07f80e6be96369c0e85` and found no P0.

Finding dispositions:

- P0 dispatch dependencies: fixed in ownership/index/phase records; R06/R07/R09 require R05+R10a and R08 requires R06+R07+R09; re-audit verified an acyclic graph.
- P0 engine requiredness: fixed with mandatory `model.execution_requirement`, required primary models, installation-independent scientific eligibility and exact V059 FAILED/PARTIAL outcomes; schema positive/negative cases passed and re-audit verified the meaning.
- P0 R03 seams: fixed with the R01 canonical hash seam plus mandatory V027/V039 actual-plan regression, and R03-owned detection Fisher/BH plus mandatory post-R05 V029 family-separation regression; re-audit verified ownership compatibility.
- Initial P1 items: fixed by freezing exact R handler names/signatures and qualifying status/ownership evidence.
- Re-audit P1: resolved by stating that R01's hardcoded map plus handler existence is registration; no separate metadata file exists. This wording clarification does not change the audited interface.

Frozen dispatch units (`packet: phase; dependencies; acceptance; risk`):

- `R01: 1; frozen baseline/contracts; V001–V010; high`
- `R02: 1; R01; V011–V020; high`
- `R03: 1; R02; V021–V030; high`
- `R04: 1; R03; V031–V040; high`
- `R05: 1; R04; V041–V050; high`
- `R10a: 1; R05; V091–V094; high`
- `R06: 2; R05+R10a; V051–V060; high`
- `R07: 2; R05+R10a; V061–V070; high`
- `R09: 2; R05+R10a; V081–V090; high`
- `R08: 2; R06+R07+R09 whole-group barrier; V071–V080; high`
- `R10b: 3; R10a+R08+R09; V095–V100; high`
- `R11: 3; R10b; V101–V110; high`
- `R12: 3; R11; V111–V120; high`

Audit checkpoints: independent freeze audit completed; every dispatch unit is high risk and remains subject to Sol diff/runtime verification; shared-interface defects return to Sol; predefined risky diffs receive a fresh Astra source/test audit; R12 private regression remains Maintainer-only. The sole legal parallel group is `R06 || R07 || R09` after accepted R05 and R10a, followed by the whole-group barrier before R08.

External gates still unresolved: `Rscript` and maintained R/scientific execution are `NOT_RUN`; private archived-study regression is `NOT_RUN`; license/ownership/public-successor-release authorization and remote description consistency remain unresolved. These do not become PASS through this freeze. `verified_commit: null`; there is no authorized commit or push.

Normative file hashes:

```text
84b846cd73af8c85f0ea4817a074d034a33f3e23f9451d9083fd1e5d04a39764  .specify/memory/constitution.md
149cc11709a5851f87a679e57fed5cb3fa79a1c0f929b436cd9acde2f375eb4f  AGENTS.md
0384a194023de91e9d1952c50eb44716e190a14bd871ba2ed5ef97c868565f6a  README.md
bf1769487738daa06b4b2aa496aba6c7f390789e5091dda6855f8d119ba1a476  docs/RELEASE_CHECKLIST.md
de9c064f65e23d0c1d3bebc99ceb52d0db9f4fece2ba82e2d99534e5ec3ac349  docs/SPEC_HARDENING_REPORT.md
5ba9e30cc1308b6b0a407a58d0f52ece9f90d5ae585963b6319cec8c4fdbc8ed  docs/adr/0001-phase-milestones.md
96dc38e5162fd4321e6334f53d271a5d98cd11276970d82159ef33a097b096aa  docs/adr/0002-primary-engine-limma.md
a99cf99edb29c762faca3df0f9994cd39eb44087a1503292f35ba3a8e4b98116  docs/adr/0003-response-descriptive-default.md
11aad7aea0e20b7f86c322541ab41e3dc67ce63ab08417dbce3baa948ea480df  docs/adr/0004-enrichment-dispatch.md
b74ff4c47ffc9127ebd3d21e4290ebc2516a5c6618a337aedfbe37ee26ec11bc  docs/adr/0005-no-silent-method-fallback.md
71a775cf58101ecf7eac1ab63c49378608c32626d2e3fc021095af5aa4c2233d  specs/001-downstream-proteomics/IMPLEMENTATION_BRIEF_R01.md
dc9190d3f9c20182d24f139e409505bb1fb5627fad78dbb38f4523753166ce8e  specs/001-downstream-proteomics/PHASES.md
1d74b5cf579ae3b7b64a199c7d3f417b52794fc2ea685c1e37865c32371b4a62  specs/001-downstream-proteomics/START_HERE.md
a1701818e40d737b83e7f15dcf4e0e3db5bb6c3f29193d98266ac4349d9936f1  specs/001-downstream-proteomics/checklists/spec-quality.md
a5e305c6f8a1ca09ee167ddfbb2cfb1b25201a398ce28e2202d91d2b728294f2  specs/001-downstream-proteomics/contracts/analysis.schema.json
492dda4017079131be4cd32f63059c2813650455a6ffe0e75786ec7213522946  specs/001-downstream-proteomics/contracts/cli-and-artifacts.md
59a4b04d59fceabdfd43a6b132ff74b3e8b2db1d3ff6dd7ebdf443d4ff2d48a6  specs/001-downstream-proteomics/contracts/example-effect-threshold.json
6d96680dc4a6cc546183315f9adb756fe31282c196050c2eb0176d7cb78d6216  specs/001-downstream-proteomics/contracts/example-independent.json
80098aac14357e15eb26b57e925827c2a7f9373e2352e4761c20f86635218889  specs/001-downstream-proteomics/contracts/example-paired.json
a44f6ec06a08168f01b8b2da08391993ede7a6aa2b6dc72e0019e661a780dcf4  specs/001-downstream-proteomics/contracts/fixtures/README.md
478ff2c26e998fd372344a92ef4b8e01e43b185d2e29edb2898c75a220cc6900  specs/001-downstream-proteomics/contracts/fixtures/features.tsv
fc331a5709031f5574b31f01200ebe73abf6b88be6898b8f6540683d5c0449ec  specs/001-downstream-proteomics/contracts/fixtures/independent-abundance.tsv
193ada9c1aebfc19535ddc1739180d7e226803bc77ac0ee18a61633bfdcbc1e9  specs/001-downstream-proteomics/contracts/fixtures/independent-observations.tsv
b6a2572d8051d12bffeac701adada507d7b7eb3f52733f0438ee358218df150a  specs/001-downstream-proteomics/contracts/fixtures/paired-abundance.tsv
0c7d75a4e1ee7975e2ff52035545304b7a7dec3bd979ac2bafa5650271a9af96  specs/001-downstream-proteomics/contracts/fixtures/paired-observations.tsv
fbac285fb477d16a368ee600fa8f76420a0274b9b513f984e6764d2253b7d6aa  specs/001-downstream-proteomics/contracts/fixtures/provenance.json
cecabd011c989489abbebbc7699aef96e4dc88d027ecf7e40f2eff4c8e1c67ff  specs/001-downstream-proteomics/contracts/run-status.schema.json
11fd29817d3a021a56407e63e86f6ce6eb56c7edb70af9f1603beb1990952a65  specs/001-downstream-proteomics/contracts/scientific-methods.md
c9f57ae6e479f0c2bb60f98f120db3bfacca129aaeeb541afa62318d0106c76f  specs/001-downstream-proteomics/contracts/semantic-validation.md
a80eee2d1949d00a56860beabc00c76119dda57422aaa072b88d5b21beaa032a  specs/001-downstream-proteomics/contracts/stage-request.schema.json
b85150b5262bf0c7d63f88e5f252fefdceac65a1e616cd6596b06157da6dc91f  specs/001-downstream-proteomics/contracts/stage-result.schema.json
8815fc92b1b76ed65de80bb15655297871f774655a3d8aa668585b5510879346  specs/001-downstream-proteomics/data-model.md
54cb00c3aef9bd7f677d4566636ce048dca98cf1a1a0bd786e617047b02f79f8  specs/001-downstream-proteomics/packet-index.md
4e35f8f1ede365d221a73c101f4b176497b91e8198dadb86e05492df570aa78e  specs/001-downstream-proteomics/packet-ownership.json
1f493e99c3089b67c7787ceb264012e4e7b228bf8bffd81aa69538325b7bf442  specs/001-downstream-proteomics/plan.md
4bb6e14fa805787975c21bef8e0bd9fdd61368dc4e52e34eaa6a508b8445d71a  specs/001-downstream-proteomics/progress.md
9d4454e8353fe608deb5a06aafbc04b8737e63ea67fa57f5766a582367e16e08  specs/001-downstream-proteomics/quickstart.md
b6573de79648dcdb76995d29d6379bcaf01222b1feb3dc8c061db688d01c0b6e  specs/001-downstream-proteomics/research.md
af6eab98413ab3df195fb1ca6000e613e54296560545ede1638e6a6f6bf25c5c  specs/001-downstream-proteomics/roadmap.json
dbd1ba755b8c67be1d6c40b900b2198117e1a4a1f545959e41f8936f096c5033  specs/001-downstream-proteomics/roadmap.md
23917ecb5de03b65ae814fc2b18b775bced401ff08e0f681481b260c9d897458  specs/001-downstream-proteomics/spec.md
17f32275072b2fa0c20f7683d0276e5ca786f5fdd2bee641bc4793d312eeca48  specs/001-downstream-proteomics/tasks.md
785324a0cb9ec372629a4e6b977e5314e5f29c6133196866f173b32052411959  specs/001-downstream-proteomics/traceability.json
76f8f2a976d0953ccf16508c16414813ac2eae5552ffa31c27b3adf45ee8d69e  specs/001-downstream-proteomics/validation-strategy.md
5a4550b35a3e2201a0d68a708d22c2a4cdb9e7858ce0b8290d4b171785f2c7d2  specs/002-foundation/plan.md
62b3d59a654df10ce43f0242c52f537171aa87bf5ec7a97a3d091504a232c882  specs/002-foundation/spec.md
34c85c6eea77d0094c68402198be74cd4f96c21ab5dbf12e5940e420878ecdd2  specs/002-foundation/tasks.md
b4daf2a41bbbcf46b0a70b940dc31623c7f162906e1c48fc98d8bdaacb244d50  specs/003-intake/plan.md
cb2b04b856a6091be4f1cf089c50f303283d0268104781bece6506079c1b8ab7  specs/003-intake/spec.md
d431af7008ab1ca7079d7cfca2c9801932101ac0e9b1c0fb7954c90cbce23ef2  specs/003-intake/tasks.md
10f5fc174c344b4dc8fc56486d3210b5c4ecec5a07299fdeac4c9ced2c209576  specs/004-preprocessing-qc/plan.md
c146c53f02456bd1288d50ad7da5239c06e90e5e69d231dbf6298f1f7adc17ce  specs/004-preprocessing-qc/spec.md
cf6fcaccb0009490008ad12a2d595eaad87477f4c294901dc02d6c942469b363  specs/004-preprocessing-qc/tasks.md
d60eb75c04b7478dc3ec5ae21f2c12e5a23e5074912936dba35e2b1741ba4d08  specs/005-design-contrasts/plan.md
0cf12338fe5f220d8a01862ec920ad742602abb49e72a799bfca81c230781626  specs/005-design-contrasts/spec.md
8bb78a29088ae1d4bd3fcdf29b6f4a51bc52feb3173c7239b2c98dff704ff3b3  specs/005-design-contrasts/tasks.md
9d0d8e61a8700f5982c93dc079b6350561504e17ab305858161b8b0d12dc6a45  specs/006-limma-inference/plan.md
4bc60e4b5efd4c28ecac4713fc8c395ed4dd43dbae86277d6a6296534a2f2137  specs/006-limma-inference/spec.md
edac343aeac1149577821b29895ad353b78f831cd45bce624cf7800e13868426  specs/006-limma-inference/tasks.md
663e230e37bfd0ecb6b9b704d6dc408697dee5836bd036dc10fef614395e1270  specs/007-assay-engines/plan.md
ed8b5e4b821a762b1c41f5a6f7193e309e927e858138c960d3b1bb3f734703a8  specs/007-assay-engines/spec.md
3141fa001e3c5a2d9a71914492d4ea23a7cdfce01fc39ac0cc3e3bd1e8b9e696  specs/007-assay-engines/tasks.md
7b170d89c6e40c689d4601319001ac0a2079488afb68135d28e83fb2708c78f5  specs/008-resources-mapping/plan.md
961af067a199f1b188ba50463535ad7427954a2808b8ddd4aa927cb5bbf7b978  specs/008-resources-mapping/spec.md
a33824a224314aec8135c3c738928743633ec67622b701b304a76fcb70ecf4b4  specs/008-resources-mapping/tasks.md
fae778462461c26d4944afa31af1d6bacbc542cee8eb325837fde14ffb95081f  specs/009-enrichment/plan.md
49bfab4e9c7a6eb33f6de54f54dce22d3e27d73f94d17a28e37af1420a33e0a0  specs/009-enrichment/spec.md
6aa15902af6dff8998c25c53316125e32d32d3a96bcb90786207f3a965d1e262  specs/009-enrichment/tasks.md
4e2b7e86f57fc6a632d8680d90b91acd24549123c07517fdea74486ab3aaad33  specs/010-treatment-response/plan.md
62e22b8a85a50df9497bcf6eb353492882e58cd3c9656e2a3e4daff3c2b60785  specs/010-treatment-response/spec.md
2ce4b32933a9ffc0b01992b911184ca0e885325b25072a656ad23c07c8948119  specs/010-treatment-response/tasks.md
8eeda26221b4eca972753911c701dbec9e84397ef63e30ec5f48d2ed593e7f0a  specs/011-reporting/plan.md
45014b7c243466b00595fa48518c42c52534b4126fdec890fce98597536db9cc  specs/011-reporting/spec.md
6d1ecede3287d315408a8c9c4ce9d5195c12d5b01347c303b495179ad69f5ee2  specs/011-reporting/tasks.md
bb3bf3396c0c8612b1fbf097f3bfe2dc11718903fecfa6dda2932de21764bf82  specs/012-validation/plan.md
3130486ef7edd5988a8077fe5c520f840c9581ef505e03338992d8af2c73d93a  specs/012-validation/spec.md
486317611a9ca5a292805a2c3661742f03dd7d4996ffdccff9c9b5cf4d3f1d66  specs/012-validation/tasks.md
a832984d8065356885514dce626114dfdcf80c5c824414b7f4d6ca80d4d205dd  specs/013-release/plan.md
9cde1cd19d2f6eb766b7f8f225499882de5854c83f2a1c98f351317913d67f70  specs/013-release/spec.md
cec3829f06edf8faaedf6d2dc1640921bdbb218cc2df6cb0a3906124da8d3699  specs/013-release/tasks.md
```
<!-- FREEZE_RECEIPT_END -->
