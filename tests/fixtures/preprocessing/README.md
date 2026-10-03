# R03 preprocessing fixtures (synthetic)

`extreme-sample-abundance.tsv` is the public independent example with
observation C4 shifted by +6 log2 on every feature (an isolated extreme
sample for V027). `constant-abundance.tsv` has eight constant features
(zero retained variable features, V030). All other R03 oracles are built
inline in `r/proteomicsCore/tests/testthat/test-preprocessing.R`,
`test-qc.R` and `tests/integration/test_qc.py`. Nothing here is study data.
