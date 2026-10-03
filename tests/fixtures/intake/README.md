# R02 intake fixtures (synthetic)

Invented, tiny protein-level inputs for V011-V020. No values come from any
study. Vendor exports under `vendor/` follow the header sets declared in
`configs/mappings/*.json`; their expected canonical rows are written out by
hand in `tests/contract/intake/test_vendor_mappings.py` before the adapter
is run.

`independent-long-shuffled.tsv` is the long form of
`configs/examples/fixtures/independent-abundance.tsv`, built by explicit
(feature_id, observation_id) keys and shuffled with a fixed seed; its two
NA cells are P07/U2 and P08/T1. Malformed-file cases (V020) and the
historical-workbook fixture (V017) are generated inside the tests because
`*.xlsx` files are git-ignored and the defects must be exact.
