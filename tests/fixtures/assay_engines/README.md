# R06 assay-engine fixtures (synthetic)

`tests/scientific/test_assay_engines.py` writes small LFQ matrices with
seeded synthetic values, dropout (missing cells, including one feature absent
from a whole group) and invented peptide-count evidence. No vendor or study
data are used; the references are direct DEqMS/proDA calls in R.
