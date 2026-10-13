"""Multi-study validation framework (packet R15d, FR-185 to FR-188; ADR 0011). Private inputs stay local."""
from .manifest import MANIFEST_SCHEMA_VERSION, SLOTS, DESIGNS, FORMATS, validate_manifest, map_columns, maxquant_filters, tmt_batch_confounding  # noqa: F401
from .runner import run_study, status_record, commit_status  # noqa: F401
