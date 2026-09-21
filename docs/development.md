# Maintained development foundation

R01 provides an installable Python package and the real proteomicsCore R I/O seam. Create a project-local environment with scripts/maintained/bootstrap.ps1 or scripts/maintained/bootstrap.sh; the scripts never install globally and never bootstrap automatically during analysis.

Schema-only validation is the only real validation mode in R01. Full validation, planning and analysis execution return typed E_CAPABILITY_NOT_IMPLEMENTED until their owning packets provide actual handlers. Missing R, R packages, or testthat is NOT_RUN/unavailable; no result is mocked.
