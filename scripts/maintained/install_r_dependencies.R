#!/usr/bin/env Rscript
# Explicit, project-local installation of the R packages the implemented
# stages use (Maintainer amendment A-2026-10-01-04).  Never run during an
# analysis; never installs globally.  Usage (from the repository root):
#   R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla scripts/maintained/install_r_dependencies.R [--with-svg]
args <- commandArgs(trailingOnly = TRUE)
lib <- Sys.getenv("R_LIBS_USER")
if (!nzchar(lib)) stop("Set R_LIBS_USER to a project-local library (for example $PWD/.r-lib).", call. = FALSE)
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .Library))
options(repos = c(CRAN = Sys.getenv("PROTEOMICS_CRAN", "https://cran.rstudio.com")), timeout = 600)
cran <- c("jsonlite", "openssl", "testthat", "statmod", "vegan", "permute", "BiocManager")
if ("--with-svg" %in% args) cran <- c(cran, "svglite")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) utils::install.packages(missing, lib = lib)
bioc <- c("limma", "impute")
missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, lib = lib, update = FALSE, ask = FALSE)
versions <- vapply(c(cran, bioc), function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "MISSING", "")
print(versions)
if (any(versions == "MISSING")) quit(status = 3L)
