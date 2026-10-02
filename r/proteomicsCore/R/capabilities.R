# Operational availability of statistical engines (packet R04, V038).
#
# Scientific eligibility is decided from the frozen design/adapter table in
# the planner without consulting installation state.  This function only
# reports whether the software for an already-eligible engine is installed;
# an absent adapter is NOT_RUN, never scientific inapplicability, and is
# never replaced by limma under another engine's name.

engine_packages <- function(engine) {
  switch(engine,
         limma = c("limma", "statmod"),
         deqms = c("limma", "DEqMS"),
         proda = c("proDA"),
         permanova = c("vegan", "permute"),
         stop(sprintf("E_ENGINE_UNKNOWN: %s", engine), call. = FALSE))
}

operational_availability <- function(engine, installed = function(pkg) requireNamespace(pkg, quietly = TRUE)) {
  packages <- engine_packages(engine)
  present <- vapply(packages, installed, logical(1))
  handler <- switch(engine, limma = "limma_stage", deqms = "assay_engine_stage", proda = "assay_engine_stage", permanova = "permanova_stage")
  implemented <- exists(handler, envir = asNamespace("proteomicsCore"), inherits = FALSE)
  list(engine = engine, packages = as.list(present), implemented = implemented,
       state = if (implemented && all(present)) "AVAILABLE" else "NOT_RUN",
       reason_code = if (!implemented) "E_CAPABILITY_NOT_IMPLEMENTED" else if (!all(present)) "E_ENGINE_NOT_AVAILABLE" else NULL)
}
