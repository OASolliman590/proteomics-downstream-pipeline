# Exploratory preranked fgsea (packet R08, SM16, V075).  Ranks are native signed
# zero-null statistics mapped through the fixed representatives; ties follow
# stable gene-ID order.  The null is gene-set randomization, not sample labels.
validate_ranks <- function(ranks) {
  if (!length(ranks) || any(!is.finite(ranks)) || anyDuplicated(names(ranks)) || any(is.na(names(ranks))))
    stop("E_FGSEA_RANK_INVALID: ranks must be finite, uniquely named genes", call. = FALSE)
  ranks[order(-ranks, names(ranks), method = "radix")]
}

run_fgsea <- function(ranks, pathways, eps, seed) {
  if (!requireNamespace("fgsea", quietly = TRUE)) stop("E_ENGINE_NOT_AVAILABLE: fgsea is not installed", call. = FALSE)
  ranks <- validate_ranks(ranks)
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  warn <- character()
  res <- withCallingHandlers(fgsea::fgsea(pathways = pathways, stats = ranks, minSize = 1, maxSize = length(ranks) - 1, eps = eps, nproc = 1),
                             warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(table = data.frame(set_id = res$pathway, es = res$ES, nes = res$NES, p_value = res$pval, native_q_value = res$padj, log2err = res$log2err, size = res$size,
                          leading_edge = vapply(res$leadingEdge, function(x) as.character(jsonlite::toJSON(x)), ""), stringsAsFactors = FALSE),
       warnings = warn)
}
