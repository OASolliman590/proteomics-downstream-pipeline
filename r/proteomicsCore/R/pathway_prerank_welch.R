# Pre-ranked enrichment on the Welch statistic (packet R15b: FR-177, SM44, SM15; ADR 0011).
#
# Ranks are the Welch t of R15a, one per gene.  Ties are refused unless the SM15 stable-ID rule breaks them
# (E_RANK_TIE); duplicated genes are refused (E_RANK_DUPLICATE).  The enrichment score uses hit weights |t| and
# miss steps 1 / (N - k).  The exact permutation P enumerates every subset of the set size, which is the
# documented named alternative when feasible; larger rankings use the existing exploratory fgsea adapter
# (run_fgsea, pathway_fgsea.R).  Every output carries the exploratory label.

# Ranks on the Welch statistic (FR-177): unique gene IDs, finite values, ties only under the SM15 stable-ID rule.
prerank_welch_ranks <- function(statistic_by_gene, tie_break = "none") {
  if (anyDuplicated(names(statistic_by_gene))) stop("E_RANK_DUPLICATE: a gene appears more than once in the ranking", call. = FALSE)
  if (any(!is.finite(statistic_by_gene))) stop("E_RANK_INVALID: ranks must be finite", call. = FALSE)
  if (anyDuplicated(unname(statistic_by_gene)) && !identical(tie_break, "sm15_stable_id"))
    stop("E_RANK_TIE: tied statistics are not broken by the SM15 stable-ID rule", call. = FALSE)
  statistic_by_gene[order(-statistic_by_gene, names(statistic_by_gene), method = "radix")]
}

# Weighted running sum (gseaParam = 1; hit weights |t|, miss step 1 / (N - k)); ES is the maximum deviation, signed.
prerank_running_es <- function(ranked_statistic, in_set) {
  n <- length(ranked_statistic); k <- sum(in_set)
  if (k == 0L || k == n) return(NA_real_)
  hit <- ifelse(in_set, abs(ranked_statistic), 0)
  hit <- hit / sum(hit[in_set])
  miss <- ifelse(in_set, 0, 1 / (n - k))
  running <- cumsum(hit - miss)
  running[which.max(abs(running))]
}

# Exact permutation P by enumerating every subset of the set size; one-sided in the direction of the observed ES.
prerank_exact <- function(ranked_statistic, set_index) {
  n <- length(ranked_statistic); k <- length(set_index)
  observed_mask <- seq_len(n) %in% set_index
  es_obs <- prerank_running_es(ranked_statistic, observed_mask)
  subsets <- utils::combn(n, k)
  es_all <- apply(subsets, 2L, function(s) prerank_running_es(ranked_statistic, seq_len(n) %in% s))
  p <- if (es_obs >= 0) mean(es_all >= es_obs - 1e-12) else mean(es_all <= es_obs + 1e-12)
  list(es = es_obs, p_value = p, subsets = ncol(subsets))
}

PRERANK_EXACT_LIMIT <- 200000

prerank_set_result <- function(ranked_statistic, set_index, set_id, seed = 20261013L) {
  n <- length(ranked_statistic); k <- length(set_index)
  if (choose(n, k) <= PRERANK_EXACT_LIMIT) {
    res <- prerank_exact(ranked_statistic, set_index)
    list(set_id = set_id, size = k, es = res$es, p_value = res$p_value, method = "prerank_welch_exact_enumeration", subsets = res$subsets)
  } else {
    fg <- run_fgsea(stats::setNames(as.numeric(ranked_statistic), names(ranked_statistic)), list(set = names(ranked_statistic)[set_index]), 0, as.integer(seed))
    list(set_id = set_id, size = k, es = fg$table$es[1], p_value = fg$table$p_value[1], method = "prerank_welch_fgsea_exploratory", subsets = NA_integer_)
  }
}

