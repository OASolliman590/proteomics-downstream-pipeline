# Exact and Monte Carlo randomization for eligible independent scores (packet R09, SM24, V088).
#
# Two-sided absolute statistic with the frozen tie rule
#   |T_perm| >= |T_obs| - (1e-12 + 1e-10 * |T_obs|).
# Exact: enumerate the complete admissible allocation space once (the observed
# allocation appears once) and report k/N.  Monte Carlo: B independent uniform
# draws with replacement from the full space, plus one forced observed
# statistic: (k+1)/(B+1), with binomial precision of the sampled tail.

randomization_extreme <- function(perm, observed) abs(perm) >= abs(observed) - (1e-12 + 1e-10 * abs(observed))

mean_difference <- function(scores, in_a) mean(scores[in_a]) - mean(scores[!in_a])

exact_label_randomization <- function(scores, group_a, limit = 2e6) {
  n <- length(scores); na <- sum(group_a)
  total <- choose(n, na)
  if (total > limit) stop(sprintf("E_SCORE_EXCHANGEABILITY_UNSUPPORTED: %g allocations exceed the exact enumeration limit; declare independent_monte_carlo", total), call. = FALSE)
  allocations <- utils::combn(n, na)
  observed <- mean_difference(scores, group_a)
  stats <- apply(allocations, 2L, function(idx) { a <- logical(n); a[idx] <- TRUE; mean_difference(scores, a) })
  k <- sum(randomization_extreme(stats, observed))
  observed_count <- sum(apply(allocations, 2L, function(idx) setequal(idx, which(group_a))))
  list(method = "exact", statistic = observed, k = k, N = ncol(allocations), p_value = k / ncol(allocations), observed_allocation_count = observed_count)
}

exact_sign_flip <- function(differences, limit = 2e6) {
  n <- length(differences)
  if (2^n > limit) stop("E_SCORE_EXCHANGEABILITY_UNSUPPORTED: sign-flip space too large for exact enumeration", call. = FALSE)
  signs <- as.matrix(expand.grid(rep(list(c(1, -1)), n)))
  observed <- mean(differences)
  stats <- as.vector(signs %*% differences) / n
  k <- sum(randomization_extreme(stats, observed))
  list(method = "exact_sign_flip", statistic = observed, k = k, N = nrow(signs), p_value = k / nrow(signs), observed_allocation_count = 1L)
}

monte_carlo_label_randomization <- function(scores, group_a, draws, seed) {
  set.seed(seed)
  n <- length(scores); na <- sum(group_a)
  observed <- mean_difference(scores, group_a)
  stats <- vapply(seq_len(draws), function(i) { a <- logical(n); a[sample.int(n, na)] <- TRUE; mean_difference(scores, a) }, numeric(1))
  k <- sum(randomization_extreme(stats, observed))
  ci <- stats::binom.test(k, draws)$conf.int
  list(method = "monte_carlo", statistic = observed, k = k, B = draws, p_value = (k + 1) / (draws + 1), forced_observed = 1L,
       sampled_tail_ci_lower = ci[1], sampled_tail_ci_upper = ci[2], seed = seed)
}
