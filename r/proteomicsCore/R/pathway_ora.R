# Exact background-aware ORA (packet R08, SM17, V076).  Universe = eligible
# tested measured mapped genes; every eligible set is tested (including k = 0)
# before BH; an empty foreground gives p = 1.
ora_test <- function(universe, foreground, sets, min_size, max_size) {
  N <- length(universe); n <- length(intersect(foreground, universe))
  rows <- lapply(names(sets), function(s) {
    members <- intersect(sets[[s]], universe); K <- length(members)
    if (K < min_size || K > max_size) return(NULL)
    k <- length(intersect(members, foreground))
    p <- if (n == 0L) 1 else stats::phyper(k - 1, K, N - K, n, lower.tail = FALSE)
    data.frame(set_id = s, N = N, K = K, n = n, k = k, p_value = p, overlap_genes = as.character(jsonlite::toJSON(sort(intersect(members, foreground), method = "radix"))),
               status = if (n == 0L) "empty_foreground" else "tested", stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
