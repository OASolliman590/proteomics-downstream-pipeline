# TMT within-plex loading and bridge normalization (packet R03, SM04, V023).
#
# Loading: subtract from every channel the difference between its observed
# log2 median and the within-plex median of channel medians.
# Bridge: b_fp is the median of a feature's normalized bridge-channel values
# in plex p, b_f the median of available b_fp across plexes; analytical
# values become y - b_fp + b_f.  A feature lacking its bridge measurement in
# a plex is ineligible for cross-plex inference (recorded, never filled).
# Bridge/reference channels are reference-only and never biological units.

.pc_tmt_columns <- function(observations) {
  for (column in c("plex_id", "channel_role")) if (!column %in% names(observations)) stop(sprintf("E_TMT_METADATA: TMT observations need column %s", column), call. = FALSE)
  roles <- observations$channel_role
  if (any(is.na(roles) | !roles %in% c("sample", "bridge", "reference"))) stop("E_TMT_METADATA: channel_role must be sample, bridge or reference", call. = FALSE)
  if (any(is.na(observations$plex_id))) stop("E_TMT_METADATA: plex_id is required for every channel", call. = FALSE)
  invisible(TRUE)
}

tmt_loading <- function(values, observations, observed = NULL) {
  .pc_tmt_columns(observations)
  usable <- .pc_obs_for_centering(values, observed)
  medians <- vapply(seq_len(ncol(values)), function(j) { v <- values[usable[, j], j]; if (!length(v)) NA_real_ else stats::median(v) }, numeric(1))
  if (any(is.na(medians))) stop("E_NORMALIZATION_COVERAGE: a TMT channel has no observed values", call. = FALSE)
  plex <- observations$plex_id
  plex_center <- stats::ave(medians, plex, FUN = stats::median)
  factors <- medians - plex_center
  list(values = sweep(values, 2L, factors, "-"),
       factors = data.frame(observation_id = colnames(values), plex_id = plex, channel_role = observations$channel_role, method = "tmt_loading",
                            center = medians, plex_median = plex_center, factor = factors, stringsAsFactors = FALSE))
}

tmt_bridge <- function(values, observations, bridge_flag) {
  plex <- observations$plex_id
  plexes <- unique(plex)
  missing_bridge <- plexes[!vapply(plexes, function(p) any(bridge_flag[plex == p]), logical(1))]
  if (length(missing_bridge)) stop(sprintf("E_TMT_BRIDGE_REQUIRED: bridge mode needs a bridge channel in every plex; missing in %s", paste(missing_bridge, collapse = ", ")), call. = FALSE)
  b_fp <- sapply(plexes, function(p) { cols <- which(plex == p & bridge_flag); apply(values[, cols, drop = FALSE], 1L, function(v) if (all(is.na(v))) NA_real_ else stats::median(v, na.rm = TRUE)) })
  b_fp <- matrix(b_fp, nrow = nrow(values), dimnames = list(rownames(values), plexes))
  b_f <- apply(b_fp, 1L, function(v) if (all(is.na(v))) NA_real_ else stats::median(v, na.rm = TRUE))
  adjusted <- values
  for (p in plexes) {
    cols <- which(plex == p)
    adjusted[, cols] <- values[, cols, drop = FALSE] - b_fp[, p] + b_f
  }
  eligible <- stats::complete.cases(b_fp)
  eligibility <- data.frame(feature_id = rownames(values), bridge_eligible = eligible,
                            reason = ifelse(eligible, "eligible", "tmt_bridge_measurement_missing_in_plex"), stringsAsFactors = FALSE)
  bridge_table <- data.frame(feature_id = rep(rownames(values), times = length(plexes)), plex_id = rep(plexes, each = nrow(values)),
                             bridge_median = as.vector(b_fp), cross_plex_median = rep(b_f, times = length(plexes)), stringsAsFactors = FALSE)
  list(values = adjusted, eligibility = eligibility, bridge = bridge_table)
}

tmt_no_bridge_check <- function(observations, group_column) {
  analytic <- observations$channel_role == "sample"
  groups_per_plex <- tapply(observations[[group_column]][analytic], observations$plex_id[analytic], function(g) length(unique(g)))
  if (any(groups_per_plex < 2L)) stop(sprintf("E_TMT_PLEX_CONFOUNDED: plex(es) %s contain a single group; plex and treatment are confounded and cannot be corrected away", paste(names(groups_per_plex)[groups_per_plex < 2L], collapse = ", ")), call. = FALSE)
  invisible(groups_per_plex)
}
