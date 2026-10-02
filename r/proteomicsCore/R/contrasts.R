# Exact contrast strategy (packet R04, SM08, V035/V037).
#
# The general exact path reparameterizes the design so that the contrast is
# the first coefficient (X* = X T^-1 with T's first row equal to c') and
# refits with identical weights/block/correlation.  Its standard error is
# exactly c'(X'V^-1X)^-1 c times the feature residual variance even with
# missing rows, weights or blocking.  limma::contrasts.fit on a full-design
# fit is retained only as a labelled approximate shortcut.

reparameterize_contrast <- function(X, contrast) {
  contrast <- as.numeric(contrast)
  if (length(contrast) != ncol(X)) stop("E_CONTRAST_COEFFICIENT: contrast length does not match design columns", call. = FALSE)
  if (all(contrast == 0)) stop("E_CONTRAST_COEFFICIENT: contrast weights are all zero", call. = FALSE)
  q <- qr(matrix(contrast, ncol = 1L))
  basis <- qr.Q(q, complete = TRUE)
  transform <- rbind(contrast, if (ncol(X) > 1L) t(basis[, -1L, drop = FALSE]) else NULL)
  X_star <- X %*% solve(transform)
  colnames(X_star) <- c("contrast", if (ncol(X) > 1L) paste0("nuisance_", seq_len(ncol(X) - 1L)) else NULL)
  rownames(X_star) <- rownames(X)
  list(X = X_star, transform = transform)
}

.pc_lmfit <- function(Y, X, weights = NULL, block = NULL, correlation = NULL) {
  args <- list(object = Y, design = X)
  if (!is.null(weights)) args$weights <- weights
  if (!is.null(block)) { args$block <- block; args$correlation <- correlation }
  suppressWarnings(do.call(limma::lmFit, args))
}

exact_contrast_fit <- function(Y, X, contrast, weights = NULL, block = NULL, correlation = NULL) {
  if (!requireNamespace("limma", quietly = TRUE)) stop("E_DEPENDENCY_UNAVAILABLE: limma is required for the exact contrast fit", call. = FALSE)
  if (!is.null(block) && is.null(correlation)) stop("E_BLOCKING_CONFLICT: a block requires the frozen consensus correlation", call. = FALSE)
  r <- reparameterize_contrast(X, contrast)
  fit <- .pc_lmfit(Y, r$X, weights, block, correlation)
  list(fit = fit, transform = r$transform, path = "contrast_as_coefficient_refit")
}

approximate_contrast_fit <- function(Y, X, contrast_matrix, weights = NULL, block = NULL, correlation = NULL) {
  fit <- .pc_lmfit(Y, X, weights, block, correlation)
  list(fit = limma::contrasts.fit(fit, contrast_matrix), path = "contrasts_fit_shortcut_approximate")
}

# Unscaled covariance among contrasts for the complete design: C (X' W V^-1 W X)^-1 C'.
contrast_unscaled_covariance <- function(X, contrasts, weights = NULL, block = NULL, correlation = NULL) {
  n <- nrow(X)
  V <- diag(n)
  if (!is.null(block)) { same <- outer(block, block, "=="); V[same] <- correlation; diag(V) <- 1 }
  root_w <- if (is.null(weights)) rep(1, n) else sqrt(weights)
  Wh <- diag(root_w, n)
  information <- t(X) %*% Wh %*% solve(V) %*% Wh %*% X
  C <- as.matrix(contrasts)
  t(C) %*% solve(information) %*% C
}
