# Response axes covariance and ratio uncertainty (packet R09, SM19/SM20, V081/V084).

# Per-feature unscaled covariance of contrasts C on the observed rows (exact GLS).
feature_contrast_covariance <- function(y, X, C, weights = NULL, block = NULL, correlation = NULL) {
  use <- !is.na(y)
  Xi <- X[use, , drop = FALSE]
  if (qr(Xi)$rank < ncol(X)) return(NULL)
  n <- sum(use); V <- diag(n)
  if (!is.null(block)) { b <- block[use]; V[outer(b, b, "==")] <- correlation; diag(V) <- 1 }
  w <- if (is.null(weights)) rep(1, n) else weights[use]
  Wh <- diag(sqrt(w), n)
  t(C) %*% solve(t(Xi) %*% Wh %*% solve(V) %*% Wh %*% Xi) %*% C
}

# Fieller confidence set for rho = a / b with Var(a), Var(b), Cov(a, b) and critical value q.
fieller_interval <- function(a, b, var_a, var_b, cov_ab, q) {
  A <- b^2 - q^2 * var_b
  B <- -2 * (a * b - q^2 * cov_ab)
  Cc <- a^2 - q^2 * var_a
  disc <- B^2 - 4 * A * Cc
  if (!all(is.finite(c(A, B, Cc)))) return(list(kind = "unavailable", lower = NA_real_, upper = NA_real_, intervals = list()))
  if (A > 0) {
    if (disc < 0) return(list(kind = "empty", lower = NA_real_, upper = NA_real_, intervals = list()))
    roots <- sort((-B + c(-1, 1) * sqrt(disc)) / (2 * A))
    return(list(kind = "bounded", lower = roots[1], upper = roots[2], intervals = list(c(roots[1], roots[2]))))
  }
  if (A < 0 && disc > 0) {
    roots <- sort((-B + c(-1, 1) * sqrt(disc)) / (2 * A))
    return(list(kind = "disjoint", lower = NA_real_, upper = NA_real_, intervals = list(c(-Inf, roots[1]), c(roots[2], Inf))))
  }
  list(kind = "unbounded", lower = -Inf, upper = Inf, intervals = list(c(-Inf, Inf)))
}

# Frozen response vocabulary (SM20): literal boundaries, no rounding.
response_class <- function(ri) {
  if (is.na(ri)) return("not_evaluable")
  if (ri < 0) return("same_direction")
  if (ri < 0.3) return("limited_return")
  if (ri < 0.8) return("partial_return")
  if (ri <= 1.2) return("near_restoration")
  "overshoot"
}

reversal_index <- function(d, t, dmin) {
  ifelse(is.finite(d) & is.finite(t) & abs(d) >= dmin & d != 0, -t / d, NA_real_)
}
