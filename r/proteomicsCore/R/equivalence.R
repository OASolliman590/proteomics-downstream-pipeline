# Residual equivalence (TOST) and formal rescue conjunction (packet R09, SM21/SM22, V085/V086).

tost <- function(r, se, df, delta, alpha) {
  if (!is.finite(delta) || delta <= 0) stop("E_EQUIVALENCE_INELIGIBLE: equivalence margin must be positive", call. = FALSE)
  t_lower <- (r + delta) / se; t_upper <- (r - delta) / se
  p_lower <- stats::pt(t_lower, df, lower.tail = FALSE)
  p_upper <- stats::pt(t_upper, df, lower.tail = TRUE)
  half <- stats::qt(1 - alpha, df) * se
  data.frame(t_lower = t_lower, t_upper = t_upper, p_lower = p_lower, p_upper = p_upper, p_value = pmax(p_lower, p_upper),
             ci_lower = r - half, ci_upper = r + half, ci_level = 1 - 2 * alpha)
}

# One-sided component tests for a frozen direction s: H0 s*d <= m_d, H0 -s*t <= m_t, plus residual TOST.
rescue_components <- function(s, d, se_d, t, se_t, r, se_r, df, disease_margin, treatment_margin, residual_margin) {
  p_disease <- stats::pt((s * d - disease_margin) / se_d, df, lower.tail = FALSE)
  p_treatment <- stats::pt((-s * t - treatment_margin) / se_t, df, lower.tail = FALSE)
  eq <- tost(r, se_r, df, residual_margin, 0.05)
  data.frame(p_disease = p_disease, p_treatment = p_treatment, p_residual_lower = eq$p_lower, p_residual_upper = eq$p_upper,
             p_value = pmax(p_disease, p_treatment, eq$p_lower, eq$p_upper))
}

intersection_union <- function(components) apply(as.matrix(components), 1L, max)
