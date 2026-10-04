# Post-DE biomarker discrimination evaluation (packet R14d, SM34-SM38, FR-149 to FR-160).
#
# Leakage-safe by construction: every learned step (declared imputation, scaling, feature filtering and selection,
# hyperparameter tuning, classifier and decision threshold) is a fit/apply pair fitted on the training portion of a
# fold only and recorded in a per-fold transform audit that is checked (E_BIOMARKER_LEAKAGE).  Folds are stratified and
# grouped by biological unit (subject in subject-blocked designs; E_BIOMARKER_GROUP_LEAKAGE).  Performance is a nested
# cross-validated estimate; the permutation test reruns the whole nested procedure, including selection and tuning,
# on unit-level label permutations.  Every result carries exactly one frozen claim label.  Results are performance
# estimates, never diagnostic or clinical-utility claims.

PD_BM_ALLOWED_TOKENS <- c("cross_validated_nested", "independently_validated", "not externally validated")

# ----------------------------------------------------------------------------- discrimination statistics
pd_auc <- function(score, y) {
  pos <- score[y == 1L]; neg <- score[y == 0L]
  r <- rank(c(pos, neg)); n1 <- length(pos); n0 <- length(neg)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Column-wise Mann-Whitney AUC (higher in the positive class).
pd_auc_columns <- function(X, y) {
  R <- apply(X, 2L, rank); n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  (colSums(R[y == 1L, , drop = FALSE]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# DeLong variance and confidence interval of one AUC (structural components; clamped to [0, 1]).
pd_delong <- function(pos, neg, level = 0.95) {
  psi <- outer(pos, neg, function(a, b) (a > b) + 0.5 * (a == b))
  v10 <- rowMeans(psi); v01 <- colMeans(psi); auc <- mean(psi)
  se <- sqrt(stats::var(v10) / length(pos) + stats::var(v01) / length(neg))
  z <- stats::qnorm(1 - (1 - level) / 2)
  list(auc = auc, se = se, lower = max(0, auc - z * se), upper = min(1, auc + z * se))
}

# Stratified bootstrap percentile CI (controls resampled before cases, as pROC does; quantile type 7).
pd_boot_auc <- function(pos, neg, B, seed, level = 0.95) {
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  aucs <- vapply(seq_len(B), function(b) { controls <- sample(neg, replace = TRUE); cases <- sample(pos, replace = TRUE)
    mean(outer(cases, controls, function(a, c) (a > c) + 0.5 * (a == c))) }, numeric(1))
  q <- stats::quantile(aucs, c((1 - level) / 2, 1 - (1 - level) / 2), names = FALSE)
  list(lower = q[1], upper = q[2], B = B, seed = seed)
}

# ----------------------------------------------------------------------------- folds and permutations on units
# Stratified, unit-grouped folds: within each unit class (sorted: "0", "1", "mixed") units are shuffled and dealt
# cyclically, continuing the deal across classes.
pd_bm_folds <- function(unit_class, k) {
  fold <- integer(length(unit_class)); offset <- 0L
  for (cls in sort(unique(unit_class))) {
    u <- which(unit_class == cls); u <- u[sample.int(length(u))]
    fold[u] <- ((seq_along(u) - 1L + offset) %% k) + 1L; offset <- offset + length(u)
  }
  fold
}

# Unit table: one row per unit with its class ("0", "1" or "mixed" when the class varies within a subject).
.bm_unit_table <- function(units, y) {
  u <- unique(units)
  cls <- vapply(u, function(x) { v <- unique(y[units == x]); if (length(v) > 1L) "mixed" else as.character(v) }, "")
  data.frame(unit = u, class = cls, stringsAsFactors = FALSE)
}

# Label permutation respecting units: class constant within unit -> permute unit labels; varying -> swap within unit.
pd_bm_permute_labels <- function(y, units) {
  ut <- .bm_unit_table(units, y)
  if (all(ut$class != "mixed")) {
    lab <- stats::setNames(as.integer(ut$class), ut$unit)[ut$unit[sample.int(nrow(ut))]]
    names(lab) <- ut$unit
    return(as.integer(lab[units]))
  }
  out <- y
  for (u in ut$unit) { w <- which(units == u); out[w] <- y[w][sample.int(length(w))] }
  out
}

.bm_stream <- function(seed, j) {
  set.seed(seed, kind = "L'Ecuyer-CMRG"); s <- .Random.seed
  for (i in seq_len(j)) s <- parallel::nextRNGStream(s)
  s
}
.bm_use_stream <- function(s) assign(".Random.seed", s, envir = globalenv())

# ----------------------------------------------------------------------------- fit/apply transforms
.bm_fit_prep <- function(X, imputation) {
  med <- if (identical(imputation, "train_median")) apply(X, 2L, stats::median, na.rm = TRUE) else NULL
  if (!is.null(med)) for (j in seq_len(ncol(X))) X[is.na(X[, j]), j] <- med[j]
  keep <- colSums(is.na(X)) == 0
  center <- colMeans(X[, keep, drop = FALSE]); sds <- apply(X[, keep, drop = FALSE], 2L, stats::sd)
  usable <- names(sds)[is.finite(sds) & sds > sqrt(.Machine$double.eps)]
  list(medians = med, features = usable, center = center[usable], scale = sds[usable], fitted_on = rownames(X))   # rows the transforms were fitted on (audit)
}
.bm_apply_prep <- function(prep, X) {
  X <- X[, prep$features, drop = FALSE]
  if (!is.null(prep$medians)) for (f in prep$features) X[is.na(X[, f]), f] <- prep$medians[[f]]
  sweep(sweep(X, 2L, prep$center, "-"), 2L, prep$scale, "/")
}

# Feature ranking/selection fitted on training data only.
.bm_rank <- function(Z, y, method, alpha) {
  if (identical(method, "top_k_auc")) { a <- pd_auc_columns(Z, y); return(names(sort(-abs(a - 0.5), method = "radix"))) }
  if (identical(method, "top_k_t")) {
    eb <- limma::eBayes(limma::lmFit(t(Z), stats::model.matrix(~ factor(y, levels = c(0L, 1L)))))
    return(rownames(eb$t)[order(-abs(eb$t[, 2]), method = "radix")])
  }
  NULL
}

.bm_selector_path <- function(Z, labels, levels, alpha) {
  fit <- glmnet::glmnet(.bm_pad(Z), factor(labels, levels = levels), family = "binomial", alpha = alpha, nlambda = 10, standardize = FALSE)
  list(fit = fit, lambda = fit$lambda)
}
.bm_pad <- function(Z) if (ncol(Z) >= 2L) Z else cbind(Z, .pad = 0)

# Classifiers: scores oriented from the fitted model's class coding (never from evaluation labels).
.bm_fit_classifier <- function(Z, labels, levels, positive, kind, grid) {
  if (identical(kind, "penalized_logistic")) {
    fit <- glmnet::glmnet(.bm_pad(Z), factor(labels, levels = levels), family = "binomial", alpha = 0, lambda = grid, standardize = FALSE)
    return(list(kind = kind, fit = fit, grid = grid, positive = positive, sign = if (identical(fit$classnames[2], positive)) 1 else -1, features = colnames(Z)))
  }
  if (startsWith(kind, "svm_")) {
    kernel <- c(svm_linear = "linear", svm_polynomial = "polynomial", svm_radial = "radial")[[kind]]
    fits <- lapply(grid, function(cost) e1071::svm(Z, factor(labels, levels = levels), kernel = kernel, cost = cost, scale = FALSE))
    return(list(kind = kind, fits = fits, grid = grid, positive = positive, features = colnames(Z)))
  }
  if (identical(kind, "random_forest")) {
    fit <- randomForest::randomForest(Z, factor(labels, levels = levels), ntree = 500)
    return(list(kind = kind, fit = fit, grid = 1, positive = positive, features = colnames(Z)))
  }
  stop(sprintf("E_BIOMARKER_CLASSIFIER: unknown classifier %s", kind), call. = FALSE)
}

# Returns a matrix of positive-class scores (rows = new observations, columns = classifier grid) and whether they are probabilities.
.bm_predict <- function(model, Z) {
  if (identical(model$kind, "penalized_logistic")) {
    p <- stats::predict(model$fit, .bm_pad(Z), s = model$grid, type = "response")
    p <- matrix(p, nrow = nrow(Z)); if (model$sign < 0) p <- 1 - p
    return(list(score = p, probability = TRUE))
  }
  if (startsWith(model$kind, "svm_")) {
    s <- vapply(model$fits, function(f) { pr <- stats::predict(f, Z, decision.values = TRUE); dv <- attr(pr, "decision.values")
      first <- strsplit(colnames(dv)[1], "/")[[1]][1]; as.numeric(dv[, 1]) * if (identical(first, model$positive)) 1 else -1 }, numeric(nrow(Z)))
    return(list(score = matrix(s, nrow = nrow(Z)), probability = FALSE))
  }
  p <- stats::predict(model$fit, Z, type = "prob")[, model$positive]
  list(score = matrix(p, ncol = 1L), probability = TRUE)
}

pd_youden <- function(score, y) {
  u <- sort(unique(score)); cand <- if (length(u) > 1L) c(u[1] - 1, (u[-1] + u[-length(u)]) / 2, u[length(u)] + 1) else c(u - 1, u + 1)
  J <- vapply(cand, function(t) mean(score[y == 1L] >= t) + mean(score[y == 0L] < t) - 1, numeric(1))
  cand[which.max(J)]
}

.bm_threshold <- function(rule, score, y) {
  if (startsWith(rule, "fixed:")) return(as.numeric(sub("^fixed:", "", rule)))
  if (identical(rule, "youden_train")) return(pd_youden(score, y))
  stop(sprintf("E_BIOMARKER_THRESHOLD_LEAKAGE: threshold rule %s is not a training-fold rule", rule), call. = FALSE)
}

# Fit the full pipeline for one configuration on training rows and score new rows.
.bm_fit_score <- function(Xtr, ytr, Xte, s, cfg, ranking = NULL, prep = NULL) {
  prep <- if (is.null(prep)) .bm_fit_prep(Xtr, s$imputation) else prep
  Ztr <- .bm_apply_prep(prep, Xtr); Zte <- .bm_apply_prep(prep, Xte)
  feats <- cfg$features(Ztr, ytr, ranking)
  labels <- ifelse(ytr == 1L, s$positive, s$negative)
  if (!length(feats)) return(list(score = matrix(0.5, nrow(Xte), length(s$classifier_grid)), features = character(), prep = prep, probability = TRUE,
                                  fitted_on = rownames(Ztr), applied_to = rownames(Zte)))
  model <- .bm_fit_classifier(Ztr[, feats, drop = FALSE], labels, s$levels, s$positive, s$classifier, s$classifier_grid)
  pr <- .bm_predict(model, Zte[, feats, drop = FALSE])
  list(score = pr$score, features = feats, prep = prep, probability = pr$probability, model = model,
       fitted_on = rownames(Ztr), applied_to = rownames(Zte))   # rows the selection/classifier saw and the rows scored (audit)
}

# Selection configurations of the grid (each returns the selected features of prepared training data).
.bm_selection_grid <- function(s, Ztr_full = NULL, ytr_full = NULL) {
  if (!is.null(s$fixed_features)) return(list(list(label = "fixed_panel", features = function(Z, y, r) intersect(s$fixed_features, colnames(Z)))))
  if (s$selection %in% c("top_k_auc", "top_k_t"))
    return(lapply(s$k, function(k) list(label = sprintf("k=%d", k), k = k, features = function(Z, y, r) { rr <- if (is.null(r)) .bm_rank(Z, y, s$selection) else r; utils::head(intersect(rr, colnames(Z)), k) })))
  alpha <- if (identical(s$selection, "lasso")) 1 else s$alpha
  lapply(seq_len(10L), function(i) list(label = sprintf("%s lambda index %d", s$selection, i), index = i, features = function(Z, y, r) {
    labels <- ifelse(y == 1L, s$positive, s$negative)
    path <- .bm_selector_path(Z, labels, s$levels, alpha); i2 <- min(i, length(path$lambda))
    b <- as.matrix(stats::coef(path$fit, s = path$lambda[i2]))[-1L, 1]; setdiff(names(b)[b != 0], ".pad") }))
}

# Inner CV on a training set: chooses the selection configuration and classifier parameter (max inner out-of-fold AUC;
# ties -> first in grid order) and the threshold from the inner out-of-fold scores.  Uses only the training rows.
.bm_tune <- function(Xtr, ytr, units_tr, s) {
  if (length(unique(ytr)) < 2L) stop("E_BIOMARKER_FOLDS: a training fold lacks one class", call. = FALSE)
  s$classifier_grid <- if (identical(s$classifier, "penalized_logistic")) {
    if (length(s$lambda)) as.numeric(s$lambda) else { prep0 <- .bm_fit_prep(Xtr, s$imputation); Z0 <- .bm_apply_prep(prep0, Xtr)
      glmnet::glmnet(.bm_pad(Z0), factor(ifelse(ytr == 1L, s$positive, s$negative), levels = s$levels), family = "binomial", alpha = 0, nlambda = 10, standardize = FALSE)$lambda }
  } else if (startsWith(s$classifier, "svm_")) (if (length(s$cost)) as.numeric(s$cost) else c(0.1, 1, 10)) else 1
  grid <- .bm_selection_grid(s)
  ut <- .bm_unit_table(units_tr, ytr)
  inner_fold <- pd_bm_folds(ut$class, s$inner_k)[match(units_tr, ut$unit)]
  inner <- array(NA_real_, c(length(ytr), length(grid), length(s$classifier_grid)))
  for (f in seq_len(s$inner_k)) {
    itr <- which(inner_fold != f); ite <- which(inner_fold == f)
    if (!length(ite) || length(unique(ytr[itr])) < 2L) next
    prep <- .bm_fit_prep(Xtr[itr, , drop = FALSE], s$imputation)
    Z <- .bm_apply_prep(prep, Xtr[itr, , drop = FALSE])
    ranking <- if (s$selection %in% c("top_k_auc", "top_k_t") && is.null(s$fixed_features)) .bm_rank(Z, ytr[itr], s$selection) else NULL
    for (g in seq_along(grid)) inner[ite, g, ] <- .bm_fit_score(Xtr[itr, , drop = FALSE], ytr[itr], Xtr[ite, , drop = FALSE], s, grid[[g]], ranking, prep)$score
  }
  aucs <- matrix(apply(inner, c(2L, 3L), function(v) { ok <- is.finite(v); if (length(unique(ytr[ok])) < 2L) NA_real_ else pd_auc(v[ok], ytr[ok]) }), nrow = length(grid))
  # ties in AUC are broken by the inner Brier score when scores are probabilities (penalised logistic, random forest), then by grid order
  brier <- matrix(apply(inner, c(2L, 3L), function(v) { ok <- is.finite(v); if (!any(ok)) NA_real_ else mean((v[ok] - ytr[ok])^2) }), nrow = length(grid))
  best <- if (all(is.na(aucs))) c(1L, 1L) else {
    cand <- which(abs(aucs - max(aucs, na.rm = TRUE)) <= 1e-12, arr.ind = TRUE)
    if (!startsWith(s$classifier, "svm_") && nrow(cand) > 1L) cand <- cand[order(brier[cand], cand[, 1], cand[, 2]), , drop = FALSE]
    cand[1, ]
  }
  scores <- inner[, best[1], best[2]]; ok <- is.finite(scores)
  list(s = s, grid = grid, best = best, threshold = .bm_threshold(s$threshold_rule, scores[ok], ytr[ok]),
       tuning = sprintf("%s; classifier parameter %s", grid[[best[1]]]$label, format(s$classifier_grid[best[2]])),
       fitted_on = rownames(Xtr), threshold_on = rownames(Xtr)[ok])   # tuning and the threshold use only the training fold's rows (audit)
}

# One outer fold: tune on the training fold, refit the whole pipeline on it, score the held-out fold, record the audit.
.bm_outer_fold <- function(X, y, units, train, test, s, rep_id, fold_id, audit) {
  if (length(intersect(units[train], units[test]))) stop("E_BIOMARKER_GROUP_LEAKAGE: a biological unit has observations in both training and test of a fold", call. = FALSE)
  Xtr <- X[train, , drop = FALSE]; ytr <- y[train]
  tune <- .bm_tune(Xtr, ytr, units[train], s)
  fit <- .bm_fit_score(Xtr, ytr, X[test, , drop = FALSE], tune$s, tune$grid[[tune$best[1]]])
  # review 2026-10-05 minor 2: the audit records the row names of the matrices each step was actually fitted on and applied to
  # (X carries the observation ids as row names), so a step wired to see held-out rows is caught by pd_bm_audit_check
  used <- list(imputation = fit$prep$fitted_on, filter_nonconstant = fit$prep$fitted_on, scaling = fit$prep$fitted_on, selection = fit$fitted_on,
               tuning = tune$fitted_on, classifier = fit$fitted_on, threshold = tune$threshold_on)
  if (!is.null(audit)) for (tr in c(if (identical(s$imputation, "train_median")) "imputation", "filter_nonconstant", "scaling", "selection", "tuning", "classifier", "threshold"))
    audit[[length(audit) + 1L]] <- list(repeat_id = rep_id, fold = fold_id, transform = tr, fitted_on = used[[tr]], applied_to = fit$applied_to,
                                        detail = switch(tr, selection = paste(fit$features, collapse = ";"), tuning = tune$tuning,
                                                        threshold = sprintf("%s = %.6g (inner out-of-fold scores of the training fold)", s$threshold_rule, tune$threshold), ""))
  list(score = fit$score[, tune$best[2]], probability = fit$probability, features = fit$features, threshold = tune$threshold, audit = audit, tuning = tune$tuning)
}

pd_bm_audit_check <- function(audit, folds) {
  for (a in audit) {
    test <- folds[[paste(a$repeat_id, a$fold)]]
    if (length(intersect(unlist(a$fitted_on), test)) || length(intersect(unlist(a$applied_to), unlist(a$fitted_on))))
      stop(sprintf("E_BIOMARKER_LEAKAGE: transform %s of repeat %d fold %d was fitted on held-out samples", a$transform, a$repeat_id, a$fold), call. = FALSE)
  }
  invisible(TRUE)
}

# The whole nested procedure for one label vector (one job = one RNG stream).
pd_bm_nested <- function(X, y, units, s, stream, record = TRUE) {
  if (is.null(rownames(X))) rownames(X) <- s$observation_ids
  if (!identical(rownames(X), as.character(s$observation_ids))) stop("E_INTEGRITY: biomarker matrix rows differ from the declared observation ids", call. = FALSE)
  .bm_use_stream(stream)
  oof <- list(); audit <- if (record) list() else NULL; per_repeat <- numeric(); selected <- list(); folds_test <- list()
  ut <- .bm_unit_table(units, y)
  for (r in seq_len(s$repeats)) {
    fold_unit <- if (identical(s$scheme, "loocv")) seq_len(nrow(ut)) else pd_bm_folds(ut$class, s$outer_k)
    fold <- fold_unit[match(units, ut$unit)]
    scores <- rep(NA_real_, length(y)); thr <- rep(NA_real_, length(y)); prob <- TRUE
    for (k in sort(unique(fold))) {
      test <- which(fold == k); train <- which(fold != k)
      res <- .bm_outer_fold(X, y, units, train, test, s, r, k, audit)
      if (record) audit <- res$audit
      scores[test] <- res$score; thr[test] <- res$threshold; prob <- res$probability
      selected[[length(selected) + 1L]] <- res$features
      folds_test[[paste(r, k)]] <- s$observation_ids[test]
      oof[[length(oof) + 1L]] <- data.frame(repeat_id = r, fold = k, observation_id = s$observation_ids[test], unit_id = units[test], label = y[test],
                                            score = res$score, threshold = res$threshold, predicted = as.integer(res$score >= res$threshold), tuning = res$tuning, stringsAsFactors = FALSE)
    }
    per_repeat <- c(per_repeat, pd_auc(scores, y))
  }
  oof <- do.call(rbind, oof)
  if (record) pd_bm_audit_check(audit, folds_test)
  list(oof = oof, per_repeat = per_repeat, pooled = pd_auc(oof$score, oof$label), selected = selected, audit = audit, probability = prob)
}

# ----------------------------------------------------------------------------- claims
pd_bm_claim_audit <- function(paths) {
  for (path in paths) {
    text <- tolower(paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
    for (token in PD_BM_ALLOWED_TOKENS) text <- gsub(token, "", text, fixed = TRUE)
    if (grepl("\\bdiagnostic\\b|clinical utility|\\bvalidated\\b", text, perl = TRUE))
      stop(sprintf("E_BIOMARKER_CLAIM: %s contains a forbidden claim (diagnostic, clinical utility or validated)", basename(path)), call. = FALSE)
  }
  invisible(TRUE)
}

.bm_locked_hash <- function(locked) {
  m <- locked$model
  classifier <- if (identical(m$kind, "penalized_logistic")) list(kind = m$kind, lambda = locked$parameter,
                                                                  coefficients = as.list(as.matrix(stats::coef(m$fit, s = locked$parameter))[, 1]), sign = m$sign)
                else list(kind = m$kind, parameter = locked$parameter, serialized_sha256 = as.character(openssl::sha256(serialize(if (is.null(m$fit)) m$fits else m$fit, NULL, version = 3))))
  value <- list(features = locked$features, medians = as.list(locked$prep$medians), center = as.list(locked$prep$center[locked$features]),
                scale = as.list(locked$prep$scale[locked$features]), classifier = classifier, threshold = locked$threshold, positive = m$positive)
  json <- as.character(jsonlite::toJSON(value, auto_unbox = TRUE, digits = NA))
  structure(as.character(openssl::sha256(charToRaw(json))), json = json)
}

# Lock the discovery model: tuning by CV on all discovery units, refit on all of them; threshold from the CV scores.
pd_bm_lock <- function(X, y, units, s, stream) {
  .bm_use_stream(stream)
  tune <- .bm_tune(X, y, units, s)
  s2 <- tune$s; s2$classifier_grid <- tune$s$classifier_grid[tune$best[2]]
  fit <- .bm_fit_score(X, y, X, s2, tune$grid[[tune$best[1]]])
  locked <- list(features = fit$features, prep = fit$prep, model = fit$model, parameter = s2$classifier_grid, threshold = tune$threshold, tuning = tune$tuning)
  locked$hash <- .bm_locked_hash(locked)
  locked
}

pd_bm_apply_locked <- function(locked, Xv) {
  missing <- setdiff(locked$features, colnames(Xv))
  if (length(missing)) stop(sprintf("E_VALIDATION_FEATURES: the validation cohort lacks locked features %s", paste(missing, collapse = ", ")), call. = FALSE)
  # review 2026-10-05 minor 1: only the locked features are needed (centring/scaling are per feature), so the
  # preparation is restricted to them and a cohort may lack features the model does not use
  f <- locked$features; prep <- locked$prep
  prep <- list(medians = if (is.null(prep$medians)) NULL else prep$medians[f], features = f, center = prep$center[f], scale = prep$scale[f])
  Z <- .bm_apply_prep(prep, Xv[, f, drop = FALSE])
  as.numeric(.bm_predict(locked$model, Z[, f, drop = FALSE])$score[, 1])
}

# ----------------------------------------------------------------------------- stage handler
post_de_biomarker_stage <- function(request) .pc_run_stage(request, function(out) {
  p <- request$parameters
  outputs <- list(); warnings <- list(); refusals <- list(); figures <- list(); written <- character()
  emit <- function(relative, id, type) outputs[[length(outputs) + 1L]] <<- .pc_output(out, relative, id, type)
  write_tsv <- function(df, relative, id, type) { .pd_check_claims(df); .pc_write_tsv(df, file.path(out, relative)); written <<- c(written, file.path(out, relative)); emit(relative, id, type) }
  write_json <- function(value, relative, id, type) { .pc_write_json(value, file.path(out, relative)); written <<- c(written, file.path(out, relative)); emit(relative, id, type) }
  add_fig <- function(x) { for (f in x$files) emit(f$relative_path, f$artifact_id, "Figure"); figures <<- c(figures, x$records) }
  refuse <- function(analysis, item, code, reason) {
    refusals[[length(refusals) + 1L]] <<- data.frame(analysis = analysis, item = item, reason_code = code, reason = reason, stringsAsFactors = FALSE)
    warnings[[length(warnings) + 1L]] <<- .pc_warning(request, code, sprintf("%s %s: %s", analysis, item, reason), "post_de/biomarker/refusals.tsv")
  }
  for (pkg in c("glmnet")) if (identical(p$classifier, "penalized_logistic") && !requireNamespace(pkg, quietly = TRUE)) stop("E_ENGINE_NOT_AVAILABLE: R package glmnet is required for the penalised-logistic classifier", call. = FALSE)
  if (startsWith(p$classifier, "svm_") && !requireNamespace("e1071", quietly = TRUE)) stop("E_ENGINE_NOT_AVAILABLE: R package e1071 is required for SVM classifiers", call. = FALSE)
  if (identical(p$classifier, "random_forest") && !requireNamespace("randomForest", quietly = TRUE)) stop("E_ENGINE_NOT_AVAILABLE: R package randomForest is required for the random-forest classifier", call. = FALSE)
  plan <- .pd_verify_inputs(request)
  prim <- .pd_primary(request); obs <- prim$obs
  keep <- obs[[p$group_column]] %in% c(p$positive, p$negative)
  obs <- obs[keep, , drop = FALSE]; Yall <- prim$Y[, obs$observation_id, drop = FALSE]
  y <- as.integer(obs[[p$group_column]] == p$positive)
  units <- .pd_units(obs, p$subject_column)
  features <- unlist(p$candidate_features)
  X <- t(Yall[features, , drop = FALSE])
  s <- list(positive = p$positive, negative = p$negative, levels = unlist(p$class_levels), imputation = p$imputation, selection = p$selection$method,
            k = as.integer(unlist(p$selection$k)), alpha = if (is.null(p$selection$alpha)) 0.5 else p$selection$alpha, classifier = p$classifier,
            lambda = unlist(p$classifier_params$lambda), cost = unlist(p$classifier_params$cost), threshold_rule = p$threshold_rule,
            scheme = p$cv$scheme, outer_k = as.integer(p$cv$k), repeats = if (identical(p$cv$scheme, "loocv")) 1L else as.integer(p$cv$repeats), inner_k = as.integer(p$cv$inner_k),
            observation_ids = obs$observation_id, fixed_features = NULL)
  seed <- as.integer(p$seed); level <- 0.95
  result_rows <- list()

  # FR-150 single-feature discrimination on biological units (in-sample; prespecified or in-sample auto direction)
  sf <- p$single_feature
  if (isTRUE(sf$enabled) && identical(sf$state, "ELIGIBLE")) {
    ut <- unique(units)
    unit_y <- vapply(ut, function(u) y[units == u][1], integer(1))
    rows <- list()
    for (f in unlist(sf$features)) {
      v <- vapply(ut, function(u) { x <- Yall[f, units == u]; if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE) }, numeric(1))
      ok <- is.finite(v); pos <- v[ok & unit_y == 1L]; neg <- v[ok & unit_y == 0L]
      if (length(pos) < 2L || length(neg) < 2L) next
      sign_dir <- if (identical(sf$expected_direction, "lower_in_positive")) -1 else 1
      d <- pd_delong(sign_dir * pos, sign_dir * neg, level)
      auto <- identical(sf$direction, "train_only") && d$auc < 0.5
      if (auto) { sign_dir <- -sign_dir; d <- pd_delong(sign_dir * pos, sign_dir * neg, level) }
      b <- if ("bootstrap" %in% unlist(sf$ci)) pd_boot_auc(sign_dir * pos, sign_dir * neg, as.integer(sf$bootstrap), seed, level) else list(lower = NA_real_, upper = NA_real_)
      rows[[length(rows) + 1L]] <- data.frame(feature_id = f, n_positive_units = length(pos), n_negative_units = length(neg), auc = d$auc,
        direction = if (sign_dir > 0) "higher_in_positive" else "lower_in_positive", direction_rule = sf$direction, direction_optimistic = identical(sf$direction, "train_only"),
        delong_lower = d$lower, delong_upper = d$upper, bootstrap_lower = b$lower, bootstrap_upper = b$upper, bootstrap_draws = if ("bootstrap" %in% unlist(sf$ci)) as.integer(sf$bootstrap) else NA_integer_,
        bootstrap_seed = seed, ci_level = level, claim_label = "in_sample", stringsAsFactors = FALSE)
    }
    sf_df <- .pd_rows(rows, c("feature_id", "auc", "claim_label"))
    write_tsv(sf_df, "single_feature_auc.tsv", "post_de_single_feature_auc", "DiscriminationResult")
  } else if (isTRUE(sf$enabled)) refuse("single_feature", "auc", sf$reason_code, sf$reason)

  # FR-152-FR-155 nested cross-validation (job 0 of the RNG streams)
  nested <- pd_bm_nested(X, y, units, s, .bm_stream(seed, 0L), record = TRUE)
  oof <- nested$oof; oof$claim_label <- "cross_validated_nested"; oof$positive_class <- p$positive
  write_tsv(oof, "oof_predictions.tsv", "post_de_oof_predictions", "CvPerformance")
  con <- file(file.path(out, "transform_audit.jsonl"), open = "wb")
  for (a in nested$audit) writeBin(charToRaw(paste0(as.character(jsonlite::toJSON(a, auto_unbox = TRUE)), "\n")), con)
  close(con); emit("transform_audit.jsonl", "post_de_transform_audit", "CvPerformance")
  q <- stats::quantile(nested$per_repeat, c(0.025, 0.975), names = FALSE)
  cv <- data.frame(analysis = "nested_cv", classifier = p$classifier, selection = p$selection$method, scheme = s$scheme, outer_k = s$outer_k, inner_k = s$inner_k, repeats = s$repeats,
                   n_positive_units = length(unique(units[y == 1L])), n_negative_units = length(unique(units[y == 0L])), grouped_by = if (is.null(p$subject_column)) "biological_unit" else p$subject_column,
                   pooled_oof_auc = nested$pooled, repeat_auc_mean = mean(nested$per_repeat), repeat_auc_sd = if (length(nested$per_repeat) > 1L) stats::sd(nested$per_repeat) else NA_real_,
                   repeat_auc_q025 = q[1], repeat_auc_q975 = q[2], seed = seed, rng = "L'Ecuyer-CMRG stream 0", claim_label = "cross_validated_nested", stringsAsFactors = FALSE)
  write_tsv(data.frame(repeat_id = seq_along(nested$per_repeat), auc = nested$per_repeat, claim_label = "cross_validated_nested"), "cv_repeats.tsv", "post_de_cv_repeats", "CvPerformance")

  # FR-154 selection stability over all outer folds
  sel <- table(factor(unlist(nested$selected), levels = features))
  stability <- data.frame(feature_id = features, selection_count = as.integer(sel), outer_folds = length(nested$selected), selection_frequency = as.numeric(sel) / length(nested$selected),
                          claim_label = "cross_validated_nested", stringsAsFactors = FALSE)
  stability <- stability[order(-stability$selection_frequency, stability$feature_id), , drop = FALSE]
  write_tsv(stability, "selection_stability.tsv", "post_de_selection_stability", "SelectionStability")

  # FR-158 calibration and threshold metrics on out-of-fold predictions
  bins <- as.integer(p$calibration_bins)
  if (isTRUE(nested$probability)) {
    b <- pmin(floor(oof$score * bins) + 1L, bins)
    cal <- do.call(rbind, lapply(seq_len(bins), function(i) data.frame(bin = i, lower = (i - 1) / bins, upper = i / bins, n = sum(b == i),
      mean_predicted = if (any(b == i)) mean(oof$score[b == i]) else NA_real_, observed_fraction = if (any(b == i)) mean(oof$label[b == i]) else NA_real_, claim_label = "cross_validated_nested")))
    brier <- mean((oof$score - oof$label)^2)
  } else { cal <- data.frame(bin = integer(), claim_label = character()); brier <- NA_real_; refuse("calibration", p$classifier, "W_CALIBRATION_NOT_PROBABILISTIC", "SVM decision values are not probabilities; Brier score and calibration are not computed") }
  .bm_use_stream(.bm_stream(seed, as.integer(p$permutation$B) + 1L))
  sens <- mean(oof$predicted[oof$label == 1L] == 1L); spec <- mean(oof$predicted[oof$label == 0L] == 0L)
  ut <- .bm_unit_table(units, y); ut <- ut[order(ut$unit, method = "radix"), , drop = FALSE]   # units in sorted order for the stratified bootstrap
  boots <- t(vapply(seq_len(as.integer(sf$bootstrap)), function(i) {
    pick <- unlist(lapply(split(ut$unit, ut$class), function(u) u[sample.int(length(u), length(u), replace = TRUE)]))
    rows <- unlist(lapply(pick, function(u) which(oof$unit_id == u)))
    c(mean(oof$predicted[rows][oof$label[rows] == 1L] == 1L), mean(oof$predicted[rows][oof$label[rows] == 0L] == 0L)) }, numeric(2)))
  cv$brier <- brier; cv$threshold_rule <- p$threshold_rule; cv$sensitivity <- sens; cv$specificity <- spec
  cv$sensitivity_lower <- stats::quantile(boots[, 1], 0.025, names = FALSE, na.rm = TRUE); cv$sensitivity_upper <- stats::quantile(boots[, 1], 0.975, names = FALSE, na.rm = TRUE)
  cv$specificity_lower <- stats::quantile(boots[, 2], 0.025, names = FALSE, na.rm = TRUE); cv$specificity_upper <- stats::quantile(boots[, 2], 0.975, names = FALSE, na.rm = TRUE)
  cv$threshold_bootstrap_draws <- as.integer(sf$bootstrap)
  write_tsv(cal, "calibration.tsv", "post_de_calibration", "CalibrationResult")

  # FR-156 permutation test of the whole nested procedure (jobs 1..B, unit-level label permutations)
  perm <- p$permutation
  null <- NULL
  if (isTRUE(perm$enabled)) {
    jobs <- seq_len(as.integer(perm$B))
    run_job <- function(b) { st <- .bm_stream(seed, b); .bm_use_stream(st); yb <- pd_bm_permute_labels(y, units)
      st2 <- .Random.seed; pd_bm_nested(X, yb, units, s, st2, record = FALSE)$pooled }
    cores <- if (.Platform$OS.type == "unix") max(1L, as.integer(p$threads)) else 1L
    aucs <- unlist(if (cores > 1L) parallel::mclapply(jobs, run_job, mc.cores = cores) else lapply(jobs, run_job))
    k <- sum(aucs >= nested$pooled - 1e-12)
    null <- data.frame(permutation = jobs, pooled_oof_auc = aucs, claim_label = "cross_validated_nested")
    write_tsv(null, "permutation_null.tsv", "post_de_permutation_null", "PermutationResult")
    blocked <- !is.null(p$subject_column)
    cv$permutation_B <- length(jobs); cv$permutation_k <- k; cv$permutation_p <- (k + 1) / (length(jobs) + 1)
    cv$permutation_scheme <- if (!blocked) "biological_units" else if (all(.bm_unit_table(units, y)$class != "mixed")) paste0("whole_subjects:", p$subject_column) else paste0("within:", p$subject_column)
    cv$permutation_scope <- "whole_procedure (selection, tuning, classifier and threshold rerun)"
    cv$compute_estimate_hours <- p$compute_estimate_hours
  } else { cv$permutation_B <- NA_integer_; cv$permutation_k <- NA_integer_; cv$permutation_p <- NA_real_ }
  write_tsv(cv, "cv_performance.tsv", "post_de_cv_performance", "CvPerformance")

  # FR-157 fixed panels: only the classifier (with scaling) in-fold; same-data/unknown panels beside a nested reselection of the same size
  panel_rows <- list(); job <- as.integer(perm$B) + 2L
  for (panel in p$fixed_panels) {
    sp <- s; sp$fixed_features <- intersect(unlist(panel$feature_ids), features)
    fixed <- pd_bm_nested(X, y, units, sp, .bm_stream(seed, job), record = FALSE); job <- job + 1L
    optimistic <- panel$provenance %in% c("same_data", "unknown")
    nested_same <- NULL
    if (optimistic) { sn <- s; sn$selection <- if (s$selection %in% c("top_k_auc", "top_k_t")) s$selection else "top_k_auc"; sn$k <- length(sp$fixed_features)
      nested_same <- pd_bm_nested(X, y, units, sn, .bm_stream(seed, job), record = FALSE); job <- job + 1L }
    panel_rows[[length(panel_rows) + 1L]] <- data.frame(panel_id = panel$id, provenance = panel$provenance, n_features = length(sp$fixed_features),
      panel_sha256 = as.character(openssl::sha256(charToRaw(paste(sort(sp$fixed_features), collapse = "\n")))), pooled_oof_auc = fixed$pooled,
      repeat_auc_mean = mean(fixed$per_repeat), selection_optimistic = optimistic, nested_reselection_auc = if (optimistic) nested_same$pooled else NA_real_,
      note = if (optimistic) "panel chosen on the same (or unknown) data: its selection was not repeated inside the folds; compare with the nested reselection estimate" else "panel declared independent of these data",
      claim_label = "fixed_panel_cv", stringsAsFactors = FALSE)
  }
  write_tsv(.pd_rows(panel_rows, c("panel_id", "claim_label")), "fixed_panels.tsv", "post_de_fixed_panels", "FixedPanelPerformance")

  # FR-159 external validation: lock the discovery model, evaluate once
  vc <- p$validation
  if (!is.null(vc)) {
    locked <- pd_bm_lock(X, y, units, s, .bm_stream(seed, job))
    vm <- .pc_matrix_from_tsv(.pc_find_input(request, "input__validation_matrix"), "numeric")
    vmeta <- .pc_read_tsv(.pc_find_input(request, "input__validation_metadata"))
    id_col <- names(vmeta)[1]
    vmeta <- vmeta[vmeta[[vc$class_column]] %in% c(p$positive, p$negative), , drop = FALSE]
    Xv <- t(vm[, vmeta[[id_col]], drop = FALSE])
    hash_before <- as.character(locked$hash)
    con <- file(file.path(out, "locked_model.json"), open = "wb"); writeBin(charToRaw(attr(locked$hash, "json")), con); close(con)
    emit("locked_model.json", "post_de_locked_model", "ValidationResult")
    sv <- pd_bm_apply_locked(locked, Xv); yv <- as.integer(vmeta[[vc$class_column]] == p$positive)
    if (!identical(as.character(.bm_locked_hash(locked)), hash_before)) stop("E_VALIDATION_RETUNED: the locked model changed during validation", call. = FALSE)
    d <- pd_delong(sv[yv == 1L], sv[yv == 0L], level)
    pred <- as.integer(sv >= locked$threshold)
    validation <- data.frame(state = "evaluated_once", locked_model_sha256 = hash_before, matrix_sha256 = vc$matrix_sha256, metadata_sha256 = vc$metadata_sha256,
      overlap_audit = "no shared subjects, biological units or observations (SM23; checked at plan time)", n_positive = sum(yv == 1L), n_negative = sum(yv == 0L),
      auc = d$auc, delong_lower = d$lower, delong_upper = d$upper, threshold = locked$threshold, sensitivity = mean(pred[yv == 1L] == 1L), specificity = mean(pred[yv == 0L] == 0L),
      locked_features = paste(locked$features, collapse = ";"), claim_label = "independently_validated", stringsAsFactors = FALSE)
  } else validation <- data.frame(state = "not externally validated", locked_model_sha256 = NA_character_, auc = NA_real_,
                                  note = "no independent validation cohort was declared; all estimates are internal (in_sample or cross_validated_nested)", claim_label = "cross_validated_nested", stringsAsFactors = FALSE)
  write_tsv(validation, "validation.tsv", "post_de_validation", "ValidationResult")

  # FR-160 figures with source tables
  fm <- unlist(p$figure_formats)
  roc <- do.call(rbind, lapply(c(Inf, sort(unique(oof$score), decreasing = TRUE)), function(t) data.frame(threshold = t, fpr = mean(oof$score[oof$label == 0L] >= t), tpr = mean(oof$score[oof$label == 1L] >= t))))
  roc$claim_label <- "cross_validated_nested"
  grid <- seq(0, 1, by = 0.05)
  .bm_use_stream(.bm_stream(seed, as.integer(perm$B) + 1000L))
  bandm <- t(vapply(seq_len(min(500L, as.integer(sf$bootstrap))), function(i) {
    pick <- unlist(lapply(split(ut$unit, ut$class), function(u) u[sample.int(length(u), length(u), replace = TRUE)]))
    rows <- unlist(lapply(pick, function(u) which(oof$unit_id == u))); sc <- oof$score[rows]; lb <- oof$label[rows]
    vapply(grid, function(g) { ts <- sort(unique(sc)); fp <- vapply(ts, function(t) mean(sc[lb == 0L] >= t), numeric(1)); tp <- vapply(ts, function(t) mean(sc[lb == 1L] >= t), numeric(1))
      ok <- fp <= g + 1e-12; if (any(ok)) max(tp[ok]) else 0 }, numeric(1)) }, numeric(length(grid))))
  band <- data.frame(fpr = grid, tpr_lower = apply(bandm, 2L, stats::quantile, 0.025, names = FALSE), tpr_upper = apply(bandm, 2L, stats::quantile, 0.975, names = FALSE), claim_label = "cross_validated_nested")
  write_tsv(roc, file.path("figure_sources", "roc_cv.tsv"), "roc_cv_source", "FigureSource")
  write_tsv(band, file.path("figure_sources", "roc_cv_band.tsv"), "roc_cv_band_source", "FigureSource")
  add_fig(.pm_devices(out, "roc_cv", fm, 5.5, 5.5, function() {
    graphics::par(mar = c(4.5, 4.5, 3, 1))
    graphics::plot(roc$fpr, roc$tpr, type = "s", xlim = c(0, 1), ylim = c(0, 1), xlab = "1 - specificity", ylab = "sensitivity", lwd = 2,
                   main = sprintf("Nested-CV out-of-fold ROC (cross_validated_nested), AUC %.2f", nested$pooled), cex.main = 0.8)
    graphics::polygon(c(band$fpr, rev(band$fpr)), c(band$tpr_lower, rev(band$tpr_upper)), col = grDevices::adjustcolor("grey60", 0.4), border = NA)
    graphics::abline(0, 1, lty = 3) }))
  rep_src <- data.frame(repeat_id = seq_along(nested$per_repeat), auc = nested$per_repeat, claim_label = "cross_validated_nested")
  write_tsv(rep_src, file.path("figure_sources", "cv_auc_repeats.tsv"), "cv_auc_repeats_source", "FigureSource")
  add_fig(.pm_devices(out, "cv_auc_repeats", fm, 5.5, 4.5, function() { graphics::par(mar = c(4.5, 4.5, 3, 1))
    graphics::hist(rep_src$auc, breaks = 15, col = "grey70", border = "white", xlim = c(0, 1), xlab = "per-repeat out-of-fold AUC", main = "Nested CV repeats (cross_validated_nested)", cex.main = 0.85) }))
  if (!is.null(null)) {
    write_tsv(null, file.path("figure_sources", "permutation_null.tsv"), "permutation_null_source", "FigureSource")
    add_fig(.pm_devices(out, "permutation_null", fm, 5.5, 4.5, function() { graphics::par(mar = c(4.5, 4.5, 3, 1))
      graphics::hist(null$pooled_oof_auc, breaks = 20, col = "grey70", border = "white", xlim = c(0, 1), xlab = "pooled AUC under permuted labels",
                     main = sprintf("Whole-procedure permutation null, P = %.3g (B = %d)", cv$permutation_p, nrow(null)), cex.main = 0.85)
      graphics::abline(v = nested$pooled, col = "#B2182B", lwd = 2) }))
  }
  if (nrow(cal)) {
    write_tsv(cal, file.path("figure_sources", "calibration.tsv"), "calibration_source", "FigureSource")
    add_fig(.pm_devices(out, "calibration", fm, 5.5, 5.5, function() { graphics::par(mar = c(4.5, 4.5, 3, 1))
      graphics::plot(cal$mean_predicted, cal$observed_fraction, xlim = c(0, 1), ylim = c(0, 1), pch = 19, cex = 0.5 + sqrt(cal$n) / 3, xlab = "mean out-of-fold probability",
                     ylab = "observed fraction positive", main = sprintf("Calibration (cross_validated_nested), Brier %.3f", brier), cex.main = 0.85); graphics::abline(0, 1, lty = 3) }))
  }
  top <- utils::head(stability, 20L)
  write_tsv(top, file.path("figure_sources", "selection_stability.tsv"), "selection_stability_source", "FigureSource")
  add_fig(.pm_devices(out, "selection_stability", fm, 6.5, 4.5, function() { graphics::par(mar = c(6, 4.5, 3, 1))
    graphics::barplot(top$selection_frequency, names.arg = top$feature_id, las = 2, cex.names = 0.6, ylim = c(0, 1), col = "grey40", ylab = "selection frequency",
                      main = "In-fold selection frequency over outer folds", cex.main = 0.85) }))

  refusals_df <- .pd_refusal_frame(refusals)
  write_tsv(cbind(refusals_df, claim_label = rep("descriptive", nrow(refusals_df))), "refusals.tsv", "post_de_biomarker_refusals", "PostDeEligibility")
  write_json(list(module = "biomarker", state = "COMPLETED", reason_code = NULL, claim_label = if (!is.null(vc)) "independently_validated" else "cross_validated_nested",
                  labels = list(single_feature = "in_sample", nested_cv = "cross_validated_nested", fixed_panels = "fixed_panel_cv", validation = if (!is.null(vc)) "independently_validated" else NULL),
                  external_validation = if (is.null(vc)) "not externally validated" else "locked model evaluated once on a disjoint declared cohort",
                  contrast = p$contrast, positive_class = p$positive, negative_class = p$negative, classifier = p$classifier, selection = p$selection, cv = p$cv,
                  seed = seed, rng = "L'Ecuyer-CMRG; stream 0 nested CV, streams 1..B permutations, B+1 threshold bootstrap, then fixed panels and the locked model",
                  permutation = perm, candidate_universe = sprintf("%d features genuinely observed in every analysed observation (label-free planning restriction)", length(features)),
                  eligibility = p$eligibility, refusals = refusals_df, input_hashes = .pd_input_hashes(request), figures = figures,
                  limitations = c("Performance estimates of the analysed biological units; they make no claim about clinical use or deployment.",
                                  "Single-feature AUCs are in-sample; auto-direction AUCs are direction-optimistic.",
                                  "Fixed panels chosen on the same data are selection-optimistic; see the nested reselection estimate.")),
             "biomarker_result.json", "post_de_biomarker_result", "DiscriminationResult")
  pd_bm_claim_audit(setdiff(written, file.path(out, "validation.tsv")))
  if (!identical(validation$claim_label[1], "independently_validated")) pd_bm_claim_audit(file.path(out, "validation.tsv"))
  write_json(list(module = "biomarker", state = "COMPLETED", reason_code = NULL, claim_label = if (!is.null(vc)) "independently_validated" else "cross_validated_nested",
                  eligibility = p$eligibility, refusals = refusals_df, input_hashes = .pd_input_hashes(request), figures = figures,
                  rule = "binary declared contrast; >= min units per class; folds must contain both classes; every learned transform fitted in-fold (audited); permutation reruns the whole procedure"),
             "eligibility.json", "post_de_biomarker_eligibility", "PostDeEligibility")
  list(outputs = outputs, warnings = warnings, message = sprintf("post-DE biomarker: pooled nested-CV AUC %.3f", nested$pooled))
})
