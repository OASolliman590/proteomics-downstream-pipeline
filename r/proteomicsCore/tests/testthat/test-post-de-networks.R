# R14e network acceptance at the R level (V162, V165).  Oracles are hand formulas for the signed adjacency and TOM,
# prcomp() eigengenes, a hand adjusted Rand index and the measured-universe invariants of the null.

ns <- asNamespace("proteomicsCore")
fn <- function(name) get(name, envir = ns)

testthat::test_that("V162 signed adjacency, TOM, eigengenes and ARI equal hand computations", {
  set.seed(3); X <- matrix(rnorm(200), 40, 5); X[, 2] <- X[, 1] + rnorm(40, 0, 0.2); C <- stats::cor(X)
  A <- fn("pd_signed_adjacency")(C, 6)
  testthat::expect_equal(A[1, 2], ((1 + C[1, 2]) / 2)^6); testthat::expect_true(all(diag(A) == 0))
  T <- fn("pd_signed_tom")(A); k <- rowSums(A)
  testthat::expect_equal(T[1, 3], (sum(A[1, ] * A[, 3]) + A[1, 3]) / (min(k[1], k[3]) + 1 - A[1, 3]), tolerance = 1e-12)
  colnames(X) <- paste0("F", 1:5)
  eg <- fn("pd_eigengenes")(X, c(F1 = 1L, F2 = 1L, F3 = 0L, F4 = 2L, F5 = 2L))
  pc <- stats::prcomp(scale(X[, 1:2]), center = FALSE)$x[, 1]
  testthat::expect_equal(abs(unname(eg$eigengenes[, "M1"])), abs(unname(pc)), tolerance = 1e-10)
  testthat::expect_gt(stats::cor(eg$eigengenes[, "M1"], rowMeans(scale(X[, 1:2]))), 0)    # sign aligned with the module mean
  a <- c(1, 1, 1, 2, 2, 3); b <- c(1, 1, 2, 2, 2, 3)
  tab <- table(a, b); n <- 6; c2 <- function(x) sum(x * (x - 1) / 2)
  e <- c2(rowSums(tab)) * c2(colSums(tab)) / c2(n)
  testthat::expect_equal(fn("pd_adjusted_rand")(a, b), (c2(tab) - e) / (0.5 * (c2(rowSums(tab)) + c2(colSums(tab))) - e))
  testthat::expect_equal(fn("pd_adjusted_rand")(a, a), 1)
})

testthat::test_that("V165 degree-preserving null draws stay in the measured universe and keep degree bins", {
  universe <- sprintf("G%02d", 1:30)
  set.seed(2); edges <- data.frame(a = sample(universe, 80, TRUE), b = sample(universe, 80, TRUE), stringsAsFactors = FALSE); edges <- edges[edges$a != edges$b, ]
  degree <- stats::setNames(vapply(universe, function(g) sum(edges$a == g | edges$b == g), numeric(1)), universe)
  bins <- fn("pd_degree_bins")(degree)
  testthat::expect_true(all(table(bins) >= 5L))
  res <- fn("pd_connectivity_null")(universe[1:6], universe, edges, degree, 50L, 4L)
  testthat::expect_length(res$null, 50L); testthat::expect_identical(res$observed, sum(edges$a %in% universe[1:6] & edges$b %in% universe[1:6]))
  testthat::expect_true(all(is.finite(res$null)) && all(res$null >= 0))
  # a set gene outside the measured universe can never be drawn: the null refuses to leave the universe
  testthat::expect_error(fn("pd_connectivity_null")(c("G01", "X99"), universe, edges, degree, 5L, 4L))
})

testthat::test_that("V164/V165 undirected edges: one row per unordered pair, best score, no self-loops, min_score applied", {
  e <- data.frame(gene_a = c("A", "B", "A", "C", "D", "E"), gene_b = c("B", "A", "A", "D", "C", "F"), score = c(500, 700, 900, 100, 400, 350), stringsAsFactors = FALSE)
  u <- fn("pd_undirected_edges")(e, 300)
  testthat::expect_identical(u$gene_a, c("A", "C", "E")); testthat::expect_identical(u$gene_b, c("B", "D", "F")); testthat::expect_identical(u$score, c(700, 400, 350))
  testthat::expect_identical(nrow(fn("pd_undirected_edges")(e, 1000)), 0L)
})

testthat::test_that("V162 scale-free fit index equals the published definition computed by hand", {
  set.seed(11); k <- stats::rexp(200, 0.2)
  breaks <- seq(min(k), max(k), length.out = 11L); cls <- findInterval(k, breaks, rightmost.closed = TRUE, left.open = TRUE); cls[cls == 0L] <- 1L
  mids <- (breaks[-1] + breaks[-11]) / 2
  dk <- vapply(1:10, function(j) if (any(cls == j)) mean(k[cls == j]) else mids[j], numeric(1))
  pk <- vapply(1:10, function(j) sum(cls == j) / length(k), numeric(1))
  x <- log10(dk); y <- log10(pk + 1e-9); slope <- sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2)
  r2 <- stats::cor(x, y)^2
  got <- fn("pd_scale_free_fit")(k)
  testthat::expect_equal(got$slope, slope, tolerance = 1e-10); testthat::expect_equal(got$r2, -sign(slope) * r2, tolerance = 1e-10)
  testthat::expect_true(is.na(fn("pd_scale_free_fit")(rep(3, 10))$r2))
})

testthat::test_that("V162 bootstrap stability resamples biological units, never features (hand loop oracle)", {
  testthat::skip_if_not_installed("dynamicTreeCut")
  set.seed(21); n <- 40L; f1 <- rnorm(n); f2 <- rnorm(n)
  X <- cbind(sapply(1:8, function(i) f1 + rnorm(n, 0, 0.3)), sapply(1:8, function(i) f2 + rnorm(n, 0, 0.3)), matrix(rnorm(n * 6), n, 6))
  colnames(X) <- sprintf("F%02d", 1:22); rownames(X) <- sprintf("u%02d", 1:n)
  build <- function(Z) fn("pd_modules")(Z, "hclust_correlation", NA, 5L, 2L, 0.5)
  labels <- build(X); units <- rownames(X)
  res <- fn("pd_module_stability")(X, labels, units, "hclust_correlation", NA, 5L, 2L, 0.5, 6L, 17L)
  jaccard <- function(lb, m) { a <- names(labels)[labels == m]; max(c(0, vapply(setdiff(unique(lb), 0L), function(bm) { g <- names(lb)[lb == bm]; length(intersect(a, g)) / length(union(a, g)) }, numeric(1)))) }
  mods <- sort(setdiff(unique(labels), 0L))
  set.seed(17, kind = "L'Ecuyer-CMRG"); unit_oracle <- t(vapply(1:6, function(b) { lb <- build(X[sample(units, n, replace = TRUE), ]); vapply(mods, function(m) jaccard(lb, m), numeric(1)) }, numeric(length(mods))))
  testthat::expect_equal(unname(res$jaccard), unname(unit_oracle), tolerance = 1e-12)
  testthat::expect_true(all(vapply(res$picks, function(p) all(p %in% units) && length(p) == n, logical(1))))   # draws are units
  testthat::expect_identical(res$failed, 0L)
  # a feature-level bootstrap with the same stream is a different procedure and does not reproduce the result
  set.seed(17, kind = "L'Ecuyer-CMRG"); feature_variant <- t(vapply(1:6, function(b) { cols <- sample(colnames(X), ncol(X), replace = TRUE); Z <- X[, cols]; colnames(Z) <- make.unique(cols)
    lb <- build(Z); lb <- lb[!duplicated(cols)]; names(lb) <- cols[!duplicated(cols)]; vapply(mods, function(m) jaccard(lb, m), numeric(1)) }, numeric(length(mods))))
  testthat::expect_false(isTRUE(all.equal(unname(res$jaccard), unname(feature_variant))))
})

testthat::test_that("V163 a single module is tested (one-row eigengene fit) and equals a direct limma fit without trend", {
  set.seed(31); g <- rep(c("A", "B"), each = 12); E <- matrix(stats::rnorm(24) + (g == "B"), 1, 24, dimnames = list("M1", sprintf("o%02d", 1:24)))
  obs <- data.frame(observation_id = colnames(E), group = g, stringsAsFactors = FALSE)
  X <- cbind("(Intercept)" = 1, group.B = as.numeric(g == "B")); rownames(X) <- colnames(E)
  contrasts <- list(list(contrast_id = "B-A", weights = list(0, 1), required_groups = list("A", "B")))
  fam <- list(family_id = "module_trait__group", hypothesis_type = "protein_zero_null", role = "secondary", adjustment = "BH", q_cutoff = 0.05, dependence_assumption = "BH",
              members = list(list(model_id = "module_trait_group", contrast_id = "B-A")))
  settings <- list(group_column = "group", subject_column = NULL, blocking_mode = "none", consensus_correlation = NULL, trend = FALSE, robust = TRUE, ci_level = 0.95)
  res <- fn("pd_fit_model")("module_trait_group", E, matrix(TRUE, 1, 24, dimnames = dimnames(E)), obs, X, contrasts, settings, list(fam), list(run_id = "r", plan_hash = "h"),
                           list(policy = "available_case", minimum_observed_per_group = 2, minimum_fraction = 0.5))
  eb <- limma::eBayes(limma::lmFit(E, X), trend = FALSE, robust = TRUE)
  testthat::expect_identical(res$rows$feature_id, "M1"); testthat::expect_identical(res$rows$eligibility, "tested")
  testthat::expect_equal(res$rows$statistic, unname(eb$t[1, 2]), tolerance = 1e-10); testthat::expect_equal(res$rows$p_value, unname(eb$p.value[1, 2]), tolerance = 1e-12)
  testthat::expect_equal(res$rows$q_value, res$rows$p_value)                                      # one member: BH leaves P unchanged
})

# ----------------------------------------------------------------------------- review 2026-10-05 (minor 5)
testthat::test_that("review minor 5: the degree-preserving null matches exact enumeration on a tiny two-bin graph", {
  # 12 universe genes; G01-G06 have degree 1-2 (bin 1), G07-G12 degree >= 3 (bin 2); exact null: an independent uniform
  # subset of each bin with the set's bin counts, enumerated here by brute force (independent of the sampler).
  universe <- sprintf("G%02d", 1:12)
  edges <- data.frame(a = c("G01", "G02", "G03", "G04", "G05", "G06", "G07", "G07", "G07", "G08", "G08", "G09", "G09", "G10", "G11", "G12", "G12", "G10"),
                      b = c("G07", "G08", "G09", "G10", "G11", "G12", "G08", "G09", "G10", "G09", "G11", "G10", "G12", "G11", "G12", "G01", "G03", "G12"), stringsAsFactors = FALSE)
  degree <- stats::setNames(vapply(universe, function(g) sum(edges$a == g | edges$b == g), numeric(1)), universe)
  bins <- fn("pd_degree_bins")(degree)
  set_genes <- c("G01", "G07", "G08")
  count <- function(g) sum(edges$a %in% g & edges$b %in% g)
  b1 <- names(bins)[bins == bins[["G01"]]]; b2 <- names(bins)[bins == bins[["G07"]]]
  testthat::expect_identical(length(unique(bins)), 2L)
  combos <- expand.grid(i = seq_along(b1), j = seq_len(ncol(utils::combn(b2, 2))))
  pairs2 <- utils::combn(b2, 2)
  exact <- table(factor(apply(combos, 1, function(r) count(c(b1[r[1]], pairs2[, r[2]]))), levels = 0:3)) / nrow(combos)
  res <- fn("pd_connectivity_null")(set_genes, universe, edges, degree, 20000L, 3L)
  mc <- table(factor(res$null, levels = 0:3)) / 20000
  testthat::expect_lt(max(abs(cumsum(mc) - cumsum(exact))), 0.015)                       # Monte Carlo within ~3 SE of the exact null
  testthat::expect_identical(res$observed, count(set_genes))
})

testthat::test_that("review minor 5: each set draws from its own derived L'Ecuyer stream", {
  universe <- sprintf("G%02d", 1:30)
  set.seed(2); edges <- data.frame(a = sample(universe, 80, TRUE), b = sample(universe, 80, TRUE), stringsAsFactors = FALSE); edges <- edges[edges$a != edges$b, ]
  degree <- stats::setNames(vapply(universe, function(g) sum(edges$a == g | edges$b == g), numeric(1)), universe)
  one <- fn("pd_connectivity_null")(universe[1:6], universe, edges, degree, 50L, 4L, stream = 1L)
  two <- fn("pd_connectivity_null")(universe[1:6], universe, edges, degree, 50L, 4L, stream = 2L)
  testthat::expect_false(identical(one$null, two$null))
  testthat::expect_identical(fn("pd_connectivity_null")(universe[1:6], universe, edges, degree, 50L, 4L, stream = 2L)$null, two$null)
  testthat::expect_identical(fn("pd_connectivity_null")(universe[1:6], universe, edges, degree, 50L, 4L)$null, one$null)   # stream 1 (the seed itself) is the default
})
