# R07 acceptance at the R level (V064, V065, V066, V068, V069).
ns <- asNamespace("proteomicsCore"); fn <- function(name) get(name, envir = ns)

mapped_fixture <- function() data.frame(feature_id = c("f3", "f1", "f2", "g1", "g2", "h2", "h1"), gene_id = c("A", "A", "A", "B", "B", "C", "C"),
                                        mapping_state = "mapped", stringsAsFactors = FALSE)

testthat::test_that("V065 representative tie stages: coverage, then median, then feature_id", {
  obs <- paste0("o", 1:4)
  values <- rbind(f3 = c(10, 10, 10, NA), f1 = c(5, 5, 5, 5), f2 = c(1, 1, 1, 1),   # A: f1/f2 full coverage; f1 higher median
                  g1 = c(4, 4, 4, 4), g2 = c(6, 6, 6, 6),                         # B: equal coverage; g2 higher median
                  h2 = c(3, 3, 3, 3), h1 = c(3, 3, 3, 3))                         # C: full tie -> h1 (smaller id)
  colnames(values) <- obs
  observed <- !is.na(values)
  r <- fn("select_representatives")(mapped_fixture(), values, observed)
  winners <- r$feature_id[r$representative_state == "representative"]
  testthat::expect_setequal(winners, c("f1", "g2", "h1"))
  for (perm in list(c(4, 3, 2, 1), c(2, 4, 1, 3))) {     # column (label/order) permutation leaves winners unchanged
    r2 <- fn("select_representatives")(mapped_fixture(), values[, perm], observed[, perm])
    testthat::expect_setequal(r2$feature_id[r2$representative_state == "representative"], winners)
  }
})

testthat::test_that("V066 multi-gene groups are excluded; same-gene groups retained", {
  features <- data.frame(feature_id = c("u", "same", "multi"), accessions = c('["a1"]', '["b1","b2"]', '["c1"]'), stringsAsFactors = FALSE)
  mapping <- data.frame(source_id = c("a1", "b1", "b2", "c1", "c1"), gene_id = c("GA", "GB", "GB", "GC", "GD"), gene_symbol = "x", taxonomy_id = "1", status = "current", stringsAsFactors = FALSE)
  m <- fn("map_features_to_genes")(features, mapping, 1)
  testthat::expect_identical(m$mapping_reason, c("unique_gene", "unique_gene", "multi_gene_group"))
  testthat::expect_true(is.na(m$gene_id[3]))
})

testthat::test_that("V068 set eligibility uses set ∩ universe with declared bounds", {
  universe <- paste0("g", 1:12)
  sets <- rbind(data.frame(set_id = "s0", set_name = "s0", gene_id = c("x1", "x2"), stringsAsFactors = FALSE),
                data.frame(set_id = "s1", set_name = "s1", gene_id = c("g1", "x3"), stringsAsFactors = FALSE),
                data.frame(set_id = "s3", set_name = "s3", gene_id = c("g1", "g2", "g3", "x4"), stringsAsFactors = FALSE),
                data.frame(set_id = "s6", set_name = "s6", gene_id = paste0("g", 1:6), stringsAsFactors = FALSE))
  sets$source_gene_id <- sets$gene_id; sets$projection_status <- "native"
  el <- fn("gene_set_eligibility")(sets, universe, 2, 5)$summary
  testthat::expect_identical(el$eligible_overlap_size, c(0L, 1L, 3L, 6L))
  testthat::expect_identical(el$set_id[el$eligible], "s3")
})

testthat::test_that("V069 ORA universe is the eligible measured mapped genes", {
  gm <- data.frame(feature_id = paste0("p", 1:8), gene_id = paste0("G", 1:8), representative_state = "representative", stringsAsFactors = FALSE)
  testthat::expect_identical(fn("ora_universe")(gm, paste0("p", 1:6)), paste0("G", 1:6))
})

testthat::test_that("V064 projection keeps evidence and excludes ambiguous orthologs", {
  sets <- data.frame(set_id = "S", set_name = "S", gene_id = c("h1", "h2"), stringsAsFactors = FALSE)
  orth <- data.frame(source_gene_id = c("h1", "h2", "h2"), source_taxonomy_id = "9606", target_gene_id = c("r1", "r2", "r3"), target_taxonomy_id = "10116",
                     evidence = c("e1", "e2", "e3"), ambiguous = c("false", "true", "true"), stringsAsFactors = FALSE)
  p <- fn("project_gene_sets")(sets, orth, 10116)
  testthat::expect_identical(p$projection_status, c("projected", "ambiguous_ortholog_excluded"))
  testthat::expect_identical(p$evidence[1], "e1"); testthat::expect_true(is.na(p$gene_id[2]))
})
