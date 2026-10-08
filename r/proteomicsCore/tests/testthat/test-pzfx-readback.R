# R16b read-back (V199-V201, FR-199 to FR-201): the CRAN package pzfx reads the committed .pzfx fixtures written by the
# Python writer (tests/fixtures/figures/pzfx). Oracle: the synthetic source TSVs read here with read.delim, not the writer.
# The TSV and the XML hold the same repr() digits, so as.numeric() of the same string gives identical doubles.
# Skipped (NOT_RUN) when the pzfx package is not installed.

fixtures <- normalizePath(file.path(testthat::test_path(), "..", "..", "..", "..", "tests", "fixtures", "figures", "pzfx"), mustWork = TRUE)
read_source <- function(name) utils::read.delim(file.path(fixtures, name), stringsAsFactors = FALSE, colClasses = "character")

testthat::test_that("V199 column tables read back equal to the source values, one column per group", {
  testthat::skip_if_not_installed("pzfx")
  src <- read_source("column_source.tsv")
  for (k in 1:3) {
    protein <- c("P1", "P2", "P3")[k]
    tab <- pzfx::read_pzfx(file.path(fixtures, "column_groups.pzfx"), table = k)
    testthat::expect_identical(names(tab), c("Control", "Acute"))
    testthat::expect_identical(nrow(tab), 6L)
    for (g in c("Control", "Acute")) {
      want <- as.numeric(src$value[src$protein == protein & src$group == g])
      testthat::expect_identical(as.numeric(tab[[g]]), want)
    }
  }
})

testthat::test_that("V199 the statistics table (D-75) reads back with its comparison rows", {
  testthat::skip_if_not_installed("pzfx")
  st <- pzfx::read_pzfx(file.path(fixtures, "column_groups.pzfx"), table = 4)
  stats <- utils::read.delim(file.path(fixtures, "stats_proteins.tsv"), stringsAsFactors = FALSE)
  testthat::expect_identical(as.character(st$ROWTITLE), stats$comparison)
  testthat::expect_identical(as.numeric(st$effect), as.numeric(stats$effect))
  testthat::expect_identical(as.numeric(st$P), as.numeric(stats$P))
  testthat::expect_identical(as.numeric(st$q), as.numeric(stats$q))
  testthat::expect_identical(as.numeric(st$n_Control), as.numeric(stats$n_Control))
})

testthat::test_that("V200 XY tables read back with one point per source row", {
  testthat::skip_if_not_installed("pzfx")
  volcano <- read_source("volcano_source.tsv")
  roc <- read_source("roc_source.tsv")
  vt <- pzfx::read_pzfx(file.path(fixtures, "xy_volcano.pzfx"), table = 1)
  testthat::expect_identical(names(vt), c("log2FC", "neg_log10P"))
  testthat::expect_identical(nrow(vt), 200L)
  testthat::expect_identical(as.numeric(vt$log2FC), as.numeric(volcano$log2FC))
  testthat::expect_identical(as.numeric(vt$neg_log10P), as.numeric(volcano$neg_log10P))
  rt <- pzfx::read_pzfx(file.path(fixtures, "xy_roc.pzfx"), table = 1)
  testthat::expect_identical(names(rt), c("one_minus_specificity", "sensitivity"))
  testthat::expect_identical(nrow(rt), 30L)
  testthat::expect_identical(as.numeric(rt$sensitivity), as.numeric(roc$sensitivity))
})

testthat::test_that("V201 grouped panel reads back: replicate rows, one data set per protein and group", {
  testthat::skip_if_not_installed("pzfx")
  src <- read_source("grouped_source.tsv")
  tab <- pzfx::read_pzfx(file.path(fixtures, "grouped_panel.pzfx"), table = 1)
  datasets <- setdiff(names(tab), "ROWTITLE")
  proteins <- unique(src$protein)
  groups <- unique(src$group)
  testthat::expect_identical(length(datasets), length(proteins) * length(groups))
  testthat::expect_identical(nrow(tab), 6L)
  testthat::expect_identical(as.character(tab$ROWTITLE), paste0("rep", 1:6))
  for (p in proteins) for (g in groups) {
    want <- as.numeric(src$value[src$protein == p & src$group == g])
    testthat::expect_identical(as.numeric(tab[[paste(p, g)]]), want)
  }
})
