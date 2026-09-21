testthat::test_that("dispatch map contains only the real foundation handler", {
 testthat::expect_true(is.function(proteomicsCore::dispatch_stage)); testthat::expect_true(is.function(proteomicsCore::io_roundtrip)); testthat::expect_false("limma_stage" %in% ls(asNamespace("proteomicsCore"),all.names=TRUE))
 expected <- c("foundation.io_roundtrip","preprocessing","design","limma","assay_engines","resources","pathways","response")
 testthat::expect_setequal(names(get("foundation_handlers", envir=asNamespace("proteomicsCore"))), expected)
})

testthat::test_that("sha256_file does not leak file connections", {
  root <- file.path(tempdir(), paste0("r01 sha leak &; ", basename(tempfile("isolated-"))))
  if (dir.exists(root)) unlink(root, recursive=TRUE, force=TRUE)
  dir.create(root, recursive=TRUE)
  on.exit(unlink(root, recursive=TRUE, force=TRUE), add=TRUE)
  script <- file.path(root, "leak check Ω.R"); path <- file.path(root, "input Ω &; payload.tsv")
  stdout_path <- file.path(root, "stdout.log"); stderr_path <- file.path(root, "stderr.log")
  writeLines("connection check", path)
  writeLines(c(
    "args <- commandArgs(trailingOnly=TRUE)",
    "if (!requireNamespace('proteomicsCore', quietly=TRUE)) quit(status=3L)",
    "digest <- get('sha256_file', envir=asNamespace('proteomicsCore'))(args[[1L]])",
    "gc()",
    "warning_names <- names(warnings()); if (is.null(warning_names)) warning_names <- character()",
    "if (any(grepl('closing unused connection', warning_names, fixed=TRUE))) { message('R01_CONNECTION_LEAK'); quit(status=17L) }",
    "quit(status=0L)"
  ), script, useBytes=TRUE)
  status <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script), shQuote(path)), stdout=stdout_path, stderr=stderr_path))
  status_attr <- attr(status, "status"); status_code <- if (is.numeric(status) && length(status)==1L) as.integer(status) else if (is.null(status_attr)) 0L else as.integer(status_attr)
  info <- paste(c(readLines(stdout_path, warn=FALSE), readLines(stderr_path, warn=FALSE)), collapse="\n")
  testthat::expect_equal(status_code, 0L, info=info)
})

testthat::test_that("foundation I/O roundtrip preserves NA and UTF-8", {
  testthat::expect_true(requireNamespace("jsonlite", quietly=TRUE)); testthat::expect_true(requireNamespace("openssl", quietly=TRUE)); testthat::expect_true(requireNamespace("proteomicsCore", quietly=TRUE))
  root <- file.path(tempdir(), "roundtrip space Ω & ; shell-safe"); if (dir.exists(root)) unlink(root, recursive=TRUE, force=TRUE); dir.create(root, recursive=TRUE); on.exit(unlink(root, recursive=TRUE, force=TRUE), add=TRUE)
  outside <- file.path(dirname(root), "R01-outside-sentinel"); if (file.exists(outside)) unlink(outside, force=TRUE)
  matrix_path <- file.path(root, "matrix.tsv"); metadata_path <- file.path(root, "metadata.tsv")
  writeLines(c("feature_id\tobs-1\tobs-2", "0001\t1.25\tNA", "0002\t2.5\t3.0"), matrix_path, useBytes=TRUE)
  writeLines(c("observation_id\tlabel\tempty_field", "obs-1\tΩ specimen\t", "obs-2\tNA\tvalue"), metadata_path, useBytes=TRUE)
  digest <- function(path) { connection <- file(path, "rb"); on.exit(close(connection), add=TRUE); value <- paste0(as.character(openssl::sha256(connection)), collapse=""); testthat::expect_type(value, "character"); testthat::expect_length(value, 1L); testthat::expect_false(is.object(value)); value }
  matrix_digest <- digest(matrix_path); metadata_digest <- digest(metadata_path)
  testthat::expect_type(matrix_digest, "character"); testthat::expect_type(metadata_digest, "character")
  request <- list(schema_version="1.2.0", run_id="r", stage_id="s", capability="foundation.io_roundtrip", plan_hash=NULL,
    inputs=list(list(artifact_id="matrix",path=matrix_path,sha256=matrix_digest), list(artifact_id="metadata",path=metadata_path,sha256=metadata_digest)),
    output_temp_dir=file.path(root,"out"), config_path=NULL, parameters=list(matrix_input_id="matrix",metadata_input_id="metadata"),
    rng=list(seed=1L,kind="L'Ecuyer-CMRG",threads=1L))
  request_path <- file.path(root,"request.json"); result_path <- file.path(root,"result.json")
  jsonlite::write_json(request, request_path, auto_unbox=TRUE, na="null")
  proteomicsCore::io_roundtrip(request_path, result_path)
  result <- jsonlite::fromJSON(result_path, simplifyVector=FALSE)
  testthat::expect_identical(result$state, "COMPLETED")
  testthat::expect_true(file.exists(file.path(root,"out","sessionInfo.txt")))
  output_matrix <- read.delim(file.path(root,"out","matrix.tsv"),header=TRUE,sep="\t",quote="",comment.char="",check.names=FALSE,na.strings="NA",colClasses="character",stringsAsFactors=FALSE)
  output_metadata <- read.delim(file.path(root,"out","metadata.tsv"),header=TRUE,sep="\t",quote="",comment.char="",check.names=FALSE,na.strings=character(),colClasses="character",stringsAsFactors=FALSE,fill=TRUE)
  testthat::expect_identical(as.character(output_matrix[[1L]]),c("0001","0002")); testthat::expect_equal(as.numeric(output_matrix[[2L]]),c(1.25,2.5)); testthat::expect_true(is.na(output_matrix[[3L]][1L])); testthat::expect_equal(as.numeric(output_matrix[[3L]][2L]),3.0)
  testthat::expect_identical(as.character(output_metadata[[1L]]),c("obs-1","obs-2")); testthat::expect_identical(output_metadata$label[1L],"Ω specimen"); testthat::expect_identical(output_metadata$label[2L],"NA"); testthat::expect_identical(output_metadata$empty_field[1L],""); testthat::expect_identical(output_metadata$empty_field[2L],"value")
  testthat::expect_true(any(vapply(result$outputs,function(item)identical(item$relative_path,"sessionInfo.txt"),logical(1))))
  testthat::expect_false(file.exists(outside))
})

testthat::test_that("an intentionally failing isolated R harness is nonzero", {
  harness_source <- normalizePath(file.path(testthat::test_path(), "..", "..", "..", "..", "scripts", "maintained", "test_r.R"), mustWork=TRUE)
  root <- file.path(tempdir(), paste0("r01 harness space &; ", basename(tempfile("isolated-"))))
  if (dir.exists(root)) unlink(root, recursive=TRUE, force=TRUE)
  dir.create(file.path(root, "scripts", "maintained"), recursive=TRUE)
  dir.create(file.path(root, "r", "proteomicsCore", "tests", "testthat"), recursive=TRUE)
  on.exit(unlink(root, recursive=TRUE, force=TRUE), add=TRUE)
  harness <- file.path(root, "scripts", "maintained", "test_r.R")
  test_file <- file.path(root, "r", "proteomicsCore", "tests", "testthat", "test-intentional-failure.R")
  stdout_path <- file.path(root, "stdout.log"); stderr_path <- file.path(root, "stderr.log")
  testthat::expect_true(file.copy(harness_source, harness, overwrite=TRUE))
  writeLines(c(
    "testthat::test_that('R01 intentional harness failure', {",
    "  testthat::expect_true(FALSE, info='R01_HARNESS_ASSERTION_FAILURE_9f3a')",
    "})"
  ), test_file, useBytes=TRUE)
  oldwd <- getwd(); on.exit(setwd(oldwd), add=TRUE); setwd(root)
  status <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(harness), "--suite", "foundation"), stdout=stdout_path, stderr=stderr_path))
  setwd(oldwd)
  status_attr <- attr(status, "status"); status_code <- if (is.numeric(status) && length(status)==1L) as.integer(status) else if (is.null(status_attr)) 0L else as.integer(status_attr)
  info <- paste(c(readLines(stdout_path, warn=FALSE), readLines(stderr_path, warn=FALSE)), collapse="\n")
  testthat::expect_true(status_code != 0L, info=info)
  testthat::expect_true(grepl("R01_HARNESS_ASSERTION_FAILURE_9f3a", info, fixed=TRUE), info=info)
})

testthat::test_that("unknown R capabilities remain NOT_RUN", {
  root <- tempfile("dispatch-"); dir.create(root)
  request <- list(schema_version="1.2.0",run_id="r",stage_id="s",capability="limma",plan_hash=NULL,inputs=list(),output_temp_dir=root,config_path=NULL,parameters=list(),rng=list(seed=1L,kind="L'Ecuyer-CMRG",threads=1L))
  request_path <- file.path(root,"request.json"); result_path <- file.path(root,"result.json")
  jsonlite::write_json(request, request_path, auto_unbox=TRUE, na="null", null="null")
  proteomicsCore::dispatch_stage(request_path,result_path)
  raw_result <- paste(readLines(result_path, warn=FALSE), collapse="\n"); testthat::expect_true(grepl('"plan_hash": null', raw_result, fixed=TRUE))
  unknown_result <- jsonlite::fromJSON(result_path); testthat::expect_identical(unknown_result$state, "NOT_RUN"); testthat::expect_true(is.null(unknown_result$plan_hash))
})

testthat::test_that("io_roundtrip rejects extra parameters and duplicate input IDs", {
  root <- tempfile("negative-io-"); dir.create(root); request_path <- file.path(root,"request.json"); result_path <- file.path(root,"result.json")
  request <- list(schema_version="1.2.0",run_id="r",stage_id="s",capability="foundation.io_roundtrip",plan_hash=NULL,inputs=list(list(artifact_id="same",path="a",sha256=paste(rep("0",64),collapse="")),list(artifact_id="same",path="b",sha256=paste(rep("0",64),collapse=""))),output_temp_dir=root,config_path=NULL,parameters=list(matrix_input_id="same",metadata_input_id="same",extra="reject"),rng=list(seed=1L,kind="L'Ecuyer-CMRG",threads=1L))
  jsonlite::write_json(request,request_path,auto_unbox=TRUE,na="null",null="null")
  testthat::expect_error(proteomicsCore::io_roundtrip(request_path,result_path),"E_CONFIG_SCHEMA")
  request$parameters <- list(matrix_input_id="same",metadata_input_id="same"); jsonlite::write_json(request,request_path,auto_unbox=TRUE,na="null",null="null")
  testthat::expect_error(proteomicsCore::io_roundtrip(request_path,result_path),"E_ID_DUPLICATE")
})
