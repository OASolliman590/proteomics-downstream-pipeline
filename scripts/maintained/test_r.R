#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE); suite="foundation"
if("--suite" %in% args){index=match("--suite",args);if(index==length(args))stop("--suite requires a value",call.=FALSE);suite=args[[index+1L]]}
if(!suite %in% c("foundation","unit"))stop("unknown R suite; only foundation/unit are implemented in R01",call.=FALSE)
if(!requireNamespace("testthat",quietly=TRUE)){message("NOT_RUN: testthat is not installed");quit(status=3L)}
if(!requireNamespace("proteomicsCore",quietly=TRUE)){message("NOT_RUN: proteomicsCore is not installed");quit(status=3L)}
result=testthat::test_dir("r/proteomicsCore/tests/testthat",reporter="summary"); if(length(result)==0L)quit(status=3L)
if(any(vapply(result,function(item)identical(item$failed,TRUE),logical(1))))quit(status=1L)
