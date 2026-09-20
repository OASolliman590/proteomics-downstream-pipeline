#!/usr/bin/env Rscript
# Capture package versions inside the same process that executes each stage.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: run_stage.R SCRIPT SESSION_OUTPUT")
status <- 0L
tryCatch(source(args[1], local = .GlobalEnv, chdir = FALSE), error = function(e) {
  message(conditionMessage(e))
  status <<- 1L
})
writeLines(capture.output(sessionInfo()), args[2])
quit(status = status)
