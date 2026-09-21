#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE); value=function(name){index=match(name,args);if(is.na(index)||index==length(args))stop(sprintf("missing %s",name),call.=FALSE);args[[index+1L]]}
request_path=value("--request"); result_path=value("--result")
if(!requireNamespace("proteomicsCore",quietly=TRUE)){message("E_CAPABILITY_NOT_AVAILABLE: installed proteomicsCore package is unavailable");quit(status=3L)}
result <- proteomicsCore::dispatch_stage(request_path,result_path)
if(!is.null(result$exit_code) && as.integer(result$exit_code)!=0L) quit(status=as.integer(result$exit_code))
