foundation_handlers <- list(
 "foundation.io_roundtrip"="io_roundtrip_impl",
 "preprocessing"="preprocess_stage",
 "design"="design_stage",
 "limma"="limma_stage",
 "assay_engines"="assay_engine_stage",
 "resources"="mapping_stage",
 "pathways"="pathway_stage",
 "response"="response_stage"
)
dispatch_stage <- function(request_path,result_path,requested=NULL) {
  request=jsonlite::fromJSON(request_path,simplifyVector=FALSE); capability=if(is.null(requested))request$capability else requested; handler_name=foundation_handlers[[capability]]
  if(is.null(handler_name)||!exists(handler_name,mode="function")) result=list(schema_version="1.2.0",run_id=request$run_id,stage_id=request$stage_id,capability=request$capability,plan_hash=request$plan_hash,state="NOT_RUN",reason_code="E_CAPABILITY_NOT_IMPLEMENTED",message=paste("no R handler for",capability),started_at=NULL,finished_at=format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC"),exit_code=3L,outputs=list(),warnings=list(),session_info_path=NULL) else result=get(handler_name,mode="function")(request)
  temp_path <- paste0(result_path, ".tmp-", Sys.getpid())
  on.exit(unlink(temp_path, force=TRUE), add=TRUE)
  required_fields <- c("schema_version","run_id","stage_id","capability","plan_hash","state","reason_code","message","started_at","finished_at","exit_code","outputs","warnings","session_info_path")
  if(!all(required_fields %in% names(result))) stop("E_CONFIG_SCHEMA: stage result is missing required fields", call.=FALSE)
  if(!identical(as.character(result$schema_version),"1.2.0") || !result$state %in% c("NOT_RUN","RUNNING","COMPLETED","INAPPLICABLE","FAILED","CANCELLED","NOT_REQUESTED") || !is.numeric(result$exit_code) || length(result$exit_code)!=1L) stop("E_CONFIG_SCHEMA: invalid stage result state or exit_code", call.=FALSE)
  if(identical(result$state,"COMPLETED") && !identical(as.integer(result$exit_code),0L)) stop("E_CONFIG_SCHEMA: completed stage requires exit_code=0", call.=FALSE)
  jsonlite::write_json(result,temp_path,auto_unbox=TRUE,pretty=TRUE,na="null",null="null")
  if (!file.rename(temp_path, result_path)) { unlink(temp_path); stop("E_INTEGRITY: atomic result promotion failed", call.=FALSE) }
  result
}
