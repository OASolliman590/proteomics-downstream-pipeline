sha256_file <- function(path) {
  if (requireNamespace("openssl", quietly = TRUE)) {
    connection <- file(path, "rb"); on.exit(close(connection), add=TRUE)
    digest <- paste0(as.character(openssl::sha256(connection)), collapse="")
    if (!is.character(digest) || length(digest) != 1L || is.object(digest)) stop("E_INTEGRITY: openssl returned a non-scalar digest", call.=FALSE)
    return(digest)
  }
  stop("E_CAPABILITY_NOT_AVAILABLE: R package openssl is required for artifact hashes", call. = FALSE)
}
read_table_preserve <- function(path) read.delim(path, header=TRUE, sep="\t", quote="", comment.char="", check.names=FALSE, na.strings=character(), stringsAsFactors=FALSE, colClasses="character", fill=TRUE, blank.lines.skip=FALSE, strip.white=FALSE)
write_table_preserve <- function(data,path) write.table(data,file=path,sep="\t",quote=FALSE,row.names=FALSE,na="NA",fileEncoding="UTF-8",qmethod="double")
io_roundtrip_impl <- function(request) {
  params=request$parameters
  expected_parameters <- c("matrix_input_id","metadata_input_id")
  if (!is.list(params) || !identical(sort(names(params)),sort(expected_parameters)) || !all(vapply(params,function(x)is.character(x)&&length(x)==1L&&nzchar(x),logical(1)))) stop("E_CONFIG_SCHEMA: io_roundtrip parameters must contain exactly matrix_input_id and metadata_input_id strings",call.=FALSE)
  inputs=request$inputs
  input_ids=vapply(inputs,function(item)as.character(item$artifact_id),character(1)); if(anyDuplicated(input_ids)) stop("E_ID_DUPLICATE: input artifact IDs must be unique",call.=FALSE)
  find_input <- function(id) { matches=vapply(inputs,function(item)identical(item$artifact_id,id),logical(1)); if(sum(matches)!=1L)stop(sprintf("E_REFERENCE_UNKNOWN: input artifact %s",id),call.=FALSE); inputs[[which(matches)[1L]]] }
  matrix_input=find_input(params$matrix_input_id); metadata_input=find_input(params$metadata_input_id)
  if (!identical(sha256_file(matrix_input$path), matrix_input$sha256)) stop("E_INTEGRITY: matrix input hash mismatch", call.=FALSE)
  if (!identical(sha256_file(metadata_input$path), metadata_input$sha256)) stop("E_INTEGRITY: metadata input hash mismatch", call.=FALSE)
  warning_messages <- character()
  capture <- function(expr) withCallingHandlers(expr, warning=function(w) { warning_messages <<- c(warning_messages, conditionMessage(w)); invokeRestart("muffleWarning") })
  matrix_data=capture(read_table_preserve(matrix_input$path)); metadata_data=capture(read_table_preserve(metadata_input$path))
  if(ncol(matrix_data)<2L || ncol(metadata_data)<1L) stop("E_CONFIG_SCHEMA: matrix and metadata tables require columns",call.=FALSE)
  matrix_data[-1L] <- lapply(matrix_data[-1L], function(column) { missing <- is.na(column) | column %in% c("NA",""); numeric <- suppressWarnings(as.numeric(column)); if(any(is.na(numeric) & !missing)) stop("E_CONFIG_SCHEMA: matrix values must be numeric",call.=FALSE); numeric[missing] <- NA_real_; numeric })
  if(ncol(metadata_data)>1L) metadata_data[-1L] <- lapply(metadata_data[-1L], function(column) { column[column=="NA"] <- NA_character_; column })
  if(anyDuplicated(matrix_data[[1L]])||anyDuplicated(metadata_data[[1L]]))stop("E_ID_DUPLICATE: matrix and metadata keys must be unique",call.=FALSE)
  if (!identical(colnames(matrix_data)[-1L], as.character(metadata_data[[1L]]))) stop("E_ID_ALIGNMENT: matrix columns and metadata keys do not align", call.=FALSE)
  output_dir=normalizePath(request$output_temp_dir,mustWork=FALSE); dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  matrix_path=file.path(output_dir,"matrix.tsv"); metadata_path=file.path(output_dir,"metadata.tsv"); capture(write_table_preserve(matrix_data,matrix_path)); capture(write_table_preserve(metadata_data,metadata_path))
  session_path=file.path(output_dir,"sessionInfo.txt"); capture.output(sessionInfo(),file=session_path)
  timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz="UTC")
  warning_records <- lapply(warning_messages, function(message) list(code="R_WARNING",stage_id=request$stage_id,severity="warning",message=message,artifact_ref=NULL))
  list(schema_version="1.2.0",run_id=request$run_id,stage_id=request$stage_id,capability=request$capability,plan_hash=request$plan_hash,state="COMPLETED",reason_code=NULL,message="foundation I/O roundtrip completed",started_at=timestamp,finished_at=timestamp,exit_code=0L,outputs=list(list(artifact_id="matrix",relative_path="matrix.tsv",sha256=sha256_file(matrix_path),result_type="matrix_tsv"),list(artifact_id="metadata",relative_path="metadata.tsv",sha256=sha256_file(metadata_path),result_type="metadata_tsv"),list(artifact_id="session_info",relative_path="sessionInfo.txt",sha256=sha256_file(session_path),result_type="session_info")),warnings=warning_records,session_info_path="sessionInfo.txt")
}
io_roundtrip <- function(request_path,result_path) dispatch_stage(request_path,result_path,requested="foundation.io_roundtrip")
