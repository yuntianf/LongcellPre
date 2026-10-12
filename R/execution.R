# Resolve executable paths before constructing shell pipelines.
tool_path <- function(program) {
  path <- if(file.exists(program)) normalizePath(program, mustWork=TRUE) else Sys.which(program)
  if(!length(path) || !nzchar(path)) stop("Required executable not found: ",program,call.=FALSE)
  unname(path)
}

checked_command <- function(program, args, stdout="", stderr="", label=program) {
  status <- system2(tool_path(program), shQuote(as.character(args)), stdout=stdout, stderr=stderr)
  if(length(status) != 1L || is.na(status) || status != 0L)
    stop(label," failed (exit ",status,"). Check its log and executable version.",call.=FALSE)
  invisible(status)
}

checked_pipeline <- function(commands, log_file, stdout="") {
  bash <- tool_path("bash")
  pipeline <- paste(vapply(commands,function(x) paste(shQuote(x),collapse=" "),character(1)),collapse=" | ")
  status <- system2(bash,c("-o","pipefail","-c",shQuote(pipeline)),stdout=stdout,stderr=log_file)
  if(status != 0L) stop("External pipeline failed (exit ",status,"); see ",log_file,call.=FALSE)
  invisible(NULL)
}

check_samtools <- function(samtools) {
  executable <- tool_path(samtools)
  version <- system2(executable,"--version",stdout=TRUE,stderr=TRUE)
  if(!is.null(attr(version,"status"))) stop("Could not run samtools --version.",call.=FALSE)
  number <- sub("^samtools ([0-9]+\\.[0-9]+).*","\\1",version[1])
  if(!grepl("^[0-9]+\\.[0-9]+$",number) || utils::compareVersion(number,"1.3") < 0)
    stop("samtools >= 1.3 is required (modern sort -o syntax and quickcheck). Found: ",version[1],call.=FALSE)
  executable
}

valid_bam <- function(path, samtools="samtools", require_index=TRUE) {
  if(!file.exists(path) || file.size(path) == 0) return(FALSE)
  if(require_index && (!file.exists(paste0(path,".bai")) || file.size(paste0(path,".bai")) == 0)) return(FALSE)
  status <- suppressWarnings(system2(tool_path(samtools),c("quickcheck",shQuote(path)),stdout=FALSE,stderr=FALSE))
  identical(status,0L)
}

# Publish validated files together, restoring previous complete outputs if a
# rename fails. Temporary files must be on the same filesystem as destinations.
publish_files <- function(staged, destinations) {
  stopifnot(length(staged) == length(destinations),all(file.exists(staged)))
  backup <- paste0(destinations,".longcell-backup-",Sys.getpid())
  if(any(file.exists(backup))) stop("Previous publication backup exists; inspect ",paste(backup,collapse=", "))
  moved <- published <- rep(FALSE,length(staged))
  success <- FALSE
  on.exit({
    if(!success) {
      unlink(destinations[published])
      for(i in which(moved)) file.rename(backup[i],destinations[i])
    } else unlink(backup[moved])
  },add=TRUE)
  for(i in seq_along(staged)) {
    if(file.exists(destinations[i])) {
      if(!file.rename(destinations[i],backup[i])) stop("Cannot stage previous output: ",destinations[i])
      moved[i] <- TRUE
    }
    if(!file.rename(staged[i],destinations[i])) stop("Cannot publish output: ",destinations[i])
    published[i] <- TRUE
  }
  success <- TRUE
  invisible(NULL)
}
