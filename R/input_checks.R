fastq_format <- function(path) {
  if(!file.exists(path) || is.na(file.size(path)) || file.size(path) == 0)
    stop("FASTQ input is missing or empty: ",path,call.=FALSE)
  con <- file(path,"rb")
  signature <- readBin(con,"raw",n=2)
  close(con)
  compressed <- identical(signature,as.raw(c(31,139)))
  if(grepl("\\.gz$",path,ignore.case=TRUE) && !compressed)
    warning("Input has a .gz suffix but contains plain FASTQ. It will be read as plain data; rename or genuinely compress it before using gzip-only tools.",call.=FALSE)
  con <- if(compressed) gzfile(path,"rt") else file(path,"rt")
  on.exit(close(con),add=TRUE)
  first <- readLines(con,n=1,warn=FALSE)
  if(length(first) != 1L || !startsWith(first,"@")) stop("Input does not begin with a FASTQ read header: ",path)
  invisible(if(compressed) "gzip" else "plain")
}

read_barcodes <- function(path, width) {
  if(!file.exists(path)) stop("Barcode whitelist does not exist: ",path)
  values <- read.table(path,header=FALSE,stringsAsFactors=FALSE,comment.char="",quote="",colClasses="character")
  if(!nrow(values)) stop("Barcode whitelist is empty.")
  values <- toupper(sub("-[0-9]+$","",values[[1]]))
  if(anyNA(values) || any(!grepl("^[ACGT]+$",values)) || any(nchar(values) != width))
    stop("Barcodes must contain exactly ",width," A/C/G/T bases (an optional Cell Ranger -number suffix is accepted). Check the selected chemistry.")
  values <- unique(values)
  message("Barcode whitelist: ",length(values)," entries; barcode length: ",width)
  if(length(values) >= 100000)
    warning("A large chemistry whitelist was supplied. Sample-specific Cell Ranger barcodes can substantially reduce runtime; do not infer cells by arbitrarily truncating the whitelist.",call.=FALSE)
  values
}
