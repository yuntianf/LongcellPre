#' Load an installed BSgenome or the same indexed FASTA used for mapping
#' @param genome_name Assembly abbreviation or full BSgenome package name.
#' @param genome_path Optional FASTA; preferred by the pipeline to avoid using
#' a different assembly for mapping and sequence reconstruction.
#' @return A reference object used by sequence extraction.
load_genome <- function(genome_name, genome_path=NULL) {
  if(!is.null(genome_path)) {
    path <- normalizePath(genome_path,mustWork=TRUE)
    if(!file.exists(paste0(path,".fai"))) Rsamtools::indexFa(path)
    index <- Rsamtools::scanFaIndex(path)
    return(structure(list(path=path,seqnames=as.character(GenomicRanges::seqnames(index)),
                          lengths=as.numeric(Biostrings::width(index))),class="longcell_fasta"))
  }
  if(length(genome_name) != 1L || is.na(genome_name) || !nzchar(genome_name))
    stop("Provide genome_path or an installed BSgenome name.")
  package <- genome_name
  if(!startsWith(package,"BSgenome.")) {
    prefixes <- c(hg="Hsapiens",mm="Mmusculus",rn="Rnorvegicus",dm="Dmelanogaster",ce="Celegans",dr="Drerio")
    organism <- prefixes[substr(genome_name,1,2)]
    if(is.na(organism)) stop("Unknown assembly abbreviation; supply the full BSgenome package name or genome_path.")
    package <- paste0("BSgenome.",organism,".UCSC.",genome_name)
  }
  if(!requireNamespace(package,quietly=TRUE))
    stop("Reference package '",package,"' is not installed. Install it with BiocManager::install('",
         package,"'), or use the mapping FASTA via genome_path. Analysis does not install packages automatically.")
  BSgenome::getBSgenome(package)
}

reference_names <- function(genome) {
  if(inherits(genome,"longcell_fasta")) genome$seqnames else names(genome)
}

reference_seq <- function(genome,chr,start,end,strand="+") {
  if(!inherits(genome,"longcell_fasta"))
    return(BSgenome::getSeq(genome,chr,start=start,end=end,strand=strand))
  id <- match(chr,genome$seqnames)
  if(length(id) != 1L || is.na(id) || !is.finite(start) || !is.finite(end) ||
     start < 1 || end < start || end > genome$lengths[id])
    stop("Reference interval outside the mapping FASTA: ",chr,":",start,"-",end)
  ranges <- GenomicRanges::GRanges(chr,IRanges::IRanges(start,end))
  seq <- Rsamtools::scanFa(genome$path,param=ranges)[[1]]
  if(strand == "-") seq <- Biostrings::reverseComplement(seq)
  seq
}
