#' @title load_genome
#' @description Load the genome from BSgenome given the genome name
#' @param genome_name The name of the genome, can be abbreviations for some commonly used genomes like "hg38"
#' @importFrom BSgenome available.genomes
#' @importFrom BSgenome getBSgenome
#' @importFrom BiocManager install
#' @return A BSgenome object
load_genome <- function(genome_name) {
  genome_list = BSgenome::available.genomes()

  if(genome_name %in% genome_list){
    genome_package = genome_name
  }
  # Handle prefix dynamically based on the genome_name and organism
  else {
    # Infer organism prefix from genome_name (a simplified heuristic)
    if (startsWith(genome_name, "hg")) {
      organism_prefix <- "BSgenome.Hsapiens.UCSC."
    } else if (startsWith(genome_name, "mm")) {
      organism_prefix <- "BSgenome.Mmusculus.UCSC."
    } else if (startsWith(genome_name, "rn")) {
      organism_prefix <- "BSgenome.Rnorvegicus.UCSC."
    } else if (startsWith(genome_name, "dm")) {
      organism_prefix <- "BSgenome.Dmelanogaster.UCSC."
    } else if (startsWith(genome_name, "ce")) {
      organism_prefix <- "BSgenome.Celegans.UCSC."
    } else if (startsWith(genome_name, "dr")) {
      organism_prefix <- "BSgenome.Drerio.UCSC."
    } else {
      stop("Unknown genome name or unsupported organism. Please use the full name for the genome_name or check if the genome exists in BSgenome::available.genomes()")
    }

    # Construct the full package name
    genome_package <- paste0(organism_prefix, genome_name)
  }

  # Check if the package is installed
  if (!requireNamespace(genome_package, quietly = TRUE)) {
    message(paste0("The genome package '", genome_package, "' is not installed. Installing now..."))

    # Install the package
    tryCatch({
      BiocManager::install(genome_package)
    }, error = function(e) {
      stop(paste("Failed to install genome package:", genome_package, "\nError:", e$message))
    })
  }

  # Load the genome using getBSgenome()
  genome <- tryCatch({
    BSgenome::getBSgenome(genome_package)
  }, error = function(e) {
    stop(paste("Failed to load genome package:", genome_package, "\nError:", e$message))
  })

  return(genome)
}
