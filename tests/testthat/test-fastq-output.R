# Substitute only reference lookup, using real Biostrings and ShortRead objects.
# Each isoform string in this fixture is its desired sequence length.
export_fixture <- function() {
  env <- new.env(parent = environment(isoformCount2ReadsForGene))
  for (name in c("isoformCount2ReadsForGene", "isoformCount2Reads",
                 "extractAnnotFromQname")) {
    fun <- get(name, envir = env)
    environment(fun) <- env
    env[[name]] <- fun
  }
  env$extractSeq <- function(genome, chr, isoforms, strand, ...) {
    Biostrings::DNAStringSet(strrep("A", as.integer(isoforms)))
  }
  env$read2bins <- function(read,...) matrix(c(1,as.numeric(read)),ncol=2)
  env
}

count_fixture <- function() {
  data.frame(cell = c("c1", "c2", "c1"), gene = c("g2", "g1", "g2"),
             isoform = c("50001", "4", "7"), count = c(1, 2, 1),
             polyA = c(0, 1, 0.5))
}

bed_fixture <- function() {
  data.frame(gene = c("g1", "g2"), chr = "chr1", strand = c("+", "-"))
}

test_that("FASTQ qualities cover reads beyond 50 kb", {
  env <- export_fixture()
  mat <- count_fixture()[rep(1, 3), ]
  mat$isoform <- c("49999", "50000", "50001")
  reads <- env$isoformCount2ReadsForGene(mat, NULL, "chr1", "+")
  expect_equal(Biostrings::width(ShortRead::sread(reads)), c(49999L, 50000L, 50001L))
  expect_equal(nchar(as.character(methods::slot(methods::slot(reads, "quality"), "quality"))), c(49999L, 50000L, 50001L))
  expect_true(methods::validObject(reads))
})

test_that("zero and fractional counts preserve rounding without invalid empty reads", {
  env <- export_fixture()
  mat <- count_fixture()[1, ]
  mat$count <- 0.5
  expect_length(env$isoformCount2ReadsForGene(mat, NULL, "chr1", "+"), 0)
  mat <- mat[rep(1, 4), ]
  mat$count <- c(0, 0.49, 0.51, 2.5)
  expect_length(env$isoformCount2ReadsForGene(mat, NULL, "chr1", "+"), 3)
  expect_length(env$isoformCount2ReadsForGene(mat[FALSE, ], NULL, "chr1", "+"), 0)
})

test_that("invalid counts and FASTQ quality characters fail clearly", {
  env <- export_fixture()
  mat <- count_fixture()[1, ]
  for (bad in c(NA_real_, Inf, -1)) {
    mat$count <- bad
    expect_error(env$isoformCount2ReadsForGene(mat, NULL, "chr1", "+"), "Counts must")
  }
  mat$count <- 1
  for (bad in c("", "!!", " ", NA_character_)) {
    expect_error(env$isoformCount2ReadsForGene(mat, NULL, "chr1", "+", quality = bad), "quality must")
  }
})

test_that("exon extraction checks relative bounds and concatenates without attached packages", {
  extractReadSeq <- extractReadSeq
  environment(extractReadSeq) <- baseenv()
  seq <- Biostrings::DNAString("ACGTACGT")
  expect_identical(as.character(extractReadSeq(seq, 100, 107,
    matrix(c(100, 102, 105, 107), ncol = 2, byrow = TRUE))), "ACGCGT")
  expect_identical(as.character(extractReadSeq(seq, 100, 107,
    matrix(c(100, 100), ncol = 2))), "A")
  expect_error(extractReadSeq(seq, 100, 107, matrix(c(99, 101), ncol = 2)), "coordinates")
  expect_error(extractReadSeq(seq, 100, 107, matrix(c(106, 108), ncol = 2)), "coordinates")
  expect_error(extractReadSeq(seq, 100, 107, matrix(c(105, 101), ncol = 2)), "coordinates")
})

test_that("sequence reconstruction concatenates exons and orients negative-strand genes", {
  testthat::local_mocked_bindings(
    getSeq = function(x, names, start, end, strand) {
      Biostrings::subseq(x, start = start, end = end)
    }, .package = "BSgenome")
  env <- new.env(parent = environment(extractSeq))
  env$read2bins <- function(read, sep = ",", split = "|") {
    parts <- strsplit(read, split, fixed = TRUE)[[1]]
    do.call(rbind, lapply(parts, function(x) as.numeric(strsplit(x, sep, fixed = TRUE)[[1]])))
  }
  extract <- extractSeq
  environment(extract) <- env
  genome <- Biostrings::DNAString("ACGTACGT")
  expect_identical(as.character(extract(genome, "chr1", c("1,2|5,7", "2,4"), "+")),
                   c("ACACG", "CGT"))
  expect_identical(as.character(extract(genome, "chr1", "1,2|5,7", "-")), "CGTGT")
})

test_that("streamed FASTQ and annotations agree with the legacy return mode", {
  env <- export_fixture()
  mat <- count_fixture()
  genome <- list(chr1 = NULL)
  fastq <- tempfile(fileext = ".fq.gz")
  legacy <- tempfile(fileext = ".fq.gz")
  annot <- tempfile(fileext = ".csv")
  on.exit(unlink(c(fastq, legacy, annot)))
  returned <- env$isoformCount2Reads(mat, genome, bed_fixture(), legacy)
  result <- env$isoformCount2Reads(mat, genome, bed_fixture(), fastq,
                                  return_reads = FALSE, annotation_file = annot)
  observed <- ShortRead::readFastq(fastq)
  expect_equal(result, 4)
  expect_identical(as.character(ShortRead::sread(observed)), as.character(ShortRead::sread(returned)))
  expect_identical(as.character(ShortRead::id(observed)), as.character(ShortRead::id(returned)))
  expect_identical(as.character(methods::slot(methods::slot(observed, "quality"), "quality")), as.character(methods::slot(methods::slot(returned, "quality"), "quality")))
  table <- read.csv(annot)
  expect_identical(table$qname, as.character(ShortRead::id(observed)))
  expect_identical(table$annot, c("c1", "c1", "c2", "c2"))
  expect_false(anyDuplicated(table$qname) > 0L)
})

test_that("failed exports leave previous outputs intact and remove temporary files", {
  env <- export_fixture()
  folder <- tempfile()
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  fastq <- file.path(folder, "out.fq.gz")
  annot <- file.path(folder, "reads.csv")
  writeLines("previous FASTQ", fastq)
  writeLines("previous CSV", annot)
  bed <- bed_fixture()
  bed$strand[1] <- "invalid"
  expect_error(env$isoformCount2Reads(count_fixture(), list(chr1 = NULL), bed, fastq,
                                    return_reads = FALSE, annotation_file = annot), "exactly one")
  expect_identical(readLines(fastq), "previous FASTQ")
  expect_identical(readLines(annot), "previous CSV")
  expect_length(list.files(folder, all.files = TRUE, pattern = "^\\.longcell-"), 0)
})

test_that("missing chromosomes and no exportable counts are reported", {
  env <- export_fixture()
  fastq <- tempfile(fileext = ".fq.gz")
  on.exit(unlink(fastq))
  mat <- count_fixture()
  mat$count <- 0
  expect_error(env$isoformCount2Reads(mat, list(chr1 = NULL), bed_fixture(), fastq), "No reads remain")
  mat <- count_fixture()
  bed <- bed_fixture()
  bed$chr[bed$gene == "g2"] <- "absent"
  expect_error(env$isoformCount2Reads(mat, list(chr1 = NULL), bed, fastq,
                                    return_reads = FALSE), "absent")
})

test_that("annotation extraction respects the requested field", {
  expect_identical(extractAnnotFromQname("cell=c1|gene=g1|polyA=0_1", "gene"), "g1")
})
