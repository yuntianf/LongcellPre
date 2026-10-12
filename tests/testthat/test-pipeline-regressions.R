test_that("gene distribution preserves complete genes, row order, and old assignment", {
  data <- data.frame(gene = rep(c("g1", "g2", "g3", "g4"), c(12, 8, 4, 2)),
                     row = seq_len(26))
  chunks <- genes_distribute(data, 2)
  expect_identical(chunks[[1]]$gene, c(rep("g1", 12), rep("g4", 2)))
  expect_identical(chunks[[2]]$gene, c(rep("g2", 8), rep("g3", 4)))
  expect_identical(sort(unlist(lapply(chunks, `[[`, "row"))), seq_len(26))
  expect_identical(genes_distribute(data, 1), list(data))
})

test_that("nested output directories can be initialized", {
  root <- tempfile()
  on.exit(unlink(root, recursive = TRUE))
  out <- file.path(root, "nested", "output")
  init(out)
  expect_true(all(dir.exists(file.path(out, c("annotation", "BarcodeMatch", "bam", "out")))))
})

test_that("regenerating collapsed FASTQ also regenerates a cached BAM", {
  env <- new.env(parent = environment(UMI_consensus_out))
  env$UMI_consensus_out <- UMI_consensus_out
  environment(env$UMI_consensus_out) <- env
  env$load_genome <- function(...) NULL
  env$isoformCount2Reads <- function(..., return_reads, annotation_file) {
    expect_false(return_reads)
    expect_match(annotation_file, "reads_annot.csv$")
  }
  mapped <- FALSE
  env$fastqMap <- function(...) { mapped <<- TRUE }
  out <- tempfile()
  dir.create(out)
  on.exit(unlink(out, recursive = TRUE))
  file.create(file.path(out, c("UMI_collapsed.bam", "UMI_collapsed.bam.bai")))
  env$UMI_consensus_out(NULL, out, "genome.fa", "hg38", NULL)
  expect_true(mapped)
})
