# Refuse reuse after input contents, settings, or package versions change.
check_run_manifest <- function(work_dir,paths,settings,force=FALSE) {
  paths <- normalizePath(paths,mustWork=TRUE)
  inputs <- data.frame(path=paths,bytes=unname(file.size(paths)),md5=unname(tools::md5sum(paths)))
  current <- list(inputs=inputs,settings=settings)
  manifest <- file.path(work_dir,"run_manifest.rds")
  cached <- any(file.exists(file.path(work_dir,c("annotation/gene_bed.rds",
                   "BarcodeMatch/BarcodeMatchIso.txt","out/iso_count.txt"))))
  if(!force && cached && (!file.exists(manifest) || !identical(readRDS(manifest),current)))
    stop("Cached outputs have different or unrecorded inputs/settings. Use a new work_dir or force_rerun=TRUE.")
  saveRDS(current,manifest)
  write.table(inputs,file.path(work_dir,"input_manifest.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  writeLines(c(paste("LongcellPre",settings$LongcellPre),paste("Native backend",settings$native_backend),
    "Export: reference-derived exon sequence; synthetic quality characters; no original-read variants or terminal poly(A) tails.",
    "BAM cell metadata is in read names and reads_annot.csv; no CB/UB tags are synthesized.",
    capture.output(sessionInfo())),file.path(work_dir,"provenance.txt"))
  invisible(current)
}
