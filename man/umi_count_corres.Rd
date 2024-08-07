% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline.R
\name{umi_count_corres}
\alias{umi_count_corres}
\title{umi_count_corres}
\usage{
umi_count_corres(
  data,
  qual,
  dir,
  gene_bed,
  gtf = NULL,
  bar = "barcode",
  gene = "gene",
  isoform = "isoform",
  polyA = "polyA",
  sim_thresh = NULL,
  split = "|",
  sep = ",",
  splice_site_thresh = 10,
  verbose = FALSE,
  bed_gene_col = "gene",
  bed_strand_col = "strand",
  to_isoform = TRUE,
  filter_only_intron = TRUE,
  mid_offset_thresh = 3,
  overlap_thresh = 0,
  gtf_gene_col = "gene",
  gtf_start_col = "start",
  gtf_end_col = "end",
  gtf_iso_col = "transname",
  cores = 1
)
}
\arguments{
\item{data}{The first output dataframe from reads_extract_bc.}

\item{qual}{The adapter distance from barcode matching to represent the data quality}

\item{gtf}{The gtf annotation, each row is an exon for an isoform}

\item{bar}{The name of the column to represent cell barcode.}

\item{gene}{The name of column which stores the gene names}

\item{isoform}{The name of the column to represent the isoform.}

\item{polyA}{The name of the column to represent the polyA status.}

\item{sim_thresh}{the minimum threshold of needleman score to connect to different UMIs.}

\item{split}{The character to seperate different exons in the isoform representation.}

\item{sep}{The character to seperate the start and end position of an exon in the isoform representation.}

\item{splice_site_thresh}{The minimum threshold of occurance  for an splice site to be preserved}

\item{verbose}{if print information during UMI clustering for each cell.}

\item{filter_only_intron}{A boolean to indicate if read only cover intron part should be preserved.}

\item{overlap_thresh}{The minimum threshold for the coverage of the annotated isoform.}

\item{gtf_gene_col}{the name of the column which stores gene name in the gtf.}

\item{gtf_start_col}{the name of the column which stores the start position of exon in the gtf.}

\item{gtf_end_col}{the name of the column which stores the end position of exon in the gtf.}

\item{gtf_iso_col}{the name of the column which stores isoform name in the gtf.}

\item{cores}{The number of cores required for parallel}
}
\description{
A wrapper function to include all steps in the UMI deduplication.
}
\details{
A wrapper function to include all steps in the UMI deduplication, including UMI cluster,
isoform correction, isoform imputation and quantification.
}
