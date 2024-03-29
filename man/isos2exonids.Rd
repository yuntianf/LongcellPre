% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exon_corres.R
\name{isos2exonids}
\alias{isos2exonids}
\title{isos2exonids}
\usage{
isos2exonids(
  data,
  gene_bed,
  gene_col = "gene",
  start = "start",
  end = "end",
  id = "id",
  mid_bias = 0,
  end_bias = 10,
  end_overlap = 10,
  nonsense_label = "N",
  split = "|",
  sep = ","
)
}
\arguments{
\item{data}{A dataframe storing reads information, the gene and isoform column are required}

\item{gene_col}{The name of the column which stores the gene name in the gene bed}

\item{start}{the name of the column which stores the start position for each exon in the gene bed.}

\item{end}{the name of the column which stores the end position for each exon in the gene bed.}

\item{id}{the name of the column which stores the exon id in the gene bed.}

\item{mid_bias}{The tolerance for the bias for mapping of exons in the middle of the isoform}

\item{end_bias}{The tolerance for the bias for mapping of exons at the end of the isoform}

\item{end_overlap}{The minimum length of the end exon}

\item{nonsense_label}{The label for exons which couldn't be mapped to gene bed}

\item{split}{The character to split the exons in the isoform}

\item{sep}{The character to split the start and end position for each exon in the isoform.}
}
\value{
A dataframe as input data, but replace the isoforms by exon id sequences
}
\description{
map the exon to the annotated exons gene bed for all genes
}
