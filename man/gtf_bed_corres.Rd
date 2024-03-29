% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exon_corres.R
\name{gtf_bed_corres}
\alias{gtf_bed_corres}
\title{gtf_bed_corres}
\usage{
gtf_bed_corres(
  gtf,
  gene_bed,
  genes = NULL,
  gtf_gene_col = "gene",
  gtf_start_col = "start",
  gtf_end_col = "end",
  transname = "transname",
  bed_gene_col = "gene",
  bed_start_col = "start",
  bed_end_col = "end",
  id = "id",
  ...
)
}
\arguments{
\item{gtf}{A dataframe of gtf annotation, each row is an exon in an isoform}

\item{gene_bed}{A dataframe of gene_bed annotation, each row is an exon part}

\item{gtf_gene_col}{The name of the column storing the gene in gtf}

\item{gtf_start_col}{The name of the column storing the exon start position in gtf}

\item{gtf_end_col}{The name of the column storing the exon end position in gtf}

\item{transname}{The name of the column storing the transcript id in gtf}

\item{bed_gene_col}{The name of the column storing the gene in gene bed}

\item{bed_start_col}{The name of the column storing the exon start position in gene bed}

\item{bed_end_col}{The name of the column storing the exon end position in gene bed}

\item{id}{The name of the column storing the exon part id in gene bed}
}
\value{
A data frame including 2 columns, one is the transcript name from gtf, the other is its exon part form.
}
\description{
mapping canonical isoforms from gtf annotation to exon parts in gene bed
}
