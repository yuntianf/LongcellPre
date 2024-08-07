% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splice_site_correct.R
\name{cluster_isoform_correct}
\alias{cluster_isoform_correct}
\title{cluster_isoform_correct}
\usage{
cluster_isoform_correct(start, mid, end, concensus, polyA, preserve_mid)
}
\arguments{
\item{start}{A numeric vector indicating start position of reads within the UMI cluster}

\item{mid}{A string vector indicating patterns of middle splicing sites}

\item{end}{A numeric vector indicating end position of reads within the UMI cluster}

\item{polyA}{A boolean vector indicating the polyA existence}

\item{preserve_mid}{A string vector indicating preserved patterns of middle splice sites}
}
\value{
A vector including the corrected start, middle splicing sites, end position, polyA existence
and the number of reads within this UMI cluster
}
\description{
Corret for wrong mapping and truncation within each UMI cluster
}
\details{
Corret for wrong mapping and truncation within each UMI cluster
}
