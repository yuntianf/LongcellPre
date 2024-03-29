% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reads_filter.R
\name{isoform_size_filter}
\alias{isoform_size_filter}
\title{isoform_size_filter}
\usage{
isoform_size_filter(isoforms, size, ratio = 0.1)
}
\arguments{
\item{isoforms}{A string vector, each string represent a corrected isoform by collapse of reads
from the same UMI}
}
\value{
A numeric vector to indicating whether the cluster is filtered out or not
}
\description{
Filter out singletons and small clusters which are not likely to be
correctly clustered, given the distribution of cluster size within a series of isoform clusters.
}
