% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paralogs_ratio.R
\name{calculate_paralogs_ratio}
\alias{calculate_paralogs_ratio}
\title{Calculate ratio for paralogs}
\usage{
calculate_paralogs_ratio(paralogs, input_genes)
}
\arguments{
\item{paralogs}{A df of paralog annotations for all protein coding genes}

\item{input_genes}{A vector of input genes}
}
\value{
A dataframe with HGNC IDs and paralog score
}
\description{
Calculate ratio for paralogs
}
