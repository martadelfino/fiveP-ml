% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/protein_complex.R
\name{fetch_protein_complex}
\alias{fetch_protein_complex}
\title{Fetch Protein Complex Data from Complex Portal}
\usage{
fetch_protein_complex(
  protein_coding_genes,
  uniprot_input_gene_symbol_results_cleaned,
  save_raw = FALSE,
  save_path = NULL
)
}
\arguments{
\item{protein_coding_genes}{The df of all protein coding genes.}

\item{uniprot_input_gene_symbol_results_cleaned}{The df of uniprot results of input genes.}
}
\value{
A dataframe with Protein Complex data of the input genes from Complex Portal Database.
}
\description{
This function fetches data from Complex Portal Database for all pc coding genes.
}
