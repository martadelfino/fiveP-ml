% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manual_check.R
\name{check_protein_complex}
\alias{check_protein_complex}
\title{Check for protein complexes}
\usage{
check_protein_complex(
  protein_coding_genes,
  uniprot_input_gene_symbol_results_cleaned
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
Check for protein complexes
}
