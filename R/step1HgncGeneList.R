#' Fetch all protein coding genes data
#'
#' This function fetches data from EBI to get all protein coding genes.
#'
#' @param save_raw Boolean for choosing whether to save the raw data or not.
#' @param save_path String for the path to save the raw data. If NULL, it will save to "data/hgnc_gene_list.csv".
#' @return A dataframe with protein coding gene data from EBI database.
#' @export
fetch_hgnc_gene_list <- function(save_raw = FALSE, save_path = NULL) {
  # Fetch hgnc gene file
  protein_coding_genes <- readr::read_delim("https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/locus_types/gene_with_protein_product.txt",
    delim = "\t",
    col_names = TRUE,
    col_select = c(hgnc_id, uniprot_ids, symbol, ensembl_gene_id, entrez_id)
  ) %>%
    as.data.frame()

  # Save raw data
  if (save_raw) {
    if (is.null(save_path)) {
      save_path <- "data/hgnc_gene_list.csv"
    }
    readr::write_csv(protein_coding_genes, save_path)
  }

  cat("\n(1/12) finished running hgnc_gene_list.R\n")
  return(protein_coding_genes)
}
