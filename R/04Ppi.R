#' Fetch PPI Data
#'
#' This function fetches PPI data from STRING for all protein coding genes.
#'
#' @param save_raw Boolean for choosing whether to save the raw data or not.
#' @param save_path String for the path to save the raw data.
#' @param save_raw Boolean for choosing whether to save the raw data or not.
#' @param save_path String for the path to save the raw data.
#' @param protein_coding_genes The df of all protein coding genes.
#' @importFrom magrittr %>%
#' @return A dataframe with data from Database.
#' @export

fetch_ppi <- function(protein_coding_genes, save_raw = FALSE, save_path = NULL) {
  # Reading the protein coding genes file --------------------------------------
  hgnc_esembl <- protein_coding_genes %>%
    dplyr::select(hgnc_id, ensembl_gene_id)
  hgnc_ensembl <- data.frame(hgnc_esembl)

  # Load STRING database -------------------------------------------------------
  string_db <- STRINGdb::STRINGdb$new(
    version = "12.0", species = 9606,
    score_threshold = 700, network_type = "full",
    input_directory = ""
  )

  # Save raw data
  if (save_raw) {
    if (is.null(save_path)) {
      save_path <- "data/string_db.rds"
    }
    saveRDS(string_db, file = "string_db.rds")
  }

  # Mapping input gene list to the STRING identifiers
  input_genes_mapped <- string_db$map(hgnc_ensembl, "ensembl_gene_id", removeUnmappedRows = TRUE)
  input_genes_mapped_vector <- input_genes_mapped %>% pull(STRING_id)

  # Get interactions
  interactions <- string_db$get_interactions(input_genes_mapped_vector)

  # remove duplicate ones
  interactions <- interactions %>%
    dplyr::distinct()

  # Clean the data -------------------------------------------------------------
  # Map back to HGNC IDs
  interactions_cleaned <- input_genes_mapped %>%
    left_join(interactions, join_by(STRING_id == from),
      relationship = "many-to-many"
    ) %>%
    dplyr::select(hgnc_id, STRING_id, to, combined_score) %>%
    dplyr::rename(protein1_hgnc_id = hgnc_id) %>%
    dplyr::rename(protein1_string_id = STRING_id) %>%
    dplyr::rename(protein2_string_id = to)

  interactions_cleaned_hgnc <- input_genes_mapped %>%
    left_join(interactions_cleaned, join_by(STRING_id == protein2_string_id),
      relationship = "many-to-many"
    ) %>%
    dplyr::select(protein1_hgnc_id, protein1_string_id, hgnc_id, STRING_id, combined_score) %>%
    dplyr::rename(protein2_hgnc_id = hgnc_id) %>%
    dplyr::rename(protein2_string_id = STRING_id) %>%
    dplyr::arrange(protein1_hgnc_id)

  cat("\n(4/12) finished running ppi.R\n")
  return(interactions_cleaned_hgnc)
}
