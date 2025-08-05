#' Fetch Pathway Data
#'
#' This function fetches pathway data from Reactome for all protein coding genes and input genes.
#'
#' @param save_raw Boolean for choosing whether to save the raw data or not.
#' @param save_path String for the path to save the raw data.
#' @param protein_coding_genes The df of all protein coding genes.
#' @param input_genes The df of input genes.
#' @importFrom magrittr %>%
#' @return A list of two dataframes with pathway data from Reactome 1. input genes pathway data, 2. All protein coding genes pathway data
#' @export
fetch_pathways <- function(protein_coding_genes, input_genes, save_raw = FALSE, save_path = NULL) {
  # Reading the protein coding genes file --------------------------------------

  hgnc_uniprot_symbol_entrez <- protein_coding_genes %>%
    dplyr::select(hgnc_id, uniprot_ids, symbol, entrez_id)

  # Removing all extra proteins identifiers (only keeping the canonical ones)
  hgnc_uniprot_symbol_entrez$uniprot_ids <- trimws(sub("\\|.*", "", hgnc_uniprot_symbol_entrez$uniprot_ids))


  # Reading the input gene list ------------------------------------------------

  input_genes <- input_genes %>%
    dplyr::select(hgnc_id)


  # Get gene symbols for the input genes ---------------------------------------

  input_genes_entrez_id <- input_genes %>%
    inner_join(hgnc_uniprot_symbol_entrez, by = "hgnc_id") %>%
    pull(entrez_id) %>%
    as.character(.)


  # Retrieving the genes in all pathways ---------------------------------------

  # Obtaining data from Reactome directly, lowest level pathways
  Uniprot2Reactome <- read_delim("https://reactome.org/download/current/UniProt2Reactome.txt",
    col_names = FALSE
  )

  # Save raw data
  if (save_raw) {
    if (is.null(save_path)) {
      save_path <- "data/Uniprot2Reactome.csv"
    }
    write.csv(Uniprot2Reactome, save_path, row.names = FALSE)
  }

  # Cleaning the Uniprot to Reactome file
  Uniprot2Reactome_cleaned <- Uniprot2Reactome %>%
    dplyr::rename(
      uniprot_ids = X1, pathway_id = X2, url = X3, pathway_name = X4,
      evidence_code = X5, species = X6
    )
  Uniprot2Reactome_cleaned <- Uniprot2Reactome_cleaned %>%
    dplyr::filter(species == "Homo sapiens")
  Uniprot2Reactome_final <- Uniprot2Reactome_cleaned %>%
    dplyr::select(uniprot_ids, pathway_id)

  # Joining the file with the protein coding genes file
  Uniprot2Reactome_final_hgnc <- Uniprot2Reactome_final %>%
    left_join(hgnc_uniprot_symbol_entrez, by = "uniprot_ids", relationship = "many-to-many")
  # Some proteins don't map to hgnc ids. this is because they are immunoglobulins
  # or other immune related proteins. These are ignored.

  Uniprot2Reactome_final_hgnc_no_na <- Uniprot2Reactome_final_hgnc %>%
    filter(!is.na(hgnc_id))

  # Filter the Uniprot2Reactome final file with the input genes
  input_genes_Uniprot2Reactome <- Uniprot2Reactome_final_hgnc_no_na %>%
    dplyr::filter(hgnc_id %in% input_genes$hgnc_id)


  cat("\n(3/12) finished running pathways.R\n")
  return(list(
    input_genes_Uniprot2Reactome = input_genes_Uniprot2Reactome,
    Uniprot2Reactome_final_hgnc_no_na = Uniprot2Reactome_final_hgnc_no_na
  ))
}
