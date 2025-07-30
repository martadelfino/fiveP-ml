#' Main function for obtaining fiveP scores for the input genes
#'
#' @param input_genes A vector of gene HGNC IDs
#' @param binary Logical, if TRUE the scores will be binary (0 or 1), if FALSE the scores will be the actual ratios
#' @return A dataframe with the fiveP scores
#' @export
get_fiveP <- function(input_genes, binary = TRUE) { # eventually I can add options to save the intermediate files too

  # Data fetching functions ----------------------------------------------------
  hgnc_gene_list <- fetch_hgnc_gene_list()
  paralogs <- fetch_paralogs(hgnc_gene_list)
  pathways <- fetch_pathways(hgnc_gene_list, input_genes)
  ppi <- fetch_ppi(hgnc_gene_list)
  uniprot <- fetch_uniprot(hgnc_gene_list, input_genes)
  protein_complex <- fetch_protein_complex(hgnc_gene_list, uniprot)
  protein_families <- fetch_protein_families(hgnc_gene_list, uniprot)

  # Data processing functions --------------------------------------------------
  paralogs_ratio <- calculate_paralogs_ratio(paralogs, input_genes)
  pathways_ratio <- calculate_pathways_ratio(pathways$input_genes_Uniprot2Reactome,
                                             pathways$Uniprot2Reactome_final_hgnc_no_na,
                                             input_genes)
  ppi_ratio <- calculate_ppi_ratio(ppi, input_genes)
  protein_complex_ratio <- calculate_protein_complex_ratio(protein_complex, input_genes)
  protein_families_ratio <- calculate_protein_families_ratio(protein_families, input_genes)

  # Merging everything ---------------------------------------------------------
  protein_coding_genes <- hgnc_gene_list %>%
    dplyr::select(hgnc_id)
  input_genes <- input_genes %>%
    dplyr::select(hgnc_id)
  protein_complexes <- protein_complex_ratio %>%
    dplyr::select(hgnc_id, ratio_input_genes_in_complexes)
  protein_families <- protein_families_ratio %>%
    dplyr::select(hgnc_id, ratio_input_genes_in_families)
  pathways <- pathways_ratio %>%
    dplyr::select(hgnc_id, ratio_input_genes_in_pathways)
  paralogs <- paralogs_ratio %>%
    dplyr::select(hgnc_id, ratio_paraloginputgenes_to_paralogs)
  ppi <- ppi_ratio %>%
    dplyr::select(hgnc_id, ratio_interactioninputgenes_to_interactions)

  list_of_dfs <- list(protein_coding_genes, protein_complexes, protein_families,
                      pathways, paralogs, ppi)

  results <- reduce(list_of_dfs, left_join, by = "hgnc_id") %>%
    dplyr::rename(protein_complex_score = ratio_input_genes_in_complexes) %>%
    dplyr::rename(protein_family_score = ratio_input_genes_in_families) %>%
    dplyr::rename(pathway_score = ratio_input_genes_in_pathways) %>%
    dplyr::rename(paralog_score = ratio_paraloginputgenes_to_paralogs) %>%
    dplyr::rename(ppi_score = ratio_interactioninputgenes_to_interactions) %>%
    arrange(desc(protein_complex_score), desc(protein_family_score),
            desc(pathway_score), desc(paralog_score), desc(ppi_score))

  if (binary) {
    results <- results %>%
      dplyr::mutate(across(ends_with("_score"), ~ ifelse(. > 0, 1, 0)))
  }

  cat("finished")
  return(results)
}
