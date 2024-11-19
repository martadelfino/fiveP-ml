#' Calculate ratio for PPI
#'
#' @param ppi A df of ppi annotations for all protein coding genes
#' @param input_genes A df of input genes
#' @return A dataframe with HGNC IDs and PPI scores
#' @export
calculate_ppi_ratio <- function(ppi, input_genes) {

  # PPI annotations ------------------------------------------------------------

  ppi <- ppi %>%
    dplyr::select(!protein1_string_id) %>%
    dplyr::select(!protein2_string_id) %>%
    dplyr::rename(hgnc_id = protein1_hgnc_id)


  # Input genes ----------------------------------------------------------------

  input_genes <- input_genes %>%
    dplyr::select(hgnc_id)


  # Calculations ---------------------------------------------------------------

  # Counting the number of interactions per gene
  ppi_count1 <- ppi %>%
    group_by(hgnc_id) %>%
    dplyr::mutate(num_of_interactions = n())

  # Check if interaction protein is an input gene
  ppi_count2 <- ppi_count1 %>%
    mutate(is_interaction_input_gene_yes_or_no = ifelse(protein2_hgnc_id %in% input_genes$hgnc_id, 1, 0))

  # Count how many interactions are input genes
  ppi_count3 <- ppi_count2 %>%
    group_by(hgnc_id) %>%
    mutate(num_input_gene_interactions = sum(is_interaction_input_gene_yes_or_no))

  # Ratio of number of input gene interactions : number of interactions
  ppi_ratio <- ppi_count3 %>%
    group_by(hgnc_id) %>%
    dplyr::mutate(ratio_interactioninputgenes_to_interactions = num_input_gene_interactions / num_of_interactions)

  ppi_ratio_final <- ppi_ratio %>%
    dplyr::select(!protein2_hgnc_id) %>%
    dplyr::distinct(hgnc_id, .keep_all = TRUE) %>%
    dplyr::select(!is_interaction_input_gene_yes_or_no) %>%
    dplyr::mutate(ratio_interactioninputgenes_to_interactions = ifelse(is.na(ratio_interactioninputgenes_to_interactions),
                                                                       0, ratio_interactioninputgenes_to_interactions))


  cat('\n(10/12) finished running ppi_ratio.R\n')
  return(ppi_ratio_final)
}

