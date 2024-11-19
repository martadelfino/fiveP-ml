#' Calculate ratio for Protein Complex
#'
#' @param complexportal A df of protein complex annotations of input genes
#' @param input_genes A vector of input genes
#' @return A dataframe with HGNC IDs and protein complex score
#' @export
calculate_protein_complex_ratio <- function(complexportal, input_genes) {

  # Input genes ----------------------------------------------------------------

  input_genes <- input_genes %>%
    dplyr::select(hgnc_id)


  # ComplexPortal protein complexes of input genes -----------------------------

  complexportal_counts <- complexportal %>%
    mutate(input_gene_yes_or_no = ifelse(hgnc_id %in% input_genes$hgnc_id, 1, 0))


  # Calculations ---------------------------------------------------------------

  # Counting the number of input genes per complex
  complexportal_counts <- complexportal_counts %>%
    group_by(complex_id) %>%
    mutate(num_genes_in_complex = n(),
           num_input_gene_per_complex = sum(input_gene_yes_or_no))

  # Counting the number of unique genes in each complex that gene is related to
  complexportal_counts_per_gene <- complexportal_counts %>%
    group_by(hgnc_id) %>%
    dplyr::mutate(
      num_complexes = n_distinct(complex_id),
      num_unique_genes_in_complexes = sum(length(unique(complexportal_counts$hgnc_id[complexportal_counts$complex_id %in% complex_id])) - 1),
      num_input_genes_in_complexes = sum(unique(complexportal_counts$hgnc_id[complexportal_counts$complex_id %in% complex_id]) %in% input_genes$hgnc_id) - (hgnc_id %in% input_genes$hgnc_id)
    ) %>%
    dplyr::mutate(ratio_input_genes_in_complexes = num_input_genes_in_complexes / num_unique_genes_in_complexes) %>%
    dplyr::select(hgnc_id, symbol, uniprot_ids, complex_id, num_complexes,
                  num_unique_genes_in_complexes, num_input_genes_in_complexes,
                  ratio_input_genes_in_complexes)

  complexportal_counts_per_gene_final <- complexportal_counts_per_gene %>%
    dplyr::select(hgnc_id, symbol, uniprot_ids, num_complexes,
                  num_unique_genes_in_complexes, num_input_genes_in_complexes,
                  ratio_input_genes_in_complexes) %>% unique() %>%
    dplyr::mutate(ratio_input_genes_in_complexes = ifelse(is.na(ratio_input_genes_in_complexes),
                                                          0, ratio_input_genes_in_complexes))


  cat('\n(11/12) finished running protein_complex_ratio.R\n')
  return(complexportal_counts_per_gene_final)
}
