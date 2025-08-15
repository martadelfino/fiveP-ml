#' Calculate ratio for Protein Families
#'
#' @param panther A df of protein family annotations of input genes
#' @param input_genes A vector of input genes
#' @importFrom magrittr %>%
#' @return A dataframe with HGNC IDs and protein family score
#' @export
calculate_protein_families_ratio <- function(panther, input_genes) {
  # Input genes ----------------------------------------------------------------
  input_genes <- input_genes %>%
    dplyr::select(hgnc_id)

  # Panther protein families of input genes ------------------------------------
  panther_counts <- panther %>%
    mutate(input_gene_yes_or_no = ifelse(hgnc_id %in% input_genes$hgnc_id, 1, 0))

  # Calculations ---------------------------------------------------------------
  # Counting number of input proteins/genes in each family
  panther_counts <- panther_counts %>%
    group_by(family_id) %>%
    mutate(
      num_genes_in_family = n(),
      num_input_gene_per_family = sum(input_gene_yes_or_no)
    )

  # Counting the number of unique genes in each pathway that gene is related to
  panther_counts_per_gene <- panther_counts %>%
    group_by(hgnc_id) %>%
    dplyr::mutate(
      num_families = n_distinct(family_id),
      num_unique_genes_in_families = sum(length(unique(panther_counts$hgnc_id[panther_counts$family_id %in% family_id])) - 1),
      num_input_genes_in_families = sum(unique(panther_counts$hgnc_id[panther_counts$family_id %in% family_id]) %in% input_genes$hgnc_id) - (hgnc_id %in% input_genes$hgnc_id)
    ) %>%
    dplyr::mutate(ratio_input_genes_in_families = num_input_genes_in_families / num_unique_genes_in_families) %>%
    dplyr::select(
      hgnc_id, uniprot_ids, family_id, num_families,
      num_unique_genes_in_families, num_input_genes_in_families,
      ratio_input_genes_in_families
    ) %>%
    arrange(hgnc_id)

  panther_counts_per_gene_final <- panther_counts_per_gene %>%
    dplyr::select(
      hgnc_id, uniprot_ids, num_families,
      num_unique_genes_in_families, num_input_genes_in_families,
      ratio_input_genes_in_families
    ) %>%
    unique() %>%
    dplyr::mutate(ratio_input_genes_in_families = ifelse(is.na(ratio_input_genes_in_families), 0, ratio_input_genes_in_families))

  cat("\n(12/12) finished running protein_families_ratio.R\n")
  return(panther_counts_per_gene_final)
}
