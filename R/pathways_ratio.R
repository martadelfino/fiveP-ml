#' Calculate ratio for Pathways
#'
#' @param input_genes_pathways A df of pathway annotations based on input genes
#' @param Uniprot2Reactome A df of pathway annotations of all protein coding genes
#' @param input_genes A df of input genes
#' @return A dataframe with HGNC IDs and pathway score
#' @export
calculate_pathways_ratio <- function(input_genes_pathways, Uniprot2Reactome, input_genes) {

  # Input genes  ---------------------------------------------------------------

  input_genes <- input_genes %>%
    dplyr::select(hgnc_id)

  # Filtering Uniprot2Reactome by selecting the pathways from input genes pathways -----

  input_genes_Uniprot2Reactome <- Uniprot2Reactome %>%
    dplyr::filter(pathway_id %in% input_genes_pathways$pathway_id)


  # Counting the number of input genes per pathway -------------------------------

  # Checking if the individual genes are input genes or not
  reactome_counts <- input_genes_Uniprot2Reactome %>%
    mutate(input_gene_yes_or_no = ifelse(hgnc_id %in% input_genes$hgnc_id, 1, 0))

  # Counting the number of input genes per pathway
  reactome_counts <- reactome_counts %>%
    group_by(pathway_id) %>%
    mutate(num_genes_in_pathway = n(),
           numb_input_gene_per_pathway = sum(input_gene_yes_or_no))


  # Creating a df for the number of unique genes in each pathway -----------------

  # Counting the number of unique genes in each pathway that gene is related to
  reactome_counts_per_gene <- reactome_counts %>%
    group_by(hgnc_id) %>%
    dplyr::mutate(
      num_pathways = n_distinct(pathway_id),
      num_unique_genes_in_pathways = sum(length(unique(reactome_counts$hgnc_id[reactome_counts$pathway_id %in% pathway_id])) - 1),
      num_input_genes_in_pathways = sum(unique(reactome_counts$hgnc_id[reactome_counts$pathway_id %in% pathway_id]) %in% input_genes$hgnc_id) - (hgnc_id %in% input_genes$hgnc_id)
    ) %>%
    dplyr::mutate(ratio_input_genes_in_pathways = num_input_genes_in_pathways / num_unique_genes_in_pathways) %>%
    dplyr::select(hgnc_id, uniprot_ids, pathway_id, num_pathways,
                  num_unique_genes_in_pathways, num_input_genes_in_pathways,
                  ratio_input_genes_in_pathways) %>% arrange(hgnc_id)

  reactome_counts_per_gene_final <- reactome_counts_per_gene %>%
    dplyr::select(hgnc_id, uniprot_ids, num_pathways,
                  num_unique_genes_in_pathways, num_input_genes_in_pathways,
                  ratio_input_genes_in_pathways) %>% unique() %>%
    dplyr::mutate(ratio_input_genes_in_pathways = ifelse(is.na(ratio_input_genes_in_pathways),
                                                         0, ratio_input_genes_in_pathways))


  cat('\n(9/12) finished running pathways_ratio.R\n')
  return(reactome_counts_per_gene_final)
}
