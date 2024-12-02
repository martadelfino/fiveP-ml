#' Fetch Paralogues Data
#'
#' This function fetches paralogues data from Ensembl for all protein coding genes.
#'
#' @param protein_coding_genes The df of all protein coding genes.
#' @return A dataframe with paralogues data from Ensembl.
#' @export
fetch_paralogues <- function(protein_coding_genes) {

  # Select what I need from the protein coding genes file
  hgnc_ensembl <- protein_coding_genes %>%
    dplyr::select(hgnc_id, ensembl_gene_id)

  ensembl_id_vector <- hgnc_ensembl %>% pull(ensembl_gene_id)


  # Querying BioBart for gene paralogues -----------------------------------------

  human <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
  paralogues <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                              "hsapiens_paralog_ensembl_gene",
                                              "hsapiens_paralog_orthology_type",
                                              "hsapiens_paralog_perc_id",
                                              "hsapiens_paralog_perc_id_r1",
                                              "version"),
                               filters = "ensembl_gene_id",
                               values = ensembl_id_vector,
                               mart = human)

  # Cleaning the paralogues file -------------------------------------------------

  # Map back the input gene ensembl IDs back to HGNC IDs
  paralogues_cleaned <- paralogues %>%
    dplyr::left_join(hgnc_ensembl, by = 'ensembl_gene_id') %>%
    dplyr::select(hgnc_id, ensembl_gene_id, hsapiens_paralog_orthology_type,
                  hsapiens_paralog_ensembl_gene, hsapiens_paralog_perc_id,
                  hsapiens_paralog_perc_id_r1) %>%
    dplyr::rename(gene1_hgnc_id = hgnc_id) %>%
    dplyr::rename(gene1_ensembl_gene_id = ensembl_gene_id)

  # Now map the paralog ensembl IDs to HGNC IDs
  paralogues_cleaned_with_paralog_hgnc <- paralogues_cleaned %>%
    left_join(hgnc_ensembl, join_by(hsapiens_paralog_ensembl_gene == ensembl_gene_id),
              relationship = "many-to-many") %>%
    dplyr::rename(paralog_hgnc_id = hgnc_id) %>%
    dplyr::rename(paralog_ensembl_gene_id = hsapiens_paralog_ensembl_gene) %>%
    dplyr::rename(paralog_perc_id = hsapiens_paralog_perc_id) %>%
    dplyr::rename(paralog_perc_id_r1 = hsapiens_paralog_perc_id_r1)

  # Cleaning it up
  paralogues_cleaned_reorg <- paralogues_cleaned_with_paralog_hgnc %>%
    dplyr::select(gene1_hgnc_id, gene1_ensembl_gene_id, paralog_hgnc_id,
                  paralog_ensembl_gene_id, hsapiens_paralog_orthology_type,
                  paralog_perc_id, paralog_perc_id_r1) %>%
    dplyr::arrange(gene1_hgnc_id)

  # Taking the mean of the two paralog percentages
  paralogues_cleaned_reorg <- paralogues_cleaned_reorg %>%
    dplyr::mutate(mean_paralog_perc = (paralog_perc_id + paralog_perc_id_r1) / 2) %>%
    dplyr::mutate(max_paralog_perc = pmax(paralog_perc_id, paralog_perc_id_r1))


  cat('\n(2/12) finished running paralogues.R\n')
  return(paralogues_cleaned_reorg)
}
