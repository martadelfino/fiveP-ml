#' Calculate ratio for Paralogues
#'
#' @param paralogues A df of paralogue annotations for all protein coding genes
#' @param input_genes A vector of input genes
#' @return A dataframe with HGNC IDs and paralogue score
#' @export
calculate_paralogues_ratio <- function(paralogues, input_genes) {

  # Input gene list ------------------------------------------------------------

  input_genes <- input_genes %>%
    dplyr::select(hgnc_id)


  # Remove anything below 30% paralogues - I've decided to take the mean -------

  paralogues_filtered <- paralogues %>%
    filter(mean_paralog_perc >= 30) %>%
    rename(hgnc_id = gene1_hgnc_id) %>%
    dplyr::select(hgnc_id, paralog_hgnc_id, mean_paralog_perc, max_paralog_perc)


  # Calculations ---------------------------------------------------------------

  # Counting the number of paralogues per gene
  paralogues_filtered_count1 <- paralogues_filtered %>%
    group_by(hgnc_id) %>%
    dplyr::mutate(num_of_paralogs = n())

  # Check if paralog is an input gene
  paralogues_filtered_count2 <- paralogues_filtered_count1 %>%
    mutate(is_paralog_input_gene_yes_or_no = ifelse(paralog_hgnc_id %in% input_genes$hgnc_id, 1, 0))

  # Count how many paralogues are input genes
  paralogues_filtered_count3 <- paralogues_filtered_count2 %>%
    group_by(hgnc_id) %>%
    mutate(num_input_gene_paralogs = sum(is_paralog_input_gene_yes_or_no))

  # Ratio of number of paralogs that are input gene : number of paralogs
  paralogues_filtered_ratio <- paralogues_filtered_count3 %>%
    group_by(hgnc_id) %>%
    dplyr::mutate(ratio_paraloginputgenes_to_paralogs = num_input_gene_paralogs / num_of_paralogs)

  paralogues_filtered_final <- paralogues_filtered_ratio %>%
    dplyr::select(!paralog_hgnc_id) %>%  dplyr::select(!mean_paralog_perc) %>%
    dplyr::select(!max_paralog_perc) %>%
    dplyr::distinct(hgnc_id, .keep_all = TRUE) %>%
    dplyr::select(!is_paralog_input_gene_yes_or_no) %>%
    dplyr::mutate(ratio_paraloginputgenes_to_paralogs = ifelse(is.na(ratio_paraloginputgenes_to_paralogs),
                                                               0, ratio_paraloginputgenes_to_paralogs))


  cat('\n(8/12) finished running paralogues_ratio.R\n')
  return(paralogues_filtered_final)
}
