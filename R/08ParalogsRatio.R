#' Calculate ratio for paralogs
#'
#' @param paralogs A df of paralog annotations for all protein coding genes
#' @param input_genes A vector of input genes
#' @importFrom magrittr %>%
#' @return A dataframe with HGNC IDs and paralog score
#' @export
calculate_paralogs_ratio <- function(paralogs, input_genes) {
  # Input gene list ------------------------------------------------------------
  input_genes <- input_genes %>%
    dplyr::select(hgnc_id)

  # Remove anything below 30% paralogs - I've decided to take the mean -------
  paralogs_filtered <- paralogs %>%
    dplyr::filter(mean_paralog_perc >= 30) %>%
    dplyr::rename(hgnc_id = gene1_hgnc_id) %>%
    dplyr::select(hgnc_id, paralog_hgnc_id, mean_paralog_perc, max_paralog_perc)


  # Calculations ---------------------------------------------------------------
  # Counting the number of paralogs per gene
  paralogs_filtered_count1 <- paralogs_filtered %>%
    group_by(hgnc_id) %>%
    dplyr::mutate(num_of_paralogs = n())

  # Check if paralog is an input gene
  paralogs_filtered_count2 <- paralogs_filtered_count1 %>%
    dplyr::mutate(is_paralog_input_gene_yes_or_no = ifelse(paralog_hgnc_id %in% input_genes$hgnc_id, 1, 0))

  # Count how many paralogs are input genes
  paralogs_filtered_count3 <- paralogs_filtered_count2 %>%
    group_by(hgnc_id) %>%
    dplyr::mutate(num_input_gene_paralogs = sum(is_paralog_input_gene_yes_or_no))

  # Ratio of number of paralogs that are input gene : number of paralogs
  paralogs_filtered_ratio <- paralogs_filtered_count3 %>%
    group_by(hgnc_id) %>%
    dplyr::mutate(ratio_paraloginputgenes_to_paralogs = num_input_gene_paralogs / num_of_paralogs)

  paralogs_filtered_final <- paralogs_filtered_ratio %>%
    dplyr::select(!paralog_hgnc_id) %>%
    dplyr::select(!mean_paralog_perc) %>%
    dplyr::select(!max_paralog_perc) %>%
    dplyr::distinct(hgnc_id, .keep_all = TRUE) %>%
    dplyr::select(!is_paralog_input_gene_yes_or_no) %>%
    dplyr::mutate(ratio_paraloginputgenes_to_paralogs = ifelse(is.na(ratio_paraloginputgenes_to_paralogs),
      0, ratio_paraloginputgenes_to_paralogs
    ))

  cat("\n(8/12) finished running paralogs_ratio.R\n")
  return(paralogs_filtered_final)
}
