#' Fetch Paralogs Data
#'
#' This function fetches paralogs data from Ensembl for all protein coding genes.
#'
#' @param protein_coding_genes The df of all protein coding genes.
#' @return A dataframe with paralogs data from Ensembl.
#' @export
fetch_paralogs <- function(protein_coding_genes, chunk_size = 100, save_raw = FALSE, save_path = NULL) {

  # Select necessary columns
  hgnc_ensembl <- protein_coding_genes %>%
    dplyr::select(hgnc_id, ensembl_gene_id)

  # Split the ensembl_gene_id vector into chunks
  ensembl_chunks <- split(
    hgnc_ensembl$ensembl_gene_id,
    ceiling(seq_along(hgnc_ensembl$ensembl_gene_id) / chunk_size)
  )

  # Connect to Ensembl Mart
  human <- biomaRt::useEnsembl(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl"
  )

  # Function to query a single chunk
  # (Ensembl uses the word paralog instead of paralog)
  query_chunk <- function(chunk) {
    biomaRt::getBM(
      attributes = c("ensembl_gene_id",
                     "hsapiens_paralog_ensembl_gene",
                     "hsapiens_paralog_orthology_type",
                     "hsapiens_paralog_perc_id",
                     "hsapiens_paralog_perc_id_r1",
                     "version"),
      filters = "ensembl_gene_id",
      values = chunk,
      mart = human
    )
  }

  # Query each chunk and handle errors
  paralogs_list <- lapply(ensembl_chunks, function(chunk) {
    tryCatch({
      query_chunk(chunk)
    }, error = function(e) {
      message("Error querying chunk: ", e$message)
      NULL
    })
  })

  # Combine results into a single dataframe
  paralogs <- do.call(rbind, paralogs_list)

  # Save raw data
  if (save_raw) {
    if (is.null(save_path)) {
      save_path <- "data/paralogs.csv"
    }
    readr::write_csv(paralogs, save_path)
  }

  # Clean the results as before
  paralogs_cleaned <- paralogs %>%
    dplyr::left_join(hgnc_ensembl, by = 'ensembl_gene_id') %>%
    dplyr::select(hgnc_id, ensembl_gene_id, hsapiens_paralog_orthology_type,
                  hsapiens_paralog_ensembl_gene, hsapiens_paralog_perc_id,
                  hsapiens_paralog_perc_id_r1) %>%
    dplyr::rename(gene1_hgnc_id = hgnc_id) %>%
    dplyr::rename(gene1_ensembl_gene_id = ensembl_gene_id)

  paralogs_cleaned_with_paralog_hgnc <- paralogs_cleaned %>%
    left_join(hgnc_ensembl, join_by(hsapiens_paralog_ensembl_gene == ensembl_gene_id),
              relationship = "many-to-many") %>%
    dplyr::rename(paralog_hgnc_id = hgnc_id) %>%
    dplyr::rename(paralog_ensembl_gene_id = hsapiens_paralog_ensembl_gene) %>%
    dplyr::rename(paralog_perc_id = hsapiens_paralog_perc_id) %>%
    dplyr::rename(paralog_perc_id_r1 = hsapiens_paralog_perc_id_r1)

  paralogs_cleaned_reorg <- paralogs_cleaned_with_paralog_hgnc %>%
    dplyr::select(gene1_hgnc_id, gene1_ensembl_gene_id, paralog_hgnc_id,
                  paralog_ensembl_gene_id, hsapiens_paralog_orthology_type,
                  paralog_perc_id, paralog_perc_id_r1) %>%
    dplyr::arrange(gene1_hgnc_id)

  paralogs_cleaned_reorg <- paralogs_cleaned_reorg %>%
    dplyr::mutate(mean_paralog_perc = (paralog_perc_id + paralog_perc_id_r1) / 2) %>%
    dplyr::mutate(max_paralog_perc = pmax(paralog_perc_id, paralog_perc_id_r1))

  cat('\n(2/12) finished running paralogs.R\n')
  return(paralogs_cleaned_reorg)
}
