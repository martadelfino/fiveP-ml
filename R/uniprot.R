#' Fetch Uniprot Data
#'
#' This function fetches protein data from Uniprot for all protein coding genes_______?.
#'
#' @param protein_coding_genes The df of all protein coding genes.
#' @param input_genes The df of input genes.
#' @return A dataframe with Uniprot data of the input genes.
#' @export
fetch_uniprot <- function(protein_coding_genes, input_genes, save_raw = FALSE, save_path = NULL) {
  batch_size = 1

  # Reading the protein coding genes file --------------------------------------
  hgnc_uniprot_symbol <- protein_coding_genes %>%
    dplyr::select(hgnc_id, uniprot_ids, symbol)


  # Checks of the protein coding genes file ------------------------------------

  # Removing all extra proteins identifiers (only keeping the canonical ones)
  hgnc_uniprot_symbol$uniprot_ids <- trimws(sub('\\|.*', "", hgnc_uniprot_symbol$uniprot_ids))


  # Reading the input gene list --------------------------------------------------

  input_genes <- input_genes %>%
    dplyr::select(hgnc_id)

  # Get gene symbols for the input genes
  input_genes_symbols <- input_genes %>% left_join(hgnc_uniprot_symbol,
                                                   by = 'hgnc_id')

  #print(input_genes_symbols)
  # Access uniprot ---------------------------------------------------------------

  # Note: Uniprot has no specific update schedule. So the entry version will be
  # important to keep track of database/entry versions.

  # Obtain gene symbols
  symbol <- dplyr::select(input_genes_symbols, `symbol`)
  #print(symbol)
  vector_symbol <- symbol %>% dplyr::pull(`symbol`)  # turn object into vector



  # Ensure input is a character vector
  if (!is.character(vector_symbol)) {
    stop("Input must be a character vector of gene names.")
  }

  # Split genes into batches
  batches <- split(vector_symbol, ceiling(seq_along(vector_symbol) / batch_size))

  # Initialize an empty list to store results
  results <- list()

  for (i in seq_along(batches)) {
    # Join genes into a query string with OR logic for the current batch
    gene_query <- paste(paste0("gene:", batches[[i]]), collapse = " OR ")

    # URL-encode the query string to handle special characters
    encoded_query <- URLencode(gene_query)

    # Construct the curl command with the specified genes and desired fields
    curl_command <- paste0(
      "curl -s -H \"Accept: text/plain; format=tsv\" \"https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+(",
      encoded_query,
      ")+AND+organism_id:9606&fields=accession,id,xref_hgnc,gene_primary,xref_complexportal,xref_panther,version\""
    )

    # Execute the curl command and capture the output
    output <- system(curl_command, intern = TRUE, ignore.stderr = TRUE)

    # Combine the output into a single string
    tsv_content <- paste(output, collapse = "\n")

    # Check if output contains valid content
    if (nchar(tsv_content) == 0) {
      warning(paste("No data returned for batch", i))
      next
    }

    # Convert the TSV content into a data frame by reading from a string
    batch_data <- read_tsv(I(tsv_content), col_types = cols(.default = "c"), show_col_types = FALSE)

    # Append the batch data to the results list
    results[[i]] <- batch_data
  }

  # Combine all batch results into a single data frame
  uniprot_input_gene_symbol_results <- do.call(rbind, results)

  # Save raw data
  if (save_raw) {
    if (is.null(save_path)) {
      save_path <- "data/uniprot_input_gene_symbol_results.csv"
    }
    readr::write_csv(uniprot_input_gene_symbol_results, save_path)
  }

  ## Clean results
  # Remove trailing ';' from the HGNC column
  uniprot_input_gene_symbol_results$HGNC <- gsub(";$", "",
                                                 uniprot_input_gene_symbol_results$HGNC)
  # Expand any rows with multiple HGNCs
  rows_to_separate <- which(sapply(strsplit(uniprot_input_gene_symbol_results$HGNC, ";"),
                                   function(x) sum(grepl("HGNC:", x)) > 1))
  # Check if rows_to_separate is empty
  if (length(rows_to_separate) > 0) {
    uniprot_input_gene_symbol_results_separated <- uniprot_input_gene_symbol_results %>%
      dplyr::slice(rows_to_separate) %>%
      tidyr::separate_rows(HGNC, sep = ";") %>%
      dplyr::mutate(HGNC = ifelse(grepl("HGNC:", HGNC), HGNC, paste0("HGNC:", HGNC)))
    #cat('print(uniprot_input_gene_symbol_results_separated)')
    #print(uniprot_input_gene_symbol_results_separated)

    # Combine the separated rows with the rest of the dataframe
    uniprot_input_gene_symbol_results_combined <- bind_rows(uniprot_input_gene_symbol_results[-rows_to_separate, ],
                                                            uniprot_input_gene_symbol_results_separated)

  } else { # continue
    uniprot_input_gene_symbol_results_combined <- uniprot_input_gene_symbol_results
  }

  # Checking duplicates
  sum(duplicated(uniprot_input_gene_symbol_results_combined$HGNC))
  duplicated_rows <- uniprot_input_gene_symbol_results_combined[duplicated(uniprot_input_gene_symbol_results_combined), ]
  # Remove duplicate rows
  uniprot_input_gene_symbol_results_combined <- distinct(uniprot_input_gene_symbol_results_combined)

  # checking duplicated HGNCs
  duplicated_rows_HGNC <- uniprot_input_gene_symbol_results_combined[duplicated(uniprot_input_gene_symbol_results_combined$HGNC), ]
  # some HGNCs just have more than one protein

  # Renaming columns
  uniprot_input_gene_symbol_results_cleaned <- uniprot_input_gene_symbol_results_combined %>%
    dplyr::rename(hgnc_id = HGNC) %>% dplyr::rename(uniprot_ids = Entry) %>%
    dplyr::rename(symbol = 'Gene Names (primary)')

  merged_df <- merge(uniprot_input_gene_symbol_results_cleaned,
                     input_genes_symbols, by = c("hgnc_id", "uniprot_ids"))

  # removing extra columns after the merge
  uniprot_input_gene_symbol_results_cleaned <- merged_df %>%
    dplyr::select(hgnc_id, uniprot_ids, symbol.x, ComplexPortal, PANTHER, `Entry version`) %>%
    dplyr::rename(symbol = symbol.x)

  cat('\n(5/12) finished running uniprot.R\n')
  return(uniprot_input_gene_symbol_results_cleaned)

}
