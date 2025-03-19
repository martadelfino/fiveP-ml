#' Fetch Protein Family Data from Uniprot
#'
#' This function fetches protein family data from Uniprot for the input genes.
#'
#' @param protein_coding_genes The df of all protein coding genes.
#' @param uniprot_input_gene_symbol_results_cleaned The df of uniprot results of input genes.
#' @return A dataframe with Protein Family data from Uniprot for the input genes.
#' @export
fetch_protein_families <- function(protein_coding_genes,
                                   uniprot_input_gene_symbol_results_cleaned, save_raw = FALSE, save_path = NULL) {

  # Creating a df of protein families data from uniprot results ----------------

  input_genes_protein_families_expanded <- uniprot_input_gene_symbol_results_cleaned %>%
    tidyr::separate_rows(PANTHER, sep = ";") %>% distinct() %>%
    filter(PANTHER != "") %>%
    dplyr::select(PANTHER) %>% distinct() %>%
    dplyr::rename(family_id = PANTHER)

  # removing extra information after the ':'
  input_genes_protein_families_expanded$family_id <- trimws(sub(":.*", "", input_genes_protein_families_expanded$family_id))


  # Querying Uniprot -----------------------------------------------------------

  batch_size = 1

  # Obtain families
  family <- dplyr::select(input_genes_protein_families_expanded, family_id)
  vector_family <- family %>% dplyr::pull(family_id)   # turning object into vector

  # Ensure input is a character vector
  if (!is.character(vector_family)) {
    stop("Input must be a character vector of gene names.")
  }

  # Split genes into batches
  batches <- split(vector_family, ceiling(seq_along(vector_family) / batch_size))

  # Initialize an empty list to store results
  results <- list()

  for (i in seq_along(batches)) {
    # Join genes into a query string with OR logic for the current batch
    family_query <- paste(paste0("xref:", batches[[i]]), collapse = " OR ")

    # URL-encode the query string to handle special characters
    encoded_query <- URLencode(family_query)

    # Construct the curl command with the specified genes and desired fields
    curl_command <- paste0(
      "curl -s -H \"Accept: text/plain; format=tsv\" \"https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+(",
      encoded_query,
      ")+AND+organism_id:9606&fields=accession,xref_hgnc,gene_primary,xref_panther,version\""
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
  uniprot_input_gene_family_results <- do.call(rbind, results)

  # Save raw data
  if (save_raw) {
    if (is.null(save_path)) {
      save_path <- "data/uniprot_input_gene_family_results.csv"
    }
    readr::write_csv(uniprot_input_gene_family_results, save_path)
  }

  # Cleaning Protein Families result file from Uniprot ---------------------------

  # Selecting and renaming required columns
  proteinfamily_genes <- uniprot_input_gene_family_results %>%
    dplyr::select(Entry, HGNC, 'Gene Names (primary)', PANTHER) %>%
    dplyr::rename(uniprot_ids = Entry) %>%
    dplyr::rename(hgnc_id = HGNC) %>%
    dplyr::rename(family_id = PANTHER) %>%
    dplyr::rename(symbol = 'Gene Names (primary)')

  # Removing trailing ;
  proteinfamily_genes$hgnc_id <- gsub(";$", "", proteinfamily_genes$hgnc_id)
  proteinfamily_genes$family_id <- gsub(";$", "", proteinfamily_genes$family_id)

  # Separating families into new rows
  proteinfamily_genes_expanded <- proteinfamily_genes %>%
    tidyr::separate_rows(family_id, sep = ";") %>%
    filter(family_id != "")

  # removing extra bits
  proteinfamily_genes_expanded$family_id <- trimws(sub(":.*", "", proteinfamily_genes_expanded$family_id))
  proteinfamily_genes_expanded <- proteinfamily_genes_expanded %>% distinct() %>%
    dplyr::select(family_id, uniprot_ids, hgnc_id, symbol) %>% # fixing order of columns
    arrange(family_id) # rearranging rows


  cat('\n(7/12) finished running protein_families.R\n')
  return(proteinfamily_genes_expanded)

}

