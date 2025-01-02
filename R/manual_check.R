#' Check for protein complexes
#'
#'
#'
#' @param protein_coding_genes The df of all protein coding genes.
#' @param uniprot_input_gene_symbol_results_cleaned The df of uniprot results of input genes.
#' @return A dataframe with Protein Complex data of the input genes from Complex Portal Database.
#' @export
check_protein_complex <- function(protein_coding_genes,
                                  uniprot_input_gene_symbol_results_cleaned) {

  # Creating a new df of the complexes from uniprot results --------------------

  input_genes_protein_complexes_expanded <- uniprot_input_gene_symbol_results_cleaned %>%
    tidyr::separate_rows(ComplexPortal, sep = ";") %>%
    distinct() %>%
    filter(ComplexPortal != "") %>%
    dplyr::select(ComplexPortal) %>% distinct() %>%
    dplyr::rename(complex_id = ComplexPortal)

  # removing extra bits
  input_genes_protein_complexes_expanded$complex_id <- trimws(sub("\\[.*?\\]", "", input_genes_protein_complexes_expanded$complex_id))


  # Querying ComplexPortal -----------------------------------------------------

  ComplexPortal <- read_delim('https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv')


  ## Cleaning ComplexPortal file
  # Note: the ComplexPortal is updated every month, so the date will be important

  # Selecting and renaming required columns
  ComplexPortal_participants <- ComplexPortal %>%
    dplyr::select(`#Complex ac`, `Expanded participant list`) %>%
    dplyr::rename(complex_id = `#Complex ac`) %>%
    dplyr::rename(uniprot_ids = `Expanded participant list`)

  # Creating a new row for each protein
  ComplexPortal_participants_separated <- ComplexPortal_participants %>%
    separate_rows(uniprot_ids, sep = '\\|')

  # Removing extra information enclosed in '[]' and after '-'
  ComplexPortal_participants_separated$uniprot_ids <- gsub("\\[.*?\\]", "", ComplexPortal_participants_separated$uniprot_ids)
  ComplexPortal_participants_separated$uniprot_ids <- sub("-.*", "", ComplexPortal_participants_separated$uniprot_ids)

  # Removing the stochiometry information in the brackets
  ComplexPortal_participants_separated$uniprot_ids <- trimws(sub("\\(.*?\\)", "", ComplexPortal_participants_separated$uniprot_ids))

  # Removing rows without a uniprot_id (these would've been other things like molecules)
  ComplexPortal_participants_separated <- subset(ComplexPortal_participants_separated, uniprot_ids != "")


  # Only keeping the complexes identified in the input gene lists --------------

  # Joining the input genes protein complexes with the participants for each complex
  input_genes_complexportal_participants <- input_genes_protein_complexes_expanded %>%
    left_join(ComplexPortal_participants_separated, by = 'complex_id' ,
              relationship = "many-to-many")

  # Adding hgnc_id

  hgnc_uniprot_symbol <- protein_coding_genes %>%
    dplyr::select(hgnc_id, uniprot_ids, symbol)

  hgnc_uniprot_symbol$uniprot_ids <- trimws(sub('\\|.*', "", hgnc_uniprot_symbol$uniprot_ids))

  input_genes_complexportal_participants_hgnc <- input_genes_complexportal_participants %>%
    left_join(hgnc_uniprot_symbol, by = 'uniprot_ids',
              relationship = "many-to-many") %>% unique()


  #rows_with_na <- input_genes_complexportal_participants_hgnc[!complete.cases(input_genes_complexportal_participants_hgnc), ]
  #print(rows_with_na). could I use the checker Pilar made?

  cat('\n(1/5) finished running check_protein_complex()\n')
  return(input_genes_complexportal_participants_hgnc)

}



#' Check protein families
#'
#'
#'
#' @param protein_coding_genes The df of all protein coding genes.
#' @param uniprot_input_gene_symbol_results_cleaned The df of uniprot results of input genes.
#' @return A dataframe with Protein Family data from Uniprot for the input genes.
#' @export
check_protein_families <- function(protein_coding_genes,
                                   uniprot_input_gene_symbol_results_cleaned) {

  # Creating a df of protein families data from uniprot results ----------------

  input_genes_protein_families_expanded <- uniprot_input_gene_symbol_results_cleaned %>%
    tidyr::separate_rows(PANTHER, sep = ";") %>% distinct() %>%
    filter(PANTHER != "") %>%
    dplyr::select(PANTHER) %>% distinct() %>%
    dplyr::rename(family_id = PANTHER)

  # removing extra information after the ':'
  input_genes_protein_families_expanded$family_id <- trimws(sub("\\:.*?", "", input_genes_protein_families_expanded$family_id))


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
    output <- system(curl_command, intern = TRUE)

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
  proteinfamily_genes_expanded$family_id <- trimws(sub("\\:.*", "", proteinfamily_genes_expanded$family_id))
  proteinfamily_genes_expanded <- proteinfamily_genes_expanded %>% distinct() %>%
    dplyr::select(family_id, uniprot_ids, hgnc_id, symbol) %>% # fixing order of columns
    arrange(family_id) # rearranging rows


  cat('\n(2/5) finished running check_protein_families()\n')
  return(proteinfamily_genes_expanded)

}


#' Check Pathways
#'
#'
#'
#' @param protein_coding_genes The df of all protein coding genes.
#' @param input_genes The df of input genes.
#' @return A df of input genes pathway data
#' @export
check_pathways <- function(protein_coding_genes, input_genes) {

  # Reading the protein coding genes file --------------------------------------

  hgnc_uniprot_symbol_entrez <- protein_coding_genes %>%
    dplyr::select(hgnc_id, uniprot_ids, symbol, entrez_id)

  # Removing all extra proteins identifiers (only keeping the canonical ones)
  hgnc_uniprot_symbol_entrez$uniprot_ids <- trimws(sub('\\|.*', "", hgnc_uniprot_symbol_entrez$uniprot_ids))


  # Reading the input gene list ------------------------------------------------

  input_genes <- input_genes %>%
    dplyr::select(hgnc_id)


  # Get gene symbols for the input genes ---------------------------------------

  input_genes_entrez_id <- input_genes %>%
    inner_join(hgnc_uniprot_symbol_entrez, by = 'hgnc_id') %>%
    pull(entrez_id) %>%
    as.character(.)


  # Retrieving the genes in all pathways ---------------------------------------

  # Obtaining data from Reactome directly, lowest level pathways
  Uniprot2Reactome <- read_delim('https://reactome.org/download/current/UniProt2Reactome.txt',
                                 col_names = FALSE)

  # Cleaning the Uniprot to Reactome file
  Uniprot2Reactome_cleaned <- Uniprot2Reactome %>%
    dplyr::rename(uniprot_ids = X1, pathway_id = X2, url = X3, pathway_name = X4,
                  evidence_code = X5, species = X6)
  Uniprot2Reactome_cleaned <- Uniprot2Reactome_cleaned %>%
    dplyr::filter(species == 'Homo sapiens')
  Uniprot2Reactome_final <- Uniprot2Reactome_cleaned %>%
    dplyr::select(uniprot_ids, pathway_id)

  # Joining the file with the protein coding genes file
  Uniprot2Reactome_final_hgnc <- Uniprot2Reactome_final %>%
    left_join(hgnc_uniprot_symbol_entrez, by = 'uniprot_ids', relationship = 'many-to-many')
  # Some proteins don't map to hgnc ids. this is because they are immunoglobulins
  # or other immune related proteins. These are ignored.

  Uniprot2Reactome_final_hgnc_no_na <- Uniprot2Reactome_final_hgnc %>%
    filter(!is.na(hgnc_id))

  # Filter the Uniprot2Reactome final file with the input genes
  input_genes_Uniprot2Reactome <- Uniprot2Reactome_final_hgnc_no_na %>%
    dplyr::filter(hgnc_id %in% input_genes$hgnc_id)


  cat('\n(3/5) finished running check_pathways()\n')
  return(input_genes_Uniprot2Reactome = input_genes_Uniprot2Reactome)

}



#' Check for Paralogues
#'
#' @param paralogues A df of paralogue annotations for all protein coding genes
#' @param input_genes A vector of input genes
#' @return A dataframe with HGNC IDs and paralogues for manual check
#' @export
check_paralogues <- function(paralogues, input_genes) {

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
    dplyr::mutate(is_paralog_input_gene_yes_or_no = ifelse(paralog_hgnc_id %in% input_genes$hgnc_id, 1, 0))

  # Count how many paralogues are input genes
  paralogues_filtered_count3 <- paralogues_filtered_count2 %>%
    group_by(hgnc_id) %>%
    dplyr::mutate(num_input_gene_paralogs = sum(is_paralog_input_gene_yes_or_no)) %>%
    dplyr::filter(is_paralog_input_gene_yes_or_no == 1)

  cat('\n(4/5) finished running check_paralogues()\n')
  return(paralogues_filtered_count3)
}



#' Check for PPI
#'
#' @param ppi A df of ppi annotations for all protein coding genes
#' @param input_genes A df of input genes
#' @return A dataframe with HGNC IDs and PPI scores for input genes
#' @export
check_ppi <- function(ppi, input_genes) {

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

  ppi_final <- ppi_count3 %>%
    dplyr::select(!protein2_hgnc_id) %>%
    dplyr::distinct(hgnc_id, .keep_all = TRUE)

  cat('\n(5/5) finished running check_ppi()\n')
  return(ppi_final)

}



#' Main function for checks
#'
#' @param protein_coding_genes A df of protein coding genes
#' @param uniprot_input_gene_symbol_results_cleaned description
#' @param paralogues description
#' @param ppi description
#' @param input_genes A df of input genes
#' @return A dataframe with HGNC IDs and protein features for input genes
#' @export
gene_check <- function(protein_coding_genes,
                       uniprot_input_gene_symbol_results_cleaned,
                       paralogues, ppi, input_genes) {

  protein_complex <- check_protein_complex(protein_coding_genes,
                                           uniprot_input_gene_symbol_results_cleaned)
  protein_families <- check_protein_families(protein_coding_genes,
                                             uniprot_input_gene_symbol_results_cleaned)
  pathways <- check_pathways(protein_coding_genes, input_genes)
  paralogues <- check_paralogues(paralogues, input_genes)
  ppi <- check_ppi(ppi, input_genes)

  dfs <- list(protein_complex, protein_families, pathways, paralogues, ppi)

  # Merge all data frames by 'id' using Reduce
  merged_df <- Reduce(function(x, y) merge(x, y, by = "hgnc_id", all = TRUE), dfs)

  merged_df <- merged_df %>% distinct()

  cat('\n finished running all gene checks.\n')
  return(merged_df)
}





