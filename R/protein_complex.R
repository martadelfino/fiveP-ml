#' Fetch Protein Complex Data from Complex Portal
#'
#' This function fetches data from Complex Portal Database for all pc coding genes.
#'
#' @param protein_coding_genes The df of all protein coding genes.
#' @param uniprot_input_gene_symbol_results_cleaned The df of uniprot results of input genes.
#' @return A dataframe with Protein Complex data of the input genes from Complex Portal Database.
#' @export
fetch_protein_complex <- function(protein_coding_genes,
                                  uniprot_input_gene_symbol_results_cleaned, save_raw = FALSE, save_path = NULL) {

  # Creating a new df of the complexes from uniprot results --------------------

  input_genes_protein_complexes_expanded <- uniprot_input_gene_symbol_results_cleaned %>%
    tidyr::separate_rows(ComplexPortal, sep = ";") %>% distinct() %>%
    filter(ComplexPortal != "") %>%
    dplyr::select(ComplexPortal) %>% distinct() %>%
    dplyr::rename(complex_id = ComplexPortal)

  # removing extra bits
  input_genes_protein_complexes_expanded$complex_id <- trimws(sub("\\[.*?\\]", "", input_genes_protein_complexes_expanded$complex_id))


  # Querying ComplexPortal -----------------------------------------------------

  ComplexPortal <- read_delim('https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv')

  # Save raw data
  if (save_raw) {
    if (is.null(save_path)) {
      save_path <- "data/ComplexPortal.csv"
    }
    readr::write_csv(ComplexPortal, save_path)
  }

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

  cat('\n(6/12) finished running protein_complex.R\n')
  return(input_genes_complexportal_participants_hgnc)

}
