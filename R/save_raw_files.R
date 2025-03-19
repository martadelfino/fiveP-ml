#' saving raw files for each P
#'
#' @param input_genes The df of input genes.
#' @param save_path path to save the raw files.
#' @return just a print statement.the files will be saved.
#' @export
save_raw_files <- function(input_genes, save_path) {

  save_raw = TRUE

  # Ensure the directory exists; if not, create it
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  # Define unique file names for each output
  hgnc_file          <- file.path(save_path, "hgnc_gene_list.csv")
  paralogues_file    <- file.path(save_path, "paralogues.csv")
  pathways_file      <- file.path(save_path, "pathways.csv")
  ppi_file           <- file.path(save_path, "ppi.csv")
  uniprot_file       <- file.path(save_path, "uniprot.csv")
  protein_complex_file  <- file.path(save_path, "protein_complex.csv")
  protein_families_file <- file.path(save_path, "protein_families.csv")

  # Data fetching functions with their respective file paths
  hgnc_gene_list    <- fetch_hgnc_gene_list(save_raw = save_raw, save_path = hgnc_file)
  paralogues        <- fetch_paralogues(hgnc_gene_list, save_raw = save_raw, save_path = paralogues_file)
  pathways          <- fetch_pathways(hgnc_gene_list, input_genes, save_raw = save_raw, save_path = pathways_file)
  ppi               <- fetch_ppi(hgnc_gene_list, save_raw = save_raw, save_path = ppi_file)
  uniprot           <- fetch_uniprot(hgnc_gene_list, input_genes, save_raw = save_raw, save_path = uniprot_file)
  protein_complex   <- fetch_protein_complex(hgnc_gene_list, uniprot, save_raw = save_raw, save_path = protein_complex_file)
  protein_families  <- fetch_protein_families(hgnc_gene_list, uniprot, save_raw = save_raw, save_path = protein_families_file)

  cat("saved raw files")
  invisible(NULL)
}
