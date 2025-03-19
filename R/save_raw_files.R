#' saving raw files for each P
#'
#' @param input_genes The df of input genes.
#' @param save_path path to save the raw files.
#' @return just a print statement.the files will be saved.
#' @export
save_raw_files <- function(input_genes, save_path) {

  save_raw = TRUE

  # Data fetching functions ----------------------------------------------------
  hgnc_gene_list <- fetch_hgnc_gene_list(save_raw = save_raw, save_path = save_path)
  paralogues <- fetch_paralogues(hgnc_gene_list, save_raw = save_raw, save_path = save_path)
  pathways <- fetch_pathways(hgnc_gene_list, input_genes, save_raw = save_raw, save_path = save_path)
  ppi <- fetch_ppi(hgnc_gene_list, save_raw = save_raw, save_path = save_path)
  uniprot <- fetch_uniprot(hgnc_gene_list, input_genes, save_raw = save_raw, save_path = save_path)
  protein_complex <- fetch_protein_complex(hgnc_gene_list, uniprot, save_raw = save_raw, save_path = save_path)
  protein_families <- fetch_protein_families(hgnc_gene_list, uniprot, save_raw = save_raw, save_path = save_path)

  cat("saved raw files")
  invisible(NULL)
}
