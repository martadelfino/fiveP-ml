# fiveP-ml

fiveP is a package created to identify protein coding genes similar to a user-defined set of input or seed genes, based on four protein annotations and paralog status. The five 'Ps' are: protein complex, protein family, pathways, protein-protein interactions (PPI), and paralog.

The package identifies similar genes by computing similarity ratios for each individual P.

The similarities are calculated with the following equations, where the higher the score (between 0 and 1), the more similar that gene is to the input genes:

#### Protein complexes, Protein families, Pathways

Using pathways as the example P: For each protein coding gene, the number of unique input genes in input gene pathways, divided by the number of unique genes in input gene pathways.

$r(g)=\frac{|S\cap P(S)|}{P(S)}$

#### Paralogs, PPI

Using PPI as the example: For each protein coding gene, the number of input genes that gene interacts with, divided by the number of genes it interacts with:

$r(g)=\frac{|P(g)\cap S|}{|P(g)|}$

#### Usage

The package can be used with one function, get_fiveP(), or with individual functions for each step of the workflow. Examples are below.

The output is always for all protein coding genes. Depending on how the package is used, the user can get ratios for all fivePs or for whatever Ps they select.

The output of the fiveP can be kept as similarity ratios, where the higher the ratio the more similar that gene is to the input genes. Or, it can be turned into a binary score, which simply states whether that gene is in the same P/interaction/paralog with any of the input genes.

For the genes without data, a 'NA' is left in its place.

### Quick start / Installation

#### Dependencies

-   devtools

-   dplyr (2.5.0)

-   purr (1.0.2)

-   tidyr (1.3.1)

-   magrittr (2.0.3)

-   readr

-   STRINGdb (2.18.0)

-   biomaRt (2.62.0)

#### Install package

`devtools::install_github("martadelfino/fiveP-ml")`

#### Load package

`library(fivePml)`

#### Examples

1.  The get_fiveP() function

The get_fiveP() function is the easiest way to get a df with all five Ps. (Please note: if you get a connection error at any point, restart the session and try again. This may be due to one of the packages or file locations.)

```{r}

input_genes <- tibble(hgnc_id = c("HGNC:19743", 
                                  "HGNC:9202",
                                  "HGNC:8653",
                                  "HGNC:6936",
                                  "HGNC:4878"))

result <- get_fiveP(input_genes, binary = TRUE) # this will give you scores of 0 or 1.
head(result)


```

2.  Running the workflow steps manually

However, it is also possible to get one P at a time. The workflow is separated into three main steps:

1.  Data fetching
2.  Data processing / ratio calculations
3.  Merging

(Please note: if you get a connection error at any point, try that function again. You may need to restart the session.)

```{r}
# 1. Data fetching functions ----------------------------------------------------
# First, all protein coding genes are obtained
hgnc_gene_list <- fetch_hgnc_gene_list()
# Then the five P are obtained. Some functions require the input genes at this step. 
paralogs <- fetch_paralogs(hgnc_gene_list)
pathways <- fetch_pathways(hgnc_gene_list, input_genes)
ppi <- fetch_ppi(hgnc_gene_list)
# (Accessing UNIPROT is required for the protein complex and protein families step)
uniprot <- fetch_uniprot(hgnc_gene_list, input_genes)
protein_complex <- fetch_protein_complex(hgnc_gene_list, uniprot)
protein_families <- fetch_protein_families(hgnc_gene_list, uniprot)

# 2. Data processing functions --------------------------------------------------
# Second, the functions for the ratio calculations. 
# Some functions require the input genes at this step.
paralogs_ratio <- calculate_paralogs_ratio(paralogs, input_genes)
pathways_ratio <- calculate_pathways_ratio(
  pathways$input_genes_Uniprot2Reactome,
  pathways$Uniprot2Reactome_final_hgnc_no_na,
  input_genes
)
ppi_ratio <- calculate_ppi_ratio(ppi, input_genes)
protein_complex_ratio <- calculate_protein_complex_ratio(protein_complex, input_genes)
protein_families_ratio <- calculate_protein_families_ratio(protein_families, input_genes)

# 3. Merging everything ---------------------------------------------------------
# Then all fiveP are meregd. Any of the fiveP can be selected.  
protein_coding_genes <- hgnc_gene_list %>%
  dplyr::select(hgnc_id)
input_genes <- input_genes %>%
  dplyr::select(hgnc_id)
protein_complexes <- protein_complex_ratio %>%
  dplyr::select(hgnc_id, ratio_input_genes_in_complexes)
protein_families <- protein_families_ratio %>%
  dplyr::select(hgnc_id, ratio_input_genes_in_families)
pathways <- pathways_ratio %>%
  dplyr::select(hgnc_id, ratio_input_genes_in_pathways)
paralogs <- paralogs_ratio %>%
  dplyr::select(hgnc_id, ratio_paraloginputgenes_to_paralogs)
ppi <- ppi_ratio %>%
  dplyr::select(hgnc_id, ratio_interactioninputgenes_to_interactions)

list_of_dfs <- list(
  protein_coding_genes, protein_complexes, protein_families,
  pathways, paralogs, ppi
)

results <- list_of_dfs %>%
  purrr::reduce(left_join, by = "hgnc_id") %>%
  dplyr::rename(protein_complex_score = ratio_input_genes_in_complexes) %>%
  dplyr::rename(protein_family_score = ratio_input_genes_in_families) %>%
  dplyr::rename(pathway_score = ratio_input_genes_in_pathways) %>%
  dplyr::rename(paralog_score = ratio_paraloginputgenes_to_paralogs) %>%
  dplyr::rename(ppi_score = ratio_interactioninputgenes_to_interactions) %>%
  arrange(
    desc(protein_complex_score), desc(protein_family_score),
    desc(pathway_score), desc(paralog_score), desc(ppi_score)
  )

# Optional: turning the ratios into binary scores -------------------------------
results_binary <- results %>%
  dplyr::mutate(across(ends_with("_score"), ~ ifelse(. > 0, 1, 0)))

```

### Other comments

License is GPL (\>= 3) because for PPI, we use STRINGdb, which has GLP 2 license.
