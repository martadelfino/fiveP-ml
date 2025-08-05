# fiveP-ml

This repo contains the scripts to obtain the '5 Ps' of input genes: protein complex, protein family, paralogs, pathways, and protein-protein interactions (PPI). The output is a data frame with all protein coding genes and a binary score (0 or 1) describing whether each gene is in the same 'P' as any of the input genes.

### Calculations for each P

Paralogue and PPI are scores are calculated differently than protein complex, protein family, and pathway scores.

#### Protein complexes, Protein families, Pathways

Using pathways as the example P: For all protein coding genes, a score is calculated by looking only at pathways of the input genes. The formula for the pathway score for each protein coding gene:

( \# of unique input genes in input gene pathways ) / ( \# of unique genes in input gene pathways )

#### Paralogs, PPI

Using PPI as the example: For all protein coding genes, a score is calculated by dividing the number of input genes by all genes in the PPI. The formula for the PPI score for each protein coding gene:

( \# of input genes that gene interacts with ) / ( \# of genes it interacts with )

#### Final results

The result is each protein coding gene is assigned a score from 0 to 1, for each 'P'. There is also the option to keep these as ratios. For the genes without data, a 'NA' is left in its place. The output is a data frame with HGNC IDs, and columns for each P.

### Quick start / Installation

Install devtools

`install.packages("devtools")`

Install dependencies

```{r}
#if (!requireNamespace("tidyverse", quietly = TRUE))
#install.packages("tidyverse")

#if (!requireNamespace("devtools", quietly = TRUE))
#install.packages("devtools")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#if (!requireNamespace("STRINGdb", quietly = TRUE)) {
#BiocManager::install("STRINGdb") }

#if (!requireNamespace("biomaRt", quietly = TRUE)) {
#BiocManager::install("biomaRt") }

```

Install package

`devtools::install_github("martadelfino/fiveP-ml")`

Load package

`library(fivePml)`

Example

```{r}

input_genes <- tibble(hgnc_id = c("HGNC:19743", 
                                  "HGNC:9202",
                                  "HGNC:8653",
                                  "HGNC:6936",
                                  "HGNC:4878"))

result <- get_fiveP(input_genes, binary = TRUE) # this will give you scores of 0 or 1.
head(result)


```

### Other comments

License is GPL (\>= 3) because for PPI, we use STRINGdb, which has GLP 2 license.
