# fiveP-ml

This repo contains the scripts to obtain the '5 Ps' of input genes: protein complex, protein family, paralogues, pathways, and protein-protein interactions (PPI). The output is a data frame with all protein coding genes and a binary score (0 or 1) describing whether each gene is in the same 'P' as any of the input genes.

### Installation

Install devtools

`install.packages("devtools")`

Install dependencies

```{r}
if (!requireNamespace("tidyverse", quietly = TRUE))
install.packages("tidyverse")

if (!requireNamespace("devtools", quietly = TRUE))
install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

if (!requireNamespace("STRINGdb", quietly = TRUE)) {
BiocManager::install("STRINGdb") }

if (!requireNamespace("biomaRt", quietly = TRUE)) {
BiocManager::install("biomaRt") }

```

Install package

`devtools::install_github("martadelfino/fiveP-ml")`

Load package

`library(fivePml)`

### Formula for each P 

Paralogue and PPI are scores are calculated differently than protein complex, protein family, and pathway scores.

#### Protein complexes, Protein families, Pathways

Using pathways as the example P: For all protein coding genes, a score is calculated by looking only at pathways of the input genes. The formula for the pathway score for each protein coding gene:

( \# of unique input genes in input gene pathways ) / ( \# of unique genes in input gene pathways )

#### Paralogues, PPI

Using PPI as the example: For all protein coding genes, a score is calculated by dividing the number of input genes by all genes in the PPI. The formula for the PPI score for each protein coding gene:

( \# of input genes that gene interacts with ) / ( \# of genes it interacts with )

#### Final results

The result is each protein coding gene is assigned a score from 0 to 1, for each 'P'. For the genes without data, a 'NA' is left in its place. The output is a data frame with HGNC IDs, and columns for each P.

(Add back option to keep the ratios).
