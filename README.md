---
editor_options: 
  markdown: 
    wrap: 72
---

# fiveP-ml

This repo contains the scripts to obtain the '5 Ps' of a number of input
genes: protein complex, protein family, paralogues, pathways, and
protein-protein interactions. The output is a data frame with all
protein coding genes and a binary score (0 or 1) describing whether each
gene is in the same 'P' as any of the input genes.

### Install devtools

`install.packages("devtools")`

### Install dependencies

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

### Install the package

`devtools::install_github("martadelfino/fiveP-ml")`

### Load the package

`library(fiveP)`
