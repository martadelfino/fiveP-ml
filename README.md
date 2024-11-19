# fiveP-R

### Install devtools

`install.packages("devtools")`

### Install the package

`devtools::install_github("martadelfino/fiveP-R")`

### Load the package

`library(fiveP)`

Example:

```{r}

# Reading the file 
gene_classes <- readr::read_delim('CodeReview_data.txt', delim = '\t',  
                                      show_col_types = FALSE)

```

```{r}

# Getting the HGNC IDs
AR <- gene_classes %>%
  dplyr::filter(ndd_ar_classes == 'positive') %>%
  dplyr::select(hgnc_id) #%>% dplyr::slice_sample(n = 5)
```

```{r}

AR_results <- get_fiveP(AR)


```

```{r}
print(AR_results)
```
