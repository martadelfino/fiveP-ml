# In order to ensure your input gene list fiveP results are significant,
# you can calculate the GO term similarity of your input gene list, versus the
# input gene list fiveP and random gene lists fiveP results. This means, once
# you have the fiveP for your input gene list, you can do the following steps:
#
# 1. Creating random gene lists, and get fiveP for each.
#
# 2. Calculate GO term enrichment for your input gene list, input gene list
#    fiveP, and random gene list fivePs, at various fiveP thresholds.
#
# 3. Calculate GO term similarity between: input gene list enriched GO terms &
#    input gene list fiveP enriched GO terms, input gene list enriched GO terms
#    & random gene list enriched GO terms. You can do this at various fiveP thresholds.
#
# 4. Plot GO semantic similarity
#
# 5. Calculate Wilcoxon Signed Rank Test for the GO similarity scores obtained at step 3.


# Requierments --------------------------------------------------------------------
#library(tidyverse)
#library(biomaRt)
#library(fivePml)
#library(STRINGdb)
#library(clusterProfiler)
#library(org.Hs.eg.db)
#library(GOSemSim)


# Step 1 --------------------------------------------------------------------------
## You've created the random gene list, you obtained fiveP for the random gene list and for
## the input gene list.

# Input gene list
ndd_ad # Just a df with one hgnc_id column

# Input gene list fiveP
ndd_ad_fivep # get_fiveP() output

# Random gene list fiveP
random_gene_list_fivep # get_fiveP() output


# Step 2 --------------------------------------------------------------------------
## Now, you need to obtain the GO term enrichment for input gene list, input fiveP gene list
## (at many fiveP thresholds), and random gene list fiveP (at many fiveP thresholds)

# defining the fiveP thresholds
filter_positive_thresholds <- function(df) {
  thresholds <- round(seq(0, 1, by = 0.05), 2)

  lapply(thresholds, function(th) {
    if (th == 0) {
      df %>%
        dplyr::filter(if_any(-hgnc_id, ~ .x > 0)) %>%
        dplyr::select(hgnc_id)
    } else if (th == 1) {
      df %>%
        dplyr::filter(if_any(-hgnc_id, ~ .x == 1)) %>%
        dplyr::select(hgnc_id)
    } else {
      df %>%
        dplyr::filter(if_any(-hgnc_id, ~ .x >= th)) %>%
        dplyr::select(hgnc_id)
    }
  }) %>%
    setNames(paste0("threshold_", thresholds))
}

# Apply it to the input gene list and random gene list
random_gene_list_fivep_thresholds <- filter_positive_thresholds(random_gene_list_fivep)
ndd_ad_fivep_thresholds <- filter_positive_thresholds(ndd_ad_fivep)

# GO Term enrichment function
go_enrichment_function <- function(hgnc_id_column) {
  vector_hgnc_id <- pull(hgnc_id_column)
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  mapped_hgnc_id <- getBM(attributes = c("hgnc_id", "entrezgene_id", "external_gene_name"),
                          filters = "hgnc_id",
                          values = vector_hgnc_id,
                          mart = ensembl)
  head(mapped_hgnc_id)
  entrez <- pull(mapped_hgnc_id, entrezgene_id)
  head(entrez)
  go_results <- enrichGO(gene = entrez,
                         OrgDb        = org.Hs.eg.db,
                         keyType      = "ENTREZID",
                         ont          = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.2,
                         readable      = TRUE)
  #head(go_results)
  return(go_results)
}

# Apply it to all thresholds of input gene list fiveP and random gene list fiveP
random_gene_list_fivep_thresholds_goenrich <- lapply(random_gene_list_fivep_thresholds, go_enrichment_function)
ndd_ad_fivep_thresholds_goenrich <- lapply(ndd_ad_fivep_thresholds, go_enrichment_function)

# Apply it to input gene list
ndd_ad_go_enrich <- go_enrichment_function(ndd_ad)


# Step 3 --------------------------------------------------------------------------
## Next, calculate GO term similarity between the GO term enrichment results. The similarity
## is computed betwee: 1. input gene list & input gene list fiveP and
## 2. input gene list & random gene list fiveP

# GO similarity function
hsGO <- godata('org.Hs.eg.db', ont="BP")
go_similarity_function <- function(go_list1, go_list2, p_value = TRUE, p_adjust = 0.05, n_count = 20, measure = "Wang", combine = "BMA") {
  # Work on the enrichment results directly
  go_df1 <- go_list1@result
  go_df2 <- go_list2@result

  # Apply p-value filtering if requested
  if (p_value) {
    go_df1 <- go_df1 %>% dplyr::filter(p.adjust < p_adjust)
    go_df2 <- go_df2 %>% dplyr::filter(p.adjust < p_adjust)
  }

  # Now select top terms
  go_ids1 <- go_df1 %>%
    dplyr::slice_max(order_by = Count, n = n_count) %>%
    dplyr::pull(ID)

  go_ids2 <- go_df2 %>%
    dplyr::slice_max(order_by = Count, n = n_count) %>%
    dplyr::pull(ID)

  # Compute similarity
  sim_score <- mgoSim(go_ids1, go_ids2, semData = hsGO,
                      measure = measure, combine = combine)

  return(sim_score)
}

# GO similarity between gene lists
random_gene_list_fivep_thresholds_gosim <- lapply(random_gene_list_fivep_thresholds_goenrich,
                                                  function(x) go_similarity_function(ndd_ad_go_enrich,x,p_value = TRUE))
ndd_ad_fivep_thresholds_gosim <- lapply(ndd_ad_fivep_thresholds_goenrich,
                                        function(x) go_similarity_function(ndd_ad_go_enrich,x,p_value = TRUE))


# Step 4 --------------------------------------------------------------------------
## Plot the GO similarity results

thresholds <- sprintf("%.2f", seq(0, 1, by = 0.05))
df <- data.frame(
  threshold_label = thresholds,
  ndd_ad_random_fivep = as.numeric(random_gene_list_fivep_thresholds_gosim),
  ndd_ad_ndd_ad_fivep = as.numeric(ndd_ad_fivep_thresholds_gosim)
)

ndd_ad_fivep_random_fivep_gosim_plot <- plot_go_similarity_1v1(df)

# Step 5 --------------------------------------------------------------------------
## Wilcoxon-signed rank test

wilcox.test(ndd_ad_ndd_ad_fivep, ndd_ad_random_fivep, paired = TRUE)

