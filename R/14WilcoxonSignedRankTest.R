#' Function for calculating Wilcoxon signed-rank test between two sets of GO similarity scores
#'
#' @param vector A numeric vector of GO similarity scores for the input gene list
#' @param list A list of numeric vectors, each representing GO similarity scores for a random gene list
#' @importFrom magrittr %>%
#' @return A data frame with Wilcoxon test results (p-values and statistics) for each comparison
#' @export
run_wilcox_tests <- function(vector, list) {
  # Run Wilcoxon test for each element in dd
  results <- lapply(seq_along(list), function(i) {
    test <- wilcox.test(vector, list[[i]], paired = TRUE)
    data.frame(
      list_index = i,
      p_value = test$p.value,
      statistic = test$statistic
    )
  })

  # Combine into a single data frame
  results_df <- do.call(rbind, results)
  return(results_df)
}
