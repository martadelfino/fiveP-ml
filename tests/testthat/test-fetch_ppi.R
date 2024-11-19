test_that("fetch_ppi returns a data frame", {
  dummy_genes <- data.frame(
    hgnc_id = c("HGNC:1", "HGNC:2"),
    ensembl_gene_id = c("ENSG000001", "ENSG000002")
  )

  result <- fetch_ppi(dummy_genes)

  # Check the result is a data frame
  expect_s3_class(result, "data.frame")

  # Check for expected columns
  expected_cols <- c("protein1_hgnc_id", "protein1_string_id", "protein2_hgnc_id", "protein2_string_id", "combined_score")
  expect_true(all(expected_cols %in% colnames(result)))

  # Check that the function handles empty input gracefully
  empty_result <- fetch_ppi(data.frame(hgnc_id = character(), ensembl_gene_id = character()))
  expect_equal(nrow(empty_result), 0)
})
