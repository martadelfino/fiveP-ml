Best Practices for Testing

Test Edge Cases:

What happens if the input dataframe is empty?
What if some required columns are missing?
What if the data contains unexpected values?


Test Expected Output:

Ensure the function returns the expected structure (e.g., data frame with correct columns).


Automate Testing:

Running devtools::check() will automatically run all tests during package checks.




A typical test file uses testthat syntax like:

r
Copy code
test_that("my_function works as expected", {
  x <- my_function(1, 2)
  expect_equal(x, 3)
})

test_that("my_function handles errors properly", {
  expect_error(my_function("a", "b"))
})
Note the main testthat functions:

test_that("description", { ... }) blocks group related expectations.
expect_equal(), expect_error(), expect_true(), etc., check that the function behaves as expected.
