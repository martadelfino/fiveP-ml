Best Practices for Testing

Test Edge Cases:

What happens if the input dataframe is empty?
What if some required columns are missing?
What if the data contains unexpected values?


Test Expected Output:

Ensure the function returns the expected structure (e.g., data frame with correct columns).


Automate Testing:

Running devtools::check() will automatically run all tests during package checks.
