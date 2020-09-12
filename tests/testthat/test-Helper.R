testdata1 <- matrix(data = as.numeric(NA), nrow = 5, ncol = 5,
                    dimnames = list(paste("Gene", seq_len(5)),
                                    paste("Cell", seq_len(5))))
testdata2 <- matrix(data = seq_len(25), nrow = 5, ncol = 5, byrow = T)
testdata3 <- matrix(data = seq_len(25), nrow = 5, ncol = 5, byrow = T,
                    dimnames = list(paste("Gene", seq_len(5)),
                                    paste("Cell", seq_len(5))))
testdata4 <- testdata3
testdata4[3,2] <- NA
testdata5 <- testdata4
storage.mode(testdata5) <- "character"
testdata6 <- as.data.frame(testdata4)

output4 <- testdata4
output4[is.na(output4)] <- 0

test_that("DataCheck_Matrix works", {

  # Class of output is a matrix
  expect_is(DataCheck_Matrix(testdata3), "matrix")
  expect_is(DataCheck_Matrix(testdata4), "matrix")
  expect_is(DataCheck_Matrix(testdata6), "matrix")

  # Throws error when not numeric
  expect_error(DataCheck_Matrix(testdata5), "Input must be numeric.\n")
  # when no dimnames
  expect_error(DataCheck_Matrix(testdata2), "Input has NULL dimnames.\n")
  # when all NULL
  expect_error(DataCheck_Matrix(testdata1), "Input has only NAs.\n")

  # Converts NA to 0
  expect_equal(DataCheck_Matrix(testdata4), output4)
})
