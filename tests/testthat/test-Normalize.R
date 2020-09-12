testdata <- matrix(data = seq_len(25), nrow = 5, ncol = 5, byrow = T,
                   dimnames = list(paste("Gene", seq_len(5)),
                                   paste("Cell", seq_len(5))))
testdata[3,2] <- NA

test_that("NormalizeRPM works", {

  # All colSums are 1M
  expect_equivalent(colSums(NormalizeRPM(testdata)), rep(10^6, ncol(testdata)))

  # Class of output is a matrix
  expect_is(NormalizeRPM(testdata), "matrix")

  # Log works
  expect_equal(NormalizeRPM(testdata, log = TRUE),
               log2(1+NormalizeRPM(testdata)))

  # Log with no pseudocount throws warning
  expect_warning(NormalizeRPM(testdata, log = TRUE, pseudo.count = 0),
                 "Using 0 pseudocount: Inf may be generated.\n")

  # Scale works
  expect_equal(NormalizeRPM(testdata, scale = 100),
               NormalizeRPM(testdata)/100)
})
