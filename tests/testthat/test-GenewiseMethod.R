testdata <- matrix(data = seq_len(25), nrow = 5, ncol = 5, byrow = TRUE,
                 dimnames = list(paste("Gene", seq_len(5)),
                                 paste("Cell", seq_len(5))))
testdata2 <- testdata
testdata2[2,] <- 0
testdata3 <- testdata
testdata3[3,sample(seq_len(ncol(testdata3)),ncol(testdata3)-1)] <- 0


test_that("MaskData works", {

  # right number of entries per gene is masked
  expect_true(all(rowSums(MaskData(testdata, mask = 0.2) == 0) == 1))

  # handles well only 0 rows
  expect_equivalent(rowSums(MaskData(testdata2, mask = 0.2) == 0),
                    c(1,ncol(testdata2),rep(1,3)))

  # right number of entries per gene is masked with a lot of zeros
  expect_equivalent(rowSums(MaskData(testdata3, mask = 0.2) == 0),
                    c(1,1,ncol(testdata2),1,1))

})

