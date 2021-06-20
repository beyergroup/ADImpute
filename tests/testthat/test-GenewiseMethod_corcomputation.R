context("GenewiseMethod functions - correlation computation")

# Gene 1: expressed everywhere and masked everywhere
# Gene 2: not expressed and not masked
# Gene 3: expressed in 4/5 and masked in 2/4
testreal <- matrix(c(rep(1,10),rep(2,10),rep(0,20),rep(c(2,2,4,4),5)),
    nrow = 3, ncol = 20, byrow = TRUE,
    dimnames = list(paste("Gene", seq_len(3)), paste("Cell", seq_len(20))))
testmask <- testreal
testmask[1,] <- 0 ; testmask[3,c(1,3,5,7,10,12,15,17,18,19)] <- 0
# correlation of 1 for gene 1, NA for gene 2 and variable for gene 3
testimp <- testreal
testimp[3,c(3,7,12)] <- 1
expected_cor <- cor.test(testreal[3,(testmask[3,] == 0) & (testimp[3,] != 0)],
    testimp[3,(testmask[3,] == 0) & (testimp[3,] != 0)])$estimate
testimp1 <- testimp; testimp1[3,2] <- 10 # should not change cor (not masked)

test_that("ComputeCorGenewise works", {

  # outputs the correct value for case gene 1
  expect_equivalent(ComputeCorGenewise(real = testreal[1,],
      masked = (testreal[1,] != 0) & (testmask[1,] == 0),
      imputed = testimp[1,], baseline = FALSE), 1)

  # outputs NA for Baseline
  expect_equivalent(ComputeCorGenewise(real = testreal[1,],
      masked = (testreal[1,] != 0) & (testmask[1,] == 0),
      imputed = testimp[1,], baseline = TRUE), NA)

  # outputs NA when no entry was masked (gene 2)
  expect_true(is.na(ComputeCorGenewise(real = testreal[2,],
      masked = (testreal[2,] != 0) & (testmask[2,] == 0),
      imputed = testimp[2,], baseline = FALSE)))

  # outputs correct result
  expect_equivalent(ComputeCorGenewise(real = testreal[3,],
      masked = (testreal[3,] != 0) & (testmask[3,] == 0),
      imputed = testimp[3,], baseline = FALSE), expected_cor)

  # does not evaluate non-masked values (cor should be the same as before)
  expect_equivalent(ComputeCorGenewise(real = testreal[3,],
      masked = (testreal[3,] != 0) & (testmask[3,] == 0),
      imputed = testimp1[3,], baseline = FALSE), expected_cor)
})


testimplist <- list("Network" = testimp, "DrImpute" = testimp1)

test_that("ComputeCor works", {

  # outputs matrix
  expect_is(ComputeCor(real = testreal, masked = testmask,
      imputed = testimplist), "matrix")
  # of ncol = length(imputed)
  expect_equal(ncol(ComputeCor(real = testreal, masked = testmask,
      imputed = testimplist)), length(testimplist))
  # and nrow = nrow(real)
  expect_equal(nrow(ComputeCor(real = testreal, masked = testmask,
      imputed = testimplist)), nrow(testreal))

  # correct output matching ComputeCorGenewise
  expected_out <- lapply(testimplist, function(m) vapply(rownames(testreal),
      function(g) ComputeCorGenewise(real = testreal[g,],
          masked = (testreal[g,] != 0) & (testmask[g,] == 0),
          imputed = m[g,], baseline = FALSE), FUN.VALUE = 1))
  expected_out <- do.call(cbind, expected_out)
  expect_identical(ComputeCor(real = testreal, masked = testmask,
      imputed = testimplist), expected_out)

  # case for 1 method only
  expect_is(ComputeCor(real = testreal, masked = testmask,
      imputed = testimplist[1]), "matrix")

  # throws error if colnames and rownames don't match between objects
  testmask1 <- testmask2 <- testmask
  rownames(testmask1) <- seq_len(nrow(testmask1))
  colnames(testmask2) <- seq_len(ncol(testmask2))
  expect_error(ComputeCor(real = testreal, masked = testmask1,
      imputed = testimplist))
  expect_error(ComputeCor(real = testreal, masked = testmask2,
      imputed = testimplist))
})

