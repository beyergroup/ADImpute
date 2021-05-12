context("GenewiseMethod functions - MSE computation")

# Gene 1: expressed everywhere and masked everywhere
# Gene 2: not expressed and not masked
# Gene 3: expressed in 4/5 and masked in 2/4
testreal <- matrix(c(rep(1,5),rep(0,5),2,2,4,0,4), nrow = 3, ncol = 5,
    byrow = TRUE, dimnames = list(paste("Gene", seq_len(3)),
        paste("Cell", seq_len(5))))
testmask <- testreal
testmask[1,] <- 0 ; testmask[3,1] <- testmask[3,3] <- 0
# MSE of 0 for gene 1, NA for gene 2 and variable for gene 3
testimp <- testreal
testimp[3,3] <- 0
expected_MSE <- mean(c((4-0)^2,0)) # baseline MSE; NA for other methods
testimp1 <- testimp; testimp1[3,2] <- 10 # should not change MSE (not masked)

test_that("ComputeMSEGenewise works", {

    # outputs the correct value for case gene 1
    expect_equal(ComputeMSEGenewise(real = testreal[1,],
        masked = (testreal[1,] != 0) & (testmask[1,] == 0),
        imputed = testimp[1,], baseline = TRUE), 0)
    expect_equal(ComputeMSEGenewise(real = testreal[1,],
        masked = (testreal[1,] != 0) & (testmask[1,] == 0),
        imputed = testimp[1,], baseline = FALSE), 0)

    # outputs NA when no entry was masked (gene 2)
    expect_true(is.na(ComputeMSEGenewise(real = testreal[2,],
        masked = (testreal[2,] != 0) & (testmask[2,] == 0),
        imputed = testimp[2,], baseline = TRUE)))
    expect_true(is.na(ComputeMSEGenewise(real = testreal[2,],
        masked = (testreal[2,] != 0) & (testmask[2,] == 0),
        imputed = testimp[2,], baseline = FALSE)))

    # separates correctly between Baseline and other methods
    # for Baseline, look at all masked values, even imputation of 0:
    expect_equal(ComputeMSEGenewise(real = testreal[3,],
        masked = (testreal[3,] != 0) & (testmask[3,] == 0),
        imputed = testimp[3,], baseline = TRUE), expected_MSE)
    # for other methods, only look at imputations > 0
    expect_equal(ComputeMSEGenewise(real = testreal[3,],
        masked = (testreal[3,] != 0) & (testmask[3,] == 0),
        imputed = testimp[3,], baseline = FALSE), 0)

    # does not evaluate non-masked values (MSE should be the same as before)
    expect_equal(ComputeMSEGenewise(real = testreal[3,],
        masked = (testreal[3,] != 0) & (testmask[3,] == 0),
        imputed = testimp1[3,], baseline = TRUE), expected_MSE)
    expect_equal(ComputeMSEGenewise(real = testreal[3,],
        masked = (testreal[3,] != 0) & (testmask[3,] == 0),
        imputed = testimp1[3,], baseline = FALSE), 0)
})


testimplist <- list("Baseline" = testimp, "DrImpute" = testimp1)

test_that("ComputeMSE works", {

    # outputs matrix
    expect_is(ComputeMSE(real = testreal, masked = testmask,
        imputed = testimplist), "matrix")
    # of ncol = length(imputed)
    expect_equal(ncol(ComputeMSE(real = testreal, masked = testmask,
        imputed = testimplist)), length(testimplist))
    # and nrow = nrow(real)
    expect_equal(nrow(ComputeMSE(real = testreal, masked = testmask,
        imputed = testimplist)), nrow(testreal))

    # correct output matching ComputeMSEGenewise
    expected_out <- lapply(testimplist, function(m) vapply(rownames(testreal),
        function(g) ComputeMSEGenewise(real = testreal[g,],
            masked = (testreal[g,] != 0) & (testmask[g,] == 0),
            imputed = m[g,], baseline = identical(testimplist$Baseline, m)),
        FUN.VALUE = 1))
    expected_out <- do.call(cbind, expected_out)
    expect_identical(ComputeMSE(real = testreal, masked = testmask,
        imputed = testimplist), expected_out)

    # case for 1 method only
    expect_is(ComputeMSE(real = testreal, masked = testmask,
        imputed = testimplist[1]), "matrix")

    # throws error if colnames and rownames don't match between objects
    testmask1 <- testmask2 <- testmask
    rownames(testmask1) <- seq_len(nrow(testmask1))
    colnames(testmask2) <- seq_len(ncol(testmask2))
    expect_error(ComputeMSE(real = testreal, masked = testmask1,
        imputed = testimplist))
    expect_error(ComputeMSE(real = testreal, masked = testmask2,
        imputed = testimplist))
})


test_that("CrossValidateImputation works", {

    mse_mat <- CrossValidateImputation(data = ADImpute::demo_data)

    # returns MSE matrix, of the right dimensions (ncol = # tested methods)
    expect_is(mse_mat, "matrix")
    expect_equal(ncol(mse_mat),3)

})


test_that("AggregateMSE works", {

    # works with array of 3rd dim 1

    # does the median correctly

    # does the mean correctly

    # returns a matrix

})


test_that("ChooseMethod works", {

    # outputs a character vector

    # picks the lowest MSE method

    # gets rid of MSEs non-NA in less than 2 methods


    # it writes when supposed to


    # outputs message
})

