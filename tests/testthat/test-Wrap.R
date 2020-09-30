methods <- EvaluateMethods(data = ADImpute::demo_data,
                            net.coef = ADImpute::demo_net, cores = 2)

test_that("EvaluateMethods works", {

    # output is a character vector
    expect_is(methods, "character")
    # entries are within appropriate methods
    expect_true(all(methods %in% c("Baseline","Network","DrImpute")))

    # test that it throws an appropriate error when no data
    expect_error(EvaluateMethods(data = NULL, net.coef = ADImpute::demo_net,
                                    cores = 2),
                 "'data' must have non-NULL value")
    # when no network coefficients
    expect_error(EvaluateMethods(data = ADImpute::demo_data, net.coef = NULL,
                                    cores = 2),
                 "'net.coef' must have non-NULL value")
    # when no valid methods
    expect_error(EvaluateMethods(data = ADImpute::demo_data,
                    net.coef = ADImpute::demo_net, cores = 2, do = NULL),
                 "Please provide appropriate imputation methods")
    # when wrong methods are passed
    expect_error(EvaluateMethods(data = ADImpute::demo_data,
                    net.coef = ADImpute::demo_net, cores = 2, do = "wrong"),
                   paste0("Please provide at least one supported method"))
})


imputation <- Impute(data = ADImpute::demo_data, method.choice = methods,
                     net.coef = ADImpute::network.coefficients, cores = 2)

test_that("Impute works", {

    # output is a list
    expect_is(imputation, "list")
    # elements are matrices
    lapply(imputation, function(m) expect_is(m, "matrix"))

    # test that it throws an appropriate error when no data
    expect_error(Impute(data = NULL, net.coef = ADImpute::demo_net,
                        method.choice = methods, cores = 2),
                 "'data' must have non-NULL value")
    # when no network coefficients
    expect_error(Impute(data = ADImpute::demo_data, net.coef = NULL,
                        method.choice = methods, cores = 2),
                 "'net.coef' must have non-NULL value")
    # when wrong methods are passed
    expect_error(Impute(data = ADImpute::demo_data,
                        net.coef = ADImpute::demo_net,
                        method.choice = methods, cores = 2, do = "wrong"),
                 "Please provide at least one supported method")
    expect_warning(Impute(data = ADImpute::demo_data,
                          net.coef = ADImpute::demo_net,
                          method.choice = methods, cores = 2, do = c("Baseline",
                                                                     "wrong")),
                   paste0("The following methods were detected as input but ",
                          "are not supported and will be ignored: 'wrong'"))
    # warning for SCRABBLE and scImpute
    expect_warning(tryCatch(Impute(data = ADImpute::demo_data,
                                   method.choice = methods,
                                   net.coef = ADImpute::network.coefficients,
                                   cores = 2, do = "scImpute")))
    expect_warning(tryCatch(Impute(data = ADImpute::demo_data,
                                   method.choice = methods,
                                   net.coef = ADImpute::network.coefficients,
                                   cores = 2, do = "SCRABBLE")))

})
