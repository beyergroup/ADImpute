methods <- EvaluateMethods(data = ADImpute::demo_data,
    net.coef = ADImpute::demo_net, cores = 2)

test_that("EvaluateMethods works with matrix input", {

    # output is a character vector
    expect_is(methods, "character")
    # entries are within appropriate methods
    expect_true(all(methods %in% c("Baseline","Network","DrImpute")))

    # test that it throws an appropriate error when no data
    expect_error(EvaluateMethods(data = NULL, net.coef = ADImpute::demo_net,
                                    cores = 2))
    # when no network coefficients
    expect_error(EvaluateMethods(data = ADImpute::demo_data, net.coef = NULL,
                                    cores = 2))
    # when no valid methods
    expect_error(EvaluateMethods(data = ADImpute::demo_data,
                    net.coef = ADImpute::demo_net, cores = 2, do = NULL),
                 "Please provide appropriate imputation methods")
    # when wrong methods are passed
    expect_error(EvaluateMethods(data = ADImpute::demo_data,
                    net.coef = ADImpute::demo_net, cores = 2, do = "wrong"),
                   paste0("Please provide at least one supported method"))
})

sce <- NormalizeRPM(sce = ADImpute::demo_sce)
sce <- EvaluateMethods(sce = sce, net.coef = ADImpute::demo_net, cores = 2)

test_that("EvaluateMethods works with SingleCellExperiment input", {

    # output is a character vector in the SingleCellExperiment
    expect_is(sce, "SingleCellExperiment")
    expect_is(SingleCellExperiment::int_elementMetadata(sce)$ADImpute$method,
        "character")
    # entries are within appropriate methods
    methods <- SingleCellExperiment::int_elementMetadata(sce)$ADImpute$method
    expect_true(all(stats::na.omit(methods) %in% c("Baseline","Network",
        "DrImpute")))

    # when no network coefficientsD
    expect_error(EvaluateMethods(sce = sce, net.coef = NULL, cores = 2))
    # when no valid methods
    expect_error(EvaluateMethods(sce = sce, net.coef = ADImpute::demo_net,
        cores = 2, do = NULL), "Please provide appropriate imputation methods")
    # when wrong methods are passed
    expect_error(EvaluateMethods(sce = sce, net.coef = ADImpute::demo_net,
        cores = 2, do = "wrong"),
        paste0("Please provide at least one supported method"))

})

imputation <- Impute(data = ADImpute::demo_data, method.choice = methods,
                     net.coef = ADImpute::network.coefficients, cores = 2)

test_that("Impute works with matrix input", {

    # output is a list
    expect_is(imputation, "list")
    # elements are matrices
    lapply(imputation, function(m) expect_is(m, "matrix"))

    # test that it throws an appropriate error when no data
    expect_error(Impute(data = NULL, net.coef = ADImpute::demo_net,
                        method.choice = methods, cores = 2))
    # when no network coefficients
    expect_error(Impute(data = ADImpute::demo_data, net.coef = NULL,
                        method.choice = methods, cores = 2))
    # when wrong methods are passed
    expect_error(Impute(data = ADImpute::demo_data,
                        net.coef = ADImpute::demo_net,
                        method.choice = methods, cores = 2, do = "wrong"))
    expect_warning(Impute(data = ADImpute::demo_data,
                          net.coef = ADImpute::demo_net,
                          method.choice = methods, cores = 2, do = c("Baseline",
                                                                     "wrong")))
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


imputation <- Impute(data = ADImpute::demo_data, do = "Baseline",
                     cores = 2, true.zero.thr = .3)

test_that("Impute works with biological zero determination", {

    # output is a list
    expect_is(imputation, "list")
    # elements are appropriate
    expect_is(imputation$imputations, "list")
    expect_is(imputation$zerofiltered, "list")
    expect_is(imputation$dropoutprobabilities, "matrix")

    # entries that are zero are less after setting biological zeros
    expect_true(sum(imputation$imputations$Baseline == 0) <
                    sum(imputation$zerofiltered$Baseline == 0))

    # works when providing labels
    expect_is(Impute(data = ADImpute::demo_data, do = "Baseline", cores = 2,
           labels = c(rep("A",20),rep("B",30)), true.zero.thr = .3), "list")
})

sce <- Impute(sce = sce, cores = 2)

test_that("Impute works with SingleCellExperiment input", {

    # output is a list
    expect_is(sce, "SingleCellExperiment")
    # elements are matrices
    lapply(SummarizedExperiment::assays(sce),
        function(m) expect_is(m, c("matrix", "dgCMatrix")))

    # when no network coefficients
    expect_error(Impute(sce = sce, net.coef = NULL, cores = 2))
    # when wrong methods are passed
    expect_error(Impute(sce = sce, net.coef = ADImpute::demo_net, cores = 2,
        do = "wrong"), "Please provide at least one supported method")
    expect_warning(Impute(sce = sce, net.coef = ADImpute::demo_net,
        cores = 2, do = c("Baseline", "wrong")),
        paste0("The following methods were detected as input but ",
            "are not supported and will be ignored: 'wrong'"))
    # warning for SCRABBLE and scImpute
    expect_warning(tryCatch(Impute(sce = sce, cores = 2, do = "scImpute",
        net.coef = ADImpute::network.coefficients)))
    expect_warning(tryCatch(Impute(sce = sce, cores = 2, do = "SCRABBLE",
        net.coef = ADImpute::network.coefficients)))

})
