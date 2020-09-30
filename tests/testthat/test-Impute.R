test_that("ImputeBaseline works", {

    out <- ImputeBaseline(data = ADImpute::demo_data)

    # test output is a matrix or dgCMatrix
    expect_is(out, c("matrix"))
    # test output dimensions are the same as in
    expect_equivalent(dim(ADImpute::demo_data), dim(out))

    # test that it throws an appropriate error when no data
    expect_error(ImputeBaseline(data = NULL),
                 "Please provide an input data matrix.")

    # test number of zeros decreases
    expect_true(sum(ADImpute::demo_data == 0) > sum(out == 0))
})


test_that("ImputeDrImpute works", {

    out <- ImputeDrImpute(data = ADImpute::demo_data)

    # test output is a matrix or dgCMatrix
    expect_is(out, c("matrix"))
    # test output dimensions are the same as in
    expect_equivalent(dim(ADImpute::demo_data), dim(out))

    # test that it throws an appropriate error when no data
    expect_error(ImputeDrImpute(data = NULL),
                 "Please provide an input data matrix.")

    # test number of zeros decreases
    expect_true(sum(ADImpute::demo_data == 0) > sum(out == 0))

})


test_that("ImputeNetwork works", {

    out <- ImputeNetwork(data = ADImpute::demo_data,
                         net.coef = ADImpute::demo_net, cores = 2)

    # test output is a matrix or dgCMatrix
    expect_is(out, c("matrix","dgCMatrix"))
    # test output dimensions are the same as in
    expect_equivalent(dim(ADImpute::demo_data), dim(out))

    # test that it throws an appropriate error when no data
    expect_error(ImputeNetwork(data = NULL, net.coef = ADImpute::demo_net,
                               cores = 2),
                 "Please provide an input data matrix.")
    # when no network coefficients
    expect_error(ImputeNetwork(data = ADImpute::demo_data, net.coef = NULL,
                               cores = 2),
                 "Please provide valid network coefficients.")

    # test number of zeros decreases
    expect_true(sum(ADImpute::demo_data == 0) >
                  sum(out == 0))
})


