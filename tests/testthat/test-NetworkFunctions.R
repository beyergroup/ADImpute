testdata <- matrix(data = seq_len(25), nrow = 5, ncol = 5, byrow = TRUE,
                   dimnames = list(paste("Gene", seq_len(5)),
                                   paste("Cell", seq_len(5))))
testnet <- matrix(data = sample(rnorm(1000), 30), nrow = 5, ncol = 6,
                  dimnames = list(paste("Gene", c(1,2,3,6,7)),
                                  c("O", paste("Gene", c(1,3,4,7,8)))))


test_that("ArrangeData works", {

  # throws error when path to coefficients is not given
  expect_error(ArrangeData(testdata),
               "Please provide a valid path for network coefficients.\n")
  # when path to coefficients is invalid
  expect_error(ArrangeData(testdata, network.path = "123/abc.zip"),
               "Please provide a valid path for network coefficients.\n")

  # returns a list
  expect_is(ArrangeData(testdata, net.coef = testnet), "list")

  # returns the proper names and classes
  expect_equal(names(ArrangeData(testdata, net.coef = testnet)),
               c("data","network","O"))
  expect_is(ArrangeData(testdata, net.coef = testnet)$O,
            "numeric")
  expect_is(ArrangeData(testdata, net.coef = testnet)$network,
            "matrix")
  expect_is(ArrangeData(testdata, net.coef = testnet)$data,
            "matrix")

  arranged <- ArrangeData(testdata, net.coef = testnet)
  # length of intercept is the same as number of network targets
  expect_equal(length(arranged$O), nrow(arranged$network))

  # data is limited to genes that are predictors or targets
  expect_true(all(rownames(arranged$data) %in% c(rownames(arranged$network),
                                                 colnames(arranged$network))))
  # network is limited to genes present in data
  expect_true(all(union(rownames(arranged$network),
                        colnames(arranged$network)) %in%
                    rownames(arranged$data)))

  # all predictors in data are maintained
  expect_true(all(intersect(colnames(testnet),rownames(data)) %in%
                    colnames(arranged$network)))
  # all targets in data are maintained
  expect_true(all(intersect(rownames(testnet),rownames(data)) %in%
                    rownames(arranged$network)))

})

testdata1 <- testdata
testdata1[3,2] <- 0

test_that("CenterData works", {

  # returns a list of centered data and centers
  expect_is(CenterData(testdata), "list")
  expect_is(CenterData(testdata)$data, "matrix")
  expect_is(CenterData(testdata)$center, "numeric")

  # returns correct centered data and centers
  expect_equivalent(CenterData(testdata)$data, testdata - rowMeans(testdata))
  expect_equivalent(CenterData(testdata)$center, rowMeans(testdata))

  # handles zeros as instructed
  expect_equivalent(CenterData(testdata1)$center[3],
                    mean(testdata1[3,testdata1[3,] != 0]))
  expect_equivalent(CenterData(testdata1)$data[3,],
                    testdata1[3,] - mean(testdata1[3,testdata1[3,] != 0]))
  expect_equivalent(CenterData(testdata1)$center[-3],
                    rowMeans(testdata)[-3])
  expect_equivalent(CenterData(testdata1)$data[-3,],
                    CenterData(testdata)$data[-3,])
})


# cell_expression <- ADImpute::demo_data[,1]
# net <- ADImpute::demo_net
# net <- net[intersect(rownames(net),names(cell_expression)),
#            intersect(colnames(net),names(cell_expression))]
# cell_expression <- cell_expression[intersect(unique(unlist(dimnames(net))),
#                                              names(cell_expression))]
#
#
# test_that("ImputeNonPredictiveDropouts works", {
#
#   expect_true(length(ImputeNonPredictiveDropouts(net, cell_expression)) ==
#                 nrow(net))
# })


