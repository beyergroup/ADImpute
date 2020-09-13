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
  expect_is(ArrangeData(testdata, network.coefficients = testnet), "list")

  # returns the proper names and classes
  expect_equal(names(ArrangeData(testdata, network.coefficients = testnet)),
               c("data","network","O"))
  expect_is(ArrangeData(testdata, network.coefficients = testnet)$O,
            "numeric")
  expect_is(ArrangeData(testdata, network.coefficients = testnet)$network,
            "matrix")
  expect_is(ArrangeData(testdata, network.coefficients = testnet)$data,
            "matrix")

  arranged <- ArrangeData(testdata, network.coefficients = testnet)
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
