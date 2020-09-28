
test_that("GetDropoutProbabilities works", {

  res = GetDropoutProbabilities(data = ADImpute::demo_data, thre = .2,
                                 ncores = 2, cell.clusters = 2)
  expect_is(res, "matrix")
  expect_true(all(na.omit(res) >= 0) & all(na.omit(res) <= 1))

  # handling only one cell cluster
  res = GetDropoutProbabilities(data = ADImpute::demo_data, thre = .2,
                                ncores = 2, cell.clusters = 1)
  expect_is(res, "matrix")
  expect_true(all(na.omit(res) >= 0) & all(na.omit(res) <= 1))

  # handling labels
  res = GetDropoutProbabilities(data = ADImpute::demo_data, thre = .2, ncores =
                                  2, labels = c(rep("A",25), rep("B",25)))
  expect_is(res, "matrix")
  expect_true(all(na.omit(res) >= 0) & all(na.omit(res) <= 1))
})

