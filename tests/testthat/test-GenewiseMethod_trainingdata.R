context("GenewiseMethod functions - data splitting and masking")

testdata <- matrix(data = seq_len(25), nrow = 5, ncol = 5, byrow = TRUE,
    dimnames = list(paste("Gene", seq_len(5)), paste("Cell", seq_len(5))))
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

  # choice of entries to mask is random
  expect_false(identical(MaskData(testdata, mask = 0.2),
                         MaskData(testdata, mask = 0.2)))

})


test_that("SplitData works", {

  # right number of samples are included in the training process
  expect_true(ncol(SplitData(testdata, ratio = .2)) == 1)
  expect_true(ncol(SplitData(testdata, ratio = .4)) == 2)

  # throws error when ratio is too small and no samples are selected
  expect_error(SplitData(testdata, ratio = .1))

  # writes when supposed to
  dir.create(dir1 <- file.path(tempdir(), "testdir"))
  cur_dir <- getwd(); setwd(dir1)
  expect_false(file.exists("training.txt"))
  expect_false(file.exists("validation.txt"))
  SplitData(testdata, write.to.file = T)
  expect_true(file.exists("training.txt"))
  # outputs correct validation when asked to
  out <- SplitData(testdata, train.only = FALSE, write.to.file = TRUE)
  expect_true(file.exists("validation.txt"))
  setwd(cur_dir); unlink(dir1, recursive = TRUE)

  # outputs matrix
  expect_is(SplitData(testdata), "matrix")

  # choice of split is random
  expect_false(identical(colnames(SplitData(testdata, ratio = .8)),
                         colnames(SplitData(testdata, ratio = .8))))

})


test_that("CreateTrainData works", {

  # outputs appropriate list
  out <- CreateTrainData(testdata)
  expect_is(out, "list")
  expect_true(length(out) == 2)
  expect_identical(names(out), c("train","mask"))
  expect_identical(dim(out[["train"]]),dim(out[["mask"]]))

})
