
test_that("CheckArguments works", {

  do = "Baseline"; write = FALSE; train.ratio = .7; sce = NULL; data = NULL;
  drop_thre = 4; labels = NULL; tr.length = NULL; bulk = NULL; mask.ratio = -1;
  folds = 2; scale = 2; pseudocount = 1; cell.clusters = NULL; cores = 1;
  type = "oops"; net.implementation = "wrong"; train.only = "yes"
  expect_error(CheckArguments(environment()))

})


testdata1 <- matrix(data = as.numeric(NA), nrow = 5, ncol = 5,
                    dimnames = list(paste("Gene", seq_len(5)),
                                    paste("Cell", seq_len(5))))
testdata2 <- matrix(data = seq_len(25), nrow = 5, ncol = 5, byrow = TRUE)
testdata3 <- matrix(data = seq_len(25), nrow = 5, ncol = 5, byrow = TRUE,
                    dimnames = list(paste("Gene", seq_len(5)),
                                    paste("Cell", seq_len(5))))
testdata4 <- testdata3
testdata4[3,2] <- NA
testdata5 <- testdata4
storage.mode(testdata5) <- "character"
testdata6 <- as.data.frame(testdata4)

output4 <- testdata4
output4[is.na(output4)] <- 0

test_that("DataCheck_Matrix works", {

  # Class of output is a matrix
  expect_is(DataCheck_Matrix(testdata3), "matrix")
  expect_is(DataCheck_Matrix(testdata4), "matrix")
  expect_is(DataCheck_Matrix(testdata6), "matrix")

  # Throws error when not numeric
  expect_error(DataCheck_Matrix(testdata5), "Input must be numeric.\n")
  # when no dimnames
  expect_error(DataCheck_Matrix(testdata2), "Input has NULL dimnames.\n")
  # when all NULL
  expect_error(DataCheck_Matrix(testdata1), "Input has only NAs.\n")

  # Converts NA to 0
  expect_equal(DataCheck_Matrix(testdata4), output4)
})


testtrlength1 <- unname(ADImpute::transcript_length)
testtrlength2 <- ADImpute::transcript_length[,2, drop = FALSE]
testtrlength3 <- ADImpute::transcript_length
testtrlength3$hgnc_symbol <- seq_len(nrow(testtrlength3))
testtrlength4 <- ADImpute::transcript_length
testtrlength4$hgnc_symbol <- rep("", nrow(testtrlength4))

test_that("DataCheck_TrLength works", {

  # Class of output is a data.frame
  expect_is(DataCheck_TrLength(ADImpute::transcript_length),"data.frame")
  expect_is(DataCheck_TrLength(as.matrix(ADImpute::transcript_length)),
            "data.frame")

  # Throws error when required colnames not present
  expect_error(DataCheck_TrLength(testtrlength1),
      paste("Transcript length data must contain the following colnames:",
          "hgnc_symbol, transcript_length\n"))
  expect_error(DataCheck_TrLength(testtrlength2),
      paste("Transcript length data must contain the following colnames:",
          "hgnc_symbol, transcript_length\n"))
  # when gene names are not character/factor
  expect_error(DataCheck_TrLength(testtrlength3),
      "hgnc_symbol column must be character/factor.\n")
  # when there are no rows
  expect_error(DataCheck_TrLength(testtrlength3[0,]),
      "Not enough rows in transcript length data.\n")
  # when all gene symbols are empty strings
  expect_error(DataCheck_TrLength(testtrlength4),
      "Not enough non-empty gene symbols in transcript length data.\n")

  # Converts character lengths to numeric
  expect_is(DataCheck_TrLength(as.matrix(
    ADImpute::transcript_length))$transcript_length, "numeric")

  # Removes empty gene symbols
  expect_false(any(as.character(DataCheck_TrLength(
    ADImpute::transcript_length)$hgnc_symbol) == ""))
})
