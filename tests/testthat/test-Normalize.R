testdata <- matrix(data = seq_len(25), nrow = 5, ncol = 5, byrow = TRUE,
                   dimnames = list(paste("Gene", seq_len(5)),
                                   paste("Cell", seq_len(5))))
testdata[3,2] <- NA

test_that("NormalizeRPM works with matrix input", {

  # All colSums are 1M
  expect_equivalent(colSums(NormalizeRPM(testdata)), rep(10^6, ncol(testdata)))

  # Class of output is a matrix
  expect_is(NormalizeRPM(testdata), "matrix")

  # Log works
  expect_equal(NormalizeRPM(testdata, log = TRUE),
               log2(1+NormalizeRPM(testdata)))

  # Log with no pseudocount throws warning
  expect_warning(NormalizeRPM(testdata, log = TRUE, pseudo.count = 0),
                 "Using 0 pseudocount: Inf may be generated.\n")

  # Scale works
  expect_equal(NormalizeRPM(testdata, scale = 100),
               NormalizeRPM(testdata)/100)
})


sce <- SingleCellExperiment::SingleCellExperiment(list(counts = testdata))

test_that("NormalizeRPM works with SingleCellExperiment input", {

  # The correct slots are created
  expect_true(all(c("cpm","normcounts") %in%
      names(SummarizedExperiment::assays(NormalizeRPM(sce = sce)))))
  expect_equal(SingleCellExperiment::cpm(NormalizeRPM(sce = sce)),
      SingleCellExperiment::normcounts(NormalizeRPM(sce = sce)))

  # All colSums are 1M
  expect_equivalent(colSums(SingleCellExperiment::cpm(NormalizeRPM(sce = sce))),
      rep(10^6, ncol(sce)))

  # Class of output is a SingleCellExperiment
  expect_is(NormalizeRPM(sce = sce), "SingleCellExperiment")

  # Log works
  expect_equal(SingleCellExperiment::logcounts(NormalizeRPM(sce = sce,
      log = TRUE)), log2(1+SingleCellExperiment::cpm(NormalizeRPM(sce = sce))))
  expect_equal(SingleCellExperiment::logcounts(NormalizeRPM(sce = sce,
      log = TRUE)), log2(1+SingleCellExperiment::cpm(NormalizeRPM(sce = sce,
      log = TRUE))))

  # Log with no pseudocount throws warning
  expect_warning(NormalizeRPM(sce = sce, log = TRUE, pseudo.count = 0),
                 "Using 0 pseudocount: Inf may be generated.\n")

  # Scale works
  expect_equal(SingleCellExperiment::normcounts(NormalizeRPM(sce = sce,
      scale = 100)),
      SingleCellExperiment::normcounts(NormalizeRPM(sce = sce))/100)
  expect_equal(SingleCellExperiment::normcounts(NormalizeRPM(sce = sce,
      scale = 100)),
      SingleCellExperiment::cpm(NormalizeRPM(sce = sce, scale = 100))/100)

})

rownames(testdata) <- sample(levels(ADImpute::transcript_length$hgnc_symbol),
                             nrow(testdata))
testtrlength <- data.frame("hgnc_symbol" = rownames(testdata),
                           "transcript_length" = rep(2,nrow(testdata)))
testoutput <- testdata[,1]/2
testoutput <- testoutput*1000000/sum(testoutput)

test_that("NormalizeTPM works with matrix input",{

  # Class of output is a matrix
  expect_is(NormalizeTPM(testdata), "matrix")

  # Log works
  expect_equal(NormalizeTPM(testdata, log = TRUE),
               log2(1+NormalizeTPM(testdata)))

  # Log with no pseudocount throws warning
  expect_warning(NormalizeTPM(testdata, log = TRUE, pseudo.count = 0),
                 "Using 0 pseudocount: Inf may be generated.\n")

  # Scale works
  expect_equal(NormalizeTPM(testdata, scale = 100),
               NormalizeTPM(testdata)/100)

  # Providing own gene length works
  expect_equal(NormalizeTPM(testdata, tr_length = testtrlength)[,1],
               testoutput)
})

sce <- SingleCellExperiment::SingleCellExperiment(list(counts = testdata))

test_that("NormalizeTPM works with SingleCellExperiment input",{

  # The correct slots are created
  expect_true(all(c("tpm","normcounts") %in%
      names(SummarizedExperiment::assays(NormalizeTPM(sce = sce)))))
  expect_equal(SingleCellExperiment::tpm(NormalizeTPM(sce = sce)),
      SingleCellExperiment::normcounts(NormalizeTPM(sce = sce)))

  # Class of output is a SingleCellExperiment
  expect_is(NormalizeTPM(sce = sce), "SingleCellExperiment")

  # Log works
  expect_equal(SingleCellExperiment::logcounts(NormalizeTPM(sce = sce,
      log = TRUE)),
      log2(1+SingleCellExperiment::normcounts(NormalizeTPM(sce = sce))))
  expect_equal(SingleCellExperiment::logcounts(NormalizeTPM(sce = sce,
      log = TRUE)), log2(1+
      SingleCellExperiment::normcounts(NormalizeTPM(sce = sce, log = TRUE))))

  # Log with no pseudocount throws warning
  expect_warning(NormalizeTPM(sce = sce, log = TRUE, pseudo.count = 0),
      "Using 0 pseudocount: Inf may be generated.\n")

  # Scale works
  expect_equal(SingleCellExperiment::normcounts(NormalizeTPM(sce = sce,
      scale = 100)),
      SingleCellExperiment::normcounts(NormalizeTPM(sce = sce))/100)
  expect_equal(SingleCellExperiment::normcounts(NormalizeTPM(sce = sce,
      scale = 100)),
      SingleCellExperiment::tpm(NormalizeTPM(sce = sce, scale = 100))/100)

  # Providing own gene length works
  expect_equal(SingleCellExperiment::normcounts(NormalizeTPM(sce = sce,
      tr_length = testtrlength))[,1], testoutput)
})
