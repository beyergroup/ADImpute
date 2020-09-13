test_that("WriteTXT works", {

  file <- tempfile()
  expect_false(file.exists(file))

  WriteTXT(iris, file)

  # Confirm file is written
  expect_true(file.exists(file))

  # Confirm file contents match expected
  expect_equivalent(iris, read.table(file, stringsAsFactors = TRUE))

  unlink(file)
})

test_that("WriteCSV works", {

  file <- tempfile()
  expect_false(file.exists(file))

  WriteCSV(iris, file)

  # Confirm file is written
  expect_true(file.exists(file))

  # Confirm file contents match expected
  expect_equivalent(iris, read.csv(file, stringsAsFactors = TRUE,
                                   row.names = 1))

  unlink(file)
})
