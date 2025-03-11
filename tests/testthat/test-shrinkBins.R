test_that("shrinkBins", {
  se <- SummarizedExperiment()
  expect_error(
    shrinkBins(iris, se),
    "Input needs to be a SummarizedExperiment",
    fixed = TRUE
  )
})
