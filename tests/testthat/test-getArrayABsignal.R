se.noranges <- SummarizedExperiment(rowRanges = GRanges())
se.bad <- SummarizedExperiment(rowRanges = GRanges("chr1:1-10"), assays = SimpleList(a = matrix()))

test_that("getArrayABsignal bad inputs", {
  expect_error(getArrayABsignal(iris), err.verifySE)
  expect_error(getArrayABsignal(se.noranges), err.verifyCoords)
  expect_error(getArrayABsignal(se.bad), err.verifyAssayNames.beta)
})

test_that(".preprocessArrays", {
  with_mocked_bindings(
    {
      expect_error(
        compartmap:::preprocessArrays(gr),
        "The minfi package must be installed for this functionality"
      )
    },
    requireNamespace = function(package, ..., quietly = FALSE) FALSE,
    .package = "base"
  )
})
