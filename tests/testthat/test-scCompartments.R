se <- SummarizedExperiment(rowRanges = GRanges("chr1:1-2"), )
se.noranges <- SummarizedExperiment(rowRanges = GRanges())
se.bad <- SummarizedExperiment(rowRanges = GRanges("chr1:1-10"), assays = SimpleList(a = matrix()))
se.rna <- SummarizedExperiment(rowRanges = GRanges("chr1:1-10"), assays = SimpleList(counts = matrix()))
se.array <- SummarizedExperiment(rowRanges = GRanges("chr1:1-10"), assays = SimpleList(Beta = matrix()))

test_that("scCompartments bad inputs", {
  expect_error(scCompartments(iris), err.verifySE)
  expect_error(scCompartments(se.noranges), err.verifyCoords)
  expect_error(scCompartments(se.bad, assay = "a"), "Supported assays are 'atac', and 'rna'")
  expect_error(scCompartments(se.array, assay = "atac"), err.verifyAssayNames.atac_counts)
  expect_error(scCompartments(se.rna, assay = "array"), "Supported assays are 'atac', and 'rna'")
})
