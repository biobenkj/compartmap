se <- SummarizedExperiment(rowRanges = GRanges("chr1:1-2"), )
se.noranges <- SummarizedExperiment(rowRanges = GRanges())
se.rna <- SummarizedExperiment(rowRanges = GRanges("chr1:1-10"), assays = SimpleList(counts = matrix()))
se.array <- SummarizedExperiment(rowRanges = GRanges("chr1:1-10"), assays = SimpleList(Beta = matrix()))
se.bad <- SummarizedExperiment(rowRanges = GRanges("chr1:1-10"), assays = SimpleList(a = matrix()))

test_that("getDenoisedCorMatrix bad inputs", {
  expect_error(getDenoisedCorMatrix(iris), err.verifySE)
  expect_error(getDenoisedCorMatrix(se.bad), err.verifyAssayNames.rna_counts)
  expect_error(getDenoisedCorMatrix(se.rna, assay = "array"), err.verifyAssayNames.beta)
  expect_error(getDenoisedCorMatrix(se.array, assay = "rna"), err.verifyAssayNames.rna_counts)
})
