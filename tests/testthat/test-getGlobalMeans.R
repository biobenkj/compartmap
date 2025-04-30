test_that("getGlobalMeans", {
  nrows <- 100
  ncols <- 6
  mat.counts <- matrix(runif(nrows * ncols, 1, 10), nrows)
  mat.beta <- matrix(runif(nrows * ncols, 0, 1), nrows)
  gr <- GRanges(
    Rle(c("chr1", "chr2", "chr3", "chr4"), c(25, 25, 25, 25)),
    IRanges(1:100, width = 1:1)
  )
  mat.rownames <- as.character(gr)

  se.rna <- SummarizedExperiment(
    assays = SimpleList(counts = mat.counts),
    rowRanges = gr
  )
  expected.rna <- matrix(rowMeans(mat.counts), ncol = 1)
  colnames(expected.rna) <- "globalMean"
  rownames(expected.rna) <- mat.rownames

  expected.rna.two_col <- matrix(rowMeans(mat.counts[, 1:2]))
  colnames(expected.rna.two_col) <- "globalMean"
  rownames(expected.rna.two_col) <- mat.rownames

  se.array <- SummarizedExperiment(
    assays = SimpleList(Beta = mat.beta),
    rowRanges = gr
  )
  expected.array <- matrix(rowMeans(flogit(mat.beta)), ncol = 1)
  colnames(expected.array) <- "globalMean"
  rownames(expected.array) <- mat.rownames

  expect_equal(
    getGlobalMeans(se.rna, assay = "rna"),
    expected.rna
  )
  expect_equal(
    getGlobalMeans(se.rna, assay = "rna"),
    expected.rna
  )
  expect_equal(
    getGlobalMeans(se.rna, assay = "atac"),
    expected.rna
  )
  expect_equal(
    getGlobalMeans(se.array, assay = "array"),
    expected.array
  )
  expect_equal(
    getGlobalMeans(se.rna, assay = "rna", targets = c(1, 2)),
    expected.rna.two_col
  )

  se.invalid <- SummarizedExperiment(
    assays = SimpleList(badname = mat.counts),
    rowRanges = gr
  )
  expect_error(
    getGlobalMeans(se.invalid, assay = "rna")
  )
})
