test_that("getBinMatrix", {
  mat <- matrix(5, nrow = 10, ncol = 10)
  mat.na <- mat
  mat.na[1, 1] <- NA

  chrs <- c("chr2", "chr2", "chr1", "chr3")
  gr <- GRanges(Rle(chrs, c(1, 1, 1, 2)), IRanges(1:5, width=1))
  gr.empty <- GRanges()

  expect_error(
    getBinMatrix(mat.na, gr),
    "Matrix must not contain NAs",
    fixed = TRUE
  )

  expect_error(
    getBinMatrix(mat, gr.empty),
    "Provided GRanges must have length equal to the matrix number of rows",
    fixed = TRUE
  )
})
