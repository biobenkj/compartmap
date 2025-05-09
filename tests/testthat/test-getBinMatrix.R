test_that("getBinMatrix", {
  mat <- matrix(5, nrow = 10, ncol = 10)
  mat.na <- mat
  mat.na[1, 1] <- NA

  chrs <- c("chr2", "chr2", "chr1", "chr3")
  gr <- GRanges(Rle(chrs, c(1, 1, 1, 2)), IRanges(1:5, width = 1))
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

test_that(".makeBins", {
  chrs <- c("chr2", "chr2", "chr1", "chr3")
  gr <- GRanges(Rle(chrs, c(1, 1, 1, 2)), IRanges(1:5, width = 1))
  compartmap:::.makeBins(1, 100, 10, "chr1")

  bins <- list(
    c(start = 1, end = 100, res = 10, chr = "chr1"),
    c(start = 2, end = 256, res = 6, chr = "chr2"),
    c(start = 5, end = 500, res = 2, chr = "chrX")
  )

  lapply(bins, function(input) {
    start <- as.numeric(input['start'])
    end <- as.numeric(input['end'])
    res <- as.numeric(input['res'])
    chr <- input['chr']
    ex.start <- seq(start, end, by = res)
    ex.end <- c(seq((start + res) - 1, end - 1, by = res), end - 1)
    ex.chr <- rep(chr, length(ex.start)) |> as.factor() |> unname()
    gdf <- as.data.frame(compartmap:::.makeBins(start, end, res, chr))

    expect_equal(gdf$start, ex.start)
    expect_equal(gdf$end, ex.end)
    expect_equal(gdf$seqnames, ex.chr)
  })
})
