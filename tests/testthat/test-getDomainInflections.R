test_that("getDomainInflections", {
  chrs <- c("chr2", "chr2", "chr1", "chr3")
  gr <- GRanges(
    Rle(chrs, c(1, 3, 2, 4)),
    IRanges(1:10, width=10:1)
  )
  expect_error(
    getDomainInflections(iris),
    "'gr' is not a GRanges object",
    fixed = TRUE
  )
  mcols(gr) <- data.frame(score = rep(1, 10), score2 = rep(5, 10))
  expect_error(
    getDomainInflections(gr, what = "na"),
    "'what' is not present in mcols(gr)",
    fixed = TRUE
  )
})
