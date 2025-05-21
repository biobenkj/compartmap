ranges <- c("chr1:1-10", "chr1:10-20", "chr1:21-30", "chr1:31-40")
gr <- GRanges(ranges)
cols <- list(
  'a' = c(1, 2, 3, 4),
  'b' = c(1, 2, 3, 4),
  'c' = c(1, 2, 3, 4),
  'd' = c(1, 2, 3, 4)
)
mcols(gr) <- cols

mat <- matrix(1:4, nrow = 4, ncol = 4)
rownames(mat) <- ranges
colnames(mat) <- letters[1:4]

test_that("granges <-> matrix", {
  # setAs
  expect_equal(as(gr, "matrix"), mat)
  expect_equal(as(mat, "GRanges"), gr)

  # setMethod
  expect_equal(granges(mat), gr)
  expect_equal(as.matrix(gr), mat)
})

test_that("granges <-> matrix methods", {
  expect_true(existsMethod("granges", "matrix"))
  expect_true(existsMethod("as.matrix", "GRanges"))
})
