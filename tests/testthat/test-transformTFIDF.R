# transformTFIDF {{{
test_that("transformTFIDF", {

  expect_error(
    transformTFIDF(iris),
    "Input needs to be a matrix",
    fixed = TRUE
  )

})
# }}}

# .binarizeMatrix {{{
test_that(".binarizeMatrix", {
  mat <- c(-20, -5, -1, -0.5, 0, 1, 0.25, 2, 10, 15) |> matrix(ncol = 5)
  mat.binary <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1) |> matrix(ncol = 5)
  expect_equal(
    compartmap:::.binarizeMatrix(mat),
    mat.binary
  )
  # TODO: handling NA, Inf
})
# }}}

# .tfidf {{{
test_that(".tfidf", {
  i <- c(1,3:8)
  j <- c(2,9,6:10)
  x <- 1:7
  mat <- sparseMatrix(i, j, x = x)                    ##  8 x 10 "dgCMatrix"
  expect_equal(
    compartmap:::.tfidf(mat, 1:8),
    mat * 1:8
  )
})
# }}}



