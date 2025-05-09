test_that("meanSmoother", {
  mat <- matrix(1:5, nrow = 5, ncol = 5)
  expect_message(
    meanSmoother(mat, k = 0),
    "Returning unsmoothed x. This is probably an error.",
    fixed = TRUE
  )
  expect_equal(
    meanSmoother(mat, k = 0),
    mat
  )

  # TODO get better error
  expect_error(
    meanSmoother(mat[1:2], k = 3),
    "length(x) >= k is not TRUE",
    fixed = TRUE
  )
})

test_that(".meanSmoother.internal", {
})
