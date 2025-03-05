test_that("getShrinkageTargets", {
  nrows <- 10
  ncols <- 25
  mat.no_colnames <- matrix(runif(nrows * ncols, 1, 10), nrows)
  mat <- mat.no_colnames
  colnames(mat) <- letters[-26]

  se.no_colnames <- SummarizedExperiment(assays=SimpleList(counts=mat.no_colnames))
  se.no_colnames.full <- se.no_colnames
  se.no_colnames.two_cols <- se.no_colnames[, 1:2]

  se <- SummarizedExperiment(assays=SimpleList(counts=mat))
  se.full <- se
  se.two_cols <- se[, 1:2]

  # subset without colnames
  expect_equal(
    getShrinkageTargets(se.no_colnames, 1:25),
    se.no_colnames.full
  )
  expect_equal(
    getShrinkageTargets(se.no_colnames, 1:2),
    se.no_colnames.two_cols
  )

  # subset with colnames and indices
  expect_equal(
    getShrinkageTargets(se, letters[-26]),
    se.full
  )
  expect_equal(
    getShrinkageTargets(se, c(1:25)),
    se.full
  )
  expect_equal(
    getShrinkageTargets(se, c(1, 2)),
    se.two_cols
  )

  # errors
  expect_error(
    getShrinkageTargets(se, c(letters[1:26]))
  )
  expect_error(
    getShrinkageTargets(se.no_colnames, c(1:26))
  )
  expect_error(
    getShrinkageTargets(se.no_colnames, c(letters[-26]))
  )
})

