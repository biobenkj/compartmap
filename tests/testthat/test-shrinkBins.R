test_that("shrinkBins", {
  se <- SummarizedExperiment()
  expect_error(
    shrinkBins(iris, se),
    "Input needs to be a SummarizedExperiment",
    fixed = TRUE
  )
})

test_that("shinkBins atac_fun", {
  vec <- runif(10, 1, 10)
  expect_equal(
    compartmap:::atac_fun(vec),
    sqrt(mean(vec)) * length(vec)
  )
})

test_that("shinkBins .ebayes", {
  vec <- 1:10
  targets <- c(4, 5)
  prior = 5
  expect_equal(
    compartmap:::.ebayes(vec, prior),
    prior + sd(vec) * (vec - prior)
  )
  expect_equal(
    compartmap:::.ebayes(vec, prior, targets),
    prior + sd(vec[targets]) * (vec - prior)
  )
})

test_that("shinkBins .jse", {
  vec <- 1:10
  prior = 5
  targets <- c(4, 5)

  m <- length(vec)
  yv.norm <- sum((vec - prior)^2)

  expect_equal(
    compartmap:::.jse(vec, prior),
    {
      sdsq <- sd(vec)^2
      c <- 1 - (((m - 3) * sdsq) / yv.norm)
      prior + c * (vec - prior)
    }
  )
  expect_equal(
    compartmap:::.jse(vec, prior, targets),
    {
      sdsq_targeted <- sd(vec[targets])^2
      c <- 1 - (((m - 3) * sdsq_targeted) / yv.norm)
      prior + c * (vec - prior)
    }
  )
})

test_that("shinkBins .shrinkRow", {
  matRow <- c((1:100), 50)
  names(matRow) <- c(paste(1:100), "globalMean")
  mean.val <- matRow[101]
  targets = c(10:90)
  expect_equal(
    compartmap:::.shrinkRow(matRow, jse = T, targets = NULL),
    compartmap:::.jse(matRow[1:100], mean.val)
  )
  expect_equal(
    compartmap:::.shrinkRow(matRow, jse = F, targets = NULL),
    compartmap:::.ebayes(matRow[1:100], mean.val)
  )
  expect_equal(
    compartmap:::.shrinkRow(matRow, jse = T, targets = 10:90),
    compartmap:::.jse(matRow[1:100], mean.val, targets)
  )
  expect_equal(
    compartmap:::.shrinkRow(matRow, jse = F, targets = 10:90),
    compartmap:::.ebayes(matRow[1:100], mean.val, targets)
  )
})
