# TODO better tests/edge cases?

.z <- compartmap:::.z
.n_approx <- compartmap:::.n_approx
.p_approx <- compartmap:::.p_approx
q <- 0.95
n1 <- n2 <- 5

test_that(".z", {
  expect_equal(.z(q), 1.96, tolerance = 3e-5)
  expect_equal(.z(1), Inf)
  expect_equal(.z(0), 0)
  expect_equal(.z(-1), -Inf)
})

test_that(".n_approx", {
  expect_equal(
    .n_approx(n1, n2, 1),
    Inf
  )
  expect_equal(
    .n_approx(5, 5, q),
    10 + .z(q)^2
  )
})

test_that(".p_approx", {
  n <- .n_approx(n1, n2, q)
  expect_equal(
    .p_approx(n1, n2, q),
    0.5
  )
})

test_that("agrestiCoulCI", {
  expect_equal(
    agrestiCoullCI(n1, n2, q),
    data.frame(
      conf.est = 0.5000000,
      conf.est.lowerCI = 0.2365931,
      conf.est.upperCI = 0.7634069
    ),
    tolerance = 1e-7
  )
})
