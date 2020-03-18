context("Test Univariate Distributions")

test_that("Density handles NA correctly", {
  expect_true(is.na(dcnorm(NA)))
  expect_false(is.na(dcnorm(0)))
  expect_true(is.na(dcnorm(c(0, NA)))[2])
})

test_that("Density constant under affine shifts", {
  expect_equal(dcnorm(0, mean = 0), dcnorm(1, mean = 1))
  expect_equal(dcnorm(0, mean = 0), dcnorm(1 + 1i, mean = 1 + 1i))
  expect_equal(dcnorm(0, mean = 1+1i), dcnorm(1 + 1i, mean = 0))
  expect_equal(dcnorm(0, mean = 0+1i), dcnorm(0, mean = 1 + 0i))
})

test_that("Density matches cmvnorm for proper case", {
  skip_if_not_installed("cmvnorm")

  skip("cmvnorm appears wrong")
  # I can't figure out what exactly cmvnorm::dcmvnorm is doing, but
  # dcmvnorm(0, mean = 0, sigma = 1) should be 1/pi, not Inf

  set.seed(125)

  X <- rcnorm(250)
  f <- Vectorize(cmvnorm::dcmvnorm, "z", SIMPLIFY = TRUE)

  expect_equal(dcnorm(X), f(X))
})

test_that("Special cases of density work", {
  expect_equal(1 / pi, dcnorm(0))
  expect_equal(1 / pi, dcnorm(1 + 1i, mean = 1 + 1i))
})

test_that("Recycling works - density", {
  set.seed(125)
  X <- rcnorm(10)

  short_mean <- 1:5
  long_mean  <- c(1:5, 1:5)

  expect_equal(dcnorm(X, short_mean),
               dcnorm(X, long_mean))

  short_X <- rcnorm(10)
  long_X  <- c(short_X, short_X)

  mean <- 1:20

  expect_equal(dcnorm(short_X, mean),
               dcnorm(long_X, mean))
})

test_that("Recycling works - RNG", {
  set.seed(150)
  short_X <- rcnorm(20, 1:10)

  set.seed(150)
  long_X  <- rcnorm(20, c(1:10, 1:10))

  expect_equal(short_X, long_X)
})
