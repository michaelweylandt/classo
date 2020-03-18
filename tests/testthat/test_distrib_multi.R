context("Test Mutlivariate Distributions")

test_that("Density errors on NA", {
  expect_error(dmvcnorm(NA))
  expect_error(dmvcnorm(0, mean = NA))
})

test_that("Density constant under affine shifts", {
  A <- c(3 + 4i, 2 + 2i)
  B <- c(5 + 12i, 6 + 13i)
  C <- c(0 + 1i, 1 + 0i)

  expect_equal(dmvcnorm(A, mean = B), dmvcnorm(A + C, mean = B + C))
  expect_equal(dmvcnorm(A, mean = A), dmvcnorm(A + C, mean = A + C))
  expect_equal(dmvcnorm(0, mean = 1+1i), dmvcnorm(1 + 1i, mean = 0))
  expect_equal(dmvcnorm(c(0, 0), mean = c(1 + 1i, 1 + 1i)),
               dmvcnorm(c(1 + 1i, 1 + 1i), mean = c(0, 0)))

  expect_equal(dmvcnorm(0, mean = 0+1i), dmvcnorm(0, mean = 1 + 0i))
})

test_that("Density matches cmvnorm for proper case", {
  skip_if_not_installed("cmvnorm")

  ## I can't figure out what cmvnorm::dcmvnorm is doing here
  ## The formula in the documentation appears right, but
  ## it doesn't match what comes out

  Z   <- c(0 + 0i, 1 + 1i, 2 + 2i)
  mu  <- c(0.5 + 0i, 1 + 1.5i, 2 + 1.5i)
  cov <- cplx_cov_spec(Sigma = diag(1:3))
  Sigma <- cov$Sigma
  Sigma_inv <- solve(Sigma)

  # Denominator - det(pi * Gamma) - must be manually calculated
  den <- prod(eigen(pi * Sigma)$values)

  expect_equal(dmvcnorm(Z, mu, cov),
               as.double(exp(-Conj(Z - mu) %*% Sigma_inv %*% (Z - mu)) / den))

  skip("cmvnorm::dcmvnorm seems off")
  expect_equal(dmvcnorm(Z, mu, cov),
               cmvnorm::dcmvnorm(Z, mu, Sigma))
})

test_that("Special cases of density work", {
  expect_equal(1 / pi, dmvcnorm(0))
  expect_equal(1 / pi^2, dmvcnorm(c(0, 0)))
  expect_equal(1 / pi^3, dmvcnorm(c(0, 0, 0)))
  expect_equal(1 / pi, dmvcnorm(1 + 1i, mean = 1 + 1i))
  expect_equal(1 / pi^2, dmvcnorm(rep(1 + 1i, 2), mean = rep(1 + 1i, 2)))
})

test_that("Recycling works - density", {
  set.seed(125)
  X <- rmvcnorm(10, mean = c(0 + 2i, 4 + 1i))

  short_mean <- cbind(1:5, 0 + 1i)
  long_mean  <- cbind(c(1:5, 1:5), 0 + 1i)

  expect_equal(dmvcnorm(X, short_mean),
               dmvcnorm(X, long_mean))

  short_X <- rmvcnorm(10, cov = cplx_cov_spec(Sigma = 3))
  long_X  <- rbind(short_X, short_X)

  mean <- matrix(1:20 + 3i, ncol = 1)

  expect_equal(dmvcnorm(short_X, mean),
               dmvcnorm(long_X, mean))
})

test_that("Recycling works - RNG", {
  set.seed(150)
  short_X <- rmvcnorm(20, cbind(1:10, 2i))

  set.seed(150)
  long_X  <- rmvcnorm(20, rbind(cbind(1:10, 2i),
                              cbind(1:10, 2i)))

  expect_equal(short_X, long_X)
})
