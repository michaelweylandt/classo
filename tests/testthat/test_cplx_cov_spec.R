context("Test Covariance Spec Class")

test_that("Scalar variance supported", {
  cov_spec <- cplx_cov_spec(1)

  expect_equal(cov_spec$Sigma, matrix(1 + 0i, 1, 1))
  expect_equal(cov_spec$R, matrix(0 + 0i, 1, 1))
  expect_equal(cov_spec$W, matrix(c(0.5, 0, 0, 0.5), 2, 2))
})

test_that("No R => Proper", {
  cov_spec <- cplx_cov_spec(Sigma = matrix(c(4, 0, 0, 2), 2, 2))

  expect_true(classo:::is_proper(cov_spec))
  expect_true(all(cov_spec$R == 0))
})

test_that("Error checking Sigma works", {
  ## Error if not positive definite
  expect_error(cplx_cov_spec(-3))
  expect_error(cplx_cov_spec(Sigma = matrix(c(4, 0, 0, -1), 2, 2)))

  ## Error if not symmetric
  expect_error(cplx_cov_spec(Sigma = matrix(c(4, 2, 1, 1), 2, 2)))

  ## Can't supply other arguments
  expect_error(cplx_cov_spec(4, W = 2))
  expect_error(cplx_cov_spec(4, A = 2))
})

test_that("Round trip (Sigma <-> W) works", {

})
