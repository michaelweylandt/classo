context("Test that helper functions for CLassoFit class work")

test_that("coef() works with original lambda", {
  set.seed(5454)

  n <- 200
  p <- 500

  beta <- complex(p); beta[1:10] <- 2 * rcnorm(10, mean = 1 + 1i)
  X <- matrix(rcnorm(n * p), ncol=p)
  y <- X %*% beta + rcnorm(n)

  nlambda <- 10

  ## With intercept
  cfit <- classo(X, y, nlambda=nlambda, intercept=TRUE, standardize=FALSE,
                 thresh=1e-14)

  c_cfit <- coef(cfit)
  expect_equal(NROW(c_cfit), p + 1)
  expect_equal(NCOL(c_cfit), nlambda)

  expect_equal(c_cfit[1,,drop=TRUE], cfit$intercept)
  expect_equal(c_cfit[-1,,drop=FALSE], cfit$coef,
               check.attributes=FALSE)

  ## Repeat without intercept
  cfit <- classo(X, y, nlambda=nlambda, intercept=FALSE, standardize=FALSE,
                 thresh=1e-14)

  c_cfit <- coef(cfit)
  expect_equal(NROW(c_cfit), p + 1)
  expect_equal(NCOL(c_cfit), nlambda)

  expect_equal(c_cfit[1,,drop=TRUE], rep(0 + 0i, nlambda))
  expect_equal(c_cfit[-1,,drop=FALSE], cfit$coef,
               check.attributes=FALSE)
})

test_that("coef() works with new lambda (exact=FALSE)", {
  set.seed(4114)

  n <- 200
  p <- 500

  beta <- complex(p); beta[1:10] <- 2
  X <- matrix(rcnorm(n * p), ncol=p)
  y <- X %*% beta + rcnorm(n)

  nlambda <- 30

  ## With intercept
  cfit <- classo(X, y, nlambda=nlambda, thresh = 1e-14,
                  intercept=TRUE, standardize=FALSE)

  ### Truly new value
  c_cfit <- coef(cfit, lambda=mean(cfit$lambda[23:24]))

  expect_true(all(abs(cfit$coef[, 23]) >= abs(c_cfit[-1])))
  expect_true(all(abs(cfit$coef[, 24]) <= abs(c_cfit[-1])))

  ### Old min endpoint
  c_cfit <- coef(cfit, lambda=min(cfit$lambda))
  expect_equal(c_cfit[-1], cfit$coef[,1], check.attributes=FALSE)

  ### Old max endpoint
  c_cfit <- coef(cfit, lambda=max(cfit$lambda))
  expect_equal(c_cfit[-1], cfit$coef[,nlambda], check.attributes=FALSE)

  ### Exact match to old internal value
  c_cfit <- coef(cfit, lambda=cfit$lambda[15])
  expect_equal(c_cfit[-1], cfit$coef[,15], check.attributes=FALSE)
})

test_that("coef() works with new lambda (exact=TRUE)", {
  set.seed(6772)

  n <- 200
  p <- 500

  beta <- complex(p); beta[1:10] <- 2
  X <- matrix(rcnorm(n * p), ncol=p)
  y <- X %*% beta + rcnorm(n)

  nlambda <- 30

  ## Without intercept
  cfit <- classo(X, y, nlambda=nlambda, thresh=1e-14, intercept=FALSE, standardize=FALSE)

  ### Truly new value
  c_cfit <- coef(cfit, lambda=mean(cfit$lambda[23:24]),
                 exact=TRUE)

  expect_true(all(abs(cfit$coef[, 23]) >= abs(c_cfit[-1])))
  expect_true(all(abs(cfit$coef[, 24]) <= abs(c_cfit[-1])))

  ### Need to refit with high precision to match original fit

  ### Old min endpoint
  c_cfit <- coef(cfit, lambda=min(cfit$lambda),
                 exact=TRUE, thresh=1e-14)
  expect_equal(c_cfit[-1], cfit$coef[,1], check.attributes=FALSE)

  ### Old max endpoint
  c_cfit <- coef(cfit, lambda=max(cfit$lambda),
                 exact=TRUE, thresh=1e-14)
  expect_equal(c_cfit[-1], cfit$coef[,nlambda], check.attributes=FALSE)

  ### Exact match to old internal value
  c_cfit <- coef(cfit, lambda=cfit$lambda[15],
                 exact=TRUE, thresh=1e-14)
  expect_equal(c_cfit[-1], cfit$coef[,15], check.attributes=FALSE)

})

test_that("predict() works -- training data", {
  set.seed(8119)

  n <- 200
  p <- 500

  beta <- complex(p); beta[1:10] <- rcnorm(10, mean = 3 + 2i)
  X <- matrix(rcnorm(n * p), ncol=p)
  y <- X %*% beta + rcnorm(n)

  nlambda <- 30

  ## With intercept
  cfit <- classo(X, y, nlambda=nlambda, thresh=1e-14, intercept=TRUE, standardize=TRUE)

  expect_equal(predict(cfit, lambda=1),
               cbind(1, X) %*% as.matrix(coef(cfit, lambda=1)),
               check.attributes=FALSE)

  expect_equal(predict(cfit),
               cbind(1, X) %*% as.matrix(coef(cfit)),
               check.attributes=FALSE)

  expect_equal(predict(cfit, type="response"), ## Doesn't change for Gaussian
               cbind(1, X) %*% as.matrix(coef(cfit)),
               check.attributes=FALSE)

  ## Without intercept + with offset
  o <- runif(n, 0.5, 1.5) + runif(n, 0.5, 1.5) * im
  cfit <- classo(X, y, nlambda=nlambda, intercept=FALSE, offset=o, thresh=1e-14)

  expect_equal(predict(cfit, lambda=1),
               X %*% as.matrix(coef(cfit, lambda=1))[-1,, drop=FALSE] + o,
               check.attributes=FALSE)

  expect_equal(predict(cfit),
               X %*% as.matrix(coef(cfit))[-1,, drop=FALSE] + o,
               check.attributes=FALSE)

  expect_equal(predict(cfit, type="response"), ## Doesn't change for Gaussian
               X %*% as.matrix(coef(cfit))[-1,, drop=FALSE] + o,
               check.attributes=FALSE)
})

test_that("predict() works -- test data", {
  set.seed(8119)

  n <- 200
  p <- 500

  beta <- complex(p); beta[1:10] <- rcnorm(10, mean = 3 + 2i)
  X <- matrix(rcnorm(n * p), ncol=p); X2 <- matrix(rcnorm(n * p), ncol=p)
  y <- X %*% beta + rcnorm(n)

  nlambda <- 30

  ## With intercept
  cfit <- classo(X, y, nlambda=nlambda, intercept=TRUE, standardize=FALSE, thresh=1e-14)

  expect_equal(predict(cfit, newx=X2, lambda=1),
               cbind(1, X2) %*% as.matrix(coef(cfit, lambda=1)),
               check.attributes=FALSE)

  expect_equal(predict(cfit, newx=X2),
               cbind(1, X2) %*% as.matrix(coef(cfit)),
               check.attributes=FALSE)

  expect_equal(predict(cfit, newx=X2, type="response"), ## Doesn't change for Gaussian
               cbind(1, X2) %*% as.matrix(coef(cfit)),
               check.attributes=FALSE)

  ## Without intercept + with offset
  o <- runif(n, 0.5, 1.5)
  cfit <- classo(X, y, nlambda=nlambda, intercept=FALSE, offset=o, thresh=1e-14)

  expect_equal(predict(cfit, lambda=1, newx=X2, offset=o),
               X2 %*% as.matrix(coef(cfit, lambda=1))[-1,, drop=FALSE] + o,
               check.attributes=FALSE)

  expect_equal(predict(cfit, newx=X2, offset=o),
               X2 %*% as.matrix(coef(cfit))[-1,, drop=FALSE] + o,
               check.attributes=FALSE)

  expect_equal(predict(cfit, newx=X2, offset=o, type="response"), ## Doesn't change for Gaussian
               X2 %*% as.matrix(coef(cfit))[-1,, drop=FALSE] + o,
               check.attributes=FALSE)

  expect_error(predict(cfit, newx=X2)) ## Offset required
})
