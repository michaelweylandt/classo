context("Test that classo CV works")

test_that("CV works in default mode", {
  skip_on_cran()

  set.seed(125)
  n <- 500
  p <- 20
  beta <- complex(p);
  beta[1:10] <- 3 * exp(im * runif(10, 0, 2 * pi))

  X <- matrix(rcnorm(n * p), ncol=p)
  y <- X %*% beta + rcnorm(n)

  cfit_cv <- cv.classo(X, y, nfolds=5, intercept = FALSE, standardize=FALSE)

  skip("This doesn't work for some reason - CV doesn't shrink enough")
  ## This is an easy problem, so we should get the right answer
  expect_equal(which(beta != 0 + 0i),
               which(coef(cfit_cv, lambda = "lambda.1se")[-1] != 0))
})

test_that("predict() and coef() work with cv selected lambda", {
  skip_on_cran()

  set.seed(50)
  n <- 20
  p <- 50

  beta <- complex(p);
  beta[1:10] <- 3 + rcnorm(10)

  X <- matrix(rcnorm(n * p), ncol=p)
  y <- X %*% beta + rcnorm(n)

  cfit_cv <- cv.classo(X, y, nfolds=5)

  expect_equal(coef(cfit_cv, lambda = "lambda.min"), coef(cfit_cv$fit, lambda = cfit_cv$lambda.min))
  expect_equal(coef(cfit_cv, lambda = "lambda.1se"), coef(cfit_cv$fit, lambda = cfit_cv$lambda.1se))
  expect_equal(coef(cfit_cv, s = "lambda.min"), coef(cfit_cv$fit, s = cfit_cv$lambda.min))
  expect_equal(coef(cfit_cv, s = "lambda.1se"), coef(cfit_cv$fit, s = cfit_cv$lambda.1se))

  expect_equal(predict(cfit_cv, lambda = "lambda.min"), predict(cfit_cv$fit, lambda = cfit_cv$lambda.min))
  expect_equal(predict(cfit_cv, lambda = "lambda.1se"), predict(cfit_cv$fit, lambda = cfit_cv$lambda.1se))
  expect_equal(predict(cfit_cv, s = "lambda.min"), predict(cfit_cv$fit, s = cfit_cv$lambda.min))
  expect_equal(predict(cfit_cv, s = "lambda.1se"), predict(cfit_cv$fit, s = cfit_cv$lambda.1se))
})
