context("Test Complex Sparse Regression (classo)")

test_that("Input validation works", {
  set.seed(123)
  ## Good inputs
  n <- 200
  p <- 500
  g <- 10

  X <- matrix(rcnorm(n * p), ncol=p)
  y <- rcnorm(n)

  weights <- runif(n, 0.5, 1.5); weights <- weights / sum(weights) * n;
  offset <- runif(n, -0.5, 0.5) + im * runif(n, -0.5, 0.5);

  lambda <- seq(0.5, 5, length.out=100)
  thresh <- 1e-2  ## Low values here speed up tests

  expect_silent(classo(X, y, weights=weights, offset=offset,
                       lambda=lambda, thresh = thresh))

  ## X, y match
  expect_error(classo(X, rep(y, 2), weights=weights, offset=offset,
                      lambda=lambda, thresh = thresh))
  expect_error(classo(X, rep(y, length.out = n - 1), weights=weights, offset=offset,
                      lambda=lambda, thresh = thresh))

  ## Weights check
  expect_error(classo(X, y, weights = rep(weights, 3),
                      offset = offset, lambda = lambda, thresh = thresh))
  expect_error(classo(X, y, weights = -1 * weights,
                      offset = offset, lambda = lambda, thresh = thresh))
  expect_error(classo(X, y, weights = 0 * weights,
                      offset = offset, lambda = lambda, thresh = thresh))
  expect_error(classo(X, y, weights = im * weights,
                      offset = offset, lambda = lambda, thresh = thresh))
  expect_warning(classo(X, y, weights = 2 * weights,
                        offset = offset, lambda = lambda, thresh = thresh))

  ## Offsets check
  expect_error(exclusive_lasso(X, y, groups=groups,
                               weights=weights, offset=rep(offset, 2),
                               lambda=lambda, family="gaussian",
                               thresh=thresh, thresh_prox=thresh_prox))

  ## Convergence thresholds check
  expect_error(classo(X, y, weights=weights, offset=offset,
                      lambda=lambda, thresh=-1 * thresh))

  ## Lambda check
  expect_error(classo(X, y, weights=weights, offset=offset,
                      nlambda=-30, thresh=thresh))
  expect_error(classo(X, y, weights=weights, offset=offset,
                      nlambda=0, thresh=thresh))
  expect_error(classo(X, y, weights=weights, offset=offset,
                      lambda=-lambda, thresh=thresh))
  expect_error(classo(X, y, weights=weights, offset=offset,
                      lambda.min.ratio=2, thresh=thresh))
  expect_warning(classo(X, y, weights=weights, offset=offset,
                        lambda=rev(lambda), thresh=thresh))
})

test_that("Dynamic defaults work", {
  set.seed(100)

  n <- 200
  p <- 500
  g <- 10

  X <- matrix(rcnorm(n * p), ncol=p)
  y <- rcnorm(n)

  cfit <- classo(X, y, thresh = 1e-2)

  expect_true(all(cfit$weights == 1))
  expect_true(all(cfit$offset == 0))
  expect_equal(length(cfit$lambda),
               formals(classo)$nlambda)
  expect_equal(min(cfit$lambda)/max(cfit$lambda),
               eval(formals(classo)$lambda.min.ratio,
                    envir=data.frame(nobs=n, nvars=p)))

})

test_that("Preserves column names", {
  set.seed(500)

  n <- 200
  p <- 500
  g <- 10

  X <- matrix(rcnorm(n * p), ncol=p);
  colnames(X) <- paste0("A", 1:p)
  y <- rcnorm(n)

  cfit <- classo(X, y, thresh = 1e-2)

  expect_equal(rownames(cfit$coef),
               colnames(X))
})

test_that("Standardization works", {
  set.seed(538)

  n <- 200
  p <- 500
  g <- 10

  X <- matrix(rcnorm(n * p), ncol=p) %*% diag(seq(from=1, to=20, length.out=p))
  X[] <- scale(X, center=TRUE, scale=FALSE)

  colnames(X) <- paste0("A", 1:p)

  beta <- rep(0, p); beta[1:g] <- 2
  y <- X %*% beta + rcnorm(n)

  X_sc <- X; X_sc[] <- scale(X_sc)

  cfit <- classo(X, y, thresh=1e-8)

  cfit_sc <- classo(X_sc, y, lambda=cfit$lambda, thresh=1e-8)

  expect_equal(cfit$intercept, cfit_sc$intercept)
  expect_equal(scale(cfit$X), cfit_sc$X, check.attributes=FALSE)
  ## If Xsc gets multiplied by attr(... "scale")
  ## then cfit_sc$coef gets divided
  expect_equal(cfit$coef,
               cfit_sc$coef / attr(scale(X), "scaled:scale"))

})

test_that("Returns prox for unitary case", {
  soft_thresh <- function(x, lambda) exp(im * Arg(x)) * pmax(Mod(x) - lambda, 0)
  ## If X is unitary then
  ##
  ## argmin_{beta} \frac{1}{2n} |y - X*\beta|_2^2 + \lambda P(\beta) \\
  ## argmin_{beta} \frac{1}{2n} |X^Hy - X^H * X*\beta|_2^2 + \lambda P(\beta) \\
  ## argmin_{beta} \frac{1}{2n} |X^Hy - \beta|_2^2 + \lambda P(\beta) \\
  ## argmin_{beta} \frac{1}{2} |X^Hy - \beta|_2^2 + n * \lambda P(\beta) \\
  ## = prox_{n * lambda * P()}(X^Hy)
  ##
  set.seed(240)

  n <- 100
  p <- 100
  g <- 10

  X <- matrix(rcnorm(n * p), ncol=p); X <- qr.Q(qr(X))
  beta <- numeric(p); beta[1:g] <- 2
  y <- X %*% beta + rcnorm(n)

  lambda <- 0.75

  cfit <- classo(X, y, lambda=lambda,
                 standardize=FALSE, intercept=FALSE)
  prox_coefs <- soft_thresh(cplx_crossprod(X, y), lambda * n)

  expect_equal(as.matrix(cfit$coef), prox_coefs, check.attributes=FALSE)
})

test_that("Returns ridge when alpha = 0", {
  set.seed(5454)

  ## Low-Dimensional
  n <- 200
  p <- 20

  beta <- rep(1, p)
  X <- matrix(rcnorm(n * p), ncol=p); X[] <- scale(X)
  y <- X %*% beta + rcnorm(n)

  nlambda <- 10

  cfit <- classo(X, y, nlambda=nlambda, alpha = 0,
                 intercept=FALSE, standardize=FALSE, thresh=1e-14)

  for(i in seq_len(nlambda)){
    expect_equal(solve(cplx_crossprod(X)/n +cfit$lambda[i] * diag(1, p, p),
                       cplx_crossprod(X, y)/n),
                 as.matrix(cfit$coef[,i, drop=FALSE]),
                 check.attributes=FALSE)
  }

  ## High-Dimensional
  n <- 200
  p <- 500
  groups <- 1:p

  beta <- numeric(p); beta[1:5] <- 3

  X <- matrix(rcnorm(n * p), ncol=p); X[] <- scale(X)
  y <- X %*% beta + rcnorm(n)

  nlambda <- 10

  cfit <- classo(X, y, nlambda=nlambda, intercept=FALSE, standardize=FALSE,
                 thresh=1e-14, alpha = 0)

  for(i in seq_len(nlambda)){
    expect_equal(solve(cplx_crossprod(X)/n + cfit$lambda[i] * diag(1, p, p),
                       cplx_crossprod(X, y)/n),
                 as.matrix(cfit$coef[, i, drop=FALSE]),
                 check.attributes=FALSE)
  }
})

test_that("Matches closed form solution", {
  set.seed(3003)

  ## Low-Dimensional
  n <- 200
  p <- 20

  beta <- runif(p, 2, 3) * exp(im * runif(p, 0, 2 *pi))

  X <- matrix(rcnorm(n * p), ncol=p); X[] <- scale(X)
  y <- X %*% beta + rcnorm(n)

  nlambda <- 20

  cfit <- classo(X, y, nlambda=nlambda, intercept=FALSE, standardize=FALSE,
                 thresh=1e-14)

  for(i in seq_len(nlambda)){
    beta_hat <- cfit$coef[, i]
    ph_beta_hat <- exp(im * Arg(beta_hat))
    supp <- which(beta_hat != 0)

    if(length(supp) == 0){
      next
    }

    X_supp <- X[, supp, drop=FALSE]

    ## Apply shrinkage
    beta_analytic <- solve(cplx_crossprod(X_supp)/n,
                           cplx_crossprod(X_supp, y)/n - cfit$lambda[i] * ph_beta_hat[supp])

    expect_equal(as.vector(beta_analytic),
                 as.vector(beta_hat[supp]))
  }

  ## High-Dimensional
  n <- 100
  p <- 200
  s <- 10

  beta <- complex(p)
  beta[1:s] <- runif(s, 2, 3) * exp(im * runif(s, 0, 2 *pi))

  X <- matrix(rcnorm(n * p), ncol=p); X[] <- scale(X)
  y <- X %*% beta + rcnorm(n)

  nlambda <- 20

  cfit <- classo(X, y, nlambda=nlambda, intercept=FALSE, standardize=FALSE,
                 thresh=1e-14)

  for(i in seq_len(nlambda)){
    beta_hat <- cfit$coef[, i]
    ph_beta_hat <- exp(im * Arg(beta_hat))
    supp <- which(beta_hat != 0)

    if(length(supp) == 0){
      next
    }

    X_supp <- X[, supp, drop=FALSE]

    ## Apply shrinkage
    beta_analytic <- solve(cplx_crossprod(X_supp)/n,
                           cplx_crossprod(X_supp, y)/n - cfit$lambda[i] * ph_beta_hat[supp])

    expect_equal(as.vector(beta_analytic),
                 as.vector(beta_hat[supp]))
  }
})

test_that("Intercepts work", {
  set.seed(1945)
  n <- 100
  p <- 200

  s <- 10

  ## Basic case
  X <- matrix(rcnorm(n * p), ncol=p)
  beta <- numeric(p); beta[1:s] <- 2 * exp(im * runif(s, 0, 2 * pi))

  y <- X %*% beta + rcnorm(n)
  ym <- matrix(y, nrow=length(y), ncol=100, byrow=FALSE)

  f <- classo(X, y, standardize=TRUE, thresh=1e-12)
  expect_equal(colMeans(ym - X %*% coef(f)[-1, ]), coef(f)[1,])
  expect_equal(colMeans(predict(f)), rep(mean(y), 100))

  f <- classo(X, y, standardize=FALSE, thresh=1e-12)
  expect_equal(colMeans(ym - X %*% coef(f)[-1, ]), coef(f)[1,])
  expect_equal(colMeans(predict(f)), rep(mean(y), 100))

  ## + Weights
  w <- runif(n, 0.5, 1.5); w <- w * n / sum(w);

  f <- classo(X, y, standardize=TRUE, weights=w, thresh=1e-12)
  expect_equal(apply(ym - X %*% coef(f)[-1, ], 2, weighted.mean, w), coef(f)[1,])
  expect_equal(apply(predict(f), 2, weighted.mean, w), rep(weighted.mean(y, w), 100))

  f <- classo(X, y, standardize=FALSE, weights=w, thresh=1e-12)
  expect_equal(apply(ym - X %*% coef(f)[-1, ], 2, weighted.mean, w), coef(f)[1,])
  expect_equal(apply(predict(f), 2, weighted.mean, w), rep(weighted.mean(y, w), 100))

  ## + Offsets
  o <- runif(n, 1.5, 2.5);

  f <- classo(X, y, standardize=TRUE, offset=o, thresh=1e-12)
  expect_equal(colMeans(ym - o - X %*% coef(f)[-1, ]), coef(f)[1,])
  expect_equal(colMeans(predict(f)), rep(mean(y), 100))

  f <- classo(X, y, standardize=FALSE, offset=o, thresh=1e-12)
  expect_equal(colMeans(ym - o - X %*% coef(f)[-1, ]), coef(f)[1,])
  expect_equal(colMeans(predict(f)), rep(mean(y), 100))

  ## + Weights + Offsets
  f <- classo(X, y, standardize=TRUE, offset=o, weights=w, thresh=1e-12)
  expect_equal(apply(ym - o - X %*% coef(f)[-1, ], 2, weighted.mean, w), coef(f)[1,])
  expect_equal(apply(predict(f), 2, weighted.mean, w), rep(weighted.mean(y, w), 100))

  f <- classo(X, y, standardize=FALSE, offset=o, weights=w, thresh=1e-12)
  expect_equal(apply(ym - o - X %*% coef(f)[-1, ], 2, weighted.mean, w), coef(f)[1,])
  expect_equal(apply(predict(f), 2, weighted.mean, w), rep(weighted.mean(y, w), 100))

  ## Different distribution of X
  X <- matrix(rcnorm(n * p, mean=2 + 3i, cov = cplx_cov_spec(W = matrix(c(4, 1, 1, 2), ncol = 2))), ncol=p)
  beta <- complex(p); beta[1:s] <- 2 * exp(im * runif(s, 0, 2 * pi))
  y <- X %*% beta + rcnorm(n) + 1
  ym <- matrix(y, nrow=length(y), ncol=100, byrow=FALSE)

  f <- classo(X, y, standardize=TRUE, offset=o, weights=w, thresh=1e-12)
  expect_equal(apply(ym - o - X %*% coef(f)[-1, ], 2, weighted.mean, w), coef(f)[1,])
  expect_equal(apply(predict(f), 2, weighted.mean, w), rep(weighted.mean(y, w), 100))

  ym <- matrix(y, nrow=length(y), ncol=10, byrow=FALSE)
  f <- classo(X, y, standardize=FALSE, offset=o, weights=w, nlambda = 10)
  expect_equal(apply(ym - o - X %*% coef(f)[-1, ], 2, weighted.mean, w), coef(f)[1,])
  expect_equal(apply(predict(f), 2, weighted.mean, w), rep(weighted.mean(y, w), 10))
})

test_that("Estimation works in low-dim + low-penalty case", {
  set.seed(672)
  n <- 500
  p <- 20

  ## Basic case
  X <- matrix(rcnorm(n * p, cov = cplx_cov_spec(Sigma = 3)), ncol=p)
  beta <- rep(2, p) * exp(im * runif(p, 0, 2*pi))

  y <- X %*% beta + (3 + 0i) ## No noise

  f <- classo(X, y, thresh=1e-12, lambda=1e-14)
  expect_equal(coef(f)[1,,drop=TRUE], 3 + 0i, check.attributes=FALSE)
  expect_equal(coef(f)[-1,,drop=TRUE], beta, check.attributes=FALSE)

  ## + Weights
  w <- runif(n, 2, 3); w <- w / sum(w) * n
  f <- classo(X, y, weights=w, thresh=1e-12, lambda=1e-14)
  expect_equal(coef(f)[1,,drop=TRUE], 3 + 0i, check.attributes=FALSE)
  expect_equal(coef(f)[-1,,drop=TRUE], beta, check.attributes=FALSE)
})
