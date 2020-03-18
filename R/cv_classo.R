#' CV for the Complex Lasso
#'
#' @rdname cv.classo
#' @export
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom stats sd weighted.mean
#' @param X The matrix of predictors (\eqn{X \in \C^{n \times p}}{X})
#' @param y The response vector (\eqn{y \in \C^n})
#' @param ... Additional arguments passed to \code{\link{classo}}.
#' @param type.measure The loss function to be used for cross-validation.
#' @param nfolds The number of folds (\eqn{K}) to be used for K-fold CV
  #' @param parallel Should CV run in parallel? If a parallel back-end for the
#'    \code{foreach} package is registered, it will be used. See the
#'    \code{foreach} documentation for details of different backends.
#' @param weights Weights applied to individual
#'     observations. If not supplied, all observations will be equally
#'     weighted. Will be re-scaled to sum to \eqn{n} if
#'     necessary. (Cf. the \code{weight} argument of
#'     \code{\link[stats]{lm}})
#' @param offset A vector of length \eqn{n} included in the linear
#'     predictor.
#' @examples
#' n <- 200
#' p <- 500
#' beta <- complex(p);
#' beta[1:10] <- 3 * exp((1:10) * im * 2*pi)
#'
#' X <- matrix(rcnorm(n * p), ncol=p)
#' y <- X %*% beta + rcnorm(n)
#'
#' cfit_cv <- cv.classo(X, y, nfolds=5)
#' print(cfit_cv)
#' plot(cfit_cv)
#'
#' # coef() and predict() work just like
#' # corresponding methods for classo()
#' # but can also specify lambda="lambda.min" or "lambda.1se"
#' coef(cfit_cv, lambda="lambda.1se")
cv.classo <- function(X, y, ...,
                      offset = rep(0, NROW(X)),
                      weights = rep(1, NROW(X)),
                      type.measure=c("mse", "deviance", "mae"),
                      nfolds=10, parallel=FALSE){

  tic <- Sys.time()

  type.measure <- match.arg(type.measure)

  ## The full data fit below will handle most input checking, but this needs to
  ## be done here to be reflected in the CV fits
  if(sum(weights) != NROW(X)){
    weights <- weights * NROW(X) / sum(weights)
    warning(sQuote("sum(weights)"), " is not equal to ", sQuote("NROW(X)."), " Renormalizing...")
  }

  fit <- classo(X = X, y = y, offset = offset, weights = weights, ...)

  lambda <- fit$lambda

  fold_ids <- split(sample(NROW(X)), rep(1:nfolds, length.out=NROW(X)))

  `%my.do%` <- if(parallel) `%dopar%` else `%do%`

  loss_func <- switch(type.measure,
                      mse = function(test_true, test_pred, w) weighted.mean(Mod(test_true - test_pred)^2, w),
                      mae = function(test_true, test_pred, w) weighted.mean(Mod(test_true - test_pred), w),
                      # For complex lasso (Gaussian), deviance = mse
                      deviance = function(test_true, test_pred, w) weighted.mean(Mod(test_true - test_pred)^2, w),
                      stop(sQuote(type.measure), "loss has not yet been implemented."))

  i <- NULL ## Hack to avoid global variable warning in foreach call below

  cv_err <- foreach(i=1:nfolds, .inorder=FALSE,
                    .packages=c("classo")) %my.do% {

                      X_tr <- X[-fold_ids[[i]], ];           X_te <- X[fold_ids[[i]], ]
                      y_tr <- y[-fold_ids[[i]]];             y_te <- y[fold_ids[[i]]]
                      offset_tr <- offset[-fold_ids[[i]]];   offset_te <- offset[fold_ids[[i]]]
                      weights_tr <- weights[-fold_ids[[i]]]; weights_te <- weights[fold_ids[[i]]]

                      my_fit <- classo(X = X_tr, y = y_tr, lambda=lambda,
                                       offset = offset_tr, weights = weights_tr, ...)

                      apply(predict(my_fit, newx=X_te, offset = offset_te), 2,
                            function(y_hat) loss_func(y_te, y_hat, weights_te))
                    }

  cv_err <- do.call(cbind, cv_err)

  cv_res <- apply(cv_err, 1,
                  function(x){
                    m <- mean(x)
                    se <- sd(x) / (length(x) - 1)
                    up <- m + se
                    lo <- m - se
                    c(m, se, up, lo)
                  })

  min_ix <- which.min(cv_res[1,])
  lambda.min <- lambda[min_ix]

  ## The "One-standard error" rule is defined as the largest lambda such that
  ## CV(lambda) <= CV(lambda_min) + SE(CV(lambda_min)) where lambda_min is the
  ## lambda giving minimal CV error

  lambda_min_plus_1se <- cv_res[3, min_ix]
  oneSE_ix <- max(which(cv_res[1,] <= lambda_min_plus_1se))
  lambda.1se <- lambda[oneSE_ix]

  r <- list(fit=fit,
            lambda=lambda,
            cvm=cv_res[1,],
            cvsd=cv_res[2,],
            cvup=cv_res[3,],
            cvlo=cv_res[4,],
            lambda.min=lambda.min,
            lambda.1se=lambda.1se,
            name=type.measure,
            time=Sys.time() - tic)

  class(r) <- c("CLassoFit_cv", class(r))

  r
}

#' @export
print.CLassoFit_cv <- function(x, ...){
  icat("Complex Lasso CV", "\n")
  icat("------------------", "\n")
  icat("\n")
  icat("Lambda (Min Rule):", x$lambda.min, "\n")
  icat("Lambda (1SE Rule):", x$lambda.1se, "\n")
  icat("\n")
  icat("Loss Function:", x$name, "\n")
  icat("\n")
  icat("Time: ", sprintf("%2.3f %s", x$time, attr(x$time, "units")), "\n")
  icat("\n")
  icat("Full Data Fit", "\n")
  icat("-------------", "\n")
  print(x$fit, indent=2)

  invisible(x)
}

#' @export
#' @importFrom stats predict
predict.CLassoFit_cv <- function(object, ...){
  dots <- list(...)
  if("s" %in% names(dots)){
    s <- dots$s
    if(s == "lambda.min"){
      s <- object$lambda.min
    }
    if(s == "lambda.1se"){
      s <- object$lambda.1se
    }
    dots$s <- s
  }
  if("lambda" %in% names(dots)){
    lambda <- dots$lambda
    if(lambda == "lambda.min"){
      lambda <- object$lambda.min
    }
    if(lambda == "lambda.1se"){
      lambda  <- object$lambda.1se
    }
    dots$lambda  <- lambda
  }

  do.call(predict, c(list(object$fit), dots))
}


#' @export
#' @importFrom stats coef
coef.CLassoFit_cv <- function(object, ...){
  dots <- list(...)
  if("s" %in% names(dots)){
    s <- dots$s
    if(s == "lambda.min"){
      s <- object$lambda.min
    }
    if(s == "lambda.1se"){
      s <- object$lambda.1se
    }
    dots$s <- s
  }
  if("lambda" %in% names(dots)){
    lambda <- dots$lambda
    if(lambda == "lambda.min"){
      lambda <- object$lambda.min
    }
    if(lambda == "lambda.1se"){
      lambda  <- object$lambda.1se
    }
    dots$lambda  <- lambda
  }

  do.call(coef, c(list(object$fit), dots))
}
