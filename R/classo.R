#' Fit the Complex Lasso (C-Lasso)
#'
#' Fit a linear model to complex-valued data via maximum penalized likelihood
#' using the (complex) lasso penalty. The regularization path is computed
#' along a grid of values for the regularization parameter (lambda).
#' The interface is intentionally similar to that of \code{\link[glmnet]{glmnet}} in
#' the package of the same name.
#'
#' The complex lasso problem is given by
#' \deqn{\frac{1}{2n}|y - X\beta|_2^2 + \lambda |\beta|_{1}}
#' where \eqn{X}, \eqn{y}, and \eqn{\beta} all are complex. Like its real counterpart,
#' this problem can be efficiently solved via coordinate descent (implemented here).
#'
#' @param X The matrix of predictors (\eqn{X \in \C^{n \times p}}{X})
#' @param y The response vector (\eqn{y \in \C^n}{y})
#' @param nlambda The number of lambda values to use in computing the
#'    regularization path. Note that the time to run is typically sublinear
#'    in the grid size due to the use of warm starts.
#' @param lambda.min.ratio The smallest value of lambda to be used, as a fraction
#'      of the largest value of lambda used. Unlike the lasso, there is no
#'      value of lambda such that the solution is wholly sparse, but we still
#'      use lambda_max from the lasso.
#' @param lambda A user-specified sequence of lambdas to use.
#' @param standardize Should \code{X} be centered and scaled before fitting?
#' @param intercept Should the fitted model have an (unpenalized) intercept term?
#' @param thresh The convergence threshold used for the coordinate-descent
#'    algorithm used to solve the penalized regression problem. We require
#'    \eqn{\|\beta^{(k)} - \beta^{(k-1)} < \text{thresh}} before stopping.
#' @param weights Weights applied to individual
#'     observations. If not supplied, all observations will be equally
#'     weighted. Will be re-scaled to sum to \eqn{n} if
#'     necessary. (Cf. the \code{weight} argument of
#'     \code{\link[stats]{lm}})
#' @param offset A vector of length \eqn{n} included in the linear
#'     predictor.
#' @param alpha The "elastic net" mixing parameter. The penalty is defined as
#'     \deqn{\frac{1-\alpha}{2}\|\beta\|_2^2 + \alpha \|\beta\|_1}
#'     so \eqn{alpha = 1} gives the lasso while \eqn{alpha = 0} gives ridge
#'     regression. (Cf. the \code{alpha} argument of \code{\link[glmnet]{glmnet}}.)
#' @include internals.R
#' @examples
#' n <- 200
#' p <- 500
#'
#' X <- matrix(rcnorm(n * p), ncol = p)
#' beta <- rep(0, p)
#' beta[1:5] <- 1:5 * exp(im * seq(0, 2*pi, length.out = 6)[-6])
#' y <- X %*% beta
#'
#' cfit <- classo(X, y)
#' @return An object of class \code{CLassoFit} containing \itemize{
#' \item \code{coef} - A matrix of estimated coefficients
#' \item \code{intercept} - A vector of estimated intercepts if \code{intercept=TRUE}
#' \item \code{X, y} - The data used to fit the model
#' \item \code{lambda} - The vector of \eqn{\lambda}{lambda} used
#' \item \code{nnz} - The number of non-zero coefficients at each value of
#'       \eqn{\lambda}{lambda}
#' }
#' @export
classo <- function(X,
                   y,
                   weights = rep(1, nobs),
                   offset = rep(0, nobs),
                   alpha = 1,
                   nlambda=100,
                   lambda.min.ratio=ifelse(nobs < nvars, 0.01, 1e-04),
                   lambda,
                   standardize=TRUE,
                   intercept=TRUE,
                   thresh=1e-07){

    tic <- Sys.time()

    ####################
    ##
    ## Input Validation
    ##
    ####################

    nobs <- NROW(X);
    nvars <- NCOL(X);

    if(length(y) != nobs){
        stop(sQuote("NROW(X)"), " and ", sQuote("length(y)"), " must match.")
    }

    if(anyNA(X) || anyNA(y)){
        stop(sQuote("classo"), " does not support missing data.")
    }

    if(!all(is.finite(X))){
        stop("All elements of ", sQuote("X"), " must be finite.")
    }

    if(!all(is.finite(y))){
        stop("All elements of ", sQuote("y"), " must be finite.")
    }


    if(missing(weights)){
        weights <- rep(1, nobs)
    }

    if(length(weights) != nobs){
        stop(sQuote("NROW(X)"), " and ", sQuote("length(weights)"), " must match.")
    }

    if(any(weights <= 0)){
        stop("Observation weights must be strictly positive.")
    }

    if(sum(weights) != nobs){
        weights <- weights * nobs / sum(weights)
        warning(sQuote("sum(weights)"), " is not equal to ", sQuote("NROW(X)."), " Renormalizing...")
    }

    if(missing(offset)){
        offset <- rep(0, nobs)
    }

    if(length(offset) != nobs){
        stop(sQuote("NROW(X)"), " and ", sQuote("length(offset)"), " must match.")
    }

    nlambda <- as.integer(nlambda)

    if((lambda.min.ratio <= 0) || (lambda.min.ratio >= 1)){
        stop(sQuote("lambda.min.ratio"), " must be in the interval (0, 1).")
    }

    if(is.null(colnames(X))){
        colnames(X) <- paste0("V", seq_len(NCOL(X)))
    }

    if(standardize){
        Xsc <- scale(X, center=TRUE, scale=TRUE)
        X_scale <- attr(Xsc, "scaled:scale", exact=TRUE)
        X_center <- attr(Xsc, "scaled:center", exact=TRUE)

        if(!all(is.finite(Xsc))){
            stop("Non-finite ", sQuote("X"), " found after standardization.")
        }
    } else {
        Xsc <- X
        X_scale <- rep(1, nvars)
        X_center <- rep(0, nvars)
    }

    if(missing(lambda)){
        ## TODO - Check this...
        lambda_max <- max(abs(h(Xsc) %*% (y - offset - weighted.mean(y, weights) * intercept)/nobs))
        lambda <- logspace(lambda.min.ratio * lambda_max, lambda_max, length.out=nlambda)
    }

    if(length(lambda) < 1){
        stop("Must solve for at least one value of lambda.")
    }

    if(any(lambda <= 0)){
        stop("All values of ", sQuote("lambda"), " must be positive.")
    }

    if(is.unsorted(lambda)){
        warning("User-supplied ", sQuote("lambda"), " is not increasing. Sorting for maximum performance.")
        lambda <- sort(lambda)
    }

    if(thresh <= 0){
        stop(sQuote("thresh"), " must be positive.")
    }

    if(length(alpha) != 1){
        stop(sQuote("alpha"), " must be a scalar.")
    }

    if((alpha > 1) | (alpha < 0)){
        stop(sQuote("alpha"), " must be between 0 and 1.")
    }

    res <- classo_cd(X = Xsc,
                     y = y,
                     o = offset,
                     w = weights,
                     alpha = alpha,
                     lambda = lambda,
                     intercept = intercept,
                     thresh = thresh)

    ## Convert intercept to R vector (arma::cvec => R 1-column matrix)
    res$intercept <- as.vector(res$intercept)

    ## Convert coefficients and intercept back to original scale
    if(standardize){
        ## To get back to the original X, we multiply by X_scale,
        ## so we divide beta to keep things on the same unit
        res$coef <- res$coef / X_scale
        if(intercept){
            ## Map back to original X (undo scale + center)
            ##
            ## We handled the scaling above, now we adjust for the
            ## centering of X: beta(X - colMeans(X)) = beta * X - beta * colMeans(X)
            ## To uncenter we add back in beta * colMeans(X), summed over all observations
            res$intercept <- res$intercept - colSums(res$coef * X_center)
        }
    }

    if(!is.null(colnames(X))){
        rownames(res$coef) <- colnames(X)
    }

    result <- list(coef=res$coef,
                   intercept=res$intercept,
                   y=y,
                   X=X,
                   offset=offset,
                   alpha=alpha,
                   iters=res$k,
                   standardize=standardize,
                   lambda=lambda,
                   nnz=apply(res$coef, 2, function(x) sum(Mod(x) != 0)),
                   time=Sys.time() - tic)

    class(result) <- c("CLassoFit", class(result))

    result
}

has_intercept <- function(x){
    !is.null(x$intercept)
}

has_offset <- function(x){
    any(x$offset != 0)
}

#' @export
print.CLassoFit <- function(x, ..., indent=0){
    icat("Complex Lasso Fit", "\n", indent=indent)
    icat("-------------------", "\n", indent=indent)
    icat("\n", indent=indent)
    icat("N: ", NROW(x$X), ". P: ", NCOL(x$X), ".\n", sep="", indent=indent)
    icat("\n", indent=indent)
    icat("Grid:", length(x$lambda), "values of lambda. \n", indent=indent)
    icat("  Miniumum:", min(x$lambda), "\n", indent=indent)
    icat("  Maximum: ", max(x$lambda), "\n", indent=indent)
    icat("  Number of selected variables:", min(x$nnz), " --> ", max(x$nnz), "\n", indent=indent)
    icat("\n", indent=indent)
    icat("Fit Options:\n", indent=indent)
    icat("  - Intercept:              ", has_intercept(x), "\n", indent=indent)
    icat("  - Standardize X:          ", x$standardize, "\n", indent=indent)
    icat("  - L2 Mixing (Elastic Net): ", round(100 * (1- x$alpha), 3), "%\n", indent=indent)
    icat("\n", indent=indent)
    icat("Time: ", sprintf("%2.3f %s", x$time, attr(x$time, "units")), "\n", indent=indent)
    icat("\n", indent=indent)

    invisible(x)
}

# Refit exclussive lasso on new lambda grid
# Used internally by predict(exact=TRUE)
#' @importFrom utils modifyList
update_fit <- function(object, lambda, ...){
    ARGS <- list(X=object$X,
                 y=object$y,
                 standardize=object$standardize,
                 intercept=has_intercept(object),
                 lambda=lambda)

    ARGS <- modifyList(ARGS, list(...))

    do.call(classo, ARGS)
}

#' @rdname predict.CLassoFit
#' @export
#' @importFrom stats predict
coef.CLassoFit <- function(object,
                           lambda=s,
                           s=NULL,
                           exact=FALSE,
                           ...){

    predict(object, lambda=lambda, type="coefficients",
            exact=exact, ...)
}

#' Make predictions using the complex lasso.
#'
#' Make predictions using the complex lasso. Similar to \code{\link[glmnet]{predict.glmnet}}.
#' \code{coef(...)} is a wrapper around \code{predict(..., type="coefficients")}.
#'
#' @rdname predict.CLassoFit
#' @export
#' @param object An \code{CLassoFit} object produced by \code{\link{classo}}.
#' @param newx New data \eqn{X \in C^{m \times p}}{X} on which to make predictions. If not
#'    supplied, predictions are made on trainng data.
#' @param s An alternate argument that may be used to supply \code{lambda}. Included for
#'    compatability with \code{\link[glmnet]{glmnet}}.
#' @param lambda The value of the regularization paramter (\eqn{lambda}) at which to
#'    return the fitted coefficients or predicted values. If not supplied, results for
#'    the entire regularization path are returned. Can be a vector.
#' @param type The type of "prediction" to return. If \code{type="link"}, returns
#'    the linear predictor. If \code{type="response"}, returns the expected
#'    value of the response. If \code{type="coefficients"}, returns the coefficients
#'    used to calculate the linear predictor. (Cf. the \code{type} argument
#'    of \code{\link[glmnet]{predict.glmnet}})
#' @param exact Should the complex lasso be re-run for provided values of \code{lambda}?
#'    If \code{FALSE}, approximate values obtained by linear interpolation on grid points
#'    are used instead. (Cf. the \code{exact} argument of \code{\link[glmnet]{predict.glmnet}})
#' @param offset An offset term used in predictions. If not supplied, all offets are
#'    taken to be zero. If the original fit was made with an offset, \code{offset} will
#'    be required.
#' @param ... Additional arguments passed to \code{\link{classo}} if
#'    \code{exact=TRUE} and ignored otherwise.
#' @examples
#' n <- 200
#' p <- 500
#' beta <- numeric(p);
#' beta[1:10] <- 3 * exp(im * runif(10, 0, 2*pi))
#'
#' X <- matrix(rcnorm(n * p), ncol=p)
#' y <- X %*% beta + rcnorm(n)
#'
#' cfit <- classo(X, y)
#' coef(cfit, lambda=1)
#' predict(cfit, lambda=1, newx = -X)
predict.CLassoFit <- function(object,
                              newx,
                              lambda=s,
                              s=NULL,
                              type=c("link", "response", "coefficients"),
                              exact=FALSE,
                              offset,
                              ...){

    type <- match.arg(type)

    ## Get coefficients first
    if(!is.null(lambda)){
        if(exact){
            object <- update_fit(object, lambda=lambda, ...)

            if(has_intercept(object)){
                int <- matrix(object$intercept, nrow=1,
                              dimnames=list("(Intercept)", NULL))
            } else {
                int <- matrix(0,
                              nrow=1,
                              ncol=length(object$lambda),
                              dimnames=list("(Intercept)", NULL))
            }

            coef <- rbind(int, object$coef)
        } else {
            if(has_intercept(object)){
                int <- matrix(object$intercept,
                              nrow=1,
                              dimnames=list("(Intercept)", NULL))
            } else {
                int <- matrix(0,
                              nrow=1,
                              ncol=length(object$lambda),
                              dimnames=list("(Intercept)", NULL))
            }

            coef <- rbind(int, object$coef)
            lambda <- clamp(lambda, range=range(object$lambda))

            coef <- lambda_interp(coef,
                                  old_lambda=object$lambda,
                                  new_lambda=lambda)
        }
    } else {
        if(has_intercept(object)){
             int <- matrix(object$intercept,
                           nrow=1,
                           dimnames=list("(Intercept)", NULL))
        } else {
             int <- matrix(0,
                           nrow=1,
                           ncol=length(object$lambda),
                           dimnames=list("(Intercept)", NULL))
        }

        coef <- rbind(int, object$coef)
    }

    if(type == "coefficients"){
        return(coef) ## Done
    }

    if(missing(newx)){
        link <- cbind(1, object$X) %*% coef + object$offset
    } else {
        if(missing(offset)){
            if(has_offset(object)){
                stop("Original fit had an offset term but ", sQuote("offset"), " not supplied.")
            } else {
                offset <- rep(0, NROW(newx))
            }
        }
        link <- cbind(1, newx) %*% coef + offset
    }

    as.matrix(link)
}

#' @export
as.data.frame.CLassoFit <- function(x, row.names = NULL, optional = FALSE, ...){
    # FIXME - Other args...
    data.frame(coef = as.vector(coef(x)),
               variable = rownames(coef(x)),
               lambda = rep(x$lambda, each = NCOL(x$X) + 1))
}
