#' Estimate Complex Graphical Models
#'
#' @export
#' @param X An n-by-p complex data matrix.
#' @param method Which estimator should be used: the Meinshausen-Buhlmann pseudo-likelihood
#'               (neighborhood selection) procecdure (\code{mb}) or the full penalized
#'               likelihood (\code{glasso})?
#' @param proper Should the data be assumed to come from a proper (circularly symmetric)
#'               Gaussian? Currently, improper graphical models are not supported.
#' @param rule Which symmetrization rule should be used in the neighborhood selection
#'             procedure? The \code{or} rule assumes an edge between X_i and X_j
#'             if the regression of X_i on X_j or the regression of X_j on X_i give
#'             non-zero coefficients. The \code{and} rule is more conservative
#'             (sparser) and only includes an edge if both regression coefficients
#'             are non-zero.
#' @param nlambda How many values of \eqn{\lambda} to use.
#' @param lambda.min.ratio The smallest value of \eqn{\lambda} to use compared to
#'                         the theoretical maximum.
#' @param lambda A user-provided set of regularization parameters.
#' @param parallel Should CV run in parallel? If a parallel back-end for the
#'    \code{foreach} package is registered, it will be used. See the
#'    \code{foreach} documentation for details of different backends.
cglasso <- function(X, method = c("mb", "glasso"),
                    proper = TRUE, rule = c("or", "and"),
                    nlambda=100, lambda.min.ratio = 0.001, lambda, parallel = TRUE){
  tic <- Sys.time()

  method <- match.arg(method)
  rule    <- match.arg(rule)

  if(!isTRUE(proper)){
    stop("Improper Gaussian graphical models not yet supported.")
  }

  if(anyNA(X)){
    stop(dQuote("classo"), " does not currently support missing values in ", dQuote("X."))
  }

  if(is_hermitian(X)){
    warning(dQuote("X"), " looks like a covariance matrix, not a data matrix.")
  }

  if(!missing(lambda)){
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
  }

  if(is.null(colnames(X))){
    colnames(X) <- paste0("V", seq_len(NCOL(X)))
  }

  if(method == "mb"){
    fit <- cglasso_mb(X, rule = rule, nlambda = nlambda,
                      lambda.min.ratio = lambda.min.ratio, lambda = lambda,
                      parallel = parallel)
  } else {
    fit <- cglasso_cglasso(X, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, lambda = lambda)
  }

  fit$X      <- X
  fit$proper <- proper
  fit$rule    <- rule
  fit$method <- method
  fit$nnz    <- apply(fit$edges, 3, function(X) sum(X != 0) / 2) # Divide by two to avoid double counting undirected edges
  fit$time <- Sys.time() - tic

  class(fit) <- c("CGLassoFit", class(fit))

  fit
}

#' @importFrom stats weighted.mean
#' @noRd
cglasso_mb <- function(X, rule, nlambda, lambda.min.ratio, lambda, parallel){
  ## We standardize X here and figure out a global lambda grid
  ## No need to unscale, so can be a bit sloppy
  X <- scale(X)

  nobs <- NROW(X)
  nvar <- NCOL(X)

  if(missing(lambda)){
    # For each pseudo-likelihood, the lambda max is X_{-j}^H %*% X_j / nobs
    # (Compare this to the general regression case: we omit offsets, weights and intercepts.)
    # Taking the max over all j, we get our "global" lambda_max. A cheap approximation
    # is to use X in lieu of X_{-j} and then the calculation reduces to X^HX / n
    # but we have to kill the diagonal or we have too high a penalty
    #
    XHX <- cplx_crossprod(X); diag(XHX) <- 0
    lambda_max <- max(Mod(XHX)/nobs)
    lambda <- logspace(lambda.min.ratio * lambda_max, lambda_max, length.out=nlambda)
  }

  ## Now do sparse regression at this lambda grid over all observations. We omit offsets
  ## and weights to keep things simple and intercepts because we are assuming a mean zero Gaussian
  ## (especialy after standardizing X above). We keep only the sparsity pattern of each regression

  `%my.do%` <- if(parallel) `%dopar%` else `%do%`

  j <- NULL ## Hack to avoid global variable warning in foreach call below

  edge_matrices <- foreach(j=seq_len(nvar), .inorder=TRUE,
                           .packages=c("classo")) %my.do% {

      cfit <- classo(X[,-j, drop = FALSE], X[, j, drop=TRUE], lambda = lambda,
                     standardize = FALSE, intercept = FALSE)

      # Matrix of edge indicators
      # Each *ROW* is a different lambda
      # Each *COLUMN* is a different predictor
      edge_indicators <- Mod(coef(cfit)[-1, , drop=FALSE]) != 0

      ## This sets things up so that each *column* is part of the adjacency matrix
      ## for a specific lambda, appropriately zero padded so that these can be combined
      ## in later steps.
      ##
      ## When we "rbind" these, each *column* will be a "vec-d out" adjacency matrix
      ## so some basic dimension shifting will give a "cube" of adjacency matrices
      if(j == 1) {
        rbind(0, edge_indicators)
      } else {
        rbind(edge_indicators[seq_len(j - 1), , drop = FALSE],
              0,
              edge_indicators[-seq_len(j - 1), , drop = FALSE])
      }

  }
  edge_matrices <- do.call(rbind, edge_matrices)
  edge_array <- array(edge_matrices,
                      dim = c(nvar, nvar, length(lambda)),
                      dimnames = list(colnames(X), colnames(X), NULL))
  sym_func <- if(rule == "and") ceiling else floor

  ## This gives us a NVAR x NVAR x NLAMBDA array, each slice of which along the third
  ## axis gives a graph
  ##
  ## The symmetrizing functions works because AND on {0, 1}^2 is FLOOR((x + y)/2)
  ## and OR is CEILING
  ##
  ## The call to aperm gives the "flipped" graph so we have the appropriate
  ## thing to match against
  edge_array <- sym_func((edge_array + aperm(edge_array, c(2, 1, 3))) / 2)

  list(edges  = edge_array,
       lambda = lambda)
}

cglasso_cglasso <- function(X, nlambda = nlambda,
                            lambda.min.ratio = lambda.min.ratio,
                            lambda = lambda){
  ## Standardize X here for consistency with neighborhood selection
  X <- scale(X)

  nobs <- NROW(X)
  nvar <- NCOL(X)

  ## Since X is centered, we get the covariance by
  Sigma_hat <- cplx_crossprod(X) / (nobs - 1)

  if(missing(lambda)){
    ## Lambda_max for the glasso problem is the maximum off-diagonal element of Sigma_hat
    ## So elementwise multiply Sigma hat with a 0/1 mask and take the maximum element
    lambda_max <- max(Mod(Sigma_hat * (1 - eye(nvar))))
    lambda <- logspace(lambda.min.ratio * lambda_max, lambda_max, length.out=nlambda)
  }

  cglasso_admm_cpp(Sigma_hat, lambda)
}


#' @export
print.CGLassoFit <- function(x, ..., indent=0){
  icat("Complex Graphical Lasso Fit", "\n", indent=indent)
  icat("-------------------", "\n", indent=indent)
  icat("\n", indent=indent)
  icat("N: ", NROW(x$X), ". P: ", NCOL(x$X), ".\n", sep="", indent=indent)
  icat("\n", indent=indent)
  icat("Grid:", length(x$lambda), "values of lambda. \n", indent=indent)
  icat("  Miniumum:", min(x$lambda), "\n", indent=indent)
  icat("  Maximum: ", max(x$lambda), "\n", indent=indent)
  icat("  Number of edges:", min(x$nnz), " --> ", max(x$nnz), "\n", indent=indent)
  icat("\n", indent=indent)
  icat("Fit Options:\n", indent=indent)
  if(x$method == "mb"){
    icat("  - Estimator: ", "Neighborhood Selection (Pseudo-Likelihood)", "\n", indent=indent)
    icat("  - Symmetrization Rule: ", toupper(x$rule), "\n", indent=indent)
  } else {
    icat("  - Estimator: ", "Full Likelihood (Graphical Lasso)", "\n", indent=indent)
  }
  icat("  - Propriety (Circular Symmetry) Enforced: ", x$proper, "\n", indent=indent)
  icat("\n", indent=indent)
  icat("Time: ", sprintf("%2.3f %s", x$time, attr(x$time, "units")), "\n", indent=indent)
  icat("\n", indent=indent)

  invisible(x)
}
