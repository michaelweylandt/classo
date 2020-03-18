# FIXME - For all mean parameters, should actually be 0+0i, but interacts oddly
#         with roxygen2. Figure out why

#' Density Function and RNG for the Univariate Complex Gaussian Distribution
#'
#' Currently, these functions are vectorized in \code{x} and \code{mean} but
#' \code{cov} must be a single \code{\link{cplx_cov_spec}}.
#'
#' @param x A complex vector of observations
#' @param mean Complex mean vector (default) zero.
#' @param cov A covariance specification of class \code{\link{cplx_cov_spec}}
#' @param log Logical. Should densities be returned on the log scale?
#' @export
#' @rdname univariate_distribution
#' @examples
#' z <- rcnorm(250, mean = 2 + 3i)
#' mean(z)
#' all.equal(dcnorm(z), dcnorm(z + 2i, mean = 2i))
dcnorm <- function(x, mean = 0, cov = cplx_cov_spec(Sigma = 1), log = FALSE){
  # Manual R-type recycling
  n_x <- length(x); n_mean = length(mean); n <- max(n_x, n_mean)

  x_recy    <- rep(x, length.out = n)
  mean_recy <- rep(mean, length.out = n)
  W      <- cov$W

  if(dim(cov) != 1){
    stop(dQuote("cov"), " must be a univariate covariance specification.")
  }

  na_mask <- is.na(x_recy) | is.na(mean_recy)

  x_recy[na_mask] <- mean_recy[na_mask] <- 0

  res <- dcnorm_cpp(x = x_recy, mean = mean_recy, W = W, log = log)
  res[na_mask] <- NA
  dim(res) <- dim(x)

  res
}

#' @rdname univariate_distribution
#' @param n The number of samples to generate
#' @export
#' @importFrom mvtnorm rmvnorm
rcnorm <- function(n, mean = 0, cov = cplx_cov_spec(Sigma = 1)){
  if(dim(cov) != 1){
    stop(dQuote("cov"), " must be a univariate covariance specification.")
  }

  W      <- cov$W
  vec_real_to_complex(rmvnorm(n, mean = c(0, 0), sigma = W)) + mean
}

#' Density Function and RNG for the Multivariate Complex Gaussian Distribution
#'
#' Currently, these functions are vectorized in \code{x} and \code{mean} but
#' \code{cov} must be a single \code{\link{cplx_cov_spec}}.
#'
#' @param x A complex vector of observations. Rows of \code{x} are separate observations.
#' @param mean Complex mean vector (default) zero. Rows of \code{mean} are separate observations.
#' @param cov A covariance specification of class \code{\link{cplx_cov_spec}}
#' @param log Logical. Should densities be returned on the log scale?
#' @export
#' @rdname multivariate_distribution
#' @examples
#' cov <- cplx_cov_spec(Sigma = diag(5))
#' Z   <- rmvcnorm(500, cov = cov)
#' den <- dmvcnorm(Z, cov = cov)
dmvcnorm <- function(x, mean = matrix(0, ncol = p),
                     cov = cplx_cov_spec(Sigma = diag(1, p)), log = FALSE){
  if(any(is.na(x))){
    stop("NAs in ", dQuote("x"), " not supported.")
  }

  if(!is.matrix(x)){
    x <- matrix(x, ncol = length(x))
  }

  p <- NCOL(x)

  if(!is.matrix(mean)){
    mean <- matrix(mean, ncol = length(mean))
  }

  if(any(is.na(mean))){
    stop("NAs in ", dQuote("mean"), " not supported.")
  }

  if(NCOL(mean) != p){
    stop("Dimensions of ", dQuote("x"), " and ", dQuote("mean"), " do not match.")
  }

  if(dim(cov) != p){
    stop("Dimensions of ", dQuote("x"), " and ", dQuote("cov"), " do not match.")
  }

  n_x <- NROW(x); n_mean <- NROW(mean); n <- max(n_x, n_mean)
  x    <- recycle_rows(x, n)
  mean <- recycle_rows(mean, n)
  W      <- cov$W

  as.vector(dmvcnorm_cpp(x = x, mean = mean, W = W, log = log))
}

#' @export
#' @param n The number of samples to generate
#' @rdname multivariate_distribution
rmvcnorm <- function(n, mean = matrix(0, ncol = p),
                     cov = cplx_cov_spec(Sigma = diag(1, NCOL(mean)))){
  if(missing(mean) && missing(cov)){
    stop("At least one of ", dQuote("mean"), " and ", dQuote("cov"), " must be supplied.")
  }

  if(!missing(cov)){
    p <- dim(cov)
  }

  if(!is.matrix(mean)){
    mean <- matrix(mean, ncol = length(mean))
  }

  W <- cov$W
  p <- dim(cov)

  mat_real_to_complex(rmvnorm(n, mean = rep(0, 2 * p), sigma = W)) + recycle_rows(mean, n)
}
