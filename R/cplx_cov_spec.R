#' Complex Covariance Specification
#'
#' The \code{cplx_cov_spec} class encodes the second order moments (covariance
#' and relation matrices) of a potentially improper Gaussian random vector.
#' To specify a covariance, three interfaces are supported, but only one can be
#' used at a time: \itemize{
#' \item \code{Sigma} (and optionally \code{R}): directly specify \eqn{\Sigma=E[ZZ^H]}
#' and optionally \eqn{R=E[ZZ^T]}.
#' \item \code{W}: specify the joint covariance of the real and imaginary parts
#' of \eqn{Z} via \eqn{W=\text{Cov}([R(Z) I(Z)])}.
#' \item {A}: specify the joint covariance of \eqn{Z} and \eqn{\overline{Z}}
#' via \eqn{A = \text{Cov}[Z \overline{Z}]}. This form is not recommended as it
#' puts somewhat tricky constraints on A.}
#'
#' @param Sigma A p-by-p complex covariance matrix (\eqn{E[ZZ^H]})
#' @param R A p-by-p relation matrix. If not provided, assumed zero (proper).
#'          (\eqn{E[ZZ^T]})
#' @param W A (2p)-by-(2p) real covariance matrix. (\eqn{E[([R(Z) I(Z)])([R(Z) I(Z)])^T]})
#' @param A A (2p)-by-(2p) augmented complex covariance matrix. (\eqn{E[([Z \overline{Z}])([Z \overline{Z}])^H]})
#' @param ... Additional arguments, currently disallowed.
#' @export
#' @examples
#' # Standard univariate complex Gaussian
#' (cov_spec <- cplx_cov_spec(1))
#'
#' # Standard multivariate complex Gaussian
#' (cov_spec <- cplx_cov_spec(diag(1, 10)))
#'
#' # Specify the joint distribution of real and imaginary parts
#' W <- matrix(c(4, 0, 0, 2), ncol = 2, nrow = 2)
#' cov_spec <- cplx_cov_spec(W = W)
#' # Because the real and imaginary parts have different variance, W is improper
#' print(cov_spec)
cplx_cov_spec <- function(Sigma, ..., R = matrix(0, NROW(Sigma), NCOL(Sigma)), W, A){
  chkDots(...)

  if(!missing(Sigma)){
    ## Simplest constructor:
    ## i) Check if Sigma is Hermitian and PD
    ## ii) Generate R (by default if needed) and check if symmetric
    ## iii) Construct W directly

    if(!missing(W)){
      stop("Cannot supply both ", dQuote("Sigma"), " and ", dQuote("W."))
    }
    if(!missing(A)){
      stop("Cannot supply both ", dQuote("Sigma"), " and ", dQuote("A."))
    }

    Sigma <- as.matrix(Sigma); storage.mode(Sigma) <- "complex"
    R     <- as.matrix(R);     storage.mode(R)     <- "complex"

    if(!is_pd(Sigma)){
      stop(dQuote("Sigma"), " must be Hermitian and positive-definite.")
    }

    # If R not specified, will be auto-generated via default argument
    if(!all(dim(Sigma) == dim(R))){
      stop("Dimensions of ", dQuote("Sigma"), " and ", dQuote("R"), " do not match.")
    }

    # Build W from Sigma and R
    p <- NCOL(Sigma)
    W <- matrix(NA_real_, nrow = 2 * p, ncol = 2 * p)
    W[1:p, 1:p] <- 0.5 * Re(Sigma + R)
    W[1:p, (p+1):(2*p)] <- 0.5 * Im(-Sigma + R)
    W[(p+1):(2*p), 1:p] <- 0.5 * Im(Sigma + R)
    W[(p+1):(2*p), (p+1):(2*p)] <- 0.5 * Re(Sigma - R)

  } else if(!missing(W)){
    if(!missing(Sigma)){
      stop("Cannot supply both ", dQuote("W"), " and ", dQuote("Sigma."))
    }

    if(!missing(R)){
      stop("Cannot supply both ", dQuote("W"), " and ", dQuote("R."))
    }

    if(!missing(A)){
      stop("Cannot supply both ", dQuote("W"), " and ", dQuote("A."))
    }

    if(!is_pd(W)){
      stop(dQuote("W"), " must be Hermitian and positive-definite.")
    }

    if(!is.numeric(W)){
      stop(dQuote("W"), " must be a real matrix.")
    }

    two_p <- NCOL(W)

    if(!is_even(two_p)){
      stop(dQuote("W"), " must have even dimensions to be a complex covariance.")
    }

    p <- two_p / 2
    Sigma <- R <- matrix(NA_complex_, p, p)

    Var_X  <- W[1:p, 1:p, drop = FALSE]
    Var_Y  <- W[(p+1):two_p, (p+1):two_p, drop = FALSE]
    Cov_XY <- W[1:p, (p+1):two_p, drop = FALSE]

    Sigma <- Var_X + Var_Y + im * (t(Cov_XY) - Cov_XY)
    R     <- Var_X - Var_Y + im * (t(Cov_XY) + Cov_XY)

  } else if(!missing(A)){
    if(!missing(Sigma)){
      stop("Cannot supply both ", dQuote("A"), " and ", dQuote("Sigma."))
    }

    if(!missing(R)){
      stop("Cannot supply both ", dQuote("A"), " and ", dQuote("R."))
    }

    if(!missing(W)){
      stop("Cannot supply both ", dQuote("A"), " and ", dQuote("W."))
    }

    if(!is_aug_cov(A)){
      stop(dQuote("A"), " must satisfy the augmented covariance conditions.")
    }

    if(!is.numeric(A)){
      stop(dQuote("A"), " must be a real matrix.")
    }

    two_p <- NCOL(A)

    if(!is_even(two_p)){
      stop(dQuote("A"), " must have even dimensions to be a complex covariance.")
    }

    p <- two_p / 2

    Sigma <- A[1:p, 1:p, drop = FALSE]
    R     <- A[1:p, (p+1):two_p, drop = FALSE]

    W <- matrix(NA_real_, nrow = two_p, ncol = two_p)
    W[1:p, 1:p] <- 0.5 * Re(Sigma + R)
    W[1:p, (p+1):(2*p)] <- 0.5 * Im(-Sigma + R)
    W[(p+1):(2*p), 1:p] <- 0.5 * Im(Sigma + R)
    W[(p+1):(2*p), (p+1):(2*p)] <- 0.5 * Re(Sigma - R)
  }

  RES <- list(Sigma = Sigma, R = R, W = W)
  class(RES) <- "cplx_cov_spec"

  if(all(R == 0)){
    class(RES) <- c("proper_cplx_cov_spec", class(RES))
  }

  return(RES)
}

#' @export
print.proper_cplx_cov_spec <- function(x, ...){
  cat("Complex Covariance Specification of a proper ", dim(x), "-vector.\n", sep = "")
  cat("\nSigma = \n")
  print(x$Sigma, ...)
  invisible(x)
}

#' @export
print.cplx_cov_spec <- function(x, ...){
  cat("Complex Covariance Specification of an improper ", dim(x), "-vector.\n", sep = "")
  cat("\nSigma = \n")
  print(x$Sigma, ...)
  cat("\n\nR = \n")
  print(x$R, ...)
  invisible(x)
}

is_proper <- function(x){
  inherits(x, "proper_cplx_cov_spec")
}

#' @export
dim.cplx_cov_spec <- function(x){
  NCOL(x$Sigma)
}
