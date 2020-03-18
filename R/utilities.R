#################################################################
#   Utilities for working with complex valued data
#   -- design is somewhat similar to the 'cmvnorm' and 'mvtnorm' packages
#   -- for now, we are assuming the circularly symmetric complex normal
#        for the multivariate case, so just
################################################################

eye <- function(n){
    diag(1, nrow=n, ncol=n)
}

#' Complex Utilities
#'
#' A collection of small utilities useful for working with complex-valued data.
#'
#' @details The following helper functions are provided:
#' \itemize{
#' \item \code{im} The imaginary unit \eqn{\sqrt{-1}}.
#' \item \code{h} The conjugate transpose function
#' \item \code{cplx_crossprod} and \code{cplx_tcrossprod} work like
#'       the \code{\link[base]{crossprod}} and \code{\link[base]{tcrossprod}} functions
#'       in the \code{base} package, but use the Hermitian (conjugate) transpose
#'       instead of the standard transpose.
#' }
#'
#' @rdname classo_utilities
#' @export
im <- complex(real=0, imaginary=1)

# Hermitian conjugate
#' @rdname classo_utilities
#' @export
#' @param X A (possibly complex) matrix
h <- function(X){ Conj(t(X)) }

## Versions of "isSymmetric" and "isHermitian" that actually do what they say
## for complex matrices

is_symmetric <- function(Z){
    if(!is.matrix(Z)){
        stop(dQuote("is_symmetric"), " only defined for matrices.")
    }
    isTRUE(all.equal(Z, t(Z)))
}

is_hermitian <- function(Z){
    if(!is.matrix(Z)){
        stop(dQuote("is_hermitian"), " only defined for matrices.")
    }
    isTRUE(all.equal(Z, h(Z)))
}

is_pd <- function(Z){
    if(!is.matrix(Z)){
        stop(dQuote("is_pd"), " only defined for matrices.")
    }
    if(is.complex(Z) && !is_hermitian(Z)){
        stop(dQuote("Z"), " must be Hermitian to check positive-definiteness.")
    }
    if(is.numeric(Z) && !is_symmetric(Z)){
        stop(dQuote("Z"), " must be symmetric to check positive-definiteness.")
    }

    all(eigen(Z, symmetric = TRUE, only.values = TRUE)$values > 0)
}

is_aug_cov <- function(Z){
  if(!is.matrix(Z)){
    stop(dQuote("is_aug_cov"), " only defined for matrices.")
  }

  two_p <- NCOL(Z)

  if(!is_even(two_p)){
    stop(dQuote("Z"), " must have even dimensions to be an augmented covariance matrix.")
  }

  if(!is_hermitian(Z)){
    stop(dQuote("Z"), " must be Hermitian to be an augmented covariance matrix.")
  }

  p <- two_p / 2;

  Block11 <- Z[1:p, 1:p]
  Block12 <- Z[1:p, (p+1):two_p]
  Block21 <- Z[(p+1):two_p, 1:p]
  Block22 <- Z[(p+1):two_p, (p+1):two_p]

  if(!isTRUE(all.equal(Block11, Conj(Block22)))){
    # Z is not an augmented covariance matrix: its top-left and bottom-right blocks must be conjugate.
    return(FALSE)
  }

  if(!isTRUE(all.equal(Block12, Conj(Block21)))){
    # Z is not an augmented covariance matrix: its top-right and bottom-left blocks must be conjugate.
    return(FALSE)
  }

  TRUE
}

## From ?is.integer
is_integer <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}

is_even <- function(x){
  is_integer(x / 2)
}

stop_if_not_scalar <- function(x){
  if(length(x) != 1L){
    cl <- sys.call(0)
    nm <- deparse(cl[[2]])
    stop(dQuote(nm), " must be a scalar.")
  }
}

recycle_rows <- function(X, n){
  X[rep(seq(1, NROW(X)), length.out = n), , drop = FALSE]
}

#' @param x A complex matrix with no non-finite elements
#' @param y A complex matrix with no non-finite elements
#' @rdname classo_utilities
#' @export
cplx_crossprod <- function(x, y = x){
  cplx_crossprod_cpp(x, y)
}

#' @rdname classo_utilities
#' @export
cplx_tcrossprod <- function(x, y = x){
  cplx_tcrossprod_cpp(x, y)
}

vec_real_to_complex <- function(vec){
  n <- length(vec)

  x_ix <- 1:(n/2); y_ix <- (n/2 + 1):n

  vec[x_ix] + im * vec[y_ix]
}


mat_real_to_complex <- function(mat){
  p <- NCOL(mat)

  x_ix <- 1:(p/2); y_ix <- (p/2+1):p

  mat[,x_ix] + im * mat[,y_ix]
}

#-- Unused utilities from earlier versions --

# is_valid_cov_mat <- function(Sigma){
#
#     if(!isHermitian(Sigma, tol=sqrt(.Machine$double.eps))){
#         return(FALSE)
#     }
#
#     return(all(eigen(Sigma, symmetric=TRUE, only.values=TRUE)$values > 0))
# }
#
# cov_mat_real_to_complex <- function(Sigma){
#     p <- NROW(Sigma)/2
#
#     x_ix <- 1:p; y_ix <- (p+1):(2*p)
#     # Following wikipedia
#     # V_{xx} + V_{yy} + i(V_{yx} - V_{xy})
#     Sigma[x_ix, x_ix] + Sigma[y_ix, y_ix] + im * (Sigma[y_ix, x_ix] - Sigma[x_ix, y_ix])
# }
#
# cov_mat_complex_to_real <- function(Sigma){
#     p <- NROW(Sigma)
#
#     V <- matrix(0, nrow=2*p, ncol=2*p)
#
#     x_ix <- 1:p; y_ix <- (p+1):(2*p)
#
#     V[x_ix, x_ix] = 1/2 * Re(Sigma)
#     V[x_ix, y_ix] = 1/2 * Im(-Sigma)
#     V[y_ix, x_ix] = 1/2 * Im(Sigma)
#     V[y_ix, y_ix]  = 1/2 * Re(Sigma)
#
#     V
# }
#
#
# vec_complex_to_real <- function(vec){
#     c(Re(vec), Im(vec))
# }
#
#
# mat_complex_to_real <- function(mat){
#     rbind(Re(mat), Im(mat))
# }
#
# complex_cov <- function(X, mean){
#     if(is.vector(X)){
#         X <- matrix(X, ncol=1)
#     }
#
#     if(missing(mean)){
#         X <- X - colMeans(X)
#         h(X) %*% X / (NROW(X) - 1)
#     } else {
#         X <- X - mean(X)
#         h(X) %*% X / NROW(X)
#     }
#
# }
#
# matrix_trace <- function(A) sum(diag(A))
#
# log_det <- function(A) log(prod(eigen(A, only.values=TRUE, symmetric=TRUE)$values))
