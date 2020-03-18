// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "classo.h"

// R's crossprod() and tcrossprod() use a transpose instead of a conjugate transpose
// even for complex matrices, so we define our own versions which make more sense
// for complex matrices.
// [[Rcpp::export]]
arma::cx_mat cplx_crossprod_cpp(const arma::cx_mat& X,
                                const arma::cx_mat& Y){
  return X.t() * Y;
}

// [[Rcpp::export]]
arma::cx_mat cplx_tcrossprod_cpp(const arma::cx_mat& X,
                                 const arma::cx_mat& Y){
  return X * Y.t();
}

// [[Rcpp::export]]
arma::vec dcnorm_cpp(const arma::cx_rowvec& x,
                     const arma::cx_rowvec& mean,
                     const arma::mat& W,
                     bool  log = false){
  // Fast computation of the (univariate) complex Gaussian PDF
  //
  // Vectorized in x and mean; not in W
  // We assume all recycling and size checking happens in R code
  //
  // Internally, this works by converting to the real and imaginary parts
  // and computing the density of the two-dimensional real Gaussian.
  //
  // Implementation of 2-variate density follows
  // https://gallery.rcpp.org/articles/dmvnorm_arma/

  arma::uword n = x.n_elem;
  arma::vec log_dens(n);

  const arma::mat chol_inv  = arma::inv(arma::trimatu(arma::chol(W)));
  // This is 1/sqrt(det(Sigma)) with the sqrt coming from chol and the  the 1/ coming from arma::inv above
  const double sqrt_log_det   = arma::sum(arma::log(chol_inv.diag()));
  // In the case of a two vector, normalizing constant is just 1/(2pi)
  const double log_norm_const = - std::log(2.0 * arma::datum::pi);
  const double log_const = log_norm_const + sqrt_log_det;

  const arma::cx_rowvec z = x - mean;
  arma::mat z_r(2, n);

  z_r.row(0) = arma::real(z);
  z_r.row(1) = arma::imag(z);

  const arma::mat z_r_whitened = chol_inv * z_r;

  for(arma::uword i = 0; i < n; i++){
    log_dens(i) = log_const - 0.5 * arma::dot(z_r_whitened.col(i), z_r_whitened.col(i));
  }

  if(log){
    return log_dens;
  } else {
    return arma::exp(log_dens);
  }
}


// [[Rcpp::export]]
arma::vec dmvcnorm_cpp(const arma::cx_mat& x,
                       const arma::cx_mat& mean,
                       const arma::mat& W,
                       bool  log = false){
  // Fast computation of the (multivariate) complex Gaussian PDF
  //
  // Vectorized along rows of x and mean; not in W
  // We assume all recycling and size checking happens in R code
  //
  // Internally, this works by converting to the real and imaginary parts
  // and computing the density of the two-dimensional real Gaussian.
  //
  // Implementation of 2-variate density follows
  // https://gallery.rcpp.org/articles/dmvnorm_arma/

  arma::uword n = x.n_rows;
  arma::uword p = x.n_cols;
  arma::vec log_dens(n);

  const arma::mat chol_inv  = arma::inv(arma::trimatu(arma::chol(W)));
  // This is 1/sqrt(det(Sigma)) with the sqrt coming from chol and the  the 1/ coming from arma::inv above
  const double sqrt_log_det   = arma::sum(arma::log(chol_inv.diag()));

  // Add 0.0 to convert to double; else -p goes badly
  const double log_norm_const = -(p + 0.0) * std::log(2.0 * arma::datum::pi);
  const double log_const = log_norm_const + sqrt_log_det;

  const arma::cx_mat z = x - mean;
  arma::mat z_r(2 * p, n);

  z_r.rows(0, p - 1)     = arma::real(z).t();
  z_r.rows(p, 2 * p - 1) = arma::imag(z).t();

  const arma::mat z_r_whitened = chol_inv * z_r;

  for(arma::uword i = 0; i < n; i++){
    log_dens(i) = log_const - 0.5 * arma::dot(z_r_whitened.col(i), z_r_whitened.col(i));
  }

  if(log){
    return log_dens;
  } else {
    return arma::exp(log_dens);
  }
}

// Used by both classo and cglasso
arma::cx_double soft_thresh(const arma::cx_double x, double lambda){
  const double r = std::abs(x);
  if(r > lambda){
    return (r - lambda) * x / r;
  } else {
    return arma::cx_double(0.0, 0.0);
  }
}
