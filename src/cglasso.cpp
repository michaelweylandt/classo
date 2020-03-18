// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "classo.h"

arma::cx_mat soft_thresh(const arma::cx_mat& X, double lambda){
  arma::uword p = X.n_rows;
  arma::cx_mat result(p, p);
  for(arma::uword i = 0; i < p; i++){
    for(arma::uword j = 0; j < p; j++){
      if(i == j){
        result(i, j) = X(i, j);
      } else {
        result(i, j) = soft_thresh(X(i, j), lambda);
      }
    }
  }

  return result;
}


// [[Rcpp::export]]
Rcpp::List cglasso_admm_cpp(const arma::cx_mat& S,
                            const arma::vec& lambda,
                            double thresh = 1e-7){

  const arma::uword p = S.n_cols; // S is the empirical covariance matrix
  const arma::uword n_lambda = lambda.n_elem;

  // Storage for results
  arma::cx_cube Theta_hat(p, p, n_lambda, arma::fill::zeros);
  arma::ucube   edges(p, p, n_lambda, arma::fill::zeros);

  // Working variables
  arma::cx_mat X = S + arma::eye(p, p); // Ensure we start at a feasible point
  arma::cx_mat Z = X;
  arma::cx_mat Z_old = Z;
  arma::cx_mat U(p, p, arma::fill::zeros);

  // Number of ADMM iterations -- used to check for interrupts
  arma::uword k = 0;

  // This approach doesn't take advantage of sparsity, so no point going in reverse
  for(arma::uword i = 0 ; i < n_lambda; i++){
    const double cur_lambda = lambda(i);

    do {
      Z_old = Z;
      arma::cx_mat eigen_vecs;
      arma::vec eigen_vals;

      eig_sym(eigen_vals, eigen_vecs, Z - U - S);
      eigen_vals = (eigen_vals + arma::sqrt(eigen_vals % eigen_vals + 4)) / 2.0;

      X = eigen_vecs * arma::diagmat(eigen_vals) * eigen_vecs.t();
      Z = soft_thresh(X + U, cur_lambda);
      U += X - Z;

      k++; // Increment loop counter

      if((k % CLASSO_CHECK_USER_INTERRUPT_RATE_ADMM) == 0){
        Rcpp::checkUserInterrupt();
      }
      if(k >= CLASSO_MAX_ITERATIONS_ADMM * n_lambda){
        Rcpp::stop("[CLasso::cglasso] Maximum number of ADMM iterations reached.");
      }

    } while( arma::norm(Z - Z_old) > thresh );

    // FIXME - We know diag(Theta_hat) has to be real and positive
    // but floating point noise is giving small imaginary components. Can we
    // enforce that here?

    Theta_hat.slice(i) = Z;

    edges.slice(i) = arma::abs(Z) > 0;
    edges.slice(i).diag().zeros();
  }

  return Rcpp::List::create(Rcpp::Named("Theta_hat") = Theta_hat,
                            Rcpp::Named("edges") = edges,
                            Rcpp::Named("lambda") = lambda);
}
