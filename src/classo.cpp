// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "classo.h"

// [[Rcpp::export]]
Rcpp::List classo_cd(const arma::cx_mat& X,
                     const arma::cx_vec& y,
                     const arma::cx_vec& o,
                     const arma::vec& w,
                     const arma::vec& lambda,
                     double alpha = 1, // Elastic Net Mixing parameter: Lasso when alpha is 1
                     double thresh=1e-7,
                     bool intercept=true){

    const arma::uword n = X.n_rows;
    const arma::uword p = X.n_cols;
    const arma::uword n_lambda = lambda.n_elem;

    // Working copies of coefficients
    arma::cx_vec r = y - o;
    arma::cx_vec beta_working(p, arma::fill::zeros);
    arma::cx_vec beta_old(p, arma::fill::zeros);
    arma::cx_double gamma = 0;
    arma::cx_double gamma_old = 0;

    arma::vec u(p);
    for(arma::uword i=0; i<p; i++){
        u(i) = arma::sum(arma::square(arma::abs(X.col(i))) % w);
    }

    // Storage for intercepts
    arma::cx_vec Gamma(n_lambda, arma::fill::zeros); // Storage for intercepts

    // R's Matrix package doesn't support complex sparse matrices, hence
    // RcppArmadillo can't handle arma::sp_cx_mat, so we use a dense matrix
    // for coefficients here. It costs a bit in memory, but saves performance
    // from constantly having to swap back and forth from dense and sparse matrices
    arma::cx_mat Beta(p, n_lambda, arma::fill::zeros);

    // Number of cd iterations -- used to check for interrupts
    arma::uword k = 0;

    // For first iteration we want to loop over all variables since
    // we haven't identified the active set yet
    bool full_loop = true;
    arma::uword full_loop_count = 0; // Number of full loops completed
                                     // We require at least CLASSO_FULL_LOOP_MIN
                                     // full loops before moving to the next value
                                     // of lambda to ensure convergence

    // Iterate from highest to smallest lambda
    // to take advantage of
    // warm starts for sparsity
    //
    // Note the non-traditional loop header: https://stackoverflow.com/a/4206815
    for(arma::uword i=n_lambda; i-- > 0;){
        const double cur_lambda = lambda(i);

        do {
            beta_old = beta_working;
            for(arma::uword j=0; j < p; j++){
                arma::cx_double beta = beta_working(j);

                if((!full_loop) && (std::abs(beta) == 0)){
                    continue;
                }

                const arma::cx_vec xj = X.col(j);

                r += xj * beta;
                beta = soft_thresh(arma::cdot(xj, w % r), n * cur_lambda * alpha) / (u(j) + n * cur_lambda * (1-alpha));
                r -= xj * beta;

                beta_working(j) = beta;
            }

            if(intercept){
                gamma_old = gamma;
                r += gamma;
                gamma = arma::sum(r % w) / (n + 0.0);
                r -= gamma;
            }

            k++; // Increment loop counter

            if((k % CLASSO_CHECK_USER_INTERRUPT_RATE_CD) == 0){
                Rcpp::checkUserInterrupt();
            }
            if(k >= CLASSO_MAX_ITERATIONS_CD * n_lambda){
                Rcpp::stop("[CLasso::classo] Maximum number of coordinate descent iterations reached.");
            }

            if(full_loop){
                full_loop_count++;
                // Only check this if not already in full loop mode
            } else if(arma::norm(beta_working - beta_old) + std::abs(gamma - gamma_old) < CLASSO_FULL_LOOP_FACTOR * thresh){
                // If it looks like we're closing in on a solution,
                // switch to looping over all variables
                full_loop = true;
            }

        } while((full_loop_count < CLASSO_FULL_LOOP_MIN) || std::abs(gamma - gamma_old) + arma::norm(beta_working - beta_old) > thresh);

        // Switch back to active set loops for next iteration
        full_loop = false;
        full_loop_count = 0;

        Gamma(i) = gamma;

        // Load sparse matrix storage
        for(arma::uword j = 0; j < p; j++){
            if((std::real(beta_working(j)) != 0) || (std::imag(beta_working(j)) != 0)){
                Beta(j, i) = beta_working(j);
            }
        }
    }

    if(intercept){
        return Rcpp::List::create(Rcpp::Named("intercept")=Gamma,
                                  Rcpp::Named("coef")=Beta,
                                  Rcpp::Named("k")=k);
    } else {
        return Rcpp::List::create(Rcpp::Named("intercept")=R_NilValue,
                                  Rcpp::Named("coef")=Beta,
                                  Rcpp::Named("k")=k);
    }
}
