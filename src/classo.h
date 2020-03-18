#include <limits>
#include <cmath>
#include <complex>

#define CLASSO_CHECK_USER_INTERRUPT_RATE_CD 500
#define CLASSO_MAX_ITERATIONS_CD 50000
#define CLASSO_FULL_LOOP_FACTOR 10
#define CLASSO_FULL_LOOP_MIN 2
#define CLASSO_CHECK_USER_INTERRUPT_RATE_ADMM 500
#define CLASSO_MAX_ITERATIONS_ADMM 50000

// [[Rcpp::depends(RcppArmadillo)]]

// We only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

arma::cx_double soft_thresh(arma::cx_double, double);

