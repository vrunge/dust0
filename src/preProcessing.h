// [[Rcpp::depends(RcppArmadillo)]]

#ifndef PREP_H
#define PREP_H

#include <RcppArmadillo.h>

using namespace Rcpp;


// --------- // MAD Estimator // --------- //
//
// Provided some data vector, computes the Median Absolute Deviation estimator
//
// Parameters:
//  - data (vector): a vector of numeric values

double madEstimator
(
    NumericVector& data
);

#endif