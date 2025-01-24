#ifndef PREP_H
#define PREP_H

#include <Rcpp.h>

using namespace Rcpp;

// --------- // MAD Estimator // --------- //
//
// Provided some data vector, computes the Median Absolute Deviation estimator
//
// Parameters:
//  - data (vector): a vector of numeric values

double sdDiff
(
    std::vector<double>& y, std::string method = "HALL"
);

#endif
