// [[Rcpp::depends(RcppArmadillo)]]

#ifndef COSTFUNCTIONS_H
#define COSTFUNCTIONS_H


#include <RcppArmadillo.h>
#include "multipleConstraints.h"
#include <forward_list>
#include <cmath>
#include <random>
#include <limits>
#include <iostream>
#include <fstream>
#include <string>
using namespace Rcpp;
using namespace arma;

// --------- // Sum // --------- //
//
// Sums the values of input vector
//
// Parameters:
//  - data (vector): a vector of numeric values

colvec colPow
(
    colvec& inputCol, int power = 2
);

// --------- // Simple test // --------- //
//
// Performs the simple lagrangian duality test to find the lowest "visible " point
// of a given cost function
// 
// Parameters:
//  - t: the current step in the OP main loop
//  - i: the current step in the pruning sub loop
//  - j: a random index drawn from the available indices, less than i
//  - valuesCumsum (vector): the cumulative sum of the data (index 0 == 0)
//  - costRecord (vector): the optimal model cost at each OP index

double test
(
    int t,
    int i,
    int j,
    const NumericVector& valuesCumsum,
    const NumericVector& costRecord
);

double test1D
(
    int t,
    int i,
    int j,
    const NumericVector& valuesCumsum,
    const NumericVector& costRecord,
    double optimalCost
);

double test1DRand
(
    int t,
    int i,
    int j,
    const NumericVector& valuesCumsum,
    const NumericVector& costRecord,
    double optimalCost
);

double testMV1C
(
    int t,
    int i,
    int j,
    const dmat& valuesCumsum,
    const NumericVector& costRecord,
    double optimalCost
);

double testMV
(
    int t,
    int i,
    MCHandler& mcHandler,
    const arma::dmat& valuesCumsum,
    const NumericVector& costRecord,
    double optimalCost
);


// --------- // modelCost // --------- //
//
// Computes the gaussian model cost in the OP algorithm at step t with cursor at i
// 
// Parameters:
//  - t: the current step in the OP main loop
//  - i: the current step in the OP sub loop
//  - penalty: the penalty value in the penalized changepoint detection model
//  - valuesCumsum (vector): the cumulative sum of the data (index 0 == 0)
//  - costRecord (vector): the optimal model cost at each OP index

double modelCost
(
    int t,
    int i,
    NumericVector& valuesCumsum
);

double modelCostMV
(
    int t,
    int i,
    dmat& valuesCumsum
);


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