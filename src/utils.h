#ifndef UTILS
#define UTILS

#include <RcppArmadillo.h>

using namespace Rcpp;

double FindBoundaryCoef(const arma::rowvec& x, const arma::rowvec& dk, const arma::rowvec& weights, bool& shrink, std::vector<unsigned int>& shrink_indices);

void updateHessian(arma::dmat& inverseHessian, const arma::rowvec& mu_diff, const arma::rowvec& grad_diff, const arma::dmat& I);

#endif
