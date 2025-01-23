#ifndef UTILS
#define UTILS

#include <RcppArmadillo.h>

using namespace Rcpp;

void clip_stepsize_to_negative_element(const arma::rowvec& mu,
                                  const arma::rowvec& direction,
                                  double& max_stepsize);

void clip_stepsize_to_negative_sum(const std::vector<int>& sign, const double& mu_sum, const arma::rowvec& direction, double& direction_sum, double& max_stepsize);

double FindBoundaryCoef(const arma::rowvec& x,
                        const arma::rowvec& dk,
                        const arma::rowvec& weights,
                        bool& shrink,
                        std::vector<unsigned int>& shrink_indices);

void updateHessian(arma::dmat& inverseHessian,
                   const arma::rowvec& mu_diff,
                   const arma::rowvec& grad_diff,
                   const arma::dmat& I);

#endif
