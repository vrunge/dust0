#ifndef FLAT_GAUSS_OP_MD
#define FLAT_GAUSS_OP_MD

#include <RcppArmadillo.h>

using namespace Rcpp;

double Cost_MD(const unsigned int& t, const unsigned int& s, const arma::dmat& cumsum);

double statistic_MD(const double& value);

#endif

