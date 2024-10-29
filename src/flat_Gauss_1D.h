#ifndef FLAT_GAUSS_OP_1D
#define FLAT_GAUSS_OP_1D

#include <Rcpp.h>

using namespace Rcpp;

double Cost_1D(const unsigned int& t, const unsigned int& s, const std::vector<double>& cumsum);

double statistic_1D(const double& value);

#endif

