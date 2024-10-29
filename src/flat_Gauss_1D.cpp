#include "flat_Gauss_1D.h"

using namespace Rcpp;

double Cost_1D(const unsigned int& t, const unsigned int& s, const std::vector<double>& cumsum)
{
  return - .5 * pow(cumsum[t] - cumsum[s], 2) / (t - s);
}

double statistic_1D(const double& value)
{
  return value;
}

