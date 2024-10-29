#include "flat_Gauss_MD.h"

using namespace Rcpp;

double CostGauss_MD(const unsigned int& t, const unsigned int& s, const arma::dmat& cumsum)
{
  double result = 0;
  for (unsigned int row = 0; row < cumsum.n_rows; row++)
    result += pow(cumsum(row, t) - cumsum(row, s), 2);

  return - .5 * result / (t - s);
}




