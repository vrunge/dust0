// [[Rcpp::depends(Rcpp)]]

#include <Rcpp.h>

using namespace Rcpp;


// --------- // MAD Estimator // --------- //
//
// Provided some data vector, computes the Median Absolute Deviation estimator
//
// Parameters:
//  - data (vector): a vector of numeric values

// [[Rcpp::export]]
double madEstimator(NumericVector& data) {

  int n = data.size();
  int medianIndex = n / 2 + n % 2;

  NumericVector sorted_data(data.begin(), data.end());
  std::sort(sorted_data.begin(), sorted_data.end());

  double median = sorted_data[medianIndex];

  double mad = 0.;
  NumericVector::iterator y = data.begin();

  do
  {
    mad += fabs(*y - median);
    ++y;
  }
  while (y != data.end());

  return mad / n;
}
