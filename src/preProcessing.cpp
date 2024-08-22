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



// Function to compute the median of a vector
double median(std::vector<double> vec) {
  size_t n = vec.size();
  if (n == 0) throw std::runtime_error("Cannot compute median of an empty vector");

  std::sort(vec.begin(), vec.end());

  if (n % 2 == 0) {
    return (vec[n / 2 - 1] + vec[n / 2]) / 2.0;
  } else {
    return vec[n / 2];
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////


std::vector<double> transformVector(const std::vector<double>& y) {
  // Check if the vector has at least two elements
  if (y.size() < 2) {
    throw std::runtime_error("Vector must contain at least two elements.");
  }

  // Vector to store the differences
  std::vector<double> differences(y.size() - 1);

  // Compute the differences
  std::adjacent_difference(y.begin(), y.end(), differences.begin());

  // Remove the first element of the result (which is the first element of y itself)
  differences.erase(differences.begin());

  // Scaling factor
  double scale = std::sqrt(2);

  // Scale the differences
  for (double& diff : differences) {
    diff /= scale;
  }
  return differences;
}



// [[Rcpp::export]]
double mad_std(const std::vector<double>& data) {
  if (data.size() < 2) throw std::runtime_error("Insufficient data to compute MAD");

  std::vector<double> dataTrans = transformVector(data);
  // Compute the median of the data
  double medianValue = median(dataTrans);

  // Compute the absolute deviations from the median
  std::vector<double> absDeviations;
  for (double value : data) {
    absDeviations.push_back(std::fabs(value - medianValue));
  }

  // Compute the median of the absolute deviations
  double madValue = median(absDeviations);

  // Scale factor to convert MAD to standard deviation
  //const double scaleFactor = 1.4826;

  // Return the scaled MAD
  return madValue;
}

