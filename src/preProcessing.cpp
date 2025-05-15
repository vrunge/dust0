#include <Rcpp.h>

using namespace Rcpp;

#include <Rcpp.h>
#include <cmath>
#include <algorithm> // for std::max_element, for std::all_of and std::isfinite

#include "preProcessing.h"


//' Calculate Standard Deviation or MAD of Differences in a Numeric Vector
//'
//' The `sdDiff` function calculates a measure of variability (standard deviation or MAD)
//' of a numeric vector.
//' It supports three methods: "HALL", "MAD", and "SD".
//'
//' @param y A numeric vector.
//' @param method A character string specifying the method to use.
//'   Options are: \code{"HALL"}, \code{"MAD"}, and \code{"SD"}.
//'   Default is \code{"HALL"}.
//'
//' @return A numeric value representing the calculated measure of variability
//'   according to the specified method.
//'   \itemize{
//'     \item \code{"HALL"}: Calculates the standard deviation using a specific weighted
//'           difference method (HALL method).
//'     \item \code{"MAD"}: Returns the MAD (Median Absolute Deviation) of the differences
//'           between consecutive elements.
//'     \item \code{"SD"}: Returns the standard deviation of the differences
//'           between consecutive elements.
//'   }
//'
//' @examples
//' y <- c(rnorm(500, mean = 1), rnorm(500,mean = 2))
//' sdDiff(y, "HALL")
//' sdDiff(y, "MAD")
//' sdDiff(y, "SD")
//'
//' @export
// [[Rcpp::export]]
double sdDiff(std::vector<double>& y, std::string method)
{
  ///////////////////////  HALL
  ///////////////////////  HALL
  ///////////////////////  HALL
  if(method == "HALL")
  {
    if(!std::all_of(y.begin(), y.end(), [](double val) { return std::isfinite(val); }) || y.size() < 5)
    {
      Rcpp::stop("y is not a numeric vector or length < 5 (the HALL method cannot be used)");
    }

    int n = y.size();
    Rcpp::NumericVector wei = {0.1942, 0.2809, 0.3832, -0.8582};
    Rcpp::NumericMatrix mat(4, n);

    // Constructing the matrix `mat`
    for(int i = 0; i < 4; i++)
    {
      for(int j = 0; j < n; j++)
      {
        mat(i, j) = wei[i] * y[j];
      }
    }

    // Adjusting the elements according to the R function
    for(int j = 0; j < n - 1; j++){mat(1, j) = mat(1, j + 1);}
    for(int j = 0; j < n - 2; j++){mat(2, j) = mat(2, j + 2);}
    for(int j = 0; j < n - 3; j++){mat(3, j) = mat(3, j + 3);}

    // Computing the result
    double sumSquares = 0.0;
    double columnSum = 0.0;
    for(int j = 0; j < n - 3; j++)
    {
      columnSum = 0.0;
      for(int i = 0; i < 4; i++){columnSum += mat(i, j);}
      sumSquares += columnSum * columnSum;
    }
    return std::sqrt(sumSquares / (n - 3));
  }
  ///////////////////////  MAD
  ///////////////////////  MAD
  ///////////////////////  MAD
  if(method == "MAD")
  {
    if(!std::all_of(y.begin(), y.end(), [](double val) { return std::isfinite(val); }) || y.size() < 2)
    {
      Rcpp::stop("y is not a numeric vector or length < 2 (the MAD method cannot be used)");
    }
    int n = y.size() - 1;

    Rcpp::NumericVector result(n);
    double sqrt2 = std::sqrt(2.0);
    for (int i = 0; i < n ; ++i) {result[i] = (y[i + 1] - y[i]) / sqrt2;}

    // Calculate the median
    std::vector<double> sorted_x(result.begin(), result.end());
    std::nth_element(sorted_x.begin(), sorted_x.begin() + n / 2, sorted_x.end());
    double median_x = sorted_x[n / 2];
    if (n % 2 == 0)
    {
      std::nth_element(sorted_x.begin(), sorted_x.begin() + n / 2 - 1, sorted_x.end());
      median_x = (median_x + sorted_x[n / 2 - 1]) / 2.0;
    }

    // Calculate the absolute deviations from the median
    std::vector<double> abs_dev(n);
    for (int i = 0; i < n; i++)
    {
      abs_dev[i] = std::abs(result[i] - median_x);
    }

    // Calculate the median of the absolute deviations
    std::nth_element(abs_dev.begin(), abs_dev.begin() + n / 2, abs_dev.end());
    double mad = abs_dev[n / 2];
    if (n % 2 == 0)
    {
      std::nth_element(abs_dev.begin(), abs_dev.begin() + n / 2 - 1, abs_dev.end());
      mad = (mad + abs_dev[n / 2 - 1]) / 2.0;
    }
    return mad * 1.4826;
  }

  ///////////////////////  SD
  ///////////////////////  SD
  ///////////////////////  SD
  if(method == "SD")
  {
    int n = y.size() - 1;
    if(!std::all_of(y.begin(), y.end(), [](double val) { return std::isfinite(val); }) || y.size() < 2)
    {
      Rcpp::stop("y is not a numeric vector or length < 2 (the SD method cannot be used)");
    }

    Rcpp::NumericVector diff_y(n);
    double sqrt2 = std::sqrt(2.0);
    for (int i = 0; i < n; ++i) {diff_y[i] = (y[i + 1] - y[i]) / sqrt2;}

    // Calculate the standard deviation
    double mean = Rcpp::mean(diff_y);
    double sum_sq_diff = 0.0;
    for (int i = 0; i < n; ++i){sum_sq_diff += std::pow(diff_y[i] - mean, 2);}

    double variance = sum_sq_diff / (n - 1);
    return std::sqrt(variance);
  }
  return(0);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

//' Data Normalization Function
//'
//' @name data_normalization_1D
//'
//' @description
//' Normalizes the input time series data `y` according to the specified `type`.
//' The normalization process depends on the statistical model type, which can be one of the following:
//' "gauss" (Gaussian/normal distribution), "exp" (exponential distribution),
//' "poisson" (Poisson distribution), "geom" (geometric distribution),
//' "bern" (Bernoulli distribution), "binom" (binomial distribution),
//' "negbin" (negative binomial distribution), or "variance"
//'
//' @param y A numeric vector representing the time series to be normalized and then segmented.
//' @param type A string specifying the model type for normalization.
//' The available options are "gauss", "exp", "poisson", "geom", "bern", "binom", "negbin", "variance".
//' The default is "gauss".
//' @return A numeric vector that is the normalized version of the input time series `y`.
//' @examples
//' # Normalize a random time series using the Gaussian model
//' normalized_y <- data_normalization_1D(rnorm(100), type = "gauss")
//'
//' # Normalize using the Poisson model
//' normalized_y <- data_normalization_1D(rpois(100, lambda = 3), type = "poisson")
//'
//' # Normalize using the Exponential model
//' normalized_y <- data_normalization_1D(rexp(100), type = "exp")
//'
//' @export
// [[Rcpp::export]]
std::vector<double> data_normalization_1D(std::vector<double>& y,
                                          std::string type = "gauss")
{
  int n = y.size();

  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  if (type == "variance")
  {
    double mean_y = std::accumulate(y.begin(), y.end(), 0.0) / n;
    for (int i = 0; i < n; ++i){y[i] = y[i] - mean_y;}
    return y;
  }

  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  if (type == "gauss")
  {
    double sdNoise = sdDiff(y);
    for (int i = 0; i < n; ++i){y[i] = y[i] / sdNoise;}
    return y;
  }

  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  if(type == "poisson")
  {
    for(int i = 0; i < n; i++){if(y[i] < 0){throw std::range_error("negative data not compatible with poisson model");}}
    double mean_y = std::accumulate(y.begin(), y.end(), 0.0) / n;
    for(int i = 0; i < n; ++i){y[i] = y[i] / mean_y;}
    return y;
  }

  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  if(type == "exp")
  {
    for(int i = 0; i < n; i++){if(y[i] < 0){throw std::range_error("negative data not compatible with exp model");}}
    double mean_y = std::accumulate(y.begin(), y.end(), 0.0) / n;
    for(int i = 0; i < n; ++i) {y[i] = y[i] / mean_y;}
    return y;
  }

  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  if(type == "binom")
  {
    for(int i = 0; i < n; i++){if(y[i] < 0){throw std::range_error("negative data not compatible with exp model");}}
    double max_y = *(std::max_element(y.begin(), y.end()));

    // Return the maximum value
    for(int i = 0; i < n; ++i) {y[i] = y[i] / max_y;}
    return y;
  }

  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  if(type == "negbin")
  {
    unsigned int windowSize = 100;
    unsigned int k = y.size() / windowSize;
    double mean = 0;
    double variance = 0;
    double disp = 0;

    for(unsigned int j = 0; j < k; j++)
    {
      mean = 0;
      variance = 0;
      for(unsigned int i = j * windowSize; i < (j + 1)*windowSize; i++){mean = mean + y[i];}
      mean = mean/windowSize;
      for(unsigned int i =  j * windowSize; i < (j + 1)*windowSize; i++){variance = variance + (y[i] - mean) * (y[i] - mean);}
      variance = variance/(windowSize - 1);
      disp = disp  + (mean * mean / (variance - mean));
    }
    disp = disp/k;
    for(size_t i = 0; i < y.size(); i++){y[i] = y[i]/disp;}
    return y;
  }

  if(type == "geom" || type == "bern")
  {
    return y;
  }

  Rcpp::stop("Unsupported type specified.");
}


//' Data Normalization Function
//'
//' @name data_normalization_MD
//'
//' @description
//' Normalizes the input time series data `y` according to the specified `type`.
//' The normalization process depends on the statistical model type, which can be one of the following:
//' "gauss" (Gaussian/normal distribution), "exp" (exponential distribution),
//' "poisson" (Poisson distribution), "geom" (geometric distribution),
//' "bern" (Bernoulli distribution), "binom" (binomial distribution),
//' "negbin" (negative binomial distribution), or "variance"
//'
//' @param y A numeric matrix representing the time series to be normalized and then segmented.
//' @param type A string specifying the model type for normalization.
//' The available options are "gauss", "exp", "poisson", "geom", "bern", "binom", "negbin", "variance".
//' The default is "gauss".
//' @return A numeric matrix that is the normalized version of the input time series `y`, normalized row by row.
//' @examples
//' # Normalize a random time series using the Gaussian model
//' normalized_y <- data_normalization_MD(matrix(rnorm(100), nrow = 2), type = "gauss")
//' @export
// [[Rcpp::export]]
NumericMatrix data_normalization_MD(NumericMatrix& y,
                                std::string type = "gauss")
{
  int nCols = y.ncol();
  int nRows = y.nrow();
  NumericMatrix y_normalized = clone(y);
  std::vector<double> currentRow(nCols);
  std::vector<double> currentRowNEW(nCols);

  // Fill the vector with the first row of the matrix
  for (int i = 0; i < nRows; i++)
  {
    for (int j = 0; j < nCols; j++)
    {
      currentRow[j] = y(i, j);
    }
    currentRowNEW = data_normalization_1D(currentRow, type);
    for (int j = 0; j < nCols; j++)
    {
      y_normalized(i, j) = currentRowNEW[j];
    }
  }
  return y_normalized;
}

