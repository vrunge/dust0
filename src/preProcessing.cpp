#include <Rcpp.h>

using namespace Rcpp;

#include <Rcpp.h>
#include <cmath>

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
double sdDiff(Rcpp::NumericVector& y, std::string method = "HALL")
{
  ///////////////////////  HALL
  ///////////////////////  HALL
  ///////////////////////  HALL
  if(method == "HALL")
  {
    if(!Rcpp::is_true(Rcpp::all(Rcpp::is_finite(y))) || y.size() < 5)
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
    if(!Rcpp::is_true(Rcpp::all(Rcpp::is_finite(y))) || y.size() < 2)
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
    if(!Rcpp::is_true(Rcpp::all(Rcpp::is_finite(y))) || y.size() < 2)
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


