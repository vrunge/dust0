#include <Rcpp.h>
using namespace Rcpp;


//' Cumulative Sum (cs1): No Copy of Input Vector
//'
//' Computes the cumulative sum of the input vector \code{x} without copying the data into a new container.
//'
//' @param x A numeric vector.
//' @return A double value (currently always returns 0, but the cumulative sum is computed).
//' @examples
//' \dontrun{
//'   x <- c(1, 2, 3, 4, 5)
//'   cs1(x)  # Returns 0, but the cumulative sum is computed internally
//' }
//' @export
// [[Rcpp::export]]
double cs1(NumericVector& x)
{
  std::vector<double> cumsum;
  cumsum = std::vector<double>(x.size(), 0.);

  cumsum[0] = x[0];
  for (size_t i = 1; i < cumsum.size(); ++i)
  {
    cumsum[i] = cumsum[i-1] + x[i];
  }
  return 0;
}

//' Cumulative Sum (cs2): Copy Input Vector to std::vector
//'
//' Computes the cumulative sum of the input vector \code{x} after copying it into a \code{std::vector}.
//'
//' @param x A numeric vector.
//' @return A double value (currently always returns 0, but the cumulative sum is computed internally).
//' @examples
//' \dontrun{
//'   x <- c(1, 2, 3, 4, 5)
//'   cs2(x)  # Returns 0, but the cumulative sum is computed internally
//' }
//' @export
// [[Rcpp::export]]
double cs2(NumericVector& x)
{
  std::vector<double> data(x.begin(), x.end());

  std::vector<double> cumsum;
  cumsum = std::vector<double>(x.size(), 0.);

  cumsum[0] = x[0];
  for (size_t i = 1; i < cumsum.size(); ++i)
  {
    cumsum[i] = cumsum[i-1] + data[i];
  }
  return 0;
}

//' Cumulative Sum (cs3): No Copy of Input Vector to NumericVector
//'
//' Computes the cumulative sum of the input vector \code{x} without copying the data into a new container.
//'
//' @param x A numeric vector.
//' @return A double value (currently always returns 0, but the cumulative sum is computed internally).
//' @examples
//' \dontrun{
//'   x <- c(1, 2, 3, 4, 5)
//'   cs3(x)  # Returns 0, but the cumulative sum is computed internally
//' }
//' @export
// [[Rcpp::export]]
double cs3(NumericVector& x)
{
  NumericVector cumsum;
  cumsum = NumericVector(x.size(), 0.);

  cumsum[0] = x[0];
  for (int i = 1; i < cumsum.size(); ++i)
  {
    cumsum[i] = cumsum[i-1] + x[i];
  }
  return 0;
}

//' Cumulative Sum (cs4): Copy Input Vector to NumericVector
//'
//' Computes the cumulative sum of the input vector \code{x} after copying it into a \code{NumericVector}.
//'
//' @param x A numeric vector.
//' @return A double value (currently always returns 0, but the cumulative sum is computed internally).
//' @examples
//' \dontrun{
//'   x <- c(1, 2, 3, 4, 5)
//'   cs4(x)  # Returns 0, but the cumulative sum is computed internally
//' }
//' @export
// [[Rcpp::export]]
double cs4(NumericVector& x)
{
  NumericVector data(x.begin(), x.end());

  NumericVector cumsum;
  cumsum = NumericVector(x.size(), 0.);

  cumsum[0] = x[0];
  for (int i = 1; i < cumsum.size(); ++i)
  {
    cumsum[i] = cumsum[i-1] + data[i];
  }
  return 0;
}

//' Cumulative Sum (cs5): Move Input Vector to NumericVector
//'
//' Computes the cumulative sum of the input vector \code{x} after moving it into a \code{NumericVector}.
//'
//' @param x A numeric vector.
//' @return A double value (currently always returns 0, but the cumulative sum is computed internally).
//' @examples
//' \dontrun{
//'   x <- c(1, 2, 3, 4, 5)
//'   cs5(x)  # Returns 0, but the cumulative sum is computed internally
//' }
//' @export
// [[Rcpp::export]]
double cs5(NumericVector& x)
{
  NumericVector data = std::move(x);

  NumericVector cumsum;
  cumsum = NumericVector(x.size(), 0.);

  cumsum[0] = x[0];
  for (int i = 1; i < cumsum.size(); ++i)
  {
    cumsum[i] = cumsum[i-1] + data[i];
  }
  return 0;
}

