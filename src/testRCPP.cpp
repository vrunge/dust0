#include <Rcpp.h>
using namespace Rcpp;


//// x sans copie +  std::vector pour CUMSUM

// [[Rcpp::export]]
double cs1(NumericVector& x)
{
  std::vector<double> cumsum;
  cumsum = std::vector<double>(x.size(), 0.);

  cumsum[0] = x[0];
  for (int i = 1; i < cumsum.size(); ++i)
  {
    cumsum[i] = cumsum[i-1] + x[i];
  }
  return 0;
}


//// x copie dans vector data +  std::vector pour CUMSUM

// [[Rcpp::export]]
double cs2(NumericVector& x)
{
  std::vector<double> data(x.begin(), x.end());

  std::vector<double> cumsum;
  cumsum = std::vector<double>(x.size(), 0.);

  cumsum[0] = x[0];
  for (int i = 1; i < cumsum.size(); ++i)
  {
    cumsum[i] = cumsum[i-1] + data[i];
  }
  return 0;
}


//// x sans copie  + NumericVector pour CUMSUM

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




//// x avec copie dans NumericVector + NumericVector pour CUMSUM

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

//// x avec copie dans NumericVector avec move + NumericVector pour CUMSUM

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


