#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<double> PELT(NumericVector data, double penalty = 0)
{
  int n = data.size();
  if(penalty == 0){penalty = 2*log(n);}

  ///
  /// DATA TRANSFORMATION
  ///
  double cumsums[n + 1];
  cumsums[0] = 0;
  for (int i = 0; i < n; i++) cumsums[i + 1] = cumsums[i] + data[i];

  ///
  /// INITIALISATION
  ///
  double costQ[n + 1]; //costQ[i] optimal cost for data y(1) to y(i-1)
  int cp[n + 1]; //cp vector cp[i] = index of the last change-point for data y(1) to y(i)
  costQ[0] = -penalty;
  cp[0] = 0;

  std::vector<double> res;
  for (int t = 0; t < (n+1); t++) {res.push_back(cumsums[t]);}

  return(res);
}
