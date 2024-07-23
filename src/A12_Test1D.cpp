// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "A12_Test1D.h"

using namespace Rcpp;

// --- // Evaluates dual function at random point // --- //
double RandomTest1D::test(int t, int i, int j, double testThreshold)
{
  return modelDual(runif(1)[0], t, i, j, testThreshold);
}

// --- // Returns exact maximum value of dual function // --- //
double DeterministicTest1D::test(int t, int i, int j, double testThreshold)
{
  return modelTest(t, i, j, testThreshold);
}