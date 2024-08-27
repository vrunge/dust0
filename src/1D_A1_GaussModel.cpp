#include <Rcpp.h>
#include <cmath>

#include "1D_A1_GaussModel.h"

using namespace Rcpp;

Gauss_1D::Gauss_1D(bool use_dual_max, bool random_constraint, Nullable<double> alpha)
  : DUST_1D(use_dual_max, random_constraint, alpha) {}

double Gauss_1D::Cost(unsigned int t, unsigned int s) const
{
  return - 0.5 * (cumsum[t] - cumsum[s]) * (cumsum[t] - cumsum[s]) / (t - s);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Gauss_1D::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  return (costRecord[s] - minCost) / (t - s)
    + point * (costRecord[s] - costRecord[r]) / (s - r)
    - 0.5 * pow((cumsum[t] - cumsum[s]) / (t - s) - point * ((cumsum[s] - cumsum[r]) / (s - r)),2) / (1 - point);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Gauss_1D::dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  // Compute the optimal point on which to evaluate the duality function

  double A = (cumsum[t] - cumsum[s])/ (t - s); // m_it
  double B = (cumsum[s] - cumsum[r])/ (s - r); // m_ji
  double AmB = std::abs(A - B);
  double B2p2C = std::sqrt(B*B + 2*(costRecord[s] - costRecord[r])/ (s - r));

  // Case 1: mu* = 0
  // deduce the following condition from the formula for mu*
  if (AmB >= B2p2C)
    return (costRecord[s] - minCost) / (t - s)  - 0.5 * A*A;

  // Case 2: mu* > 0
    return (costRecord[s] - minCost) / (t - s) + 0.5 * (AmB - B2p2C)*(AmB - B2p2C) - 0.5*A*A;
}


