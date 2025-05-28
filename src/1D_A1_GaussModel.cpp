#include <Rcpp.h>
#include <cmath>

#include "1D_A1_GaussModel.h"

using namespace Rcpp;

Gauss_1D::Gauss_1D(int dual_max_type, int constraints_type, Nullable<int> nbLoops)
  : DUST_1D(dual_max_type, constraints_type, nbLoops) {}

double Gauss_1D::Cost(unsigned int t, unsigned int s) const
{
  return - 0.5 * (cumsum[t] - cumsum[s]) * (cumsum[t] - cumsum[s]) / (t - s);
}

double Gauss_1D::statistic(double& data) const
{return(data);}




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Gauss_1D::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  return (costRecord[s] - minCost) / (t - s)
    + point * (costRecord[s] - costRecord[r]) / (s - r)
    - 0.5 * std::pow((cumsum[t] - cumsum[s]) / (t - s) - point * ((cumsum[s] - cumsum[r]) / (s - r)), 2) / (1 - point);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Gauss_1D::dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  // Compute the optimal point on which to evaluate the duality function

  double A = (cumsum[t] - cumsum[s])/ (t - s); // m_it
  double B = (cumsum[s] - cumsum[r])/ (s - r); // m_ji
  double absAmB = std::abs(A - B);
  double sqrtB2p2C = std::sqrt(B*B + 2*(costRecord[s] - costRecord[r])/ (s - r));

  // Case 1: mu* = 0
  // deduce the following condition from the formula for mu*
  if (absAmB >= sqrtB2p2C)
    return (costRecord[s] - minCost) / (t - s)  - 0.5 * A*A;

  // Case 2: mu* > 0
    return (costRecord[s] - minCost) / (t - s) - 0.5*A*A + 0.5 * (absAmB - sqrtB2p2C)*(absAmB - sqrtB2p2C);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Gauss_1D::muMax(double a, double b) const
{
  return 1;
}

bool Gauss_1D::isBoundary(double a) const
{
  return false;
}



double Gauss_1D::Dstar(double x) const
{
  return 0.5 * x * x;
}

double Gauss_1D::DstarPrime(double x) const
{
  return x;
}

double Gauss_1D::DstarPrimeInv(double x) const
{
  return x;
}

double Gauss_1D::DstarSecond(double x) const
{
  return 1.0;
}

std::string Gauss_1D::get_model() const { return "gauss"; }
