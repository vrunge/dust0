#include <Rcpp.h>
#include <algorithm> // for std::min
#include <cmath>

#include <random>

#include "1D_A2_PoissonModel.h"

using namespace Rcpp;

Poisson_1D::Poisson_1D(int dual_max, bool random_constraint, Nullable<double> alpha, Nullable<int> nbLoops)
  : DUST_1D(dual_max, random_constraint, alpha, nbLoops) {}

double Poisson_1D::Cost(unsigned int t, unsigned int s) const
{
  double res = 0.0;
  double diff = cumsum[t] - cumsum[s];
  if(diff != 0.0) {
    double inv_diff = 1.0 / (t - s);
    double log_term = std::log(diff * inv_diff);
    res = diff * (1.0 - log_term);
  }
  return res;
}

double Poisson_1D::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  double objectiveMean = (cumsum[t] - cumsum[s]) / (t - s); // m_it
  double constraintMean = (cumsum[s] - cumsum[r]) / (s - r); // m_ji

  ///
  /// point in the right interval
  if(constraintMean != 0){point = point * std::min(1.0, objectiveMean/constraintMean);}
  ///

  return (costRecord[s] - minCost) / (t - s)
    + point * (costRecord[s] - costRecord[r]) / (s - r)
    - (objectiveMean - point * constraintMean) * (std::log((objectiveMean - point * constraintMean) / (1 - point)) - 1);
}


double Poisson_1D::dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  return (-std::numeric_limits<double>::infinity());
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Poisson_1D::muMax(double a, double b) const
{
  double res = 1;
  if(b != 0){res = std::min(1.0, a/b);}
  return res;
}

double Poisson_1D::Dstar(double x) const
{
  return (x * (std::log(x) - 1.0));
}


double Poisson_1D::DstarPrime(double x) const
{
  return std::log(x);
}

double Poisson_1D::DstarSecond(double x) const
{
  return (1.0/x);
}




