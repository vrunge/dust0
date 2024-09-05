#include <Rcpp.h>
#include <cmath>

#include "1D_A3_ExpModel.h"

using namespace Rcpp;

Exp_1D::Exp_1D(int dual_max, bool random_constraint, Nullable<double> alpha, Nullable<int> nbLoops)
  : DUST_1D(dual_max, random_constraint, alpha, nbLoops) {}

double Exp_1D::Cost(unsigned int t, unsigned int s) const
{
  double delta_t = t - s;
  double diff_cumsum = cumsum[t] - cumsum[s];
  return delta_t * (1.0 + std::log(diff_cumsum / delta_t));
}

double Exp_1D::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  double objectiveMean = (cumsum[t] - cumsum[s]) / (t - s); // m_it
  double constraintMean = (cumsum[s] - cumsum[r]) / (s - r); // m_ji

  ///
  /// point in the right interval:
  if(constraintMean != 0){point = point * std::min(1.0, objectiveMean/constraintMean);}
  ///
  ///

  return (costRecord[s] - minCost) / (t - s)
  + point * (costRecord[s] - costRecord[r]) / (s - r)
  + (1 - point) * (std::log((objectiveMean - point * constraintMean) / (1 - point)) + 1);
}


double Exp_1D::dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  return (-std::numeric_limits<double>::infinity());
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Exp_1D::muMax(double a, double b) const
{
  double res = 1;
  if(b != 0){res = std::min(1.0, a/b);}
  return res;
}

double Exp_1D::Dstar(double x) const
{
  return (-std::log(x) - 1.0);
}


double Exp_1D::DstarPrime(double x) const
{
  return -1.0/x;
}

double Exp_1D::DstarSecond(double x) const
{
  return 1.0/std::pow(x,2);
}


