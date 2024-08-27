#include <Rcpp.h>
#include <algorithm> // for std::min
#include <cmath>

#include <random>

#include "1D_A2_PoissonModel.h"

using namespace Rcpp;

Poisson_1D::Poisson_1D(bool use_dual_max, bool random_constraint, Nullable<double> alpha)
  : DUST_1D(use_dual_max, random_constraint, alpha) {}

double Poisson_1D::Cost(unsigned int t, unsigned int s) const
{
  double res = 0;
  if(cumsum[t] - cumsum[s] != 0)
    {res = (cumsum[t] - cumsum[s])*(1 - std::log((cumsum[t] - cumsum[s])/double(t - s)));}
  return res;
}

double Poisson_1D::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  double objectiveMean = (cumsum[t] - cumsum[s]) / (t - s); // m_it
  double constraintMean = (cumsum[s] - cumsum[r]) / (s - r); // m_ji

  ///
  /// point in the right interval:
  /// TO DO: IMPROVE with exception objectiveMean = 0
  point = point * std::min(1.0, constraintMean/objectiveMean);
  ///
  ///

  return (costRecord[s] - minCost) / (t - s)
    + point * (costRecord[s] - costRecord[r]) / (s - r)
    - (objectiveMean - point * constraintMean) * (std::log((objectiveMean - point * constraintMean) / (1 - point)) - 1);
}

double Poisson_1D::dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  double max_val = Poisson_1D::dualEval(0.4, minCost, t, s, r);
  double max_val2 = Poisson_1D::dualEval(0.6, minCost, t, s, r);

  if (max_val2 > max_val)
  {
    max_val = max_val2;
    double max_val3 = Poisson_1D::dualEval(0.8, minCost, t, s, r);
    if (max_val3 > max_val){max_val = max_val3;}
  }
  else
  {
    double max_val3 = Poisson_1D::dualEval(0.2, minCost, t, s, r);
    if (max_val3 > max_val){max_val = max_val3;}
  }
  return max_val;
}
