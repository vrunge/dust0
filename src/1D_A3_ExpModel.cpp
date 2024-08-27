#include <Rcpp.h>
#include <cmath>

#include "1D_A3_ExpModel.h"

using namespace Rcpp;

Exp_1D::Exp_1D(bool use_dual_max, bool random_constraint, Nullable<double> alpha)
  : DUST_1D(use_dual_max, random_constraint, alpha) {}

double Exp_1D::Cost(unsigned int t, unsigned int s) const
{
  return (t - s)*(1 + std::log((cumsum[t] - cumsum[s])/(t - s)));
}

double Exp_1D::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  double objectiveMean = (cumsum[t] - cumsum[s]) / (t - s); // m_it
  double constraintMean = (cumsum[s] - cumsum[r]) / (s - r); // m_ji

  ///
  /// point in the right interval:
  /// TO DO: IMPROVE with exception objectiveMean = 0
  point = point * std::min(1.0, objectiveMean/constraintMean);
  ///
  ///

  return (costRecord[s] - minCost) / (t - s)
  + point * (costRecord[s] - costRecord[r]) / (s - r)
  + (1 - point) * (std::log((objectiveMean - point * constraintMean) / (1 - point) + 1));
}

double Exp_1D::dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  double max_val = Exp_1D::dualEval(0.4, minCost, t, s, r);
  double max_val2 = Exp_1D::dualEval(0.6, minCost, t, s, r);

  if (max_val2 > max_val)
  {
    max_val = max_val2;
    double max_val3 = Exp_1D::dualEval(0.8, minCost, t, s, r);
    if (max_val3 > max_val){max_val = max_val3;}
  }
  else
  {
    double max_val3 = Exp_1D::dualEval(0.2, minCost, t, s, r);
    if (max_val3 > max_val){max_val = max_val3;}
  }
  return max_val;
}
