#include <Rcpp.h>

#include <cmath>

#include "1D_A5_BernModel.h"

using namespace Rcpp;

Bern_1D::Bern_1D(bool use_dual_max, bool random_constraint, Nullable<double> alpha)
  : DUST_1D(use_dual_max, random_constraint, alpha) {}

double Bern_1D::Cost(unsigned int t, unsigned int s) const
{
  double res = 0;
  double m = (cumsum[t] - cumsum[s]) / (t - s);
  if(m != 0 && m != 1)
  {res = - double(t - s) * (m * std::log(m) + (1 - m) *  std::log(1 - m));}
  return res;
}

double Bern_1D::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  return -std::numeric_limits<double>::infinity();
  double objectiveMean = (cumsum[t] - cumsum[s]) / (t - s); // m_it
  double constraintMean = (cumsum[s] - cumsum[r]) / (s - r); // m_ji

  ///
  /// point in the right interval:
  /// TO DO: IMPROVE with exception objectiveMean = 0
  point = point * std::min(objectiveMean/constraintMean, (1 - objectiveMean)/(1 - constraintMean));
  if(constraintMean == 0 || constraintMean == 1){point = 0;}
  ///
  ///
  double R = (objectiveMean - point * constraintMean) / (1 - point);

  return (costRecord[s] - minCost) / (t - s)
  + point * (costRecord[s] - costRecord[r]) / (s - r)
  - (1 - point) * (R * log(R) + (1 - R) * log(1 - R));
}

double Bern_1D::dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  double max_val = Bern_1D::dualEval(0.4, minCost, t, s, r);
  double max_val2 = Bern_1D::dualEval(0.6, minCost, t, s, r);

  if (max_val2 > max_val)
  {
    max_val = max_val2;
    double max_val3 = Bern_1D::dualEval(0.8, minCost, t, s, r);
    if (max_val3 > max_val){max_val = max_val3;}
  }
  else
  {
    double max_val3 = Bern_1D::dualEval(0.2, minCost, t, s, r);
    if (max_val3 > max_val){max_val = max_val3;}
  }
  return max_val;
}

