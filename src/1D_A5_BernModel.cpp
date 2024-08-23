#include <Rcpp.h>

#include "1D_A5_BernModel.h"

using namespace Rcpp;

Bern_1D::Bern_1D(bool use_dual_max, bool random_constraint, Nullable<double> alpha)
  : DUST_1D(use_dual_max, random_constraint, alpha) {}

double Bern_1D::Cost(int t, int s) const
{
  double m = (cumsum[t] - cumsum[s])/(t - s);
  return -(t - s) * log(1 - m) - (cumsum[t] - cumsum[s]) * log(m / (1 - m));
}

double Bern_1D::dualEval(double point, double minCost, int t, int s, int r) const
{
  int objectiveLength = t - s;
  double objectiveMean = (cumsum[t] - cumsum[s]) / objectiveLength; // m_it
  int constraintLength = s - r;
  double constraintMean = (cumsum[s] - cumsum[r]) / constraintLength; // m_ji

  ///
  /// point in the right interval:
  /// TO DO: IMPROVE with exception objectiveMean = 0
  point = point * std::min(constraintMean/objectiveMean, (constraintMean - 1)/(objectiveMean - 1));
  ///
  ///
  double R = (constraintMean - point * objectiveMean) / (1 - point);

  return (costRecord[s] - minCost) / objectiveLength
  + point * (costRecord[s] - costRecord[r]) / constraintLength
  + (1 - point) * ((1 - R) * log(1 - R) + R * log(R));
}

double Bern_1D::dualMax(double minCost, int t, int s, int r) const
{
  return - std::numeric_limits<double>::infinity();
}

