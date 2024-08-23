#include <Rcpp.h>

#include "1D_A3_ExpModel.h"

using namespace Rcpp;

Exp_1D::Exp_1D(bool use_dual_max, bool random_constraint, Nullable<double> alpha)
  : DUST_1D(use_dual_max, random_constraint, alpha) {}

double Exp_1D::Cost(int t, int s) const
{
  return (t - s)*(1 + log((cumsum[t] - cumsum[s])/(t - s)));
}

double Exp_1D::dualEval(double point, double minCost, int t, int s, int r) const
{
  int objectiveLength = t - s;
  double objectiveMean = (cumsum[t] - cumsum[s]) / objectiveLength; // m_it
  int constraintLength = s - r;
  double constraintMean = (cumsum[s] - cumsum[r]) / constraintLength; // m_ji

  ///
  /// point in the right interval:
  /// TO DO: IMPROVE with exception objectiveMean = 0
  point = point * std::min(1.0, objectiveMean/constraintMean);
  ///
  ///
  double R = (objectiveMean - point * constraintMean) / (1 - point);

  return (costRecord[s] - minCost) / objectiveLength
  + point * (costRecord[s] - costRecord[r]) / constraintLength
  + (1 - point) * (log(R) + 1);
}

double Exp_1D::dualMax(double minCost, int t, int s, int r) const
{
  return - std::numeric_limits<double>::infinity();
}
