#include <Rcpp.h>
#include <algorithm> // for std::min

#include "1D_A2_PoissonModel.h"

using namespace Rcpp;

Poisson_1D::Poisson_1D(bool use_dual_max, bool random_constraint, Nullable<double> alpha)
  : DUST_1D(use_dual_max, random_constraint, alpha) {}

double Poisson_1D::Cost(int t, int s) const
{
  return (cumsum[t] - cumsum[s])*(1 - log((cumsum[t] - cumsum[s])/(t - s)));
}

double Poisson_1D::dualEval(double point, double minCost, int t, int s, int r) const
{
  int objectiveLength = t - s;
  double objectiveMean = (cumsum[t] - cumsum[s]) / objectiveLength; // m_it
  int constraintLength = s - r;
  double constraintMean = (cumsum[s] - cumsum[r]) / constraintLength; // m_ji

  ///
  /// point in the right interval:
  /// TO DO: IMPROVE with exception objectiveMean = 0
  point = point * std::min(1.0, constraintMean/objectiveMean);
  ///
  ///

  return (costRecord[s] - minCost) / objectiveLength
    + point * (costRecord[s] - costRecord[r]) / constraintLength
    -(constraintMean - point * objectiveMean)* (log((constraintMean - point * objectiveMean) / (1 - point)) - 1);
}

double Poisson_1D::dualMax(double minCost, int t, int s, int r) const
{
  return - std::numeric_limits<double>::infinity();
}
