#include <Rcpp.h>
#include <algorithm> // for std::min

#include "1D_A2_PoissonModel.h"

using namespace Rcpp;

Poisson_1D::Poisson_1D(bool use_dual_max, bool random_constraint, Nullable<double> alpha)
  : DUST_1D(use_dual_max, random_constraint, alpha) {}

double Poisson_1D::Cost(int t, int s) const
{
  double res = 0;
  if(cumsum[t] - cumsum[s] != 0)
    {res = (cumsum[t] - cumsum[s])*(1 - log((cumsum[t] - cumsum[s])/(t - s)));}
  return res;
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
    -(objectiveMean - point * constraintMean) * (log((objectiveMean - point * constraintMean) / (1 - point)) - 1);
}

double Poisson_1D::dualMax(double minCost, int t, int s, int r) const
{
  return - std::numeric_limits<double>::infinity();
}
