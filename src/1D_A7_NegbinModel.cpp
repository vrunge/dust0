#include <Rcpp.h>

#include "1D_A7_NegbinModel.h"

using namespace Rcpp;

Negbin_1D::Negbin_1D(bool use_dual_max, bool random_constraint, Nullable<double> alpha)
  : DUST_1D(use_dual_max, random_constraint, alpha) {}

double Negbin_1D::Cost(unsigned int t, unsigned int s) const
{
  double res = 0;
  double m = (cumsum[t] - cumsum[s])/(t - s);
  if(m != 0)
  {res = (t - s) * log(1 + m) - (cumsum[t] - cumsum[s]) * log(m / (1 + m));}
  return res;
}

double Negbin_1D::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  unsigned int objectiveLength = t - s;
  double objectiveMean = (cumsum[t] - cumsum[s]) / objectiveLength; // m_it
  unsigned int constraintLength = s - r;
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
  + (1 - point) * (R * log(R) - (1 + R) * log(1 + R));
}

double Negbin_1D::dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  return - std::numeric_limits<double>::infinity();
}

