#include <Rcpp.h>
#include <cmath>

#include "1D_A6_BinomModel.h"

using namespace Rcpp;

Binom_1D::Binom_1D(bool use_dual_max, bool random_constraint, Nullable<double> alpha)
  : DUST_1D(use_dual_max, random_constraint, alpha) {}

double Binom_1D::Cost(unsigned int t, unsigned int s) const
{
  double res = 0;
  double m = (cumsum[t] - cumsum[s]) / (t - s);
  if(m != 0 && m != 1)
  {res = - double(t - s) * (m * std::log(m) + (1 - m) *  std::log(1 - m));}
  return res;
}

double Binom_1D::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  unsigned int objectiveLength = t - s;
  double objectiveMean = (cumsum[t] - cumsum[s]) / objectiveLength; // m_it
  unsigned int constraintLength = s - r;
  double constraintMean = (cumsum[s] - cumsum[r]) / constraintLength; // m_ji

  ///
  /// point in the right interval:
  /// TO DO: IMPROVE with exception objectiveMean = 0
  point = point * std::min(objectiveMean/constraintMean, (objectiveMean - 1)/(constraintMean - 1));
  ///
  ///
  double R = (objectiveMean - point * constraintMean) / (1 - point);

  return (costRecord[s] - minCost) / objectiveLength
  + point * (costRecord[s] - costRecord[r]) / constraintLength
  + (1 - point) * ((1 - R) * std::log(1 - R) + R * std::log(R));
}

double Binom_1D::dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  double max_val = Binom_1D::dualEval(0.4, minCost, t, s, r);
  double max_val2 = Binom_1D::dualEval(0.6, minCost, t, s, r);

  if (max_val2 > max_val)
  {
    max_val = max_val2;
    double max_val3 = Binom_1D::dualEval(0.8, minCost, t, s, r);
    if (max_val3 > max_val){max_val = max_val3;}
  }
  else
  {
    double max_val3 = Binom_1D::dualEval(0.2, minCost, t, s, r);
    if (max_val3 > max_val){max_val = max_val3;}
  }
  return max_val;
}

