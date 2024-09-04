#include <Rcpp.h>
#include <cmath>

#include "1D_A6_BinomModel.h"

using namespace Rcpp;

Binom_1D::Binom_1D(int dual_max, bool random_constraint, Nullable<double> alpha, Nullable<int> nbLoops)
  : DUST_1D(dual_max, random_constraint, alpha, nbLoops) {}

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
  double objectiveMean = (cumsum[t] - cumsum[s]) / (t - s); // m_it
  double constraintMean = (cumsum[s] - cumsum[r]) / (s - r); // m_ji

  ///
  /// point in the right interval:
  if(constraintMean != 0 && constraintMean != 1){point = point * std::min(objectiveMean/constraintMean, (1 - objectiveMean)/(1 - constraintMean));}
  else{
    if(constraintMean == 0){point = point * (1 - objectiveMean);}else{point = point * objectiveMean;}
  }
  ///
  ///
  double R = (objectiveMean - point * constraintMean) / (1 - point);

  return (costRecord[s] - minCost) / (t - s)
  + point * (costRecord[s] - costRecord[r]) / (s - r)
  - (1 - point) * (R * std::log(R) + (1 - R) * std::log(1 - R));
}


double Binom_1D::dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  return (-std::numeric_limits<double>::infinity());
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Binom_1D::Dstar(double x) const
{
  return 0;
}


double Binom_1D::DstarPrime(double x) const
{
  return 0;
}


double Binom_1D::DstarSecond(double x) const
{
  return 0;
}


