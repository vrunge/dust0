#include <Rcpp.h>

#include "1D_A1_GaussModel.h"

using namespace Rcpp;

Gauss_1D::Gauss_1D(bool use_dual_max, bool random_constraint, Nullable<double> alpha)
  : DUST_1D(use_dual_max, random_constraint, alpha) {}

double Gauss_1D::Cost(unsigned int t, unsigned int s) const
{
  return - (cumsum[t] - cumsum[s]) * (cumsum[t] - cumsum[s]) / (t - s);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Gauss_1D::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  return -(minCost - costRecord[s]) / (t - s)
    + point * (costRecord[s] - costRecord[r]) / (s - r)
    - pow((cumsum[t] - cumsum[s]) / (t - s) - point * (cumsum[s] - cumsum[r]) / (s - r), 2) / (1 - point);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Gauss_1D::dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  // Compute the optimal point on which to evaluate the duality function
  //
  // Denoting y_it = y[(s+1):t]; y_ji = y[(r+1):s]
  // Denoting m_it = mean(y_it); m_ji = mean(y_ji)
  // Duality function: D(mu) = Qi + (t - s) ( mu (Qi-Qj)/(s-r) - (m_it - mu m_ji)/(1-mu))
  //
  // Formula: mu* = max(0, 1 - abs(m_it - m_ji)/sqrt((Qi-Qj)/(s-r) + m_ji^2))
  // Formula: mu* > 0, d* = D(mu*) = Qi - (t-s) m_it^2 + (t-s) (sqrt((Qi-Qj)/(s-r) + m_ji^2) - abs(m_it - m_ji))^2
  // Formula: (d* - Qt) / (t - s) = (Qi-Qt)/(t-s) - m_it^2 + (sqrt((Qi-Qj)/(s-r) + m_ji^2) - abs(m_it - m_ji))^2
  // Pruning happens if (d* - Qt) / (t - s) > 0

  double objectiveMean = (cumsum[t] - cumsum[s]) / (t - s); // m_it
  double constraintMean = (cumsum[s] - cumsum[r]) / (s - r); // m_ji

  double costJI = (costRecord[s] - costRecord[r]) / (s - r) + constraintMean*constraintMean; // Qi - Qj / s - r
  double costIT = (minCost - costRecord[s]) / (t - s) + pow(objectiveMean, 2);
  double meanGap = fabs(objectiveMean - constraintMean);

  // Case 1: mu* = 0
  // deduce the following condition from the formula for mu*
  if (costJI <= pow(meanGap, 2))
    return - costIT;

  // Case 2: mu* > 0
  return - costIT + pow(meanGap - sqrt(costJI), 2);
}
