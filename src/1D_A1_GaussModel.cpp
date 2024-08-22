#include <Rcpp.h>

#include "1D_A1_GaussModel.h"

using namespace Rcpp;

Gauss_1D::Gauss_1D(bool use_dual_max, bool random_constraint, Nullable<double> alpha)
  : DUST_1D(use_dual_max, random_constraint, alpha) {}

double Gauss_1D::Cost(int t, int s) const
{
  return - pow(cumsum[t] - cumsum[s], 2) / (t - s);
}

double Gauss_1D::dualEval(double point, double minCost, int t, int s, int r) const
{
  int objectiveLength = t - s;
  double objectiveMean = (cumsum[t] - cumsum[s]) / objectiveLength; // m_it
  
  int constraintLength = s - r;
  double constraintMean = (cumsum[s] - cumsum[r]) / constraintLength; // m_ji
  
  return (costRecord[s] - minCost) / objectiveLength - pow(objectiveMean, 2)
    + point * ((costRecord[s] - costRecord[r]) / constraintLength + pow(constraintMean, 2))
    - pow(constraintMean - point * objectiveMean, 2) / (1 - point);
}

double Gauss_1D::dualMax(double minCost, int t, int s, int r) const
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
  
  // Rcout << "t = " << t << "; s = " << s << "; r = " << r << std::endl;
  
  int objectiveLength = t - s;
  double objectiveMean = (cumsum[t] - cumsum[s]) / objectiveLength; // m_it
  
  int constraintLength = s - r;
  double constraintMean = (cumsum[s] - cumsum[r]) / constraintLength; // m_ji
  double sqGapMean = pow(constraintMean, 2);
  
  double costJI = (costRecord[s] - costRecord[r]) / constraintLength + sqGapMean; // Qi - Qj / s - r
  double costIT = (minCost - costRecord[s]) / objectiveLength + pow(objectiveMean, 2);
  double meanGap = fabs(objectiveMean - constraintMean);
  
  // Case 1: mu* = 0
  // deduce the following condition from the formula for mu*
  if (costJI <= pow(meanGap, 2))
    return - costIT;
  
  // Case 2: mu* > 0
  return - costIT + pow(meanGap - sqrt(costJI), 2);
}
