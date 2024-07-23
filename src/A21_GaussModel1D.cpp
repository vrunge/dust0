// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "A21_GaussModel1D.h"

using namespace Rcpp;

double GaussModel1D::modelCost(int t, int i)
{
  return - pow(valuesCumsum[t] - valuesCumsum[i], 2) / (t - i);
}

double GaussModel1D::modelDual(double point, int t, int i, int j, double testThreshold)
{
  int objectiveLength = t - i;
  double objectiveMean = (valuesCumsum[t] - valuesCumsum[i]) / objectiveLength; // m_it
  
  int constraintLength = i - j;
  double constraintMean = (valuesCumsum[i] - valuesCumsum[j]) / constraintLength; // m_ji
  
  return (costRecord[i] - testThreshold) / objectiveLength - pow(objectiveMean, 2)
    + point * ((costRecord[i] - costRecord[j]) / constraintLength + pow(constraintMean, 2))
    - pow(constraintMean - point * objectiveMean, 2) / (1 - point);
}

double GaussModel1D::modelTest(int t, int i, int j, double testThreshold)
{
  
  // Compute the optimal point on which to evaluate the duality function
  //
  // Denoting y_it = y[(i+1):t]; y_ji = y[(j+1):i]
  // Denoting m_it = mean(y_it); m_ji = mean(y_ji)
  // Duality function: D(mu) = Qi + (t - i) ( mu (Qi-Qj)/(i-j) - (m_it - mu m_ji)/(1-mu))
  //
  // Formula: mu* = max(0, 1 - abs(m_it - m_ji)/sqrt((Qi-Qj)/(i-j) + m_ji^2))
  // Formula: mu* > 0, d* = D(mu*) = Qi - (t-i) m_it^2 + (t-i) (sqrt((Qi-Qj)/(i-j) + m_ji^2) - abs(m_it - m_ji))^2
  // Formula: (d* - Qt) / (t - i) = (Qi-Qt)/(t-i) - m_it^2 + (sqrt((Qi-Qj)/(i-j) + m_ji^2) - abs(m_it - m_ji))^2
  // Pruning happens if (d* - Qt) / (t - i) > 0
  
  int objectiveLength = t - i;
  double objectiveMean = (valuesCumsum[t] - valuesCumsum[i]) / objectiveLength; // m_it
  
  int constraintLength = i - j;
  double constraintMean = (valuesCumsum[i] - valuesCumsum[j]) / constraintLength; // m_ji
  double sqGapMean = pow(constraintMean, 2);
  
  double costJI = (costRecord[i] - costRecord[j]) / constraintLength + sqGapMean; // Qi - Qj / i - j
  double costIT = (testThreshold - costRecord[i]) / objectiveLength + pow(objectiveMean, 2);
  double meanGap = fabs(objectiveMean - constraintMean);
  
  // Case 1: mu* = 0
  // deduce the following condition from the formula for mu*
  if (costJI <= pow(meanGap, 2))
    return - costIT;
  
  // Case 2: mu* > 0
  return - costIT + pow(meanGap - sqrt(costJI), 2);
}
