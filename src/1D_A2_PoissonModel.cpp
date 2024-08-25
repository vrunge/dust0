#include <Rcpp.h>
#include <algorithm> // for std::min

#include <random>

#include "1D_A2_PoissonModel.h"

using namespace Rcpp;

Poisson_1D::Poisson_1D(bool use_dual_max, bool random_constraint, Nullable<double> alpha)
  : DUST_1D(use_dual_max, random_constraint, alpha) {}

double Poisson_1D::Cost(unsigned int t, unsigned int s) const
{
  double res = 0;
  if(cumsum[t] - cumsum[s] != 0)
    {res = (cumsum[t] - cumsum[s])*(1 - log((cumsum[t] - cumsum[s])/(t - s)));}
  return res;
}

double Poisson_1D::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  unsigned int objectiveLength = t - s;
  double objectiveMean = (cumsum[t] - cumsum[s]) / objectiveLength; // m_it
  unsigned int constraintLength = s - r;
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

double Poisson_1D::dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const
{

  double max_val = Poisson_1D::dualEval(0.3, minCost, t, s, r);
  double max_val2 = Poisson_1D::dualEval(0.5, minCost, t, s, r);
  double max_val3 = Poisson_1D::dualEval(0.7, minCost, t, s, r);

  if (max_val2 > max_val){max_val = max_val2;}
  if (max_val3 > max_val){max_val = max_val3;}
  //double current_val = 0;
  //for (int i = 1; i < 3; ++i)
    //{
    //current_val =  Poisson_1D::dualEval(i/10, minCost, t, s, r);
    //  if (current_val > max_val){max_val = current_val;}
  //}
  return max_val;
  //return - std::numeric_limits<double>::infinity();
}
