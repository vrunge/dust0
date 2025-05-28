#include <Rcpp.h>
#include <cmath>

#include "1D_A8_VarianceModel.h"

using namespace Rcpp;

Variance_1D::Variance_1D(int dual_max_type, int constraints_type, Nullable<int> nbLoops)
  : DUST_1D(dual_max_type, constraints_type, nbLoops) {}

double Variance_1D::Cost(unsigned int t, unsigned int s) const
{
  double delta_t = 1.0*(t - s);
  double diff_cumsum = cumsum[t] - cumsum[s];
  if(diff_cumsum <= 0){diff_cumsum = 1e-100;} /// choice  1e-100 to avoid -Inf
  return 0.5 * delta_t * (1.0 + std::log(std::abs(diff_cumsum / delta_t)));
}

double Variance_1D::statistic(double& data) const
{return(data * data);}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Variance_1D::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  double objectiveMean = (cumsum[t] - cumsum[s]) / (t - s); // m_it
  double constraintMean = (cumsum[s] - cumsum[r]) / (s - r); // m_ji

  ///
  /// point in the right interval:
  if(constraintMean != 0){point = point * std::min(1.0, objectiveMean/constraintMean);}
  ///
  ///

  return (costRecord[s] - minCost) / (t - s)
  + point * (costRecord[s] - costRecord[r]) / (s - r)
  + (1 - point) * 0.5 * (std::log((objectiveMean - point * constraintMean) / (1 - point)) + 1);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Variance_1D::dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  return (-std::numeric_limits<double>::infinity());
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Variance_1D::muMax(double a, double b) const
{
  double res = 1;
  if(b != 0){res = std::min(1.0, a/b);}
  return res;
}


bool Variance_1D::isBoundary(double a) const
{
  return false;
}

double Variance_1D::Dstar(double x) const
{
  return -0.5 * (std::log(x) + 1.0);
}


double Variance_1D::DstarPrime(double x) const
{
  return -0.5 /x;
}

double Variance_1D::DstarPrimeInv(double x) const
{
  return -0.5 /x;
}

double Variance_1D::DstarSecond(double x) const
{
  return 0.5/std::pow(x,2);
}

std::string Variance_1D::get_model() const { return "variance"; }
