#include <Rcpp.h>
#include <algorithm> // for std::min
#include <cmath>

#include <random>

#include "1D_A2_PoissonModel.h"

using namespace Rcpp;

Poisson_1D::Poisson_1D(int dual_max_type, int constraints_type, Nullable<int> nbLoops)
  : DUST_1D(dual_max_type, constraints_type, nbLoops) {}

double Poisson_1D::Cost(unsigned int t, unsigned int s) const
{
  double res = 0.0;
  double diff = cumsum[t] - cumsum[s];
  if(diff > 0.0) /// should be != 0, we use > for robustness
  {
    res = diff * (1.0 - std::log(diff / (t - s)));
  }
  return res;
}

double Poisson_1D::statistic(double& data) const
{return(data);}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Poisson_1D::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  double objectiveMean = (cumsum[t] - cumsum[s]) / (t - s); // m_it
  double constraintMean = (cumsum[s] - cumsum[r]) / (s - r); // m_ji

  ///
  /// point in the right interval
  if(constraintMean != 0){point = point * std::min(1.0, objectiveMean/constraintMean);}
  ///

  return (costRecord[s] - minCost) / (t - s)
    + point * (costRecord[s] - costRecord[r]) / (s - r)
    - (objectiveMean - point * constraintMean) * (std::log((objectiveMean - point * constraintMean) / (1 - point)) - 1);
}


double Poisson_1D::dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  return (-std::numeric_limits<double>::infinity());
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Poisson_1D::muMax(double a, double b) const
{
  if (b != 0) return std::min(1., a/b);
  return 1.;
}

double Poisson_1D::Dstar(double x) const
{
  return (x * (std::log(x) - 1.0));
}


double Poisson_1D::DstarPrime(double x) const
{
  return std::log(x);
}

double Poisson_1D::DstarSecond(double x) const
{
  return pow(x, -1);
}

std::string Poisson_1D::get_model() const { return "poisson"; }
