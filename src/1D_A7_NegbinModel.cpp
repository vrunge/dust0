#include <Rcpp.h>
#include <cmath>

#include "1D_A7_NegbinModel.h"

using namespace Rcpp;

Negbin_1D::Negbin_1D(int dual_max, bool random_constraint, Nullable<int> nbLoops)
  : DUST_1D(dual_max, random_constraint, nbLoops) {}

double Negbin_1D::Cost(unsigned int t, unsigned int s) const
{
  double res = 0;
  double m = (cumsum[t] - cumsum[s])/(t - s);
  if(m > 0)
  {res = double(t - s) * std::log(1 + m) - (cumsum[t] - cumsum[s]) * std::log(m / (1 + m));}
  return res;
}

double Negbin_1D::statistic(double& data) const
{return(data);}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


double Negbin_1D::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  double objectiveMean = (cumsum[t] - cumsum[s]) / (t - s); // m_it
  double constraintMean = (cumsum[s] - cumsum[r]) / (s - r); // m_ji

  ///
  /// point in the right interval:
  if(constraintMean != 0){point = point * std::min(1.0, objectiveMean/constraintMean);}
  ///
  ///
  double R = (objectiveMean - point * constraintMean) / (1 - point);

  return (costRecord[s] - minCost) / (t - s)
  + point * (costRecord[s] - costRecord[r]) / (s - r)
  + (1 - point) * ((1 + R) * std::log(1 + R) - R * std::log(R));
}


double Negbin_1D::dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  return (-std::numeric_limits<double>::infinity());
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Negbin_1D::muMax(double a, double b) const
{
  double res = 1;
  if(b != 0){res = std::min(1.0, a/b);}
  return res;
}

double Negbin_1D::Dstar(double x) const
{
  return x*std::log(x) - (1.0 + x)*std::log(1.0 + x);
}


double Negbin_1D::DstarPrime(double x) const
{
  return std::log(x) - std::log(1.0 + x);
}


double Negbin_1D::DstarSecond(double x) const
{
  return 1.0/x + 1.0/(1.0 + x);
}


