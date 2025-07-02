#include <Rcpp.h>
#include <cmath>

#include "1D_A7_NegbinModel.h"

using namespace Rcpp;

Negbin_1D::Negbin_1D(int dual_max_type, int constraints_type, Nullable<int> nbLoops)
  : DUST_1D(dual_max_type, constraints_type, nbLoops) {}

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


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Negbin_1D::dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  /// -Dstar(a - (b-a)*x) - (c - (d-c)*x)
  double a = (cumsum[t] - cumsum[s]) / (t - s);
  double b = (cumsum[s] - cumsum[r]) / (s - r);
  double c = (minCost - costRecord[s]) / (t - s);
  double d = (costRecord[s] - costRecord[r]) / (s - r);
  double mu_max = muMax(a, b);

  double ratio = -(d - c)/(b - a);
  double val = DstarPrimeInv(ratio);

  double x = 1.0/(b - a) * (a - val);
  double x_max = mu_max/(1 - mu_max);
  if(x < 0){x = 0;}
  if(x > x_max)
  {
    x = 0;
  }
  return -Dstar(a - (b-a)*x) - (c - (d-c)*x);
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

bool Negbin_1D::isBoundary(double a) const
{
  if(a == 0){return true;}
  return false;
}

double Negbin_1D::Dstar(double x) const
{
  return x*std::log(x) - (1.0 + x)*std::log(1.0 + x);
}


double Negbin_1D::DstarPrime(double x) const
{
  return std::log(x) - std::log(1.0 + x);
}

double Negbin_1D::DstarPrimeInv(double x) const
{
  return std::exp(x) / (1 - std::exp(x));
}


double Negbin_1D::DstarSecond(double x) const
{
  return 1.0/x + 1.0/(1.0 + x);
}

std::string Negbin_1D::get_model() const { return "negbin"; }
