#include <cmath>

#include "MD_A7_NegbinModel.h"

using namespace Rcpp;

Negbin_MD::Negbin_MD(int dual_max, bool random_constraint, Nullable<double> alpha, Nullable<int> nbLoops)
  : DUST_MD(dual_max, random_constraint, alpha, nbLoops) {}

double Negbin_MD::Cost(const unsigned int& t, const unsigned int& s) const
{
  double res = 0;
  double delta = t - s;
  double inv_delta = pow(t - s, -1);
  double diff;
  double ratio;
  for (unsigned int row = 0; row < d; row++)
  {
    diff = cumsum(row, t) - cumsum(row, s);
    ratio = diff * inv_delta;
    if (ratio == 0 || ratio == 1)
      continue;
    res += delta * std::log(1 + ratio) - diff * std::log(ratio / (1 + ratio));
  }
  return res;
}

double Negbin_MD::statistic(const double& data) const
{return(data);}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Negbin_MD::muMax(const double& a, const double& b) const
{
  double res = 1;
  if(b != 0){res = std::min(1.0, a/b);}
  return res;
}

double Negbin_MD::Dstar(const double& x) const
{
  return x*std::log(x) - (1.0 + x)*std::log(1.0 + x);
}


double Negbin_MD::DstarPrime(const double& x) const
{
  return std::log(x) - std::log(1.0 + x);
}


double Negbin_MD::DstarSecond(const double& x) const
{
  return 1.0/x + 1.0/(1.0 + x);
}


