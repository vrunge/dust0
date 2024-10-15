#include <cmath>

#include "MD_A8_VarianceModel.h"

using namespace Rcpp;

Variance_MD::Variance_MD(int dual_max, bool random_constraint, Nullable<double> alpha, Nullable<int> nbLoops)
  : DUST_MD(dual_max, random_constraint, alpha, nbLoops) {}

double Variance_MD::Cost(const unsigned int& t, const unsigned int& s) const
{
  double res = 0;
  double delta = t - s;
  double inv_delta = pow(t - s, -1);
  double diff;
  for (unsigned int row = 0; row < d; row++)
  {
    diff = cumsum(row, t) - cumsum(row, s);
    res += 0.5 * delta * (1.0 + std::log(diff * inv_delta));
  }
  return res;
}

double Variance_MD::statistic(const double& data) const
{return(data * data);}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Variance_MD::muMax(const double& a, const double& b) const
{
  double res = 1;
  if(b != 0){res = std::min(1.0, a/b);}
  return res;
}

double Variance_MD::Dstar(const double& x) const
{
  return -0.5 * (std::log(x) + 1.0);
}


double Variance_MD::DstarPrime(const double& x) const
{
  return -0.5 /x;
}

double Variance_MD::DstarSecond(const double& x) const
{
  return 0.5/std::pow(x,2);
}


