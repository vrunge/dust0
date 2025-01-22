#include <cmath>

#include "MD_A3_ExpModel.h"

using namespace Rcpp;

Exp_MD::Exp_MD(int dual_max, bool random_constraint, Nullable<int> nbLoops)
  : DUST_MD(dual_max, random_constraint, nbLoops) {}

double Exp_MD::Cost(const unsigned int& t, const unsigned int& s) const
{
  double diff;
  double res = 0;
  double inv_delta = pow(t - s, -1);
  for (unsigned int row = 0; row < d; row++)
  {
    diff = cumsum(row, t) - cumsum(row, s);
    if(diff <= 0){diff = 1e-100;} /// choice  1e-100 to avoid -Inf
    res += (1.0 + std::log(diff * inv_delta));
  }
  return res * (t - s);
}

double Exp_MD::statistic(const double& data) const
{return(data);}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Exp_MD::muMax(const double& a, const double& b) const
{
  double res = 1;
  if(b != 0){res = std::min(1.0, a/b);}
  return res;
}

double Exp_MD::Dstar(const double& x) const
{
  return (-std::log(x) - 1.0);
}


double Exp_MD::DstarPrime(const double& x) const
{
  return -1.0/x;
}

double Exp_MD::DstarSecond(const double& x) const
{
  return 1.0/std::pow(x,2);
}

std::string Exp_MD::get_model() const { return "exp"; }


