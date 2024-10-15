#include <cmath>

#include "MD_A4_GeomModel.h"

using namespace Rcpp;

Geom_MD::Geom_MD(int dual_max, bool random_constraint, Nullable<double> alpha, Nullable<int> nbLoops)
  : DUST_MD(dual_max, random_constraint, alpha, nbLoops) {}

double Geom_MD::Cost(const unsigned int& t, const unsigned int& s) const
{
  double res = 0;
  double ratio;
  double diff;
  double delta = t - s;
  double inv_delta = pow(delta, -1);
  for (unsigned int row = 0; row < d; row++)
  {
    diff = cumsum(row, t) - cumsum(row, s);
    ratio = diff * inv_delta;
    if (diff != 1)
      continue;
    res += delta * std::log(ratio - 1) - diff * std::log((ratio - 1) / ratio);
  }
return res;
}



double Geom_MD::statistic(const double& data) const
{return(data);}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Geom_MD::muMax(const double& a, const double& b) const
{
  double res = 1;
  if(b != 1){res = std::min(1.0, (a-1)/(b-1));}
  return res;
}

double Geom_MD::Dstar(const double& x) const
{
  return (x - 1.0)*std::log(x - 1.0) - x*std::log(x);
}


double Geom_MD::DstarPrime(const double& x) const
{
  return std::log(x - 1.0) - std::log(x);
}

double Geom_MD::DstarSecond(const double& x) const
{
  return 1.0/(x - 1.0) - 1.0/x;
}





