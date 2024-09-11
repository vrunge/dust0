#include <Rcpp.h>
#include <cmath>

#include "1D_A4_GeomModel.h"

using namespace Rcpp;

Geom_1D::Geom_1D(int dual_max, bool random_constraint, Nullable<double> alpha, Nullable<int> nbLoops)
  : DUST_1D(dual_max, random_constraint, alpha, nbLoops) {}

double Geom_1D::Cost(unsigned int t, unsigned int s) const
{
  double res = 0;
  double m = (cumsum[t] - cumsum[s])/(t - s);
  if(m != 1)
  {res = (t - s) * std::log(m - 1) - (cumsum[t] - cumsum[s]) * std::log((m - 1) / m);}
  return res;
}

double Geom_1D::statistic(double& data) const
{return(data);}

double Geom_1D::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  double objectiveMean = (cumsum[t] - cumsum[s]) / (t - s); // m_it
  double constraintMean = (cumsum[s] - cumsum[r]) / (s - r); // m_ji

  ///
  /// point in the right interval:
  if(constraintMean != 1){point = point * std::min(1.0, (objectiveMean - 1)/(constraintMean - 1));}
  ///
  ///
  double R = (objectiveMean - point * constraintMean) / (1 - point);

  return (costRecord[s] - minCost) / (t - s)
  + point * (costRecord[s] - costRecord[r]) / (s - r)
  + (1 - point) * (R * std::log(R) - (R - 1) * std::log(R - 1));
}

double Geom_1D::dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  return (-std::numeric_limits<double>::infinity());
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Geom_1D::muMax(double a, double b) const
{
  double res = 1;
  if(b != 1){res = std::min(1.0, (a-1)/(b-1));}
  return res;
}

double Geom_1D::Dstar(double x) const
{
  return (x - 1.0)*std::log(x - 1.0) - x*std::log(x);
}


double Geom_1D::DstarPrime(double x) const
{
  return std::log(x - 1.0) - std::log(x);
}

double Geom_1D::DstarSecond(double x) const
{
  return 1.0/(x - 1.0) - 1.0/x;
}





