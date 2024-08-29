#include <Rcpp.h>
#include <cmath>

#include "1D_A7_NegbinModel.h"

using namespace Rcpp;

Negbin_1D::Negbin_1D(bool use_dual_max, bool random_constraint, Nullable<double> alpha, Nullable<int> nbLoops)
  : DUST_1D(use_dual_max, random_constraint, alpha, nbLoops) {}

double Negbin_1D::Cost(unsigned int t, unsigned int s) const
{
  double res = 0;
  double m = (cumsum[t] - cumsum[s])/(t - s);
  if(m != 0)
  {res = double(t - s) * std::log(1 + m) - (cumsum[t] - cumsum[s]) * std::log(m / (1 + m));}
  return res;
}

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
  const double phi = (1 + sqrt(5)) / 2;  // Golden ratio
  double a = 0.0;
  double b = 1.0;
  double c = 1 - 1/phi;
  double d = 1/phi;

  double fc = Negbin_1D::dualEval(c, minCost, t, s, r);
  double fd = Negbin_1D::dualEval(d, minCost, t, s, r);
  double max_val = std::max(fc, fd);

  for (int i = 0; i < nb_Loops; i++)
  {
    if (fc > fd)
    {
      b = d;
      d = c;
      fd = fc;
      c = b - (b - a) / phi;
      fc = Negbin_1D::dualEval(c, minCost, t, s, r);
    }
    else
    {
      a = c;
      c = d;
      fc = fd;
      d = a + (b - a) / phi;
      fd = Negbin_1D::dualEval(d, minCost, t, s, r);
    }
    max_val = std::max(max_val, std::max(fc, fd));
    if(max_val > 0){break;}
  }
  return max_val;
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Negbin_1D::Dstar(double x) const
{
  return 0;
}


double Negbin_1D::DstarPrime(double x) const
{
  return 0;
}


double Negbin_1D::DstarSecond(double x) const
{
  return 0;
}


