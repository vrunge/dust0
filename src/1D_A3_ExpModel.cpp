#include <Rcpp.h>
#include <cmath>

#include "1D_A3_ExpModel.h"

using namespace Rcpp;

Exp_1D::Exp_1D(bool use_dual_max, bool random_constraint, Nullable<double> alpha, Nullable<int> nbLoops)
  : DUST_1D(use_dual_max, random_constraint, alpha, nbLoops) {}

double Exp_1D::Cost(unsigned int t, unsigned int s) const
{
  double delta_t = t - s;
  double diff_cumsum = cumsum[t] - cumsum[s];
  return delta_t * (1.0 + std::log(diff_cumsum / delta_t));
}

double Exp_1D::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const
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
  + (1 - point) * (std::log((objectiveMean - point * constraintMean) / (1 - point)) + 1);
}

double Exp_1D::dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  const double phi = (1 + sqrt(5)) / 2;  // Golden ratio
  double a = 0.0;
  double b = 1.0;
  double c = 1 - 1/phi;
  double d = 1/phi;

  double fc = Exp_1D::dualEval(c, minCost, t, s, r);
  double fd = Exp_1D::dualEval(d, minCost, t, s, r);
  double max_val = std::max(fc, fd);

  for (int i = 0; i < nb_Loops; i++)
  {
    if (fc > fd)
    {
      b = d;
      d = c;
      fd = fc;
      c = b - (b - a) / phi;
      fc = Exp_1D::dualEval(c, minCost, t, s, r);
    }
    else
    {
      a = c;
      c = d;
      fc = fd;
      d = a + (b - a) / phi;
      fd = Exp_1D::dualEval(d, minCost, t, s, r);
    }
    max_val = std::max(max_val, std::max(fc, fd));
    if(max_val > 0){break;}
  }
  return max_val;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Exp_1D::Dstar(double x) const
{
  return 0;
}


double Exp_1D::DstarPrime(double x) const
{
  return 0;
}

double Exp_1D::DstarSecond(double x) const
{
  return 0;
}


