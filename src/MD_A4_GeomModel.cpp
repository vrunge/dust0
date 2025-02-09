#include <cmath>

#include "MD_A4_GeomModel.h"

using namespace Rcpp;

Geom_MD::Geom_MD(int dual_max_type, int constraints_type, Nullable<unsigned> nbLoops)
  : DUST_MD(dual_max_type, constraints_type, nbLoops) {}

double Geom_MD::Cost(const unsigned int& t, const unsigned int& s) const
{
  double res = 0;
  double ratio;
  double diff;
  double delta = double(t - s);
  double inv_delta = pow(delta, -1);
  for (unsigned int row = 0; row < dim; row++)
  {
    diff = cumsum(row, t) - cumsum(row, s);
    ratio = diff * inv_delta;
    if (ratio <= 1)
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

std::array<double, 2> Geom_MD::muInterval(const arma::colvec& a, const arma::colvec& b, double& c, double& d) const
{
  std::array<double, 2> interval = {0, std::numeric_limits<double>::infinity() };

  for (unsigned int i = 0; i < a.n_elem; ++i)
  {
    if (a[i] > 0 && b[i] < 0){interval[1] = std::min(interval[1], -a[i]/b[i]);}
    else if (a[i] < 0 && b[i] > 0) {interval[0] = std::max(interval[0], -a[i]/b[i]);}
    if (a[i]-c > 0 && b[i]-d < 0){interval[1] = std::min(interval[1], -(a[i]-c)/(b[i]-d));}
    else if (a[i]-c < 0 && b[i]-d > 0) {interval[0] = std::max(interval[0], -(a[i]-c)/(b[i]-d));}
  }

  if (c > 0 && d < 0)
  {interval[1] = std::min(interval[1], -c / d);}
  else if (c < 0 && d > 0){interval[0] = std::max(interval[0], -c / d);}

  return(interval);
}


void Geom_MD::clipStepSizeModel(const double& m_elem, const arma::rowvec& constraint_means, const double& mu_sum, const arma::rowvec& direction, const double& direction_sum, double& max_stepsize) const
{
  double dot_product = arma::dot(constraint_means, direction);
  double sum_diff = dot_product - direction_sum;

  if (dot_product < 0)
  {
    double new_stepsize = -(1 + mu_sum) * m_elem / dot_product;
    if (new_stepsize < max_stepsize)
      max_stepsize = new_stepsize;
  }

  if (sum_diff > 0)
  {
    double new_stepsize = (1 + mu_sum) * (1 - m_elem) / sum_diff;
    if (new_stepsize < max_stepsize)
      max_stepsize = new_stepsize;
  }
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Geom_MD::dual1D_Eval(double& point, const arma::colvec& a, const arma::colvec& b, double& c, double& d, double& e, double& f) const
{
  return(-std::numeric_limits<double>::infinity());
}

////////////////////////////////////////////////////////////////////////////////

double Geom_MD::dual1D_Max(double& argmax, arma::colvec& a, arma::colvec& b, double& c, double& d, double& e, double& f) const
{
  double Max = -std::numeric_limits<double>::infinity();

  return(Max);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

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

std::string Geom_MD::get_model() const { return "geom"; }
