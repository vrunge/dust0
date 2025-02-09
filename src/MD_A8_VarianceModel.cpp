#include <cmath>

#include "MD_A8_VarianceModel.h"

using namespace Rcpp;

Variance_MD::Variance_MD(int dual_max_type, int constraints_type, Nullable<unsigned> nbLoops)
  : DUST_MD(dual_max_type, constraints_type, nbLoops) {}

double Variance_MD::Cost(const unsigned int& t, const unsigned int& s) const
{
  double res = 0;
  double delta = t - s;
  double inv_delta = pow(t - s, -1);
  double diff;
  for (unsigned int row = 0; row < dim; row++)
  {
    diff = cumsum(row, t) - cumsum(row, s);
    if(diff <= 0){diff = 1e-100;} /// choice  1e-100 to avoid -Inf
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

std::array<double, 2> Variance_MD::muInterval(const arma::colvec& a, const arma::colvec& b, double& c, double& d) const
{
  std::array<double, 2> interval = {0, std::numeric_limits<double>::infinity() };

  for (unsigned int i = 0; i < a.n_elem; ++i)
  {
    if (a[i] > 0 && b[i] < 0){interval[1] = std::min(interval[1], -a[i]/b[i]);}
    else if (a[i] < 0 && b[i] > 0) {interval[0] = std::max(interval[0], -a[i]/b[i]);}
  }
  if (c > 0 && d < 0)
  {interval[1] = std::min(interval[1], -c / d);}
  else if (c < 0 && d > 0){interval[0] = std::max(interval[0], -c / d);}

  return(interval);
}


void Variance_MD::clipStepSizeModel(const double& m_elem, const arma::rowvec& constraint_means, const double& mu_sum, const arma::rowvec& direction, const double& direction_sum, double& max_stepsize) const
{
  double dot_product = arma::dot(constraint_means, direction);
  if (dot_product > 0)
    return;

  double new_stepsize = -(1 + mu_sum) * m_elem / dot_product;
  if (new_stepsize < max_stepsize)
    max_stepsize = new_stepsize;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Variance_MD::dual1D_Eval(double& point, const arma::colvec& a, const arma::colvec& b, double& c, double& d, double& e, double& f) const
{
  return(-std::numeric_limits<double>::infinity());
}

////////////////////////////////////////////////////////////////////////////////



double Variance_MD::dual1D_Max(double& argmax, arma::colvec& a, arma::colvec& b, double& c, double& d, double& e, double& f) const
{
  double Max = -std::numeric_limits<double>::infinity();

  return(Max);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


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

std::string Variance_MD::get_model() const { return "variance"; }
