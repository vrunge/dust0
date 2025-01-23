#include <cmath>

#include "MD_A8_VarianceModel.h"

using namespace Rcpp;

Variance_MD::Variance_MD(int dual_max, bool random_constraint, Nullable<int> nbLoops)
  : DUST_MD(dual_max, random_constraint, nbLoops) {}

double Variance_MD::Cost(const unsigned int& t, const unsigned int& s) const
{
  double res = 0;
  double delta = t - s;
  double inv_delta = pow(t - s, -1);
  double diff;
  for (unsigned int row = 0; row < d; row++)
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

double Variance_MD::dualMax(arma::colvec& a, arma::colvec& b, double& c, double& d, double& e, double& f) const
{
  return 0.;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Variance_MD::muMax(const double& a, const double& b) const
{
  double res = 1;
  if(b != 0){res = std::min(1.0, a/b);}
  return res;
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
