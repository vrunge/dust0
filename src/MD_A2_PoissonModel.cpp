#include "MD_A2_PoissonModel.h"

using namespace Rcpp;

Poisson_MD::Poisson_MD(int dual_max_type, int constraints_type, Nullable<unsigned> nbLoops)
  : DUST_MD(dual_max_type, constraints_type, nbLoops) {}

double Poisson_MD::Cost(const unsigned int& t, const unsigned int& s) const
{
  double diff;
  double res = 0;
  double inv_delta = pow(t - s, -1);
  for (unsigned int row = 0; row < dim; row++)
  {
    diff = cumsum(row, t) - cumsum(row, s);
    if (diff <= 0)
      continue;
    res += diff * (1. - std::log(diff * inv_delta));
  }
  return res;
}

double Poisson_MD::statistic(const double& value) const
{
  return value;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Poisson_MD::muMax(const double& a, const double& b) const
{
  if (b != 0) return std::min(1., a/b);
  return 1.;
}


std::array<double, 2> Poisson_MD::muInterval(const arma::colvec& a, const arma::colvec& b, double& c, double& d) const
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


void Poisson_MD::clipStepSizeModel(const double& m_elem, const arma::rowvec& constraint_means, const double& mu_sum, const arma::rowvec& direction, const double& direction_sum, double& max_stepsize) const
{
  // DIRECTION_SUM SHOULD BE POSITIVE HERE
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

double Poisson_MD::dual1D_Eval(double& point, const arma::colvec& a, const arma::colvec& b, double& c, double& d, double& e, double& f) const
{
  return(-std::numeric_limits<double>::infinity());
}

////////////////////////////////////////////////////////////////////////////////

double Poisson_MD::dual1D_Max(double& argmax, arma::colvec& a, arma::colvec& b, double& c, double& d, double& e, double& f) const
{
  double Max = -std::numeric_limits<double>::infinity();

  return(Max);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


double Poisson_MD::Dstar(const double& x) const
{
  return (x * (std::log(x) - 1.0));
}

double Poisson_MD::DstarPrime(const double& x) const
{
  return std::log(x);
}

double Poisson_MD::DstarSecond(const double& x) const
{
  return pow(x, -1);
}

std::string Poisson_MD::get_model() const { return "poisson"; }


