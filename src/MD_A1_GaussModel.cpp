#include "MD_A1_GaussModel.h"

using namespace Rcpp;

Gauss_MD::Gauss_MD(int dual_max, bool random_constraint, Nullable<int> nbLoops)
  : DUST_MD(dual_max, random_constraint, nbLoops) {}

double Gauss_MD::Cost(const unsigned int& t, const unsigned int& s) const
{
  double result = 0;
  for (unsigned int row = 0; row < d; row++)
    result += pow(cumsum(row, t) - cumsum(row, s), 2);

  return - .5 * result / (t - s);
}

double Gauss_MD::statistic(const double& value) const
{
  return value;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Gauss_MD::muMax(const double& a, const double& b) const
{
  return 1.;
}

void Gauss_MD::clipStepSizeModel(const double& m_elem, const arma::rowvec& constraint_means, const double& mu_sum, const arma::rowvec& direction, const double& direction_sum, double& max_stepsize) const
{
  // DIRECTION_SUM MUST BE POSITIVE HERE
}

double Gauss_MD::Dstar(const double& x) const
{
  return .5 * pow(x, 2);
}

double Gauss_MD::DstarPrime(const double& x) const
{
  return x;
}

double Gauss_MD::DstarSecond(const double& x) const
{
  return 1.;
}

std::string Gauss_MD::get_model() const { return "gauss"; }





