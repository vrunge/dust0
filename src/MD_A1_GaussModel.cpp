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

double Gauss_MD::dual1DMax(arma::colvec& a, arma::colvec& b, double& c, double& d, double& e, double& f) const
{
  /// Find mu_l and mu_r = dualInterval
  //std::array<double, 2>  mu_interval = muInterval(a,b,c,d);

  // function -0.5 sum (a_i + X b_i)^2/(c + x d) -e X - f
  double A2 = arma::dot(a,a);
  double B2 = arma::dot(b,b);
  double AB = arma::dot(a,b);
  // roots of Ax^2 + Bx + C ?
  double A = d*B2 + 2*e*d*d;
  double B = 2*c*B2 + 4*c*d*e;
  double C = 2*c*AB - d*A2 + 2*c*c*e;

  double delta = B*B - 4*A*C;

  return 0.;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Gauss_MD::muMax(const double& a, const double& b) const
{
  return 1.;
}


std::array<double, 2> Gauss_MD::muInterval(const arma::colvec& a, const arma::colvec& b, double& c, double& d) const
{
  std::array<double, 2> interval = {0, std::numeric_limits<double>::infinity() };
  if (c > 0 && d < 0)
  {
    interval[1] = -c / d;
  }
  else if (c < 0 && d > 0)
  {
    interval[0] = -c / d;
  }
  return(interval);
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



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





