#include <Rcpp.h>
#include <algorithm> // for std::min
#include <cmath>

#include <random>

#include "1D_A2_PoissonModel.h"

using namespace Rcpp;

Poisson_1D::Poisson_1D(int dual_max, bool random_constraint, Nullable<double> alpha, Nullable<int> nbLoops)
  : DUST_1D(dual_max, random_constraint, alpha, nbLoops) {}

double Poisson_1D::Cost(unsigned int t, unsigned int s) const
{
  double res = 0.0;
  double diff = cumsum[t] - cumsum[s];
  if(diff != 0.0) {
    double inv_diff = 1.0 / (t - s);
    double log_term = std::log(diff * inv_diff);
    res = diff * (1.0 - log_term);
  }
  return res;
}

double Poisson_1D::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  double objectiveMean = (cumsum[t] - cumsum[s]) / (t - s); // m_it
  double constraintMean = (cumsum[s] - cumsum[r]) / (s - r); // m_ji

  ///
  /// point in the right interval
  if(constraintMean != 0){point = point * std::min(1.0, objectiveMean/constraintMean);}
  ///

  return (costRecord[s] - minCost) / (t - s)
    + point * (costRecord[s] - costRecord[r]) / (s - r)
    - (objectiveMean - point * constraintMean) * (std::log((objectiveMean - point * constraintMean) / (1 - point)) - 1);
}


double Poisson_1D::dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  return (-std::numeric_limits<double>::infinity());
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Poisson_1D::Dstar(double x) const
{
  return 0;
}


double Poisson_1D::DstarPrime(double x) const
{
  return 0;
}

double Poisson_1D::DstarSecond(double x) const
{
  return 0;
}





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// double a = (cumsum[t] - cumsum[s]) / (t - s); // m_it
// double b = (cumsum[s] - cumsum[r]) / (s - r); // m_ji
// double C = (costRecord[s] - costRecord[r]) / (s - r);

// std::cout << "a. " << a << " ---"  << "b. " << b << " ---" << "C. " << C << " ---"<< std::endl;
// double f_prime;
// double f_second;
// double mu_new;

// double mu = 0.5;
// double m;
// std::cout << "NEWWWWW ";
// for (int i = 0; i < 6; ++i)
// {
//   m =  (a - mu*b) / (1 - mu);
//   f_prime = Poisson_1D::Dstar(m) - ((a-b)/(1-mu)) * Poisson_1D::DstarPrime(m) + C;
//   f_second = Poisson_1D::DstarPrime(m)*(1 - (a-b)/pow(1-mu,2))- pow(a-b,2)/pow(1-mu,3) * Poisson_1D::DstarSecond(m);

  //   mu_new = mu - f_prime / f_second;
  //   mu_new = std::max(0.0, std::min(1.0, mu_new));
  // Check for pruning
  // if (Poisson_1D::dualEval(mu_new, minCost, t, s, r) > 0) {break;}

  //   std::cout << "f_prime. " << mu << " &&& " << m << " --- " << f_prime  << " +++ " << f_second << std::endl;

  //   mu = mu_new;
  // }

// return Poisson_1D::dualEval(mu, minCost, t, s, r);
