#include "utils.h"

using namespace Rcpp;


// ######### // FINDBOUNDARY // ######### //
// -> clip direction to restrict optimization within the simplex / distorted simplex (i.e. if mu_max not 1)

void clip_stepsize_to_negative_element(const arma::rowvec &mu, const arma::rowvec &direction, double &max_stepsize)
{
  unsigned size = mu.n_elem;
  for (unsigned i = 0; i < size; i++)
  {
    double distance_to_zero;
    if (mu(i) > 0)
    {
      if (direction(i) > 0)
      {
        continue;
      }
      double distance_to_zero = -mu(i) / direction(i);
      if (distance_to_zero < max_stepsize)
      {
        max_stepsize = distance_to_zero;
      }
      continue;
    }
    if (direction(i) < 0)
    {
      continue;
    }
    distance_to_zero = -mu(i) / direction(i);
    if (distance_to_zero < max_stepsize)
    {
      max_stepsize = distance_to_zero;
    }
  }
}

void clip_stepsize_to_negative_sum(const std::vector<int>& sign, const double& mu_sum, const arma::rowvec& direction, double& direction_sum, double& max_stepsize)
{
  for (unsigned i = 0; i < direction.n_elem; i++)
  {
    direction_sum += sign[i] * direction(i);
  }

  if (direction_sum < 0)
  {
    double new_stepsize = -(1 + mu_sum) / direction_sum;
    if (new_stepsize < max_stepsize)
    {
      max_stepsize = new_stepsize;
    }
  }
}

double FindBoundaryCoef(const arma::rowvec &x,
                        const arma::rowvec &direction,
                        const arma::rowvec &weights,
                        bool &shrink,
                        std::vector<unsigned int> &shrink_indices)
{
  double max_coef = 1.; // record current min
  double coef; // store current value
  double sum_x = 0.; // for sum condition (sum < 1)
  double sum_d = 0.;
  for (unsigned int col = 0; col < x.size(); col++)
  {
    if (x(col) == 0. && direction(col) < 0.)
    { 
      shrink = true; shrink_indices.push_back(col);
    }; // if point already on boundary, do not push further

    double x_norm = x(col) * weights(col);
    double d_norm = direction(col) * weights(col);
    sum_x += x_norm;
    sum_d += d_norm;

    coef = x_norm / d_norm; // coef for reaching boundary i
    if (coef > 0.)
    {
      if (coef < max_coef) 
      {
        max_coef = coef; // update closest boundary
      }
    }
  }

  if (sum_d == 0.)
  { 
    return 0.;
  }

  coef = (1 - 1e-9) * (1 + sum_x) / sum_d; // sum boundary (boundary defined by sum(x / mu_max) == 1)
  // coef = (1 - 1e-9) * (1 - sum_x) / sum_d; // sum boundary (boundary defined by sum(x / mu_max) == 1)
  if (coef > 0)
  {
    if (coef < max_coef)
    {
      max_coef = coef;
    }
  }

  return max_coef;
  // if (min_coef == std::numeric_limits<double>::infinity()) return 0.;
};

/////////////// ///////////////
/////////////// ///////////////
/////////////// ///////////////

void updateHessian(arma::dmat &inverseHessian,
                   const arma::rowvec &mu_diff,
                   const arma::rowvec &grad_diff,
                   const arma::dmat &I)
{
  double rho = pow(arma::dot(grad_diff, mu_diff), -1);
  arma::dmat mult = I - mu_diff.t() * (grad_diff * rho);
  inverseHessian = mult * inverseHessian * mult.t() + mu_diff.t() * (mu_diff * rho);
}
