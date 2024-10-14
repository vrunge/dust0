#include "utils.h"

using namespace Rcpp;



// ######### // FINDBOUNDARYK // ######### //
// -> clip direction to restrict optimization within the simplex / distorted simplex (i.e. if mu_max not 1)
double FindBoundaryCoef(const arma::rowvec& x, const arma::rowvec& dk, const arma::rowvec& weights, bool& shrink, std::vector<unsigned int>& shrink_indices)
{
  unsigned int size = x.n_elem;
  double min_coef = std::numeric_limits<double>::infinity(); // record current min
  double coef; // store current value
  double x_norm;
  double d_norm;
  double sum_x = 0; // for sum condition (sum < 1)
  double sum_d = 0;
  for (unsigned int col = 0; col < size; col++)
  {
    if (x(col) == 0 && dk(col) < 0) { shrink = true; shrink_indices.push_back(col); }; // if point already on boundary, do not push further
    x_norm = x(col) * weights(col);
    d_norm = dk(col) * weights(col);
    sum_x += x_norm;
    sum_d += d_norm;
    coef = - x_norm / d_norm; // coef for reaching boundary i
    if (coef > 0)
    {
      if (coef < min_coef) min_coef = coef; // update closest boundary
    }
  }

  if (sum_d == 0) return 0.;

  coef = (1 - .1) * (1 - sum_x) / sum_d; // sum boundary (boundary defined by sum(x / mu_max) == 1)
  // coef = (1 - 1e-9) * (1 - sum_x) / sum_d; // sum boundary (boundary defined by sum(x / mu_max) == 1)
  if (coef > 0)
  {
    if (coef < min_coef) min_coef = coef;
  }

  if (min_coef < 1) return min_coef;
  if (min_coef == std::numeric_limits<double>::infinity()) return 0.;
  return 1.;
};

void updateHessian(arma::dmat& inverseHessian, const arma::rowvec& mu_diff, const arma::rowvec& grad_diff, const arma::dmat& I)
{
  double rho = pow(arma::dot(grad_diff, mu_diff), -1);
  arma::dmat mult = I - mu_diff.t() * (grad_diff * rho);
  inverseHessian = mult * inverseHessian * mult.t() + mu_diff.t() * (mu_diff * rho);
}
