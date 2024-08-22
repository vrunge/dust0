#include <Rcpp.h>

// --- // Models // --- //
#include "1D_A1_GaussModel.h"

using namespace Rcpp;


// ---------------------------- //
// --- //////////////////// --- //
// --- // Object factory // --- //
// --- //////////////////// --- //
// ---------------------------- //

DUST_1D *newModule1D(const std::string& model, const std::string& method, Nullable<double> alpha)
{
  bool use_dual_max;
  bool random_constraint;
  if (method == "all.random")
  {
    use_dual_max = false;
    random_constraint = true;
  }
  else if (method == "half.random")
  {
    use_dual_max = true;
    random_constraint = true;
  }
  else
  {
    use_dual_max = true;
    random_constraint = false;
  }

  if (model == "gauss")
    return new Gauss_1D(use_dual_max, random_constraint, alpha);
  if (model == "poisson")
    return new Gauss_1D(use_dual_max, random_constraint, alpha);
  if (model == "exp")
    return new Gauss_1D(use_dual_max, random_constraint, alpha);
  if (model == "geom")
    return new Gauss_1D(use_dual_max, random_constraint, alpha);
  if (model == "bern")
    return new Gauss_1D(use_dual_max, random_constraint, alpha);
  if (model == "binom")
    return new Gauss_1D(use_dual_max, random_constraint, alpha);
  if (model == "negbin")
    return new Gauss_1D(use_dual_max, random_constraint, alpha);
  return nullptr;
}



// --------------------------------- //
// --- ///////////////////////// --- //
// --- // Exposing the module // --- //
// --- ///////////////////////// --- //
// --------------------------------- //

RCPP_MODULE(DUSTMODULE1D)
{
  class_<DUST_1D>("DUST_1D")

    .factory<const std::string&, const std::string&, Nullable<double>>(newModule1D)

    .method("init_raw", &DUST_1D::init)
    .method("compute", &DUST_1D::compute)
    .method("get_partition", &DUST_1D::get_partition)
    .method("quick_raw", &DUST_1D::quick)
  ;
}
