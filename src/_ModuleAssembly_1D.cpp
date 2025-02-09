#include <Rcpp.h>

// --- // Models // --- //
#include "1D_A1_GaussModel.h"
#include "1D_A2_PoissonModel.h"
#include "1D_A3_ExpModel.h"
#include "1D_A4_GeomModel.h"
#include "1D_A5_BernModel.h"
#include "1D_A6_BinomModel.h"
#include "1D_A7_NegbinModel.h"
#include "1D_A8_VarianceModel.h"

using namespace Rcpp;

// ---------------------------- //
// --- //////////////////// --- //
// --- // Object factory // --- //
// --- //////////////////// --- //
// ---------------------------- //

DUST_1D *newModule1D(const std::string& model,
                     const std::string& method,
                     Nullable<int> nbLoops)
{
  ///////////////////  method separation into 2 ///////////////////
  std::vector<std::string> method_INFO;
  size_t pos = method.find('_');  // Find the position of the underscore

  if (pos != std::string::npos)
  {
    method_INFO.push_back(method.substr(0, pos));        // First part before the underscore
    method_INFO.push_back(method.substr(pos + 1));       // Second part after the underscore
  }
  else
  {
    method_INFO.push_back(method);
    method_INFO.push_back(method);
  }

  ///////////////////  DEFAULT CHOICE  /////////////////////////////////
  ///////////////////  DEFAULT CHOICE  = best choice ///////////////////
  ///////////////////  DEFAULT CHOICE  /////////////////////////////////
  /// DEFAULT FASTEST CHOICE
  /// DEFAULT FASTEST CHOICE
  int dual_max_type = 4; /// quasi newton
  int constraints_type = 0; /// deterministic choice for constraint
  if(model == "gauss"){dual_max_type = 1;} /// exact max eval
  /// DEFAULT FASTEST CHOICE
  /// DEFAULT FASTEST CHOICE

  if (method_INFO[0] == "randIndex"){constraints_type = 0;}
  else if (method_INFO[0] == "detIndex"){constraints_type = 1;}

  if (method_INFO[1] == "Eval0"){dual_max_type = 0;} //algo0
  else if (method_INFO[1] == "Eval1"){dual_max_type = 1;} //algo1
  else if (method_INFO[1] == "Eval2"){dual_max_type = 2;} //algo2
  else if (method_INFO[1] == "Eval3"){dual_max_type = 3;} //algo3
  else if (method_INFO[1] == "Eval4"){dual_max_type = 4;} //algo4
  else if (method_INFO[1] == "Eval5"){dual_max_type = 5;} //algo5
  else if (method_INFO[1] == "Eval6"){dual_max_type = 6;} //algo6

  if (model == "gauss")
    return new Gauss_1D(dual_max_type, constraints_type, nbLoops);
  if (model == "poisson")
    return new Poisson_1D(dual_max_type, constraints_type, nbLoops);
  if (model == "exp")
    return new Exp_1D(dual_max_type, constraints_type, nbLoops);
  if (model == "geom")
    return new Geom_1D(dual_max_type, constraints_type, nbLoops);
  if (model == "bern")
    return new Bern_1D(dual_max_type, constraints_type, nbLoops);
  if (model == "binom")
    return new Binom_1D(dual_max_type, constraints_type, nbLoops);
  if (model == "negbin")
    return new Negbin_1D(dual_max_type, constraints_type, nbLoops);
  if (model == "variance")
    return new Variance_1D(dual_max_type, constraints_type, nbLoops);
  return nullptr;
}


// --------------------------------- //
// --- ///////////////////////// --- //
// --- // Exposing the module // --- //
// --- ///////////////////////// --- //
// --------------------------------- //

//' @title MyModule: Exposing DUST_1D to R
//'
//' @name DUST_1D
//'
//' @description
//' This module exposes the \code{DUST_1D} C++ class to R, allowing you to create
//' instances of \code{DUST_1D} and call its methods directly from R.
//'
//' @export
RCPP_MODULE(DUSTMODULE1D)
{
  class_<DUST_1D>("DUST_1D")

    .factory<const std::string&, const std::string&, Nullable<int>>(newModule1D)

    .method("append_c", &DUST_1D::append)
    .method("update_partition", &DUST_1D::update_partition)
    .method("get_partition", &DUST_1D::get_partition)
    .method("get_info", &DUST_1D::get_info)
    .method("one_dust", &DUST_1D::one_dust)
  ;
}

