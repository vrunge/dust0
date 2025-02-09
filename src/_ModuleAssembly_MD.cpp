// --- // Models // --- //
#include "MD_A_DUST.h"
#include "MD_A1_GaussModel.h"
#include "MD_A2_PoissonModel.h"
#include "MD_A3_ExpModel.h"
#include "MD_A4_GeomModel.h"
#include "MD_A5_BernModel.h"
#include "MD_A6_BinomModel.h"
#include "MD_A7_NegbinModel.h"
#include "MD_A8_VarianceModel.h"

using namespace Rcpp;

// ---------------------------- //
// --- //////////////////// --- //
// --- // Object factory // --- //
// --- //////////////////// --- //
// ---------------------------- //

DUST_MD *newModuleMD(const std::string& model,
                     const std::string& method,
                     Nullable<unsigned> nbLoops)
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
  /// FASTEST CHOICE
  /// FASTEST CHOICE
  int dual_max_type = 4; /// quasi newton
  int constraints_type = 0; /// deterministic choice
  if(model == "gauss"){dual_max_type = 1;}
  /// FASTEST CHOICE
  /// FASTEST CHOICE

  if (method_INFO[0] == "randIndex"){constraints_type = 0;} // rand index choice
  else if (method_INFO[0] == "detIndex"){constraints_type = 1;} // det index choice

  if (method_INFO[1] == "Eval0") {dual_max_type = 0;} //algo0
  else if (method_INFO[1] == "Eval1"){dual_max_type = 1;} //algo1
  else if (method_INFO[1] == "Eval2"){dual_max_type = 2;} //algo2
  else if (method_INFO[1] == "Eval3"){dual_max_type = 3;} //algo3
  else if (method_INFO[1] == "Eval4"){dual_max_type = 4;} //algo4
  else if (method_INFO[1] == "Eval5"){dual_max_type = 5;} //algo5
  else if (method_INFO[1] == "Eval6"){dual_max_type = 6;} //algo6

  if (model == "gauss")
    return new Gauss_MD(dual_max_type, constraints_type, nbLoops);
  if (model == "poisson")
    return new Poisson_MD(dual_max_type, constraints_type, nbLoops);
  if (model == "exp")
    return new Exp_MD(dual_max_type, constraints_type, nbLoops);
  if (model == "geom")
    return new Geom_MD(dual_max_type, constraints_type, nbLoops);
  if (model == "bern")
    return new Bern_MD(dual_max_type, constraints_type, nbLoops);
  if (model == "binom")
    return new Binom_MD(dual_max_type, constraints_type, nbLoops);
  if (model == "negbin")
    return new Negbin_MD(dual_max_type, constraints_type, nbLoops);
  if (model == "variance")
    return new Variance_MD(dual_max_type, constraints_type, nbLoops);
  return nullptr;
}


// --------------------------------- //
// --- ///////////////////////// --- //
// --- // Exposing the module // --- //
// --- ///////////////////////// --- //
// --------------------------------- //

//' @title MyModule: Exposing DUST_MD to R
//'
//' @name DUST_MD
//'
//' @description
//' This module exposes the \code{DUST_MD} C++ class to R, allowing you to create
//' instances of \code{DUST_MD} and call its methods directly from R.
//'
//' @export
RCPP_MODULE(DUSTMODULEMD)
{
  class_<DUST_MD>("DUST_MD")

  .factory<const std::string&, const std::string&, Nullable<unsigned>>(newModuleMD)

  .method("append_c", &DUST_MD::append)
  .method("update_partition", &DUST_MD::update_partition)
  .method("get_partition", &DUST_MD::get_partition)
  .method("get_info", &DUST_MD::get_info)
  .method("one_dust", &DUST_MD::one_dust)
  ;
}




