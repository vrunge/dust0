// --- // Models // --- //
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
                     Nullable<int> nbLoops)
{
  ///////////////////  method separation into 2 ///////////////////
  std::vector<std::string> indices_max;
  size_t pos = method.find('_');  // Find the position of the underscore

  if (pos != std::string::npos)
  {
    indices_max.push_back(method.substr(0, pos));        // First part before the underscore
    indices_max.push_back(method.substr(pos + 1));       // Second part after the underscore
  }
  else
  {
    indices_max.push_back(method);
    indices_max.push_back(method);
  }

  ///////////////////  DEFAULT CHOICE  /////////////////////////////////
  ///////////////////  DEFAULT CHOICE  = best choice ///////////////////
  ///////////////////  DEFAULT CHOICE  /////////////////////////////////
  /// FASTEST CHOICE
  /// FASTEST CHOICE
  int dual_max = 4; /// quasi newton
  bool random_constraint = false; /// deterministic choice
  if(model == "gauss"){dual_max = 1;}

  if (indices_max[0] == "randIndex"){random_constraint = true;}
  else if (indices_max[0] == "detIndex"){random_constraint = false;}

  if (indices_max[1] == "Eval0"){dual_max = 0;} //algo0
  else if (indices_max[1] == "Eval1"){dual_max = 1;} //algo1
  else if (indices_max[1] == "Eval2"){dual_max = 2;} //algo2
  else if (indices_max[1] == "Eval3"){dual_max = 3;} //algo3
  else if (indices_max[1] == "Eval4"){dual_max = 4;} //algo4
  else if (indices_max[1] == "Eval5"){dual_max = 5;} //algo5
  else if (indices_max[1] == "Eval6"){dual_max = 6;} //algo6

  if (model == "gauss")
    return new Gauss_MD(dual_max, random_constraint, nbLoops);
  if (model == "poisson")
    return new Poisson_MD(dual_max, random_constraint, nbLoops);
  if (model == "exp")
    return new Exp_MD(dual_max, random_constraint, nbLoops);
  if (model == "geom")
    return new Geom_MD(dual_max, random_constraint, nbLoops);
  if (model == "bern")
    return new Bern_MD(dual_max, random_constraint, nbLoops);
  if (model == "binom")
    return new Binom_MD(dual_max, random_constraint, nbLoops);
  if (model == "negbin")
    return new Negbin_MD(dual_max, random_constraint, nbLoops);
  if (model == "variance")
    return new Variance_MD(dual_max, random_constraint, nbLoops);
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

  .factory<const std::string&, const std::string&, Nullable<int>>(newModuleMD)

  .method("prepare", &DUST_MD::prepare)
  .method("compute", &DUST_MD::compute)
  .method("get_partition", &DUST_MD::get_partition)
  .method("one_dust", &DUST_MD::one_dust)
  ;
}




