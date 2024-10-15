// --- // Models // --- //
#include "MD_A1_GaussModel.h"
#include "MD_A2_PoissonModel.h"

using namespace Rcpp;

// ---------------------------- //
// --- //////////////////// --- //
// --- // Object factory // --- //
// --- //////////////////// --- //
// ---------------------------- //

DUST_MD *newModuleMD(const std::string& model,
                     const std::string& method,
                     Nullable<double> alpha,
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
  int dual_max = 2;
  bool random_constraint = false;
  if(model == "gauss"){dual_max = 1;}

  if (indices_max[0] == "randIndex"){random_constraint = true;}
  if (indices_max[0] == "detIndex"){random_constraint = false;}


  if (indices_max[1] == "Eval0"){dual_max = 0;} //algo0
  else if (indices_max[1] == "Eval1"){dual_max = 1;} //algo1
  else if (indices_max[1] == "Eval2"){dual_max = 2;} //algo2
  else if (indices_max[1] == "Eval3"){dual_max = 3;} //algo3
  else if (indices_max[1] == "Eval4"){dual_max = 4;} //algo4
  else if (indices_max[1] == "Eval5"){dual_max = 5;} //algo5

  if (model == "gauss")
    return new Gauss_MD(dual_max, random_constraint, alpha, nbLoops);
  if (model == "poisson")
    return new Poisson_MD(dual_max, random_constraint, alpha, nbLoops);
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

  .factory<const std::string&, const std::string&, Nullable<double>, Nullable<int>>(newModuleMD)

  .method("init_raw", &DUST_MD::init)
  .method("compute", &DUST_MD::compute)
  .method("get_partition", &DUST_MD::get_partition)
  .method("quick_raw", &DUST_MD::quick)
  ;
}

