#include <Rcpp.h>

// --- // Models // --- //
#include "2D_DUSTmeanVar.h"
#include "2D_DUSTreg.h"

using namespace Rcpp;
using namespace std;
// ---------------------------- //
// --- //////////////////// --- //
// --- // Object factory // --- //
// --- //////////////////// --- //
// ---------------------------- //

DUST_meanVar *newModuleMeanVar(const std::string& method,
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

  ///////////////////  DEFAULT CHOICE  ///////////////////
  int constraint_indices = 11;
  int dual_max_type = 2;

  if (method_INFO[0] == "randIndex"){constraint_indices = 10;}
  if (method_INFO[0] == "detIndex"){constraint_indices = 11;}
  if (method_INFO[0] == "rand2Index"){constraint_indices = 20;}
  if (method_INFO[0] == "det2Index"){constraint_indices = 21;}

  if (method_INFO[1] == "Eval0"){dual_max_type = 0;} //algo0
  if (method_INFO[1] == "Eval1"){dual_max_type = 1;} //algo1
  if (method_INFO[1] == "Eval2"){dual_max_type = 2;} //algo2
  if (method_INFO[1] == "Eval3"){dual_max_type = 3;} //algo3
  if (method_INFO[1] == "Eval4"){dual_max_type = 4;} //algo4
  if (method_INFO[1] == "Eval5"){dual_max_type = 5;} //algo5
  if (method_INFO[1] == "Eval6"){dual_max_type = 6;} //algo6

  return new DUST_meanVar(dual_max_type, constraint_indices, nbLoops);
}


// --------------------------------- //
// --- ///////////////////////// --- //
// --- // Exposing the module // --- //
// --- ///////////////////////// --- //
// --------------------------------- //

//' @title MyModule: Exposing DUST_meanVar to R
//'
//' @name DUST_meanVar
//'
//' @description
//' This module exposes the \code{DUST_meanVar} C++ class to R, allowing you to create
//' instances of \code{DUST_meanVar} and call its methods directly from R.
//'
//' @export
RCPP_MODULE(DUSTMODULEMeanVar)
{
  class_<DUST_meanVar>("DUST_meanVar")

  .factory<const std::string&, Nullable<int>>(newModuleMeanVar)

  .method("init_raw", &DUST_meanVar::init)
  .method("compute", &DUST_meanVar::compute)
  .method("get_partition", &DUST_meanVar::get_partition)
  .method("quick_raw", &DUST_meanVar::quick)
  ;
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


DUST_reg *newModuleReg(const std::string& method,
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

  ///////////////////  DEFAULT CHOICE  ///////////////////
  int constraint_indices = 11;
  int dual_max_type = 2;

  if (method_INFO[0] == "randIndex"){constraint_indices = 10;}
  if (method_INFO[0] == "detIndex"){constraint_indices = 11;}
  if (method_INFO[0] == "rand2Index"){constraint_indices = 20;}
  if (method_INFO[0] == "det2Index"){constraint_indices = 21;}

  if (method_INFO[1] == "Eval0"){dual_max_type = 0;} //algo0
  if (method_INFO[1] == "Eval1"){dual_max_type = 1;} //algo1
  if (method_INFO[1] == "Eval2"){dual_max_type = 2;} //algo2
  if (method_INFO[1] == "Eval3"){dual_max_type = 3;} //algo3
  if (method_INFO[1] == "Eval4"){dual_max_type = 4;} //algo4
  if (method_INFO[1] == "Eval5"){dual_max_type = 5;} //algo5
  if (method_INFO[1] == "Eval6"){dual_max_type = 6;} //algo6

  return new DUST_reg(dual_max_type, constraint_indices, nbLoops);
}


// --------------------------------- //
// --- ///////////////////////// --- //
// --- // Exposing the module // --- //
// --- ///////////////////////// --- //
// --------------------------------- //

//' @title MyModule: Exposing DUST_reg to R
//'
//' @name DUST_reg
//'
//' @description
//' This module exposes the \code{DUST_reg} C++ class to R, allowing you to create
//' instances of \code{DUST_reg} and call its methods directly from R.
//'
//' @export
RCPP_MODULE(DUSTMODULEreg)
{
  class_<DUST_reg>("DUST_reg")

  .factory<const std::string&, Nullable<int>>(newModuleReg)

  .method("init_raw", &DUST_reg::init)
  .method("compute", &DUST_reg::compute)
  .method("get_partition", &DUST_reg::get_partition)
  .method("quick_raw", &DUST_reg::quick)
  ;
}



