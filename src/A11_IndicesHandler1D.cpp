// [[Rcpp::depends(Rcpp)]]

#include <Rcpp.h>
#include "A11_IndicesHandler1D.h"

using namespace Rcpp;

// --- // init_handler creates the correct IndicesModule object // --- //
// --- // and feeds it to parent class                          // --- //

void RandomConstraint1D::init_handler()
{
  set_method(std::make_unique<I_Random>(n, alpha));
}

void DeterministicConstraint1D::init_handler()
{
  set_method(std::make_unique<I_Deterministic>());
}
