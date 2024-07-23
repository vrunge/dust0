// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

// --- // Algorithm file // --- //
#include "A1__Algorithm1D.h"

// --- // Model files // --- //
#include "A21_GaussModel1D.h"

using namespace Rcpp;


// ---------------------------- //
// --- //////////////////// --- //
// --- // Gaussian model // --- //
// --- //////////////////// --- //
// ---------------------------- //

// --- // Random evaluation of dual, random constraint in duality problem // --- //
class GaussAllRandom1D
  : public AllRandom1D,
    public GaussModel1D
{
public:
  GaussAllRandom1D(double alpha_) : AllRandom1D(alpha_) {}
};

// --- // Max value of dual, random constraint in duality problem // --- //
class GaussHalfRandom1D
  : public HalfRandom1D,
    public GaussModel1D
{
public:
  GaussHalfRandom1D(double alpha_) : HalfRandom1D(alpha_) {}
};

// --- // Max value of dual, constraint is largest index smaller than tested index // --- //
class GaussDeterministic1D
  : public Deterministic1D,
    public GaussModel1D
{};


// ---------------------------- //
// --- //////////////////// --- //
// --- // Object factory // --- //
// --- //////////////////// --- //
// ---------------------------- //

Skeleton1D *newModule1D(const std::string& model, const std::string& method, double alpha)
{
  if (model == "gauss")
  {
    if (method == "fastest" || method == "deterministic")
    {
      return new GaussDeterministic1D;
    }
    else if (method == "all.random")
    {
      return new GaussAllRandom1D(alpha);
    }
    else if (method == "half.random")
    {
      return new GaussHalfRandom1D(alpha);
    }
  }
  return nullptr;
}


// --------------------------------- //
// --- ///////////////////////// --- //
// --- // Exposing the module // --- //
// --- ///////////////////////// --- //
// --------------------------------- //


RCPP_MODULE(DUSTMODULE1D)
{
  class_<Skeleton1D>("Skeleton1D")

    .constructor()

    .factory<const std::string&, const std::string&, double>(newModule1D)

    .method("fit_raw", &Skeleton1D::fit)
    .method("compute", &Skeleton1D::compute)
    .method("get_partition", &Skeleton1D::get_partition)
    .method("quick_raw", &Skeleton1D::quick)

  ;
}
