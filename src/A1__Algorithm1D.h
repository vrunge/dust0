// [[Rcpp::depends(RcppArmadillo)]]

#ifndef ALGORITHM1D_H
#define ALGORITHM1D_H

#include <RcppArmadillo.h>

#include "A___Skeleton1D.h"
#include "A11_IndicesHandler1D.h"
#include "A12_Test1D.h"

using namespace Rcpp;


// --- // Define variations of the DUST algorithm           // --- //
// --- // by combining an indices_handler and a test method // --- //

class AllRandom1D :
  public virtual RandomConstraint1D,
  public virtual RandomTest1D
{
public:
  AllRandom1D(double alpha_);
};

class HalfRandom1D :
  public virtual RandomConstraint1D,
  public virtual DeterministicTest1D
{
public:
  HalfRandom1D(double alpha_);
};

class Deterministic1D :
  public virtual DeterministicConstraint1D,
  public virtual DeterministicTest1D
{};

#endif