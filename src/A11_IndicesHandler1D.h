// [[Rcpp::depends(RcppArmadillo)]]

#ifndef HANDLER1D_H
#define HANDLER1D_H

#include <RcppArmadillo.h>
#include "A___Skeleton1D.h"

using namespace Rcpp;


class RandomConstraint1D : public virtual Skeleton1D {
public:
  void init_handler() override;
  double alpha;
};

class DeterministicConstraint1D : public virtual Skeleton1D {
public:
  void init_handler() override;
};

#endif
