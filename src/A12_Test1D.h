// [[Rcpp::depends(Rcpp)]]

#ifndef TEST1D_H
#define TEST1D_H

#include <Rcpp.h>
#include "A___Skeleton1D.h"

using namespace Rcpp;


class RandomTest1D : public virtual Skeleton1D {
public:
  double test(int t, int i, int j, double testThreshold) override;
};

class DeterministicTest1D : public virtual Skeleton1D {
public:
  double test(int t, int i, int j, double testThreshold) override;
};

#endif
