// [[Rcpp::depends(Rcpp)]]

#ifndef GAUSS1D_H
#define GAUSS1D_H

#include <Rcpp.h>
#include "A___Skeleton1D.h"

using namespace Rcpp;


class GaussModel1D : public virtual Skeleton1D {
public:
  double modelCost(int t, int i) override;
  double modelTest(int t, int i, int j, double testThreshold) override;
  double modelDual(double point, int t, int i, int j, double testThreshold) override;
};

#endif
