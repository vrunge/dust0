#ifndef Negbin_1D_H
#define Negbin_1D_H

#include <Rcpp.h>

#include "1D_A_DUST.h"

using namespace Rcpp;

class Negbin_1D : public DUST_1D {
public:
  Negbin_1D(int dual_max, bool random_constraint, Nullable<double> alpha = Nullable<double>(), Nullable<int> nbLoops = Nullable<int>());
protected:
  double Cost(unsigned int t, unsigned int s) const override;
  double dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const override;
  double dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const override;

  double muMax(double a, double b) const override;

  double Dstar(double x) const override;
  double DstarPrime(double x) const override;
  double DstarSecond(double x) const override;
};

#endif
