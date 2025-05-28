#ifndef Gauss_1D_H
#define Gauss_1D_H

#include <Rcpp.h>

#include "1D_A_DUST.h"

using namespace Rcpp;

class Gauss_1D : public DUST_1D {
public:
  Gauss_1D(int dual_max_type, int constraints_type, Nullable<int> nbLoops = Nullable<int>());
protected:
  double Cost(unsigned int t, unsigned int s) const override;
  double statistic(double& data) const override;

  double dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const override;
  double dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const override;

  double muMax(double a, double b) const override;
  bool isBoundary(double a) const override;

  double Dstar(double x) const override;
  double DstarPrime(double x) const override;
  double DstarPrimeInv(double x) const override;
  double DstarSecond(double x) const override;

  std::string get_model() const override;
};

#endif
