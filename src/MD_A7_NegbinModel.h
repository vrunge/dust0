#ifndef Negbin_MD_H
#define Negbin_MD_H

#include "MD_A_DUST.h"

using namespace Rcpp;

class Negbin_MD : public DUST_MD {
public:
  Negbin_MD(int dual_max, bool random_constraint, Nullable<int> nbLoops = Nullable<int>());
protected:
  double Cost(const unsigned int& t, const unsigned int& s) const override;
  double statistic(const double& value) const override;

  double muMax(const double& a, const double& b) const override;

  double dualMax(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> r) override { return 0; }
  double dualEval(std::vector<unsigned int> point, const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> r) override { return 0; }
  double dualMax(arma::colvec& a, arma::colvec& b, double& c, double& d, double& e, double& f) const override;

  double Dstar(const double& x) const override;
  double DstarPrime(const double& x) const override;
  double DstarSecond(const double& x) const override;

  std::string get_model() const override;
};

#endif

