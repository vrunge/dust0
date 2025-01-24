#ifndef Gauss_MD_H
#define Gauss_MD_H

#include "MD_A_DUST.h"

using namespace Rcpp;

class Gauss_MD : public DUST_MD {
public:
  Gauss_MD(int dual_max, bool random_constraint, Nullable<int> nbLoops = Nullable<int>());
protected:
  double Cost(const unsigned int& t, const unsigned int& s) const override;
  double statistic(const double& value) const override;

  double muMax(const double& a, const double& b) const override;
  std::array<double, 2> muInterval(const arma::colvec& a, const arma::colvec& b, double& c, double& d) const override;
  void clipStepSizeModel(const double& m_elem, const arma::rowvec& constraint_means, const double& mu_sum, const arma::rowvec& direction, const double& direction_sum, double& max_stepsize) const override;

  double dualMax(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> r) override { return 0; }
  double dualEval(std::vector<unsigned int> point, const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> r) override { return 0; }
  virtual std::array<double, 2> dual1D_ArgmaxMax(arma::colvec& a, arma::colvec& b, double& c, double& d, double& e, double& f) const override;

  double Dstar(const double& x) const override;
  double DstarPrime(const double& x) const override;
  double DstarSecond(const double& x) const override;

  std::string get_model() const override;
};

#endif
