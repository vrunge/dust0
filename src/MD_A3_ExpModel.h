#ifndef Exp_MD_H
#define Exp_MD_H

#include "MD_A_DUST.h"

using namespace Rcpp;

class Exp_MD : public DUST_MD {
public:
  Exp_MD(int dual_max_type, int constraints_type, Nullable<unsigned> nbLoops = Nullable<unsigned>());
  protected:
  double Cost(const unsigned int& t, const unsigned int& s) const override;
  double statistic(const double& value) const override;

  double muMax(const double& a, const double& b) const override;
  std::array<double, 2> muInterval(const arma::colvec& a, const arma::colvec& b, double& c, double& d) const override;
  void clipStepSizeModel(const double& m_elem, const arma::rowvec& constraint_means, const double& mu_sum, const arma::rowvec& direction, const double& direction_sum, double& max_stepsize) const override;

  virtual double dual1D_Eval(double& point, const arma::colvec& a, const arma::colvec& b, double& c, double& d, double& e, double& f) const override;
  virtual double dual1D_Max(double& argmax, arma::colvec& a, arma::colvec& b, double& c, double& d, double& e, double& f) const override;

  double Dstar(const double& x) const override;
  double DstarPrime(const double& x) const override;
  double DstarSecond(const double& x) const override;

  std::string get_model() const override;
};
#endif

