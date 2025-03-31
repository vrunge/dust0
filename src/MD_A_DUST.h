#ifndef DUST_MD_H
#define DUST_MD_H

#include <RcppArmadillo.h>
#include <random>
#include <forward_list>
// #include <fstream>

#include "MD_B_Indices.h"
#include "utils.h"
#include <array>


using namespace Rcpp;

class DUST_MD
{
public:
  DUST_MD(int dual_max_type,
          int constraints_type,
          Nullable<unsigned> nbLoops = Nullable<unsigned>());

  virtual ~DUST_MD();

  // --- // Setup // --- //
  // append is accessible by user
  void append(const arma::dmat& inData,
            Nullable<double> inPenalty = Nullable<double>(),
            Nullable<unsigned int> inNbL = Nullable<unsigned int>(),
            Nullable<unsigned int> inNbR = Nullable<unsigned int>());

  // --- // Main computation accessible by user
  void update_partition(); ///

  // --- // Result retrieval // --- //
  // get_partition is accessible by user
  List get_partition();
  List get_info();

  // --- // Wrapper method for quick use of the class // --- //
  // one_dust is accessible by user
  List one_dust(const arma::dmat& inData,
             Nullable<double> inPenalty = Nullable<double>(),
             Nullable<unsigned int> inNbL = Nullable<unsigned int>(),
             Nullable<unsigned int> inNbR = Nullable<unsigned int>());

  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////

protected:
  unsigned int n; // number of observations
  unsigned int dim; // data dimension
  unsigned int nb_max; // = nb_l + nb_r
  unsigned int nb_l; // number of constraints before the current index to prune
  unsigned int nb_r; // number of constraints after the current index to prune

  const double phi = (1 + sqrt(5)) / 2;  // Golden ratio
  const double m1 = 0.01;  // Armijo
  arma::dmat cumsum;
  std::vector<double> costRecord;
  unsigned nb_Loops; // number of loops in optimization step (For dual max)


  ///////////////////////////////////////////////////////////////////////////
  ///
  /// VIRTUAL START
  ///
  /// computation on data
  virtual double Cost(const unsigned int& t, const unsigned int& s) const = 0;
  virtual double statistic(const double& value) const = 0;

  /// dual domain bounds
  virtual double muMax(const double& a, const double& b) const = 0;
  virtual std::array<double, 2> muInterval(const arma::colvec& a, const arma::colvec& b, double& c, double& d) const = 0;
  virtual void clipStepSizeModel(const double& m_elem, const arma::rowvec& constraint_means, const double& mu_sum, const arma::rowvec& direction, const double& direction_sum, double& max_stepsize) const = 0;

  /// dual evaluation
  virtual double dual1D_Eval(double& point, const arma::colvec& a, const arma::colvec& b, double& c, double& d, double& e, double& f) const = 0;
  virtual double dual1D_Max(double& argmax, arma::colvec& a, arma::colvec& b, double& c, double& d, double& e, double& f) const = 0;

  /// dual non linear function
  virtual double Dstar(const double& x) const = 0;
  virtual double DstarPrime(const double& x) const = 0;
  virtual double DstarSecond(const double& x) const = 0;

  virtual std::string get_model() const = 0;
  ///
  /// VIRTUAL END
  ///
  ///////////////////////////////////////////////////////////////////////////

  // arma::colvec objectiveMean;
  // arma::dmat constraintMean;
  // arma::rowvec linearTerm;
  // double constantTerm;
  // arma::rowvec mu;
  double dual_Eval();
  double dual_Eval(double &nonLinear);
  void grad_Eval(const double nonLinear);
  void update_dual_parameters_l(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int>& l);
  void update_dual_parameters_l_r(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int>& l, std::vector<unsigned int>& r);

  //////////// RANDOM NUMBER GENERATOR ////////////

  std::minstd_rand0 engine;  // Random number engine
  std::uniform_real_distribution<double> dist;  // Uniform distribution [0, 1)

  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////

private:
  // --- // Test and Indices init // --- //
  void pruning_method();

  ////////////////////////////////////
  // --- // MAX DUAL METHODS // --- //
  ////////////////////////////////////
  // 0: random eval
  // 1: random direction, exact evaluation in this direction
  // 2: barycenter evaluation
  // 3: coordinate descent
  // 4: Quasi-Newton
  // 5: PELT
  // 6: OP
  bool dualMaxAlgo0(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> l, std::vector<unsigned int> r);
  bool dualMaxAlgo1(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> l, std::vector<unsigned int> r);
  bool dualMaxAlgo2(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> l, std::vector<unsigned int> r);
  bool dualMaxAlgo3(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> l, std::vector<unsigned int> r);
  bool dualMaxAlgo4(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> l, std::vector<unsigned int> r);
  bool dualMaxAlgo42(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> l, std::vector<unsigned int> r);
  bool dualMaxAlgo5(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> l, std::vector<unsigned int> r);
  bool dualMaxAlgo6(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> l, std::vector<unsigned int> r);
  bool dualMaxAlgo7(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> l, std::vector<unsigned int> r);

  bool (DUST_MD::*current_test)(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> l, std::vector<unsigned int> r);

  // --- // Result processing // --- //
  std::forward_list<unsigned int> backtrack_changepoints();

  // --- // Private fields // --- //
  int dual_max_type;
  int constraints_type;

  Indices_MD* indices;
  std::vector<int> nb_indices;
  double penalty;

  std::vector<int> changepointRecord;

  ///////////////////////////////////////
  ///
  /// PARAMETERS defining the DUAL FUNCTION
  ///
  /// link to constraints (row)
  /// vector size nb_l_r (= number of constraints)
  ///
  arma::rowvec mu; /// for the optimization step
  arma::rowvec mu_sign;
  arma::rowvec mu_max; /// constraints on mu (case no right constraint)
  /////  /////  /////  ///// link to dimension (column)
  /// vector size d = dimension
  arma::colvec objectiveMean; // mean between s and t
  arma::dmat constraintMean; // mean between all r and s   /// size d x nb_l_r
  arma::rowvec linearTerm; /// sum (pm 1)mu Delta Q
  double constantTerm;
  ///
  ///
  ///
  ///////////////////////////////////////


  arma::colvec m_value; // m(mu)
  arma::rowvec inv_max;  /// 1/max
  arma::rowvec grad;   /// grad in the optimization, of the dual
  arma::colvec nonLinearGrad; // Dstar'(m(mu))

  arma::dmat Identity; /// size nb_l_r x nb_l_r
  arma::dmat inverseHessian; // for quasi newton optimizing /// size nb_l_r x nb_l_r

};

#endif
