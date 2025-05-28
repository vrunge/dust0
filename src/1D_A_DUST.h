#ifndef DUST_1D_H
#define DUST_1D_H

#include <Rcpp.h>
#include <random> /// FOR RANDOM NUMBER IN DUAL EVAL

#include "1D_B_Indices.h"

using namespace Rcpp;


class DUST_1D
{

public:
  DUST_1D(int dual_max_type,
          int constraints_type,
          Nullable<int> nbLoops = Nullable<int>());

  virtual ~DUST_1D();

  // --- // Setup // --- //
  // fit is accessible by user
  void append(std::vector<double>& inData, Nullable<double> inPenalty = Nullable<double>());

  // --- // Main computation // --- //
  void update_partition();

  // --- // Result retrieval // --- //
  // get_partition is accessible by user
  List get_partition();

  List get_info();

  // --- // Wrapper method for quick use of the class // --- //
  // quick is accessible by user
  List one_dust(std::vector<double>& inData, Nullable<double> inPenalty = Nullable<double>());

  ////////////////////////////////
  ////////////////////////////////
  ////////////////////////////////

protected:
  const double phi = (1 + sqrt(5)) / 2;  // Golden ratio
  const double m1 = 0.01;  // Armijo
  std::vector<double> cumsum;
  std::vector<double> costRecord;
  int nb_Loops; // number of loops in optimization step (For dual max)

  virtual double Cost(unsigned int t, unsigned int s) const = 0;
  virtual double statistic(double& data) const = 0;

  virtual double dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const = 0;
  virtual double dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const = 0;

  virtual double muMax(double a, double b) const = 0;
  virtual bool isBoundary(double a) const = 0;

  virtual double Dstar(double x) const = 0;
  virtual double DstarPrime(double x) const = 0;
  virtual double DstarPrimeInv(double x) const = 0;
  virtual double DstarSecond(double x) const = 0;

  virtual std::string get_model() const = 0;

  //////////// RANDOM NUMBER GENERATOR ////////////

  std::minstd_rand0 engine;  // Random number engine
  std::uniform_real_distribution<double> dist;  // Uniform distribution [0, 1)

  ////////////////////////////////
  ////////////////////////////////
  ////////////////////////////////

private:
  // --- // Test and Indices init // --- //
  void pruning_method();

  bool isSpecialCase(double objectiveMean, double constraintMean);
  bool specialCasePruning(double objectiveMean, double constraintMean,
                          double linearTerm, double constantTerm, double mu_max);

  // --- // MAX DUAL METHODS // --- //
  // --- //   // --- //   // --- //   // --- //
  // 0: random eval
  // 1: exact eval
  // 2: golden-section search
  // 3: binary search. At each step, we evaluate the tangent line to the current point at its max to stop the search at early step (when possible)
  // 4: Quasi-Newton
  // 5: PELT
  // 6: OP
  bool dualMaxAlgo0(double minCost, unsigned int t, unsigned int s, unsigned int r);
  bool dualMaxAlgo1(double minCost, unsigned int t, unsigned int s, unsigned int r);
  bool dualMaxAlgo2(double minCost, unsigned int t, unsigned int s, unsigned int r);
  bool dualMaxAlgo3(double minCost, unsigned int t, unsigned int s, unsigned int r);
  bool dualMaxAlgo4(double minCost, unsigned int t, unsigned int s, unsigned int r);
  bool dualMaxAlgo5(double minCost, unsigned int t, unsigned int s, unsigned int r);
  bool dualMaxAlgo6(double minCost, unsigned int t, unsigned int s, unsigned int r);

  bool (DUST_1D::*current_test)(double minCost, unsigned int t, unsigned int s, unsigned int r);

  // --- // Result processing // --- //
  std::forward_list<unsigned int> backtrack_changepoints();

  // --- // Private fields // --- //
  int dual_max_type;
  int constraints_type;

  Indices_1D* indices;
  std::vector<int> nb_indices;

  unsigned int n; // number of observations
  double penalty;

  std::vector<int> changepointRecord;

};

#endif
