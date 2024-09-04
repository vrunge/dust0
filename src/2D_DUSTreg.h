#ifndef DUST_reg_H
#define DUST_reg_H

#include <Rcpp.h>
#include <random> /// FOR RANDOM NUMBER IN DUAL EVAL

#include "1D_B_Indices.h"

using namespace Rcpp;


class DUST_reg
{

  ////////////////////////////////
  ////////////////////////////////
  ////////////////////////////////

public:
  DUST_reg(int dual_max,
           bool random_constraint,
           Nullable<double> alpha = Nullable<double>(),
           Nullable<int> nbLoops = Nullable<int>());

  ~DUST_reg();

  // --- // Setup // --- //
  // fit is accessible by user
  void init(DataFrame& inData, Nullable<double> inPenalty = Nullable<double>());

  // --- // Main computation // --- //
  void compute();

  // --- // Result retrieval // --- //
  // get_partition is accessible by user
  List get_partition();

  // --- // Wrapper method for quick use of the class // --- //
  // quick is accessible by user
  List quick(DataFrame& inData, Nullable<double> inPenalty = Nullable<double>());

  ////////////////////////////////
  ////////////////////////////////
  ////////////////////////////////

  private:

  const double phi = (1 + sqrt(5)) / 2;  // Golden ratio
  std::vector<double> A;
  std::vector<double> B;
  std::vector<double> C;
  std::vector<double> D;
  std::vector<double> E;
  std::vector<double> F;

  std::vector<double> costRecord;
  int nb_Loops; // number of loops in optimization step (For dual max)

  double Cost(unsigned int t, unsigned int s) const;
  double dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const;
  double dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const;

  double Dstar(double x) const;
  double DstarPrime(double x) const;
  double DstarSecond(double x) const;

  //////////// RANDOM NUMBER GENERATOR ////////////

  std::minstd_rand0 engine;  // Random number engine
  std::uniform_real_distribution<double> dist;  // Uniform distribution [0, 1)

  ////////////////////////////////
  ////////////////////////////////
  ////////////////////////////////

  // --- // Test and Indices init // --- //
  void init_method();

  // --- // MAX DUAL METHODS // --- //
  // --- //   // --- //   // --- //   // --- //
  double dualMaxAlgo0(double minCost, unsigned int t, unsigned int s, unsigned int r);
  double dualMaxAlgo1(double minCost, unsigned int t, unsigned int s, unsigned int r);
  double dualMaxAlgo2(double minCost, unsigned int t, unsigned int s, unsigned int r);
  double dualMaxAlgo3(double minCost, unsigned int t, unsigned int s, unsigned int r);
  double dualMaxAlgo4(double minCost, unsigned int t, unsigned int s, unsigned int r);
  double dualMaxAlgo5(double minCost, unsigned int t, unsigned int s, unsigned int r);

  double (DUST_reg::*current_test)(double minCost, unsigned int t, unsigned int s, unsigned int r);

  // --- // Result processing // --- //
  std::forward_list<unsigned int> backtrack_changepoints();

  // --- // Private fields // --- //
  int dual_max;
  bool random_constraint;
  double alpha;

  Indices* indices;
  std::vector<int> nb_indices;

  unsigned int n; // number of observations
  NumericVector data_x;
  NumericVector data_y;
  double penalty;

  std::vector<int> changepointRecord;

};

#endif
