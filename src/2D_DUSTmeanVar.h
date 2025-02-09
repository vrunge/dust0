#ifndef DUST_meanVar_H
#define DUST_meanVar_H

#include <Rcpp.h>
#include <random> /// FOR RANDOM NUMBER IN DUAL EVAL

#include "2D_B_Indices2.h"

using namespace Rcpp;


class DUST_meanVar
{

  ////////////////////////////////
  ////////////////////////////////
  ////////////////////////////////

public:
  DUST_meanVar(int dual_max_type,
               int constraint_indices,
               Nullable<int> nbLoops = Nullable<int>());

  ~DUST_meanVar();

  // --- // Setup // --- //
  // fit is accessible by user
  void init(std::vector<double>& inData, Nullable<double> inPenalty = Nullable<double>());

  // --- // Main computation // --- //
  void compute(std::vector<double>& inData);

  // --- // Result retrieval // --- //
  // get_partition is accessible by user
  List get_partition();

  // --- // Wrapper method for quick use of the class // --- //
  // quick is accessible by user
  List quick(std::vector<double>& inData, Nullable<double> inPenalty = Nullable<double>());

  ////////////////////////////////
  ////////////////////////////////
  ////////////////////////////////

  private:
  const double phi = (1 + sqrt(5)) / 2;  // Golden ratio
  const double m1 = 0.01;  // Armijo

  std::vector<double> cumsum;
  std::vector<double> cumsum2;
  std::vector<double> costRecord;
  int nb_Loops; // number of loops in optimization step (For dual max)

  double Cost(unsigned int t, unsigned int s) const;
  double dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const;
  double dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const;

  double muMax(double a, double b, double a2, double b2) const;


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
  bool dualMaxAlgo0(double minCost, unsigned int t, unsigned int s, unsigned int r);
  bool dualMaxAlgo1(double minCost, unsigned int t, unsigned int s, unsigned int r);
  bool dualMaxAlgo2(double minCost, unsigned int t, unsigned int s, unsigned int r);
  bool dualMaxAlgo3(double minCost, unsigned int t, unsigned int s, unsigned int r);
  bool dualMaxAlgo4(double minCost, unsigned int t, unsigned int s, unsigned int r);
  bool dualMaxAlgo5(double minCost, unsigned int t, unsigned int s, unsigned int r);
  bool dualMaxAlgo6(double minCost, unsigned int t, unsigned int s, unsigned int r);

  bool (DUST_meanVar::*current_test)(double minCost, unsigned int t, unsigned int s, unsigned int r);

  // --- // Result processing // --- //
  std::forward_list<unsigned int> backtrack_changepoints();

  // --- // Private fields // --- //
  int dual_max_type;
  int constraint_indices;

  Indices_2D2* indices;
  std::vector<int> nb_indices;

  unsigned int n; // number of observations
  double penalty;

  std::vector<int> changepointRecord;
};

#endif
