#ifndef DUST_1D_H
#define DUST_1D_H

#include <Rcpp.h>
#include <random> /// FOR RANDOM NUMBER IN DUAL EVAL

#include "1D_B_Indices.h"

using namespace Rcpp;


class DUST_1D
{

  ////////////////////////////////
  ////////////////////////////////
  ////////////////////////////////

public:
  DUST_1D(int dual_max,
          bool random_constraint,
          Nullable<double> alpha = Nullable<double>(),
          Nullable<int> nbLoops = Nullable<int>());

  virtual ~DUST_1D();

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

protected:
  const double phi = (1 + sqrt(5)) / 2;  // Golden ratio
  std::vector<double> cumsum;
  std::vector<double> costRecord;
  int nb_Loops; // number of loops in optimization step (For dual max)

  virtual double Cost(unsigned int t, unsigned int s) const = 0;
  virtual double dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const = 0;
  virtual double dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const = 0;

  virtual double Dstar(double x) const = 0;
  virtual double DstarPrime(double x) const = 0;
  virtual double DstarSecond(double x) const = 0;

  //////////// RANDOM NUMBER GENERATOR ////////////

  std::minstd_rand0 engine;  // Random number engine
  std::uniform_real_distribution<double> dist;  // Uniform distribution [0, 1)

  ////////////////////////////////
  ////////////////////////////////
  ////////////////////////////////

private:
  // --- // Test and Indices init // --- //
  void init_method();

  // --- // MAX DUAL METHODS // --- //
  // --- //   // --- //   // --- //   // --- //
  double dualMaxAlgo0(double minCost, unsigned int t, unsigned int s, unsigned int r);
  double dualMaxAlgo1(double minCost, unsigned int t, unsigned int s, unsigned int r);
  double dualMaxAlgo2(double minCost, unsigned int t, unsigned int s, unsigned int r);
  double dualMaxAlgo3(double minCost, unsigned int t, unsigned int s, unsigned int r);

  double (DUST_1D::*current_test)(double minCost, unsigned int t, unsigned int s, unsigned int r);

  // --- // Result processing // --- //
  std::forward_list<unsigned int> backtrack_changepoints();

  // --- // Private fields // --- //
  int dual_max;
  bool random_constraint;
  double alpha;

  Indices* indices;
  std::vector<int> nb_indices;

  unsigned int n; // number of observations
  double penalty;

  std::vector<int> changepointRecord;

};

#endif
