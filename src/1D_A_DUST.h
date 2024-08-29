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
  DUST_1D(bool use_dual_max,
          bool random_constraint,
          Nullable<double> alpha = Nullable<double>(),
          Nullable<int> nbLoops = Nullable<int>());

  virtual ~DUST_1D();

  // --- // Setup // --- //
  // fit is accessible by user
  void init(NumericVector& inData, Nullable<double> inPenalty = Nullable<double>());

  // --- // Main computation // --- //
  void compute();

  // --- // Result retrieval // --- //
  // get_partition is accessible by user
  List get_partition();

  // --- // Wrapper method for quick use of the class // --- //
  // quick is accessible by user
  List quick(NumericVector& inData, Nullable<double> inPenalty = Nullable<double>());

  ////////////////////////////////
  ////////////////////////////////
  ////////////////////////////////

protected:
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

  // --- // Test handling // --- //
  double exact_test(double minCost, unsigned int t, unsigned int s, unsigned int r);
  double random_test(double minCost, unsigned int t, unsigned int s, unsigned int r);

  double (DUST_1D::*current_test)(double minCost, unsigned int t, unsigned int s, unsigned int r);

  // --- // Result processing // --- //
  std::forward_list<unsigned int> backtrack_changepoints();

  // --- // Private fields // --- //
  bool use_dual_max;
  bool random_constraint;
  double alpha;

  Indices* indices;

  unsigned int n; // number of observations
  NumericVector data;
  double penalty;

  IntegerVector changepointRecord;

};

#endif
