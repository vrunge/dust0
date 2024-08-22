#ifndef DUST_1D_H
#define DUST_1D_H

#include <Rcpp.h>
#include "1D_B_Indices.h"

using namespace Rcpp;


class DUST_1D {
public:
  DUST_1D(bool use_dual_max, bool random_constraint, Nullable<double> alpha = Nullable<double>());
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

protected:
  NumericVector cumsum;
  NumericVector costRecord;
  
  virtual double Cost(int t, int s) const = 0;
  virtual double dualEval(double point, double minCost, int t, int s, int r) const = 0;
  virtual double dualMax(double minCost, int t, int s, int r) const = 0;
  
private:
  // --- // Test and Indices init // --- //
  void init_method();
  
  // --- // Test handling // --- //
  double exact_test(double minCost, int t, int s, int r);
  double random_test(double minCost, int t, int s, int r);
  double (DUST_1D::*current_test)(double minCost, int t, int s, int r);
  
  // --- // Result processing // --- //
  std::forward_list<int> backtrack_changepoints();
  
  // --- // Private fields // --- //
  bool use_dual_max;
  bool random_constraint;
  double alpha;
  
  Indices* indices;
  
  int n; // number of observations
  NumericVector data;
  double penalty;
  
  IntegerVector changepointRecord;
};

#endif
