// [[Rcpp::depends(Rcpp)]]

#ifndef SKELETON1D_H
#define SKELETON1D_H

#include <Rcpp.h>
#include "B___IndicesModule.h"

using namespace Rcpp;


class Skeleton1D {
protected:
  std::unique_ptr<IndicesModule> validIndices;

public:

  Skeleton1D();
  virtual ~Skeleton1D();

  // --- // Setup // --- //
  // fit is accessible by user
  void set_method(std::unique_ptr<IndicesModule> validIndices_);
  void fit(NumericVector inData, Nullable<double> inPenalty = Nullable<double>());

  // --- // Main computation // --- //
  void compute();

  // --- // Method specific test // --- //
  virtual void init_handler();
  virtual double test(int t, int i, int j, double testThreshold);

  // --- // Model specific computations // --- //
  virtual double modelCost(int t, int i);
  virtual double modelTest(int t, int i, int j, double testThreshold);
  virtual double modelDual(double point, int t, int i, int j, double testThreshold);

  // --- // Result processing // --- //
  void backtrack_changepoints();

  // --- // Result retrieval // --- //
  // get_partition is accessible by user
  List get_partition();

  // --- // Wrapper method for quick use of the class // --- //
  // quick is accessible by user
  List quick(NumericVector inData, Nullable<double> inPenalty = Nullable<double>());

  // ------------------------------------ //
  // ------------------------------------ //
  // ------------------------------------ //

  // --- // Public fields // --- //
  NumericVector valuesCumsum;
  NumericVector costRecord;
  int n;

private:
  double penalty;

  NumericVector data;

  std::forward_list<int> changepoints;
  IntegerVector changepointsForward; // trouvefr nomo cmmun avec Q, "online", "record" ?
};

#endif
