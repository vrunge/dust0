// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "costFunctions.h"
#include "RandomList.h"
#include <forward_list>
#include <cmath>
#include <random>
#include <limits>
using namespace Rcpp;


// --------- // DUST // --------- //
//
// Provided some data vector, uses the DUST algorithm to return its optimal par-
// titioning
//
// Parameters:
//  - data (vector): a vector of numeric values
//  - penalty: the value of the penalty in the penalized changepoint detection
//    model
//  - alpha: value controlling the size of the random vector used to select the
//    constraint index. the smaller the value, the longer the random vector,
//    inducing a behaviour closer to "true" randomness
//
// Rd file [TODO]


// [[Rcpp::export]]
List dust1DRandConstraint(NumericVector data, double penalty = 0, double alpha = 1e-9) {
  const int n = data.size();

  if(penalty == 0)
    penalty = 2 * log(n);


  // Initialize incremented vectors

  IntegerVector changepointsForward(n + 1, 0); // changepointsForward records the optimal last change point at each OP step
  NumericVector valuesCumsum(n + 1, 0.0), costRecord(n + 1, - penalty); // valuesCumsum stores the cumsum of the data and costRecord contains the optimal model cost at each OP step


  // Initialize OP step values

  double lastCost; // temporarily stores the cost for the model with last changepoint at some i, then keeps the cost of the model with last changepoint at the first possible index in the t-th OP step ...
  // ... storing it allows pruning of the first available index
  int optimalChangepoint; // stores the optimal last changepoint for the current OP step
  double optimalCost; // stores the cost for the model with optimal last changepoint for the current OP step

  // Initialize pruning step values and vectors

  int i;
  double testValue; // the value to be checked vs. the test threshold = optimalCost

  RandomList validIndices(n, alpha); // the available indices (decreasing)
  validIndices.add(0);
  validIndices.add(1);


  // First OP step (t = 1)

  valuesCumsum[1] = data[0];
  costRecord[1] = - pow(data[0], 2);
  changepointsForward[1] = 0;


  // Main loop

  for (int t = 2; t <= n; t++)
  {
    // update valuesCumsum
    valuesCumsum[t] =
      valuesCumsum[t - 1] + data[t - 1];

    // OP step
    validIndices.reset();
    optimalCost = std::numeric_limits<double>::infinity();
    do
    {
      i = validIndices.read();
      lastCost = costRecord[i] + modelCost(t, i, valuesCumsum);
      if (lastCost < optimalCost)
      {
        optimalCost = lastCost;
        optimalChangepoint = i;
      }
      validIndices.next();
    }
    while(validIndices.check());
    // END (OP step)

    // OP update
    optimalCost += penalty;
    costRecord[t] = optimalCost;
    changepointsForward[t] = optimalChangepoint;

    // if (t % 5) {
    //   validIndices.add(t);
    //   continue;
    // }

    // DUST step
    validIndices.reset_prune();

    // DUST loop
    do
    {
      testValue = test(t, validIndices.read(), *validIndices.draw(), valuesCumsum, costRecord); // compute test value
      if (testValue > optimalCost) // prune as needs pruning
      {
        // remove the pruned index and its pointer
        // removing the elements increments the cursors i and pointerIt, while before stands still
        validIndices.prune();
      }
      else
      {
        // increment all cursors
        validIndices.next_prune();
      }
    }
    while (validIndices.check_prune()); // exit the loop if we may not draw a valid constraint index
    // END (DUST loop)

    // Prune the last index (analoguous with a null (mu* = 0) duality simple test)
    if (lastCost > optimalCost) {
      validIndices.prune();
    }

    // update the available indices
    validIndices.add(t);
  }

  // Backtrack des changepoints
  std::forward_list<int> changepoints {n};
  for (int newChangepoint = changepointsForward[n]; newChangepoint != 0; newChangepoint = changepointsForward[newChangepoint])
  {
    changepoints.push_front(newChangepoint);
  }

  // Output
  List output;
  output["changepoints"] = changepoints;
  output["lastIndexSet"] = validIndices.get_list();
  output["costQ"] = costRecord;

  return output;
}
