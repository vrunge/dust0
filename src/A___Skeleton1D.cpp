// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "A___Skeleton1D.h"
#include "preProcessing.h"

using namespace Rcpp;


// --- // Constructor // --- //
Skeleton1D::Skeleton1D() : validIndices(nullptr) {}

// --- // Default destructor // --- //
Skeleton1D::~Skeleton1D() {}

// --- // Indices-handler-specific method // --- //
void Skeleton1D::init_handler(){}

// --- // Test-method-specific method // --- //
double Skeleton1D::test(int t, int i, int j, double testThreshold) {return -1;}

// --- // Model-specific methods // --- //
double Skeleton1D::modelCost(int t, int i) {return -1;}
double Skeleton1D::modelTest(int t, int i, int j, double testThreshold) {return -1;}
double Skeleton1D::modelDual(double point, int t, int i, int j, double testThreshold) {return -1;}

// --- // Takes indices handler pointer and embeds it into protected validIndces // --- //
void Skeleton1D::set_method(std::unique_ptr<IndicesModule> validIndices_)
{
  validIndices = std::move(validIndices_);
}


// --- // Fits the data, i. e. initializes all data-dependent vectors // --- //
void Skeleton1D::fit(NumericVector inData, Nullable<double> inPenalty)
{
  data = inData;
  n = data.size();
  
  if (inPenalty.isNull()) {
    penalty = 2 * pow(madEstimator(data), 2) * log(n);
  }
  else
  {
    penalty = as<double>(inPenalty);
  }
  
  changepointsForward = IntegerVector(n + 1, 0);
  valuesCumsum = NumericVector(n + 1, 0.);
  costRecord = NumericVector(n + 1, -penalty);
  changepoints = std::forward_list<int>{n};
  
  init_handler();
  
  validIndices->add(0);
  validIndices->add(1);
}

// --- // Computes the optimal partition of the data // --- //
void Skeleton1D::compute()
{
  // Initialize OP step value
  double lastCost; // temporarily stores the cost for the model with last changepoint at some i, then keeps the cost of the model with last changepoint at the first possible index in the t-th OP step ...
  // ... storing it allows pruning of the first available index
  int optimalChangepoint; // stores the optimal last changepoint for the current OP step
  double optimalCost; // stores the cost for the model with optimal last changepoint for the current OP step
  
  
  // Initialize pruning step values and vectors
  int i;
  double testValue; // the value to be checked vs. the test threshold = optimalCost
  
  
  // First OP step (t = 1)
  valuesCumsum[1] = data[0];
  costRecord[1] = modelCost(1, 0);
  changepointsForward[1] = 0;
  
  
  // Main loop
  for (int t = 2; t <= n; t++)
  {
    // update valuesCumsum
    valuesCumsum[t] =
      valuesCumsum[t - 1] + data[t - 1];
    
    // OP step
    validIndices->reset();
    optimalCost = std::numeric_limits<double>::infinity();
    do
    {
      i = validIndices->get_current();
      lastCost = costRecord[i] + modelCost(t, i);
      if (lastCost < optimalCost)
      {
        optimalCost = lastCost;
        optimalChangepoint = i;
      }
      validIndices->next();
    }
    while(validIndices->check());
    // END (OP step)
    
    // OP update
    optimalCost += penalty;
    costRecord[t] = optimalCost;
    changepointsForward[t] = optimalChangepoint;
    
    // DUST step
    validIndices->reset_prune();

    // DUST loop
    while (validIndices->check_prune())
    {
      testValue = test(t, validIndices->get_current(), validIndices->get_constraint(), optimalCost); // compute test value
      if (testValue > 0) // prune as needs pruning
      {
        // remove the pruned index and its pointer
        // removing the elements increments the cursors i and pointerIt, while before stands still
        validIndices->prune_current();
      }
      else
      {
        // increment all cursors
        validIndices->next_prune();
      }
      validIndices->new_constraint();
    }
    // while (validIndices->check_prune()); // exit the loop if we may not draw a valid constraint index
    // END (DUST loop)

    // Prune the last index (analoguous with a null (mu* = 0) duality simple test)
    if (lastCost > optimalCost) {
      validIndices->prune_current();
    }
    
    // update the available indices
    validIndices->add(t);
  }
}

// --- // Builds changepoints // --- //
void Skeleton1D::backtrack_changepoints()
{
  for (int newChangepoint = changepointsForward[n]; newChangepoint != 0; newChangepoint = changepointsForward[newChangepoint])
  {
    changepoints.push_front(newChangepoint);
  }
}

// --- // Retrieves optimal partition // --- //
List Skeleton1D::get_partition()
{
  backtrack_changepoints();
  
  return List::create(
    _["changepoints"] = this->changepoints,
    _["lastIndexSet"] = this->validIndices->get_list(),
    _["costQ"] = this->costRecord
  );
}

// --- // Wrapper method for quickly computing               // --- //
// --- // and retrieving the optimal partition of input data // --- //
List Skeleton1D::quick(NumericVector inData, Nullable<double> inPenalty)
{
  fit(inData, inPenalty);
  compute();
  return get_partition();
}
