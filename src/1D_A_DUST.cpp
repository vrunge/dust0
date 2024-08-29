#include <Rcpp.h>
#include <cmath>

#include <random> /// FOR RANDOM NUMBER IN DUAL EVAL

#include "1D_A_DUST.h"
#include "preProcessing.h"

using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// --- // Constructor // --- //
DUST_1D::DUST_1D(bool use_dual_max, bool random_constraint, Nullable<double> alpha_, Nullable<int> nbLoops)
  : use_dual_max(use_dual_max),
    random_constraint(random_constraint),
    indices(nullptr)
{
  if(alpha_.isNull())
  {
    alpha = 1e-9;
  }
  else
  {
    alpha = as<double>(alpha_);
  }
  if(nbLoops.isNull())
  {
    nb_Loops = 10;
  }
  else
  {
    nb_Loops = as<int>(nbLoops);
  }
}

DUST_1D::~DUST_1D()
{
  delete indices;
}

void DUST_1D::init_method()
{
  delete indices;
  if(random_constraint)
  {
    indices = new RandomIndices(n, alpha);
  }
  else
  {
    indices = new DeterministicIndices;
  }

  if(use_dual_max)
  {
    current_test = &DUST_1D::exact_test;
  }
  else
  {
    current_test = &DUST_1D::random_test;
  }

  engine.seed(std::random_device{}());
  dist = std::uniform_real_distribution<double>(0.0, 1.0);

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// --- // Fits the data, i. e. initializes all data-dependent vectors // --- //
void DUST_1D::init(NumericVector& inData, Nullable<double> inPenalty)
{
  data = std::move(inData);
  n = data.size();

  if (inPenalty.isNull())
  {
    penalty = 2 * pow(sdDiff(data), 2) * std::log(n); //to do
  }
  else
  {
    penalty = as<double>(inPenalty);
  }

  changepointRecord = IntegerVector(n + 1, 0);

  //cumsum = NumericVector(n + 1, 0.);
  //costRecord = NumericVector(n + 1, -penalty);
  cumsum = std::vector<double>(n + 1, 0.);
  costRecord = std::vector<double>(n + 1, -penalty);


  init_method();

  indices->add(0);
  indices->add(1);

}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// --- // Algorithm-specific method // --- //
void DUST_1D::compute()
{
  // Initialize OP step value
  double lastCost; // temporarily stores the cost for the model with last changepoint at some i
                   // then keeps the cost of the model with last changepoint at the first possible index in the t-th OP step ...
                   // ... storing it allows pruning of the first available index
  double minCost;
  unsigned int argMin; // stores the optimal last changepoint for the current OP step

  // First OP step (t = 1)
  unsigned int t = 1;
  unsigned int s = 0;
  cumsum[1] = data[0];
  costRecord[1] = Cost(t, s);
  changepointRecord[1] = 0;

  // Main loop
  for (t = 2; t <= n; t++)
  {
    // update cumsum
    cumsum[t] =
      cumsum[t - 1] + data[t - 1];

    // OP step
    indices->reset();
    minCost = std::numeric_limits<double>::infinity();
    do
    {
      s = indices->get_current();
      lastCost = costRecord[s] + Cost(t, s);
      if (lastCost < minCost)
      {
        minCost = lastCost;
        argMin = s;
      }
      indices->next();
    }
    while(indices->check());
    // END (OP step)

    // OP update
    minCost += penalty;
    costRecord[t] = minCost;
    changepointRecord[t] = argMin;

    // DUST step
    indices->reset_prune();

    // DUST loop
    while (indices->check_prune())
    {
      if ((this->*current_test)(minCost, t, indices->get_current(), indices->get_constraint()) > 0) // prune as needs pruning
      {
        // remove the pruned index and its pointer
        // removing the elements increments the cursors i and pointerIt, while before stands still
        indices->prune_current();
      }
      else
      {
        // increment all cursors
        indices->next_prune();
      }
    }
    // END (DUST loop)

    // Prune the last index (analoguous with a null (mu* = 0) duality simple test)
    if (lastCost > minCost)
    {
      indices->prune_last();
    }

    // update the available indices
    indices->add(t);
  }
}

// --- // Test methods // --- //
double DUST_1D::exact_test(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  return dualMax(minCost, t, s, r);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


double DUST_1D::random_test(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  return dualEval(dist(engine), minCost, t, s, r);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


// --- // Builds changepoints // --- //
std::forward_list<unsigned int> DUST_1D::backtrack_changepoints()
{
  std::forward_list<unsigned int> changepoints {n};
  for (int newChangepoint = changepointRecord[n]; newChangepoint != 0; newChangepoint = changepointRecord[newChangepoint])
  {
    changepoints.push_front(newChangepoint);
  }
  return changepoints;
}

// --- // Retrieves optimal partition // --- //
List DUST_1D::get_partition()
{
  costRecord.erase(costRecord.begin()); ///// REMOVE FIRST ELEMENT /////
  return List::create(
    _["changepoints"] = backtrack_changepoints(),
    _["lastIndexSet"] = indices->get_list(),
    _["costQ"] = costRecord
  );
}

// --- // Wrapper method for quickly computing               // --- //
// --- // and retrieving the optimal partition of input data // --- //
List DUST_1D::quick(NumericVector& inData, Nullable<double> inPenalty)
{
  init(inData, inPenalty);
  compute();
  return get_partition();
}
