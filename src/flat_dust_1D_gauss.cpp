#include <Rcpp.h>
#include <cmath>

#include "1D_B_Indices.h"
#include "preProcessing.h"

using namespace Rcpp;

double Cost(const unsigned int& t, const unsigned int& s, const std::vector<double>& cumsum)
{
  return - 0.5 * (cumsum[t] - cumsum[s]) * (cumsum[t] - cumsum[s]) / (t - s);
}


double DualMax(const double& minCost, const unsigned int& t, const unsigned int& s, const unsigned int& r, const std::vector<double>& cumsum, const std::vector<double>& costRecord)
{
  // Compute the optimal point on which to evaluate the duality function

  double A = (cumsum[t] - cumsum[s])/ (t - s); // m_it
  double B = (cumsum[s] - cumsum[r])/ (s - r); // m_ji
  double absAmB = std::abs(A - B);
  double sqrtB2p2C = std::sqrt(B*B + 2*(costRecord[s] - costRecord[r])/ (s - r));

  // Case 1: mu* = 0
  // deduce the following condition from the formula for mu*
  if (absAmB >= sqrtB2p2C)
    return (costRecord[s] - minCost) / (t - s)  - 0.5 * A*A;

  // Case 2: mu* > 0
  return (costRecord[s] - minCost) / (t - s) - 0.5*A*A + 0.5 * (absAmB - sqrtB2p2C)*(absAmB - sqrtB2p2C);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


List flat_DUST_1D(const std::vector<double>& inData, Nullable<double> inPenalty = R_NilValue)
{
  DeterministicIndices_1D indices;

  const unsigned int n = inData.size();

  double penalty;
  if (inPenalty.isNull()){penalty = 2 * std::log(n);}else{penalty = as<double>(inPenalty);}

  std::vector<int> changepointRecord;
  changepointRecord.reserve(n + 1);
  changepointRecord.push_back(0);

  std::vector<double> cumsum;
  cumsum.reserve(n + 1);
  cumsum.push_back(0);

  std::vector<double> costRecord;
  costRecord.reserve(n + 1);
  costRecord.push_back(-penalty);

  indices.add(0);

  // Initialize OP step value
  double lastCost; // temporarily stores the cost for the model with last changepoint at some i
  // then keeps the cost of the model with last changepoint at the first possible index in the t-th OP step ...
  // ... storing it allows pruning of the first available index
  double minCost;
  unsigned int argMin = 0; // stores the optimal last changepoint for the current OP step

  // First OP step (t = 1)
  unsigned int t = 1;
  unsigned int s = 0;
  cumsum.push_back(inData[0]);
  costRecord.push_back(Cost(t, s, cumsum));
  changepointRecord.push_back(0);

  // Main loop
  for (t = 2; t <= n; t++)
  {
    // update cumsum
    cumsum.push_back(cumsum.back() + inData[t - 1]);

    // OP step
    indices.reset();
    minCost = std::numeric_limits<double>::infinity();
    do
    {
      s = indices.get_current();
      lastCost = costRecord[s] + Cost(t, s, cumsum);
      if (lastCost < minCost)
      {
        minCost = lastCost;
        argMin = s;
      }
      indices.next();
    }
    while(indices.check());
    // END (OP step)

    // OP update
    minCost += penalty;
    costRecord.push_back(minCost);
    changepointRecord.push_back(argMin);

    // DUST step
    indices.reset_prune();

    // DUST loop
    while (indices.check_prune())
    {
      if (DualMax(minCost, t, indices.get_current(), indices.get_constraint(), cumsum, costRecord) > 0) // prune as needs pruning
      {
        // remove the pruned index and its pointer
        // removing the elements increments the cursors i and pointersCurrent, while before stands still
        indices.prune_current();
      }
      else
      {
        // increment all cursors
        indices.next_prune();
      }
    }
    // END (DUST loop)

    // Prune the last index (analoguous with a null (mu* = 0) duality simple test)
    if (lastCost > minCost)
    {
      indices.prune_last();
    }

    // update the available indices
    indices.add(t);
  }

  std::forward_list<unsigned int> changepoints {n};
  for (int newChangepoint = changepointRecord[n]; newChangepoint != 0; newChangepoint = changepointRecord[newChangepoint])
  {
    changepoints.push_front(newChangepoint);
  }

  costRecord.erase(costRecord.begin()); ///// REMOVE FIRST ELEMENT /////
  indices.remove_first(); ///// REMOVE FIRST ELEMENT /////

  return List::create(
    _["changepoints"] = changepoints,
    _["lastIndexSet"] = indices.get_list(),
    _["costQ"] = costRecord
  );
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

List flat2_DUST_1D(const std::vector<double>& inData, Nullable<double> inPenalty = R_NilValue)
{
  DeterministicIndices_1D indices;

  const unsigned int n = inData.size();

  double penalty;
  if (inPenalty.isNull()){penalty = 2 * std::log(n);}else{penalty = as<double>(inPenalty);}

  std::vector<int> changepointRecord(n + 1);
  changepointRecord[0] = 0;

  std::vector<double> cumsum(n + 1);
  cumsum[0] = 0;

  std::vector<double> costRecord(n + 1);
  costRecord[0] = -penalty;

  indices.add(0);

  // Initialize OP step value
  double lastCost; // temporarily stores the cost for the model with last changepoint at some i
  // then keeps the cost of the model with last changepoint at the first possible index in the t-th OP step ...
  // ... storing it allows pruning of the first available index
  double minCost;
  unsigned int argMin = 0; // stores the optimal last changepoint for the current OP step

  // First OP step (t = 1)
  unsigned int t = 1;
  unsigned int s = 0;
  cumsum[1] = inData[0];
  costRecord[1] = Cost(t, s, cumsum);
  changepointRecord[1] = 0;

  // Main loop
  for (t = 2; t <= n; t++)
  {
    // update cumsum
    cumsum[t] = cumsum[t - 1] + inData[t - 1];

    // OP step
    indices.reset();
    minCost = std::numeric_limits<double>::infinity();
    do
    {
      s = indices.get_current();
      lastCost = costRecord[s] + Cost(t, s, cumsum);
      if (lastCost < minCost)
      {
        minCost = lastCost;
        argMin = s;
      }
      indices.next();
    }
    while(indices.check());
    // END (OP step)

    // OP update
    minCost += penalty;
    costRecord[t] = minCost;
    changepointRecord[t] = argMin;

    // DUST step
    indices.reset_prune();

    // DUST loop
    while (indices.check_prune())
    {
      if (DualMax(minCost, t, indices.get_current(), indices.get_constraint(), cumsum, costRecord) > 0) // prune as needs pruning
      {
        // remove the pruned index and its pointer
        // removing the elements increments the cursors i and pointersCurrent, while before stands still
        indices.prune_current();
      }
      else
      {
        // increment all cursors
        indices.next_prune();
      }
    }
    // END (DUST loop)

    // Prune the last index (analoguous with a null (mu* = 0) duality simple test)
    if (lastCost > minCost)
    {
      indices.prune_last();
    }

    // update the available indices
    indices.add(t);
  }

  std::forward_list<unsigned int> changepoints {n};
  for (int newChangepoint = changepointRecord[n]; newChangepoint != 0; newChangepoint = changepointRecord[newChangepoint])
  {
    changepoints.push_front(newChangepoint);
  }

  costRecord.erase(costRecord.begin()); ///// REMOVE FIRST ELEMENT /////
  indices.remove_first(); ///// REMOVE FIRST ELEMENT /////

  return List::create(
    _["changepoints"] = changepoints,
    _["lastIndexSet"] = indices.get_list(),
    _["costQ"] = costRecord
  );
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


RCPP_MODULE(FLATDUST1D)
{
  function("dust.gauss.1D", flat_DUST_1D);
}

RCPP_MODULE(FLAT2DUST1D)
{
  function("dust2.gauss.1D", flat2_DUST_1D);
}


