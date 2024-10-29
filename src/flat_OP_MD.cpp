#include "flat_Gauss_MD.h"

#include <forward_list>

using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// --- // Fits the data, i. e. initializes all data-dependent vectors // --- //
//[[Rcpp::export]]
List flat_OP_MD(const arma::dmat& inData, Nullable<double> inPenalty = R_NilValue)
{
  unsigned int n = inData.n_cols;
  unsigned int d = inData.n_rows;

  double penalty;
  if (inPenalty.isNull())
  {
    penalty = 2 * d * std::log(n); //to do
  }
  else
  {
    penalty = as<double>(inPenalty);
  }

  std::vector<int> changepointRecord;
  changepointRecord.reserve(n + 1);
  changepointRecord.push_back(0);

  std::vector<double> costRecord;
  costRecord.reserve(n + 1);
  costRecord.push_back(-penalty);

  arma::dmat cumsum(d, n + 1);

  // Initialize OP step value
  double lastCost;
  // then keeps the cost of the model with last changepoint at the first possible index in the t-th OP step ...
  // ... storing it allows pruning of the first available index
  double minCost;
  unsigned int argMin; // stores the optimal last changepoint for the current OP step

  // First OP step (t = 1)
  unsigned int t = 1;
  unsigned int s = 0;

  for (unsigned int row = 0; row < d; row++)
    cumsum(row, 1) = statistic_MD(inData(row));

  costRecord.push_back(Cost_MD(t, s, cumsum));
  changepointRecord.push_back(0);

  // Main loop
  for (t = 2; t <= n; t++)
  {
    // update cumsum
    auto col_prev = cumsum.col(t - 1);
    for (unsigned int row = 0; row < d; row++)
      cumsum(row, t) = col_prev(row) + statistic_MD(inData(row, t - 1));

    // OP step
    minCost = std::numeric_limits<double>::infinity();
    for (unsigned int s = 0; s < t; s++)
    {
      // Rcout << "t = " << t << "; s = " << s << std::endl;
      lastCost = costRecord[s] + Cost_MD(t, s, cumsum);
      if (lastCost < minCost)
      {
        minCost = lastCost;
        argMin = s;
      }
    }
    // END (OP step)

    // OP update
    costRecord.push_back(minCost + penalty);
    changepointRecord.push_back(argMin);
  }

  std::forward_list<unsigned int> changepoints {n};
  for (int newChangepoint = changepointRecord[n]; newChangepoint != 0; newChangepoint = changepointRecord[newChangepoint])
  {
    changepoints.push_front(newChangepoint);
  }

  changepoints.push_front(0);

  double totalCost = 0;
  auto it = changepoints.begin();      // Iterator to the current element
  auto next_it = std::next(it); // Iterator to the next element

  while (next_it != changepoints.end())
  {
    totalCost += Cost_MD(*next_it, *it, cumsum);  // Apply cost function to consecutive elements
    it = next_it;  // Move to the next element
    ++next_it;     // Advance the next iterator
  }

  changepoints.pop_front();

  costRecord.erase(costRecord.begin()); ///// REMOVE FIRST ELEMENT /////

  return List::create(
    _["changepoints"] = changepoints,
    _["segmentation_cost"] = totalCost,
    _["costQ"] = costRecord
  );
}

////
//// IDEA : propose a new quick method with a loop of "compute (K)"
//// solving the K fixed (number of change) problem.
//// new 3 functions : ini, comute and get_partition
////






