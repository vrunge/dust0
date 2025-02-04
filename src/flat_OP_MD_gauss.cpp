#include "flat_Gauss_MD.h"

#include <forward_list>

using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// --- // Fits the data, i. e. initializes all data-dependent vectors // --- //

//' Optimal Partitioning in MD with Flat Model
//'
//' Computes the optimal partitioning of multi-dimensional data using a flat model with an optional penalty parameter.
//'
//' @param inData A numeric matrix representing the input data for partitioning.
//' @param inPenalty An optional numeric penalty parameter to control the number of partitions. Defaults to \code{NULL}, indicating a default penalty is used.
//' @return A list containing the results of the optimal partitioning, including identified change points and other model details.
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
  unsigned int argMin = 0; // stores the optimal last changepoint for the current OP step

  // First OP step (t = 1)
  unsigned int t = 1;
  unsigned int s = 0;

  for (unsigned int row = 0; row < d; row++)
    cumsum(row, 1) = inData(row);

  costRecord.push_back(CostGauss_MD(t, s, cumsum));
  changepointRecord.push_back(0);

  // Main loop
  for (t = 2; t <= n; t++)
  {
    // update cumsum
    auto col_prev = cumsum.col(t - 1);
    for (unsigned int row = 0; row < d; row++)
      cumsum(row, t) = col_prev(row) + inData(row, t - 1);

    // OP step
    minCost = std::numeric_limits<double>::infinity();
    for (unsigned int s = 0; s < t; s++)
    {
      lastCost = costRecord[s] + CostGauss_MD(t, s, cumsum);
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
  changepoints.pop_front();
  costRecord.erase(costRecord.begin()); ///// REMOVE FIRST ELEMENT /////

  return List::create(
    _["changepoints"] = changepoints,
    _["costQ"] = costRecord
  );
}

////
//// IDEA : propose a new quick method with a loop of "compute (K)"
//// solving the K fixed (number of change) problem.
//// new 3 functions : ini, comute and get_partition
////

RCPP_MODULE(FLATOPMD)
{
  function("OP.gauss.MD", flat_OP_MD);
}






