#include <Rcpp.h>
#include <forward_list>

using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// --- // Fits the data, i. e. initializes all data-dependent vectors // --- //

//' Optimal Partitioning in 1D with Flat Model
//'
//' Computes the optimal partitioning of one-dimensional data using a flat model with an optional penalty parameter.
//'
//' @param inData A numeric vector representing the input data for partitioning.
//' @param inPenalty An optional numeric penalty parameter to control the number of partitions. Defaults to \code{NULL}, indicating a default penalty is used.
//' @return A list containing the results of the optimal partitioning, including identified change points and other model details.
//' @examples
//' data <- rnorm(100)
//' result <- flat_OP_1D(data, inPenalty = 1.0)
//[[Rcpp::export]]
List flat_OP_1D(const std::vector<double>& inData, Nullable<double> inPenalty = R_NilValue)
{
  unsigned int n = inData.size();

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
  costRecord.push_back(- .5 * pow(cumsum[t] - cumsum[s], 2) / (t - s));
  changepointRecord.push_back(0);

  // Main loop
  for (t = 2; t <= n; t++)
  {
    // update cumsum
    cumsum.push_back(cumsum.back() + inData[t - 1]);

    // OP step
    minCost = std::numeric_limits<double>::infinity();
    for (unsigned int s = 0; s < t; s++)
    {
      lastCost = costRecord[s] - .5 * pow(cumsum[t] - cumsum[s], 2) / (t - s);
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


RCPP_MODULE(FLATOP1D)
{
  function("OP.gauss.1D", flat_OP_1D);
}




