#include "MD2_A_DUST.h"

using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// --- // Constructor // --- //
DUST_MD2::DUST_MD2(int dual_max, bool random_constraint, Nullable<double> alpha_, Nullable<int> nbLoops)
  : dual_max(dual_max),
    random_constraint(random_constraint),
    indices(nullptr)
{
  if(alpha_.isNull()){alpha = 1e-9;}else{alpha = as<double>(alpha_);}
  if(nbLoops.isNull()){nb_Loops = 10;}else{nb_Loops = as<int>(nbLoops);}
}


////////////////////////////////////////////////////////////////////////////////

DUST_MD2::~DUST_MD2()
{
  delete indices;
}

////////////////////////////////////////////////////////////////////////////////

void DUST_MD2::init_method()
{
  delete indices;

  /// /// ///
  /// /// /// index METHOD
  /// /// ///
  if(random_constraint)
  {
    indices = new RandomIndices_MD2(nb_l, nb_r);
  }
  else
  {
    indices = new VariableIndices_MD2(nb_l, nb_r);
  }

  /// /// ///
  /// /// /// dual_max METHOD
  /// /// ///
  if(dual_max == 0){current_test = &DUST_MD2::dualMaxAlgo0;}
  if(dual_max == 1){current_test = &DUST_MD2::dualMaxAlgo1;}
  if(dual_max == 2){current_test = &DUST_MD2::dualMaxAlgo2;}
  if(dual_max == 3){current_test = &DUST_MD2::dualMaxAlgo3;}
  if(dual_max == 4){current_test = &DUST_MD2::dualMaxAlgo4;}
  if(dual_max == 5){current_test = &DUST_MD2::dualMaxAlgo5;}
  if(dual_max == 6){current_test = &DUST_MD2::dualMaxAlgo6;}
  /// /// /// INIT RANDOM GENERATOR
  /// /// ///
  engine.seed(std::random_device{}());
  dist = std::uniform_real_distribution<double>(0., 1.);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

bool DUST_MD2::dualMaxAlgo0(const double& minCost, const unsigned int& t,
                            const unsigned int& s,
                            std::vector<unsigned int> r,
                            std::vector<unsigned int> r2)
{
  // Draw a random point and evaluate the corresponding dual value
  unsigned int r_size = r.size();

  double constantTerm = (costRecord[s] - minCost) / (t - s); // Dst // !!! CAPTURED IN OPTIM !!! //

  auto col_t = cumsum.col(t);
  auto col_s = cumsum.col(s);
  for (unsigned int row = 0; row < d; row++)
    objectiveMean(row) = (col_t(row) - col_s(row)) / (t - s);

  /// resize the elements:
  linearTerm.resize(r_size);
  constraintMean.resize(d, r_size);
  mu_max.resize(r_size);
  inv_max.resize(r_size);

  double mean_sum = std::accumulate(objectiveMean.begin(), objectiveMean.end(), 0.0)/d;

  // Compute constraintMean matrix and mu_max
  unsigned int j = 0;
  for (auto k: r)
  {
    double constraint_mean_sum = 0;
    linearTerm(j) = (costRecord[s] - costRecord[k]) / (s - k);
    auto col_k = cumsum.col(k);
    for (unsigned int row = 0; row < d; row++)
    {
      constraintMean(row, j) = (col_s(row) - col_k(row)) / (s - k);
      constraint_mean_sum += constraintMean(row, j);
    }

    mu_max(j) = muMax(mean_sum, constraint_mean_sum / d);
    inv_max(j) = pow(mu_max(j), -1);
    j++;
  }

  // Uniform random vector u that sums to 1, scaled depending on mu_max
  std::vector<double> u;
  u.reserve(r_size + 2);
  u.push_back(0.);

  for (unsigned int i = 1; i < r_size + 1; i++)
    u.push_back(dist(engine));

  // mu
  mu.resize(r_size);

  double mu_sum = 0;
  std::sort(u.begin() + 1, u.end());
  for (unsigned int i = 0; i < r_size; i++)
  {
    mu(i) = mu_max(i) * (u[i + 1] - u[i]); /// change if constraint s < r
    mu_sum += mu(i); /// change if constraint s < r
  }
  double inv_sum = pow(mu_sum, -1); ///???

  double linDot = 0;
  for (unsigned int col = 0; col < r_size; col++)
    linDot += mu(col) * linearTerm(col);

  inv_sum = pow(1 - mu_sum, -1);
  double nonLinear = 0;
  for (unsigned int row = 0; row < d; row++)
    nonLinear += Dstar(inv_sum * (objectiveMean(row) - arma::dot(mu, constraintMean.row(row))));

  return constantTerm + linDot - (1 - mu_sum) * nonLinear > 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool DUST_MD2::dualMaxAlgo1(const double& minCost, const unsigned int& t,
                            const unsigned int& s,
                            std::vector<unsigned int> r,
                            std::vector<unsigned int> r2)
{
  // BARYCENTRE test
  unsigned int r_size = r.size();

  double constantTerm = (costRecord[s] - minCost) / (t - s); // Dst // !!! CAPTURED IN OPTIM !!! //

  auto col_t = cumsum.col(t);
  auto col_s = cumsum.col(s);
  for (unsigned int row = 0; row < d; row++)
    objectiveMean(row) = (col_t(row) - col_s(row)) / (t - s);

  linearTerm.resize(r_size);
  constraintMean.resize(d, r_size);

  mu_max.resize(r_size);

  double mean_sum = std::accumulate(objectiveMean.begin(), objectiveMean.end(), 0.0)/d;

  unsigned int j = 0;
  for (auto k: r)
  {
    double constraint_mean_sum = 0;
    linearTerm(j) = (costRecord[s] - costRecord[k]) / (s - k);
    auto col_k = cumsum.col(k);
    for (unsigned int row = 0; row < d; row++)
    {
      constraintMean(row, j) = (col_s(row) - col_k(row)) / (s - k);
      constraint_mean_sum += constraintMean(row, j);
    }

    mu_max(j) = muMax(mean_sum, constraint_mean_sum / d);
    j++;
  }

  mu.resize(r_size);
  double mu_sum = 0;
  double x = pow(r_size + 1, -1);
  for (unsigned int i = 0; i < r_size; i++)
  {
    mu(i) = mu_max(i) * x;
    mu_sum += mu(i);
  }
  double inv_sum = pow(1 - mu_sum, -1);

  double linDot = arma::dot(mu, linearTerm);

  double nonLinear = 0;
  for (unsigned int row = 0; row < d; row++)
    nonLinear += Dstar(inv_sum * (objectiveMean(row) - arma::dot(mu, constraintMean.row(row))));

  return constantTerm + linDot - (1 - mu_sum) * nonLinear > 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool DUST_MD2::dualMaxAlgo2(const double& minCost, const unsigned int& t,
                            const unsigned int& s,
                            std::vector<unsigned int> r,
                            std::vector<unsigned int> r2)
{return(false);
}

bool DUST_MD2::dualMaxAlgo3(const double& minCost, const unsigned int& t,
                            const unsigned int& s,
                            std::vector<unsigned int> r,
                            std::vector<unsigned int> r2)
{return(false);
}

bool DUST_MD2::dualMaxAlgo4(const double& minCost, const unsigned int& t,
                            const unsigned int& s,
                            std::vector<unsigned int> r,
                            std::vector<unsigned int> r2)
{return(false);
}

bool DUST_MD2::dualMaxAlgo5(const double& minCost, const unsigned int& t,
                            const unsigned int& s,
                            std::vector<unsigned int> r,
                            std::vector<unsigned int> r2)
{return(false);
}

bool DUST_MD2::dualMaxAlgo6(const double& minCost, const unsigned int& t,
                            const unsigned int& s,
                            std::vector<unsigned int> r,
                            std::vector<unsigned int> r2)
{return(false);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// --- // Fits the data, i. e. initializes all data-dependent vectors // --- //
void DUST_MD2::init(const arma::dmat& inData,
                    Nullable<double> inPenalty,
                    Nullable<unsigned int> inNbL,
                    Nullable<unsigned int> inNbR)
{
  n = inData.n_cols;
  d = inData.n_rows;

  if (inPenalty.isNull())
  {
    penalty = 2 * d * std::log(n);
  }
  else
  {
    penalty = as<double>(inPenalty);
  }

  /// read the number of constraints + default choice
  if (inNbR.isNull()){nb_r = 1;}else{nb_r = std::min(d, as<unsigned int>(inNbR));}
  if (inNbL.isNull()){nb_l = d - 1;}else{nb_l = std::min(d, as<unsigned int>(inNbL));}
  nb_max = nb_l + nb_r;


  changepointRecord = std::vector<int>(n + 1, 0);
  nb_indices = std::vector<int>(n, 0);

  cumsum = arma::dmat(d, n + 1);
  costRecord = std::vector<double>(n + 1, -penalty);

  init_method();

  indices->set_size(n);
  indices->add(0);
  indices->add(1);

  mu             = arma::rowvec(nb_max);
  mu_max         = arma::rowvec(nb_max);
  inv_max        = arma::rowvec(nb_max);
  grad           = arma::rowvec(nb_max);

  linearTerm     = arma::rowvec(nb_max);

  m_value        = arma::colvec(d);
  objectiveMean  = arma::colvec(d);
  constraintMean = arma::dmat(d, nb_max);
  nonLinearGrad  = arma::colvec(d);

  Identity       = arma::dmat(nb_max, nb_max, arma::fill::eye);
  inverseHessian = arma::dmat(nb_max, nb_max);
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// --- // Algorithm-specific method // --- //
void DUST_MD2::compute(const arma::dmat& inData)
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

  for (unsigned int row = 0; row < d; row++)
    cumsum(row, 1) = statistic(inData(row));

  costRecord[1] = Cost(t, s);
  changepointRecord[1] = 0;

  int nbt = 2;
  nb_indices[0] = 1;

  // Main loop
  for (t = 2; t <= n; t++)
  {
    // update cumsum
    auto col_prev = cumsum.col(t - 1);
    for (unsigned int row = 0; row < d; row++)
      cumsum(row, t) = col_prev(row) + statistic(inData(row, t - 1));

    // OP step
    indices->reset();
    minCost = std::numeric_limits<double>::infinity();
    do
    {
      s = *(indices->current);
      // Rcout << "t = " << t << "; s = " << s << std::endl;
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
      if ((this->*current_test)(minCost, t,
                                *(indices->current),
                                  indices->get_constraints_l(),
                                  indices->get_constraints_r())) // prune as needs pruning
      {
        // remove the pruned index and its pointer
        // removing the elements increments the cursors i and pointerIt, while before stands still
        indices->prune_current();
        nbt--;
      }
      else
      {
        // increment all cursors
        indices->next_prune();
      }
    }
    // END (DUST loop)

    // Prune the last index (analoguous with a null (mu* = 0) duality simple test)
    // if (lastCost > minCost)
    // {
    //   indices->prune_last();
    //   nbt--;
    // }

    // update the available indices
    indices->add(t);
    nb_indices[t - 1] = nbt;
    nbt++;
  }

}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


// --- // Builds changepoints // --- //
std::forward_list<unsigned int> DUST_MD2::backtrack_changepoints()
{
  std::forward_list<unsigned int> changepoints {n};
  for (int newChangepoint = changepointRecord[n]; newChangepoint != 0; newChangepoint = changepointRecord[newChangepoint])
  {
    changepoints.push_front(newChangepoint);
  }
  return changepoints;
}


// --- // Retrieves optimal partition // --- //
List DUST_MD2::get_partition()
{
  costRecord.erase(costRecord.begin()); ///// REMOVE FIRST ELEMENT /////
  indices->remove_last(); ///// REMOVE FIRST ELEMENT /////

  std::forward_list<unsigned int> chpts = backtrack_changepoints();
  std::vector<unsigned int> lastIndexSet = indices->get_list();
  std::reverse(lastIndexSet.begin(), lastIndexSet.end());

  return List::create(
    _["changepoints"] = chpts,
    _["lastIndexSet"] = lastIndexSet,
    _["nb"] = nb_indices,
    _["costQ"] = costRecord
  );
}


// --- // Wrapper method for quickly computing               // --- //
// --- // and retrieving the optimal partition of input data // --- //
List DUST_MD2::quick(const arma::dmat& inData,
                     Nullable<double> inPenalty,
                     Nullable<unsigned int> inNbL,
                     Nullable<unsigned int> inNbR)
{
  init(inData, inPenalty, inNbL, inNbR);
  compute(inData);
  return get_partition();
}

////
//// IDEA : propose a new quick method with a loop of "compute (K)"
//// solving the K fixed (number of change) problem.
//// new 3 functions : init, compute and get_partition
////






