#include "MD_A_DUST.h"

using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// --- // Constructor // --- //
DUST_MD::DUST_MD(int dual_max, bool random_constraint, Nullable<double> alpha_, Nullable<int> nbLoops)
  : dual_max(dual_max),
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

DUST_MD::~DUST_MD()
{
  delete indices;
}

void DUST_MD::init_method()
{
  delete indices;

  /// /// ///
  /// /// /// index METHOD
  /// /// ///
  if(random_constraint)
  {
    indices = new RandomIndices_MD(nb_max);
  }
  else
  {
    indices = new VariableIndices_MD(nb_max);
  }

  /// /// ///
  /// /// /// dual_max METHOD
  /// /// ///
  if(dual_max == 0){current_test = &DUST_MD::dualMaxAlgo0;}
  if(dual_max == 1){current_test = &DUST_MD::dualMaxAlgo1;}
  if(dual_max == 2){current_test = &DUST_MD::dualMaxAlgo2;}
  if(dual_max == 3){current_test = &DUST_MD::dualMaxAlgo3;}
  if(dual_max == 4){current_test = &DUST_MD::dualMaxAlgo4;}
  if(dual_max == 5){current_test = &DUST_MD::dualMaxAlgo5;}
  if(dual_max == 6){current_test = &DUST_MD::dualMaxAlgo6;}
  /// /// /// INIT RANDOM GENERATOR
  /// /// ///
  engine.seed(std::random_device{}());
  dist = std::uniform_real_distribution<double>(0., 1.);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

bool DUST_MD::dualMaxAlgo0(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> r)
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


bool DUST_MD::dualMaxAlgo1(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> r)
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


bool DUST_MD::dualMaxAlgo2(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> r)
{return(false);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool DUST_MD::dualMaxAlgo3(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> r)
{return(false);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool DUST_MD::dualMaxAlgo4(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> r)
{

  // ############################################## //
  // ############################################## //
  // ######### // FIRST STEP AT (0, 0) // ######### //
  // ############################################## //
  // ############################################## //


  // ######### // PELT TEST // ######### //
  // Formula: Dst - D*(Sst)              //

  double constantTerm = (costRecord[s] - minCost) / (t - s); // Dst // !!! CAPTURED IN OPTIM !!! //

  auto col_t = cumsum.col(t);
  auto col_s = cumsum.col(s);

  double nonLinear = 0; // D*(Sst) // !!! UPDATED IN OPTIM !!! //
  for (unsigned int row = 0; row < d; row++)
  {
    objectiveMean(row) = (col_t(row) - col_s(row)) / (t - s);
    nonLinear += Dstar(objectiveMean(row));
  }

  double test_value = constantTerm - nonLinear; // !!! UPDATED IN OPTIM !!! //

  if (test_value > 0) { return true; } // PELT test
  //
  //
  // ######### // TANGENT HYPERPLANE TEST // ######### //
  // Formula: D(mu) + (mu_tan - mu) * grad(mu)         //
  // mu_tan is the highest point on the tangent hyper- //
  // plane at mu. values 0 except at                   //
  // i* = argmax(grad(mu)), tan_mu[i*] = mu_max[i*]    //
  // if formula <= 0, then pruning is impossible       //


  // First compute the constraint-related objects
  unsigned int r_size = r.size(); // !!! UPDATED IN OPTIM !!! //

  linearTerm.resize(r_size);
  constraintMean.resize(d, r_size);

  mu_max.resize(r_size);
  inv_max.resize(r_size);

  double mean_sum = std::accumulate(objectiveMean.begin(), objectiveMean.end(), 0.0)/d;


  // Initialize the constraint mean matrix, which contains the mean-vectors associated with each constraint as columns
  // + Initialize mu_max, which is based on the value of the value of each mean-vector compared to the value of the "objective" mean-vector
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

  // Coordinates of the point that maximizes the hyperplane tangent defined at the current mu (init)
  arma::rowvec tangent_max(r_size); // !!! SHRUNK IN OPTIM !!! //
  double grad_max = -std::numeric_limits<double>::infinity(); // !!! UPDATED IN OPTIM !!! //
  unsigned int grad_argmax = 0; // !!! UPDATED IN OPTIM !!! //

  for (unsigned int row = 0; row < d; row++)
    nonLinearGrad(row) = DstarPrime(objectiveMean(row));

  grad.resize(r_size);

  // Grad and hyperplane tangent compuation (formula: linearTerm + nonLinear * t(1_p) + t(nonLinearGrad) * (Srs - Sst * t(1_p)))
  for (unsigned int col = 0; col < r_size; col++)
  {
    double dot_product = 0;
    auto col_k = constraintMean.col(col);
    for (unsigned int row = 0; row < d; row++)
    {
      dot_product += nonLinearGrad(row) * (objectiveMean(row) - col_k(row));
    }
    grad(col) = linearTerm(col) + nonLinear - dot_product;

    double norm = grad(col) * inv_max(col);
    if (norm > grad_max)
    {
      grad_max = norm;
      grad_argmax = col;
    }
  }

  if (grad_max > 0)
  {
    tangent_max(grad_argmax) = mu_max(grad_argmax); // maximum value on the hyperplane is at the corner of the simplex corresponding to the largest grad value
    if (test_value + arma::dot(tangent_max, grad) <= 0) { return false; } // check if the tangent hyperplane ever reaches 0 on the simplex triangle (from concavity property of the dual)
  }
  else return false; // if grad is fully negative, then no improvement can be made on the test value


  // ############################################## //
  // ############################################## //
  // ######### // OPTIM INITIALIZATION // ######### //
  // ############################################## //
  // ############################################## //
  // Initialize all dynamic objects used in the op- //
  // timization recursion.                          //
  // ############################################## //

  mu.resize(r_size);
  for (unsigned int col = 0; col < r_size; col++)
    mu(col) = 0;

  double mu_sum = 0; // !!! UPDATED IN OPTIM !!! //
  double inv_sum = 0; // !!! UPDATED IN OPTIM !!! //

  inverseHessian.resize(r_size, r_size);
  arma::dmat I = Identity.submat(0, 0, r_size - 1, r_size - 1);

  // Initialize inverseHessian as minus identity
  for (unsigned int row = 0; row < r_size; row++)
    for(unsigned int col = 0; col < r_size; col++)
    {
      if (row == col) inverseHessian(row, col) = -1;
      else inverseHessian(row, col) = 0;
    }

  // // ######################################### //
  // // ######################################### //
  // // ######### // OPTIM RECURSION // ######### //
  // // ######################################### //
  // // ######################################### //

  std::function<bool(std::vector<unsigned int>)> optim = [&] (std::vector<unsigned int> zero_index)
  {
    if (zero_index.size() > 0)
    {
      // Shrink all relevant objects
      r_size -= zero_index.size();
      if (r_size <= 0) { return false; }

      grad_argmax = 0;
      for (auto k = zero_index.rbegin(); k != zero_index.rend(); ++k)
      {
        r.erase(r.begin() + *k);

        linearTerm.shed_col(*k);
        constraintMean.shed_col(*k);

        tangent_max.shed_col(*k);

        mu.shed_col(*k);
        grad.shed_col(*k);

        mu_max.shed_col(*k);
        inv_max.shed_col(*k);

        inverseHessian.shed_col(*k);
        inverseHessian.shed_row(*k);

      }
      I = Identity.submat(0, 0, r_size - 1, r_size - 1);
    }

    bool shrink = false;
    std::vector<unsigned int> shrink_indices;
    shrink_indices.reserve(r_size);

    arma::rowvec direction(r_size); // direction and intensity of the update
    double direction_scale; // scaling applied to direction when direction pushes past boundaries
    arma::rowvec mu_diff(r_size); // (mu+) - mu
    arma::rowvec grad_diff = -grad; // (g+) - g; initialized as -g to avoid storing 2 values of gradient.

    auto updateTestValue = [&] ()
    {
      // Armijo step-size:
      // evaluate D at mu + dk
      // test based on gradient value and m1
      // if false, scale dk by half
      // repeat until test is true

      // 1. Update mu and D(mu).
      double linDot = 0; // mu.dot(linearTerm);
      for (unsigned int col = 0; col < r_size; col++)
      {
        mu_diff(col) = direction(col) * direction_scale;

        mu(col) += mu_diff(col);
        mu_sum += mu_diff(col);

        linDot += mu(col) * linearTerm(col);
      }

      inv_sum = pow(1 - mu_sum, -1);
      nonLinear = 0;
      for (unsigned int row = 0; row < d; row++)
      {
        m_value(row) = inv_sum * (objectiveMean(row) - arma::dot(mu, constraintMean.row(row)));
        nonLinear += Dstar(m_value(row));
      }

      double new_test = constantTerm + linDot - (1 - mu_sum) * nonLinear;

      // 2. define gradient condition
      arma::rowvec gradCondition(r_size);
      for (unsigned int col = 0; col < r_size; col++)
        gradCondition(col) = m1 * grad(col);

      // 3. Scale dk until test is valid
      unsigned int iter = 0;
      while(new_test < test_value + arma::dot(mu_diff, gradCondition) && iter < 100)
      {
        for(unsigned int col = 0; col < r_size; col++)
        {
          mu_diff(col) *= .5;
          mu(col) -= mu_diff(col);
          linDot -= mu_diff(col) * linearTerm(col);
          mu_sum -= mu_diff(col);
        }

        inv_sum = pow(1 - mu_sum, -1);
        nonLinear = 0;
        for (unsigned int row = 0; row < d; row++)
        {
          m_value(row) = inv_sum * (objectiveMean(row) - arma::dot(mu, constraintMean.row(row)));
          nonLinear += Dstar(m_value(row));
        }

        new_test = constantTerm + linDot - (1 - mu_sum) * nonLinear; // update values

        iter++;
      }
      test_value = new_test;

      // Project mu onto the interior simplex
      // If any null value, then shrink
      for (unsigned int col = 0; col < r_size; col++)
      {
        if(mu(col) < 0) { mu_sum -= mu(col); mu(col) = 0; shrink = true; shrink_indices.push_back(col); }
        else if(mu(col) == 0) { shrink = true; shrink_indices.push_back(col); }
      }
    };

    auto updateGrad = [&] ()
    {
      // Initialize tangent max search.
      grad_max = -std::numeric_limits<double>::infinity();
      grad_argmax = 0;

      for (unsigned int row = 0; row < d; row++)
        nonLinearGrad(row) = DstarPrime(m_value(row));

      // Update grad value.
      for (unsigned int col = 0; col < r_size; col++)
      {
        double dot_product = 0;
        auto col_k = constraintMean.col(col);
        for (unsigned int row = 0; row < d; row++)
          dot_product += nonLinearGrad(row) * (col_k(row) - m_value(row));
        grad(col) = linearTerm(col) + nonLinear + dot_product;
        grad_diff(col) += grad(col);
        if (grad_diff(col) == 0)
        {
          if (grad(col) > 0) { grad_diff(col) = -1e-16; }
          else { grad_diff(col) = 1e-16; }
        }

        double norm = grad(col) * inv_max(col);
        if (norm > grad_max)
        {
          grad_max = norm;
          grad_argmax = col;
        }
      }

      // maximum value on the hyperplane is at the corner of the simplex corresponding to the largest grad value. if no value positive, then it is at (0,...)
      if (grad_max > 0)
        tangent_max(grad_argmax) = mu_max(grad_argmax);
    };

    // Optim loop
    unsigned int iter = 0;
    do
    {
      // Update direction and intensity then clip it to stay within the bounds of the positive simplex
      direction = (-grad) * inverseHessian;
      direction_scale = FindBoundaryCoef(mu, direction, inv_max, shrink, shrink_indices); // may trigger shrink if pushing past boundary
      if (shrink)
      {
        return optim(shrink_indices);
      }
      if (direction_scale == 0) { return false; } // stops optimization if no movement is produced

      // update mu and D(mu) + check for shrink
      updateTestValue();

      if(test_value > 0) { return true; } // success, index s is pruned

      // update grad and tangent max location
      updateGrad();

      // compute tangent max value
      double dot_product = 0;
      for (unsigned int col = 0; col < r_size; col++)
        dot_product += grad(col) * (tangent_max(col) - mu(col));

      if (test_value + dot_product <= 0) { return false; } // check if the tangent hyperplane ever reaches 0 on the simplex triangle (from concavity property of the dual)

      tangent_max(grad_argmax) = 0;

      // trigger shrinking, which reduces the dimension of the solution space by one
      if (shrink) return optim(shrink_indices);

      // Update inverse hessian estimation
      updateHessian(inverseHessian, mu_diff, grad_diff, I);
      grad_diff = -grad;
      iter++;
    } while (iter < 50);

    // Rcout << "Reached iter limit" << std::endl;

    return false;
  };

  // return false;
  return optim(std::vector<unsigned int>());
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

bool DUST_MD::dualMaxAlgo5(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> r)
{
  // ############################################## //
  // ############################################## //
  // ######### // FIRST STEP AT (0, 0) // ######### //
  // ############################################## //
  // ############################################## //


  // ######### // PELT TEST // ######### //
  // Formula: Dst - D*(Sst)              //

  double constantTerm = (costRecord[s] - minCost) / (t - s); // Dst // !!! CAPTURED IN OPTIM !!! //

  auto col_t = cumsum.col(t);
  auto col_s = cumsum.col(s);

  double nonLinear = 0; // D*(Sst) // !!! UPDATED IN OPTIM !!! //
  for (unsigned int row = 0; row < d; row++)
  {
    objectiveMean(row) = (col_t(row) - col_s(row)) / (t - s);
    nonLinear += Dstar(objectiveMean(row));
  }

  double test_value = constantTerm - nonLinear; // !!! UPDATED IN OPTIM !!! //

  return test_value > 0; // PELT test
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool DUST_MD::dualMaxAlgo6(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> r)
{return(false);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// --- // Fits the data, i. e. initializes all data-dependent vectors // --- //
void DUST_MD::init(const arma::dmat& inData, Nullable<double> inPenalty, Nullable<unsigned int> inNbMax)
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
  /// if not precised, the number of constraints =
  ///  = the dimension of the parameter space  = number of time series
  if (inNbMax.isNull())
  {
    nb_max = d;
  }
  else
  {
    nb_max = std::min(d, as<unsigned int>(inNbMax));
  }

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
void DUST_MD::compute(const arma::dmat& inData)
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
      if ((this->*current_test)(minCost, t, *(indices->current), indices->get_constraints())) // prune as needs pruning
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
std::forward_list<unsigned int> DUST_MD::backtrack_changepoints()
{
  std::forward_list<unsigned int> changepoints {n};
  for (int newChangepoint = changepointRecord[n]; newChangepoint != 0; newChangepoint = changepointRecord[newChangepoint])
  {
    changepoints.push_front(newChangepoint);
  }
  return changepoints;
}


// --- // Retrieves optimal partition // --- //
List DUST_MD::get_partition()
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
List DUST_MD::quick(const arma::dmat& inData, Nullable<double> inPenalty, Nullable<unsigned int> inNbMax)
{
  init(inData, inPenalty, inNbMax);
  compute(inData);
  return get_partition();
}

////
//// IDEA : propose a new quick method with a loop of "compute (K)"
//// solving the K fixed (number of change) problem.
//// new 3 functions : ini, comute and get_partition
////






