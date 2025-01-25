#include "MD_A_DUST.h"
#include "fstream"

using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// --- // Constructor // --- //
DUST_MD::DUST_MD(int dual_max, bool random_constraint, Nullable<int> nbLoops)
  : n(0),
    d(0),
    dual_max(dual_max),
    random_constraint(random_constraint),
    indices(nullptr)
{
  if(nbLoops.isNull()){nb_Loops = 10;}else{nb_Loops = as<int>(nbLoops);}
}


////////////////////////////////////////////////////////////////////////////////

DUST_MD::~DUST_MD()
{
  delete indices;
}

////////////////////////////////////////////////////////////////////////////////

void DUST_MD::pruning_method()
{
  delete indices;

  /// /// ///
  /// /// /// index METHOD
  /// /// ///
  if(random_constraint){indices = new RandomIndices_MD(nb_l, nb_r);}
  else{indices = new DeterministicIndices_MD(nb_l, nb_r);}

  /// /// ///
  /// /// /// dual_max METHOD
  /// /// ///
  if(dual_max == 0){current_test = &DUST_MD::dualMaxAlgo0;}
  if(dual_max == 1){current_test = &DUST_MD::dualMaxAlgo1;}
  if(dual_max == 2){current_test = &DUST_MD::dualMaxAlgo2;}
  if(dual_max == 3){current_test = &DUST_MD::dualMaxAlgo3;}
  if(dual_max == 4){current_test = nb_r > 0 ? &DUST_MD::dualMaxAlgo42 : &DUST_MD::dualMaxAlgo4;}
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


// --- // Fits the data, i. e. initializes all data-dependent vectors // --- //
void DUST_MD::append(const arma::dmat& inData,
                     Nullable<double> inPenalty,
                     Nullable<unsigned int> inNbL,
                     Nullable<unsigned int> inNbR)
{
  bool first_execution = (n == 0);

  n += inData.n_cols;

  if (!first_execution)
  {
    if (d != inData.n_rows)
      throw std::invalid_argument("new data has invalid n_rows. got " + std::to_string(inData.n_rows) + ", expected " + std::to_string(d));
  }
  else { d = inData.n_rows; }

  if (inPenalty.isNull()){penalty = 2 * d * std::log(n);}else{penalty = as<double>(inPenalty);}

  changepointRecord.reserve(n + 1);
  nb_indices.reserve(n);
  costRecord.reserve(n + 1);

  cumsum.resize(d, n + 1);

  if (first_execution)
  {
    changepointRecord.push_back(0);
    nb_indices.push_back(1);
    costRecord.push_back(-penalty);
    for (unsigned row = 0; row < d; row++)
      cumsum(row, 0) = 0.;

    // read the number of constraints + default choice
    if (inNbL.isNull()){nb_l = d - 1;}else{nb_l = std::min(d, as<unsigned int>(inNbL));}
    if (inNbR.isNull()){nb_r = 1;}else{nb_r = std::min(d - nb_l, as<unsigned int>(inNbR));}
    nb_max = nb_l + nb_r;

    pruning_method();
    indices->set_init_size(n);
    indices->add(0);

    // Initialize optim objects
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
  else{ indices->set_init_size(n); }

  // store the data as cumsum
  unsigned current_filled_cols = cumsum.n_cols - inData.n_cols - 1; // -1 for correct indexing first column
  for (unsigned int data_col = 0; data_col < inData.n_cols; data_col++)
  {
    unsigned cumsum_col = data_col + current_filled_cols;
    for (unsigned int row = 0; row < d; row++)
    {
      cumsum(row, cumsum_col + 1) = cumsum(row, cumsum_col) + statistic(inData(row, data_col));
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// --- // Algorithm-specific method // --- //
void DUST_MD::update_partition()
{
  std::ofstream log_file;
  log_file.close();
  log_file.open("output.log", std::ofstream::out | std::ofstream::app);
  log_file.clear();
  auto cout_buf = Rcout.rdbuf(log_file.rdbuf());
  Rcout << "correctly setup log_file" << std::endl;

  int nbt = nb_indices.back();

  // Main loop
  // Main loop
  for (unsigned t = indices->get_first() + 1; t <= n; t++)
  {
    // OP step
    // OP step
    indices->reset();
    double lastCost;
    double minCost = std::numeric_limits<double>::infinity();
    unsigned argMin;
    do
    {
      unsigned s = *indices->current;
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
    // END (OP step)

    // OP update
    minCost += penalty;
    costRecord.push_back(minCost);
    changepointRecord.push_back(argMin);

    // DUST step
    // DUST step
    indices->reset_prune();

    // DUST loop
    while (indices->check())
    {
      auto l = indices->get_constraints_l();
      auto r = indices->get_constraints_r();
      Rcout << l.size() << ", " << r.size() << std::endl;
      if ((this->*current_test)(minCost, t,
           *(indices->current),
           indices->get_constraints_l(),
           indices->get_constraints_r())) // prune as needs pruning
      {
        // remove the pruned index
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
    // END (DUST loop)

    // Prune the last index (analogous with a null (mu* = 0) duality simple test)
    // if (lastCost > minCost)
    // {
    //   indices->prune_last();
    //   nbt--;
    // }

    indices->add(t);
    nb_indices.push_back(nbt);
    nbt++;
  }

  Rcout.rdbuf(cout_buf);
  log_file.close();
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

List DUST_MD::get_info()
{
  return List::create(
    _["data_statistic"] = cumsum,
    _["data_dimensions"] = std::vector<unsigned> { d, n },
    _["current_penalty"] = penalty,
    _["model"] = get_model(),
    _["pruning_algo"] = dual_max,
    _["pruning_random_constraint"] = random_constraint,
    _["pruning_nb_constraints"] = std::vector<unsigned> { nb_l, nb_r },
    _["pruning_nb_loops"] = nb_Loops
  );
}

// --- // Retrieves optimal partition // --- //
List DUST_MD::get_partition()
{
  // costRecord.erase(costRecord.begin()); ///// REMOVE FIRST ELEMENT /////
  // indices->remove_last(); ///// REMOVE FIRST ELEMENT /////

  std::forward_list<unsigned int> chpts = backtrack_changepoints();
  std::vector<unsigned int> lastIndexSet = indices->get_list();
  std::reverse(lastIndexSet.begin(), lastIndexSet.end());

  return List::create(
    _["changepoints"] = chpts,
    _["lastIndexSet"] = lastIndexSet,
    _["nb"] = std::vector<unsigned>(nb_indices.begin() + 1, nb_indices.end()),
    _["costQ"] = std::vector<double>(costRecord.begin() + 1, costRecord.end())
  );
}


// --- // Wrapper method for quickly computing               // --- //
// --- // and retrieving the optimal partition of input data // --- //
List DUST_MD::one_dust(const arma::dmat& inData,
                       Nullable<double> inPenalty,
                       Nullable<unsigned int> inNbL,
                       Nullable<unsigned int> inNbR)
{
  append(inData, inPenalty, inNbL, inNbR);
  update_partition();
  return get_partition();
}

// arma::colvec objectiveMean;
// arma::dmat constraintMean;
// arma::rowvec linearTerm;
// double constantTerm;
// arma::rowvec mu;
double dual_Eval()
{

  return(-std::numeric_limits<double>::infinity());
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool DUST_MD::dualMaxAlgo0(const double& minCost, const unsigned int& t,
                            const unsigned int& s,
                            std::vector<unsigned int> l,
                            std::vector<unsigned int> r)
{
  // Draw a random point and evaluate the corresponding dual value
  unsigned int r_size = l.size();

  constantTerm =  - (minCost - costRecord[s]) / (t - s); // Dst // !!! CAPTURED IN OPTIM !!! //

  for (unsigned int row = 0; row < d; row++)
    objectiveMean(row) = (cumsum(t,row) - cumsum(s,row)) / (t - s);

  /// resize the elements:
  linearTerm.resize(r_size);
  constraintMean.resize(d, r_size);
  mu_max.resize(r_size);
  inv_max.resize(r_size);

  double mean_sum = std::accumulate(objectiveMean.begin(), objectiveMean.end(), 0.0)/d;

  // Compute constraintMean matrix and mu_max
  unsigned int j = 0;
  for (auto k: l)
  {
    double constraint_mean_sum = 0;
    linearTerm(j) = (costRecord[s] - costRecord[k]) / (s - k);

    for (unsigned int row = 0; row < d; row++)
    {
      constraintMean(row, j) = (cumsum(s,row) - cumsum(k,row)) / (s - k);
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
    mu(i) = mu_max(i) * (u[i + 1] - u[i]); /// change if constraint s < l
    mu_sum += mu(i); /// change if constraint s < l
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

// BARYCENTRE test
// BARYCENTRE test
// BARYCENTRE test

bool DUST_MD::dualMaxAlgo1(const double& minCost, const unsigned int& t,
                            const unsigned int& s,
                            std::vector<unsigned int> l,
                            std::vector<unsigned int> r)
{
  unsigned int r_size = l.size();
  unsigned int r2_size = r.size();

  ///////
  /////// constantTerm AND objectiveMean
  ///////
  constantTerm = - (minCost - costRecord[s]) / (t - s); // Dst // !!! CAPTURED IN OPTIM !!! //

  auto col_t = cumsum.col(t);
  auto col_s = cumsum.col(s);
  for (unsigned int row = 0; row < d; row++)
    objectiveMean(row) = (col_t(row) - col_s(row)) / (t - s);

  ///////
  /////// resize vectors. r_size + r2_size constraints
  ///////
  linearTerm.resize(r_size + r2_size);
  constraintMean.resize(d, r_size + r2_size);
  mu_max.resize(r_size + r2_size);

  /// /d ????
  double mean_sum = std::accumulate(objectiveMean.begin(), objectiveMean.end(), 0.0)/d;

  ///////
  /////// constraintMean, constraint_mean_sum
  ///////
  unsigned int j = 0;
  for (auto k: l) ///////// WITH l
  {
    double constraint_mean_sum = 0;
    linearTerm(j) = (costRecord[s] - costRecord[k]) / (s - k);
    auto col_k = cumsum.col(k);
    for (unsigned int row = 0; row < d; row++)
    {
      constraintMean(row, j) = (col_s(row) - col_k(row)) / (s - k);
      constraint_mean_sum += constraintMean(row, j);
    }

    mu_max(j) = muMax(mean_sum, constraint_mean_sum / d); /// WHY? //// to be reviewed
    j++;
  }
  for (auto k: r) ///////// SAME WITH r
  {
    double constraint_mean_sum = 0;
    linearTerm(j) = (costRecord[s] - costRecord[k]) / (s - k);
    auto col_k = cumsum.col(k);
    for (unsigned int row = 0; row < d; row++)
    {
      constraintMean(row, j) = (col_s(row) - col_k(row)) / (s - k);
      constraint_mean_sum += constraintMean(row, j);
    }

    mu_max(j) = muMax(mean_sum, constraint_mean_sum / d); /// WHY? //// to be reviewed... (no mu_max)
    j++;
  }


  ///////
  /////// mu
  ///////
  mu.resize(r_size + r2_size);
  double mu_sum = 0;
  double x = pow(r_size + r2_size + 1, -1); ///

  for (unsigned int i = 0; i < r_size; i++) ///////// WITH l
  {
    mu(i) = mu_max(i) * x;
    mu_sum += mu(i);
  }
  for (unsigned int i = 0; i < r2_size; i++) ///////// WITH r
  {
    mu(r_size + i) = - mu_max(r_size + i) * x; /// with a minus here
    mu_sum = mu(r_size + i);
  }

  ///////
  /////// dual value in mu
  ///////
  double inv_sum = pow(1 - mu_sum, -1);
  double linDot = arma::dot(mu, linearTerm);
  double nonLinear = 0;
  for (unsigned int row = 0; row < d; row++)
    nonLinear += Dstar(inv_sum * (objectiveMean(row) - arma::dot(mu, constraintMean.row(row))));

  return(constantTerm + linDot - (1 - mu_sum) * nonLinear > 0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool DUST_MD::dualMaxAlgo2(const double& minCost, const unsigned int& t,
                            const unsigned int& s,
                            std::vector<unsigned int> l,
                            std::vector<unsigned int> r)
{
  return(false);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// COORDINATE DESCENT
//// COORDINATE DESCENT
//// COORDINATE DESCENT

bool DUST_MD::dualMaxAlgo3(const double& minCost,
                           const unsigned int& t,
                           const unsigned int& s,
                           std::vector<unsigned int> l,
                           std::vector<unsigned int> r)
{
  unsigned int l_size = l.size();
  unsigned int r_size = r.size();
  unsigned int nb_l_r = l_size + r_size; /// not equal to nb_max (if no size for getting constraints)

  ///////  ///////  ///////  ///////  ///////  ///////
  /////// DUAL parameter CONSTRUCTION
  ///////  ///////  ///////  ///////  ///////  ///////
  ///////
  ///////  AAAAA SIGNS: vector of -1s and 1s depending on whether constraint is to the left or right resp.
  ///////
  std::vector<int> sign;
  sign.reserve(nb_l_r); //// number of constraints
  for (unsigned int i = 0; i < nb_l_r; i++)
    sign.push_back(i < l_size ? -1 : 1);
  ///////
  /////// BBBBB: constantTerm AND objectiveMean
  ///////
  constantTerm = - (minCost - costRecord[s]) / (t - s);
  for (unsigned int row = 0; row < d; row++){objectiveMean(row) = (cumsum(row, t) - cumsum(row, s)) / (t - s);}
  ///////
  /////// CCCCC: linearTerm AND constraintMean
  ///////
  linearTerm.resize(nb_l_r);
  constraintMean.resize(d, nb_l_r);
  for (unsigned int j = 0; j < l_size ; j++)
  {
    linearTerm(j) = (costRecord[s] - costRecord[l[j]]) / (s - l[j]);
    for (unsigned int row = 0; row < d; row++)
      {constraintMean(row, j) = (cumsum(row,s) - cumsum(row, l[j])) / (s - l[j]);}
  }
  for (unsigned int j = 0; j < r_size ; j++)
  {
    linearTerm(l_size + j) = (costRecord[s] - costRecord[r[j]]) / (s - r[j]);
    for (unsigned int row = 0; row < d; row++)
      {constraintMean(row, l_size + j) = (cumsum(row,s) - cumsum(row,r[j])) / (s - r[j]);}
  }
  ///////  ///////  ///////  ///////  ///////  ///////
  /////// DUAL parameter CONSTRUCTION END
  ///////  ///////  ///////  ///////  ///////  ///////

  ///////
  /////// mu INITIAL = 0
  ///////
  mu.resize(nb_l_r);
  for (unsigned int i = 0; i < nb_l_r; i++){mu(i) = 0;}
  arma::colvec a = arma::colvec(d, arma::fill::zeros);
  arma::colvec b = arma::colvec(d, arma::fill::zeros);
  double c = 1;
  double d;
  double e;
  double f = constantTerm;
  /// the index to consider in the dual for mu


  ///////  ///////  ///////  ///////  ///////  ///////
  /////// COORDINATE DESCENT MAIN LOOP
  ///////  ///////  ///////  ///////  ///////  ///////
  unsigned int k_dual;
  std::array<double, 2> argmaxMax = { 0, -std::numeric_limits<double>::infinity() };

  for (unsigned int k = 0; k < nb_Loops; k++)
  {
    //// THE index for mu to optimize
    k_dual  = nb_Loops % nb_l_r;
    //// UPDATE the coefficients a,b,c,d,e,f
    /// a and b
    for(unsigned int j = 0; j < nb_l_r; j++)
    {
      if(j != k_dual){for(unsigned int row = 0; row < d; row++){a(row) = a(row) + sign[j]*mu(j)*constraintMean(row,j);}}
    }
    for(unsigned int row = 0; row < d; row++){b(row) = sign[k_dual]*constraintMean(row,k_dual);}

    /// c, d, e, f
    for(unsigned int j = 0; j < nb_l_r; j++){c += sign[j]*mu(j);}
    c -= sign[k_dual]*mu(k_dual);
    d = sign[k_dual];
    e = sign[k_dual]*linearTerm[k_dual];
    for(unsigned int j = 0; j < nb_l_r; j++)
    {
      if(j != k_dual){f += sign[j]*mu(j)*linearTerm[j];}
    }

    argmaxMax = dual1D_ArgmaxMax(a, b, c, d, e, f);  //// FIND ARGMAX AND MAX
    if(argmaxMax[1] > 0){return true;}  /// early return and/or update mu
    mu(k_dual) = argmaxMax[0];          //// UPDATE MU
  }
  return(false);
}





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool DUST_MD::dualMaxAlgo4(const double& minCost, const unsigned int& t,
                           const unsigned int& s,
                           std::vector<unsigned int> l,
                           std::vector<unsigned int> r)
{

  // ############################################## //
  // ############################################## //
  // ######### // FIRST STEP AT (0, 0) // ######### //
  // ############################################## //
  // ############################################## //


  // ######### // PELT TEST // ######### //
  // Formula: Dst - D*(Sst)              //

  constantTerm = - (minCost - costRecord[s]) / (t - s); // Dst // !!! CAPTURED IN OPTIM !!! //

  arma::subview_col col_t = cumsum.col(t);
  arma::subview_col col_s = cumsum.col(s);

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
  unsigned l_size = l.size(); // !!! UPDATED IN OPTIM !!! //

  linearTerm.resize(l_size);
  constraintMean.resize(d, l_size);

  mu_max.resize(l_size);
  inv_max.resize(l_size);

  double mean_sum = std::accumulate(objectiveMean.begin(), objectiveMean.end(), 0.0) / d;

  // Initialize the constraint mean matrix, which contains the mean-vectors associated with each constraint as columns
  // + Initialize mu_max, which is based on the value of the value of each mean-vector compared to the value of the "objective" mean-vector
  unsigned int j = 0;
  for (auto k: l)
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
  arma::rowvec tangent_max(l_size); // !!! SHRUNK IN OPTIM !!! //
  double grad_max = -std::numeric_limits<double>::infinity(); // !!! UPDATED IN OPTIM !!! //
  unsigned int grad_argmax = 0; // !!! UPDATED IN OPTIM !!! //

  for (unsigned int row = 0; row < d; row++)
    nonLinearGrad(row) = DstarPrime(objectiveMean(row));

  grad.resize(l_size);

  // Grad and hyperplane tangent compuation (formula: linearTerm + nonLinear * t(1_p) + t(nonLinearGrad) * (Srs - Sst * t(1_p)))
  for (unsigned int col = 0; col < l_size; col++)
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

  mu.resize(l_size);
  for (unsigned int col = 0; col < l_size; col++)
    mu(col) = 0;

  double mu_sum = 0; // !!! UPDATED IN OPTIM !!! //
  double inv_sum = 0; // !!! UPDATED IN OPTIM !!! //

  inverseHessian.resize(l_size, l_size);
  arma::dmat I = Identity.submat(0, 0, l_size - 1, l_size - 1);

  // Initialize inverseHessian as minus identity
  for (unsigned int row = 0; row < l_size; row++)
    for(unsigned int col = 0; col < l_size; col++)
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
        l_size -= zero_index.size();
        if (l_size <= 0) { return false; }

        grad_argmax = 0;
        for (auto k = zero_index.rbegin(); k != zero_index.rend(); ++k)
        {
          l.erase(l.begin() + *k);

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
        I = Identity.submat(0, 0, l_size - 1, l_size - 1);
      }

      bool shrink = false;
      std::vector<unsigned int> shrink_indices;
      shrink_indices.reserve(l_size);

      arma::rowvec direction(l_size); // direction and intensity of the update
      double direction_scale; // scaling applied to direction when direction pushes past boundaries
      arma::rowvec mu_diff(l_size); // (mu+) - mu
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
        for (unsigned int col = 0; col < l_size; col++)
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
        arma::rowvec gradCondition(l_size);
        for (unsigned int col = 0; col < l_size; col++)
          gradCondition(col) = m1 * grad(col);

        // 3. Scale dk until test is valid
        unsigned int iter = 0;
        while(new_test < test_value + arma::dot(mu_diff, gradCondition) && iter < 100)
        {
          for(unsigned int col = 0; col < l_size; col++)
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
        for (unsigned int col = 0; col < l_size; col++)
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
        for (unsigned int col = 0; col < l_size; col++)
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
        for (unsigned int col = 0; col < l_size; col++)
          dot_product += grad(col) * (tangent_max(col) - mu(col));

        if (test_value + dot_product <= 0) { return false; } // check if the tangent hyperplane ever reaches 0 on the simplex triangle (from concavity property of the dual)

        tangent_max(grad_argmax) = 0;

        // trigger shrinking, which reduces the dimension of the solution space by one
        if (shrink) return optim(shrink_indices);

        // Update inverse hessian estimation
        updateHessian(inverseHessian, mu_diff, grad_diff, I);
        grad_diff = -grad;
        iter++;
      } while (iter < nb_Loops);

      // Rcout << "Reached iter limit" << std::endl;

      return false;
    };

    // return false;
    return optim(std::vector<unsigned int>());
}

bool DUST_MD::dualMaxAlgo42(const double& minCost, const unsigned int& t,
                           const unsigned int& s,
                           std::vector<unsigned int> l,
                           std::vector<unsigned int> r)
{

  // ############################################## //
  // ############################################## //
  // ######### // FIRST STEP AT (0, 0) // ######### //
  // ############################################## //
  // ############################################## //

  Rcout << "testing" << std::endl;


  // ######### // PELT TEST // ######### //
  // Formula: Dst - D*(Sst)              //

  constantTerm = - (minCost - costRecord[s]) / (t - s); // Dst // !!! CAPTURED IN OPTIM !!! //

  auto col_t = cumsum.col(t);
  auto col_s = cumsum.col(s);

  double nonLinear = 0; // D*(Sst) // !!! UPDATED IN OPTIM !!! //
  for (unsigned int row = 0; row < d; row++)
  {
    objectiveMean(row) = (col_t(row) - col_s(row)) / (t - s);
    nonLinear += Dstar(objectiveMean(row));
  }

  double test_value = constantTerm - nonLinear; // !!! UPDATED IN OPTIM !!! //

  Rcout << "Pelt test" << std::endl;
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
  const unsigned l_size = l.size(); // !!! UPDATED IN OPTIM !!! //
  const unsigned r_size = r.size();
  const unsigned total_size = l_size + r_size;

  std::vector<unsigned> a;
  a.reserve(total_size);
  std::copy(l.begin(), l.end(), std::back_inserter(a));
  std::copy(r.begin(), r.end(), std::back_inserter(a));

  // vector of -1s and 1s depending on whether constraint is to the left or right resp.
  std::vector<int> sign;
  sign.reserve(total_size);
  for (auto i = 0; i < total_size; i++)
    sign.push_back(i < l_size ? -1 : 1);

  linearTerm.resize(total_size);
  constraintMean.resize(d, total_size);

  // Initialize the constraint mean matrix, which contains the mean-vectors associated with each constraint as columns
  // + Initialize mu_max, which is based on the value of the value of each mean-vector compared to the value of the "objective" mean-vector
  unsigned int j = 0;
  for (auto k: a)
  {
    linearTerm(j) = (costRecord[s] - sign[j] * costRecord[k]) / (s - k);
    auto col_k = cumsum.col(k);
    for (unsigned int row = 0; row < d; row++)
    {
      constraintMean(row, j) = sign[j] * (col_s(row) - col_k(row)) / (s - k);
    }
    j++;
  }

  Rcout << "Linear and constraint mean" << std::endl;

  for (unsigned int row = 0; row < d; row++)
    nonLinearGrad(row) = DstarPrime(objectiveMean(row));

  grad.resize(total_size);

  bool neg_grad = true; // !!! UPDATED IN OPTIM !!! //
  // Grad and hyperplane tangent compuation (formula: linearTerm + nonLinear * t(1_p) + t(nonLinearGrad) * (Srs - Sst * t(1_p)))
  for (unsigned int col = 0; col < total_size; col++)
  {
    double dot_product = 0;
    auto col_k = constraintMean.col(col);
    for (unsigned int row = 0; row < d; row++)
    {
      dot_product += nonLinearGrad(row) * (col_k(row) - objectiveMean(row));
    }
    grad(col) = - sign[col] * (nonLinear + dot_product + linearTerm(col));


    if (grad(col) > 0)
    {
      neg_grad = false;
    }
  }

  Rcout << "grad" << grad << std::endl;

  if (neg_grad)
  {
    return false;
  }

  // ############################################## //
  // ############################################## //
  // ######### // OPTIM INITIALIZATION // ######### //
  // ############################################## //
  // ############################################## //
  // Initialize all dynamic objects used in the op- //
  // timization recursion.                          //
  // ############################################## //

  mu.resize(total_size);
  for (unsigned int col = 0; col < total_size; col++)
    mu(col) = 0;

  double mu_sum = 0; // !!! UPDATED IN OPTIM !!! //
  double inv_sum = 0; // !!! UPDATED IN OPTIM !!! //

  inverseHessian.resize(total_size, total_size);
  arma::dmat I = Identity.submat(0, 0, total_size - 1, total_size - 1);

  // Initialize inverseHessian as minus identity
  for (unsigned int row = 0; row < total_size; row++)
    for(unsigned int col = 0; col < total_size; col++)
      inverseHessian(row, col) = row == col ? -1 : 0;

    // // ######################################### //
    // // ######################################### //
    // // ######### // OPTIM RECURSION // ######### //
    // // ######################################### //
    // // ######################################### //


    arma::rowvec direction(total_size); // direction and intensity of the update
    double direction_scale; // scaling applied to direction when direction pushes past boundaries
    arma::rowvec mu_diff(l_size); // (mu+) - mu
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
      for (unsigned int col = 0; col < total_size; col++)
      {
        mu_diff(col) = direction(col) * direction_scale;

        mu(col) += mu_diff(col);
        mu_sum += sign[col] * mu_diff(col);

        linDot += mu(col) * linearTerm(col);
      }

      inv_sum = pow(1 - mu_sum, -1);
      nonLinear = 0;
      for (unsigned int row = 0; row < d; row++)
      {
        m_value(row) = inv_sum * (objectiveMean(row) + arma::dot(mu, constraintMean.row(row)));
        nonLinear += Dstar(m_value(row));
      }

      double new_test = -(1 + mu_sum) * nonLinear - linDot + constantTerm;

      // 2. define gradient condition
      arma::rowvec gradCondition(l_size);
      for (unsigned int col = 0; col < l_size; col++)
        gradCondition(col) = m1 * grad(col);

      // 3. Scale dk until test is valid
      unsigned int iter = 0;
      while(new_test < test_value + arma::dot(mu_diff, gradCondition) && iter < 100)
      {
        for(unsigned int col = 0; col < l_size; col++)
        {
          mu_diff(col) *= .5;
          mu(col) -= mu_diff(col);
          linDot -= mu_diff(col) * linearTerm(col);
          mu_sum -= sign[col] * mu_diff(col);
        }

        inv_sum = pow(1 - mu_sum, -1);
        nonLinear = 0;
        for (unsigned int row = 0; row < d; row++)
        {
          m_value(row) = inv_sum * (objectiveMean(row) + arma::dot(mu, constraintMean.row(row)));
          nonLinear += Dstar(m_value(row));
        }

        new_test = -(1 + mu_sum) * nonLinear - linDot + constantTerm; // update values

        iter++;
      }
      test_value = new_test;

      // Project mu onto the interior simplex
      // If any null value, then shrink
      for (unsigned int col = 0; col < l_size; col++)
      {
        if(mu(col) < 0) { mu_sum -= sign[col] * mu(col); mu(col) = 0; }
      }
    };

    auto updateGrad = [&] ()
    {
      for (unsigned int row = 0; row < d; row++)
        nonLinearGrad(row) = DstarPrime(m_value(row));

      for (unsigned int col = 0; col < l_size; col++)
      {
        double dot_product = 0;
        auto col_k = constraintMean.col(col);
        for (unsigned int row = 0; row < d; row++)
        {
          dot_product += nonLinearGrad(row) * (col_k(row) - objectiveMean(row));
        }
        grad(col) = -sign[col] * (nonLinear + dot_product + linearTerm(col));


        if (grad(col) > 0)
        {
          neg_grad = false;
        }
      }

      // Update grad value.
      for (unsigned int col = 0; col < l_size; col++)
      {
        double dot_product = 0;
        auto col_k = constraintMean.col(col);
        for (unsigned int row = 0; row < d; row++)
          dot_product += nonLinearGrad(row) * (col_k(row) - m_value(row));
        grad(col) = -sign[col] * (nonLinear + dot_product + linearTerm(col));
        grad_diff(col) += grad(col);
        if (grad_diff(col) == 0)
        {
          if (grad(col) > 0) { grad_diff(col) = -1e-16; }
          else { grad_diff(col) = 1e-16; }
        }
      }
    };

    // Optim loop
    unsigned int iter = 0;
    do
    {
      // Update direction and intensity then clip it to stay within the bounds of the positive simplex
      direction = (-grad) * inverseHessian;
      direction_scale = 1.;

      Rcout << "direction " << direction << std::endl;

      double direction_sum ;
      clip_stepsize_to_negative_element(mu, direction, direction_scale);
      clip_stepsize_to_negative_sum(sign, mu_sum, direction, direction_sum, direction_scale);

      for (auto row = 0; row < d; row++)
      {
        clipStepSizeModel(m_value(row), constraintMean.row(row), mu_sum, direction, direction_sum, direction_scale);
      }

      Rcout << "direction scale " << direction_scale << std::endl;

      if (direction_scale <= 0) { return false; } // stops optimization if no movement is produced

      // update mu and D(mu) + check for shrink
      updateTestValue();

      Rcout << "mu " << mu << std::endl;

      if(test_value > 0) { return true; } // success, index s is pruned

      // update grad and tangent max location
      updateGrad();

      Rcout << "grad " << grad << std::endl;

      // Update inverse hessian estimation
      updateHessian(inverseHessian, mu_diff, grad_diff, I);

      Rcout << "Inverse hessian " << inverseHessian << std::endl;

      grad_diff = -grad;
      iter++;
    } while (iter < nb_Loops);

    // return false;
    return false;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

bool DUST_MD::dualMaxAlgo5(const double& minCost, const unsigned int& t,
                           const unsigned int& s,
                           std::vector<unsigned int> l,
                           std::vector<unsigned int> r)
{
  // ############################################## //
  // ############################################## //
  // ######### // FIRST STEP AT (0, 0) // ######### //
  // ############################################## //
  // ############################################## //


  // ######### // PELT TEST // ######### //
  // Formula: Dst - D*(Sst)              //

  constantTerm = - (minCost - costRecord[s]) / (t - s); // Dst // !!! CAPTURED IN OPTIM !!! //

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


bool DUST_MD::dualMaxAlgo6(const double& minCost, const unsigned int& t,
                           const unsigned int& s,
                           std::vector<unsigned int> l,
                           std::vector<unsigned int> r)
{
  return(false);
}

bool DUST_MD::dualMaxAlgo7(const double& minCost, const unsigned int& t,
                           const unsigned int& s,
                           std::vector<unsigned int> l,
                           std::vector<unsigned int> r)
{
  return false;
}



