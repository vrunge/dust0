#include "MD_A_DUST.h"
#include "fstream"

using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// --- // Constructor // --- //
DUST_MD::DUST_MD(int dual_max, bool random_constraint, Nullable<unsigned> nbLoops)
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

  /// Exception
  /// Exception
  if (!first_execution)
  {
    if (d != inData.n_rows)
      throw std::invalid_argument("new data has invalid n_rows. got " + std::to_string(inData.n_rows) + ", expected " + std::to_string(d));
  }
  else { d = inData.n_rows; }

  /// PENALTY value default
  if (inPenalty.isNull()){penalty = 2 * d * std::log(n);}else{penalty = as<double>(inPenalty);}

  /// RETURN
  changepointRecord.reserve(n + 1);
  nb_indices.reserve(n);
  costRecord.reserve(n + 1);

  /// CUMSUM => statistics
  cumsum.resize(d, n + 1);

  if (first_execution)
  {
    changepointRecord.push_back(0);
    nb_indices.push_back(1);
    costRecord.push_back(-penalty);
    for (unsigned row = 0; row < d; row++)
      cumsum(row, 0) = 0.;

    // read the number of constraints + default choice
    if (inNbL.isNull())
    {
      if (inNbR.isNotNull())
      {
        nb_r = std::min(d, as<unsigned>(inNbR));
        nb_l = d - nb_r;
      }
      else 
      {
        nb_r = 0;
        nb_l = d;
      }
    }
    else
    {
      nb_l = std::min(d, as<unsigned>(inNbL));
      nb_r = d - nb_l;
    }
    nb_max = nb_l + nb_r;
    Rcout << "nb_l = " << nb_l << "; nb_r = " << nb_r << "; nb_max = " << nb_max << std::endl;

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
  else{indices -> set_init_size(n);}

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
  /////
  ///// output
  /////
  //std::ofstream log_file;
  //log_file.close();
  //log_file.open("output.log", std::ofstream::out | std::ofstream::app);
  //log_file.clear();
  //auto cout_buf = Rcout.rdbuf(log_file.rdbuf());
  //Rcout << "correctly setup log_file" << std::endl;

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
    unsigned argMin = 0;
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
      ///// OUTPUTLOG
      //Rcout << "fetching constraints" << std::endl;
      if ((this->*current_test)(minCost,
                                t,
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
        indices -> next_prune();
      }
    }
    // END (DUST loop)
    // END (DUST loop)

    // Prune the last index (analogous with (mu* = 0) duality simple test in mu = zero)
    // if (lastCost > minCost)
    // {
    //   indices->prune_last();
    //   nbt--;
    // }


    indices->add(t);
    nb_indices.push_back(nbt);
    nbt++;
  }

  /////  OUTPUT
  //Rcout.rdbuf(cout_buf);
  //log_file.close();
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


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// arma::colvec objectiveMean;
// arma::dmat constraintMean;
// arma::rowvec linearTerm;
// double constantTerm;
// arma::rowvec mu;
//
// COMPUTE :  - [(1 + mu_sum) * nonLinear + linDot + constantTerm]
// sign in front of mu already included
//
double DUST_MD::dual_Eval()
{
  double mu_sum = 0;
  for (unsigned int i = 0; i < mu.n_elem; i++){mu_sum += mu(i);}
  double coeff = pow(1 + mu_sum, -1);

  double Linear = arma::dot(mu, linearTerm);
  double nonLinear = 0;
  for (unsigned int i = 0; i < d; i++)
    nonLinear += Dstar(coeff * (objectiveMean(i) + arma::dot(mu, constraintMean.row(i))));

  return( - ((1 + mu_sum) * nonLinear + Linear + constantTerm));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double DUST_MD::dual_Eval(double &nonLinear)
{
  double mu_sum = std::accumulate(mu.begin(), mu.end(), 0.);
  double coeff = pow(1+ mu_sum, -1);

  double Linear = arma::dot(mu, linearTerm);
  nonLinear = 0;
  for (unsigned int i = 0; i < d; i++)
    nonLinear += Dstar(coeff * (objectiveMean(i) + arma::dot(mu, constraintMean.row(i))));

  return( - ((1 + mu_sum) * nonLinear + Linear + constantTerm));
}

void DUST_MD::grad_Eval(const double nonLinear)
{  
  for (unsigned col = 0; col < mu.size(); col++)
  {
    double dot_product = 0.;
    for (unsigned row = 0; row < m_value.size(); row++)
    {
      dot_product += DstarPrime(m_value(row)) * (constraintMean(row, col) - m_value(row));
    }
    grad(col) = nonLinear + dot_product + linearTerm(col);
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void DUST_MD::update_dual_parameters_l(const double& minCost,
                                       const unsigned int& t,
                                       const unsigned int& s,
                                       std::vector<unsigned int>& l)
{
  unsigned size_l = l.size();

  /// RESIZE
  mu.resize(size_l);
  mu_max.resize(size_l);
  mu_max.fill(1.0);
  linearTerm.resize(size_l);
  constraintMean.resize(d, size_l);

  /// UDDATE DUAL FUNCTION parameters
  constantTerm =  (minCost - costRecord[s]) / (t - s);
  for (unsigned int row = 0; row < d; row++)
  { 
    objectiveMean(row) = (cumsum(row, t) - cumsum(row, s)) / (t - s);
  }

  for (unsigned int j = 0; j < size_l; j++)
  {
    linearTerm(j) = (costRecord[s] - costRecord[l[j]]) / (s - l[j]);
    for (unsigned int row = 0; row < d; row++)
    {
      constraintMean(row, j) = (cumsum(row,s) - cumsum(row, l[j])) / (s - l[j]);
      mu_max(j) = std::min(mu_max(j), muMax(objectiveMean(row), constraintMean(row, j)));
    }
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void DUST_MD::update_dual_parameters_l_r(const double& minCost,
                                        const unsigned int& t,
                                        const unsigned int& s,
                                        std::vector<unsigned int>& l,
                                        std::vector<unsigned int>& r)
{
  unsigned size_l = l.size();
  unsigned size_r = r.size();

  /// RESIZE
  mu.resize(size_l + size_r);
  mu_max.resize(size_l + size_r);
  mu_max.fill(1.0);
  linearTerm.resize(size_l + size_r);
  constraintMean.resize(d, size_l + size_r);

  /// UDDATE DUAL FUNCTION parameters
  constantTerm =  (minCost - costRecord[s]) / (t - s);
  for (unsigned int row = 0; row < d; row++){objectiveMean(row) = (cumsum(row, t) - cumsum(row, s)) / (t - s);}
  ///////
  for (unsigned int j = 0; j < size_l; j++)
  {
    linearTerm(j) = (costRecord[s] - costRecord[l[j]]) / (s - l[j]);
    for (unsigned int row = 0; row < d; row++)
    {
      constraintMean(row, j) = (cumsum(row,s) - cumsum(row, l[j])) / (s - l[j]);
    }
  }
  ///////
  /// BIG DANGER = difference of two unsigned INT
  /// BIG DANGER = difference of two unsigned INT
  /// BIG DANGER = difference of two unsigned INT
  for (unsigned int j = 0; j < size_r; j++)
  {
    linearTerm(size_l + j) = -(costRecord[s] - costRecord[r[j]]) / (r[j] - s);
    for (unsigned int row = 0; row < d; row++)
    {
      constraintMean(row, size_l + j) = -(cumsum(row,s) - cumsum(row, r[j])) / (r[j] - s);
    }
  }

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Draw a random point and evaluate the corresponding dual value

bool DUST_MD::dualMaxAlgo0(const double& minCost,
                           const unsigned int& t,
                           const unsigned int& s,
                           std::vector<unsigned int> l,
                           std::vector<unsigned int> r)
{
  update_dual_parameters_l(minCost, t, s, l);

  /// CHOOSE ONE MU RANDOMLY. Only using left constraints
  // Random vector u in the SIMPLEX with boundary mu(i) = mu_max(i) on the i-th axis
  std::vector<double> u;
  u.reserve(l.size());
  for (unsigned int i = 0; i < l.size(); i++){u.push_back(dist(engine));}

  double sum_all = 0;
  double scaling_factor = dist(engine); // Uniform random number in [0, 1]
  for (unsigned int i = 0; i < l.size(); i++){sum_all = sum_all + u[i]/mu_max(i);}
  for (unsigned int i = 0; i < l.size(); i++){mu(i) = -scaling_factor/sum_all*u[i];} // minus for r < s constraints

  return(dual_Eval() > 0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////
////////// random direction in l constraints + 1D dual optimization
//////////

bool DUST_MD::dualMaxAlgo1(const double& minCost, const unsigned int& t,
                           const unsigned int& s,
                           std::vector<unsigned int> l,
                           std::vector<unsigned int> r)
{
  update_dual_parameters_l_r(minCost, t, s, l, r);
  /// CHOOSE ONE MU RANDOMLY
  // Random direction u in the positive orthant with sign included:
  /// (-1) mu for left constraint
  /// (+1) mu for right constraint
  std::vector<double> u;
  u.reserve(l.size() + r.size());

  Rcout << "ALL INDEX "<< std::endl;
  for (unsigned int i = 0; i < l.size(); i++){Rcout << "l: " << l[i]<< std::endl;}
  Rcout << "s: " << s << std::endl;
  for (unsigned int i = 0; i < r.size(); i++){Rcout << "r: "  << r[i]<< std::endl;}
  Rcout << t<< std::endl;

  Rcout << "ALL INDEX "<< std::endl;
  for (unsigned int i = 0; i < nb_l; i++){Rcout << "l: " << l[i]<< std::endl;}
  Rcout << "s: " << s << std::endl;
  for (unsigned int i = 0; i < nb_r; i++){Rcout << "r: "  << r[i]<< std::endl;}
  Rcout << t<< std::endl;


  /// push_back => START by the end
  for (unsigned int i = 0; i < l.size(); i++){u.push_back(-dist(engine));}
  for (unsigned int i = 0; i < r.size(); i++){u.push_back(dist(engine));}

  //
  // build the 1D dual
  //
  arma::colvec b(d);
  for (unsigned int row = 0; row < d; row++)
  {
    b(row) = 0;
    for (unsigned int i = 0; i < u.size(); i++)
    {
      b(row) += u[i] * constraintMean(row, i);
    }
  }
  double D = 0;
  double e = 0;
  for (unsigned int i = 0; i < u.size(); i++)
  {
    D += u[i];
    e += u[i]*linearTerm(i);
  }
  double argmax = 0; // initial value for finding max
  //
  // optimize the 1D dual with dual1D_Max(argmax,a,b,c,d,e,f)
  // a = objectiveMean
  double c = 1;
  // f = constantTerm
  return(dual1D_Max(argmax,objectiveMean,b,c,D,e,constantTerm) > 0);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// BARYCENTRE test
// BARYCENTRE test
// BARYCENTRE test

bool DUST_MD::dualMaxAlgo2(const double& minCost,
                           const unsigned int& t,
                           const unsigned int& s,
                           std::vector<unsigned int> l,
                           std::vector<unsigned int> r)
{

  update_dual_parameters_l_r(minCost, t, s, l, r);

  ///////
  /////// mu
  ///////
  mu.resize(l.size() + r.size());
  double mu_sum = 0;
  double x = pow(l.size() + r.size() + 1, -1); ///

  for (unsigned int i = 0; i < l.size(); i++) ///////// WITH l
  {
    mu(i) = mu_max(i) * x;
    mu_sum += mu(i);
  }
  for (unsigned int i = 0; i < r.size(); i++) ///////// WITH r
  {
    mu(l.size() + i) = - mu_max(l.size() + i) * x; /// with a minus here
    mu_sum = mu(l.size() + i);
  }

  //Rcout << "test"  <<  std::endl;
  //for (unsigned int i = 0; i < nb_l; i++){Rcout << mu(i) << " -- ";}
  //Rcout << std::endl;
  //Rcout << nb_l << std::endl;
  //Rcout << dual_Eval() << std::endl;

  return(dual_Eval() > 0);
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
  //Rcout << "DUST_MD DUST_MD DUST_MD DUST_MD DUST_MD " << std::endl;

  update_dual_parameters_l_r(minCost, t, s, l, r);
  unsigned int nb_l_r = l.size() + r.size();

  std::vector<int> sign;
  sign.reserve(nb_l_r); //// number of constraints
  for (unsigned int i = 0; i < nb_l_r; i++)
    sign.push_back(i < l.size() ? -1 : 1);


  ///////
  /////// 1D dual parameters
  ///////
  arma::colvec a = arma::colvec(d, arma::fill::zeros);
  arma::colvec b = arma::colvec(d, arma::fill::zeros);
  double c = 1;
  double D;
  double e;
  double f = constantTerm;


  ///////
  /////// mu INITIAL = 0
  ///////
  for (unsigned int i = 0; i < nb_l_r; i++){mu(i) = 0;}


  ///////  ///////  ///////  ///////  ///////  ///////
  /////// COORDINATE DESCENT MAIN LOOP
  ///////  ///////  ///////  ///////  ///////  ///////
  unsigned int k_dual;
  double Max = -std::numeric_limits<double>::infinity();
  double argmax = 0;

  for (unsigned int k = 0; k < nb_Loops; k++)
  {
    for(unsigned int row = 0; row < d; row++)
    {
      a(row) = objectiveMean(row);
    }
    c = 1;
    f = constantTerm;


    //// THE index for mu to optimize
    k_dual  = k % nb_l_r;
    //// UPDATE the coefficients a,b,c,d,e,f
    /// a and b
    for(unsigned int j = 0; j < nb_l_r; j++)
    {
      if(j != k_dual){for(unsigned int row = 0; row < d; row++){a(row) = a(row) + sign[j]*mu(j)*constraintMean(row,j);}}
    }
    for(unsigned int row = 0; row < d; row++){b(row) = sign[k_dual]*constraintMean(row,k_dual);}

    /// c, d, e, f
    for(unsigned int j = 0; j < nb_l_r; j++)
    {
      if(j != k_dual){c += sign[j]*mu(j);}
    }
    D = sign[k_dual];
    e = sign[k_dual]*linearTerm[k_dual];
    for(unsigned int j = 0; j < nb_l_r; j++)
    {
      if(j != k_dual){f += sign[j]*mu(j)*linearTerm[j];}
    }

    Max = dual1D_Max(argmax, a, b, c, D, e, f);  //// FIND ARGMAX AND MAX
    //Rcout << "dual1D_Max " << std::endl;
    //Rcout << Max << std::endl;
    //Rcout << argmax << std::endl;
    //Rcout << f << " / " << constantTerm << std::endl;
    if(Max > 0){return true;}  /// early return and/or update mu
    mu(k_dual) = argmax;          //// UPDATE MU
  }
  return(false);
}





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool DUST_MD::dualMaxAlgo4(const double &minCost, const unsigned int &t,
                           const unsigned int &s,
                           std::vector<unsigned int> l,
                           std::vector<unsigned int> r)
{

  // ############################################## //
  // ############################################## //
  // ######### // FIRST STEP AT (0, 0) // ######### //
  // ############################################## //
  // ############################################## //


  // ######### // PELT TEST // ######### //
  // Formula: Dst - D*(Sst) //


  Rcout << std::endl << "t: " << t << "; s: " << s << std::endl;
  Rcout << "l: " << std::endl;
  for (auto lk: l)
    Rcout << lk << "; ";
  Rcout << "r: " << std::endl;
  for (auto rk: r)
    Rcout << rk << "; ";

  update_dual_parameters_l(minCost, t, s, l);
  Rcout << "updated parameters..." << std::endl;

  double nonLinear = 0; // D*(Sst) // !!! UPDATED IN OPTIM !!! //
  for (unsigned int row = 0; row < d; row++)
  {
    nonLinear += Dstar(objectiveMean(row));
  }
  Rcout << "computed nonLinear ..." << std::endl;

  double test_value = -(constantTerm + nonLinear); // !!! UPDATED IN OPTIM !!! //

  if (test_value > 0) { return true; } // PELT test
  return false;
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

  // Coordinates of the point that maximizes the hyperplane tangent defined at the current mu (init)
  arma::rowvec tangent_max(l_size); // !!! SHRUNK IN OPTIM !!! //
  double grad_max = -std::numeric_limits<double>::infinity(); // !!! UPDATED IN OPTIM !!! //
  unsigned int grad_argmax = 0; // !!! UPDATED IN OPTIM !!! //

  for (unsigned int row = 0; row < d; row++)
    nonLinearGrad(row) = DstarPrime(objectiveMean(row));

  mu.fill(0.);
  grad.resize(l_size);
  grad_Eval(nonLinear);

  // Grad and hyperplane tangent compuation (formula: linearTerm + nonLinear * t(1_p) + t(nonLinearGrad) * (Srs - Sst * t(1_p)))
  for (unsigned int col = 0; col < l_size; col++)
  {
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
    if (test_value + arma::dot(tangent_max, grad) <= 0)
    {
      return false;
    } // check if the tangent hyperplane ever reaches 0 on the simplex triangle (from concavity property of the dual)
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

  inverseHessian.resize(l_size, l_size);
  arma::dmat I = Identity.submat(0, 0, l_size - 1, l_size - 1);

  // Initialize inverseHessian as minus identity
  for (unsigned int row = 0; row < l_size; row++)
  {
    for(unsigned int col = 0; col < l_size; col++)
    {
      if (row == col) inverseHessian(row, col) = -1;
      else inverseHessian(row, col) = 0;
    }
  }

  // // ######################################### //
  // // ######################################### //
  // // ######### // OPTIM RECURSION // ######### //
  // // ######################################### //
  // // ######################################### //

  std::function<bool(std::vector<unsigned int>)> optim = [&] (std::vector<unsigned int> &&zero_index)
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
      for (unsigned int col = 0; col < l_size; col++)
      {
        mu_diff(col) = direction(col) * direction_scale;
        mu(col) -= mu_diff(col);
      }

      double new_test = dual_Eval(nonLinear);

      // 2. define gradient condition
      arma::rowvec gradCondition = m1 * grad;

      // 3. Scale dk until test is valid
      unsigned iter = 0;
      while(new_test < test_value + arma::dot(mu_diff, gradCondition) && ++iter < 100)
      {
        for(unsigned col = 0; col < l_size; col++)
        {
          mu_diff(col) *= .5;
          mu(col) += mu_diff(col);
        }
        new_test = dual_Eval(nonLinear);
      }
      test_value = new_test;

      // Project mu onto the interior simplex
      // If any null value, then shrink
      for (unsigned int col = 0; col < l_size; col++)
      {
        if(mu(col) > 0) { mu(col) = 0.; shrink = true; shrink_indices.push_back(col); }
        else if(mu(col) == 0.) { shrink = true; shrink_indices.push_back(col); }
      }
    };

    auto updateGrad = [&] ()
    {
      // Initialize tangent max search.
      grad_max = -std::numeric_limits<double>::infinity();
      grad_argmax = 0;

      grad_Eval(nonLinear);

      // Update grad value.
      for (unsigned int col = 0; col < l_size; col++)
      {
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
      {
        tangent_max(grad_argmax) = mu_max(grad_argmax);
      }
    };

    // Optim loop
    unsigned iter = 0;
    do
    {
      // Update direction and intensity then clip it to stay within the bounds of the positive simplex
      direction = (-grad) * inverseHessian;
      direction_scale = FindBoundaryCoef(mu, direction, inv_max, shrink, shrink_indices); // may trigger shrink if pushing past boundary
      if (shrink)
      {
        return optim(std::move(shrink_indices));
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
        dot_product += grad(col) * (tangent_max(col) + mu(col));

      if (test_value + dot_product <= 0) { return false; } // check if the tangent hyperplane ever reaches 0 on the simplex triangle (from concavity property of the dual)

      tangent_max(grad_argmax) = 0;

      // trigger shrinking, which reduces the dimension of the solution space by one
      if (shrink) return optim(std::move(shrink_indices));

      // Update inverse hessian estimation
      updateHessian(inverseHessian, mu_diff, grad_diff, I);
      grad_diff = -grad;
    } while (++iter < nb_Loops);

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
  for (unsigned i = 0; i < total_size; i++)
  {
    sign.push_back(i < l_size ? -1 : 1);
  }

  linearTerm.resize(total_size);
  constraintMean.resize(d, total_size);

  // Initialize the constraint mean matrix, which contains the mean-vectors associated with each constraint as columns
  // + Initialize mu_max, which is based on the value of the value of each mean-vector compared to the value of the "objective" mean-vector
  unsigned int j = 0;
  for (auto k: a)
  {
    linearTerm(j) = sign[j] * (costRecord[s] - costRecord[k]) / (s - k);
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
  {
    for(unsigned int col = 0; col < total_size; col++)
      inverseHessian(row, col) = row == col ? -1 : 0;
  }

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

    inv_sum = pow(1 + mu_sum, -1);
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
    while(new_test < test_value + arma::dot(mu_diff, gradCondition) && iter < nb_Loops)
    {
      for(unsigned int col = 0; col < l_size; col++)
      {
        mu_diff(col) *= .5;
        mu(col) -= mu_diff(col);
        linDot -= mu_diff(col) * linearTerm(col);
        mu_sum -= sign[col] * mu_diff(col);
      }

      inv_sum = pow(1 + mu_sum, -1);
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
        dot_product += nonLinearGrad(row) * (constraintMean(row, col) - objectiveMean(row));
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

    for (unsigned row = 0; row < d; row++)
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



