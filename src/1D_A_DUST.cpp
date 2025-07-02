#include <Rcpp.h>
#include <cmath>
#include <random> /// FOR RANDOM NUMBER IN DUAL EVAL
#include "1D_A_DUST.h"
#include "preProcessing.h"

#include <fstream> /////// T0 WRITE  EMMELINE
#include <iostream>  /////// T0 WRITE  EMMELINE

using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// --- // Constructor // --- //
DUST_1D::DUST_1D(int dual_max_type, int constraints_type, Nullable<int> nbLoops)
  : dual_max_type(dual_max_type),
    constraints_type(constraints_type),
    indices(nullptr),
    n(0)
{
  if(nbLoops.isNull()){nb_Loops = 10;}else{nb_Loops = as<int>(nbLoops);}
}

////////////////////////////////////////////////////////////////////////////////

DUST_1D::~DUST_1D()
{
  delete indices;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void DUST_1D::pruning_method()
{
  delete indices;

  /// /// ///
  /// /// /// index METHOD
  /// /// ///
  if(constraints_type == 0){indices = new RandomIndices_1D();}
  else{indices = new DeterministicIndices_1D;}

  /// /// ///
  /// /// /// dual_max_type METHOD
  /// /// ///
  if(dual_max_type == 0){current_test = &DUST_1D::dualMaxAlgo0;}
  if(dual_max_type == 1){current_test = &DUST_1D::dualMaxAlgo1;}
  if(dual_max_type == 2){current_test = &DUST_1D::dualMaxAlgo2;}
  if(dual_max_type == 3){current_test = &DUST_1D::dualMaxAlgo3;}
  if(dual_max_type == 4){current_test = &DUST_1D::dualMaxAlgo4;}
  if(dual_max_type == 5){current_test = &DUST_1D::dualMaxAlgo5;}
  if(dual_max_type == 6){current_test = &DUST_1D::dualMaxAlgo6;}

  /// /// ///
  /// /// /// INIT RANDOM GENERATOR
  /// /// ///
  engine.seed(std::random_device{}());
  dist = std::uniform_real_distribution<double>(0.0, 1.0);
}


// --- // Fits the data, i. e. initializes all data-dependent vectors // --- //
// --- // Fits the data, i. e. initializes all data-dependent vectors // --- //
// --- // Fits the data, i. e. initializes all data-dependent vectors // --- //

void DUST_1D::append(std::vector<double>& inData, Nullable<double> inPenalty)
{
  bool first_execution = (n == 0);

  // update total size
  n += inData.size();

  // Reserve memory for all vectors
  cumsum.reserve(n + 1);
  changepointRecord.reserve(n + 1);
  costRecord.reserve(n + 1);
  nb_indices.reserve(n + 1);

  if (first_execution)
  {
    // On first execution: initialize all vectors
    if (inPenalty.isNull()) { penalty = 2 * std::log(n); } else { penalty = as<double>(inPenalty); }

    cumsum.push_back(0);
    costRecord.push_back(-penalty);
    changepointRecord.push_back(0);
    nb_indices.push_back(1);

    pruning_method();
    indices->add(0);
  }

  // store data by updating the cumsum
  for (auto current_datum = inData.begin(); current_datum != inData.end(); ++current_datum)
    cumsum.push_back(cumsum.back() + statistic(*current_datum));
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// --- // Algorithm-specific method // --- //
void DUST_1D::update_partition()
{
  // Initialize OP step value
  double lastCost;
  unsigned int nbt = nb_indices.back();

  // Main loop
  for (unsigned t = indices->get_first() + 1; t <= n; t++)
  {
    // OP step
    // OP step
    indices->reset();
    double minCost = std::numeric_limits<double>::infinity();
    unsigned argMin = 0;
    do
    {
      unsigned int s = indices->get_current();
      lastCost = costRecord[s] + Cost(t, s); // without the penalty beta
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

    // OP update minCost and save values
    // minCost = Q_t. Here + beta to get Q_t = Q_i + C(y_it) + beta
    minCost += penalty;
    costRecord.push_back(minCost);
    changepointRecord.push_back(argMin);

    // DUST step
    // DUST step
    indices->reset_prune();

    // DUST loop
    // DUST loop
    while (indices->check_prune()) // is true, while we are not on the smallest index
    {
      // prune as needs pruning
      if ((this->*current_test)(minCost, t, indices->get_current(), indices->get_constraint()))
      {
        // remove the pruned index and its pointer
        // removing the elements increments the cursors i and pointersCurrent, while before stands still
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

    // Prune the last index (analogous with a "mu* = 0" duality simple test)
    // this is the smallest available index = PELT RULE
    if (lastCost > minCost)
    {
      indices->prune_last();
      nbt--;
    }

    // update the available indices
    indices->add(t);
    nb_indices.push_back(nbt);
    nbt++;
  }
}


// --- // Wrapper method for quickly computing               // --- //
// --- // and retrieving the optimal partition of input data // --- //
List DUST_1D::one_dust(std::vector<double>& inData, Nullable<double> inPenalty)
{
  append(inData, inPenalty);
  update_partition();
  return get_partition();
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// --- // Builds changepoints // --- //
// --- // Builds changepoints // --- //
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
// --- // Retrieves optimal partition // --- //
// --- // Retrieves optimal partition // --- //
List DUST_1D::get_partition()
{
  std::forward_list<unsigned int> chpts = backtrack_changepoints();

  return List::create(
    _["changepoints"] = chpts,
    _["lastIndexSet"] = indices->get_list(),
    _["nb"] = std::vector<unsigned>(nb_indices.begin() + 1, nb_indices.end()),
    _["costQ"] = std::vector<double>(costRecord.begin() + 1, costRecord.end())
  );
}

// --- // get_info object // --- //
// --- // get_info object // --- //
// --- // get_info object // --- //
List DUST_1D::get_info()
{
  return List::create(
    _["data_statistic"] = cumsum,
    _["data_length"] = n,
    _["current_penalty"] = penalty,
    _["model"] = get_model(),
    _["pruning_algo"] = dual_max_type,
    _["pruning_constraints_type"] = constraints_type,
    _["pruning_nb_loops"] = nb_Loops
  );
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// SPECIAL CASE BEFORE DUAL EVALUATION / OPTIMIZATION

bool DUST_1D::isSpecialCase(double objectiveMean, double constraintMean)
{
  if(isBoundary(objectiveMean) == true){return true;}
  if(objectiveMean == constraintMean){return true;}
  return false;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// if special_cases == true => return true;

bool DUST_1D::specialCasePruning(double objectiveMean, double constraintMean,
                                 double linearTerm, double constantTerm, double mu_max)
{
  ///
  /// case "objectiveMean = 0"
  ///
  if(isBoundary(objectiveMean) == true)
  {
    if(isBoundary(constraintMean) == false) /// case mu_max = 0
    {
      if (constantTerm > 0) {return true;}
      return false;
    }
    else /// case linear
    {
      if (constantTerm > 0){return true;} /// left bound
      /// case bern / binom : mu_max != 1
      if (mu_max * linearTerm + constantTerm > 0){return true;} /// right bound
      return false;
    }
  }

  ///
  /// case isBoundary(objectiveMean) == false AND isBoundary(constraintMean) == true
  ///  SEE MU_MAX function
  ///

  ///
  /// case "objectiveMean = constraintMean", here not on the boundary
  /// ALWAYS mu_max = 1
  ///
  if(objectiveMean == constraintMean)
  {
    if (- Dstar(objectiveMean) + constantTerm > 0){return true;} /// left bound
    if (linearTerm + constantTerm > 0){return true;} /// right bound
  }
  return false;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// RANDOM EVAL

bool DUST_1D::dualMaxAlgo0(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  return (dualEval(dist(engine), minCost, t, s, r) > 0);
}

/// EXACT EVAL // only available with Gaussian cost (fixed variance)

bool DUST_1D::dualMaxAlgo1(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  double a = (cumsum[t] - cumsum[s]) / (t - s);
  double b = (cumsum[s] - cumsum[r]) / (s - r);
  double c = (costRecord[s] - costRecord[r]) / (s - r);
  double d = (costRecord[s] - minCost) / (t - s);
  double mu_max = muMax(a, b);

  ///
  /// case dual domain = one point OR dual = linear function
  ///
  if(isSpecialCase(a,b) == true){return specialCasePruning(a,b,c,d,mu_max);}
  return (dualMax(minCost, t, s, r) > 0);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// Golden-section search

bool DUST_1D::dualMaxAlgo2(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  double a = (cumsum[t] - cumsum[s]) / (t - s);
  double b = (cumsum[s] - cumsum[r]) / (s - r);

  double lt = 0.0;
  double rt = muMax(a, b);
  double c = (1 - 1/phi) * rt;
  double d = 1/phi;

  double linear = (costRecord[s] - costRecord[r]) / (s - r);
  double cst = (costRecord[s] - minCost) / (t - s);

  double fc = - (1.0 - c) * Dstar((a - c*b)/(1 - c)) + c * linear + cst;
  double fd = - (1.0 - d) * Dstar((a - d*b)/(1 - d)) + d * linear + cst;
  if(fc > 0 || fd > 0){return true;}
  double max_val = std::max(fc, fd);

  for (int i = 0; i < nb_Loops; i++)
  {
    if (fc > fd)
    {
      rt = d;
      d = c;
      fd = fc;
      c = rt - (rt - lt) / phi;
      fc = - (1.0 - c) * Dstar((a - c*b)/(1 - c)) + c * linear + cst;
    }
    else
    {
      lt = c;
      c = d;
      fc = fd;
      d = lt + (rt - lt) / phi;
      fd = - (1.0 - d) * Dstar((a - d*b)/(1 - d)) + d * linear + cst;
    }
    max_val = std::max(max_val, std::max(fc, fd));
    if(max_val > 0){return true;}
  }
  return false;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// binary search. At each step, we evaluate the tangent line to the current point
// at its max to stop the search at early step (when possible)

bool DUST_1D::dualMaxAlgo3(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  double a = (cumsum[t] - cumsum[s]) / (t - s);
  double b = (cumsum[s] - cumsum[r]) / (s - r);
  double c = (costRecord[s] - costRecord[r]) / (s - r);
  double d = (costRecord[s] - minCost) / (t - s);
  double mu_max = muMax(a, b);

  ///
  /// case dual domain = one point OR dual = linear function
  ///
  if(isSpecialCase(a,b) == true){return specialCasePruning(a,b,c,d,mu_max);}

  ///
  /// PELT TEST
  ///
  double nonLinear = Dstar(a);
  double test_value = - nonLinear + d; // dual value in 0
  if (test_value > 0) {return true;} // PELT test


  double meanGap = a - b;
  double grad = nonLinear - meanGap * DstarPrime(a) + c;
  if (test_value + mu_max * grad <= 0) {return false;} //check value in mu = mu_max

  /// iterative algo
  /// iterative algo
  /// iterative algo
  double lt = 0; // left interval value
  double rt = mu_max; // right interval value
  double mu = 0.5 * mu_max; // center interval value
  double m;

  for (int i = 0; i < nb_Loops; ++i)
  {
    m = (a - mu * b) / (1 - mu);
    test_value = - (1 - mu) * Dstar(m) + mu * c + d; // test value
    if (test_value > 0) {return true;} // pruning
    grad = Dstar(m) - (a - b) / (1 - mu) *DstarPrime(m) + c; // gradiant value
    if(grad > 0)
    {
      if (test_value + (rt - mu) * grad <= 0) {return false;} //check value in rt
      lt = mu; // left interval value
      mu = mu + 0.5 * (rt - lt); // center interval value
    }
    else
    {
      if (test_value + (lt - mu) * grad <= 0) {return false;} //check value in lt
      rt = mu; // left interval value
      mu = lt + 0.5 * (rt - lt); // center interval value
    }
  }
  return false;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


bool DUST_1D::dualMaxAlgo4(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  ///
  /// TEST PELT in MU = 0
  ///
  double objectiveMean = (cumsum[t] - cumsum[s]) / (t - s);
  double constantTerm = (costRecord[s] - minCost) / (t - s);
  double nonLinear = Dstar(objectiveMean);
  double test_value = - nonLinear + constantTerm;

  if (test_value > 0) {return true;} // PELT test

  ///
  /// mu_max (stop ALGO if mu_max == 0) !!!
  ///
  double constraintMean = (cumsum[s] - cumsum[r]) / (s - r);
  double mu_max = muMax(objectiveMean, constraintMean);
  if (mu_max == 0) {return false;}

  ///
  /// INIT Quasi Newton algo in zero
  ///
  double linearTerm = (costRecord[s] - costRecord[r]) / (s - r);
  double meanGap = objectiveMean - constraintMean;
  double grad = nonLinear - meanGap * DstarPrime(objectiveMean) + linearTerm;

  /// if grad < 0 the max is in zero, no pelt pruning at this step, return false
  if (grad < 0) {return false;}

  // the duality function is concave, meaning we can check if the tangent at mu ever reaches the desired value.
  if (test_value + mu_max * grad <= 0) {return false;}
  double mu = 0;
  double direction;
  double mu_diff;
  double m_value;
  double grad_diff = -grad; // stores grad difference between two steps, initialized each step as g_k, then once g_{k+1} is computed, y += g_{k+1}
  double inverseHessian = -1;

  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  auto updateDirection = [&] () // update and clip direction
  {
    direction = - inverseHessian * grad;
    if (mu + direction > mu_max) { direction = mu_max - 1e-9 - mu; } // clip ascent direction, as dual may increase beyond the bound
    else if (mu + direction < 0) { direction = -mu + 1e-9; }
  };

  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  auto getM = [&] (double point) // for use in Dstar, DstarPrime
  {
    return (objectiveMean - point * constraintMean) / (1 - point);
  };

  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  auto updateTestValue = [&] ()
  {
    double gradCondition = m1 * grad;

    // Initialize all values
    mu_diff = direction;
    mu += mu_diff;
    m_value = getM(mu);
    nonLinear = Dstar(m_value);
    double new_test = - (1 - mu) * nonLinear + mu * linearTerm + constantTerm;

    int i = 0;
    while(new_test < test_value + mu_diff * gradCondition)
    {
      mu_diff *= .5; // shrink if unsuitable stepsize
      mu -= mu_diff; // relay shrinking
      m_value = getM(mu); // update values
      nonLinear = Dstar(m_value); // update values
      new_test = - (1 - mu) * nonLinear + mu * linearTerm + constantTerm; // update values
      i++;
      if (i == 10) { break; }
    }
    test_value = new_test;
  };

  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  auto updateGrad = [&] ()
  {
    grad = nonLinear - meanGap * DstarPrime(m_value) / (1 - mu) + linearTerm;
  };

  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  auto updateHessian = [&] () // uses BFGS formula (in 1D, the formula simplifies to s / y)
  {
    grad_diff += grad;
    if(grad_diff == 0) { grad_diff = 1e-9; }

    inverseHessian = mu_diff / grad_diff;
    grad_diff = - grad; // setting up y for next step
  };

  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  int i = 0;
  do
  {
    updateDirection();
    updateTestValue();
    if(test_value > 0) {return true;} // index s is pruned
    updateGrad();

    if (grad > 0) // the duality function is concave, meaning we can check if the tangent at mu ever reaches the desired value.
    {
      if (test_value + (mu_max - mu) * grad <= 0) {return false;}
    }
    else if (test_value - mu * grad <= 0) {return false;}
    updateHessian();
    i++;
  } while (i < nb_Loops);
  return false;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


bool DUST_1D::dualMaxAlgo5(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  double objectiveMean = (cumsum[t] - cumsum[s]) / (t - s);
  double constantTerm = (costRecord[s] - minCost) / (t - s);
  double nonLinear = Dstar(objectiveMean);
  double test_value = - nonLinear + constantTerm;

  if (test_value > 0) {return true;} // PELT test
  return false;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


bool DUST_1D::dualMaxAlgo6(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  //Rcout << "START START START" << std::endl;

  // 1. Get model string from get_info()
  Rcpp::List info = get_info();
  std::string model = Rcpp::as<std::string>(info["model"]);
  std::string filename = "dataset_1D_" + model + ".csv";
  std::ofstream file(filename, std::ios::app);

  if (file.tellp() == 0)
  {
    file << "r,s,t,objectiveMean,constraintMean,linearTerm,constantTerm,mu,muMax,pruning,Qt,S1t,ab\n"; // Replace with actual column names
  }

  ///
  ///
  ///
  double constantTerm = (costRecord[s] - minCost) / (t - s);
  double linearTerm = (costRecord[s] - costRecord[r]) / (s - r);
  double nonLinear;

  double objectiveMean = (cumsum[t] - cumsum[s]) / (t - s);
  double constraintMean = (cumsum[s] - cumsum[r]) / (s - r);
  double mu_max = muMax(objectiveMean, constraintMean);

  double test_value;
  // Emmeline: ON AJOUTE
  file << r << "," << s << "," << t << "," << objectiveMean << "," << constraintMean << "," << linearTerm << "," << constantTerm << ",";  // Emmeline

  //// SPECTIAL CASES
  //// SPECTIAL CASES
  //// SPECTIAL CASES
  if(isBoundary(objectiveMean) == true)
  {
    if(isBoundary(constraintMean) == false) /// case mu_max = 0
    {
      /// PELT PELT PELT
      /// PELT PELT PELT
      /// PELT PELT PELT
      file << 0 << ","<<  mu_max << "," << (constantTerm > 0) << "," << minCost <<  "," << cumsum[t] << "," << "f-nf" <<"\n";  // Emmeline
      file.close();  // Emmeline
      if (constantTerm > 0) {return true;}
      return false;
    }
    else /// case linear
    {
      /// linear
      /// linear
      /// linear
      linearTerm = (costRecord[s] - costRecord[r]) / (s - r);
      if (constantTerm > 0)
      {
        file << 0 << ","<<  mu_max << "," << true << "," << minCost <<  "," << cumsum[t] << "," << "f-f" <<"\n";  // Emmeline
        file.close();  // Emmeline
        return true;
      }
      /// what is mu_max here ? verify
      if (mu_max * linearTerm + constantTerm > 0)
      {
        file << mu_max << ","<<  mu_max << "," << true << "," << minCost <<  "," << cumsum[t] << "," << "f-f" <<"\n";  // Emmeline
        file.close();  // Emmeline
        return true;
      }
      file << -1 << ","<<  mu_max << "," << false << "," << minCost <<  "," << cumsum[t] << "," << "f-f" <<"\n";  // Emmeline
      file.close();  // Emmeline
      return false;
    }
  }
  if(objectiveMean == constraintMean)
  {
    /// linear test with Dstar != 0
    /// linear test with Dstar != 0
    /// linear test with Dstar != 0
    linearTerm = (costRecord[s] - costRecord[r]) / (s - r);
    nonLinear = Dstar(objectiveMean);
    if (- nonLinear + constantTerm > 0)
    {
      file << 0 << ","<<  mu_max << "," << true << "," << minCost <<  "," << cumsum[t] << "," << "a=b" <<"\n";  // Emmeline
      file.close();  // Emmeline
      return true;
    }
    if (-(1 - mu_max) * nonLinear + mu_max * linearTerm + constantTerm > 0)
    {
      file << mu_max << ","<<  mu_max << "," << true << "," << minCost <<  "," << cumsum[t] << "," << "a=b" <<"\n";  // Emmeline
      file.close();  // Emmeline
      return true;
    }
    file << -1 << ","<<  mu_max << "," << false << "," << minCost <<  "," << cumsum[t] << "," << "a=b" <<"\n";  // Emmeline
    file.close();  // Emmeline
    return false;
  }

  /// case objectiveMean != Boundary
  /// PELT PELT PELT
  /// PELT PELT PELT
  /// PELT PELT PELT
  nonLinear = Dstar(objectiveMean);
  test_value = - nonLinear + constantTerm;
  if (test_value > 0)
  {
    file << 0 << ","<<  mu_max << "," << true << "," << minCost <<  "," << cumsum[t] << "," << "a!=b" <<"\n";  // Emmeline
    file.close();  // Emmeline
    return true;
  }

  ///
  /// INIT Quasi Newton algo in zero
  ///
  linearTerm = (costRecord[s] - costRecord[r]) / (s - r);
  double meanGap = objectiveMean - constraintMean;
  double grad = nonLinear - meanGap * DstarPrime(objectiveMean) + linearTerm;

  //Rcout << "PROBLEM : " << objectiveMean << std::endl;
  //Rcout << "constraintMean : " << constraintMean << std::endl;
  //Rcout << "mu_max : " << mu_max << std::endl;


  // the duality function is concave, meaning we can check if the tangent at mu ever reaches the desired value.
  double mu = 0;
  double direction;
  double mu_diff;
  double m_value;
  double grad_diff = -grad; // stores grad difference between two steps, initialized each step as g_k, then once g_{k+1} is computed, y += g_{k+1}
  double inverseHessian = -1;

  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  auto updateDirection = [&] () // update and clip direction
  {
    direction = - inverseHessian * grad;
    if (mu + direction > mu_max) { direction = mu_max - 1e-9 - mu; } // clip ascent direction, as dual may increase beyond the bound
    else if (mu + direction < 0) { direction = -mu + 1e-9; }
  };

  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  auto getM = [&] (double point) // for use in Dstar, DstarPrime
  {
    return (objectiveMean - point * constraintMean) / (1 - point);
  };

  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  auto updateTestValue = [&] ()
  {
    double gradCondition = m1 * grad;

    // Initialize all values
    mu_diff = direction;
    mu += mu_diff;
    m_value = getM(mu);
    nonLinear = Dstar(m_value);
    double new_test = - (1 - mu) * nonLinear + mu * linearTerm + constantTerm;

    int i = 0;
    while(new_test < test_value + mu_diff * gradCondition)
    {
      mu_diff *= .5; // shrink if unsuitable stepsize
      mu -= mu_diff; // relay shrinking
      m_value = getM(mu); // update values
      nonLinear = Dstar(m_value); // update values
      new_test = - (1 - mu) * nonLinear + mu * linearTerm + constantTerm; // update values
      i++;
      if (i == 10) { break; }
    }
    test_value = new_test;
  };


  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  auto updateGrad = [&] ()
  {
    grad = nonLinear - meanGap * DstarPrime(m_value) / (1 - mu) + linearTerm;
  };

  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  auto updateHessian = [&] () // uses BFGS formula (in 1D, the formula simplifies to s / y)
  {
    grad_diff += grad;
    if(grad_diff == 0) { grad_diff = 1e-9; }

    inverseHessian = mu_diff / grad_diff;
    grad_diff = - grad; // setting up y for next step
  };

  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  int i = 0;
  //Rcout << "start" << std::endl;
  //Rcout << mu_max << std::endl;
  do
  {
    //Rcout << mu << std::endl;
    updateDirection();
    updateTestValue();
    updateGrad();
    updateHessian();
    i++;
  } while (i < 100);

  file << mu << ","<<  mu_max << "," << (test_value > 0) << "," << minCost <<  "," << cumsum[t] << "," << "a!=b" <<"\n";  // Emmeline
  file.close();  // Emmeline


  if(test_value > 0) {return true;} else {return false;}
}




