#include <Rcpp.h>
#include <cmath>

#include <random> /// FOR RANDOM NUMBER IN DUAL EVAL
#include <limits>

#include "2D_DUSTmeanVar.h"
#include "preProcessing.h"

using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// --- // Constructor // --- //
DUST_meanVar::DUST_meanVar(int dual_max_type, int constraint_indices, Nullable<int> nbLoops)
  : dual_max_type(dual_max_type),
    constraint_indices(constraint_indices),
    indices(nullptr)
{
  if(nbLoops.isNull()){nb_Loops = 10;}else{nb_Loops = as<int>(nbLoops);}
}

DUST_meanVar::~DUST_meanVar()
{
  delete indices;
}

void DUST_meanVar::init_method()
{
  delete indices;

  /// /// ///
  /// /// /// index METHOD
  /// /// ///
  if(constraint_indices == 10){indices = new RandomIndices_2D;} ///(n, alpha) ??? only in 1D ?
  if(constraint_indices == 20){indices = new RandomIndices_2D2;}

  if(constraint_indices == 11){indices = new DeterministicIndices_2D;}
  if(constraint_indices == 21){indices = new DeterministicIndices_2D2;}

  /// /// ///
  /// /// /// dual_max_type METHOD
  /// /// ///
  if(dual_max_type == 0){current_test = &DUST_meanVar::dualMaxAlgo0;}
  if(dual_max_type == 1){current_test = &DUST_meanVar::dualMaxAlgo1;}
  if(dual_max_type == 2){current_test = &DUST_meanVar::dualMaxAlgo2;}
  if(dual_max_type == 3){current_test = &DUST_meanVar::dualMaxAlgo3;}
  if(dual_max_type == 4){current_test = &DUST_meanVar::dualMaxAlgo4;}
  if(dual_max_type == 5){current_test = &DUST_meanVar::dualMaxAlgo5;}

  /// /// ///
  /// /// /// INIT RANDOM GENERATOR
  /// /// ///
  engine.seed(std::random_device{}());
  dist = std::uniform_real_distribution<double>(0.0, 1.0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool DUST_meanVar::dualMaxAlgo0(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  if(s + 1 == t){return false;}
  //if(r + 1 == s){return false;} // => Vb = 0
  return (dualEval(dist(engine), minCost, t, s, r) > 0);
}

bool DUST_meanVar::dualMaxAlgo1(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  return false;
}

bool DUST_meanVar::dualMaxAlgo2(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  if(s + 1 == t){return false;}
 // if(r + 1 == s){return false;} // => Vb = 0
  double Mt = (cumsum[t] - cumsum[s]) / (t - s);
  double Mt2 = (cumsum2[t] - cumsum2[s]) / (t - s);
  double Ms = (cumsum[s] - cumsum[r]) / (s - r);
  double Ms2 = (cumsum2[s] - cumsum2[r]) / (s - r);

  double lt = 0.0;
  double rt = muMax(Mt, Ms, Mt2, Ms2);
  double c = (1 - 1/phi) * rt;
  double d = 1/phi;

  double linear = (costRecord[s] - costRecord[r])/(s - r);
  double cst = (costRecord[s] - minCost)/(t - s);

  double A = (Mt2 - c *  Ms2)/(1 - c);
  double B = (Mt - c *  Ms)/(1 - c);
  double fc =  0.5 * (1 - c) * (1 + std::log(A - B*B)) + c * linear + cst;

  A = (Mt2 - d * Ms2)/(1 - d);
  B = (Mt - d * Ms)/(1 - d);
  double fd = 0.5 * (1 - d) * (1 + std::log(A - B*B)) + d * linear + cst;
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
      A = (Mt2 - c * Ms2)/(1 - c);
      B = (Mt - c * Ms)/(1 - c);
      fc =  0.5 * (1 - c) * (1 + std::log(A - B*B)) + c * linear + cst;
    }
    else
    {
      lt = c;
      c = d;
      fc = fd;
      d = lt + (rt - lt) / phi;
      A = (Mt2 - d * Ms2)/(1 - d);
      B = (Mt - d * Ms)/(1 - d);
      fd = 0.5 * (1 - d) * (1 + std::log(A - B*B)) + d * linear + cst;
    }
    max_val = std::max(max_val, std::max(fc, fd));
    if(max_val > 0){return true;}
  }
  return false;
}



bool DUST_meanVar::dualMaxAlgo3(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  return false;
}



bool DUST_meanVar::dualMaxAlgo4(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  if(s + 1 == t){return false;}
  //if(r + 1 == s){return false;} // => Vb = 0

  double a = (cumsum[t] - cumsum[s]) / (t - s);
  double a2 = (cumsum2[t] - cumsum2[s]) / (t - s);

  double constantTerm = (costRecord[s] - minCost) / (t - s);
  double nonLinear = 0.5 * (1 + std::log(a2 - a*a));
  double test_value = nonLinear + constantTerm;  //dual in mu = 0

  if (test_value > 0) {return true;} // PELT test (eval dual in 0)

  /// /// ///

  double b = (cumsum[s] - cumsum[r]) / (s - r);
  double b2 = (cumsum2[s] - cumsum2[r]) / (s - r);

  double linearTerm = (costRecord[s] - costRecord[r]) / (s - r);
  double term = a2 - b2 - 2 * a * (a - b);

  double grad = - nonLinear + term *  0.5 / (a2 - a*a) + linearTerm; //dual prime in mu = 0
  double mu_max = muMax(a, b, a2, b2);

  /////////
  /////////

  // the duality function is concave, meaning we can check if the tangent at mu ever reaches the desired value.
  if (test_value + mu_max * grad <= 0) {return false;}
  double mu = 0;
  double direction;
  double mu_diff;
  double m_value;
  double m_value2;
  double grad_diff = -grad; // stores grad difference between two steps, initialized each step as g_k, then once g_{k+1} is computed, y += g_{k+1}
  double inverseHessian = -1;

  auto updateDirection = [&] () // update and clip direction
  {
    direction = - inverseHessian * grad;
    if (mu + direction > mu_max) { direction = mu_max - 1e-9 - mu; } // clip ascent direction, as dual may increase beyond the bound
    else if (mu + direction < 0) { direction = -mu + 1e-9; }
  };

  auto updateTestValue = [&] ()
  {
    double gradCondition = m1 * grad;
    // Initialize all values
    mu_diff = direction;
    mu += mu_diff;
    m_value = pow(1 - mu, -1) * (a - mu * b);
    m_value2 = pow(1 - mu, -1) * (a2 - mu * b2);
    nonLinear = 0.5 * (1 + std::log(m_value2 - m_value*m_value));
    double new_test = (1 - mu) * nonLinear + mu * linearTerm + constantTerm; ///eval dual at mu

    int i = 0;
    while(new_test < test_value + mu_diff * gradCondition)
    {
      mu_diff *= .5; // shrink if unsuitable stepsize
      mu -= mu_diff; // relay shrinking
      m_value = pow(1 - mu, -1) * (a - mu * b);
      m_value2 = pow(1 - mu, -1) * (a2 - mu * b2);
      nonLinear = 0.5 * (1 + std::log(m_value2 - m_value*m_value)); // update values
      new_test = (1 - mu) * nonLinear + mu * linearTerm + constantTerm; // update values
      i++;
      if (i == 10) { break; }
   }
    test_value = new_test;
  };

  auto updateGrad = [&] ()
  {
    grad = - nonLinear +
      ((a2 - b2)*pow(1 - mu, -1) - 2*(a - mu*b)*(a - b)*pow(1 - mu, -2))*  0.5 / (m_value2 - m_value*m_value) + linearTerm;
  };

  auto updateHessian = [&] () // uses BFGS formula (in 1D, the formula simplifies to s / y)
  {
    grad_diff += grad;
    if(grad_diff == 0) { grad_diff = 1e-9; }

    inverseHessian = mu_diff / grad_diff;
    grad_diff = - grad; // setting up y for next step
  };

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
  } while (i < 100);
  return false;
}



bool DUST_meanVar::dualMaxAlgo5(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  if(s + 1 == t){return false;}
  //if(r + 1 == s){return false;} // => Vb = 0

  double a = (cumsum[t] - cumsum[s]) / (t - s);
  double a2 = (cumsum2[t] - cumsum2[s]) / (t - s);

  double constantTerm = (costRecord[s] - minCost) / (t - s);
  double nonLinear = 0.5 * (1 + std::log(a2 - a*a));
  double test_value = nonLinear + constantTerm;  //dual in mu = 0

  if (test_value > 0) {return true;} // PELT test (eval dual in 0)
  return false;
}


bool DUST_meanVar::dualMaxAlgo6(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  return false;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// --- // Fits the data, i. e. initializes all data-dependent vectors // --- //
void DUST_meanVar::init(std::vector<double>& inData, Nullable<double> inPenalty)
{
  n = inData.size();

  if (inPenalty.isNull())
  {
    penalty = 2 * pow(sdDiff(inData), 2) * std::log(n); //to do
  }
  else
  {
    penalty = as<double>(inPenalty);
  }

  changepointRecord = std::vector<int>(n + 1, 0);
  nb_indices = std::vector<int>(n, 0);

  cumsum = std::vector<double>(n + 1, 0.);
  cumsum2 = std::vector<double>(n + 1, 0.);
  costRecord = std::vector<double>(n + 1, -penalty);

  init_method();

  indices->add(0);
  indices->add(1);

}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// --- // Algorithm-specific method // --- //
void DUST_meanVar::compute(std::vector<double>& inData)
{
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
  cumsum2[1] = inData[0] * inData[0];
  costRecord[1] = Cost(t, s);
  changepointRecord[1] = 0;

  int nbt = 2;
  nb_indices[0] = 1;

  // Main loop
  for (t = 2; t <= n; t++)
  {
    // update cumsum and cumsum2
    cumsum[t] = cumsum[t - 1] + inData[t - 1];
    cumsum2[t] = cumsum2[t - 1] + inData[t - 1] * inData[t - 1];

    // OP step
    indices->reset();
    minCost = std::numeric_limits<double>::infinity();
    do
    {
      s = *(indices->current);
      lastCost = costRecord[s] + Cost(t, s);
      if (lastCost <= minCost) ////// <=
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
    while (indices->check())
    {
      if ((this->*current_test)(minCost, t, *(indices->current), indices->get_constraint_l())) // prune as needs pruning
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

    // Prune the last index (analoguous with a null (mu* = 0) duality simple test)
    if (lastCost > minCost)
    {
      // TO DO
      // TO DO
      // TO DO
      // Discuss with Simon. What it does exactly?
      //indices->prune_last();
    }

    // update the available indices
    indices->add(t);
    nb_indices[t - 1] = nbt;
    nbt++;
  }
}



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


// --- // Builds changepoints // --- //
std::forward_list<unsigned int> DUST_meanVar::backtrack_changepoints()
{
  std::forward_list<unsigned int> changepoints {n};
  for (int newChangepoint = changepointRecord[n]; newChangepoint != 0; newChangepoint = changepointRecord[newChangepoint])
  {
    changepoints.push_front(newChangepoint);
  }
  return changepoints;
}


// --- // Retrieves optimal partition // --- //
List DUST_meanVar::get_partition()
{
  costRecord.erase(costRecord.begin()); ///// REMOVE FIRST ELEMENT /////
  indices->remove_last(); ///// REMOVE FIRST ELEMENT /////

  std::forward_list<unsigned int> chpts = backtrack_changepoints();

  return List::create(
    _["changepoints"] = chpts,
    _["lastIndexSet"] = indices->get_list(),
    _["nb"] = nb_indices,
    _["costQ"] = costRecord
  );
}

// --- // Wrapper method for quickly computing               // --- //
// --- // and retrieving the optimal partition of input data // --- //
List DUST_meanVar::quick(std::vector<double>& inData, Nullable<double> inPenalty)
{
  init(inData, inPenalty);
  compute(inData);
  return get_partition();
}


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


double DUST_meanVar::Cost(unsigned int t, unsigned int s) const
{

  if(s + 1 == t){return(std::numeric_limits<double>::infinity());} // infinite cost segment if one data point only
  double m = (cumsum[t] - cumsum[s]) / (t - s);
  double var = (cumsum2[t] - cumsum2[s]) / (t - s) - m * m;
  //if(var <= 0){return(-std::numeric_limits<double>::infinity());}
  return 0.5 * (t - s) * (1 + std::log(var));
}


double DUST_meanVar::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  if(s + 1 == t){return(-std::numeric_limits<double>::infinity());}
  //if(r + 1 == s){return(-std::numeric_limits<double>::infinity());} // => Vb = 0
  double Mt = (cumsum[t] - cumsum[s]) / (t - s);
  double Mt2 = (cumsum2[t] - cumsum2[s]) / (t - s);
  double Ms = (cumsum[s] - cumsum[r]) / (s - r);
  double Ms2 = (cumsum2[s] - cumsum2[r]) / (s - r);

  // Compute variance terms
  double Va = Mt2 - std::pow(Mt, 2);
  double Vb = Ms2 - std::pow(Ms, 2);

  /// pruning if same mean and same variance
  //if(Mt == Ms && Va == Vb){return(std::numeric_limits<double>::infinity());}

  double u = (Va + Vb) * (1 + std::pow((Mt - Ms) / std::sqrt(Va + Vb), 2));

  if(Vb > 0){point = point * ((u - std::sqrt(std::pow(u, 2) - 4.0 * Va * Vb)) / (2.0 * Vb));}
  else{point = point * (Va / (Va + pow(Mt - Ms, 2)));}

  //std::cout << point << " ";
  double A = (Mt2 - point *  Ms2)/(1 - point);
  double B = (Mt - point *  Ms)/(1 - point);

  return (costRecord[s] - minCost) / (t - s)
  + point * (costRecord[s] - costRecord[r]) / (s - r)
  + 0.5 * (1 - point) * (1 + std::log(A - B*B));
}


/////////////////////////////////////////////////////////////

double DUST_meanVar::muMax(double a, double b, double a2, double b2) const
{
  double Va = a2 - std::pow(a, 2);
  double Vb = b2 - std::pow(b, 2);
  double u = (Va + Vb) * (1 + std::pow((a - b) / std::sqrt(Va + Vb), 2));

  if(Vb > 0){return((u - std::sqrt(std::pow(u, 2) - 4.0 * Va * Vb)) / (2.0 * Vb));}
  else{return(Va / (Va + pow(a - b, 2)));}
}








