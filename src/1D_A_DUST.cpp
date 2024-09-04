#include <Rcpp.h>
#include <cmath>

#include <random> /// FOR RANDOM NUMBER IN DUAL EVAL

#include "1D_A_DUST.h"
#include "preProcessing.h"

using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// --- // Constructor // --- //
DUST_1D::DUST_1D(int dual_max, bool random_constraint, Nullable<double> alpha_, Nullable<int> nbLoops)
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

DUST_1D::~DUST_1D()
{
  delete indices;
}

void DUST_1D::init_method()
{
  delete indices;

  /// /// ///
  /// /// /// index METHOD
  /// /// ///
  if(random_constraint)
  {
    indices = new RandomIndices(n, alpha);
  }
  else
  {
    indices = new DeterministicIndices;
  }

  /// /// ///
  /// /// /// dual_max METHOD
  /// /// ///
  if(dual_max == 0)
  {
    current_test = &DUST_1D::dualMaxAlgo0;
  }
  if(dual_max == 1)
  {
    current_test = &DUST_1D::dualMaxAlgo1;
  }
  if(dual_max == 2)
  {
    current_test = &DUST_1D::dualMaxAlgo2;
  }
  if(dual_max == 3)
  {
    current_test = &DUST_1D::dualMaxAlgo3;
  }
  if(dual_max == 4)
  {
    current_test = &DUST_1D::dualMaxAlgo4;
  }
  if(dual_max == 5)
  {
    current_test = &DUST_1D::dualMaxAlgo5;
  }
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double DUST_1D::dualMaxAlgo0(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  return dualEval(dist(engine), minCost, t, s, r);
}

double DUST_1D::dualMaxAlgo1(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  return dualMax(minCost, t, s, r);
}

double DUST_1D::dualMaxAlgo2(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  double a = 0.0;
  double b = 1.0;
  double c = 1 - 1/phi;
  double d = 1/phi;

  double fc = dualEval(c, minCost, t, s, r);
  double fd = dualEval(d, minCost, t, s, r);
  double max_val = std::max(fc, fd);

  for (int i = 0; i < nb_Loops; i++)
  {
    if (fc > fd)
    {
      b = d;
      d = c;
      fd = fc;
      c = b - (b - a) / phi;
      fc = dualEval(c, minCost, t, s, r);
    }
    else
    {
      a = c;
      c = d;
      fc = fd;
      d = a + (b - a) / phi;
      fd = dualEval(d, minCost, t, s, r);
    }
    max_val = std::max(max_val, std::max(fc, fd));
    if(max_val > 0){break;}
  }
  return max_val;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double DUST_1D::dualMaxAlgo3(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  double a = (cumsum[t] - cumsum[s]) / (t - s); // m_it
  double b = (cumsum[s] - cumsum[r]) / (s - r); // m_ji
  double C = (costRecord[s] - costRecord[r]) / (s - r);

  ////DANGER
  double mu = 0.5*std::min(1.0, a/b);
  double m = (a - mu*b) / (1.0 - mu);
  return -(1.0 - mu) * Dstar(m) + mu * C - (minCost - costRecord[s]) / (t - s);
}

double DUST_1D::dualMaxAlgo4(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  return (-std::numeric_limits<double>::infinity());
}

double DUST_1D::dualMaxAlgo5(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  return (-std::numeric_limits<double>::infinity());
}


/*
double DUST_1D::dualMaxAlgo3(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  double a = (cumsum[t] - cumsum[s]) / (t - s); // m_it
  double b = (cumsum[s] - cumsum[r]) / (s - r); // m_ji
  double C = (costRecord[s] - costRecord[r]) / (s - r);

  //std::cout << "c " << C << std::endl;
  double f_prime;
  double f_second;
  double mu_new;

  ////DANGER
  double mu = 0.01;

  double m;
  for (int i = 0; i < nb_Loops; ++i)
  {

    m =  (a - mu*b) / (1.0 - mu);
    //std::cout << "dualEvaldualEval" << -(1.0 - mu) * Dstar(m) + mu * C - (minCost - costRecord[s]) / (t - s) << std::endl;

    f_prime = Dstar(m) - ((a-b)/(1.0 - mu)) * DstarPrime(m) + C;
    f_second = - (std::pow(a-b,2)/std::pow(1.0-mu,3)) * DstarSecond(m);

    mu_new = mu - (f_prime / f_second);
    mu_new = std::max(0.0, std::min(1.0, mu_new));
    // Check for pruning
    //if(dualEval(mu_new, minCost, t, s, r) > 0) {break;}
    //std::cout << "mu " << mu << " munew " << mu_new << " ab " << a << " " << b << std::endl;
    mu = mu_new;
    //std::cout <<  " Dstar(m)" << Dstar(m)  << " DstarPrime(m) " << ((a-b)/(1.0 - mu)) * DstarPrime(m)  << " DstarSecond(m) " << DstarSecond(m) << std::endl;
    //std::cout <<  " --- " << f_prime  << " --- " << f_second << " +++ " << f_prime / f_second << std::endl;
  }
  //std::cout << std::endl;
  return -(1.0 - mu) * Dstar(m) + mu * C - (minCost - costRecord[s]) / (t - s);
  //return(-std::numeric_limits<double>::infinity());
}

*/


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// --- // Fits the data, i. e. initializes all data-dependent vectors // --- //
void DUST_1D::init(std::vector<double>& inData, Nullable<double> inPenalty)
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
  costRecord = std::vector<double>(n + 1, -penalty);


  init_method();

  indices->add(0);
  indices->add(1);

}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// --- // Algorithm-specific method // --- //
void DUST_1D::compute(std::vector<double>& inData)
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
  cumsum[1] = inData[0];
  costRecord[1] = Cost(t, s);
  changepointRecord[1] = 0;

  int nbt = 2;
  nb_indices[0] = 1;

  // Main loop
  for (t = 2; t <= n; t++)
  {
    // update cumsum
    cumsum[t] =
      cumsum[t - 1] + inData[t - 1];

    // OP step
    indices->reset();
    minCost = std::numeric_limits<double>::infinity();
    do
    {
      s = indices->get_current();
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
      if ((this->*current_test)(minCost, t, indices->get_current(), indices->get_constraint()) > 0) // prune as needs pruning
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
    if (lastCost > minCost)
    {
      indices->prune_last();
      nbt--;
    }

    // update the available indices
    indices->add(t);
    nbt++;
    nb_indices[t - 1] = nbt;
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
List DUST_1D::get_partition()
{
  costRecord.erase(costRecord.begin()); ///// REMOVE FIRST ELEMENT /////
  return List::create(
    _["changepoints"] = backtrack_changepoints(),
    _["lastIndexSet"] = indices->get_list(),
    _["nb"] = nb_indices,
    _["costQ"] = costRecord
  );
}

// --- // Wrapper method for quickly computing               // --- //
// --- // and retrieving the optimal partition of input data // --- //
List DUST_1D::quick(std::vector<double>& inData, Nullable<double> inPenalty)
{
  init(inData, inPenalty);
  compute(inData);
  return get_partition();
}
