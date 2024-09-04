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
DUST_meanVar::DUST_meanVar(int dual_max, bool random_constraint, Nullable<double> alpha_, Nullable<int> nbLoops)
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
    current_test = &DUST_meanVar::dualMaxAlgo0;
  }
  if(dual_max == 1)
  {
    current_test = &DUST_meanVar::dualMaxAlgo1;
  }
  if(dual_max == 2)
  {
    current_test = &DUST_meanVar::dualMaxAlgo2;
  }
  if(dual_max == 3)
  {
    current_test = &DUST_meanVar::dualMaxAlgo3;
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

double DUST_meanVar::dualMaxAlgo0(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  return dualEval(dist(engine), minCost, t, s, r);
}

double DUST_meanVar::dualMaxAlgo1(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  return (-std::numeric_limits<double>::infinity());
}

double DUST_meanVar::dualMaxAlgo2(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  if(s + 1 == t){return(-std::numeric_limits<double>::infinity());}
  if(r + 1 == s){return(-std::numeric_limits<double>::infinity());}

  double a = 0.0;
  double b = 1.0;
  double c = 1 - 1/phi;
  double d = 1/phi;

  double fc = DUST_meanVar::dualEval(c, minCost, t, s, r);
  double fd = DUST_meanVar::dualEval(d, minCost, t, s, r);
  double max_val = std::max(fc, fd);

  for (int i = 0; i < nb_Loops; i++)
  {
    if (fc > fd)
    {
      b = d;
      d = c;
      fd = fc;
      c = b - (b - a) / phi;
      fc = DUST_meanVar::dualEval(c, minCost, t, s, r);
    }
    else
    {
      a = c;
      c = d;
      fc = fd;
      d = a + (b - a) / phi;
      fd = DUST_meanVar::dualEval(d, minCost, t, s, r);
    }
    max_val = std::max(max_val, std::max(fc, fd));
    if(max_val > 0){break;}
  }
  return max_val;
}



double DUST_meanVar::dualMaxAlgo3(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  return (-std::numeric_limits<double>::infinity());
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
  unsigned int argMin; // stores the optimal last changepoint for the current OP step

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
      // TO DO
      // TO DO
      // TO DO
      // Discuss with Simon. What it does exactly?
      //indices->prune_last();
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
  return List::create(
    _["changepoints"] = backtrack_changepoints(),
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

  if(s + 1 == t){return(std::numeric_limits<double>::infinity());}
  double m = (cumsum[t] - cumsum[s]) / (t - s);
  return 0.5 * (t - s) * (1 + std::log((cumsum2[t] - cumsum2[s]) / (t - s) - m * m));
}


double DUST_meanVar::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  if(s + 1 == t){return(-std::numeric_limits<double>::infinity());}
  if(r + 1 == s){return(-std::numeric_limits<double>::infinity());}
  double Mt = (cumsum[t] - cumsum[s]) / (t - s);
  double Mt2 = (cumsum2[t] - cumsum2[s]) / (t - s);
  double Ms = (cumsum[s] - cumsum[r]) / (s - r);
  double Ms2 = (cumsum2[s] - cumsum2[r]) / (s - r);

  // Compute variance terms
  double Va = Mt2 - std::pow(Mt, 2);
  double Vb = Ms2 - std::pow(Ms, 2);

  double u = (Va + Vb) * (1 + std::pow((Mt - Ms) / std::sqrt(Va + Vb), 2));
  point = point * (u - std::sqrt(std::pow(u, 2) - 4.0 * Va * Vb)) / (2.0 * Vb);

  //std::cout << point << " ";
  double A = (Mt2 - point *  Ms2)/(1 - point);
  double B = (Mt - point *  Ms)/(1 - point);

  return (costRecord[s] - minCost) / (t - s)
  + point * (costRecord[s] - costRecord[r]) / (s - r)
  + 0.5 * (1 - point) * (1 + std::log(A - B*B));
}


/////////////////////////////////////////////////////////////

double DUST_meanVar::Dstar(double x) const
{
  return 0;
}

double DUST_meanVar::DstarPrime(double x) const
{
  return 0;
}

double DUST_meanVar::DstarSecond(double x) const
{
  return 0;
}




