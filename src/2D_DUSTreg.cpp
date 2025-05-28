#include <Rcpp.h>
#include <cmath>

#include <random> /// FOR RANDOM NUMBER IN DUAL EVAL
#include <limits>

#include "2D_DUSTreg.h"
#include "preProcessing.h"

using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// --- // Constructor // --- //
DUST_reg::DUST_reg(int dual_max_type, int constraint_indices, Nullable<int> nbLoops)
  : dual_max_type(dual_max_type),
    constraint_indices(constraint_indices),
    indices(nullptr)
{
  if(nbLoops.isNull()){nb_Loops = 10;}else{nb_Loops = as<int>(nbLoops);}
}

DUST_reg::~DUST_reg()
{
  delete indices;
}

void DUST_reg::init_method()
{
  delete indices;

  /// /// ///
  /// /// /// index METHOD
  /// /// ///
  if(constraint_indices == 10){indices = new RandomIndices_2D;}
  if(constraint_indices == 20){indices = new RandomIndices_2D2;}

  if(constraint_indices == 11){indices = new DeterministicIndices_2D;}
  if(constraint_indices == 21){indices = new DeterministicIndices_2D2;}


  /// /// ///
  /// /// /// dual_max_type METHOD
  /// /// ///
  if(dual_max_type == 0){current_test = &DUST_reg::dualMaxAlgo0;}
  if(dual_max_type == 1){current_test = &DUST_reg::dualMaxAlgo1;}
  if(dual_max_type == 2){current_test = &DUST_reg::dualMaxAlgo2;}
  if(dual_max_type == 3){current_test = &DUST_reg::dualMaxAlgo3;}
  if(dual_max_type == 4){current_test = &DUST_reg::dualMaxAlgo4;}
  if(dual_max_type == 5){current_test = &DUST_reg::dualMaxAlgo5;}

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

bool DUST_reg::dualMaxAlgo0(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  return (dualEval(dist(engine), minCost, t, s, r) > 0);
}

bool DUST_reg::dualMaxAlgo1(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  return false;
}

bool DUST_reg::dualMaxAlgo2(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  if(s + 1 == t){return false;}
  if(r + 1 == s){return false;}

  double a = 0.0;
  double b = 1.0;
  double c = 1 - 1/phi;
  double d = 1/phi;

  double fc = dualEval(c, minCost, t, s, r);
  double fd = dualEval(d, minCost, t, s, r);
  if(fc > 0 || fd > 0){return true;}
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
    if(max_val > 0){return true;}
  }
  return false;
}



bool DUST_reg::dualMaxAlgo3(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  return false;
}

bool DUST_reg::dualMaxAlgo4(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  return false;
}

bool DUST_reg::dualMaxAlgo5(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  return false;
}

bool DUST_reg::dualMaxAlgo6(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  return false;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// --- // Fits the data, i. e. initializes all data-dependent vectors // --- //
void DUST_reg::init(DataFrame& inData, Nullable<double> inPenalty)
{
  data_x = inData["x"];
  data_y = inData["y"];

  n = data_y.size();
  penalty = as<double>(inPenalty);

  A = std::vector<double>(n + 1, 0.);
  B = std::vector<double>(n + 1, 0.);
  C = std::vector<double>(n + 1, 0.);
  D = std::vector<double>(n + 1, 0.);
  E = std::vector<double>(n + 1, 0.);
  F = std::vector<double>(n + 1, 0.);

  changepointRecord = std::vector<int>(n + 1, 0);
  nb_indices = std::vector<int>(n, 0);

  costRecord = std::vector<double>(n + 1, -penalty);

  init_method();

  indices->add(0);
  indices->add(1);

}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// --- // Algorithm-specific method // --- //
void DUST_reg::compute()
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
  A[1] = data_x[0] * data_x[0];
  B[1] = data_x[0];
  C[1] = 1;
  D[1] = - data_x[0] * data_y[0];
  E[1] = - data_y[0];
  F[1] = data_y[0] * data_y[0];

  int nbt = 2;
  nb_indices[0] = 1;

  costRecord[1] = Cost(t, s);
  changepointRecord[1] = 0;

  // Main loop
  for (t = 2; t <= n; t++)
  {
    // update cumsum and cumsum2
    A[t] = A[t - 1] + data_x[t - 1] * data_x[t - 1];
    B[t] = B[t - 1] + data_x[t - 1];
    C[t] = t;
    D[t] = D[t - 1] - data_x[t - 1] * data_y[t - 1];
    E[t] = E[t - 1] - data_y[t - 1];
    F[t] = F[t - 1] + data_y[t - 1] * data_y[t - 1];


    // OP step
    indices->reset();
    minCost = std::numeric_limits<double>::infinity();
    do
    {
      s = *(indices->current);
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
std::forward_list<unsigned int> DUST_reg::backtrack_changepoints()
{
  std::forward_list<unsigned int> changepoints {n};
  for (int newChangepoint = changepointRecord[n]; newChangepoint != 0; newChangepoint = changepointRecord[newChangepoint])
  {
    changepoints.push_front(newChangepoint);
  }
  return changepoints;
}



// --- // Retrieves optimal partition // --- //
List DUST_reg::get_partition()
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
List DUST_reg::quick(DataFrame& inData, Nullable<double> inPenalty)
{
  init(inData, inPenalty);
  compute();
  return get_partition();
}


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


double DUST_reg::Cost(unsigned int t, unsigned int s) const
{
  if(s + 1 == t){return(std::numeric_limits<double>::infinity());}

  double Adiff = A[t] - A[s];
  double Bdiff = B[t] - B[s];
  double Cdiff = C[t] - C[s];
  double Ddiff = D[t] - D[s];
  double Ediff = E[t] - E[s];
  double Fdiff = F[t] - F[s];

  double num = 2.0 * Bdiff * Ddiff * Ediff - Adiff * Ediff * Ediff - Cdiff * Ddiff * Ddiff;
  double denom = Adiff * Cdiff - Bdiff * Bdiff;

  return num / denom + Fdiff;
}


double DUST_reg::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  if(s + 1 == t){return(-std::numeric_limits<double>::infinity());}
  if(r + 1 == s){return(-std::numeric_limits<double>::infinity());}
  double Mt = (B[t] - B[s]) / (t - s);
  double Mt2 = (A[t] - A[s]) / (t - s);
  double Ms = (B[s] - B[r]) / (s - r);
  double Ms2 = (A[s] - A[r]) / (s - r);

  // Compute variance terms
  double Va = Mt2 - std::pow(Mt, 2);
  double Vb = Ms2 - std::pow(Ms, 2);
  double u = (Va + Vb) * (1 + std::pow((Mt - Ms) / std::sqrt(Va + Vb), 2));
  point = point * ((t - s)/ (s - r)) * (u - std::sqrt(std::pow(u, 2) - 4.0 * Va * Vb)) / (2.0 * Vb);

  double Adiff = A[t] - A[s] - point * (A[s] - A[r]);
  double Bdiff = B[t] - B[s] - point * (B[s] - B[r]);
  double Cdiff = C[t] - C[s] - point * (C[s] - C[r]);
  double Ddiff = D[t] - D[s] - point * (D[s] - D[r]);
  double Ediff = E[t] - E[s] - point * (E[s] - E[r]);
  double Fdiff = F[t] - F[s] - point * (F[s] - F[r]);

  double num = 2.0 * Bdiff * Ddiff * Ediff - Adiff * Ediff * Ediff - Cdiff * Ddiff * Ddiff;
  double denom = Adiff * Cdiff - Bdiff * Bdiff;

  return (costRecord[s] - minCost)
    + point * (costRecord[s] - costRecord[r])
    + num / denom + Fdiff;
}


double DUST_reg::muMax(double a, double b) const
{
  return 0;
}


double DUST_reg::Dstar(double x) const
{
  return 0;
}

double DUST_reg::DstarPrime(double x) const
{
  return 0;
}

double DUST_reg::DstarSecond(double x) const
{
  return 0;
}




