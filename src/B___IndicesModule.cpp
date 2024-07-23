// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "B___IndicesModule.h"

using namespace Rcpp;

// --------------------------- //
// --- /////////////////// --- //
// --- // Parent Module // --- //
// --- /////////////////// --- //
// --------------------------- //

IndicesModule::~IndicesModule() {}

void IndicesModule::reset()
{
  current = list.begin();
}

void IndicesModule::next()
{
  ++current;
}

bool IndicesModule::check()
{
  return current != list.end();
}

int IndicesModule::get_current()
{
  return *current;
}

std::forward_list<int> IndicesModule::get_list()
{
  return list;
}


// --------------------------- //
// --- /////////////////// --- //
// --- // Random Module // --- //
// --- /////////////////// --- //
// --------------------------- //

I_Random::I_Random(int size, double alpha) {
  // --- // Generate pseudo-random vector // --- //
  double k = std::max(2., ceil(pow(size, .2)));
  randomU = Rcpp::runif(log(alpha) / log(1 - 1/k));
  u = randomU.begin();
}

void I_Random::add(int value)
{
  list.push_front(value);
  pointers.push_back(&list.front());
  nb++;
}

void I_Random::reset_prune()
{
  current = list.begin();
  before = list.before_begin();
  pointersCurrent = pointers.rbegin();
  
  nbC = nb;
  new_constraint();
}

void I_Random::next_prune()
{
  before = current;
  ++current;
  ++pointersCurrent;
}

// --- // If no constraint can be selected, exit loop // --- //
bool I_Random::check_prune()
{
  return nbC > 0;
}

void I_Random::prune_current()
{
  current = list.erase_after(before);
  pointersCurrent = std::vector<int*>::reverse_iterator(pointers.erase(std::next(pointersCurrent).base()));
  nb--;
}

// --- // Select new constraint // --- //
// Optimisation possible, car dans le cas random on sélectionne une nouvelle contrainte avant de vérifier qu'elle sera utilisée
void I_Random::new_constraint()
{
  nbC--;
  constraint = pointers[floor(nbC * (*u))];
  
  ++u;
  if (u == randomU.end())
  {
    u = randomU.begin();
  }
}

int I_Random::get_constraint()
{
  return *constraint;
}


// ---------------------------------- //
// --- ////////////////////////// --- //
// --- // Deterministic Module // --- //
// --- ////////////////////////// --- //
// ---------------------------------- //

void I_Deterministic::add(int value)
{
  list.push_front(value);
}

void I_Deterministic::reset_prune()
{
  before = list.before_begin();
  current = std::next(before);
  constraint = std::next(current);
}

void I_Deterministic::next_prune()
{
  before = current;
  current = constraint;
}

bool I_Deterministic::check_prune()
{
  return constraint != list.end();
}

void I_Deterministic::prune_current()
{
  current = list.erase_after(before);
}

void I_Deterministic::new_constraint()
{
  ++constraint;
}

int I_Deterministic::get_constraint()
{
  return *constraint;
}
