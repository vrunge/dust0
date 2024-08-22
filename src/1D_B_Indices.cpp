#include <Rcpp.h>

#include "1D_B_Indices.h"

using namespace Rcpp;

// --------------------------- //
// --- /////////////////// --- //
// --- // Parent Module // --- //
// --- /////////////////// --- //
// --------------------------- //

Indices::~Indices() {}

void Indices::reset()
{
  current = list.begin();
}

void Indices::next()
{
  ++current;
}

bool Indices::check()
{
  return current != list.end();
}

int Indices::get_current()
{
  return *current;
}

std::forward_list<int> Indices::get_list()
{
  return list;
}


// --------------------------- //
// --- /////////////////// --- //
// --- // Random Module // --- //
// --- /////////////////// --- //
// --------------------------- //

RandomIndices::RandomIndices(int size, double alpha) {
  // --- // Generate pseudo-random vector // --- //
  double k = std::max(2., ceil(pow(size, .2)));
  randomU = Rcpp::runif(log(alpha) / log(1 - 1/k));
  u = randomU.begin();
}

void RandomIndices::add(int value)
{
  list.push_front(value);
  pointers.push_back(&list.front());
  nb++;
}

void RandomIndices::reset_prune()
{
  current = list.begin();
  before = list.before_begin();
  pointersCurrent = pointers.rbegin();
  
  nbC = nb - 1;
}

void RandomIndices::next_prune()
{
  before = current;
  ++current;
  ++pointersCurrent;
  new_constraint();
}

void RandomIndices::prune_current()
{
  current = list.erase_after(before);
  pointersCurrent = std::vector<int*>::reverse_iterator(pointers.erase(std::next(pointersCurrent).base()));
  nb--;
  new_constraint();
}

// --- // If no constraint can be selected, exit loop // --- //
bool RandomIndices::check_prune()
{
  return nbC > 0;
}

void RandomIndices::prune_last()
{
  list.erase_after(before);
  pointers.erase(std::next(pointersCurrent).base());
  nb--;
}

// --- // Select new constraint // --- //
// Optimisation possible, car dans le cas random on sélectionne une nouvelle contrainte avant de vérifier qu'elle sera utilisée
void RandomIndices::new_constraint()
{
  nbC--;
}

int RandomIndices::get_constraint()
{
  constraint = pointers[floor(nbC * (*u))];
  
  ++u;
  if (u == randomU.end())
  {
    u = randomU.begin();
  }
  
  return *constraint;
}


// ---------------------------------- //
// --- ////////////////////////// --- //
// --- // Deterministic Module // --- //
// --- ////////////////////////// --- //
// ---------------------------------- //

void DeterministicIndices::add(int value)
{
  list.push_front(value);
}

void DeterministicIndices::reset_prune()
{
  before = list.before_begin();
  current = std::next(before);
  constraint = std::next(current);
}

void DeterministicIndices::next_prune()
{
  before = current;
  current = constraint;
  new_constraint();
}

void DeterministicIndices::prune_current()
{
  current = list.erase_after(before);
  new_constraint();
}

bool DeterministicIndices::check_prune()
{
  return constraint != list.end();
}

void DeterministicIndices::prune_last()
{
  list.erase_after(before);
}

void DeterministicIndices::new_constraint()
{
  ++constraint;
}

int DeterministicIndices::get_constraint()
{
  return *constraint;
}
